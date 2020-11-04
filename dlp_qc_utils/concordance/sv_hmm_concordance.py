import numpy as np 
from single_cell.utils.csvutils import IrregularCsvInput, CsvOutput, CsvInput
import statistics
import pandas as pd
import pyranges as pr


def get_cna_changepoints(cna_data, change_region_length = 500000):
    cn_data_summary = cna_data.groupby(['chr', 'start', 'end'])
    cn_data_summary = cn_data_summary[['copy', 'state']]
    cn_data_summary = cn_data_summary.aggregate(
        {'copy': statistics.median, 'state': statistics.median}
    )
    cn_data_summary["cn_change"]=cn_data_summary.groupby('chr').state.transform(pd.Series.diff)
    cn_data_summary.reset_index(inplace=True)
    cn_data_summary = cn_data_summary[(cn_data_summary.cn_change!=0)].dropna()
    #pandas diff() only includes the higher side of a place of change: i.e.: 1,1,2,2 --> NaN, 0, changepoint, 0
    #include the lower side to get: 1,1,2,2 --> NaN, changepoint, changepoint, 0.
    #accomplish by expanding the lower boundary of every changepoint by the changepoint size. 
    #shouldnt matter if start becomes < 0, as breakpoints can be < 0.
    cn_data_summary["start"] = cn_data_summary.start - change_region_length
    print(cn_data_summary)

    return cn_data_summary



def anno_bkp(row, cn, side="1"):
    if side == "1":
        chrom = str(row.chromosome_1)
        pos = row.position_1
    else:
        chrom = str(row.chromosome_2)
        pos = row.position_2        
    pos = int(pos)
    matches = cn[(cn.start < pos) 
        & (cn.end > pos) 
        & (cn.chr == chrom)]
    matches = matches.reset_index()
    if matches.empty:
        row["chr_pos_{}".format(side)] = np.NaN
        row["start_pos_{}".format(side)] = np.NaN
        row["end_pos_{}".format(side)] = np.NaN
        row["state_pos_{}".format(side)] = np.NaN
        row["cn_change_pos_{}".format(side)] = np.NaN
    else:
        row["chr_pos_{}".format(side)] = matches.chr.tolist()[0]
        row["start_pos_{}".format(side)] = matches.start.tolist()[0]
        row["end_pos_{}".format(side)] = matches.end.tolist()[0]
        row["state_pos_{}".format(side)] = matches.state.tolist()[0]
        row["cn_change_pos_{}".format(side)] = matches.cn_change.tolist()[0]
    row["pos_{}_in_cn_change_region".format(side)] = not matches.empty
    return row


def match_changepoints(changepoints, sv):
    sv = sv.apply(lambda row: anno_bkp(row, changepoints, side="2"), axis=1)
    sv = sv.reset_index()
    sv = sv.apply(lambda row: anno_bkp(row, changepoints), axis=1)
    return sv


def flatten_breakpoints(sv):
    out = pd.DataFrame({"chromosome":[], "position":[], "type":[], "prediction_id":[]})
    out["chromosome"] = sv.chromosome_1.tolist() + sv.chromosome_2.tolist()
    out["position"] = sv.position_1.tolist() + sv.position_2.tolist()
    out["rearrangement_type"] = sv.type.tolist() + sv.type.tolist()
    out["prediction_id"] = sv.index.tolist() + sv.index.tolist()

    return out


def get_concordance(sv, changepoints, csv, verbose=True):

    sv = CsvInput(sv).read_csv()

    essential_columns = ["chromosome_1", "chromosome_2", "position_1", "position_2", "type"]

    assert all(col in sv.columns for col in essential_columns)

    sv = sv[["chromosome_1", "chromosome_2", "position_1", "position_2", "type"]]

    matched = match_changepoints(changepoints, sv)

    if csv:
        matched.to_csv(csv, sep="\t", index=False)


    if verbose:
        n_concordant = len(get_concordant_svs(matched))
        n_total = len(matched)
        print("HMM-concordant: {}/{}, {}%".format(n_concordant, n_total, (n_concordant/n_total) * 100))

    return matched


def get_concordant_svs(matched):
    return  matched[(matched.pos_1_in_cn_change_region == True) & (matched.pos_2_in_cn_change_region  == True)]


def get_discordant_svs(matched):
    return matched[(matched.pos_1_in_cn_change_region == False) | (matched.pos_2_in_cn_change_region  == False)]


def split_concordance(matched, csv=False):

    if csv:
        matched = pd.read_csv(matched, sep="\t")

    concordant = get_concordant_svs(matched)
    discordant = get_discordant_svs(matched)

    n_concrdant = len(concordant)
    n_total = len(matched)
    
    return concordant, discordant
    

def check_changepoints_against_matches(concordant, changepoints, verbose=True):
    
    matching_changepoints = pd.DataFrame({"chrom": [], "start":[], "end":[]})
    matching_changepoints["chrom"] = concordant.chr_pos_1 + concordant.chr_pos_2
    matching_changepoints["start"] = concordant.start_pos_1 + concordant.start_pos_2
    matching_changepoints["end"] = concordant.end_pos_1 + concordant.end_pos_2

    matching_changepoints = pd.unique(matching_changepoints.values.ravel())

    if verbose:
        n_matching_changepoints = len(matching_changepoints)
        n_total= len(changepoints)
        print("matching changepoints: {}/{}, {}%".format(n_matching_changepoints, n_total, (n_matching_changepoints/n_total) * 100))


    return match_changepoints

def get_length(r):
    pass

def bin(r):
    if r.chromosome_1 != r.chromosome_2:
        return "interchromosomal"
    else:
        l =  r.end - r.start    
        if l > 100000:
            return "large"
        else:
            return "small"

def bin_breakpoints(sv):
    sv["length"] = sv.apply(lambda r: bin(r), axis=1)


