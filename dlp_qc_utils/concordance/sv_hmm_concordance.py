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


def flatten_matched_breakpoints(matched, csv=False):
    if csv:
        matched = read_matched_data(matched)

    flat = pd.DataFrame({"chromosome":[], "position":[], "type":[], "prediction_id":[],
        "hmm_chrom":[], "hmm_start":[], "hmm_end":[]}
    )

    flat["chromosome"] = matched.chromosome_1.tolist() + matched.chromosome_2.tolist()
    flat["position"] = matched.position_1.tolist() + matched.position_2.tolist()
    flat["rearrangement_type"] = matched.type.tolist() + matched.type.tolist()
    flat["prediction_id"] = matched.index.tolist() + matched.index.tolist()
    flat["hmm_chrom"] = matched.chr_pos_1.tolist() + matched.chr_pos_2.tolist()
    flat["hmm_start"] = matched.start_pos_1.tolist() + matched.start_pos_2.tolist()
    flat["hmm_end"] = matched.end_pos_1.tolist() + matched.end_pos_2.tolist()

    return flat


def get_concordance(sv, changepoints, csv=True, verbose=True):

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


def get_concordant_svs(matched, csv=False):
    if csv:
        matched = read_matched_data(matched)
    return  matched[(matched.pos_1_in_cn_change_region == True) & (matched.pos_2_in_cn_change_region  == True)]


def get_discordant_svs(matched, csv=False):
    if csv:
        matched = read_matched_data(matched)
    return matched[(matched.pos_1_in_cn_change_region == False) | (matched.pos_2_in_cn_change_region  == False)]


def split_concordance(matched, csv=False):
    if csv:
        matched = read_matched_data(matched)

    concordant = get_concordant_svs(matched)
    discordant = get_discordant_svs(matched)

    n_concrdant = len(concordant)
    n_total = len(matched)
    
    return concordant, discordant
    
def read_matched_data(matched):
    matched = pd.read_csv(matched, sep="\t")
    return matched.astype({"chromosome_1":"str", "chromosome_2": "str", "chr_pos_1": "str", "chr_pos_2": "str"})


def check_changepoints_against_matches(matched, changepoints, csv=False, verbose=True):
    if csv:
        matched = read_matched_data(matched)

    matches = get_concordant_svs(matched)

    matches =  flatten_matched_breakpoints(matches)[["hmm_chrom", "hmm_start", "hmm_end"]]
    # print(matches)

    matches = matches.apply(lambda r: r.tolist(), axis=1)
    matches = [list(x) for x in set(tuple(x) for x in matches.tolist())]

    if verbose:
        n_matching_changepoints = len(matches)
        n_total= len(changepoints)
        print(len(matched))
        # print(n_matching_changepoints/n_total)


    return match_changepoints

def get_length(r):
    pass

def bin(r):
    if r.chromosome_1 != r.chromosome_2:
        return "interchromosomal"
    else:
        l =  r.position_2 - r.position_1   
        if l <= 1000:
            return "<= 1000" 
        if l > 1000 and l <= 10000:
            return "1kb-10kb" 
        if l > 10000 and l <= 100000:
            return "10kb - 100kb" 
        if l > 100000 and l <= 1000000:
            return "100kb-1mb" 
        if l > 1000000 and l <= 10000000:
            return "1mb-10mb" 
        if l > 10000000:
            return "> 10mb" 
        else:
            print(l)
            return "?"


def bin_breakpoints(sv):
    sv["length"] = sv.apply(lambda r: bin(r), axis=1)


