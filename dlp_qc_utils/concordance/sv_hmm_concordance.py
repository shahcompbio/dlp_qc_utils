import numpy as np 
from single_cell.utils.csvutils import IrregularCsvInput, CsvOutput, CsvInput
import statistics
import pandas as pd
import pyranges as pr


def get_cna_changepoints(cna_data, change_region_length = 500000):
    print(len(cna_data))
    cn_data_summary = cna_data.groupby(['chr', 'start', 'end'])
    cn_data_summary = cn_data_summary[['copy', 'state']]
    cn_data_summary = cn_data_summary.aggregate(
        {'copy': statistics.median, 'state': statistics.median}
    )
    print(len(cn_data_summary))
    cn_data_summary["cn_change"]=cn_data_summary.groupby('chr').state.transform(pd.Series.diff)
    cn_data_summary.reset_index(inplace=True)
    cn_data_summary = cn_data_summary[(cn_data_summary.cn_change!=0)].dropna()
    #pandas diff() only includes the higher side of a place of change: i.e.: 1,1,2,2 --> NaN, 0, changepoint, 0
    #include the lower side to get: 1,1,2,2 --> NaN, changepoint, changepoint, 0.
    #accomplish by expanding the lower boundary of every changepoint by the changepoint size. 
    #shouldnt matter if start becomes < 0, as breakpoints can be < 0.
    cn_data_summary["start"] = cn_data_summary.start - change_region_length

    return cn_data_summary

def get_local_breakpoints(changepoint, breakpoints, size_filter=True, stranded=False, csv=True):
    if csv:
        breakpoints = read_parsed_breakpoints(breakpoints)
    if size_filter:
        breakpoints = size_filter_svs(breakpoints)

    chrom = changepoint.chr.tolist()[0]
    start = changepoint.start.tolist()[0]
    end = changepoint.end.tolist()[0]
    matching = breakpoints[(breakpoints.chromosome_1 == chrom) | (breakpoints.chromosome_2 == chrom)]
    # print(matching)
    matching_1 = matching[(matching.position_1 >= start) & (matching.position_1 <= end)]
    matching_2 = matching[(matching.position_2 >= start) & (matching.position_2 <= end)]
    matching = pd.concat([matching_1, matching_2])
    if stranded:
        if changepoints.cn_change > 0:
            return matching[(matching.strand_1 == "-") | (matching.strand_1 == "_")]
        if changepoints.cn_change < 0:
            return matching[(matching.strand_1 == "+") | (matching.strand_1 == "+")]
    else:
        return matching


def anno_bkp(row, cn, side="1"):
    
    if side == "1":
        chrom = str(row.chromosome_1)
        pos = int(row.position_1)
        strand = str(row.strand_1)
    else:
        chrom = str(row.chromosome_2)
        pos = int(row.position_2)     
        strand = str(row.strand_2)

    if strand == "-":
        #look for -
        cn = cn[cn.cn_change > 0]
    else:
        cn = cn[cn.cn_change < 0]

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

# def anno_bkp(row, cn, side="1"):


#     if side == "1":
#         chrom = str(row.chromosome_1)
#         pos = int(row.position_1)
#     else:
#         chrom = str(row.chromosome_2)
#         pos = int(row.position_2)     

#     matches = cn[(cn.start < pos) 
#         & (cn.end > pos) 
#         & (cn.chr == chrom)]
#     matches = matches.reset_index()

#     if matches.empty:
#         row["chr_pos_{}".format(side)] = np.NaN
#         row["start_pos_{}".format(side)] = np.NaN
#         row["end_pos_{}".format(side)] = np.NaN
#         row["state_pos_{}".format(side)] = np.NaN
#         row["cn_change_pos_{}".format(side)] = np.NaN
#     else:
#         row["chr_pos_{}".format(side)] = matches.chr.tolist()[0]
#         row["start_pos_{}".format(side)] = matches.start.tolist()[0]
#         row["end_pos_{}".format(side)] = matches.end.tolist()[0]
#         row["state_pos_{}".format(side)] = matches.state.tolist()[0]
#         row["cn_change_pos_{}".format(side)] = matches.cn_change.tolist()[0]
#     row["pos_{}_in_cn_change_region".format(side)] = not matches.empty
#     return row


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


def filt_row(row, res):
    if row.chromosome_1 != row.chromosome_2:
        return True

    if abs(row.position_1 - row.position_2) >= res:
        return True
    
    return False
    



def size_filter_svs(sv, res=500000):
    return sv[sv.apply(lambda row: filt_row(row, res), axis=1)]

def read_parsed_breakpoints(sv):
    # try:
    #     sv = CsvInput(sv).read_csv()
    # except: 
    sv = pd.read_csv(sv)
    
    sv = sv.astype({"chromosome_1": "str", "chromosome_2":"str", "position_1":"int", "position_2":"int"})
        
    essential_columns = ["chromosome_1", "chromosome_2", "position_1", "position_2", "strand_1", "strand_2", "type"]

    assert all(col in sv.columns for col in essential_columns)

    sv = sv[["chromosome_1", "chromosome_2", "position_1", "position_2", "strand_1", "strand_2", "type"]]
    return sv

def get_concordance(sv, changepoints, csv=True, verbose=True):

    sv = read_parsed_breakpoints(sv)

    sv = size_filter_svs(sv)

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


