from dlp_qc_utils.concordance import sv_hmm_concordance
from dlp_qc_utils.plots import p


def test_check_changepoints_against_matches():
    changepoints =pd.read_csv("/juno/work/shah/abramsd/CODE/SV_bench/changepoints.csv")
    matched = pd.read_csv("/juno/work/shah/abramsd/CODE/SV_bench/notebooks/cons_primary", sep="\t")
    conc = get_concordant_svs(matched)
    check_changepoints_against_matches(conc, changepoints)


def test_get_cna_changepoints():