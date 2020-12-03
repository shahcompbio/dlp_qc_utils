import wgs_analysis.plots.rearrangement
from  scgenome.cnplot import plot_cell_cn_profile
import seaborn

def plot(ax, sv_hmm_concordance, changepoints, title=""):

    if len(sv_hmm_concordance) > 0:

        wgs_analysis.plots.rearrangement.chromosome_type_plot(
            ax, sv_hmm_concordance,rearrangement_types= sv_hmm_concordance.rearrangement_type.unique(),
            full_genome=True, bin_size=1000000)
        wgs_analysis.plots.rearrangement.chromosome_type_plot_legend(ax, rearrangement_types=sv_hmm_concordance.rearrangement_type.unique())
    
    ax2 = ax.twinx()

    plot_cell_cn_profile(
            ax2, changepoints, 'copy', 'state'
        )

    ax2.yaxis.tick_right()
    seaborn.despine(ax=ax2, offset=0, trim=True)

    left_color = 'tab:blue'
    ax.tick_params(axis='y', labelcolor=left_color)

    
    right_color = 'tab:red'
    ax2.set_ylabel('copy number') 
    ax2.tick_params(axis='y', labelcolor=right_color)
    ax.set_title(title)

    def bin(r):
#     print(r.position_2, r.position_1, abs(r.position_2 - r.position_1) )
    if r.chromosome_1 != r.chromosome_2:
        return "interchromosomal"
    else:
        l =  abs(r.position_2 - r.position_1) 
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
        
def type_plot(matched, ax1, ax2, fig, title="", csv=True):
    
    fig.suptitle(title, fontsize=14)
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)
    
    if csv:
        matched = read_matched_data(matched)
    print(matched)
    matches = get_concordant_svs(matched)
    non_matches = get_discordant_svs(matched)
#     print(matches)
    if len(matches) > 0:
        matches["bin"] = matches.apply(lambda row: bin(row), axis=1)
        matches=matches.groupby(["bin", "type"]).size().reset_index(name="count")
        matches.pivot("bin", "type", "count").plot(kind='bar', ax=ax1)

    if len(non_matches) > 0:
        non_matches["bin"] = non_matches.apply(lambda row: bin(row), axis=1)
        non_matches=non_matches.groupby(["bin", "type"]).size().reset_index(name="count")
        non_matches.pivot("bin", "type", "count").plot(kind='bar', ax=ax2)

    ax1.set_ylabel("count")
    ax1.set_title("concordant")
    ax2.set_title("discordant")