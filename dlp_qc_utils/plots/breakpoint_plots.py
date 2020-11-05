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