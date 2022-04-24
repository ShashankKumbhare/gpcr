
####################################################################################################################################
################################## gpcr_package/plots/scatter_plot_mod  ############################################################
####################################################################################################################################

################################################################
from ..__constants__ import *
from ..__auxil__ import *
################################################################


####################################################################################################################################
################################## scatter_proteins_residue_count_vs_len ###########################################################
# >>
def scatter_proteins_residue_count_vs_len( df                  = GPCR_DF,
                                           residue             = "S",
                                           figsize             = FIGSIZE_SCATTER,
                                           xlim                = None,
                                           ylim                = None,
                                           title_enabled       = True,
                                           title               = None,
                                           title_add           = "",
                                           label_xaxis_enabled = True,
                                           label_xaxis         = "",
                                           label_yaxis_enabled = True,
                                           label_yaxis         = "",
                                           legend_enabled      = True,
                                           legend_location     = "upper right",
                                           scale_factor        = 1,
                                           axes                = None ):
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    # Getting count of residue >>
    count_aa = [ aa_seq.count(residue) for aa_seq in df["seq"] ]
    
    # Plotting scatter plot >>
    fig.axes.scatter( df["seq_len"], count_aa, s = 1.5, color = COLOR_PLOT_SCATTER )
    
    # Setting plot elements >>
    if not title:       title       = "Scatter plot for " + AA.ABBR_LONG_NAMES[residue] + " count vs protein sequence length\n "
    if not label_xaxis: label_xaxis = "Sequence length"
    if not label_yaxis: label_yaxis = "Count"
    set_plot_elements( fig                 = fig,
                       figsize             = figsize,
                       xlim                = xlim,
                       ylim                = ylim,
                       title_enabled       = title_enabled,
                       title               = title,
                       title_add           = title_add,
                       label_xaxis_enabled = label_xaxis_enabled,
                       label_xaxis         = label_xaxis,
                       label_yaxis_enabled = label_yaxis_enabled,
                       label_yaxis         = label_yaxis,
                       legend_enabled      = legend_enabled,
                       legend_location     = legend_location,
                       scale_factor        = scale_factor )
    
    # Show plot >>
    if not axes: plt.show()

    return fig.figure

# <<
################################## scatter_proteins_residue_count_vs_len ###########################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/plots/scatter_plot_mod  ############################################################
####################################################################################################################################
