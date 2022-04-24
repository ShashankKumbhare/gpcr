
####################################################################################################################################
################################## gpcr_package/plots/box_plot_mod #################################################################
####################################################################################################################################

################################################################
from ..__constants__ import *
from ..__auxil__ import *
from ..__data__ import *
################################################################


####################################################################################################################################
################################## box_proteins_residue_loc ########################################################################
# >>
def box_proteins_residue_loc( df                 = GPCR_DF,
                              residues           = AA_ABBREVATIONS,
                              normalized         = True,
                              mean               = False,
                              # bins               = None,
                              # enable_target      = True,
                              # xlim               = False,
                              title_enabled      = True,
                              title              = None,
                              figsize            = None,
                              enable_xaxis_label = True,
                              enable_yaxis_label = True,
                              # enable_legend      = True,
                              axes               = None ):
    
    # Setting plot size >>
    if not figsize: figsize = FIGSIZE_BOX[0]*len(residues), FIGSIZE_BOX[1]
    fig = get_figure(figsize = figsize, axes = axes)
    
    # getting locations of residue locations >>
    box_data = [np.nan] * len(residues)
    for i, residue in enumerate(residues):
        locs_aa     = [ get_locations_AA(aa_seq, residue, normalized) for aa_seq in df["seq"] ]
        if mean:
            mean_locs_aa = np.array( [ stat.mean(loc_aa) for loc_aa in locs_aa ] )
            box_data[i]  = mean_locs_aa[~np.isnan(mean_locs_aa)]
        else:
            locs_aa_all  = np.array( list(flatten_list(locs_aa)) )
            box_data[i]  = locs_aa_all[~np.isnan(locs_aa_all)]
        # fig.axes.boxplot( box_data[i], widths = 0.5 )
        fig.axes.boxplot( box_data[i], positions = [i], widths = 0.5 )
    
    # Calculating mean of residue locations >>
    # mean_locs_residue = dict()
    # for amino in aminos:
    #     mean_locs_residue[amino] = get_mean(df, amino, normalized)
    #     mean_locs_residue[amino] = np.array(mean_locs_residue[amino])
    #     mean_locs_residue[amino] = mean_locs_residue[amino][~np.isnan(mean_locs_residue[amino])]
    
    # Plotting box plot >>
    # fig.axes.boxplot( box_data, positions = range(len(mean_loc_aminos)), widths = 0.5 )
    # fig.axes.boxplot( box_data, widths = 0.5 )
    
    # plt.plot([], c='#D7191C', label='Apples')
    # plt.plot([], c='#2C7BB6', label='Oranges')
    
    # fig.axes.boxplot( box_data, positions = range(len(mean_loc_aminos)), widths = 0.5 )
    # fig.axes.hist( box_data, bins = bins, alpha = 1, density = True, color = "lightseagreen")
    
    # Setting plot title, axes-lables, legend, axes-ticks >>
    len_type  = "Normalized" if normalized else "Absolute"
    mean_type = " mean " if mean else " individual "
    location  = "locations for each sequence" if mean else "locations for all the sequences"
    if title_enabled:
        # if not title: title = AA.ABBR_LONG_NAMES[residues] + mean_type + location + " \n "
        if not title: title = "Amono acid" + mean_type + location + " \n "
        fig.axes.set_title(title, fontsize = set_fontsize(figsize, "title"))
    if enable_xaxis_label:
        xlabel = "Amino acids"
        fig.axes.set_xlabel(xlabel, fontsize = set_fontsize(figsize, "axes_label"))
    if enable_yaxis_label:
        ylabel = len_type + " location in GPCR C-Tail"
        fig.axes.set_ylabel(ylabel, fontsize = set_fontsize(figsize, "axes_label"))
    # fig.axes.tick_params(axis = 'x', labelsize = set_fontsize(figsize, "xtick"))
    # fig.axes.tick_params(axis = 'y', labelsize = set_fontsize(figsize, "ytick"))
    
    if not axes: plt.show()

    return fig.figure

# <<
################################## box_proteins_residue_loc ########################################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/plots/box_plot_mod #################################################################
####################################################################################################################################
