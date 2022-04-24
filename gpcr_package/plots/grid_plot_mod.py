
####################################################################################################################################
################################## gpcr_package/plots/grid #########################################################################
####################################################################################################################################

################################################################
from ..__constants__ import *
from ..__auxil__ import *
from ..__data__ import *
from .hist_plot_mod import hist_proteins_residue_loc, hist_proteins_residue_relative_count
from .scatter_plot_mod import scatter_proteins_residue_count_vs_len
################################################################


####################################################################################################################################
################################## grid_hist_proteins_residue_loc ##################################################################
# >>
def grid_hist_proteins_residue_loc( df                    = GPCR_DF,
                                    residues              = AA_ABBREVATIONS,
                                    normalized            = True,
                                    take_mean             = True,
                                    take_log              = True,
                                    threshold_len_min     = "default",
                                    threshold_len_max     = "default",
                                    centerline_enabled    = True,
                                    enable_target         = True,
                                    bins                  = None,
                                    xlim                  = None,
                                    ylim                  = None,
                                    title_enabled         = True,
                                    title                 = "",
                                    title_add             = "",
                                    subfigsize            = SUBFIGSIZE_GRID,
                                    sharex                = "all",
                                    sharey                = "all",
                                    subplot_title_enabled = True,
                                    subplot_titles        = (),
                                    subplot_title_add     = (),
                                    label_xaxis_enabled   = True,
                                    label_xaxis           = "",
                                    label_yaxis_enabled   = True,
                                    label_yaxis           = "",
                                    legend_enabled        = True,
                                    legend_location       = "upper right",
                                    scale_factor          = 1,
                                    show                  = True):
    
    # Setting plot size >>
    ncols                = ceil( np.sqrt(len(residues)) )
    nrows                = len(residues) // ncols + bool(len(residues) % ncols)
    fig                  = Struct()
    figsize              = (ncols*subfigsize[0], nrows*subfigsize[1])
    fig.figure, fig.axes = plt.subplots( nrows, ncols, sharex = sharex, sharey = sharey, figsize = figsize )
    
    # Plotting subplots >>
    if not subplot_titles:    subplot_titles    = [ AA.ABBR_LONG_NAMES[residue] for residue in residues ]
    if not subplot_title_add: subplot_title_add = [""] * len(subplot_titles)
    k = -1
    for i in range(nrows):
        for j in range(ncols):
            k = k + 1
            if k < len(residues):
                flag_xaxis_label = True if label_xaxis_enabled and i == nrows - 1 else False
                flag_yaxis_label = True if label_yaxis_enabled and j == 0         else False
                if len(residues) == 1:
                    current_axis = fig.axes
                else:
                    current_axis = fig.axes[j] if nrows == 1 else fig.axes[i,j]
                hist_proteins_residue_loc( df                  = df,
                                           residue             = residues[k],
                                           normalized          = normalized,
                                           take_mean           = take_mean,
                                           take_log            = take_log,
                                           threshold_len_min   = threshold_len_min,
                                           threshold_len_max   = threshold_len_max,
                                           centerline_enabled  = centerline_enabled,
                                           enable_target       = enable_target,
                                           bins                = bins,
                                           xlim                = xlim,
                                           ylim                = ylim,
                                           title_enabled       = subplot_title_enabled,
                                           title               = subplot_titles[k],
                                           title_add           = subplot_title_add[k],
                                           label_xaxis_enabled = flag_xaxis_label,
                                           label_xaxis         = label_xaxis,
                                           label_yaxis_enabled = flag_yaxis_label,
                                           label_yaxis         = label_yaxis,
                                           legend_enabled      = legend_enabled,
                                           legend_location     = legend_location,
                                           scale_factor        = scale_factor,
                                           axes                = current_axis )
                
                if sharex: current_axis.xaxis.set_tick_params(labelbottom = True)
                if sharey: current_axis.yaxis.set_tick_params(labelleft   = True)
    
    # Setting plot title, axis lables >>
    if title_enabled:
        mean_type = " mean "                      if take_mean else " individual "
        location  = "locations for each sequence" if take_mean else "locations for all the sequences"
        log       = " (log)"                      if take_log  else ""
        # dataset   = ": GPCR" if df.equals(GPCR_DF) else ( ": DisProt" if df.equals(DISPROT_DF) else ( ": D2P2" if df.equals(D2P2_DF) else "" ) )
        dataset   = ": GPCR" if df.equals(GPCR_DF) else ( ": DisProt" if df.equals(DISPROT_DF) else ( ": D2P2" if df.equals(D2P2_DF) else ( ": GPCRDB" if df.equals(GPCRDB_DF) else ": xxxx" ) ) )
        title     = (title if title else "Histogram for AAs" + log + mean_type + location) + dataset + title_add
        fig.figure.suptitle(title, fontsize = set_fontsize(figsize, "gridtitle", scale_factor))
    # if enable_xaxis_label:
    #     print("")
    #     fig.figure.set_xlabel(len_type + " location in GPCR C-Tail", fontsize = set_fontsize(figsize, "axes_label"))
    # if enable_yaxis_label:
    #     ylabel = "Density" if normalized else "Count"
    #     print("")
    #     fig.figure.set_ylabel(ylabel, fontsize = set_fontsize(figsize, "axes_label"))
    
    plt.show() if show else plt.close()
    
    return fig.figure

# <<
################################## grid_hist_proteins_residue_loc ##################################################################
####################################################################################################################################


####################################################################################################################################
################################## grid_scatter_proteins_residue_count_vs_len ######################################################
# >>
def grid_scatter_proteins_residue_count_vs_len( df                    = GPCR_DF,
                                                residues              = AA_ABBREVATIONS,
                                                xlim                  = None,
                                                ylim                  = None,
                                                title_enabled         = True,
                                                title                 = "",
                                                subplot_title_enabled = True,
                                                subplot_titles        = None,
                                                sharex                = "all",
                                                sharey                = "all",
                                                label_xaxis_enabled   = True,
                                                label_xaxis           = "",
                                                label_yaxis_enabled   = True,
                                                label_yaxis           = "",
                                                legend_enabled        = True,
                                                legend_location       = "upper right",
                                                subfigsize            = SUBFIGSIZE_GRID,
                                                show                  = True):
    
    # Setting plot size >>
    ncols                = ceil( np.sqrt(len(residues)) )
    nrows                = len(residues) // ncols + bool(len(residues) % ncols)
    fig                  = Struct()
    figsize              = (ncols*subfigsize[0], nrows*subfigsize[1])
    fig.figure, fig.axes = plt.subplots( nrows, ncols, sharex = sharex, sharey = sharey, figsize = figsize )
    
    # Plotting subplots >>
    if not subplot_titles:
        subplot_titles = [ AA.ABBR_LONG_NAMES[residue] for residue in residues ]
    k = -1
    for i in range(nrows):
        for j in range(ncols):
            k = k + 1
            if k < len(residues):
                flag_xaxis_label = True if label_xaxis_enabled and i == nrows - 1 else False
                flag_yaxis_label = True if label_yaxis_enabled and j == 0         else False
                if len(residues) == 1:
                    current_axis = fig.axes
                else:
                    current_axis = fig.axes[j] if nrows == 1 else fig.axes[i,j]
                    
                scatter_proteins_residue_count_vs_len( df                  = df,
                                                       residue             = residues[k],
                                                       xlim                = xlim,
                                                       ylim                = ylim,
                                                       title_enabled       = subplot_title_enabled,
                                                       title               = subplot_titles[k],
                                                       label_xaxis_enabled = flag_xaxis_label,
                                                       label_xaxis         = label_xaxis,
                                                       label_yaxis_enabled = flag_yaxis_label,
                                                       label_yaxis         = label_yaxis,
                                                       legend_enabled      = legend_enabled,
                                                       legend_location     = legend_location,
                                                       axes                = current_axis )
                
                if sharex: current_axis.xaxis.set_tick_params(labelbottom = True)
                if sharey: current_axis.yaxis.set_tick_params(labelleft   = True)
    
    # Setting plot title, axis lables >>
    if title_enabled:
        if not title: title = "Scatter plots for amino acid" + " count vs protein sequence length \n"
        fig.figure.suptitle(title, fontsize = set_fontsize(figsize, "gridtitle"))
    
    plt.show() if show else plt.close()
    
    return fig.figure

# <<
################################## grid_scatter_proteins_residue_count_vs_len ######################################################
####################################################################################################################################


####################################################################################################################################
################################## grid_hist_proteins_residue_relative_count #######################################################
# >>
def grid_hist_proteins_residue_relative_count(  df                    = GPCR_DF,
                                                residues              = AA_ABBREVATIONS,
                                                relative_count        = True,
                                                take_log              = False,
                                                threshold_len_min     = "default",
                                                threshold_len_max     = "default",
                                                target_type           = "mean",
                                                target_enabled        = True,
                                                targets_x           = (),
                                                targets_x_colors    = (),
                                                targets_x_labels    = (),
                                                targets_y           = (),
                                                targets_y_colors    = (),
                                                targets_y_labels    = (),
                                                confidence_interval   = (0.1, 0.9),
                                                bins                  = "auto",
                                                xlim                  = None,
                                                ylim                  = None,
                                                title_enabled         = True,
                                                title                 = "",
                                                title_add             = "",
                                                subfigsize            = SUBFIGSIZE_GRID,
                                                sharex                = "all",
                                                sharey                = "all",
                                                subplot_title_enabled = True,
                                                subplot_titles        = (),
                                                subplot_title_add     = (),
                                                label_xaxis_enabled   = True,
                                                label_xaxis           = "",
                                                label_yaxis_enabled   = True,
                                                label_yaxis           = "",
                                                legend_enabled        = True,
                                                legend_location       = "upper right",
                                                scale_factor          = 1,
                                                show                  = True ):
    
    # Setting plot size >>
    ncols                = ceil( np.sqrt(len(residues)) )
    nrows                = len(residues) // ncols + bool(len(residues) % ncols)
    fig                  = Struct()
    figsize              = (ncols*subfigsize[0], nrows*subfigsize[1])
    fig.figure, fig.axes = plt.subplots( nrows, ncols, sharex = sharex, sharey = sharey, figsize = figsize )
    
    # Plotting subplots >>
    target_return     = np.array( [np.nan] * len(residues) )
    left_conf_return  = np.array( [np.nan] * len(residues) )
    right_conf_return = np.array( [np.nan] * len(residues) )
    if not subplot_titles:    subplot_titles    = [ AA.ABBR_LONG_NAMES[residue] for residue in residues ]
    if not subplot_title_add: subplot_title_add = [""]     * len(residues)
    if not targets_x:         targets_x         = [np.nan] * len(residues)
    # if not targets_x_colors:  targets_x_colors  = [np.nan] * len(residues)
    if not targets_x_colors:  targets_x_colors  = get_random_color( len(residues) )
    if not targets_x_labels:  targets_x_labels  = [""]     * len(residues)
    if not targets_y:         targets_y         = [np.nan] * len(residues)
    # if not targets_y_colors:  targets_y_colors  = [np.nan] * len(residues)
    if not targets_y_colors:  targets_y_colors  = get_random_color( len(residues) )
    if not targets_y_labels:  targets_y_labels  = [""]     * len(residues)
    k = -1
    for i in range(nrows):
        for j in range(ncols):
            k = k + 1
            if k < len(residues):
                flag_xaxis_label = True if label_xaxis_enabled and i == nrows - 1 else False
                flag_yaxis_label = True if label_yaxis_enabled and j == 0         else False
                if len(residues) == 1:
                    current_axis = fig.axes
                else:
                    current_axis = fig.axes[j] if nrows == 1 else fig.axes[i,j]
                _, target_return[k], left_conf_return[k], right_conf_return[k], _ = hist_proteins_residue_relative_count(
                                                                                df                  = df,
                                                                                residue             = residues[k],
                                                                                relative_count      = relative_count,
                                                                                take_log            = take_log,
                                                                                threshold_len_min   = threshold_len_min,
                                                                                threshold_len_max   = threshold_len_max,
                                                                                target_type         = target_type,
                                                                                target_enabled      = target_enabled,
                                                                                targets_x           = targets_x[k],
                                                                                targets_x_colors    = targets_x_colors[k],
                                                                                targets_x_labels    = targets_x_labels[k],
                                                                                targets_y           = targets_y[k],
                                                                                targets_y_colors    = targets_y_colors[k],
                                                                                targets_y_labels    = targets_y_labels[k],
                                                                                confidence_interval = confidence_interval,
                                                                                bins                = bins,
                                                                                xlim                = xlim,
                                                                                ylim                = ylim,
                                                                                title_enabled       = subplot_title_enabled,
                                                                                title               = subplot_titles[k],
                                                                                title_add           = subplot_title_add[k],
                                                                                label_xaxis_enabled = flag_xaxis_label,
                                                                                label_xaxis         = label_xaxis,
                                                                                label_yaxis_enabled = flag_yaxis_label,
                                                                                label_yaxis         = label_yaxis,
                                                                                legend_enabled      = legend_enabled,
                                                                                legend_location     = legend_location,
                                                                                scale_factor        = scale_factor,
                                                                                axes                = current_axis )
                if sharex: current_axis.xaxis.set_tick_params(labelbottom = True)
                if sharey: current_axis.yaxis.set_tick_params(labelleft   = True)
    
    # Setting plot title, axis lables >>
    if title_enabled:
        log      = " (log)"    if take_log           else ""
        relative = " relative" if relative_count     else ""
        # dataset  = ": GPCR"    if df.equals(GPCR_DF) else ( ": DisProt" if df.equals(DISPROT_DF) else ( ": D2P2" if df.equals(D2P2_DF) else "" ) )
        dataset  = ": GPCR"    if df.equals(GPCR_DF) else ( ": DisProt" if df.equals(DISPROT_DF) else ( ": D2P2" if df.equals(D2P2_DF) else ( ": GPCRDB" if df.equals(GPCRDB_DF) else ": xxxx" ) ) )
        title    = (title if title else "Histogram for AAs" + log + relative + " count") + dataset + title_add
        # title    = (title if title else "Histogram for AAs" + log + relative + " count") + df.name + title_add
        fig.figure.suptitle(title, fontsize = set_fontsize(figsize, "gridtitle"))
    
    # Show plot >>
    plt.show() if show else plt.close()
    
    return fig.figure, target_return, left_conf_return, right_conf_return

# <<
################################## grid_hist_proteins_residue_relative_count #######################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/plots/grid #########################################################################
####################################################################################################################################
