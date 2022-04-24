
####################################################################################################################################
################################## gpcr_package/plots/line_plot_mod ################################################################
####################################################################################################################################

################################################################
from ..__constants__ import *
from ..__auxil__ import *
from ..__data__ import *
################################################################


####################################################################################################################################
################################## line_compare_points #############################################################################
# >>
def line_compare_points( x                   = None,
                         y                   = None,
                         xerr                = (),
                         yerr                = (),
                         residues            = AA.ABBREVATIONS,
                         linewidth           = 1,
                         errbar_alpha        = 1,
                         errbar_linewidth    = 0.5,
                         errbar_color        = None,
                         errorbar_capsize    = None,
                         errorbar_capthick   = None,
                         scatter_marker      = "o",
                         scatter_marker_size = None,
                         scatter_color       = None,
                         target_enabled      = True,
                         targets_x           = (),
                         targets_x_labels    = (),
                         targets_y           = (),
                         targets_y_labels    = (),
                         figsize             = FIGSIZE_LINE,
                         xlim                = None,
                         ylim                = None,
                         title_enabled       = True,
                         title               = "",
                         title_add           = "",
                         label_xaxis_enabled = True,
                         label_xaxis         = "x",
                         label_yaxis_enabled = True,
                         label_yaxis         = "y",
                         legend_enabled      = True,
                         legend_location     = "upper right",
                         scale_factor        = 1,
                         axes                = None,
                         show                = True ):
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    # Plotting errorbar, points, line
    scatter_colors = [scatter_color] * len(residues) if scatter_color else get_random_color(len(residues))
    errbar_colors  = [errbar_color]  * len(residues) if errbar_color  else scatter_colors
    
    for i, residue in enumerate(residues):
        if not list(xerr):
            xerr_i = None
        else:
            if type(xerr) == tuple:
                xerr_i = ( [xerr[0][i]], [xerr[1][i]] )
            else:
                xerr_i = xerr[i]
        if not list(yerr):
            yerr_i = None
        else:
            if type(yerr) == tuple:
                yerr_i = ( [yerr[0][i]], [yerr[1][i]] )
            else:
                yerr_i = yerr[i]
        
        plt.errorbar( x[i],
                      y[i],
                      xerr      = xerr_i,
                      yerr      = yerr_i,
                      color     = scatter_colors[i],
                      ecolor    = errbar_colors[i],
                      alpha     = errbar_alpha,
                      marker    = scatter_marker,
                      ms        = scatter_marker_size,
                      capsize   = errorbar_capsize,
                      capthick  = errorbar_capthick,
                      linewidth = errbar_linewidth,
                      linestyle = "",
                      label     = AA.ABBR_LONG_NAMES[residue] )
        annotation_residue = plt.annotate(  residue,
                                            (x[i] + 0.001, y[i] + 0.001) )
        annotation_residue.set_fontsize( set_fontsize(figsize, "annotation", scale_factor_extended = scale_factor) )
        
    plt.plot(x, x, color = COLOR_PLOT_LINE, linewidth = linewidth, label = "x vs x")
    
    # Plotting targets lines >>
    if target_enabled:
        plot_x_targets(fig = fig, figsize = figsize, targets_x = targets_x, targets_x_labels = targets_x_labels)
        plot_y_targets(fig = fig, figsize = figsize, targets_y = targets_y, targets_y_labels = targets_y_labels)
    
    # Setting plot elements >>
    if not title: title = label_yaxis + " vs " + label_xaxis
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
    plt.show() if show and not axes else plt.close()
    
    return fig.figure
# <<
################################## line_compare_points #############################################################################
####################################################################################################################################


####################################################################################################################################
################################## line_hydrophobicity_in_sequences ################################################################
# >>
def line_hydrophobicity_in_sequences(   df                  = GPCR_DF,
                                        sequences           = (),
                                        normalized          = False,
                                        hydrophobicity_aas  = AA.HYDRO_PHOBICITY_PH2,
                                        linewidth           = 5,
                                        savefig             = False,
                                        color_min           = "#0000FF", # blue
                                        color_zero          = "#27E900", # green
                                        color_max           = "#FF0000", # red
                                        dataset             = "",
                                        threshold_len_min   = "default",
                                        threshold_len_max   = "default",
                                        figsize             = (30, 32),
                                        xlim                = None,
                                        ylim                = None,
                                        title_enabled       = True,
                                        title               = "",
                                        title_add           = "",
                                        label_xaxis_enabled = True,
                                        label_xaxis         = "",
                                        label_yaxis_enabled = True,
                                        label_yaxis         = "",
                                        legend_enabled      = True,
                                        legend_location     = "upper right",
                                        scale_factor        = 0.5,
                                        axes                = None,
                                        show                = True ):
    
    # Filtering data >>
    if not type(sequences) == pd.core.series.Series:
        if not sequences:
            # if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else "xxxx" ) )
            if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else ( "GPCRDB" if df.equals(GPCRDB_DF) else "xxxx" ) ) )
            if threshold_len_min or threshold_len_max:
                df                  = apply_threshold_len(df, threshold_len_min, threshold_len_max)
            sequences               = list(df["seq"])
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    # Plotting horizontal lines for each C-tail sequence >>
    min_hydrophobicity = min(AA_HYDRO_PHOBICITY_PH2 + AA_HYDRO_PHOBICITY_PH7)
    max_hydrophobicity = max(AA_HYDRO_PHOBICITY_PH2 + AA_HYDRO_PHOBICITY_PH7)
    pbar = ProgressBar()
    j    = -1
    if normalized: locs_aa_for_average = [] # << need this variable only for normalized type plot for its average plot
    for seq in pbar(sequences):
        j = j + 1
        dist_between_two_aas = 1 / len(seq) if normalized else 1
        for aa in AA_ABBREVATIONS: # << Plotting colored line for each type of aa found at their respective locations
            locs_aa = np.array( get_locations_AA( seq, aa, normalized = normalized, indexing = 0 ) )
            if normalized: locs_aa_for_average.append( list(locs_aa[~np.isnan(locs_aa)]) )
            # print(locs_aa_for_average)
            for loc_aa in locs_aa[~np.isnan(locs_aa)]:
                xmin              = loc_aa
                xmax              = loc_aa + dist_between_two_aas
                hydrophobicity_aa = hydrophobicity_aas[aa]
                color_aa          = color_fader( c1  = color_min,
                                                 c2  = color_zero,
                                                 c3  = color_max,
                                                 mix = hydrophobicity_aa / (2*max_hydrophobicity if hydrophobicity_aa > 0 else 2*abs(min_hydrophobicity)) + 0.5 ) # >> need value between 0 and 1
                fig.axes.hlines(y = j, xmin = xmin, xmax = xmax, color = color_aa, linewidth = linewidth)
    if normalized:
        locs_aa_for_average = list( set( sum(locs_aa_for_average, []) ) )
        locs_aa_for_average.sort()
    
    # Inverting y-axis >>
    plt.gca().invert_yaxis()
    
    # Setting plot elements >>
    len_type = "Normalized " if normalized else ""
    if not title:       title       = "Hydrophobicity in sequences" + ", " + dataset
    if not label_xaxis: label_xaxis = len_type + "C-tail length"
    if not label_yaxis: label_yaxis = "Sr. no."
    set_plot_elements(  fig                 = fig,
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
    
    # Saving moment figure >>
    if savefig:
        output_dir     = "output/__line_plots__/"
        new_dir_prefix = "line_" + datetime.today().strftime("%Y%m%d_%H%M%S_%f") + "_hydrophobicity_in_sequences" + "_" + dataset
        new_dir        = new_dir_prefix + ""
        new_dir_path   = output_dir + new_dir
        os.mkdir(new_dir_path)
        filename       = new_dir_prefix + ""
        fig.figure.savefig(new_dir_path + "/" + filename, dpi = 80)
    
    # Show colorbar >>
    cfig = get_figure(figsize = (figsize[0], 1), axes = axes)
    for i in range(min_hydrophobicity, max_hydrophobicity+1):
        _ = cfig.axes.scatter(i, 0, s = 200, color = color_fader( c1  = color_min,
                                                                  c2  = color_zero,
                                                                  c3  = color_max,
                                                                  mix = i / (2*max_hydrophobicity if i > 0 else 2*abs(min_hydrophobicity)) + 0.5 ) )
    # cfig.axes.set_xticks(list(hydrophobicity_aas.values()))
    cfig.axes.set_xlim( (-70, 110) )
    
    # Adding another axes at top for aa labels >>
    axes_top = cfig.axes.twiny()
    axes_top.set_xlim(cfig.axes.get_xlim())
    axes_top.set_xticks(list(hydrophobicity_aas.values()))
    axes_top.set_xticklabels(AA_ABBREVATIONS)
    axes_top.tick_params(axis = 'x', labelsize = set_fontsize(figsize, "axes_ticks", scale_factor_extended = scale_factor * 0.8), rotation = 0)
    set_plot_elements(  fig          = cfig,
                        figsize      = figsize,
                        scale_factor = scale_factor * 0.8 )
    cfig.axes.tick_params(labelleft  = False)
    display(cfig.figure)
    
    # Show plot >>
    plt.show(fig.figure.number) if show and not axes else plt.close()
    
    # Showing average >>
    avgFig  = get_figure(figsize = (figsize[0], 1), axes = axes)
    len_max = max( [ len(seq) for seq in sequences ] )
    if not normalized:
        for pos in range(len_max):
            aas_pos                    = [seq[pos] if len(seq) > pos else "" for seq in sequences ]
            hydrophobicity_aas_pos     = [ hydrophobicity_aas[aa] if aa else np.nan for aa in aas_pos ]
            hydrophobicity_aas_pos_avg = np.nanmean(np.array(hydrophobicity_aas_pos))
            xmin                  = pos
            dist_between_two_aas  = 1 / len(seq) if normalized else 1
            xmax                  = pos + dist_between_two_aas
            color_aa              = color_fader( c1  = color_min,
                                                 c2  = color_zero,
                                                 c3  = color_max,
                                                 mix = hydrophobicity_aas_pos_avg / (2*max_hydrophobicity if hydrophobicity_aas_pos_avg > 0 else 2*abs(min_hydrophobicity)) + 0.5 ) # >> need value between 0 and 1
            avgFig.axes.hlines(y = j, xmin= xmin, xmax = xmax, color = color_aa, linewidth = linewidth)
    else:
        ############################################################################################################################
        # OLD >>
        # len_max                = max( [ len(seq) for seq in sequences ] )
        # pos_len_max_normalized = [ pos / len_max for pos in range(len_max) ]
        # pos_normalized         = [[]] * len(sequences)
        # for i, seq in enumerate(sequences):
        #     pos_normalized[i] = [ pos / len(seq) for pos in range(len(seq)) ]
        # aas_pos    = [ [""]*len(sequences) ] * len_max
        # aas_pos[0] = [ seq[0] for seq in sequences ]
        # pos_normalized_check = [ pos_norm[1] for pos_norm in pos_normalized ]
        # counter_check        = [ 1 for _ in sequences ]
        # out    = [0] * len_max
        # out[0] = aas_pos[0]
        # for pos in range(1, len_max):
        #     for i, seq in enumerate(sequences):
        #         if pos_normalized_check[i] <= pos_len_max_normalized[pos]:
        #             aas_pos[pos][i]         = seq[counter_check[i]]
        #             pos_normalized_check[i] = pos_normalized[i][counter_check[i] + 1] if counter_check[i] + 1 < len(seq) else 1
        #             counter_check[i]        = counter_check[i] + 1
        #         else:
        #             aas_pos[pos][i]         = aas_pos[pos-1][i]
        #         out[pos] = tuple(aas_pos[pos])
        # out = [ list(item) for item in out ]
        # weights_sequences    = np.array([ len(seq) / len_max for seq in sequences ])
        # sum_weight           = sum(weights_sequences)
        # aas_pos = out
        # for pos in range(len_max):
        #     hydrophobicity_aas_pos     = [ hydrophobicity_aas[aa] for aa in aas_pos[pos] ]
        #     hydrophobicity_aas_pos_avg = sum( weights_sequences * hydrophobicity_aas_pos ) / sum_weight
        #     xmin                  = pos / len_max
        #     dist_between_two_aas  = 1 / len_max
        #     xmax                  = xmin + dist_between_two_aas
        #     color_aa              = color_fader( c1  = color_min,
        #                                          c2  = color_zero,
        #                                          c3  = color_max,
        #                                          mix = hydrophobicity_aas_pos_avg / (2*max_hydrophobicity if hydrophobicity_aas_pos_avg > 0 else 2*abs(min_hydrophobicity)) + 0.5 ) # >> need value between 0 and 1
        #     avgFig.axes.hlines(y = j, xmin= xmin, xmax = xmax, color = color_aa, linewidth = linewidth)
        ############################################################################################################################
        # NEW >>
        # Getting new sequences >>
        locs_aas_sequences = [ list( np.array( range( len(seq) ) ) / len(seq) ) for seq in sequences ]
        positions_all      = list( set( sum(locs_aas_sequences, []) ) )
        positions_all.sort()
        # print(positions_all)
        sequences_new      = [[]] * len(sequences)
        for k, seq in enumerate(sequences):
            sequences_new[k] = [""] * ( len(positions_all))
            j          = 0
            aa_current = seq[j]
            aa_next    = seq[j+1]
            loc_end    = (len(seq) - 1) / len(seq)
            for i, pos in enumerate(positions_all):
                if pos > loc_end or pos < locs_aas_sequences[k][j+1]:
                    sequences_new[k][i] = aa_current
                else:
                    sequences_new[k][i] = aa_next
                    j          = j + 1
                    aa_current = seq[j]
                    if j < len(seq)-1: aa_next = seq[j+1]
            sequences_new[k] = ''.join(sequences_new[k])
        
        # Calculating average >>
        len_max                    = max( [ len(seq) for seq in sequences ] )
        weights_sequences          = np.array([ len(seq) / len_max for seq in sequences ])
        sum_weight                 = sum(weights_sequences)
        hydrophobicity_aas_pos_avg = [0] * len(positions_all)
        for i, pos in enumerate(positions_all):
            hydrophobicity_aas_pos        = [ hydrophobicity_aas[aa] for aa in [ aa[i] for aa in sequences_new ] ]
            scale_factor_hydrophobicity   = 1.7
            hydrophobicity_aas_pos_avg[i] = sum( weights_sequences * hydrophobicity_aas_pos ) / sum_weight * scale_factor_hydrophobicity
        
        # Plotting average >>
        for i, pos in enumerate(positions_all):
            # hydrophobicity_aas_pos     = [ hydrophobicity_aas[aa] for aa in [ aa[i] for aa in sequences_new ] ]
            # hydrophobicity_aas_pos_avg = sum( weights_sequences * hydrophobicity_aas_pos ) / sum_weight
            xmin = pos
            xmax = positions_all[i+1] if i < (len(positions_all) - 1) else 1
            # print(xmin, xmax)
            color_aa             = color_fader( c1  = color_min,
                                                c2  = color_zero,
                                                c3  = color_max,
                                                mix = hydrophobicity_aas_pos_avg[i] / (2*max_hydrophobicity if hydrophobicity_aas_pos_avg[i] > 0 else 2*abs(min_hydrophobicity)) + 0.5 ) # >> need value between 0 and 1
            avgFig.axes.hlines(y = j, xmin= xmin, xmax = xmax, color = color_aa, linewidth = linewidth)
        
        ############################################################################################################################
    avgFig.axes.tick_params(axis = 'x', labelsize = set_fontsize(figsize, "axes_ticks", scale_factor_extended = scale_factor * 0.8), rotation = 0)
    set_plot_elements(  fig          = avgFig,
                        title        = "Avegare Hydrophobicity",
                        figsize      = figsize,
                        scale_factor = scale_factor * 0.8 )
    avgFig.axes.tick_params(labelleft  = False)
    
    return fig.figure, hydrophobicity_aas_pos_avg
# <<
################################## line_hydrophobicity_in_sequences ################################################################
####################################################################################################################################


####################################################################################################################################
################################## line_charges_in_sequences #######################################################################
# >>
def line_charges_in_sequences(  df                  = GPCR_DF,
                                sequences           = (),
                                normalized          = False,
                                charges_aas         = AA.CHARGES_PH7,
                                linewidth           = 5,
                                savefig             = False,
                                color_min           = "#0000FF", # blue
                                color_zero          = "#27E900", # green
                                color_max           = "#FF0000", # red
                                dataset             = "",
                                threshold_len_min   = "default",
                                threshold_len_max   = "default",
                                figsize             = (30, 32),
                                xlim                = None,
                                ylim                = None,
                                title_enabled       = True,
                                title               = "",
                                title_add           = "",
                                label_xaxis_enabled = True,
                                label_xaxis         = "",
                                label_yaxis_enabled = True,
                                label_yaxis         = "",
                                legend_enabled      = True,
                                legend_location     = "upper right",
                                scale_factor        = 0.5,
                                axes                = None,
                                show                = True ):
    
    # Filtering data >>
    if not type(sequences) == pd.core.series.Series:
        if not sequences:
            # if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else "xxxx" ) )
            if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else ( "GPCRDB" if df.equals(GPCRDB_DF) else "xxxx" ) ) )
            if threshold_len_min or threshold_len_max:
                df                  = apply_threshold_len(df, threshold_len_min, threshold_len_max)
            sequences               = list(df["seq"])
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    # Plotting horizontal lines for each C-tail sequence >>
    min_charge = -1
    max_charge = 1
    pbar = ProgressBar()
    j    = -1
    for seq in pbar(sequences):
        j = j + 1
        dist_between_two_aas = 1 / len(seq) if normalized else 1
        for aa in AA_ABBREVATIONS: # << Plotting colored line for each type of aa found at their respective locations
            locs_aa = np.array( get_locations_AA( seq, aa, normalized = normalized, indexing = 0 ) )
            for loc_aa in locs_aa[~np.isnan(locs_aa)]:
                xmin              = loc_aa
                xmax              = loc_aa + dist_between_two_aas
                charges_aa        = charges_aas[aa]
                color_aa          = color_fader( c1  = color_min,
                                                 c2  = color_zero,
                                                 c3  = color_max,
                                                 mix = charges_aa / (2*max_charge if charges_aa > 0 else 2*abs(min_charge)) + 0.5 ) # >> need value between 0 and 1
                fig.axes.hlines(y = j, xmin = xmin, xmax = xmax, color = color_aa, linewidth = linewidth)
    
    # Inverting y-axis >>
    plt.gca().invert_yaxis()
    
    # Setting plot elements >>
    len_type = "Normalized " if normalized else ""
    if not title:       title       = "Charges in sequences" + ", " + dataset
    if not label_xaxis: label_xaxis = len_type + "C-tail length"
    if not label_yaxis: label_yaxis = "Sr. no."
    set_plot_elements(  fig                 = fig,
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
    
    # Saving moment figure >>
    if savefig:
        output_dir     = "output/__line_plots__/"
        new_dir_prefix = "line_" + datetime.today().strftime("%Y%m%d_%H%M%S_%f") + "_charges_in_sequences" + "_" + dataset
        new_dir        = new_dir_prefix + ""
        new_dir_path   = output_dir + new_dir
        os.mkdir(new_dir_path)
        filename       = new_dir_prefix + ""
        fig.figure.savefig(new_dir_path + "/" + filename, dpi = 80)
    
    # Show colorbar >>
    cfig = get_figure(figsize = (figsize[0], 1), axes = axes)
    for i in [-1, 0, 1]:
        cfig.axes.scatter(i, 0, s = 200, color = color_fader( c1  = color_min,
                                                              c2  = color_zero,
                                                              c3  = color_max,
                                                              mix = i / (2*max_charge if i > 0 else 2*abs(min_charge)) + 0.5 ) )
    cfig.axes.set_xticks( [-1, 0, 1] )
    
    # Adding another axes at top for aa labels >>
    axes_top = cfig.axes.twiny()
    axes_top.set_xlim(cfig.axes.get_xlim())
    axes_top.set_xticks( [-1, 0, 1] )
    tick_lable_negative = str( [aa for aa in AA_ABBREVATIONS if charges_aas[aa] == -1] )
    tick_lable_zero     = str( [aa for aa in AA_ABBREVATIONS if charges_aas[aa] == 0] )
    tick_lable_positive = str( [aa for aa in AA_ABBREVATIONS if charges_aas[aa] == 1] )
    tick_lables_charges = [ tick_lable_negative, tick_lable_zero, tick_lable_positive ]
    axes_top.set_xticklabels( tick_lables_charges )
    axes_top.tick_params(axis = 'x', labelsize = set_fontsize(figsize, "axes_ticks", scale_factor_extended = scale_factor), rotation = 0)
    set_plot_elements(  fig          = cfig,
                        figsize      = figsize,
                        scale_factor = scale_factor )
    cfig.axes.tick_params(labelleft  = False)
    display(cfig.figure)
    
    # Show plot >>
    plt.show(fig.figure.number) if show and not axes else plt.close()
    
    return fig.figure
# <<
################################## line_charges_in_sequences #######################################################################
####################################################################################################################################


####################################################################################################################################
################################## help_code_weighted_avg_hydrophobicity ###########################################################
# sequences = ["ABCDABCDE",
#              "IJKLMNOP",
#              "UVWX" ]
#
# df        = gpcr.GPCRDB_DF
# df        = gpcr.apply_threshold_len(df, 20, None)
# sequences = list(df["seq"])
# sequences
#
# len_max = max( [ len(seq) for seq in sequences ] )
#
# pos_len_max_normalized = [ (pos) / len_max for pos in range(len_max) ]
# pos_len_max_normalized
#
# pos_normalized = [[]] * len(sequences)
# for i, seq in enumerate(sequences):
#     pos_normalized[i] = [ (pos) / len(seq) for pos in range(len(seq)) ]
#
# aas_pos    = [ [""]*len(sequences) ] * len_max
# aas_pos
# aas_pos[0] = [ seq[0] for seq in sequences ]
# print(f"pos = {0}")
# print(aas_pos[0])
# pos_normalized_check = [ pos_norm[1] for pos_norm in pos_normalized ]
# # print(f"{pos_normalized_check}")
# counter_check = [ 1 for _ in sequences ]
# counter_check
# for pos in range(1, len_max):
#     print(f"pos = {pos}")
#     for i, seq in enumerate(sequences):
#         # print(f"i = {i}")
#         if pos_normalized_check[i] <= pos_len_max_normalized[pos]:
#             aas_pos[pos][i]         = seq[counter_check[i]]
#             # print(aas_pos[pos])
#             # counter_check[i]        = counter_check[i] + 1
#             # print( f"i = {i}, counter_check[i] = {counter_check[i]}" )
#             # print( len(pos_normalized[i]) )
#             pos_normalized_check[i] = pos_normalized[i][counter_check[i] + 1] if counter_check[i] + 1 < len(seq) else 1
#             counter_check[i]        = counter_check[i] + 1
#             # print(f"Change: {pos_normalized_check}")
#         else:
#             aas_pos[pos][i]         = aas_pos[pos-1][i]
#             # print(f"No change{pos_normalized_check}")
#             # print(aas_pos[pos])
#     print(aas_pos[pos])
################################## help_code_avg_hydrophobicity ####################################################################
####################################################################################################################################


####################################################################################################################################
################################## line_residue_in_ctail ###########################################################################
# >>
def line_residue_in_ctail(  df                  = GPCRDB_DF,
                            sequences           = (),
                            residue             = 'S',
                            normalized          = False,
                            # aas_acidic          = AA_ACIDIC,
                            # aas_hydrophobic     = AA_HYDROPHIBOC,
                            linewidth           = None,
                            # savefig             = False,
                            dataset             = "",
                            threshold_len_min   = "default",
                            threshold_len_max   = "default",
                            figsize             = (30, 20),
                            xlim                = None,
                            ylim                = None,
                            title_enabled       = True,
                            title               = "",
                            title_add           = "",
                            label_xaxis_enabled = True,
                            label_xaxis         = "",
                            label_yaxis_enabled = True,
                            label_yaxis         = "",
                            legend_enabled      = True,
                            legend_location     = "upper right",
                            scale_factor        = 0.5,
                            axes                = None,
                            show                = True ):

    # Filtering data >>
    if not type(sequences) == pd.core.series.Series:
        if not sequences:
            # if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else "xxxx" ) )
            if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else ( "GPCRDB" if df.equals(GPCRDB_DF) else "xxxx" ) ) )
            if threshold_len_min or threshold_len_max:
                df                  = apply_threshold_len(df, threshold_len_min, threshold_len_max)
            sequences               = list(df["seq"])
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    # Plotting horizontal line for residue for each C-tail sequence >>
    color_residue    = "black"
    color_no_residue = "dodgerblue" # "lightseagreen"
    pbar = ProgressBar()
    j    = -1
    for seq in pbar(sequences):
        j = j + 1
        dist_between_two_aas = 1 / len(seq) if normalized else 1
        locs_aa              = get_locations_AA( seq, residue, normalized = normalized, indexing = 0 )
        if locs_aa == [np.nan]:    # << Plotting `color_no_residue` lines if residue is not found
            fig.axes.hlines( y         = j,
                             xmin      = 0,
                             xmax      = 1 if normalized else len(seq),
                             color     = color_no_residue,
                             alpha     = 1,
                             linewidth = linewidth)
        else:
            for loc_aa in locs_aa: # << Plotting "color_residue" line for each residue found at their respective locations
                xmin              = loc_aa
                xmax              = loc_aa + dist_between_two_aas
                # alpha = 0.5 if include_mean else 1
                fig.axes.hlines( y         = j,
                                 xmin      = xmin,
                                 xmax      = xmax,
                                 color     = color_residue,
                                 linewidth = linewidth) # , alpha = alpha)
        end_of_seq = len(seq) if not normalized else 1
        fig.axes.scatter( x = end_of_seq, y = j, s = 1, color = "black", marker = "|" )

    # Setting plot elements >>
    len_type = "Normalized " if normalized else ""
    if not title:       title       = "'" + residue + "'" + " residue in " + dataset + " C-tail sequences"
    if not label_xaxis: label_xaxis = len_type + "C-tail length"
    if not label_yaxis: label_yaxis = "Sr. no."
    if not xlim: xlim = (0, None)
    if not ylim: ylim = (0, len(sequences))
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
    
    # Inverting y-axis >>
    plt.gca().invert_yaxis()
    
    # Show plot >>
    plt.show() if show and not axes else plt.close()
    
    return fig.figure

# <<
################################## line_residue_in_ctail ###########################################################################
####################################################################################################################################


####################################################################################################################################
################################## line_motifs_in_ctail ############################################################################
# >>
def line_motifs_in_ctail(   df                  = GPCRDB_DF,
                            sequences           = (),
                            motif               = 'pdz',
                            motif_class         = (1, 2, 3),
                            normalized          = True,
                            aas_acidic          = AA_ACIDIC,
                            aas_hydrophobic     = AA_HYDROPHIBOC,
                            linewidth           = None,
                            dataset             = "",
                            threshold_len_min   = "default",
                            threshold_len_max   = "default",
                            figsize             = (30, 20),
                            xlim                = None,
                            ylim                = None,
                            title_enabled       = True,
                            title               = "",
                            title_add           = "",
                            label_xaxis_enabled = True,
                            label_xaxis         = "",
                            label_yaxis_enabled = True,
                            label_yaxis         = "",
                            legend_enabled      = True,
                            legend_location     = "upper right",
                            scale_factor        = 0.5,
                            axes                = None,
                            show                = True ):

    # Filtering data >>
    if not type(sequences) == pd.core.series.Series:
        if not sequences:
            # if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else "xxxx" ) )
            if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else ( "GPCRDB" if df.equals(GPCRDB_DF) else "xxxx" ) ) )
            if threshold_len_min or threshold_len_max:
                df                  = apply_threshold_len(df, threshold_len_min, threshold_len_max)
            sequences               = list(df["seq"])
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    # Plotting horizontal line for residue for each C-tail sequence >>
    color_class_I_motif   = "red"
    color_class_II_motif  = "green"
    color_class_III_motif = "blue"
    color_no_residue      = "grey" # "lightseagreen"
    pbar = ProgressBar()
    j    = -1
    for seq in pbar(sequences):
        j = j + 1
        dist_between_two_aas = 1 / len(seq) if normalized else 1
        locs_motif_all = []
        # Plotting location of Class I `pdz` motifs >>
        if 1 in motif_class:
            locs_1st = get_locations_AA( seq, ["S", "T"], normalized = False, indexing = 0 )
            locs_3rd = []
            for loc in locs_1st:
                if (loc < len(seq)-2) and (seq[int(loc)+2] in aas_hydrophobic):
                    locs_3rd.append(loc+2)
            locs_motif_I   = np.array(locs_3rd) - 1
            locs_motif_I   = locs_motif_I / (len(seq) if normalized else 1)
            # locs_motif_all = locs_motif_all.append(locs_motif_I)
            locs_motif_all = np.concatenate([locs_motif_all, locs_motif_I])
            for loc_aa in locs_motif_I:
                xmin              = loc_aa
                xmax              = loc_aa + dist_between_two_aas
                fig.axes.hlines( y         = j,
                                 xmin      = xmin,
                                 xmax      = xmax,
                                 color     = color_class_I_motif,
                                 linewidth = linewidth)
        # Plotting location of Class II `pdz` motifs >>
        if 2 in motif_class:
            locs_1st = get_locations_AA( seq, aas_hydrophobic, normalized = False, indexing = 0 )
            locs_3rd = []
            for loc in locs_1st:
                if (loc < len(seq)-2) and (seq[int(loc)+2] in aas_hydrophobic):
                    locs_3rd.append(loc+2)
            locs_motif_II  = np.array(locs_3rd) - 1
            locs_motif_II  = locs_motif_II / (len(seq) if normalized else 1)
            # locs_motif_all = locs_motif_all.append(locs_motif_II)
            locs_motif_all = np.concatenate([locs_motif_all, locs_motif_II])
            for loc_aa in locs_motif_II:
                xmin              = loc_aa
                xmax              = loc_aa + dist_between_two_aas
                fig.axes.hlines( y         = j,
                                 xmin      = xmin,
                                 xmax      = xmax,
                                 color     = color_class_II_motif,
                                 linewidth = linewidth)
        # Plotting location of Class III `pdz` motifs >>
        if 3 in motif_class:
            locs_1st = get_locations_AA( seq, aas_acidic, normalized = False, indexing = 0 )
            locs_3rd = []
            for loc in locs_1st:
                if (loc < len(seq)-2) and (seq[int(loc)+2] in aas_acidic):
                    locs_3rd.append(loc+2)
            locs_motif_III = np.array(locs_3rd) - 1
            locs_motif_III = locs_motif_III / (len(seq) if normalized else 1)
            # locs_motif_all = locs_motif_all.append(locs_motif_III)
            locs_motif_all = np.concatenate([locs_motif_all, locs_motif_III])
            for loc_aa in locs_motif_III:
                xmin              = loc_aa
                xmax              = loc_aa + dist_between_two_aas
                fig.axes.hlines( y         = j,
                                 xmin      = xmin,
                                 xmax      = xmax,
                                 color     = color_class_III_motif,
                                 linewidth = linewidth)
        # Plotting `color_no_residue` lines if "S" / "T" is not found >>
        if locs_motif_all.size == 0:
            fig.axes.hlines( y         = j,
                             xmin      = 0,
                             xmax      = 1 if normalized else len(seq),
                             color     = color_no_residue,
                             alpha     = 0.25,
                             linewidth = linewidth)
        # Terminating sequences with a marker >>
        end_of_seq = len(seq) if not normalized else 1
        fig.axes.scatter( x = end_of_seq, y = j, s = 1, color = "black", marker = "|" )

    # Setting plot elements >>
    len_type = "Normalized " if normalized else ""
    if not title:       title       = "'" + motif + "'" + " binding motif in " + dataset + " C-tail sequences"
    if not label_xaxis: label_xaxis = len_type + "C-tail length"
    if not label_yaxis: label_yaxis = "Sr. no."
    if not xlim: xlim = (0, None)
    if not ylim: ylim = (0, len(sequences))
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
    
    # Inverting y-axis >>
    plt.gca().invert_yaxis()
    
    # Show plot >>
    plt.show() if show and not axes else plt.close()
    
    return fig.figure

# <<
################################## line_motifs_in_ctail ############################################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/plots/line_plot_mod ################################################################
####################################################################################################################################
