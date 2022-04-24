
####################################################################################################################################
################################## gpcr_package/plots/heatmap_plot_mod #############################################################
####################################################################################################################################

################################################################
from ..__constants__ import *
from ..__auxil__ import *
from ..__data__ import *
################################################################


####################################################################################################################################
################################## heatmap_1 #######################################################################################
# >>
def heatmap_1( df                  = GPCR_DF,
               sequences           = (),
               lag                 = 1,
               moments             = 1,       # "all": "all implies [1, 1.5, 2, 3, 4]"
               join_sequences      = False,
               get_figlist         = False,
               savefig             = False,
               show_sum            = False,
               clim                = (None, None),
               dataset             = "",
               threshold_len_min   = "default",
               threshold_len_max   = "default",
               figsize             = FIGSIZE_HEATMAP,
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
               scale_factor        = 0.85,
               cmap                = None,  # COLOR_COLORBAR_CMAP,
               axes                = None,
               show                = True ):
    
    # Filtering data >>
    if not type(sequences) == pd.core.series.Series:
        if not sequences:
            # if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else "xxxx" ) )
            if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else ( "GPCRDB" if df.equals(GPCRDB_DF) else "xxxx" ) ) )
            df                      = apply_threshold_len(df, threshold_len_min, threshold_len_max)
            sequences               = list(df["seq"])
    
    # Computing relative frequencies >>
    aas                     = AA_ABBREVATIONS
    aa_pair_freqs_dicts     = np.array([dict()] * len(sequences))
    aa_pair_freqs_rel_dicts = np.array([dict()] * len(sequences))
    sum_of_count            = np.array([dict()] * len(sequences))
    if join_sequences:
        for i, seq in enumerate(sequences):
            aa_pair_freqs_dicts[i]     = { aa_start : { aa_next : count_overlap(seq, aa_start, aa_next, lag) + count_overlap(seq, aa_next, aa_start, lag)                        for aa_next in aas } for aa_start in aas }
    else:
        for i, seq in enumerate(sequences):
            aa_pair_freqs_dicts[i]     = { aa_start : { aa_next : count_overlap(seq, aa_start, aa_next, lag) + count_overlap(seq, aa_next, aa_start, lag)                        for aa_next in aas } for aa_start in aas }
            sum_of_count[i]            = { aa_start : sum( aa_pair_freqs_dicts[i][aa_start].values() )                                                                                                for aa_start in aas }
            aa_pair_freqs_rel_dicts[i] = { aa_start : { aa_next : aa_pair_freqs_dicts[i][aa_start][aa_next] / sum_of_count[i][aa_start] if sum_of_count[i][aa_start] != 0 else np.nan for aa_next in aas } for aa_start in aas }
    
    # Calculating statistics for all sequences >>
    if join_sequences:
        sum_lens_sequences             = sum([ len(seq) for seq in sequences ] )
        no_of_sequences                = len(sequences)
        moment_aa_pair_freqs_rel_dicts = { aa_start : { aa_next : sum( [ dct[aa_start][aa_next]         for dct in aa_pair_freqs_dicts ] ) / (sum_lens_sequences - no_of_sequences * lag) for aa_next in aas } for aa_start in aas }
        moments                        = 1
    else:
        moments                        = [1, 1.5, 2, 3, 4] if moments == "all" or get_figlist else (moments if type(moments) == list else [moments])
        moment_aa_pair_freqs_rel_dicts = [{}] * len(moments)
        moment_str                     = [""] * len(moments)
        # clims                          = [np.nan] * len(moments)
        for i, moment in enumerate(moments):
            if   moment == 1:   take_moment = np.nanmean;                                 moment_str[i] = "mean" # ;     clims[i] = (0, 0.18) if not clim else clim
            elif moment == 1.5: take_moment = np.nanstd;                                  moment_str[i] = "stddev" # ;   clims[i] = (0, 0.25)
            elif moment == 2:   take_moment = np.nanvar;                                  moment_str[i] = "var" # ;      clims[i] = (0, 0.06)
            elif moment == 3:   take_moment = lambda x: skew(x, nan_policy = "omit");     moment_str[i] = "skew" # ;     clims[i] = (0, 16)
            elif moment == 4:   take_moment = lambda x: kurtosis(x, nan_policy = "omit"); moment_str[i] = "kurtosis" # ; clims[i] = (0, 290)
            else:               print("Please enter correct moment. Allowed values are: 1, 1.5, 2, 3, 4"); return None, None, None, None
            moment_aa_pair_freqs_rel_dicts[i] = { aa_start : { aa_next : take_moment( [ dct[aa_start][aa_next] for dct in aa_pair_freqs_rel_dicts ] ) for aa_next in aas } for aa_start in aas }
    
    # Making dataframe >>
    moment_aa_pair_freqs_rel_dfs = [ dict_to_df( dct = dct ) for dct in moment_aa_pair_freqs_rel_dicts ]
    if show_sum:
        for df in moment_aa_pair_freqs_rel_dfs: df["Sum"] = df.sum(axis = 1)
    
    # Plotting heatmap >>
    figs         = [ get_figure(figsize = figsize, axes = axes) for _ in range(len(moment_aa_pair_freqs_rel_dfs)) ]
    title_prefix = (title if title else "Co-occurrence Prob(AA-pair: lag ") + str(lag)+ "): "
    for i, (fig, df, moment) in enumerate(zip(figs, moment_aa_pair_freqs_rel_dfs, moments)):
        heatmap_data = df
        fig.axes.matshow( heatmap_data, cmap = cmap, vmin = clim[0] if moment == 1 else None, vmax = clim[1] if moment == 1 else None )
        fig.axes.set_xticks(np.arange(0, len(heatmap_data.columns),                     1))
        fig.axes.set_yticks(np.arange(0, len(heatmap_data.columns[:len(heatmap_data)]), 1))
        fig.axes.set_xticklabels(heatmap_data.columns,                     rotation = None)
        fig.axes.set_yticklabels(heatmap_data.columns[:len(heatmap_data)], rotation = None)
        fig.axes.tick_params(axis = 'x', labelsize = None, labelbottom = True, labeltop   = True)
        fig.axes.tick_params(axis = 'y', labelsize = None, labelleft   = True, labelright = True)
        # Adding colorbar >>
        get_colorbar_heatmap(fig = fig)
        # Setting plot elements >>
        title        = title_prefix + moment_str[i] + ", " + dataset
        newline      = "\n" if title or title_add else ""
        set_plot_elements( fig                 = fig,
                           figsize             = figsize,
                           xlim                = xlim,
                           ylim                = ylim,
                           title_enabled       = title_enabled,
                           title               = title,
                           title_add           = title_add + newline,
                           label_xaxis_enabled = label_xaxis_enabled,
                           label_xaxis         = label_xaxis,
                           label_yaxis_enabled = label_yaxis_enabled,
                           label_yaxis         = label_yaxis,
                           legend_enabled      = legend_enabled,
                           legend_location     = legend_location,
                           scale_factor        = scale_factor )
    
    # Saving moment figure >>
    if savefig:
        output_dir     = "output/__heatmap_plots__/"
        new_dir_prefix = "heatmap_" + datetime.today().strftime("%Y%m%d_%H%M%S_%f") + "_co_occur_prob" + "_" + dataset + "_lag_" + str(lag)
        no_of_figures  = "_(" + str( len(moments) + (len(sequences) if get_figlist else 0) ) + "_figures)"
        new_dir        = new_dir_prefix + no_of_figures
        new_dir_path   = output_dir + new_dir
        os.mkdir(new_dir_path)
        for i, (fig, moment) in enumerate(zip(figs, moments)):
            filename = new_dir_prefix + "_"*(len(moments)-i) + moment_str[i]
            fig.figure.savefig(new_dir_path + "/" + filename, dpi = 80)
    
    # Adding plots for individual sequences >>
    if not join_sequences:
        if get_figlist:
            fig_list              = np.array([Struct()] * len(sequences))
            aa_pair_freqs_rel_dfs = np.array([dict()] * len(sequences))
            pbar = ProgressBar()
            i    = -1
            for aa_pair_freqs_rel_dict in pbar(aa_pair_freqs_rel_dicts):
                i = i + 1
                fig_list[i]              = get_figure(figsize = figsize, axes = axes)
                # Making dataframe >>
                aa_pair_freqs_rel_dfs[i] = dict_to_df( dct = aa_pair_freqs_rel_dict )
                if show_sum:  aa_pair_freqs_rel_dfs[i]["Sum"] = aa_pair_freqs_rel_dfs[i].sum(axis = 1)
                # Plotting heatmap >>
                heatmap_data             = aa_pair_freqs_rel_dfs[i]
                fig_list[i].axes.matshow( heatmap_data, cmap = cmap, vmin = 0, vmax = 1 )
                fig_list[i].axes.set_xticks(np.arange(0, len(heatmap_data.columns),                     1))
                fig_list[i].axes.set_yticks(np.arange(0, len(heatmap_data.columns[:len(heatmap_data)]), 1))
                fig_list[i].axes.set_xticklabels(heatmap_data.columns,                     rotation = None)
                fig_list[i].axes.set_yticklabels(heatmap_data.columns[:len(heatmap_data)], rotation = None)
                fig_list[i].axes.tick_params(axis = 'x', labelsize = None, labelbottom = True, labeltop   = True)
                fig_list[i].axes.tick_params(axis = 'y', labelsize = None, labelleft   = True, labelright = True)
                set_plot_elements( fig                 = fig_list[i],
                                   figsize             = figsize,
                                   xlim                = xlim,
                                   ylim                = ylim,
                                   title_enabled       = title_enabled,
                                   title               = title_prefix + " " + dataset + f", seq no. {str(i)}",
                                   title_add           = title_add + newline,
                                   label_xaxis_enabled = label_xaxis_enabled,
                                   label_xaxis         = label_xaxis,
                                   label_yaxis_enabled = label_yaxis_enabled,
                                   label_yaxis         = label_yaxis,
                                   legend_enabled      = legend_enabled,
                                   legend_location     = legend_location,
                                   scale_factor        = scale_factor )
                # Adding colorbar >>
                get_colorbar_heatmap(fig = fig_list[i])
                # Saving individual figures >>
                if savefig:
                    filename = new_dir_prefix + "_seq_" + str(i)
                    fig_list[i].figure.savefig(new_dir_path + "/" + filename, dpi = 80)
                plt.close(fig_list[i].figure.number)
        else:
            fig_list              = []
            aa_pair_freqs_rel_dfs = []
    else:
        fig_list              = []
        aa_pair_freqs_rel_dfs = []
    
    # Show plot >>
    plt.figure(fig.figure.number)
    plt.show() if show and not axes else plt.close()
    
    return [fig.figure for fig in figs], [fig.figure for fig in fig_list], moment_aa_pair_freqs_rel_dfs, aa_pair_freqs_rel_dfs
# <<
################################## heatmap_1 #######################################################################################
####################################################################################################################################


####################################################################################################################################
################################## heatmap_2 #######################################################################################
# >>
def heatmap_2( df                  = GPCR_DF,
               sequences           = (),
               lag                 = 1,
               moment              = 1,
               join_sequences      = False,
               threshold_len_min   = "default",
               threshold_len_max   = "default",
               get_figlist         = False,
               savefig             = False,
               figsize             = FIGSIZE_HEATMAP,
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
               scale_factor        = 1,
               cmap                = None,  # COLOR_COLORBAR_CMAP,
               axes                = None,
               show                = True ):
    
    # Filtering data >>
    dataset  = "GPCR"    if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else "" ) )
    if not type(sequences) == pd.core.series.Series:
        if not sequences:
            df = apply_threshold_len(df, threshold_len_min, threshold_len_max)
            # sequences = df["seq"]
            sequences = list(df["seq"])
            # sequences = np.array(df["seq"])
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    # Computing relative frequencies >>
    aas                     = AA_ABBREVATIONS
    aa_pair_freqs_dicts     = np.array([dict()] * len(sequences))
    aa_pair_freqs_rel_dicts = np.array([dict()] * len(sequences))
    if join_sequences:
        for i, seq in enumerate(sequences):
            aa_pair_freqs_dicts[i]     = { aa_start : { aa_next : count_overlap(seq, aa_start, aa_next, lag)                    for aa_next in aas } for aa_start in aas }
    else:
        for i, seq in enumerate(sequences):
            aa_pair_freqs_rel_dicts[i] = { aa_start : { aa_next : count_overlap(seq, aa_start, aa_next, lag) / (len(seq) - lag) for aa_next in aas } for aa_start in aas }
    
    # Calculating statistics for all sequences >>
    if join_sequences:
        sum_lens_sequences            = sum([ len(seq) for seq in sequences ] )
        no_of_sequences               = len(sequences)
        moment_aa_pair_freqs_rel_dict = { aa_start : { aa_next : sum( [ dct[aa_start][aa_next]         for dct in aa_pair_freqs_dicts ] ) / (sum_lens_sequences - no_of_sequences * lag) for aa_next in aas } for aa_start in aas }
        moment                        = 1
    else:
        if   moment == 1:   take_moment = stat.mean
        elif moment == 1.5: take_moment = np.std
        elif moment == 2:   take_moment = np.var
        elif moment == 3:   take_moment = skew
        elif moment == 4:   take_moment = kurtosis
        else:               print("Please enter correct moment. Allowed values are: 1, 1.5, 2, 3, 4"); return None, None, None, None
        moment_aa_pair_freqs_rel_dict = { aa_start : { aa_next : take_moment( [ dct[aa_start][aa_next] for dct in aa_pair_freqs_rel_dicts ] ) for aa_next in aas } for aa_start in aas }
    
    # Making dataframe >>
    def get_dataframe(dct):
        dataframe = pd.DataFrame()
        for key_1, dict_inside in dct.items():
            dataframe[key_1] = dict_inside.values()
        dataframe.index = dct.keys()
        return dataframe
    moment_aa_pair_freqs_rel_df = get_dataframe( dct = moment_aa_pair_freqs_rel_dict )
    
    # Plotting heatmap >>
    heatmap_data = moment_aa_pair_freqs_rel_df
    fig.axes.matshow( heatmap_data, cmap = cmap )
    fig.axes.set_xticks(np.arange(0, len(heatmap_data.columns), 1))
    fig.axes.set_yticks(np.arange(0, len(heatmap_data.columns), 1))
    fig.axes.set_xticklabels(heatmap_data.columns, rotation = None)
    fig.axes.set_yticklabels(heatmap_data.columns, rotation = None)
    fig.axes.tick_params(axis = 'x', labelsize = None, labelbottom = True, labeltop   = True)
    fig.axes.tick_params(axis = 'y', labelsize = None, labelleft   = True, labelright = True)
    
    # Adding colorbar >>
    get_colorbar_heatmap(fig = fig)
    
    # Setting plot elements >>
    title_prefix = (title if title else "Co-occurrence Prob(AA-pair: lag ") + str(lag)+ "): "
    title        = title_prefix + "moment " + str(moment) + ", " + dataset
    newline      = "\n" if title or title_add else ""
    set_plot_elements( fig                 = fig,
                       figsize             = figsize,
                       xlim                = xlim,
                       ylim                = ylim,
                       title_enabled       = title_enabled,
                       title               = title,
                       title_add           = title_add + newline,
                       label_xaxis_enabled = label_xaxis_enabled,
                       label_xaxis         = label_xaxis,
                       label_yaxis_enabled = label_yaxis_enabled,
                       label_yaxis         = label_yaxis,
                       legend_enabled      = legend_enabled,
                       legend_location     = legend_location,
                       scale_factor        = scale_factor )
    
    # Saving figure >>
    output_dir      = "output/heatmap/"
    new_dir_prefix  = "heatmap_" + datetime.today().strftime('%Y%m%d_%H%M%S_%f') + "_lag_" + str(lag)
    new_dir         = new_dir_prefix + "_moment_" + str(moment) + "_" + dataset
    dirtype         = "figlist" if get_figlist else ""
    new_dir_path    = output_dir + new_dir + ("_" + dirtype if dirtype else "")
    filename_prefix = new_dir_prefix
    if savefig:
        os.mkdir(new_dir_path)
        filename = filename_prefix + "_moment_" + str(moment) + "_" + dataset
        fig.figure.savefig(new_dir_path + "/" + filename, dpi = 80)
    
    # Adding plots for individual sequences >>
    if not join_sequences:
        if get_figlist:
            fig_list              = np.array([Struct()] * len(sequences))
            aa_pair_freqs_rel_dfs = np.array([dict()] * len(sequences))
            pbar = ProgressBar()
            i    = -1
            for aa_pair_freqs_rel_dict in pbar(aa_pair_freqs_rel_dicts):
                i = i + 1
                fig_list[i]              = get_figure(figsize = figsize, axes = axes)
                # Making dataframe >>
                aa_pair_freqs_rel_dfs[i] = get_dataframe( dct = aa_pair_freqs_rel_dict )
                # Plotting heatmap >>
                heatmap_data             = aa_pair_freqs_rel_dfs[i]
                fig_list[i].axes.matshow( heatmap_data, cmap = cmap )
                fig_list[i].axes.set_xticks(np.arange(0, len(heatmap_data.columns), 1))
                fig_list[i].axes.set_yticks(np.arange(0, len(heatmap_data.columns), 1))
                fig_list[i].axes.set_xticklabels(heatmap_data.columns, rotation = None)
                fig_list[i].axes.set_yticklabels(heatmap_data.columns, rotation = None)
                fig_list[i].axes.tick_params(axis = 'x', labelsize = None, labelbottom = True, labeltop   = True)
                fig_list[i].axes.tick_params(axis = 'y', labelsize = None, labelleft   = True, labelright = True)
                set_plot_elements( fig                 = fig_list[i],
                                   figsize             = figsize,
                                   xlim                = xlim,
                                   ylim                = ylim,
                                   title_enabled       = title_enabled,
                                   title               = title_prefix + " " + dataset + f", seq no. {str(i)}",
                                   title_add           = title_add + newline,
                                   label_xaxis_enabled = label_xaxis_enabled,
                                   label_xaxis         = label_xaxis,
                                   label_yaxis_enabled = label_yaxis_enabled,
                                   label_yaxis         = label_yaxis,
                                   legend_enabled      = legend_enabled,
                                   legend_location     = legend_location,
                                   scale_factor        = scale_factor )
                # Adding colorbar >>
                get_colorbar_heatmap(fig = fig_list[i])
                # Saving figure >>
                if savefig:
                    filename = filename_prefix +  "_" + dataset + "_seq_" + str(i)
                    fig_list[i].figure.savefig(new_dir_path + "/" + filename, dpi = 80)
                plt.close(fig_list[i].figure.number)
        else:
            fig_list              = []
            aa_pair_freqs_rel_dfs = []
    else:
        fig_list              = []
        aa_pair_freqs_rel_dfs = []
    
    # Show plot >>
    plt.figure(fig.figure.number)
    plt.show() if show and not axes else plt.close()
    
    return fig.figure, [fig.figure for fig in fig_list], moment_aa_pair_freqs_rel_df, aa_pair_freqs_rel_dfs
# <<
################################## heatmap_2 #######################################################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/plots/heatmap_plot_mod #############################################################
####################################################################################################################################
