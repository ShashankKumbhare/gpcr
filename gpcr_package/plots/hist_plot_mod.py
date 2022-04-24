
####################################################################################################################################
################################## gpcr_package/plots/hist_plot_mod ################################################################
####################################################################################################################################

################################################################
from ..__constants__ import *
from ..__auxil__ import *
################################################################


####################################################################################################################################
################################## hist_quick ######################################################################################
# >>
def hist_quick( hist_data           = "",
                take_log            = False,
                target_enabled      = True,
                targets_x           = (),
                targets_x_labels    = (),
                bins                = None,
                figsize             = FIGSIZE_HIST_DEFAULT,
                xlim                = None,
                ylim                = None,
                title_enabled       = True,
                title               = "",
                title_add           = "",
                label_xaxis_enabled = True,
                label_xaxis         = "Data",
                label_yaxis_enabled = True,
                label_yaxis         = "Count",
                legend_enabled      = True,
                legend_location     = "upper right",
                scale_factor        = 1,
                axes                = None,
                show                = True ):
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    # Plotting histograms >>
    if take_log: hist_data = np.log(hist_data)
    hist_data_filtered = [x for x in hist_data if ~np.isnan(x)]
    hist_data_nan      = [x for x in hist_data if np.isnan(x)]
    bins               = set_bins(hist_data, bins, xlim)
    fig.axes.hist( hist_data_filtered, bins = bins, color = COLOR_PLOT_HIST_DEFAULT )
    
    # Printing stats for all C-tail lengths >>
    count_total    = len(hist_data);                  print(f"Total count: {count_total}")
    count_nan      = len(hist_data_nan);              print(f"Nan count:   {count_nan}")
    count          = len(hist_data_filtered);         print(f"Hist Count:  {count}")
    minimum        = min(hist_data_filtered);         print(f"Min:         {minimum}")
    mean           = stat.mean(hist_data_filtered);   print(f"Mean:        {mean:.4f}")
    median         = stat.median(hist_data_filtered); print(f"Median:      {median:.4f}")
    maximum        = max(hist_data_filtered);         print(f"Max:         {maximum}")
    skew_values    = skew(hist_data_filtered);        print(f"Skew:        {skew_values}")
    kurtosis_value = kurtosis(hist_data_filtered);    print(f"kurtosis:    {kurtosis_value}")
    
    # Plotting targets lines >>
    if target_enabled:
        fig.axes.axvline( mean,   color = COLOR_MEAN_LINE,   label = f"Mean  : {mean:.4f}",   linewidth = set_linewidth(figsize) )
        fig.axes.axvline( median, color = COLOR_MEDIAN_LINE, label = f"Median: {median:.4f}", linewidth = set_linewidth(figsize) )
        plot_x_targets(fig = fig, figsize = figsize, targets_x = targets_x, targets_x_labels = targets_x_labels)
    
    # Setting plot elements >>
    if not title: title = label_yaxis + " vs " + label_xaxis + title_add
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
################################## hist_quick ######################################################################################
####################################################################################################################################


####################################################################################################################################
################################## hist_gpcr_Ctail_lens ############################################################################
# >>
def hist_gpcr_Ctail_lens( df                   = GPCR_DF,
                          gpcr_classes         = GPCR_CLASSES,
                          # take_log             = True,
                          # enable_target        = True,
                          bins                 = None,
                          figsize              = (14, 8),
                          xlim                 = None,
                          ylim                 = None,
                          title_enabled        = True,
                          title                = "",
                          title_add            = "",
                          label_xaxis_enabled  = True,
                          label_xaxis          = "",
                          label_yaxis_enabled  = True,
                          label_yaxis          = "",
                          legend_enabled       = True,
                          legend_location      = "upper right",
                          scale_factor         = 1,
                          axes                 = None ):
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    linewidth_vline = 1.5
    
    # Plotting histograms for all C-tail lengths >>
    all_C_tail_lens = df["seq_len"]
    bins = set_bins(all_C_tail_lens, bins, xlim)
    y, x, _         = fig.axes.hist( all_C_tail_lens, bins = bins, alpha = 0.2, color = "grey", label = "All classes combined" )
    
    # Printing stats for all C-tail lengths >>
    mean_C_tail_len_all = np.mean(all_C_tail_lens)
    print("All GPCRs: Count = " + str(len(all_C_tail_lens)) + ", Mean C-tail length = " + str(mean_C_tail_len_all))
    plt.axvline( mean_C_tail_len_all, color = "black", label = "Mean length for all", linewidth = linewidth_vline )
    
    # Plotting histograms for C-tail lengths of individual classes >>
    C_tail_lens = df[ df["class"].isin(gpcr_classes)]["seq_len" ]
    colors = set_colors_for_classes(gpcr_classes)
    for i, gpcr_class in enumerate(gpcr_classes):
        C_tail_len = df[df["class"] == gpcr_class]["seq_len"]
        fig.axes.hist( C_tail_len,
                       bins  = bins,
                       color = colors[i],
                       label = "Class: " + gpcr_class )
        # Printing stats for C-tail lengths of individual classes >>
        mean_C_tail_len_i = np.mean(C_tail_len)
        print("Class " + gpcr_class + " GPCRs: Count = " + str(len(C_tail_len)) + ", Mean C-tail length = " + str(mean_C_tail_len_i))
        plt.axvline( mean_C_tail_len_i, color = colors[i], label = "Mean length for " + gpcr_class, linewidth = linewidth_vline )
    
    # Printing stats for requested C-tail lengths >>
    mean_C_tail_len = np.mean(C_tail_lens)
    print("Class " + ", ".join(gpcr_classes) + " GPCRs: Count = " + str(len(C_tail_lens)) + ", Mean C-tail length = " + str(mean_C_tail_len))
    plt.axvline( mean_C_tail_len, color = "red", label = "Mean length for " + ",".join(gpcr_classes), linewidth = 2 )
    
    # Setting plot elements >>
    if not title:       title       = "Histogram of GPCR C-tail length"
    if not label_xaxis: label_xaxis = "C-tail length"
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
################################## hist_gpcr_Ctail_lens ############################################################################
####################################################################################################################################


####################################################################################################################################
################################## hist_proteins_residue_loc #######################################################################
# >>
def hist_proteins_residue_loc( df                  = GPCR_DF,
                               residue             = "S",
                               normalized          = True,
                               take_mean           = True,
                               take_log            = True,
                               threshold_len_min   = "default",
                               threshold_len_max   = "default",
                               centerline_enabled  = True,
                               enable_target       = True,
                               bins                = None,
                               figsize             = FIGSIZE_HIST,
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
                               axes                = None ):
    
    # Filtering data >>
    df = apply_threshold_len(df, threshold_len_min, threshold_len_max)
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    if residue == "":
        plt.hist( [] )
        return None
    
    # Getting locations of residue locations >>
    locs_aa     = [ get_locations_AA(aa_seq, residue, normalized) for aa_seq in df["seq"] ]
    if take_mean:
        mean_locs_aa = np.array( [ stat.mean(loc_aa) for loc_aa in locs_aa ] )
        hist_data    = mean_locs_aa[~np.isnan(mean_locs_aa)]
    else:
        locs_aa_all  = np.array( list(flatten_list(locs_aa)) )
        hist_data    = locs_aa_all[~np.isnan(locs_aa_all)]
    
    # Taking log of histogram data if requested >>
    if take_log:
        hist_data    = np.log(hist_data)
        # Plotting fit gaussian for histogram >>
        mean_fit, std_fit = norm.fit(hist_data)
        xmin_fit          = mean_fit - 3*std_fit
        xmax_fit          = mean_fit + 3*std_fit
        x_fit             = np.linspace(xmin_fit, xmax_fit, 100)
        y_fit             = norm.pdf(x_fit, mean_fit, std_fit)
        fig.axes.plot(x_fit, y_fit ,   color = COLOR_FIT_GAUSSIAN, linewidth = 1.2)                               # << fit gaussian
        fig.axes.axvline(x = mean_fit, color = COLOR_MEAN_LINE, linewidth = 1.2, label = f"Mean: {mean_fit:.2f}") # << mean for gaussian
        left_conf          = mean_fit - std_fit
        fig.axes.axvline(x = left_conf,  linestyle = "--",   color = COLOR_CONFI_LINE, label = f"Std Dev: {std_fit:.2f}") # << left confidence interval, 1 std
        right_conf         = mean_fit + std_fit
        fig.axes.axvline(x = right_conf, linestyle = "--",   color = COLOR_CONFI_LINE)                            # << right confidence interval, 1 std
        fig.axes.axvspan(left_conf, right_conf, alpha = 0.2, color = COLOR_CONFI_LINE)                            # << color area between confidence interval
    else:
        # Plotting mean & standard error >>
        if enable_target:
            mean_overall = stat.mean(hist_data)
            std_dev      = np.std(hist_data)
            std_err      = std_dev / np.sqrt(len(hist_data))
            fig.axes.axvline(x = mean_overall, color = "red", linewidth = 1.2, label = f"Mean: {mean_overall:.2f}") # << mean of data
            left_conf          = mean_overall - std_err
            fig.axes.axvline(x = left_conf,  linestyle = "--",   color = COLOR_STD_ERR_LINE, label = f"Std Err: {std_err:.2f}") # << left confidence interval, 1 std
            right_conf         = mean_overall + std_err
            fig.axes.axvline(x = right_conf, linestyle = "--",   color = COLOR_STD_ERR_LINE)                    # << right confidence interval, 1 std
            fig.axes.axvspan(left_conf, right_conf, alpha = 0.2, color = COLOR_STD_ERR_LINE)                    # << color area between confidence interval
    
    # Plotting histogram >>
    bins = set_bins(hist_data, bins, xlim)
    fig.axes.hist(hist_data, bins = bins, alpha = 1, density = False, color = COLOR_PLOT_HIST_DEFAULT)
    
    # Plotting center line >>
    if centerline_enabled and normalized and not take_log:
        norm_center = 0.5
        if take_log: norm_center = np.log(norm_center)
        fig.axes.axvline(x = norm_center, color = "black", alpha = 1, label = "center line")
    
    # # Plotting verticle target line for mean location >>
    # if enable_target:
    #     mean_overall = stat.mean(hist_data[~np.isnan(hist_data)])
    #     fig.axes.axvline(x = mean_overall, color = "red", alpha = 1, label = f"Mean: {mean_overall:.2f}")
    
    # Setting plot elements >>
    len_type  = "Normalized" if normalized else "Absolute"
    mean_type = " mean " if take_mean else " individual "
    location  = "locations for each sequence" if take_mean else "locations for all the sequences"
    log       = " (log)" if take_log else ""
    if not title:       title = "Histogram for '" + AA.ABBR_LONG_NAMES[residue] + "'" + mean_type + location
    if not label_xaxis: label_xaxis = len_type + log + " location"
    # if not label_yaxis: label_yaxis = "Density" if normalized else "Count"
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
################################## hist_proteins_residue_loc #######################################################################
####################################################################################################################################


####################################################################################################################################
################################## hist_proteins_residue_relative_count ############################################################
# >>
def hist_proteins_residue_relative_count( df                  = GPCR_DF,
                                          residue             = "S",
                                          relative_count      = True,
                                          take_log            = False,
                                          threshold_len_min   = "default",
                                          threshold_len_max   = "default",
                                          target_type         = "mean",
                                          target_enabled      = True,
                                          targets_x           = (),
                                          targets_x_colors    = (),
                                          targets_x_labels    = (),
                                          targets_y           = (),
                                          targets_y_colors    = (),
                                          targets_y_labels    = (),
                                          confidence_interval = (0.1, 0.9),
                                          bins                = "auto",
                                          figsize             = FIGSIZE_HIST,
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
                                          axes                = None ):
    
    # Filtering data >>
    df = apply_threshold_len(df, threshold_len_min, threshold_len_max)
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    if residue == "":
        plt.hist( [] )
        return None
    
    # Getting count of residue >>
    if relative_count:
        count_aa = [ aa_seq.count(residue) / len(aa_seq) * 100 for aa_seq in df["seq"] ]
        # count_aa = [ aa_seq.count(residue) / len(aa_seq) for aa_seq in df["seq"] if aa_seq.count(residue) !=0 ]
    else:
        count_aa = [ aa_seq.count(residue) for aa_seq in df["seq"] ]
    hist_data = count_aa
    
    # Taking log of histogram data if requested >>
    if take_log:
        hist_data    = np.log(hist_data)
        # Plotting fit gaussian for histogram >>
        mean_fit, std_fit = norm.fit(hist_data)
        left_conf          = mean_fit - std_fit
        right_conf         = mean_fit + std_fit
        if target_enabled:
            xmin_fit          = mean_fit - 3*std_fit
            xmax_fit          = mean_fit + 3*std_fit
            x_fit             = np.linspace(xmin_fit, xmax_fit, 100)
            y_fit             = norm.pdf(x_fit, mean_fit, std_fit)
            fig.axes.plot(x_fit, y_fit ,   color = COLOR_FIT_GAUSSIAN, linewidth = 1.2)                                       # << fit gaussian
            fig.axes.axvline(x = mean_fit, color = COLOR_MEAN_LINE, linewidth = 1.2, label = f"Mean: {mean_fit:.2f}")         # << mean for gaussian
            fig.axes.axvline(x = left_conf,  linestyle = "--",   color = COLOR_CONFI_LINE, label = f"Std Dev: {std_fit:.2f}") # << left confidence interval, 1 std
            fig.axes.axvline(x = right_conf, linestyle = "--",   color = COLOR_CONFI_LINE)                                    # << right confidence interval, 1 std
            fig.axes.axvspan(left_conf, right_conf, alpha = 0.2, color = COLOR_CONFI_LINE)                                    # << color area between confidence interval
        target_center     = mean_fit
        target_left  = left_conf
        target_right = right_conf
    else:
        if target_type == "median":
            target_center = np.quantile(hist_data, 0.5)
            target_left   = np.quantile(hist_data, confidence_interval[0])
            target_right  = np.quantile(hist_data, confidence_interval[1])
            # color_span    = COLOR_QUANTILE_LINE
            color_center  = COLOR_MEDIAN_LINE
            # color_left    = color_span
            # color_right   = color_span
            # label_span    = f"Conf interval: {(confidence_interval[1] - confidence_interval[0])*100:.1f}%"
            label_center  = f"Median: {target_center:.3f}"
            # label_left    = f"Left conf: {confidence_interval[0]*100:.1f}%"
            # label_right   = f"Right conf: {confidence_interval[1]*100:.1f}%"
        elif target_type == "mean":
            target_center = np.mean(hist_data)
            std_err       = np.std(hist_data) / np.sqrt( len(hist_data) )
            if not axes: print( f"Std Err: {std_err}")
            target_left   = target_center - std_err
            target_right  = target_center + std_err
            # color_span    = COLOR_STD_ERR_LINE
            color_center  = COLOR_MEAN_LINE
            # color_left    = color_span
            # color_right   = color_span
            # label_span    = None # f"2 Std Err"
            label_center  = f"Mean: {target_center:.3f}"
            # label_left    = f"+1 Std Err: {target_left:.3f}"
            # label_right   = f"-1 Std Err: {target_right:.3f}"
        if target_enabled:
            # fig.axes.axvspan(target_left, target_right, alpha = 0.15, color = color_span, label = label_span)  # << area between left and right
            fig.axes.axvline(x = target_center, color = color_center, linewidth = 1.2,    label = label_center) # << center
            # fig.axes.axvline(x = target_left,   color = color_left,   linestyle = "--",   label = label_left)   # << left
            # fig.axes.axvline(x = target_right,  color = color_right,  linestyle = "--",   label = label_right)  # << right
            plot_x_targets(fig = fig, figsize = figsize, targets_x = targets_x, targets_x_colors = targets_x_colors, targets_x_labels = targets_x_labels)
            plot_y_targets(fig = fig, figsize = figsize, targets_y = targets_y, targets_y_colors = targets_y_colors, targets_y_labels = targets_y_labels)
    
    # Plotting histogram >>
    bins = set_bins(hist_data, bins, xlim)
    fig.axes.hist( hist_data, bins = bins, alpha = 1, density = False, color = COLOR_PLOT_HIST_DEFAULT)
    
    # Setting plot elements >>
    log      = " (log)" if take_log else ""
    # relative = " relative" if relative_count else ""
    count    = " relative count %" if relative_count else " count"
    # title    = (title if title else "Histogram for '" + AA.ABBR_LONG_NAMES[residue] + "'" + log + relative +" count")
    title    = (title if title else "Histogram for '" + AA.ABBR_LONG_NAMES[residue] + "'" + log + count)
    # if not label_xaxis: label_xaxis = "Amino acid" + log + relative + " count"
    if not label_xaxis: label_xaxis = "Amino acid" + log + count
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

    return fig.figure, target_center, target_left, target_right, hist_data

# <<
################################## hist_proteins_residue_relative_count ############################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/plots/hist_plot_mod ################################################################
####################################################################################################################################
