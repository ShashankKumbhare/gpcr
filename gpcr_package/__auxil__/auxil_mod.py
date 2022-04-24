
####################################################################################################################################
################################## gpcr_package/__auxil__/auxil_mod ################################################################
####################################################################################################################################

################################################################
from ..__constants__ import *
################################################################


####################################################################################################################################
################################## get_random_color ################################################################################
# >>
def get_random_color(num_of_colors):
    import random
    color = [ "#"+''.join([ random.choice('0123456789ABCDEF') for _ in range(6) ]) for _ in range(num_of_colors)]
    return color
# <<
################################## get_random_color ################################################################################
####################################################################################################################################


####################################################################################################################################
################################## get_random_color ################################################################################
# >>
def color_fader(c1 = "red", c2 = "green", c3 = None, mix = 0): # Fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1 = np.array(clr.to_rgb(c1))
    c2 = np.array(clr.to_rgb(c2))
    if not c3:
        A = c1; B = c2; MIX = mix
    else:
        c3 = np.array(clr.to_rgb(c3))
        if   mix <  0.5: A = c1;  B = c2; MIX = 2 * mix
        elif mix >= 0.5: A = c2;  B = c3; MIX = 2 * (mix - 0.5)
    return clr.to_hex( (1-MIX)*A + MIX*B )
# <<
################################## get_random_color ################################################################################
####################################################################################################################################


####################################################################################################################################
################################## set_colors_for_classes ##########################################################################
# >>
def set_colors_for_classes(gpcr_classes = GPCR_CLASSES):
    colors = [""] * len(gpcr_classes)
    for i, gpcr_classe in enumerate(gpcr_classes):
        if gpcr_classe == "A" : colors[i] = COLOR_GPCR_CLASSES[0]
        if gpcr_classe == "B" : colors[i] = COLOR_GPCR_CLASSES[1]
        if gpcr_classe == "C" : colors[i] = COLOR_GPCR_CLASSES[2]
        if gpcr_classe == "F" : colors[i] = COLOR_GPCR_CLASSES[3]
        if gpcr_classe == "AD": colors[i] = COLOR_GPCR_CLASSES[4]
    return colors
# <<
################################## set_colors_for_classes ##########################################################################
####################################################################################################################################


####################################################################################################################################
################################## fasta_to_pandas_df ##############################################################################
# >>
def fasta_to_pandas_df(fasta_file):
    print("Please wait...")
    iter_fasta_file                       = SeqIO.parse(fasta_file, "fasta")       # >> returns an iterator
    iter_fasta_file, copy_iter_fasta_file = it.tee(iter_fasta_file)   # >> using `copy_iter_fasta_file` to get the keys from fasta file
    keys = [ key for key in vars( next(copy_iter_fasta_file) ) ]      # >> getting keys from fasta file
    iter_fasta_file, copy_iter_fasta_file = it.tee(iter_fasta_file)   # >> using `copy_iter_fasta_file` to get the length of records
    len_records                           = len(list( copy_iter_fasta_file ))

    df               = pd.DataFrame( data = [], index = [], columns = [] )
    iters_fasta_file = it.tee(iter_fasta_file, len(keys))             # >> making copies of `iter_fasta_file` for each column of dataframe
    columns = keys                                                    # >> setting column names as keys
    print("Adding columns...")
    for i, key in enumerate(keys):
        print(f"Adding columns {i}...")
        # Removing `_` from the 1st position in keys if present. Ex. `_seq` --> `seq` >>
        if key.startswith("_"): columns[i] = "". join( [ char for j, char in enumerate(key) if j != 0 ] )  # >> modifying column names

        df[columns[i]] = [np.nan] * len_records
        # >> adding columns to dataframe
        if columns[i] == "seq":
            df[columns[i]] = [ str(getattr(record, key)) for record in iters_fasta_file[i] ]
        else:
            df[columns[i]] = [ getattr(record, key) if not(getattr(record, key) == [] or getattr(record, key) == {}) else "" for record in iters_fasta_file[i] ]
        print(f"Columns {i} added....")
    # Adding a column `seq_len` for lengths of sequences >>
    df["seq_len"] = [ len(seq) for seq in df["seq"] ]
    print("Done!")
    return df
# <<
################################## fasta_to_pandas_df ##############################################################################
####################################################################################################################################


####################################################################################################################################
################################## save_variables_to_file ##########################################################################
# >>
def save_variables_to_file(variables, file_pkl):
    with open(file_pkl, 'wb') as f:
        pickle.dump(variables, f)
# <<
################################## save_variables_to_file ##########################################################################
####################################################################################################################################


####################################################################################################################################
################################## load_variables_from_file ########################################################################
# >>
def load_variables_from_file(file_pkl):
    with open(file_pkl, 'rb') as f:
        variables = pickle.load(f)
    return variables
# <<
################################## load_variables_from_file ########################################################################
####################################################################################################################################


####################################################################################################################################
################################## apply_threshold_len #############################################################################
# >>
def apply_threshold_len(df, threshold_len_min = "default", threshold_len_max = "default"):
    if threshold_len_min == "default": threshold_len_min = df.threshold_len_min
    if threshold_len_max == "default": threshold_len_max = df.threshold_len_max
    if threshold_len_min:              df                = df[ df["seq_len"] >= threshold_len_min ]
    if threshold_len_max:              df                = df[ df["seq_len"] <  threshold_len_max ]
    df.index = list( range( len(df) ) )
    return df
# <<
################################## apply_threshold_len #############################################################################
####################################################################################################################################


####################################################################################################################################
################################## get_location_AA #################################################################################
# >>
def get_locations_AA(aa_seq, req_aa, normalized = False, indexing = 1 ):
    norm_len_factor = len(aa_seq) if normalized else 1
    if type(req_aa) == list:
        if indexing == 0: locs = [ i       / norm_len_factor for i, aa in enumerate(aa_seq) if aa in req_aa ]
        if indexing == 1: locs = [ (i + 1) / norm_len_factor for i, aa in enumerate(aa_seq) if aa in req_aa ]
    else:
        if indexing == 0: locs = [ i       / norm_len_factor for i, aa in enumerate(aa_seq) if aa == req_aa ]
        if indexing == 1: locs = [ (i + 1) / norm_len_factor for i, aa in enumerate(aa_seq) if aa == req_aa ]
    
    if len(locs) == 0: locs = [np.nan]
    return locs
# <<
################################## get_location_AA #################################################################################
####################################################################################################################################


####################################################################################################################################
################################## flatten_list ####################################################################################
# >>
def flatten_list(some_list):
    for element in some_list:
        if type(element) in (tuple, list):
            for item in flatten_list(element):
                yield item
        else:
            yield element
# <<
################################## flatten_list ####################################################################################
####################################################################################################################################


####################################################################################################################################
################################## set_fontsize ####################################################################################
# >>
def set_fontsize(figsize, string_type = "none", scale_factor_extended = 1 ):
    if string_type   == "title" or string_type == "gridtitle":
        if string_type   == "title":
            scale_factor = SCALE_FACTOR_FONT_TITLE
        elif string_type == "gridtitle":
            scale_factor = SCALE_FACTOR_FONT_TITLE_GRID
        else:
            scale_factor = SCALE_FACTOR_FONT_TITLE
        fontsize = min(figsize) * scale_factor
    else:
        if string_type   == "axes_label":
            scale_factor = SCALE_FACTOR_FONT_AXES_LABEL
        elif string_type == "legend":
            scale_factor = SCALE_FACTOR_FONT_LEGEND
        elif string_type == "axes_ticks" or "annotation":
            scale_factor = SCALE_FACTOR_FONT_AXES_TICKS
        else:
            scale_factor = SCALE_FACTOR_FONT_DEFAULT
        fontsize = min(figsize) * scale_factor
    return fontsize * scale_factor_extended
# <<
################################## set_fontsize ####################################################################################
####################################################################################################################################


####################################################################################################################################
################################## set_linewidth ###################################################################################
# >>
def set_linewidth(figsize, scale_factor = SCALE_FACTOR_LINE_DEFAULT):
    linewidth = min(figsize) * scale_factor
    return linewidth
# <<
################################## set_linewidth ###################################################################################
####################################################################################################################################


####################################################################################################################################
################################## set_bins ########################################################################################
# >>
def set_bins(hist_data, bins, xlim):
    if not bins:
        if xlim:
            xlim_bins_min = xlim[0]
            xlim_bins_max = xlim[1]
        else:
            xlim_bins_min = min(hist_data)
            xlim_bins_max = max(hist_data)
        bins = np.linspace(xlim_bins_min, xlim_bins_max, num = BINS_PLOT_HIST + 1)
    elif type(bins) == np.ndarray:
        bins = list(bins)
    return bins
# <<
################################## set_bins ########################################################################################
####################################################################################################################################


####################################################################################################################################
################################## get_figure ######################################################################################
# >>
def get_figure(figsize = None, axes = None):
    fig = Struct()
    if axes:
        fig.axes   = axes
        fig.figure = axes.figure
        # fig.figure = axes[0].figure
    else:
        fig.figure = plt.figure(figsize = figsize)
        fig.axes   = plt.gca()
    return fig
# <<
################################## get_figure ######################################################################################
####################################################################################################################################


####################################################################################################################################
################################## set_plot_elements ###############################################################################
# >>
def set_plot_elements( fig                 = None,
                       figsize             = FIGSIZE_DEFAULT,
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
                       legend_location     = None,
                       scale_factor        = 1 ):
    figsize = tuple([scale_factor*figsize[0], scale_factor*figsize[1]])
    if title_enabled:       fig.axes.set_title(title + title_add, fontsize = set_fontsize(figsize, "title"))
    if label_xaxis_enabled: fig.axes.set_xlabel(label_xaxis,      fontsize = set_fontsize(figsize, "axes_label"))
    if label_yaxis_enabled: fig.axes.set_ylabel(label_yaxis,      fontsize = set_fontsize(figsize, "axes_label"))
    if legend_enabled and fig.axes.get_legend_handles_labels()[0]:
        fig.axes.legend(loc = legend_location, fontsize = set_fontsize(figsize, "legend"))
    if xlim: fig.axes.set_xlim(xlim)
    if ylim: fig.axes.set_ylim(ylim)
    fig.axes.tick_params(axis = 'x', labelsize = set_fontsize(figsize, "axes_ticks"), rotation = 0)
    fig.axes.tick_params(axis = 'y', labelsize = set_fontsize(figsize, "axes_ticks"))
# <<
################################## set_plot_elements ###############################################################################
####################################################################################################################################


####################################################################################################################################
################################## plot_x_targets ##################################################################################
# >>
def plot_x_targets( fig              = None,
                    figsize          = FIGSIZE_DEFAULT,
                    targets_x        = (),
                    targets_x_colors = (),
                    targets_x_labels = () ):
    if type(targets_x)        == int or type(targets_x)== float: targets_x        = [targets_x]
    if type(targets_x_colors) == str:                            targets_x_colors = [targets_x_colors]
    if type(targets_x_labels) == str:                            targets_x_labels = [targets_x_labels]
    if not targets_x_colors: targets_x_colors = get_random_color( len(targets_x) )
    for target_x, target_x_color, target_x_label in it.zip_longest(targets_x, targets_x_colors, targets_x_labels):
        fig.axes.axvline( target_x, color = target_x_color, label = target_x_label, linewidth = set_linewidth(figsize) )
# <<
################################## plot_x_targets ##################################################################################
####################################################################################################################################


####################################################################################################################################
################################## plot_y_targets ##################################################################################
# >>
def plot_y_targets( fig              = None,
                    figsize          = FIGSIZE_DEFAULT,
                    targets_y        = (),
                    targets_y_colors = (),
                    targets_y_labels = () ):
    if type(targets_y)        == int or type(targets_y)== float: targets_y        = [targets_y]
    if type(targets_y_colors) == str:                            targets_y_colors = [targets_y_colors]
    if type(targets_y_labels) == str:                            targets_y_labels = [targets_y_labels]
    if not targets_y_colors: targets_y_colors = get_random_color( len(targets_y) )
    for target_y, target_y_color, target_y_label in it.zip_longest(targets_y, targets_y_colors, targets_y_labels):
        fig.axes.axhline( target_y, color = target_y_color, label = target_y_label, linewidth = set_linewidth(figsize) )
# <<
################################## plot_y_targets ##################################################################################
####################################################################################################################################


####################################################################################################################################
################################## get_colorbar_heatmap ############################################################################
# >>
def get_colorbar_heatmap(fig):
# def get_colorbar_heatmap(fig, vmin = None, vmax = None):
    PCM = None
    for child in fig.axes.get_children():
        if isinstance(child, plt.cm.ScalarMappable):
            PCM = child
            break
    fig.figure.colorbar(PCM, ax = fig.axes)
    # fig.axes.set_clim(vmin = vmin, vmax = vmax)
# <<
################################## get_colorbar_heatmap ############################################################################
####################################################################################################################################


####################################################################################################################################
################################## count_overlap ###################################################################################
# >>
def count_overlap(  seq        = None,
                    char_start = "",
                    char_next  = "",
                    lag        = 1 ):
    if len(seq) > lag:
        count = sum([ 1 for i in range(len(seq) - lag) if seq[i] == char_start and (seq[i + lag] == char_next if char_next else True) ])
    else:
        count = 0
    return count
# <<
################################## count_overlap ###################################################################################
####################################################################################################################################


####################################################################################################################################
################################## dict_to_df ######################################################################################
# >>
def dict_to_df(dct):
        dataframe = pd.DataFrame()
        for key_1, dict_inside in dct.items():
            dataframe[key_1] = dict_inside.values()
        dataframe.index = dct.keys()
        return dataframe.transpose()
# <<
################################## dict_to_df ######################################################################################
####################################################################################################################################



####################################################################################################################################
################################## gpcr_package/__auxil__/auxil_mod ################################################################
####################################################################################################################################
