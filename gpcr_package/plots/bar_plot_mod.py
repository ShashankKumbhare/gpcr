
####################################################################################################################################
################################## gpcr_package/plots/bar_plot_mod #################################################################
####################################################################################################################################

################################################################
from ..__constants__ import *
from ..__auxil__ import *
from ..__data__ import *
################################################################


####################################################################################################################################
################################## bar_gpcr_classes ################################################################################
# >>
def bar_gpcr_classes(   df                  = GPCR_DF,
                        gpcr_classes        = GPCR_CLASSES,
                        threshold           = GPCR_DF.threshold,
                        stacked             = True,
                        normalized          = "class", # can be True or False
                        dataset             = "",
                        target_enabled      = True,
                        targets_x           = (),
                        targets_x_labels    = (),
                        targets_y           = (),
                        targets_y_labels    = (),
                        figsize             = FIGSIZE_BAR,
                        xlim                = None,
                        ylim                = None,
                        title_enabled       = True,
                        title               = "",
                        title_add           = "",
                        label_xaxis_enabled = True,
                        label_xaxis         = "GPCR class",
                        label_yaxis_enabled = True,
                        label_yaxis         = "Seq count",
                        legend_enabled      = True,
                        legend_location     = "upper right",
                        scale_factor        = 1,
                        axes                = None,
                        show                = True ):
    
    # Getting a figure >>
    fig = get_figure(figsize = figsize, axes = axes)
    
    # Recognising dataset >>
    if not dataset: dataset = "GPCR" if df.equals(GPCR_DF) else ( "DisProt" if df.equals(DISPROT_DF) else ( "D2P2" if df.equals(D2P2_DF) else ( "GPCRDB" if df.equals(GPCRDB_DF) else "xxxx" ) ) )
    
    # Counting all, short & long tails for gpcr classes >>
    lens_classes_all_tails_gpcrs  = [ len( df[ df["class"] == gpcr_class ] ) for gpcr_class in gpcr_classes ]
    lens_classes_short_tail_gpcrs = [ len( df[ np.logical_and(df["class"] == gpcr_class, df["seq_len"] <  threshold) ] ) for gpcr_class in gpcr_classes ]
    lens_classes_long_tail_gpcrs  = [ len( df[ np.logical_and(df["class"] == gpcr_class, df["seq_len"] >= threshold) ] ) for gpcr_class in gpcr_classes ]
    
    # Plotting bar plot >>
    x_indices = np.arange( 1 + len(gpcr_classes) )
    width     = 0.25
    
    height_all_tails   = np.array( [sum(lens_classes_all_tails_gpcrs)]  + lens_classes_all_tails_gpcrs )
    height_short_tails = np.array( [sum(lens_classes_short_tail_gpcrs)] + lens_classes_short_tail_gpcrs )
    height_long_tails  = np.array( [sum(lens_classes_long_tail_gpcrs)]  + lens_classes_long_tail_gpcrs )
    
    if normalized == "class":
        height_short_tails = ( height_short_tails / height_all_tails ) * 100
        height_long_tails  = ( height_long_tails  / height_all_tails )  * 100
        height_all_tails   = ( height_all_tails   / height_all_tails )   * 100
    elif normalized:
        height_all_tails   = ( height_all_tails   / len(df) ) * 100
        height_short_tails = ( height_short_tails / len(df) ) * 100
        height_long_tails  = ( height_long_tails  / len(df) ) * 100
        
    if stacked:
        plt.bar(    x      = x_indices,
                    height = height_long_tails + height_short_tails,
                    width  = width,
                    color  = COLOR_LT,
                    label  = 'Long tails'    )
        plt.bar(    x      = x_indices,
                    height = height_short_tails,
                    width  = width,
                    color  = COLOR_ST,
                    label  = 'Short tails'    )
    else:
        plt.bar(    x      = x_indices - width,
                    height = height_all_tails,
                    width  = width,
                    color  = COLOR_AT,
                    label  = 'All tails'    )
        plt.bar(    x      = x_indices,
                    height = height_short_tails,
                    width  = width,
                    color  = COLOR_ST,
                    label  = 'Short tails'    )
        plt.bar(    x      = x_indices + width,
                    height = height_long_tails,
                    width  = width,
                    color  = COLOR_LT,
                    label  = 'Long tails'    )
    
    # Adding x-labels >>
    plt.xticks(ticks = x_indices, labels = ["All"] + gpcr_classes)
    
    # Change no. of tics on y-axis >>
    plt.locator_params(axis = "y", nbins = 15)
    
    # Plotting targets lines >>
    if target_enabled:
        plot_x_targets(fig = fig, figsize = figsize, targets_x = targets_x, targets_x_labels = targets_x_labels)
        plot_y_targets(fig = fig, figsize = figsize, targets_y = targets_y, targets_y_labels = targets_y_labels)
    
    # Setting plot elements >>
    if not title: title = f"{dataset} classes distribution, threshold = {threshold}"
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
                       label_yaxis         = ( "% " if normalized else "" ) + label_yaxis,
                       legend_enabled      = legend_enabled,
                       legend_location     = legend_location,
                       scale_factor        = scale_factor )
    
    # Show plot >>
    plt.show() if show and not axes else plt.close()
    
    return fig.figure
    
# <<
################################## bar_gpcr_classes ################################################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/plots/bar_plot_mod #################################################################
####################################################################################################################################
