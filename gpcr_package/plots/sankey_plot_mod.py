
####################################################################################################################################
################################## gpcr_package/plots/sankey #######################################################################
####################################################################################################################################

################################################################
from ..__constants__ import *
from ..__auxil__ import *
################################################################


####################################################################################################################################
################################## sankey_gpcr_classes #############################################################################
# >>
def sankey_gpcr_classes(df = GPCR_DF, gpcr_classes = GPCR_CLASSES, threshold = GPCR_THRESHOLD):

    len_gpcr_classes = len(gpcr_classes)

    lens_classes_gpcrs            = [ len( df[ df["class"] == gpcr_class ] ) for gpcr_class in gpcr_classes ]
    lens_classes_short_tail_gpcrs = [ len( df[ np.logical_and(df["class"] == gpcr_class, df["seq_len"] <  threshold) ] ) for gpcr_class in gpcr_classes ]
    lens_classes_long_tail_gpcrs  = [ len( df[ np.logical_and(df["class"] == gpcr_class, df["seq_len"] >= threshold) ] ) for gpcr_class in gpcr_classes ]
    
    labels        = ["All Classes"] + gpcr_classes + ST_LT_GPCRs
    it_i_labels   = iter( range( len(labels) ) )
    
    first_index   = [next(it_i_labels)]
    class_indices = [i for i in it.islice(it_i_labels, len_gpcr_classes)]
    ST_index      = [next(it_i_labels)]
    LT_index      = [next(it_i_labels)]
    
    colors = [COLOR_GPCR_CLASSES_ALL] + set_colors_for_classes(gpcr_classes) + COLOR_ST_LT
    fig = go.Figure( data = [ go.Sankey( node = dict( pad       = 15,
                                                      thickness = 20,
                                                      line      = dict(color =  "black", width = 0.5),
                                                      label     = labels,
                                                      color     = colors ),
                                         
                                         link = dict( source    = first_index*len_gpcr_classes + class_indices + class_indices,                             # >> list of label indices
                                                      target    = class_indices + ST_index*len_gpcr_classes + LT_index*len_gpcr_classes,                    # >> list of label indices
                                                      value     = lens_classes_gpcrs + lens_classes_short_tail_gpcrs + lens_classes_long_tail_gpcrs) ) ] )  # >> height for each node

    fig.update_layout( title_text = "Short and Long tail distribution of GPCR classes for Threshold = " + str(threshold),
                       font_size  = 15 )
    fig.show()

# <<
################################## sankey_gpcr_classes #############################################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/plots/sankey #######################################################################
####################################################################################################################################
