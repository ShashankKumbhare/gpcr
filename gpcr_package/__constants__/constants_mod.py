
####################################################################################################################################
################################## gpcr_package/__constants__/constants_mod ########################################################
####################################################################################################################################

################################################################
from ..__dependencies__ import *
################################################################

class Struct: pass

####################################################################################################################################
################################## GPCR related constants ##########################################################################
# >>
GPCR_DF                   = pd.read_excel("gpcr_package/__data__/GPCR/GPCR Tail Analysis_compiled_long tails.xlsx")
GPCR_DF.columns           = ["gene", "class", "name", "seq_len", "seq"]
GPCR_DF                   = GPCR_DF.drop_duplicates(subset='seq')
GPCR_DF                   = GPCR_DF.sort_values('seq_len')
GPCR_DF.index             = list( range( len(GPCR_DF) ) )
# GPCR_DF.name              = "GPCR"
GPCR_DF.threshold_len_min = 20
GPCR_DF.threshold         = 110
GPCR_DF.threshold_len_max = None
LEN_GPCR_DF               = len(GPCR_DF)

AA_LONG_NAMES          = [ 'Arginine', 'Lysine', 'Histidine', 'Aspartic acid', 'Glutamic acid', 'Asparagine', 'Glutamine', 'Cysteine', 'Serine', 'Threonine', 'Tyrosine', 'Glycine', 'Alanine', 'Proline', 'Valine', 'Leucine', 'Isoleucine', 'Methionine', 'Tryptophan', 'Phenylalanine' ]
AA_SHORT_NAMES         = [ 'Arg',      'Lys',    'His',       'Asp',           'Glu',           'Asn',        'Gln',       'Cys',      'Ser',    'Thr',       'Tyr',      'Gly',     'Ala',     'Pro',     'Val',    'Leu',     'Ile',        'Met',        'Trp',        'Phe' ]
AA_ABBREVATIONS        = [ 'R',        'K',      'H',         'D',             'E',             'N',          'Q',         'C',        'S',      'T',         'Y',        'G',       'A',       'P',       'V',      'L',       'I',          'M',          'W',          'F' ]
AA_HYDRO_PHOBICITY_PH2 = [ -26,        -37,      -42,         -18,             8,               -41,          8,           52,         -7,       13,          49,         0,         47,        -46,       79,       100,       100,          74,           84,           92 ]
AA_HYDRO_PHOBICITY_PH7 = [ -14,        -23,      8,           -55,             -31,             -28,          -31,         49,         -5,       13,          63,         0,         41,        -46,       76,       97,        99,           74,           97,           100 ]
AA_CHARGES_PH7         = [ 1,           1,       1,           -1,              -1,              0,            0,           0,          0,        0,           0,          0,          0,         0,         0,        0,         0,            0,            0,            0 ]
AA_ABBRE_SHORT_NAMES   = [ aa_abbrevation + ": " + aa_short_name for aa_abbrevation, aa_short_name in zip(AA_ABBREVATIONS, AA_SHORT_NAMES) ]
AA_ABBRE_LONG_NAMES    = [ aa_abbrevation + ": " + aa_long_name  for aa_abbrevation, aa_long_name  in zip(AA_ABBREVATIONS, AA_LONG_NAMES) ]
AA_ACIDIC              = [ 'D', 'Q' ]
AA_HYDROPHIBOC         = [ 'V', 'L', 'I', 'M', 'F', 'W', 'C' ]
AA_HYDROPHIBOC_1       = [ 'V', 'L', 'I', 'M', 'F', 'W', 'A', 'Y' ] # https://www.youtube.com/watch?v=cL2_e83v3js&t=13s
AA_HYDROPHIBOC_2       = [ 'V', 'L', 'I', 'M', 'F', 'W', 'A', 'G', 'P' ] # https://www2.chem.wisc.edu/deptfiles/genchem/netorial/modules/biomolecules/modules/protein1/prot13.htm#:~:text=Hydrophobic%20Amino%20Acids&text=The%20nine%20amino%20acids%20that,is%20the%20structure%20of%20valine.
AA_HYDROPHIBOC_3       = [ 'V', 'L', 'I']
AA_HYDROPHIBOC_4       = [ 'V', 'L', 'I', 'M', 'F', 'W', 'Y']          # Paper: Lee, PDZ domains and their binding partners structure

AA                     = Struct()
AA.SHORT_NAMES         = dict(zip(AA_ABBREVATIONS, AA_SHORT_NAMES))
AA.LONG_NAMES          = dict(zip(AA_ABBREVATIONS, AA_LONG_NAMES))
AA.ABBREVATIONS        = dict(zip(AA_ABBREVATIONS, AA_ABBREVATIONS))
AA.ABBR_SHORT_NAMES    = dict(zip(AA_ABBREVATIONS, AA_ABBRE_SHORT_NAMES))
AA.ABBR_LONG_NAMES     = dict(zip(AA_ABBREVATIONS, AA_ABBRE_LONG_NAMES))
AA.HYDRO_PHOBICITY_PH2 = dict(zip(AA_ABBREVATIONS, AA_HYDRO_PHOBICITY_PH2))
AA.HYDRO_PHOBICITY_PH7 = dict(zip(AA_ABBREVATIONS, AA_HYDRO_PHOBICITY_PH7))
AA.CHARGES_PH7         = dict(zip(AA_ABBREVATIONS, AA_CHARGES_PH7))

GPCR_CLASSES           = ["A",       "B",       "C",       "F",       "AD"]
COLOR_GPCR_CLASSES     = ["#31A354", "#3182BD", "#20B2AA", "#FF8745", "#E978AA"]
COLOR_GPCR_CLASSES_ALL = "#074457"

ST_LT                  = ["Short Tail", "Long Tail"]
ST_LT_GPCRs            = ["Short Tail GPCRs", "Long Tail GPCRs"]
COLOR_ST_LT            = ["#A1A1A1", "#3C3C3C"]

GPCR_THRESHOLD                        = GPCR_DF.threshold
LENS_CLASSES_ALL_GPCRS                = [ len( GPCR_DF[ GPCR_DF["class"] == gpcr_class ] ) for gpcr_class in GPCR_CLASSES ]
LENS_CLASSES_SHORT_TAIL_ALL_GPCRS     = [ len( GPCR_DF[ np.logical_and(GPCR_DF["class"] == gpcr_class, GPCR_DF["seq_len"] <  GPCR_THRESHOLD) ] ) for gpcr_class in GPCR_CLASSES ]
LENS_CLASSES_LONG_TAIL_ALL_GPCRS      = [ len( GPCR_DF[ np.logical_and(GPCR_DF["class"] == gpcr_class, GPCR_DF["seq_len"] >= GPCR_THRESHOLD) ] ) for gpcr_class in GPCR_CLASSES ]
# <<
################################## GPCR related constants ##########################################################################
####################################################################################################################################


####################################################################################################################################
################################## plots related constants #########################################################################
# >>

# Bin size >>
BINS_PLOT_HIST       = 50

# Figure sizes >>
FIGSIZE_HIST_DEFAULT = (8, 6)
FIGSIZE_HIST         = (6, 6)
FIGSIZE_LINE         = (10, 8)
FIGSIZE_SCATTER      = (6, 6)
FIGSIZE_BOX          = (1, 6)
FIGSIZE_HEATMAP      = (11, 8)
FIGSIZE_DEFAULT      = (10, 6)
FIGSIZE_BAR          = (8, 6)
SUBFIGSIZE_GRID      = (4, 4)

# Fontsize Scale factor >>
SCALE_FACTOR_FONT_TITLE      = 3
SCALE_FACTOR_FONT_TITLE_GRID = 1.5
SCALE_FACTOR_FONT_TITLE_BOX  = 1
SCALE_FACTOR_FONT_LEGEND     = 2
SCALE_FACTOR_FONT_AXES_LABEL = 2.5
SCALE_FACTOR_FONT_AXES_TICKS = 2
SCALE_FACTOR_FONT_DEFAULT    = 2.5

# Linewidth Scale factor >>
SCALE_FACTOR_LINE_DEFAULT    = 1/5

# Colors >>
COLOR_AT                = "#4D4B4D"
COLOR_ST                = "#1796DB"
COLOR_LT                = "#EBB240"
COLOR_PLOT_HIST_DEFAULT = "lightseagreen"
COLOR_PLOT_SCATTER      = "#15517B"
COLOR_PLOT_LINE         = "#65A7D0"
COLOR_FIT_GAUSSIAN      = "red"
COLOR_MEAN_LINE         = "red"
COLOR_STD_ERR_LINE      = "red"
COLOR_MEDIAN_LINE       = "darkviolet"
COLOR_QUANTILE_LINE     = "darkviolet"
COLOR_CONFI_LINE        = "#247390"
COLOR_COLORBAR_CMAP           = clr.LinearSegmentedColormap.from_list('custom blue', ["#02FDFF","#D701D9"], N = 256)
COLOR_COLORBAR_CMAP_RED_BLUE  = clr.LinearSegmentedColormap.from_list('custom blue', ["#0005F1","#FB0101"], N = 256)
COLOR_SEQ_COMPARISION_SEQ_TITLE              = "yellow"
COLOR_SEQ_COMPARISION_MATCH                  = "green"
COLOR_SEQ_COMPARISION_MISMATCH               = "red"
COLOR_SEQ_COMPARISION_LENGTH_COMPARISION     = "magenta"
COLOR_SEQ_COMPARISION_LENGTH_NON_COMPARISION = "white"
COLOR_SEQ_COMPARISION_MATCH_RATIO            = "blue"

# Divider line >>
divider_line = colored(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", "blue")

# <<
################################## plots related constants #########################################################################
####################################################################################################################################



####################################################################################################################################
################################## gpcr_package/__constants__/constants_mod ########################################################
####################################################################################################################################
