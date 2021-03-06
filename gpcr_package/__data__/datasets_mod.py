
####################################################################################################################################
################################## gpcr_package/__data__/datasets_mod ##############################################################
####################################################################################################################################

################################################################
from ..__auxil__.auxil_mod import load_variables_from_file as _load_variables_from_file
from ..__dependencies__.dependencies_mod import pd as _pd
################################################################


####################################################################################################################################
################################## GPCR Dataset ####################################################################################
# >>
# GPCR_DF                   = _pd.read_excel("gpcr_package/__data__/GPCR/GPCR Tail Analysis_compiled_long tails.xlsx")
# GPCR_DF.columns           = ['gene', 'class', 'name', 'seq_len', 'seq']
# GPCR_DF.threshold         = 133
# GPCR_DF.threshold_len_max = None
# GPCR_DF.threshold_len_min = 20
# <<
################################## GPCR Dataset ####################################################################################
####################################################################################################################################


####################################################################################################################################
################################## GPCRDB Dataset ##################################################################################
# >>
GPCRDB_DF                   = _pd.read_excel("gpcr_package/__data__/GPCRDB/GPCRDB.xlsx")
GPCRDB_DF.columns           = ["gene", "class", "name", "seq_len", "helix_8", "seq"]
GPCRDB_DF                   = GPCRDB_DF.drop_duplicates(subset='seq')
GPCRDB_DF                   = GPCRDB_DF.sort_values('seq_len')
GPCRDB_DF.index             = list( range( len(GPCRDB_DF) ) )
GPCRDB_DF.threshold_len_min = 20
GPCRDB_DF.threshold         = 110
GPCRDB_DF.threshold_len_max = None
# <<
################################## GPCRDB Dataset ##################################################################################
####################################################################################################################################


####################################################################################################################################
################################## DisProt Dataset #################################################################################
# >>
DISPROT_DF                   = _load_variables_from_file("gpcr_package\__data__\DisProt\DISPROT_DF.pkl")
DISPROT_DF                   = DISPROT_DF.drop_duplicates(subset='seq')
DISPROT_DF                   = DISPROT_DF.sort_values('seq_len')
DISPROT_DF.index             = list( range( len(DISPROT_DF) ) )
DISPROT_DF.name              = "DisProt"
DISPROT_DF.threshold         = 250
DISPROT_DF.threshold_len_min = 20
DISPROT_DF.threshold_len_max = 600 # 250
# <<
################################## DisProt Dataset #################################################################################
####################################################################################################################################


####################################################################################################################################
################################## D2P2 Dataset ####################################################################################
# >>
D2P2_DF                   = _load_variables_from_file("gpcr_package\__data__\D2P2\D2P2_DF.pkl")
D2P2_DF.name              = "D2P2"
D2P2_DF.threshold         = 300
D2P2_DF.threshold_len_min = 20
D2P2_DF.threshold_len_max = 300
# <<
################################## D2P2 Dataset ####################################################################################
####################################################################################################################################


####################################################################################################################################
################################## MobiDB Dataset ##################################################################################
# >>
# MOBIDB_DF = load_variables_from_file("gpcr_package\__data__\MobiDB\MOBIDB_DF.pkl")
# <<
################################## MobiDB Dataset ##################################################################################
####################################################################################################################################


####################################################################################################################################
################################## Codons Dataset ##################################################################################
# >>
_df_temp                = _pd.read_excel("gpcr_package/__data__/Codons/codons_homo_sapiens.xlsx")
_codons_list            = [ list(codons.split(", ")) for codons in _df_temp["codons"] ]
_df_temp["codons"]      = _codons_list
_freqs_str_list         = [ freqs.split(",") if type(freqs) == str else [freqs] for freqs in _df_temp["freq_codons"] ]
_freqs_float_list       = [ [ float(freq_str) / 1000 for freq_str in freqs_str ] for freqs_str in _freqs_str_list ]
_df_temp["freq_codons"] = _freqs_float_list
_df_temp["freq_aa"]     = [ sum(freq_codons)  for freq_codons in _df_temp["freq_codons"] ]
_df_temp.sum_of_freq    = sum( _df_temp["freq_aa"] )
CODONS_HOMO_SAPIENS_DF  = _df_temp
# <<
################################## MobiDB Dataset ##################################################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/__data__/datasets_mod ##############################################################
####################################################################################################################################
