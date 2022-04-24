

####################################################################################################################################
################################## Disprot Dataset #################################################################################
# >>
# Fasta file to dataframe >>
DISPROT_DF = gpcr.fasta_to_pandas_df("gpcr_package\__data__\DisProt\disprot_homo_sapiens.fasta")
# Fasta file to dataframe >>
gpcr.save_variables_to_file(DISPROT_DF, "gpcr_package\__data__\DisProt\DISPROT_DF.pkl")
# Loading dataframe from .pkl file >>
DISPROT_DF = gpcr.load_variables_from_file("gpcr_package\__data__\DisProt\DISPROT_DF.pkl")
# <<
################################## Disprot Dataset #################################################################################
####################################################################################################################################


####################################################################################################################################
################################## D2P2 Dataset ####################################################################################
# >>
# Fasta file to dataframe >>
D2P2_DF = gpcr.fasta_to_pandas_df("gpcr_package\__data__\D2P2\d2p2_homo_sapiens.fasta")
# Saving dataframe to .pkl file >>
gpcr.save_variables_to_file(D2P2_DF, "gpcr_package\__data__\D2P2\D2P2_DF.pkl")
# Loading dataframe from .pkl file >>
D2P2_DF = gpcr.load_variables_from_file("gpcr_package\__data__\D2P2\D2P2_DF.pkl")
# <<
################################## D2P2 Dataset ####################################################################################
####################################################################################################################################


####################################################################################################################################
################################## MobiDB Dataset ##################################################################################
# >>
# Fasta file to dataframe >>
MOBIDB_DF = gpcr.fasta_to_pandas_df("gpcr_package\__data__\MobiDB\mobidb_homo_sapiens.fasta")
# Fasta file to dataframe >>
gpcr.save_variables_to_file(MOBIDB_DF, "gpcr_package\__data__\MobiDB\MOBIDB_DF.pkl")
# Loading dataframe from .pkl file >>
MOBIDB_DF = gpcr.load_variables_from_file("gpcr_package\__data__\MobiDB\MOBIDB_DF.pkl")
# <<
################################### MobiDB Dataset #################################################################################
####################################################################################################################################
