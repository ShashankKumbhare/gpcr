
####################################################################################################################################
################################## gpcr_package/plots/sequence_comparision_mod #####################################################
####################################################################################################################################

################################################################
from ..__constants__ import *
from ..__auxil__ import *
################################################################


####################################################################################################################################
################################## print_comparision_seqs ##########################################################################
# >>
def print_seqs_comparision( seq_1               = "1ST_TEST_SEQ",
                            seq_2               = "2ND_TEST_SEQUENCE",
                            start_seq_1         = 0,
                            start_seq_2         = 0,
                            length_comparision  = None,
                            index_seq_1         = 1,
                            index_seq_2         = 2 ):
    
    # Setting default length for comparision: which is length of the smaller sequence >>
    if length_comparision is None: length_comparision = min([len(seq_1) - start_seq_1, len(seq_2) - start_seq_2])
    
    # AAs before comparision length >>
    out_seq_1 = colored( "".join([" "] * ((start_seq_2 - start_seq_1) if start_seq_2 > start_seq_1 else 0)) + seq_1[0: start_seq_1], COLOR_SEQ_COMPARISION_LENGTH_NON_COMPARISION)
    out_seq_2 = colored( "".join([" "] * ((start_seq_1 - start_seq_2) if start_seq_1 > start_seq_2 else 0)) + seq_2[0: start_seq_2], COLOR_SEQ_COMPARISION_LENGTH_NON_COMPARISION)
    
    # AAs in comparision length >>
    for i in range(length_comparision):
        if seq_1[start_seq_1 + i] == seq_2[start_seq_2 + i]: color = COLOR_SEQ_COMPARISION_MATCH    # >> green color for match
        else:                                                color = COLOR_SEQ_COMPARISION_MISMATCH # >> red color for mis-match
        out_seq_1 = out_seq_1 + colored(seq_1[start_seq_1 + i], color)
        out_seq_2 = out_seq_2 + colored(seq_2[start_seq_2 + i], color)
        
    # AAs after comparision length >>
    out_seq_1 = out_seq_1 + colored(seq_1[start_seq_1 + length_comparision:], COLOR_SEQ_COMPARISION_LENGTH_NON_COMPARISION)
    out_seq_2 = out_seq_2 + colored(seq_2[start_seq_2 + length_comparision:], COLOR_SEQ_COMPARISION_LENGTH_NON_COMPARISION)
    
    # Printing comparision >>
    aa_pos_1 = colored(f"(aa_loc = {str(start_seq_1).ljust(4)} - {str(start_seq_1 + length_comparision - 1).ljust(4)})", COLOR_SEQ_COMPARISION_LENGTH_COMPARISION)
    aa_pos_2 = colored(f"(aa_loc = {str(start_seq_2).ljust(4)} - {str(start_seq_2 + length_comparision - 1).ljust(4)})", COLOR_SEQ_COMPARISION_LENGTH_COMPARISION)
    print( f"{colored(f'Seq {str(index_seq_1).ljust(5)}', COLOR_SEQ_COMPARISION_SEQ_TITLE)} {aa_pos_1} : " + out_seq_1)
    print( f"{colored(f'Seq {str(index_seq_2).ljust(5)}', COLOR_SEQ_COMPARISION_SEQ_TITLE)} {aa_pos_2} : " + out_seq_2)

    return None
# <<
################################## print_comparision_seqs ##########################################################################
####################################################################################################################################


####################################################################################################################################
################################## print_details_comparision_seqs ##################################################################
def print_details_comparision_seqs( match,
                                    mismatch,
                                    length_comparision,
                                    length_non_comparision,
                                    match_ratio_comparision_length,
                                    match_ratio_overall ):
    
    print( f"Match                            : {colored(f'{match}',                                COLOR_SEQ_COMPARISION_MATCH)}" )
    print( f"Mis-match                        : {colored(f'{mismatch}',                             COLOR_SEQ_COMPARISION_MISMATCH)}" )
    print( f"Comparision length               : {colored(f'{length_comparision}',                   COLOR_SEQ_COMPARISION_LENGTH_COMPARISION)}")
    print( f"Non comparision length           : {colored(f'{length_non_comparision}',               COLOR_SEQ_COMPARISION_LENGTH_NON_COMPARISION)}")
    print( f"Comparision length match-ratio   : {colored(f'{match_ratio_comparision_length:0.10f}', COLOR_SEQ_COMPARISION_MATCH_RATIO)}     [ Comparision length match-ratio =  Match  /  Comparision length" )
    print( f"Overall length match-ratio       : {colored(f'{match_ratio_overall           :0.10f}', COLOR_SEQ_COMPARISION_MATCH_RATIO)}     [ Overall length match-ratio     = 2*Match / ( len(seq_1) + len(seq_2) ) ]" )
    print(divider_line)
    print(divider_line)

################################## print_details_comparision_seqs ##################################################################
####################################################################################################################################


####################################################################################################################################
################################## get_match_ratio #################################################################################
def get_match_ratio (   seq_1               = "1ST_TEST_SEQ",
                        seq_2               = "2ND_TEST_SEQUENCE",
                        start_seq_1_list    = 0, # >> can be a list
                        start_seq_2_list    = 0, # >> can be a list
                        length_comparision  = None,
                        index_seq_1         = 1,
                        index_seq_2         = 2,
                        show_details        = False,
                        iscycle             = False
                    ):
    
    if type(start_seq_1_list) == int: start_seq_1_list = [start_seq_1_list]
    if type(start_seq_2_list) == int: start_seq_2_list = [start_seq_2_list]
    length_comparision_is_provided = True if length_comparision else False
    
    # Creating a dict for comparision >>
    match_score_dict                                   = dict()
    match_score_dict["index_seq_1"]                    = []
    match_score_dict["index_seq_2"]                    = []
    match_score_dict["start_seq_1"]                    = []
    match_score_dict["start_seq_2"]                    = []
    match_score_dict["match"]                          = []
    match_score_dict["mismatch"]                       = []
    match_score_dict["length_comparision"]             = []
    match_score_dict["length_non_comparision"]         = []
    match_score_dict["match_ratio_comparision_length"] = []
    match_score_dict["match_ratio_overall"]            = []
    match_score_dict["seq_1"]                          = []
    match_score_dict["seq_2"]                          = []
    
    if show_details:
        print(divider_line)
        print(divider_line)
    for start_seq_1 in start_seq_1_list:
        for start_seq_2 in start_seq_2_list:
            
            # Setting default length for comparision: which is length of the smaller sequence >>
            if not length_comparision_is_provided: length_comparision = min([len(seq_1) - start_seq_1, len(seq_2) - start_seq_2])
            
            # Calculating match and mismatch >>
            match    = 0
            mismatch = 0
            for i in range(length_comparision):
                if seq_1[start_seq_1 + i] == seq_2[start_seq_2 + i]:
                    match    = match + 1
                else:
                    mismatch = mismatch + 1
            
            # Calculating match ratios >>
            length_non_comparision         = len(seq_1) + len(seq_2) - 2 * length_comparision
            match_ratio_comparision_length = match / length_comparision
            match_ratio_overall            = 2 * match / ( len(seq_1) + len(seq_2) )
            if show_details:
                print_seqs_comparision( seq_1               = seq_1,
                                        seq_2               = seq_2,
                                        start_seq_1         = start_seq_1,
                                        start_seq_2         = start_seq_2,
                                        length_comparision  = length_comparision,
                                        index_seq_1         = index_seq_1,
                                        index_seq_2         = index_seq_2 )
                print_details_comparision_seqs( match,
                                                mismatch,
                                                length_comparision,
                                                length_non_comparision,
                                                match_ratio_comparision_length,
                                                match_ratio_overall )
            
            # Append to output dict >>
            match_score_dict["index_seq_1"].append(index_seq_1)
            match_score_dict["index_seq_2"].append(index_seq_2)
            match_score_dict["start_seq_1"].append(start_seq_1)
            match_score_dict["start_seq_2"].append(start_seq_2)
            match_score_dict["match"].append(match)
            match_score_dict["mismatch"].append(mismatch)
            match_score_dict["length_comparision"].append(length_comparision)
            match_score_dict["length_non_comparision"].append(length_non_comparision)
            match_score_dict["match_ratio_comparision_length"].append(match_ratio_comparision_length)
            match_score_dict["match_ratio_overall"].append(match_ratio_overall)
            match_score_dict["seq_1"].append(seq_1)
            match_score_dict["seq_2"].append(seq_2)
            
    # Making a dataframe from comparision dict >>
    match_score_df = pd.DataFrame.from_dict(match_score_dict)
    
    # Printing the comparision for optimal overall length match-ratio >>
    optimal_match_ratio_overall_df = match_score_df[ match_score_df['match_ratio_overall'] == max(match_score_df['match_ratio_overall']) ]
    if not iscycle and show_details:
        cprint("                    Optimal overall length match-ratio", "red")
        print(divider_line)
        print(divider_line)
        for index, row_optimal_match_ratio_overall in optimal_match_ratio_overall_df.iterrows():
            print_seqs_comparision( seq_1               = seq_1,
                                    seq_2               = seq_2,
                                    start_seq_1         = int(row_optimal_match_ratio_overall["start_seq_1"]),
                                    start_seq_2         = int(row_optimal_match_ratio_overall["start_seq_2"]),
                                    length_comparision  = int(row_optimal_match_ratio_overall["length_comparision"]),
                                    index_seq_1         = int(row_optimal_match_ratio_overall["index_seq_1"]),
                                    index_seq_2         = int(row_optimal_match_ratio_overall["index_seq_2"]) )
            print_details_comparision_seqs( match                           = row_optimal_match_ratio_overall["match"],
                                            mismatch                        = row_optimal_match_ratio_overall["mismatch"],
                                            length_comparision              = row_optimal_match_ratio_overall["length_comparision"],
                                            length_non_comparision          = row_optimal_match_ratio_overall["length_non_comparision"],
                                            match_ratio_comparision_length  = row_optimal_match_ratio_overall["match_ratio_comparision_length"],
                                            match_ratio_overall             = row_optimal_match_ratio_overall["match_ratio_overall"] )
        
    return match_score_df, optimal_match_ratio_overall_df
################################## get_match_ratio #################################################################################
####################################################################################################################################


####################################################################################################################################
################################## get_match_ratio_cycle ###########################################################################
def get_match_ratio_cycle(  seq_1               = "1ST_TEST_SEQ",
                            seq_2               = "2ND_TEST_SEQUENCE",
                            # length_comparision  = None,
                            index_seq_1         = 1,
                            index_seq_2         = 2,
                            show_details        = False
                         ):
    
    match_score_half_cycle_1_df, _ = get_match_ratio(   seq_1               = seq_1,
                                                        seq_2               = seq_2,
                                                        # length_comparision  = length_comparision,
                                                        start_seq_1_list    = 0,
                                                        start_seq_2_list    = range(len(seq_2) - 1, -1, -1),
                                                        index_seq_1         = index_seq_1,
                                                        index_seq_2         = index_seq_2,
                                                        show_details        = show_details,
                                                        iscycle             = True   )
    
    match_score_half_cycle_2_df, _ = get_match_ratio(   seq_1               = seq_1,
                                                        seq_2               = seq_2,
                                                        # length_comparision  = length_comparision,
                                                        start_seq_1_list    = range(1, len(seq_1)),
                                                        start_seq_2_list    = 0,
                                                        index_seq_1         = index_seq_1,
                                                        index_seq_2         = index_seq_2,
                                                        show_details        = show_details,
                                                        iscycle             = True   )
    
    match_score_cycle    = [match_score_half_cycle_1_df, match_score_half_cycle_2_df]
    match_score_cycle_df = pd.concat(match_score_cycle)
    
    # Printing the comparision for optimal overall length match-ratio >>
    optimal_match_ratio_overall_df = match_score_cycle_df[ match_score_cycle_df['match_ratio_overall'] == max(match_score_cycle_df['match_ratio_overall']) ]
    if show_details:
        cprint("                    Optimal overall length match-ratio", "green")
        print(divider_line)
        print(divider_line)
        for index, row_optimal_match_ratio_overall in optimal_match_ratio_overall_df.iterrows():
            print_seqs_comparision( seq_1               = seq_1,
                                    seq_2               = seq_2,
                                    start_seq_1         = int(row_optimal_match_ratio_overall["start_seq_1"]),
                                    start_seq_2         = int(row_optimal_match_ratio_overall["start_seq_2"]),
                                    length_comparision  = int(row_optimal_match_ratio_overall["length_comparision"]),
                                    index_seq_1         = int(row_optimal_match_ratio_overall["index_seq_1"]),
                                    index_seq_2         = int(row_optimal_match_ratio_overall["index_seq_2"]) )
            print_details_comparision_seqs( match                           = row_optimal_match_ratio_overall["match"],
                                            mismatch                        = row_optimal_match_ratio_overall["mismatch"],
                                            length_comparision              = row_optimal_match_ratio_overall["length_comparision"],
                                            length_non_comparision          = row_optimal_match_ratio_overall["length_non_comparision"],
                                            match_ratio_comparision_length  = row_optimal_match_ratio_overall["match_ratio_comparision_length"],
                                            match_ratio_overall             = row_optimal_match_ratio_overall["match_ratio_overall"] )
    
    return match_score_cycle_df, optimal_match_ratio_overall_df

################################## get_match_ratio_cycle ###########################################################################
####################################################################################################################################


####################################################################################################################################
################################## reduce_dataset_based_on_match_ratio #############################################################
def reduce_dataset_based_on_match_ratio(    df                                  = GPCR_DF,
                                            sequences                           = (),
                                            threshold_len_min                   = "default",
                                            threshold_len_max                   = "default",
                                            limit_match_score                   = 0.8,
                                            num_of_next_seqs_to_compare_with    = 50,
                                            show_details                        = False
                                         ):
    
    # Filtering data >>
    if not type(sequences) == pd.core.series.Series:
        if not sequences:
            df                      = apply_threshold_len(df, threshold_len_min, threshold_len_max)
            sequences               = list(df["seq"])
    
    # Calculating optimal match ratios for different positions for each sequence with the next `num_of_next_seqs_to_compare_with` no. of sequences >>
    # optimal_match_ratios_overall_df = [np.nan] * (len(sequences) - 1)
    n = num_of_next_seqs_to_compare_with
    optimal_match_ratios_overall_df = [np.nan for i in range((len(sequences) - 1)) for _ in ( range(i+1, i+1+n) if (i+1+n) < len(sequences) else range(i+1, len(sequences)) ) ]
    pbar = ProgressBar()
    df_k = -1
    for i in pbar(range(len(sequences) - 1)):
        for j in ( range(i+1, i+1+n) if (i+1+n) < len(sequences) else range(i+1, len(sequences)) ):
            # _, optimal_match_ratios_overall_df[i] = get_match_ratio_cycle(  seq_1               = sequences[i],
            #                                                                 seq_2               = sequences[i+1],
            #                                                                 # length_comparision  = None,
            #                                                                 index_seq_1         = i,
            #                                                                 index_seq_2         = i+1,
            #                                                                 show_details        = False
            #                                                               )
            df_k = df_k + 1
            _, optimal_match_ratios_overall_df[df_k] = get_match_ratio_cycle(   seq_1               = sequences[i],
                                                                                seq_2               = sequences[j],
                                                                                # length_comparision  = None,
                                                                                index_seq_1         = i,
                                                                                index_seq_2         = j,
                                                                                show_details        = False
                                                                              )
    
    # Clubbingg all optimal match ratios together >>
    optimal_match_ratios_overall_df = pd.concat(optimal_match_ratios_overall_df)
    
    # Removed dataset >>
    removed_dataset_df = optimal_match_ratios_overall_df[ optimal_match_ratios_overall_df["match_ratio_overall"] >= limit_match_score ]
    
    # Reducing dataset >>
    reduced_dataset_df = df.drop(removed_dataset_df["index_seq_1"])
    
    # Printing matched comparision >>
    if show_details:
        print(divider_line)
        print(divider_line)
        cprint("                    Matched sequences", "green")
        print(divider_line)
        print(divider_line)
        for index, row_removed_dataset in removed_dataset_df.iterrows():
            print_seqs_comparision( seq_1               = row_removed_dataset["seq_1"],
                                    seq_2               = row_removed_dataset["seq_2"],
                                    start_seq_1         = int(row_removed_dataset["start_seq_1"]),
                                    start_seq_2         = int(row_removed_dataset["start_seq_2"]),
                                    length_comparision  = int(row_removed_dataset["length_comparision"]),
                                    index_seq_1         = int(row_removed_dataset["index_seq_1"]),
                                    index_seq_2         = int(row_removed_dataset["index_seq_2"]) )
            print_details_comparision_seqs( match                           = row_removed_dataset["match"],
                                            mismatch                        = row_removed_dataset["mismatch"],
                                            length_comparision              = row_removed_dataset["length_comparision"],
                                            length_non_comparision          = row_removed_dataset["length_non_comparision"],
                                            match_ratio_comparision_length  = row_removed_dataset["match_ratio_comparision_length"],
                                            match_ratio_overall             = row_removed_dataset["match_ratio_overall"] )
    
    return reduced_dataset_df, removed_dataset_df, optimal_match_ratios_overall_df, df

################################## reduce_dataset_based_on_match_ratio #############################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/plots/sequence_comparision_mod #####################################################
####################################################################################################################################
