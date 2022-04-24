
####################################################################################################################################
################################## gpcr_package/__auxil__/auxil_mod ################################################################
####################################################################################################################################

################################################################
from ..__dependencies__ import it, SeqIO, pd, np, pickle
################################################################

####################################################################################################################################
################################## auxil_func ######################################################################################
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
################################## auxil_func ######################################################################################
####################################################################################################################################


####################################################################################################################################
################################## auxil_func ######################################################################################
# >>
def save_variables_to_file(variables, file_pkl):
    with open(file_pkl, 'wb') as f:
        pickle.dump(variables, f)
# <<
################################## auxil_func ######################################################################################
####################################################################################################################################


####################################################################################################################################
################################## auxil_func ######################################################################################
# >>
def load_variables_from_file(file_pkl):
    with open(file_pkl, 'rb') as f:
        variables = pickle.load(f)
    return variables
# <<
################################## auxil_func ######################################################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/__auxil__/auxil_mod ################################################################
####################################################################################################################################
