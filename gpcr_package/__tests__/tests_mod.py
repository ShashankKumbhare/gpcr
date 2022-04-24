
####################################################################################################################################
################################## gpcr_package/__tests__/tests_mod ################################################################
####################################################################################################################################

################################################################
from ..__info__ import *
################################################################


####################################################################################################################################
################################## test_dependencies ###############################################################################
# >>
def test_dependencies(dependencies_mod):
    print("")
    print(">>>>>>>>>>>>>>>>> Testing dependencies...")
    print("`" + package_name + "` depends on the following packages:")
    print( [ dep for dep in dir(dependencies_mod) if not(dep.startswith('__', 0)) ] )
    print("<<<<<<<<<<<<<<<<< Dependencies test successful!")
    print("")
# <<
################################## test_dependencies ###############################################################################
####################################################################################################################################


####################################################################################################################################
################################## test_constants ##################################################################################
# >>
def test_constants(constants_mod):
    print("")
    print(">>>>>>>>>>>>>>>>> Testing constants...")
    print("`" + package_name + "` uses the following constants:")
    print( [ constant for constant in dir(constants_mod) if not(constant.startswith('__', 0)) ] )
    print("<<<<<<<<<<<<<<<<< Constants test successful!")
    print("")
# <<
################################## test_constants ##################################################################################
####################################################################################################################################


####################################################################################################################################
################################## test_auxils #####################################################################################
# >>
def test_auxils(auxil_list):
    print("")
    print(">>>>>>>>>>>>>>>>> Testing auxils...")
    print("`" + package_name + "` uses the following auxiliary functions:")
    print( [ auxil.__name__ for auxil in auxil_list ] )
    print("<<<<<<<<<<<<<<<<< Auxils test successful!")
    print("")
# <<
################################## test_auxils #####################################################################################
####################################################################################################################################


####################################################################################################################################
################################## test_datasets ###################################################################################
# >>
def test_datasets(datasets_mod):
    print("")
    print(">>>>>>>>>>>>>>>>> Testing datasets...")
    print("`" + package_name + "` uses the following datasets:")
    print( [ dataset for dataset in dir(datasets_mod) if not(dataset.startswith('__', 0) or dataset.startswith('_', 0)) ] )
    print("<<<<<<<<<<<<<<<<< Datasets test successful!")
    print("")
# <<
################################## test_datasets ###################################################################################
####################################################################################################################################


####################################################################################################################################
################################## test_plots ######################################################################################
# >>
def test_plots(plot_list):
    print("")
    print(">>>>>>>>>>>>>>>>> Testing plots...")
    print("`" + package_name + "` provides the following plots:")
    # print( [ plot for plot in plot_list ] )
    for plot in plot_list:
        print(plot)
    print("<<<<<<<<<<<<<<<<< Plots test successful!")
    print("")
# <<
################################## test_plots ######################################################################################
####################################################################################################################################


####################################################################################################################################
################################## test_sequence_analysis ##########################################################################
# >>
def test_sequence_analysis(sequence_analysis_tools):
    print("")
    print(">>>>>>>>>>>>>>>>> Testing sequence analysis tools...")
    print("`" + package_name + "` uses the following sequence analysis tools:")
    for sequence_analysis_tool in sequence_analysis_tools:
        print(sequence_analysis_tool)
    print("<<<<<<<<<<<<<<<<< Datasets test successful!")
    print("")
# <<
################################## test_sequence_analysis ##########################################################################
####################################################################################################################################


####################################################################################################################################
################################## gpcr_package/__tests__/tests_mod ################################################################
####################################################################################################################################