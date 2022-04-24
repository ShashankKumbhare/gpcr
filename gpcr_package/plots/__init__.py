
####################################################################################################################################
################################## gpcr_package/plots/___init__ ####################################################################
####################################################################################################################################

from .hist_plot_mod     import hist_quick, hist_gpcr_Ctail_lens, hist_proteins_residue_loc, hist_proteins_residue_relative_count
from .sankey_plot_mod   import sankey_gpcr_classes
from .grid_plot_mod     import grid_hist_proteins_residue_loc, grid_scatter_proteins_residue_count_vs_len, grid_hist_proteins_residue_relative_count
from .box_plot_mod      import box_proteins_residue_loc
from .scatter_plot_mod  import scatter_proteins_residue_count_vs_len
from .line_plot_mod     import line_compare_points, line_hydrophobicity_in_sequences, line_charges_in_sequences, line_residue_in_ctail, line_motifs_in_ctail
from .heatmap_plot_mod  import heatmap_1, heatmap_2
from .bar_plot_mod      import bar_gpcr_classes

plot_list = [ # histogram plots >>
              hist_quick, hist_gpcr_Ctail_lens, hist_proteins_residue_loc, hist_proteins_residue_relative_count,
              # sankey plots >>
              sankey_gpcr_classes,
              # grid plots >>
              grid_hist_proteins_residue_loc, grid_scatter_proteins_residue_count_vs_len, grid_hist_proteins_residue_relative_count,
              # box plots >>
              box_proteins_residue_loc,
              # scatter plots >>
              scatter_proteins_residue_count_vs_len,
              # line plots >>
              line_compare_points, line_hydrophobicity_in_sequences, line_charges_in_sequences, line_residue_in_ctail, line_motifs_in_ctail,
              # heatmap plots >>
              heatmap_1, heatmap_2,
              # bar plots >>
              bar_gpcr_classes
             ]

####################################################################################################################################
################################## gpcr_package/plots/___init__ ####################################################################
####################################################################################################################################
