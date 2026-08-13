import json
import matplotlib as mpl
import os


metrics_colors_dictionary = {"ofml"        : "viridis_r",
                              "ofml_score"  : "#6A33E0",
                              "omega_trunc" : "#FA5E32",
                              "omega_synon" : "#89E4A2",
                              "omega_miss"  : "#FABE4A",
                              "o3d_score"   : "#6DBDCC",
                              "o3d_cluster" : "skyblue",
                              "o3d_prob"    : "darkgray",
                              "frameshift"  : "#E4ACF4",
                              "inframe"     : "C5",
                              "hv_lines"    : "lightgray", # horizontal and vertical lines,
                              "hv_lines_needle" : "gray",
                              "needle_obs"  : "#003366",
                              "omega_miss_tert" : "#f5840c",
                              "omega_synon_tert": "#378c12",
                              "nonsense" : "#FA5E32",
                              "synonymous" : "#89E4A2",
                              "missense"  : "#FABE4A",
                              #"nonsense"    : "#FB8E6F",  
                              # "synonymous"  : "#ACECBD", 
                              #"missense"    : "#FBD180", 
                              "indel"       : "#ECC4F7", 
                              "splicing"    : "#A1C5DF",
                              }


####
# plotting configurations
####

plots_general_config = {

                        # fonsizes
                        "ylabel_fontsize": 6,
                        "xlabel_fontsize": 6,
                        "xylabel_fontsize": 6,
                        "title_fontsize": 7,
                        "xyticks_fontsize": 5,
                        "xticks_fontsize": 5,
                        "yticks_fontsize": 5,
                        "legend_fontsize": 5,
                        "annots_fontsize": 5,

                        "dot_size_scplot": 15,
                        "dot_size_coeffplot": 5,
                        "dot_sizebelow_coeffplot": 40,
                        "dot_color_coeffplot": "#D3D3D3",
                        "dot_colorabove_coeffplot": "#D62728",
                        "dot_colorbelow_coeffplot": "#f29c9e",
                        "dot_edgethres_coeffplot": 0.2,
                        "dot_edgewidth_coeffplot": 0.5
                        }


mpl.rcParams.update({
   'font.family': 'Arial',            # Enforce Arial
   'pdf.fonttype': 42,                # TrueType for PDF
   'ps.fonttype': 42,                 # TrueType for PS/EPS
   'svg.fonttype': 'none',            # Keep text as editable text
})


mpl.rcParams.update({
   'axes.titlesize'    : plots_general_config["title_fontsize"],       # Title font size
   'axes.labelsize'    : plots_general_config["xylabel_fontsize"],     # X and Y axis labels
   'xtick.labelsize'   : plots_general_config["xyticks_fontsize"],     # X tick labels
   'ytick.labelsize'   : plots_general_config["xyticks_fontsize"],     # Y tick labels
   'legend.fontsize'   : plots_general_config["legend_fontsize"],      # Legend text
   'figure.titlesize'  : plots_general_config["title_fontsize"],       # Figure suptitle (if used)
})