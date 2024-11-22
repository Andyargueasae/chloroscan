import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
import matplotlib
import math
import os
from os import listdir
import statistics
import seaborn as sns
from pathlib import Path
import sys
from visualization_utils import *
from plotnine import ggplot, geom_point, aes, stat_smooth, facet_wrap, geom_boxplot, geom_bar, coord_flip, geom_violin
import plotnine as gg


cross_reference_tsv = snakemake.input[0]
figure_output_dir = Path(snakemake.output[0])
batch_name = snakemake.params[0]
refine_bins_dir = snakemake.params['refine_bins_dir']

# In order to avoid confusions, split batch name by "/".
# Because we want to have a valid single dirname.
batch_name = batch_name.split("/")[-1]

figure_output_dir.mkdir(parents=True, exist_ok=True)

time_stamp_file = ".snakemake_timestamp"

refine_bins_list = [i for i in listdir(refine_bins_dir) if i != time_stamp_file]
# Safety check if there were nothing.
if (os.stat(cross_reference_tsv).st_size == 0):
    print("Empty cross reference tsv, then nothing will be visualized.")
    sys.exit(0)
elif (len(refine_bins_list)==0):
    print("Empty bins directory after binning, process shut down.")
    sys.exit(0)

# Visualization rule will take the spreadsheet created from the last rule, and plot the statistical 
# graphics. Three graphics will be demonstrated:
#   1. contig-length distribution (a boxplot, or a histogram), to show how skewed the distribution is.
#   2. The contig length - marker number scatterplot, to show how markers distribute over this plastid assembly;
#   3. The taxonomy summarization of the assembly: how many taxa are there?
#   4. The Bin Purity plot: see how each bin's contig taxonomy homogenity?

summary_dataframe = pd.read_csv(cross_reference_tsv, sep="\t", index_col=0)

assembly_length_array = np.array(summary_dataframe['contig length'])
assembly_depth_array = np.array(summary_dataframe['contig depth'])
assembly_bin_array = list(summary_dataframe['Contig2Bin'])
contig_taxon = summary_dataframe['Taxon per Contig']
contig_2_bin = summary_dataframe['Contig2Bin']

# colors = ['black', 'red', 'green', 'blue', 'brown','cyan','grey','goldenrod','lime','violet','indigo','coral','olive','azure', 'light pink', 'dark blue']
# color_iter = iter(colors)

# output the violin plot.
summary_dataframe['batch depth per contig']=compute_batch_depth(assembly_depth_array)

# depth_violin_plot(summary_dataframe, assembly_bin_array)
binned_contigs_per_batch=summary_dataframe.loc[summary_dataframe["Contig2Bin"].notna()]
simplified_binname = MAG_name_simplify(binned_contigs_per_batch, batch_name)
binned_contigs_per_batch['Bin name simplified'] = simplified_binname

#store log10 contig depth.
log10_depth = np.log10(np.array(binned_contigs_per_batch['batch depth per contig']))
binned_contigs_per_batch['log10 depth per contig'] = log10_depth
violin_plot_by_bin=(
                # The plot firstly needs ggplot attribute, to select cols.
                ggplot(binned_contigs_per_batch, aes(x="Bin name simplified", y="log10 depth per contig", color="Bin name simplified"))
                + geom_violin()
                + geom_point()
                + gg.coord_flip()
                + gg.theme(figure_size=(15, 13), 
                           axis_title_x=gg.element_text(size=14), 
                           axis_title_y=gg.element_text(size=14),
                           axis_text_x=gg.element_text(size=12),
                           axis_text_y=gg.element_text(size=12)) # Theme is responsible for fontsize, legend, figsize...
).draw(show=False)

violin_plot_by_bin.savefig(f"{figure_output_dir}/LogDepth_Violin.png", dpi=300)

#2. sequence depth, gc contents and marker gene.

plt.rcParams["figure.autolayout"] = True
plt.rcParams["figure.figsize"] = [25, 36]
plt.rcParams.update({'font.size': 24})

fig1, axsc = plt.subplots(nrows=2, ncols=1)
sns.set_theme(font_scale=3)
# store count of marker genes.
markers_per_contig = binned_contigs_per_batch['markers on the contig']
marker_count = marker_counter(markers_per_contig)
binned_contigs_per_batch['marker count per contig'] = marker_count

Scatter_GC_log_depth=(
                        ggplot(binned_contigs_per_batch, aes(x="GC contents", y="log10 depth per contig", color="Bin name simplified", size = "marker count per contig"))
                            + geom_point(alpha=0.5)
                            + gg.labs(color="Bins", size="marker count per contig")
                            + gg.theme(legend_title=gg.element_text(size=12, weight='bold'), 
                                    legend_text=gg.element_text(size=8))
                            + gg.ggtitle("dot-scaled GC x LogDepth scatterplot")
                    ).draw(show=False) # Add size and colors here.
                            # Add point scale. 

Scatter_GC_log_depth.savefig(f"{figure_output_dir}/Scatter_GCLogDepth.png", dpi=300)

# The second plot is the bar chart of the completeness and purity.
bins = [i for i in list(set(binned_contigs_per_batch['Contig2Bin']))]
# works at the x-label.
bins_quality_dict = dict()
# works as the legend.
bins_quality_dict['Bin Completeness'] = []
bins_quality_dict['Bin Purity'] = []
bin_xticks = []

# Need to consider scMAG problem: we need a cool shit to process their name. 
for i in bins:
    bin_base_names = i.split("_")
    prefix = "bin"
    for j in range(len(bin_base_names)):
        if "C" in bin_base_names[j]:
            try: 
                individual_comp = int(bin_base_names[j].replace("C", ""))
                bins_quality_dict['Bin Completeness'].append(individual_comp)
                identifier = bin_base_names[j-1]
            except ValueError:
                continue
        if "P" in bin_base_names[j]:
            individual_pur = int(bin_base_names[j].replace("P", ""))
            bins_quality_dict['Bin Purity'].append(individual_pur)
        else:
            continue
    bin_name = prefix + "." + f"{identifier}_C{individual_comp}_P{individual_pur}"
    bin_xticks.append(bin_name)

width = 0.2
x=np.arange(len(bins))
multiplier=0

for attribute, measurement in bins_quality_dict.items():
    # Add scatterplot of completeness against purity.

    offset = width * multiplier
    bin_record = axsc[1].bar(x+offset, measurement, width, label=attribute)
    axsc[1].bar_label(bin_record, padding=3)
    multiplier += 1
axsc[0].scatter(bins_quality_dict["Bin Completeness"], bins_quality_dict["Bin Purity"])
axsc[0].set_xlabel("Bin Completeness (%)")
axsc[0].set_ylabel("Bin Purity (%)")

axsc[1].set_ylabel("percentage")
axsc[1].set_xlabel("bin name")
axsc[1].set_xticks(x + width, bin_xticks, rotation=30)
axsc[1].legend(loc='lower left', ncols=len(bins))
axsc[1].set_ylim(50, 100)

plt.savefig("{}/{}_bins_scatter_bar_chart.png".format(figure_output_dir, batch_name), dpi="figure")

#3. taxonomy pie chart.
bin_set = set(contig_2_bin)
bin_set
for i in bin_set:
    if not pd.isna(i):
        contig_composition = summary_dataframe.loc[summary_dataframe['Contig2Bin'] == i]
        
        # prepare the graph.
        MAGfig, MAGax = plt.subplots(figsize=(18,18))
        sns.set_theme(font_scale=2)
        plt.style.use('_mpl-gallery-nogrid')
        contig_composition_taxon_prop_bylength, taxon_labels = count_and_prop(contig_composition)[0], count_and_prop(contig_composition)[1]
        colors = plt.get_cmap('jet')(np.linspace(0.2, 0.7, len(set(contig_composition_taxon_prop_bylength))+10))
        MAGax.pie(contig_composition_taxon_prop_bylength, startangle=90, colors=colors, radius=2, wedgeprops={"linewidth": 1, "edgecolor": "white"}, center=(4,4), frame=True, autopct='%1.2f%%')
        plt.axis('equal')
        plt.legend(taxon_labels, loc="best", fontsize="29")
        plt.title("{} bin {} taxon composition".format(batch_name, i), fontsize=30)
        plt.show()
        # Hide xticks and yticks.
        plt.xticks([])
        plt.yticks([])
        plt.savefig("{}/{}_{}_taxonomy_composition.png".format(figure_output_dir, batch_name, i), dpi="figure")
    else:
        continue