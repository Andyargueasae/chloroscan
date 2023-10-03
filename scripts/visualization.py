import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
import math
import os
import statistics
import seaborn as sns

cross_reference_tsv = snakemake.input[0]
figure_output_dir = snakemake.output[0]
batch_name = snakemake.params[0]

os.mkdir(figure_output_dir)
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


colors = ['black', 'red', 'green', 'blue', 'brown','cyan','grey','goldenrod','lime','violet','indigo','coral','olive','azure', 'light pink', 'dark blue']
color_iter = iter(colors)
def compute_batch_depth(depth_Array):
    batch_depth = []
    for i in depth_Array:
        if "\t" in str(i):
            depth_list = [float(j) for j in i.split("\t")]
            batch_depth_per_contig = sum(depth_list)
            batch_depth.append(batch_depth_per_contig)
        else:
            depth_batch = float(i)
            batch_depth.append(depth_batch)
    return batch_depth


def count_and_prop(individual_bin_df):
    bin_comp_dict = dict()
    length_array_per_bin = list(individual_bin_df['contig length'])
    taxon_array_per_bin = list(individual_bin_df['Taxon per Contig'])

    # print(taxon_set)
    for i in range(0, len(taxon_array_per_bin)):
        if taxon_array_per_bin[i] in bin_comp_dict.keys():
            bin_comp_dict[taxon_array_per_bin[i]] += int(length_array_per_bin[i])
        else:
            bin_comp_dict[taxon_array_per_bin[i]] = int(length_array_per_bin[i]) # initialize dict value for new key.

    prop_list = bin_comp_dict.values()
    label_list = bin_comp_dict.keys()
    return prop_list, label_list 

def Seaborn_Violin_plot(sample_binned_df):
    fig, ax1 = plt.subplots(figsize=(32, 18))
    sns.set(font_scale=3)
    sample_binned_df['Contig2Bin'] = pd.Categorical(sample_binned_df['Contig2Bin'], [i for i in set(sample_binned_df['Contig2Bin'])])
    groups = sample_binned_df.groupby(['Contig2Bin'])
    sns.violinplot(x="batch depth per contig",y="Contig2Bin", hue="Contig2Bin", hue_order=sorted(list(set(sample_binned_df['Contig2Bin']))), data=sample_binned_df,ax=ax1)
    for artist in ax1.lines:
        artist.set_zorder(10)
    for artist in ax1.findobj(PathCollection):
        artist.set_zorder(11)
    sns.stripplot(x="batch depth per contig",y = "Contig2Bin", data=sample_binned_df, orient="h", jitter=True, size=6, color="purple", ax=ax1)
    # sns.move_legend(ax1, "upper left", bbox_to_anchor=(1, 1))
    plt.legend(loc="best")
    plt.savefig("{}/{}_Contig_Depth_Violin_Plot.png".format(figure_output_dir, batch_name), dpi="figure")
    return
    
# output the violin plot.
summary_dataframe['batch depth per contig']=compute_batch_depth(assembly_depth_array)
# depth_violin_plot(summary_dataframe, assembly_bin_array)
binned_contigs_per_batch=summary_dataframe.loc[summary_dataframe["Contig2Bin"].notna()]
Seaborn_Violin_plot(binned_contigs_per_batch)

#2. sequence depth, gc contents and marker gene.
def marker_counter(marker_array):
    marker_count_array = []
    for i in marker_array:
        # print(i)
        if pd.isna(i):
            marker_count_array.append(0)
        else:
            marker_list = [genes for genes in i.split(",") if genes != ""]
            # print(marker_list)
            marker_count_array.append(len(marker_list))
    return marker_count_array

def compute_batch_depth(depth_Array):
    batch_depth = []
    for i in depth_Array:
        str_i = str(i)
        if "\t" in str_i:
            depth_list = [float(j) for j in i.split("\t")]
            batch_depth_per_contig = sum(depth_list)
            batch_depth.append(batch_depth_per_contig)
        else:
            depth_batch = i
            batch_depth.append(depth_batch)
    return batch_depth

fig, axsc = plt.subplots(figsize=(32,18))
sns.set(font_scale=3)
# store count of marker genes.
markers_per_contig = binned_contigs_per_batch['markers on the contig']
marker_count = marker_counter(markers_per_contig)
binned_contigs_per_batch['marker count per contig'] = marker_count

#store log10 contig depth.
log10_depth = np.log10(np.array(binned_contigs_per_batch['batch depth per contig']))
binned_contigs_per_batch['log10 depth per contig'] = log10_depth

sns.scatterplot(data=binned_contigs_per_batch, x="GC contents", y="log10 depth per contig", hue="Contig2Bin", hue_order=sorted(list(set(binned_contigs_per_batch['Contig2Bin']))), alpha=0.7, s=(np.array(5*5*binned_contigs_per_batch['marker count per contig'])+80), marker="o")
plt.legend(loc="best")
# sns.move_legend(axsc, "upper left", bbox_to_anchor=(1, 1))
plt.title("GC vs. sequence depth, dot size scaled by marker count")

plt.savefig("{}/{}_binned_contigs_scatter.png".format(figure_output_dir, batch_name), dpi="figure")

#3. taxonomy pie chart.
bin_set = set(contig_2_bin)
bin_set
for i in bin_set:
#     print(i)
    if not pd.isna(i):
        # print(i)
        contig_composition = summary_dataframe.loc[summary_dataframe['Contig2Bin'] == i]
        # print(contig_composition)
        
        # prepare the graph.
        MAGfig, MAGax = plt.subplots(figsize=(18,18))
        sns.set(font_scale=2)
        plt.style.use('_mpl-gallery-nogrid')
        contig_composition_taxon_prop_bylength, taxon_labels = count_and_prop(contig_composition)[0], count_and_prop(contig_composition)[1]
        colors = plt.get_cmap('jet')(np.linspace(0.2, 0.7, len(set(contig_composition_taxon_prop_bylength))+10))
        MAGax.pie(contig_composition_taxon_prop_bylength, startangle=90, colors=colors, radius=2, wedgeprops={"linewidth": 1, "edgecolor": "white"}, center=(4,4), frame=True, autopct='%1.2f%%')
        plt.axis('equal')
        plt.legend(taxon_labels, loc="best")
        plt.title("{} bin {} taxon composition".format(batch_name, i))
        plt.show()
        plt.savefig("{}/{}_{}_taxonomy_composition.png".format(figure_output_dir, batch_name, i), dpi="figure")
    else:
        continue

del summary_dataframe
# Visualization is done, delete dataframe.