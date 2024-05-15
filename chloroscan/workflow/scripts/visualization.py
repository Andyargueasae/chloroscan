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
    fig, ax1 = plt.subplots(figsize=(35, 12))
    matplotlib.rcParams.update({'font.size': 48})
    colors = ['yellow', 'red', 'green', 'blue', 'brown','cyan','grey','goldenrod','lime','violet','indigo','coral','olive','azure', 'light pink', 'dark blue']
    color_iter = iter(colors)
    sample_binned_df['Contig2Bin'] = pd.Categorical(sample_binned_df['Contig2Bin'], [i for i in set(sample_binned_df['Contig2Bin'])])
    log10_depth_sample = np.log10(np.array(sample_binned_df['batch depth per contig']))
    sample_binned_df['log10 depth per contig'] = log10_depth_sample
    groups = sample_binned_df.groupby("Contig2Bin")
    index = 0
    for group_name, group in groups:
        sns.violinplot(x="log10 depth per contig",y="Contig2Bin", data=group, orient="h", ax=ax1, color = colors[index])
        sns.stripplot(x="log10 depth per contig", y="Contig2Bin", data=group, orient="h", jitter=True, zorder=1, size=6, color="purple")
        index += 1
    plt.yticks(fontsize=14)
    plt.xticks(fontsize=18)
    plt.xlabel("log10 depth per contig", fontsize=20)
    plt.title("Violin Plot of batch depth for each bin")
    violin_plot_path = figure_output_dir/f"{batch_name}_Contig_Depth_Violin_Plot.png"
    plt.savefig(str(violin_plot_path), dpi="figure")
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



plt.rcParams["figure.autolayout"] = True
plt.rcParams["figure.figsize"] = [25, 36]
plt.rcParams.update({'font.size': 24})

fig1, axsc = plt.subplots(nrows=2, ncols=1)
sns.set(font_scale=3)
# store count of marker genes.
markers_per_contig = binned_contigs_per_batch['markers on the contig']
marker_count = marker_counter(markers_per_contig)
binned_contigs_per_batch['marker count per contig'] = marker_count

#store log10 contig depth.
log10_depth = np.log10(np.array(binned_contigs_per_batch['batch depth per contig']))
binned_contigs_per_batch['log10 depth per contig'] = log10_depth
groups = binned_contigs_per_batch.groupby("Contig2Bin")
for bin_name, bin_group in groups:
    group_label = bin_name.replace("_", " ")
    print(group_label)
    sns.scatterplot(data=bin_group, x="GC contents", y="log10 depth per contig", alpha=0.7, s=(np.array(5*5*bin_group['marker count per contig'])+120), marker="o", ax = axsc[0], label = group_label )
plt.legend(bbox_to_anchor=(1.02, 1), loc='lower left', borderaxespad=0)
axsc[0].set_title("GC vs. sequence depth, dot size scaled by marker count", fontsize=40)
axsc[0].set_xlabel("GC content", fontsize=45)
axsc[0].set_ylabel(r"$Sequence Depth (log_{10})$", fontsize=45)

# The second plot is the bar chart of the completeness and purity.
bins = [i for i in list(set(binned_contigs_per_batch['Contig2Bin']))]
# works at the x-label.
bins_quality_dict = dict()
# works as the legend.
bins_quality_dict['Bin Completeness'] = []
bins_quality_dict['Bin Purity'] = []
bin_xticks = []
for i in bins:
    bin_base_names = i.split("_")
    prefix = "bin"
    for j in range(len(bin_base_names)):
        if "C" in bin_base_names[j]:
            individual_comp = int(bin_base_names[j].replace("C", ""))
            bins_quality_dict['Bin Completeness'].append(individual_comp)
            identifier = bin_base_names[j-1]
        if "P" in bin_base_names[j]:
            individual_pur = int(bin_base_names[j].replace("P", ""))
            bins_quality_dict['Bin Purity'].append(individual_pur)
        else:
            continue
    bin_name = prefix + "." + f"{identifier}_C{individual_comp}_P{individual_pur}"
    # breakpoint()
    bin_xticks.append(bin_name)

width = 0.2
x=np.arange(len(bins))
multiplier=0

for attribute, measurement in bins_quality_dict.items():
    offset = width * multiplier
    bin_record = axsc[1].bar(x+offset, measurement, width, label=attribute)
    axsc[1].bar_label(bin_record, padding=3)
    multiplier += 1

axsc[1].set_ylabel("percentage")
axsc[1].set_xlabel("bin name")
axsc[1].set_xticks(x + width, bin_xticks, rotation=30)
axsc[1].legend(loc='lower left', ncols=len(bins))
axsc[1].set_ylim(50, 100)

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
        plt.legend(taxon_labels, loc="best", fontsize="29")
        plt.title("{} bin {} taxon composition".format(batch_name, i), fontsize=30)
        plt.show()
        # Hide xticks and yticks.
        plt.xticks([])
        plt.yticks([])
        plt.savefig("{}/{}_{}_taxonomy_composition.png".format(figure_output_dir, batch_name, i), dpi="figure")
    else:
        continue

del summary_dataframe
# Visualization is done, delete dataframe.