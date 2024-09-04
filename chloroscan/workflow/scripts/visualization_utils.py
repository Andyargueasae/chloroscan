import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import numpy as np

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

def MAG_name_simplify(DF, BATCH_NAME):
    bins = list(DF.loc[:, "Contig2Bin"])
    bins_simplified = []
    for bin in bins:
        bin_feature = "_".join([BATCH_NAME, *bin.split("_")[-3:]])
        bins_simplified.append(bin_feature)
    return bins_simplified

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
            depth_batch = float(i)
            batch_depth.append(depth_batch)
    return batch_depth
