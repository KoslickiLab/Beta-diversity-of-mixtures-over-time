# This will be a command line script implementing the approach in Analysis.py
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage
import numpy as np
from matplotlib.patches import Patch
from scipy.stats import shapiro
from scipy.stats import ttest_ind, mannwhitneyu
import argparse
from src import add_significance_bracket, earth_movers_distance, start_box_plot, add_stats_to_plot

# Make the CLI. Want the input to be two CSV files: count_data and metadata. Output will be the pairwise distances.
# Optionally, can generate plots and statistics too.
parser = argparse.ArgumentParser(description='Calculate pairwise distances between samples and mixtures.')
# Want the count data argument to also be invoked as -c
parser.add_argument('-c', '--count_data', type=str, help='Path to the count data CSV file.', required=True)
# Want the metadata argument to also be invoked as -m
parser.add_argument('-m', '--metadata', type=str, help='Path to the metadata CSV file.', required=True)
# Want the output argument to also be invoked as -o
parser.add_argument('-o', '--output_prefix', type=str, help='Prefix for the output (eg. -o /a/b/c/prefix will '
                                                            'output /a/b/c/prefix_pariwise_distances.csv',
                    required=True)
# Want the plot argument to also be invoked as -p. Should be a flag, no arguments needed.
parser.add_argument('-p', '--plot', action='store_true', help='Generate plots of the pairwise distances.')

# parse the args
args = parser.parse_args()
count_data_file = args.count_data
metadata_file = args.metadata
output_file_prefix = args.output_prefix
plot_flag = args.plot
# sanity check all the arguments, check if files exist and I have read access and write access to the output path
# Simultanously, load the data
try:
    df = pd.read_csv(count_data_file, sep=',', header=0, index_col=0)
except FileNotFoundError:
    raise FileNotFoundError("Count data file not found or do not have read access to it.")
try:
    meta = pd.read_csv(metadata_file, sep=',', header=0, index_col=0)
except FileNotFoundError:
    raise FileNotFoundError("Metadata file not found or do not have read access to it.")
try:
    with open(output_file_prefix, 'w'):
        pass
    # remove the file, we just wanted to check if we have write access to the directory
    os.remove(output_file_prefix)
except FileNotFoundError:
    raise FileNotFoundError("Output file path not found or do not have write access to it.")

if 'Time' not in meta.columns:
    raise ValueError("Metadata file must have a 'Time' column.")
# check if the time column is numeric (specifically, integers)
if not meta['Time'].apply(lambda x: str(x).isdigit()).all():
    raise ValueError("Time column must be numeric. Eg. 1, 2, 5, 10, etc.")

# define the metrics
metrics = ['euclidean', 'cosine', 'jensenshannon', earth_movers_distance, 'jaccard']
metrics_pretty = {'euclidean': 'Euclidean', 'cosine': 'Cosine', 'jensenshannon': 'Jensen-Shannon',
                  'earth_movers_distance': 'Earth Mover\'s Distance', 'jaccard': 'Jaccard'}


if plot_flag:
    for metric in metrics:
        # Compute pairwise distances
        # in the case of the jaccard, I need to set everything to binary in the df, but not overwrite the original df
        if metric == 'jaccard':
            binary_df = df.copy()
            binary_df[binary_df > 0] = 1
            dists = pdist(binary_df, metric=metric)
        else:
            dists = pdist(df, metric=metric)
        square_dists = squareform(dists)

        # Perform hierarchical clustering using the distance matrix
        row_linkage = linkage(dists, method='average')

        # Create a color map for each metadata category and its corresponding legend
        category_colors = {}
        legends = []

        for category in meta.columns:
            # Create a color palette, you can adjust this as needed
            unique_categories = meta[category].unique()
            palette = sns.color_palette("hls", len(unique_categories))
            lut = dict(zip(unique_categories, palette))

            # Map the metadata to colors
            category_colors[category] = meta[category].map(lut)

            # Create legend handles for each category
            handles = [mpatches.Patch(color=lut[name], label=name) for name in lut]
            legends.append(handles)

        # Combine the color mappings into a DataFrame
        row_colors = pd.DataFrame(category_colors)

        # Create a DataFrame from the square distance matrix
        distance_df = pd.DataFrame(square_dists, index=df.index, columns=df.index)

        # Generate clustermap
        g = sns.clustermap(distance_df, row_linkage=row_linkage, col_linkage=row_linkage, figsize=(20, 16),
                           cmap="viridis", row_colors=row_colors, col_colors=row_colors)

        # The Axes object of the heatmap
        ax_heatmap = g.ax_heatmap

        # Add legends for each category
        for i, legend_handles in enumerate(legends):
            legend = ax_heatmap.legend(handles=legend_handles,
                                       loc='upper left',
                                       bbox_to_anchor=(-0.3 + 0.1 * i, 1.4),
                                       borderaxespad=0.)
            ax_heatmap.add_artist(legend)
        plt.title(f'Clustering of {metrics_pretty[metric] if type(metric) is str else metrics_pretty[metric.__name__]} '
                  f'distances')
        heatmap_suffix = f"_heatmap_dendrogram" \
                         f"_{metrics_pretty[metric] if type(metric) is str else metrics_pretty[metric.__name__]}.png"
        g.savefig(os.path.abspath(output_file_prefix + heatmap_suffix), dpi=300)

# Box plots of distances over time decorated with statistics

# Filter and reorder meta
meta = meta.loc[df.index]
meta = meta.reindex(df.index)

#metrics = ['euclidean', 'cosine', 'jensenshannon', earth_movers_distance]
#metrics_pretty = {'euclidean': 'Euclidean', 'cosine': 'Cosine', 'jensenshannon': 'Jensen-Shannon',
# 'earth_movers_distance': 'Earth Mover\'s Distance'}

raw_distances_by_metric = {}
for metric in metrics:
    # Initialize data structures to store the raw distances and averages
    raw_distances = {'Time': [], 'env_1_mix': [], 'env_2_mix': []}
    average_distances = {'Time': [], 'env_1_mix': [], 'env_2_mix': []}
    pretty_name = metrics_pretty[metric] if type(metric) is str else metrics_pretty[metric.__name__]

    # Calculate distances for each time point
    for time_point in np.sort(meta['Time'].unique()):
        time_samples = meta[meta['Time'] == time_point]
        time_point_samples = df.loc[time_samples.index]
        # dist_matrix = pdist(time_point_samples, metric='euclidean')
        if metric == 'jaccard':
            # make a copy of the df
            binary_df = time_point_samples.copy()
            binary_df[binary_df > 0] = 1
            dist_matrix = pdist(binary_df, metric=metric)
        else:
            dist_matrix = pdist(time_point_samples, metric=metric)
        dist_matrix_square = squareform(dist_matrix)

        # Identifying samples
        try:
            samples_A = time_point_samples[time_samples['Environment'] == 'env1']
            samples_B = time_point_samples[time_samples['Environment'] == 'env2']
            samples_AxB = time_point_samples[time_samples['Environment'] == 'mix']
        except KeyError:
            raise ValueError("Metadata file must have an 'Environment' column with values 'env1', 'env2', and 'mix'.")
        if len(samples_A) == 0 or len(samples_B) == 0 or len(samples_AxB) == 0:
            raise ValueError("Metadata file must have an 'Environment' column with values 'env1', 'env2', and 'mix' "
                             "for each time point.")

        # Convert index values to numerical indices
        numerical_indices_A = time_point_samples.index.get_indexer(samples_A.index)
        numerical_indices_B = time_point_samples.index.get_indexer(samples_B.index)
        numerical_indices_AxB = time_point_samples.index.get_indexer(samples_AxB.index)

        # Store raw distances and calculate averages
        raw_distances['Time'].append(time_point)
        dist_env_1_mix = dist_matrix_square[np.ix_(numerical_indices_A, numerical_indices_AxB)].flatten()
        dist_env_2_mix = dist_matrix_square[np.ix_(numerical_indices_B, numerical_indices_AxB)].flatten()
        raw_distances['env_1_mix'].append(dist_env_1_mix)
        raw_distances['env_2_mix'].append(dist_env_2_mix)
        average_distances['Time'].append(time_point)
        average_distances['env_1_mix'].append(np.mean(dist_env_1_mix))
        average_distances['env_2_mix'].append(np.mean(dist_env_2_mix))

    raw_distances_by_metric[pretty_name] = raw_distances

    positions_env_1_mix, positions_env_2_mix = start_box_plot(raw_distances, average_distances, pretty_name)
    add_stats_to_plot(raw_distances, positions_env_1_mix, positions_env_2_mix)
    # save the figure, if it's asked for
    if plot_flag:
        plt.savefig(os.path.abspath(output_file_prefix + f"_boxplot_{pretty_name}.png"), dpi=300)


    # Export to CSV the raw distances and the average distances
    #raw_distances_df = pd.DataFrame(raw_distances)
    # redo the raw distances so that it looks like
    #	Time	Pair	Values
    #1	1	env_1_mix	0.49478978
    #2	1	env_1_mix	0.46863569
    #3	1	env_1_mix	0.5195925
    #4	1	env_1_mix	0.47284923
    #5	1	env_1_mix	0.48184655
    #6	1	env_1_mix	0.4928068
    #7	1	env_1_mix	0.49931512
    #8	1	env_1_mix	0.47241578
    #9	1	env_1_mix	0.517094
    #10	1	env_1_mix	0.4779161
    #11	1	env_1_mix	0.48163473
    #12	1	env_1_mix	0.49896055
    raw_distances_df = pd.DataFrame()
    # need three columns: Time, Pair, Values
    for index, time_point in enumerate(raw_distances['Time']):
        env_1_mix = raw_distances['env_1_mix'][index]
        env_2_mix = raw_distances['env_2_mix'][index]
        for value in env_1_mix:
            raw_distances_df = raw_distances_df._append({'Time': time_point, 'Environment': 'env_1_mix', 'Values': value},
                                                        ignore_index=True)
        for value in env_2_mix:
            raw_distances_df = raw_distances_df._append({'Time': time_point, 'Environment': 'env_2_mix', 'Values': value},
                                                        ignore_index=True)

    average_distances_df = pd.DataFrame(average_distances)

    # file names
    distances_file = os.path.abspath(output_file_prefix + f"_raw_distances_{pretty_name}.csv")
    ave_distances_file = os.path.join(output_file_prefix + f"_average_distances_{pretty_name}.csv")

    # export
    raw_distances_df.to_csv(distances_file)
    average_distances_df.to_csv(ave_distances_file)
