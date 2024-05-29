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


# Function to add significance brackets
def add_significance_bracket(ax, x1, x2, y, h, text, color='black'):
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.5, c=color)
    ax.text((x1 + x2) * .5, y + h, text, ha='center', va='bottom', color=color)


def earth_movers_distance(u, v):
    """
    If I equip the probability simplex with the ground metric of the L1 norm,
    it turns the EMD into a simple L1 norm. So that's how we'll do it.
    Parameters:
    - vec1: A list or numpy array representing the first distribution.
    - vec2: A list or numpy array representing the second distribution.
    Returns:
    - The Earth Mover's Distance between vec1 and vec2.
    """
    # compute L1 norm between u and v
    emd = np.linalg.norm(np.array(u) - np.array(v), ord=1)
    return emd


def start_box_plot(raw_distances, average_distances, pretty_name):
    # Creating the box plots
    plt.figure(figsize=(12, 8))

    # Improved color scheme
    colors = ['orange', 'teal']

    # Box plot positions
    n_time_points = len(raw_distances['Time'])
    positions_env_1_mix = [x * 2 for x in range(n_time_points)]  # Positions for env_1_mix box plots
    positions_env_2_mix = [x * 2 + 1 for x in range(n_time_points)]  # Positions for env_2_mix box plots

    # Create box plots
    bp_env_1_mix = plt.boxplot(raw_distances['env_1_mix'], positions=positions_env_1_mix, widths=0.6, patch_artist=True,
                           boxprops=dict(facecolor=colors[0]))
    bp_env_2_mix = plt.boxplot(raw_distances['env_2_mix'], positions=positions_env_2_mix, widths=0.6, patch_artist=True,
                           boxprops=dict(facecolor=colors[1]))

    # Add lines connecting averages
    plt.plot(positions_env_1_mix, average_distances['env_1_mix'], color=colors[0], marker='o', label='Average Env 1 vs '
                                                                                                     'mix')
    plt.plot(positions_env_2_mix, average_distances['env_2_mix'], color=colors[1], marker='o', label='Average Env 2 vs '
                                                                                                     'mix')

    # Adding titles and labels
    plt.title(f'Distribution of {pretty_name} Distances Over Time')
    plt.xlabel('Time Points')
    plt.ylabel('Distance')
    plt.xticks([r * 2 + 0.5 for r in range(n_time_points)],
               raw_distances['Time'])  # Set x-ticks to be between the pairs of box plots

    # Adding a legend
    legend_elements = [Patch(facecolor=colors[0], label='Env 1 vs mix'),
                       Patch(facecolor=colors[1], label='Env 2 vs mix')]
    plt.legend(handles=legend_elements, title='Sample Type Pairs')
    return positions_env_1_mix, positions_env_2_mix


def add_stats_to_plot(raw_distances, positions_env_1_mix, positions_env_2_mix):
    # Add significance brackets
    p_values = []

    for i, time_point in enumerate(raw_distances['Time']):
        dist_env_1_mix = raw_distances['env_1_mix'][i]
        dist_env_2_mix = raw_distances['env_2_mix'][i]

        # Perform t-test on the two independent samples
        #t_stat, p_val = ttest_ind(dist_env_1_mix, dist_env_2_mix)
        t_stat, p_val = mannwhitneyu(dist_env_1_mix, dist_env_2_mix)
        p_values.append(p_val)

    # Adjust y value for significance bracket based on your plot's max value
    y_max = np.max([np.max(raw_distances['env_1_mix']), np.max(raw_distances['env_2_mix'])])
    y_value_for_bracket = y_max * 1.1  # for example, 10% above the max value

    # Plot significance on the existing boxplot
    for i, p_val in enumerate(p_values):
        if p_val < 0.05:  # Choose your significance level here
            star_text = ('****' if p_val < 0.0001 else  # four stars for p < 0.0001
                         '***' if p_val < 0.001 else  # three stars for p < 0.001
                         '**' if p_val < 0.01 else  # two stars for p < 0.01
                         '*' if p_val < 0.05 else  # single star for p < 0.05
                         '')  # no star if p is not significant

            add_significance_bracket(plt.gca(), positions_env_1_mix[i], positions_env_2_mix[i], y_value_for_bracket,
                                     0.05 * y_max,
                                     star_text)
