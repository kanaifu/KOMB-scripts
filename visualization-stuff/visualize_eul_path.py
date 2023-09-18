import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

alien_ids = [0, 2, 376, 377, 378, 380, 60, 600, 630, 632, 637, 7, 773]
alien_ids.sort()


def extract(name):
    file = open(name, 'r')
    lines = file.readlines()
    lst = []
    ct = 0
    under = []
    for line in lines:
        if line[0] == 'A' or line[0] == 'T' or line[0] == 'C' or line[0] == 'G':
            if len(line) - 1 >= 80:
                lst.append(len(line) - 1)
                if len(line) - 1 <= 300:
                    under.append(len(line) - 1)
            else:
                ct += 1
    file.close()
    return lst, ct, under


def visualize(data, labels, name):
    # Determine the number of tuples and biomes
    num_tuples = len(data)
    num_tools = len(data[0])

    # Define a color palette for each row
    row_colors = ['skyblue', 'orange', 'green', 'red']

    # Create a figure and axes for subplots
    fig, axes = plt.subplots(num_tools, num_tuples, figsize=(6 * num_tuples, 6 * num_tools))

    # Iterate over each tuple and create a separate histogram for each biome
    for i, tools in enumerate(data):
        for j, tool in enumerate(tools):
            # Set the current subplot
            ax = axes[j][i]

            # Create a histogram for the current biome with the corresponding row color
            sns.histplot(data=tool, kde=True, color=row_colors[j], ax=ax)

            ax.set_xlabel(f'Comparison of Tools - Day {alien_ids[i]}')
            ax.set_ylabel('Count')
            ax.set_title(labels[i][j])

            # Remove x-axis label from all but the bottom subplots
            if j != num_tools - 1:
                ax.set_xlabel('')

            # Remove y-axis label from all but the leftmost subplots
            if i != 0:
                ax.set_ylabel('')

            # Add legend to the top-left subplot
            if i == 0 and j == 0:
                ax.legend(labels=['Comparison'])

            # Adjust the layout to avoid overlapping
            plt.tight_layout()
            print('*')

    # Add a common y-axis label
    fig.text(0.06, 0.5, 'Count', va='center', rotation='vertical')

    plt.savefig('/home/Users/mt102/data/plots/' + name)
    plt.close()
    print(name + ' saved.')


def histograms(k, freq):
    data = []
    data_short = []

    descr = []
    descr_short = []

    for id in alien_ids:
        ggcat, ggcat_count, ggcat_short = extract(
            f"/home/Users/mt102/shared/ggcat-res/human_gut_longitudinal/alien.{id}.{k}.f{freq}")

        euler, euler_count, euler_short = extract(
            f"/home/Users/mt102/shared/ggcat-res/human_gut_longitudinal/alien.{id}.{k}.f{freq}.euler")

        path, path_count, path_short = extract(
            f"/home/Users/mt102/shared/ggcat-res/human_gut_longitudinal/alien.{id}.{k}.f{freq}.path")

        abyss, abyss_count, abyss_short = extract(
            f"/rdf/tt40/ns58/KOMB-backup/human_gut_longitudinal_samples/alien-11-{id}-0/final.unitigs.fa")

        descr.append([f'GGCAT\ncount: {len(ggcat) // 1000}K removed: {ggcat_count // 1000}K',
                      f'GGCAT(eulertigs)\ncount: {len(euler) // 1000}K removed: {euler_count // 1000}K',
                      f'GGCAT(pathtigs)\ncount: {len(path) // 1000}K removed: {path_count // 1000}K',
                      f'Abyss\ncount: {len(abyss) // 1000}K removed: {abyss_count // 1000}K'])
        descr_short.append([f'GGCAT\ncount: {len(ggcat_short) // 1000}K removed: {ggcat_count // 1000}K',
                            f'GGCAT(eulertigs)\ncount: {len(euler_short) // 1000}K removed: {euler_count // 1000}K',
                            f'GGCAT(pathtigs)\ncount: {len(path_short) // 1000}K removed: {path_count // 1000}K',
                            f'Abyss\ncount: {len(abyss_short) // 1000}K removed: {abyss_count // 1000}K'])

        data.append((ggcat, euler, path, abyss))
        data_short.append((ggcat_short, euler_short, path_short, abyss_short))

    visualize(data, descr, f'alien-hist-eul-path-k{k}-f{freq}.png')

    visualize(data_short, descr_short, f'alien-hist-eul-path-u300-k{k}-f{freq}.png')


freqs = [1, 2, 5, 10]  # [1, 2]

lens = [35]  # [21, 25, 31, 35, 41, 45, 51]

for k in lens:
    for freq in freqs:
        if k == 35:
            histograms(k, freq)

