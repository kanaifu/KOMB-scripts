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


def visualize(data, labels, do_abyss, title, name):
    # Create a figure and axes for subplots
    val = 15
    if do_abyss:
        val = 20
    fig, axes = plt.subplots(len(data), 1, figsize=(val, 6 * len(data)))

    # Iterate over each tuple and create a DataFrame and violin plot
    for i, lists in enumerate(data):

        # Create a DataFrame for the current tuple
        if do_abyss:
            df = pd.DataFrame({
                'tool': ['Cuttle\n' + labels[i][0]] * len(lists[0]) + ['GGCAT\n' + labels[i][1]] * len(lists[1]) + [
                    'Abyss\n' + labels[i][2]] * len(lists[2]),
                'class': lists[0] + lists[1] + lists[2]
            })
        else:
            df = pd.DataFrame({
                'tool': ['Cuttle\n' + labels[i][0]] * len(lists[0]) + ['GGCAT\n' + labels[i][1]] * len(lists[1]),
                'class': lists[0] + lists[1]
            })
        # Set the current subplot
        ax = axes[i]

        # Create the violin plot using seaborn
        sns.violinplot(data=df, x='tool', y='class', scale='count', ax=ax)
        sns.swarmplot(data=df, x='tool', y='class', color='k', alpha=0.6, ax=ax)

        ax.set_xlabel('Tools')
        ax.set_ylabel('Lengths')
        ax.set_title(f'Comparison of Tools - Day {alien_ids[i]}')

    # Adjust the layout to avoid overlapping
    plt.tight_layout()
    plt.title(title)
    plt.savefig('/home/Users/mt102/data/plots/' + name)
    plt.close()
    print(name + ' saved.')


def violins(k, freq, do_abyss):
    data = []
    data_short = []

    descr = []
    descr_short = []

    for id in alien_ids:

        cuttle, cuttle_count, cuttle_short = extract(
            f"/home/Users/mt102/shared/cuttle-res/human_gut_longitudinal/alien.{id}.{k}.f{freq}.fa")
        ggcat, ggcat_count, ggcat_short = extract(
            f"/home/Users/mt102/shared/ggcat-res/human_gut_longitudinal/alien.{id}.{k}.f{freq}")

        if do_abyss:
            abyss, abyss_count, abyss_short = extract(
                f"/rdf/tt40/ns58/KOMB-backup/human_gut_longitudinal_samples/alien-11-{id}-0/final.unitigs.fa")

            descr.append([f'count: {len(cuttle) // 1000}K removed: {cuttle_count // 1000}K',
                          f'count: {len(ggcat) // 1000}K removed: {ggcat_count // 1000}K',
                          f'count: {len(abyss) // 1000}K removed: {abyss_count // 1000}K'])
            descr_short.append([f'count: {len(cuttle_short) // 1000}K removed: {cuttle_count // 1000}K',
                                f'count: {len(ggcat_short) // 1000}K removed: {ggcat_count // 1000}K',
                                f'count: {len(abyss_short) // 1000}K removed: {abyss_count // 1000}K'])

            data.append((cuttle, ggcat, abyss))
            data_short.append((cuttle_short, ggcat_short, abyss_short))

        else:
            descr.append([f'count: {len(cuttle) // 1000}K removed: {cuttle_count // 1000}K',
                          f'count: {len(ggcat) // 1000}K removed: {ggcat_count // 1000}K'])
            descr_short.append([f'count: {len(cuttle_short) // 1000}K removed: {cuttle_count // 1000}K',
                                f'count: {len(ggcat_short) // 1000}K removed: {ggcat_count // 1000}K'])

            data.append((cuttle, ggcat))
            data_short.append((cuttle_short, ggcat_short))

    visualize(data, descr, do_abyss, f'Alien-11 unitig distribution, k = {k}, freq = {freq}', f'alien-k{k}-f{freq}.png')

    visualize(data_short, descr_short, do_abyss,
              f'Alien-11 unitig distribution, k = {k}, freq = {freq}, unitig len <= 300',
              f'alien-u300-k{k}-f{freq}.png')


freqs = [15, 20, 25]  # [1, 2]

lens = [35]  # [21, 25, 31, 35, 41, 45, 51]

for k in lens:
    for freq in freqs:
        if k == 35:
            violins(k, freq, True)
        else:
            violins(k, freq, False)

