import numpy as np
import pandas as pd
from statistics import mean
from statistics import stdev
import math
from collections import Counter
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from collections import defaultdict
import glob
import pickle

CB91_Blue = '#2CBDFE'
CB91_Green = '#47DBCD'
CB91_Pink = '#F3A0F2'
CB91_Purple = '#9D2EC5'
CB91_Violet = '#661D98'
CB91_Amber = '#F5B14C'

main_dir = "/home/Users/ns58/KOMB2/PRJNA878603-Output"
info_path = "/home/Users/ns58/KOMB2/Data/PRJNA878603_SraRunTable.txt"

bucket = {}  # (gender, disease length) : [SRA IDs]


with open(info_path, "r") as file:
    lines = file.readlines()
    iterable = iter(lines)
    next(iterable)  # skip first entry
    for line in iterable:
        data = line[:-1].split(',')  # 0: SRA ID, 24: gender, 29: disease_length
        if len(data) < 30:
            print("ERROR")
            continue
        if data[29] not in bucket:
            bucket[data[29]] = []
        bucket[data[29]].append(data[0])

dpi = 300

corenesses = {}

datasets = []
for names in [bucket['Control'], bucket['Short'], bucket['Long']]:
    tables = []
    links = [f"{main_dir}/{dn}-k51-l150-c2" for dn in names]

    for fname, flink in zip(names, links):
        # Data Load
        data = np.loadtxt(flink + "/kcore.tsv", dtype=int)
        data = np.delete(data, 1, axis=1)
        data = np.delete(data, 0, axis=1)

        data_coreA = np.loadtxt(f"{flink}/CoreA_anomaly.txt")
        data_coreA = np.delete(data_coreA, 0, axis=1)
        mean_value = np.mean(data_coreA)
        std_value = np.std(data_coreA)

        # Step 2: Create a new array with values 1 or 0 based on the condition.
        coreA_new = np.where(data_coreA >= mean_value + (3 * std_value), 1, 0)

        # Now 'new_arr' contains 1 where the element is greater than or equal to the mean + 3std, and 0 otherwise.
        print("Mean:", mean_value)
        print("Standard Deviation:", std_value)

        data = np.concatenate((data, coreA_new.reshape(-1, 1)), axis=1)

        tables.append(data)

    datasets.append(np.concatenate(tables, axis=0))

# _cmap = sns.color_palette("vlag", as_cmap=True)


def make_plots(data, ax, sname):
    coreA_metric = pd.DataFrame(data=data, columns=["Coreness", "Degree", "CoreA"], copy=True)

    degree_core_count_data = np.zeros_like(data)
    degree_core_dict = defaultdict(int)
    for x in data:
        if x[1] < 0:
            x[1] = 0
        degree_core_dict[(x[1], x[0])] += 1
    for k, t in enumerate(degree_core_dict):
        i, j = t
        degree_core_count_data[k] = i, j, degree_core_dict[(i, j)]
        if i == 0 and j == 0:
            degree_core_count_data[k] = 0, 0, 0
    degree_core_count_data = degree_core_count_data[~np.all(degree_core_count_data == 0, axis=1)]

    # subsample to make plotting faster
    # coreA_sample = coreA_metric.sample(n=1000000, random_state=42)

    # Plotting
    # Degree - core - CoreA
    ax.plot([0, 10000], [0, 10000], color='k', lw=3, zorder=-1)
    scatter = ax.scatter(coreA_metric["Degree"], coreA_metric["Coreness"], c=coreA_metric["CoreA"], cmap="coolwarm",
                         vmin=0, vmax=1)
    ax.set_facecolor("white")
    ax.tick_params(labelsize=13)
    ax.set_xlabel("Degree", fontsize=15)
    ax.set_ylabel("Coreness", fontsize=15)
    ax.set_xlim((-50, 3000))
    ax.set_ylim((-50, 1000))
    ax.set_title(sname, fontsize=17, pad=11)

    return scatter


gs = gridspec.GridSpec(1, 3)
fig = plt.figure(figsize=(23, 10))
ax1, ax2, ax3 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[0, 2])
samples = ["Control", "ME/CFS (Short)", "ME/CFS (Long)"]
axes = [ax1, ax2, ax3]

for d, a, name in zip(datasets, axes, ['Control', 'ME/CFS (Short)', 'ME/CFS (Long)']):
    make_plots(d, a, name)

plt.tight_layout()
plt.subplots_adjust(wspace=0.5, top=0.8)

norm = mpl.colors.Normalize(vmin=0, vmax=1)

# Add a common legend to the right of all subplots
legend_elements = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=25, label='Normal'),
                   plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=25, label='Mean + 3 std')]
fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.95, 0.5), title='Classification', fontsize=25)

# Adjust layout for the legend
plt.subplots_adjust(right=0.85)

ax1.text(-0.13, 1.25, "B", transform=ax1.transAxes,
         size=23, weight='bold')
# ax2.text(-0.13, 1.07, "F", transform=ax2.transAxes,
# size=23, weight='bold')
# ax3.text(-0.13, 1.07, "G", transform=ax3.transAxes,
# size=23, weight='bold')
# ax4.text(-0.13, 1.07, "H", transform=ax4.transAxes,
# size=23, weight='bold')

# plt.tight_layout()
plt.suptitle("Anomaly Profiles of PRJNA samples", fontsize=23)
plt.savefig("KOMB_PRJNA_Coreness_Degree_cmap_mean+3std.png", dpi=300)

