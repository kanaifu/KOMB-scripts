import numpy as np
import pandas as pd
from statistics import mean
from statistics import stdev
import math
from collections import Counter
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
        data = np.loadtxt(flink + "/kcore.tsv", dtype=int) # #VID    Name    Coreness        Degree
        corenesses[fname] = data[:, 2].tolist()

        print(corenesses[fname][:5])

sizes = [len(bucket['Control']), len(bucket['Short']), len(bucket['Long'])]
sample = "PRJNA"
threshold = 100

fig, axs = plt.subplots(figsize=(32, 10))

pos = [1 for i in range(sizes[0])] + [1.5 for i in range(sizes[1])] + [2 for i in range(sizes[2])]

data = []

for lst in [bucket['Control'], bucket['Short'], bucket['Long']]:
    for fname in lst:
        data.append(corenesses[fname])

violins = axs.violinplot(data, pos, points=300, widths=0.3, showextrema=False)

for pc in violins['bodies'][:sizes[0]:]:
    pc.set_facecolor(CB91_Purple)
for pc in violins['bodies'][sizes[0]:sizes[0]+sizes[1]:]:
    pc.set_facecolor(CB91_Green)
for pc in violins['bodies'][sizes[0] + sizes[1]::]:
    pc.set_facecolor(CB91_Amber)

axs.set_xticks([1, 1.5, 2])
axs.set_xticklabels(["Control", "ME/CFS (Short)", "ME/CFS (Long)"])
axs.tick_params(axis="both", labelsize=35)
axs.set_ylabel("Shell number", fontsize=40, labelpad=25)
axs.set_xlabel("Diseased Type", fontsize=40, labelpad=25)

# plt.legend(bbox_to_anchor=(1.01, 1), loc='upper left')
fig.tight_layout()
plt.ylim((0, threshold))
axs.text(-0.05, 1.15, "A", transform=axs.transAxes,
         size=46, weight='bold')
plt.title("KOMB Profiles of ME/CFS study samples for 3 different diseased levels", fontsize=43, pad=30)
fig.savefig("shell-size-{}-violin-thr-{}.png".format(sample, threshold),
            bbox_inches="tight", dpi=300)

