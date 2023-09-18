import matplotlib.pyplot as plt
from matplotlib_venn import venn3_unweighted

def extract_ids(fname):
    file = open(fname, 'r')
    ids = []
    for line in file.readlines():
        line = line[:-1]
        lst = line.split('\t')
        # extract info from each line
        level, id, name = lst[3], int(lst[4]), lst[-1].strip()
        if level == 'G':
            ids.append(id)

    file.close()
    return set(ids)


dir = "/home/Users/mt102/shared/KOMB/kraken2-trusses"
info_path = "/home/Users/ns58/KOMB2/Data/PRJNA878603_SraRunTable.txt"

bucket = {}  # disease length : [SRA IDs]

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

dataset_names = [('control', bucket['Control']), ('short', bucket['Short']), ('long', bucket['Long'])]
taxids = {}

for group, fnames in dataset_names:
    fnames = [f"{dir}/{dn}-k51-l150-c2/truss.kraken2_report.txt" for dn in fnames]
    taxids[group] = set([])
    for fname in fnames:
        taxids[group].update(extract_ids(fname))


fig, ax = plt.subplots(figsize=(5, 5))
ids = ['100', '010', '001',
       '110', '101', '011', '111']

v = venn3_unweighted([taxids['control'], taxids['short'], taxids['long']],
                     set_labels=["Control", "ME/CFS (Short)", "ME/CFS (Long)"])

p = v.get_patch_by_id('100')
p.set_color('#117733')

p = v.get_patch_by_id('010')
p.set_color('#cc6677')

p = v.get_patch_by_id('001')
p.set_color('#ddcc77')

p = v.get_patch_by_id('110')
p.set_color('#994455')

p = v.get_patch_by_id('101')
p.set_color('#999933')

p = v.get_patch_by_id('011')
p.set_color('#ee8866')

p = v.get_patch_by_id('111')
p.set_color('#888888')

for i in ids:
    l = v.get_label_by_id(i)
    l.set_fontsize(16)

for l in v.set_labels:
    l.set_fontsize(18)

plt.tight_layout()
fig.savefig("PRJNA878603-MECFS-venn-G-truss-condition-dpi300.png", dpi=300)

