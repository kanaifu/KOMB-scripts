import matplotlib.pyplot as plt
from upsetplot import from_contents, UpSet

def extract_ids(fname):
    file = open(fname, 'r')
    ids = []
    for line in file.readlines():
        line = line[:-1]
        lst = line.split('\t')
        # extract info from each line
        level, id, name = lst[3], int(lst[4]), lst[-1].strip()
        if level == 'S':
            ids.append(id)

    file.close()
    return set(ids)


dir = "/home/Users/ns58/KOMB2/PRJNA878603-kraken2"
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
    fnames = [f"{dir}/{dn}.kraken2_report.txt" for dn in fnames]
    taxids[group] = set([])
    for fname in fnames:
        taxids[group].update(extract_ids(fname))

# Prepare the data for the UpSet plot
data = from_contents(taxids)

# Create and display the UpSet plot
upset = UpSet(data)
upset.plot()
plt.savefig("try3.png")
