import pandas as pd

info_path = "/home/Users/ns58/KOMB2/Data/PRJNA878603_SraRunTable.txt"

bucket = {}  # disease length : [SRA IDs]
sample_to_condition = {} # SRA ID : disease length

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
        sample_to_condition[data[0]] = data[29]
        
samples = bucket['Control'] + bucket['Long'] + bucket['Short']
sample_paths = [f"/home/Users/ns58/KOMB2/PRJNA878603-Output/{sample}-k51-l150-c2" for sample in samples]

sam1_paths = [f"/home/Users/mt102/shared/KOMB-revision/SAM-alignments/{sample}_1.sam" for sample in samples]
sam2_paths = [f"/home/Users/mt102/shared/KOMB-revision/SAM-alignments/{sample}_2.sam" for sample in samples]
for sample, sam1_path, sam2_path in zip(samples, sam1_paths, sam2_paths):
    mapped_counts = {}
    with open(sam1_path, 'r') as file:
        for line in file:
            if not line.startswith('@'):
                fields = line.split('\t')
                unitig_id = fields[2]  # Assuming the sequence ID is in the third column (column 2)

                # Increment the count for the sequence ID in the dictionary
                if unitig_id in mapped_counts:
                    mapped_counts[unitig_id] += 1
                else:
                    mapped_counts[unitig_id] = 1
    
    with open(sam2_path, 'r') as file:
        for line in file:
            if not line.startswith('@'):
                fields = line.split('\t')
                unitig_id = fields[2]  # Assuming the sequence ID is in the third column (column 2)

                # Increment the count for the sequence ID in the dictionary
                if unitig_id in mapped_counts:
                    mapped_counts[unitig_id] += 1
                else:
                    mapped_counts[unitig_id] = 1
                    
    df = pd.DataFrame(list(mapped_counts.items()), columns=['Unitig Name', 'Read Count'])
    df.to_csv(f"/home/Users/mt102/shared/KOMB-revision/mapping-summary/{sample}_summary.tsv", sep='\t')
