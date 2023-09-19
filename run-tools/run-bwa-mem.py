import os

info_path = "/home/Users/ns58/KOMB2/Data/PRJNA878603_SraRunTable.txt"

bucket = {}  # disease length : [SRA IDs]
sample_to_condition = {}  # SRA ID : disease length

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
unitig_paths = [f"/home/Users/ns58/KOMB2/PRJNA878603-Output/{sample}-k51-l150-c2/unitigs.l150.fasta" for sample in samples]
read1_paths = [f"/home/Users/ns58/KOMB2/Data/PRJNA878603-ME_CFS-NovaSeq/{sample}_1.fastq" for sample in samples]
read2_paths = [f"/home/Users/ns58/KOMB2/Data/PRJNA878603-ME_CFS-NovaSeq/{sample}_2.fastq" for sample in samples]
output1_paths = [f"/home/Users/mt102/shared/KOMB-revision/SAM-alignments/{sample}_1.sam" for sample in samples]
output2_paths = [f"/home/Users/mt102/shared/KOMB-revision/SAM-alignments/{sample}_2.sam" for sample in samples]

for sample, unitig_path, read1_path, read2_path, output1_path, output2_path in zip(samples, unitig_paths, read1_paths,
                                                                                   read2_paths, output1_paths, output2_paths):
    os.system(f"bwa mem -a -k 18 -t 60 {unitig_path} {read1_path} -o {output1_path}")
    os.system(f"bwa mem -a -k 18 -t 60 {unitig_path} {read2_path} -o {output2_path}")
