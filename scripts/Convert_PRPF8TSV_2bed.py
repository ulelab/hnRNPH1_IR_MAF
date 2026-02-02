import re

infile = "/scratch/prj/ppn_rnp_networks/users/mike.jones/data/rbpnet/PRPF8_eCLIP_RBPbinding_Prediction.tsv"
outfile = "PRPF8_eCLIP_RBPbinding_Prediction.bed"

with open(infile, 'r') as fin, open(outfile, 'w') as fout:
    lines = fin.readlines()

    for i in range(0, len(lines), 6):
        header = lines[i].strip()
        profile_target = lines[i+4].strip().split()

        # Extract metadata from header
        # Example: >PSMC6::14:52718988-52720910(+)
        match = re.match(r'^>([^:]+)::([^:]+):(\d+)-(\d+)\(([+-])\)', header)
        if not match:
            print(f"❌ Header format error on line {i+1}: {header}")
            continue

        gene, chrom, start, end, strand = match.groups()
        start = int(start)
        end = int(end)

        # Reverse profile if on minus strand
        if strand == '-':
            profile_target = profile_target[::-1]
            for j, score in enumerate(profile_target):
                pos_end = end - j
                pos_start = pos_end - 1
                bed_fields = [f"chr{chrom}", str(pos_start), str(pos_end), gene, score, strand]
                fout.write('\t'.join(bed_fields) + '\n')
        else:
            for j, score in enumerate(profile_target):
                pos_start = start + j
                pos_end = pos_start + 1
                bed_fields = [f"chr{chrom}", str(pos_start), str(pos_end), gene, score, strand]
                fout.write('\t'.join(bed_fields) + '\n')

print(f"✅ BED6 file written: {outfile}")
