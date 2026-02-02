import re
import pandas as pd

# Input settings
input_path = "/scratch/prj/ppn_rnp_networks/users/mike.jones/data/splice/CrypticSS_All_Inferences.tsv"
output_path = "/scratch/prj/ppn_rnp_networks/users/mike.jones/data/splice/CrypticSS_All.filtered.bedgraph"
logfile = "/scratch/prj/ppn_rnp_networks/users/mike.jones/software/header_parse_failures.log"

# Format: choose "matrix" (CrypticSS) or "fasta" (RBPNet)
input_format = "matrix"  # <-- change accordingly

# Logging
log = open(logfile, "w")
bedgraph_records = []

with open(input_path, "r") as infile:
    lines = infile.readlines()

    for idx, line in enumerate(lines):
        line = line.strip()
        if line == "":
            continue

        if input_format == "fasta" and not line.startswith(">"):
            continue
        if input_format == "matrix" and (line.startswith("ID") or line.startswith("#")):
            continue

        if input_format == "fasta":
            header = line.strip()[1:]
            profile_target = lines[idx + 4].strip().split()

            match = re.match(r'^([^:]+)::([^:]+):(\d+)-(\d+)\(([+-])\)', header)
            if not match:
                log.write(f"Header parse error on line {idx+1}: {header}\n")
                continue

            gene, chrom, start, end, strand = match.groups()
            start = int(start)
            end = int(end)

            for j, score in enumerate(profile_target):
                score = float(score)
                if score < 1e-6:
                    continue

                if strand == "+":
                    pos_start = start + j
                    signed_score = score
                else:
                    pos_start = end - (j + 1)
                    signed_score = -score

                bedgraph_records.append([f"{chrom}", pos_start, pos_start+1, signed_score])

        elif input_format == "matrix":
            fields = line.strip().split("\t")
            intron_id = fields[0]
            prediction_blocks = fields[1:]

            try:
                gene, loc = intron_id.split("::")
                chrom_region, strand = loc.strip(")").split("(")
                chrom, coords = chrom_region.split(":")
                start_coord, end_coord = map(int, coords.split("-"))
            except Exception as e:
                log.write(f"Header parse error on line {idx+1}: {intron_id} ({e})\n")
                continue

            for pred in prediction_blocks:
                parts = pred.strip().split(",")
                if len(parts) != 3:
                    continue

                acceptor, donor, pos = parts

                try:
                    pos = int(pos)
                    donor = float(donor)
                    if donor < 1e-6:
                        continue

                    if strand == "+":
                        pos_start = start_coord + (pos - 1)
                        signed_score = donor
                    else:
                        pos_start = end_coord - pos
                        signed_score = -donor

                    bedgraph_records.append([f"{chrom}", pos_start, pos_start+1, signed_score])

                except Exception as e:
                    log.write(f"Prediction parse error: {pred} for {intron_id} ({e})\n")

log.close()

# Convert to pandas for sorting and collapsing
df = pd.DataFrame(bedgraph_records, columns=["chrom", "start", "end", "score"])
df = df.sort_values(by=["chrom", "start", "end"])

# Collapse duplicates by averaging signal
collapsed = df.groupby(["chrom", "start", "end"], as_index=False).agg({"score": "mean"})

# Output final bedGraph
collapsed.to_csv(output_path, sep='\t', header=False, index=False, float_format='%.6f')

print(f"✅ Strand-aware bedGraph written: {output_path}")

