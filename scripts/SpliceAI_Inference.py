import argparse
import os
import subprocess
import tempfile

import numpy as np
from keras.models import load_model  # type: ignore
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode  # type: ignore


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bed", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--out", required=True)
    return parser.parse_args()


def load_spliceai_models():
    paths = (f"models/spliceai{x}.h5" for x in range(1, 6))
    return [load_model(resource_filename("spliceai", p), compile=False) for p in paths]


def read_fasta_sequences(fasta_path):
    # Parse FASTA robustly (bedtools FASTA may or may not be single-line)
    seqs = []
    current = []
    in_seq = False
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_seq:
                    seqs.append("".join(current).upper())
                current = []
                in_seq = True
                continue
            current.append(line)
        if current:
            seqs.append("".join(current).upper())
        elif in_seq:
            # Handle the case where a FASTA header exists but the sequence body is empty.
            seqs.append("")
    return seqs


def max_spliceai_score(sequence, models, context=10000):
    sequence_clean = sequence.replace("-", "").replace(" ", "").replace("*", "").upper()
    if len(sequence_clean) == 0:
        return -1.0
    sequence_no_n = sequence_clean.replace("N", "").replace("n", "")
    if len(sequence_no_n) == 0:
        return -1.0

    padded = "N" * (context // 2) + sequence_clean + "N" * (context // 2)
    x = one_hot_encode(padded)[None, :]
    y = np.mean([m.predict(x, verbose=0) for m in models], axis=0)[0]

    if len(sequence_clean) >= context:
        start = context // 2
        end = start + len(sequence_clean)
        region = y[start:end, :]
    else:
        region = y[: len(sequence_clean), :]

    # Donor-only, ignore acceptor scores
    donor_max = float(np.max(region[:, 2]))
    return donor_max


def main():
    args = parse_args()
    models = load_spliceai_models()

    with tempfile.TemporaryDirectory() as tmpdir:
        slopped_bed = os.path.join(tmpdir, "slop24.bed")
        seq_fa = os.path.join(tmpdir, "slop24.fa")

        subprocess.run(
            ["bedtools", "slop", "-i", args.bed, "-g", args.genome, "-b", "24"],
            stdout=open(slopped_bed, "w"),
        )

        subprocess.run(
            [
                "bedtools",
                "getfasta",
                "-fi",
                args.fasta,
                "-bed",
                slopped_bed,
                "-s",
                "-name",
            ],
            stdout=open(seq_fa, "w"),
        )

        seqs = read_fasta_sequences(seq_fa)
        scores = [max_spliceai_score(seq, models) for seq in seqs]

    with open(args.bed) as fin, open(args.out, "w") as fout:
        for line, score in zip(fin, scores):
            cols = line.rstrip("\n").split("\t")
            while len(cols) < 5:
                cols.append(".")
            cols[4] = f"{score:.6f}"
            fout.write("\t".join(cols) + "\n")


if __name__ == "__main__":
    main()

