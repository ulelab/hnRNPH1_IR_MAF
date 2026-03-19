import argparse
import os
import subprocess
import tempfile

import numpy as np
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode


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
    seqs = []
    with open(fasta_path) as f:
        lines = [line.strip() for line in f if line.strip()]
    for i in range(0, len(lines), 2):
        seqs.append(lines[i + 1].upper())
    return seqs


def max_spliceai_score(sequence, models, context=10000):
    padded = "N" * (context // 2) + sequence + "N" * (context // 2)
    x = one_hot_encode(padded)[None, :]
    y = np.mean([m.predict(x, verbose=0) for m in models], axis=0)[0]

    if len(sequence) >= context:
        start = context // 2
        end = start + len(sequence)
        region = y[start:end, :]
    else:
        region = y[: len(sequence), :]

    acceptor_max = float(np.max(region[:, 1]))
    donor_max = float(np.max(region[:, 2]))
    return max(acceptor_max, donor_max)


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

