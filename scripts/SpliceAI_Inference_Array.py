import os
import sys
import numpy as np
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai.utils import one_hot_encode

# Silence TensorFlow output and progress bars
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
import tensorflow as tf
tf.get_logger().setLevel('ERROR')

# Paths
input_tsv = '/scratch/prj/ppn_rnp_networks/users/mike.jones/data/splice/batches.tsv'
input_dir = '/scratch/prj/ppn_rnp_networks/users/mike.jones/data/splice/split_batches'
output_dir = '/scratch/prj/ppn_rnp_networks/users/mike.jones/software/spliceAI/inf'

# Load models
context = 10000
paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(resource_filename('spliceai', x), compile=False) for x in paths]
predict_fn = [model.predict for model in models]

# Get SLURM array index
array_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 1))
with open(input_tsv) as f:
    fasta_file = f.readlines()[array_id - 1].strip()

fasta_path = os.path.join(input_dir, fasta_file)

# Read FASTA and run predictions
output_lines = []

with open(fasta_path) as f:
    lines = f.readlines()

for i in range(0, len(lines), 2):
    header = lines[i].strip()[1:]  # Remove '>'
    sequence = lines[i+1].strip().upper()

    padded_seq = 'N' * (context // 2) + sequence + 'N' * (context // 2)
    x = one_hot_encode(padded_seq)[None, :]

    y = np.mean([fn(x, verbose=0) for fn in predict_fn], axis=0)
    scores = y[0]  # shape (len + context, 4)

    # Context trimming logic
    if len(sequence) >= context:
        start = context // 2
        end = start + len(sequence)
        acceptor_probs = scores[start:end, 1]
        donor_probs = scores[start:end, 2]
    else:
        acceptor_probs = scores[:len(sequence), 1]
        donor_probs = scores[:len(sequence), 2]

    values = [f"{acc:.6f},{don:.6f},{pos}" for pos, (acc, don) in enumerate(zip(acceptor_probs, donor_probs), start=1)]
    output_lines.append(f"{header}\t" + "\t".join(values))

# Write to file
basename = os.path.splitext(fasta_file)[0]
outfile = os.path.join(output_dir, f"{basename}_triplets.tsv")
with open(outfile, 'w') as f:
    f.write("ID\t[Acceptor,Donor,Position]...\n")
    f.write("\n".join(output_lines))

