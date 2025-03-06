import numpy as np
from Bio import SeqIO
from collections import defaultdict

###
### reads, weights, qualities
###
### When unique_modus = True:
### - Identical reads are summarized.
### - The weight represents the number of identical reads.
### - The quality scores are averaged for identical reads.
### When unique_modus = False:
### - All reads are kept separately, even if they are identical.
### - Each read has a weight of 1.
### - Quality scores are kept as they are, without averaging.

def load_and_process_reads(fname_fasta, fname_qualities, alphabet, unique_modus):
    if unique_modus:
        return process_reads_unique(fname_fasta, fname_qualities, alphabet)
    else:
        return process_reads_non_unique(fname_fasta, fname_qualities, alphabet)

def process_reads_non_unique(fname_fasta, fname_qualities, alphabet):
    with open(fname_qualities, "rb") as f:
        qualities = np.load(f, allow_pickle=True)

    alphabet_dict = {nucleotide: i for i, nucleotide in enumerate(alphabet)}
    reads_seq_binary = []
    read_mapping = []

    for idx, seq in enumerate(SeqIO.parse(fname_fasta, "fasta")):
        seq_str = str(seq.seq)
        binary_seq = np.eye(len(alphabet))[np.array([alphabet_dict.get(n, -1) for n in seq_str])]
        reads_seq_binary.append(binary_seq)
        read_mapping.append((seq.id, 1, []))  # (id, weight, identical_reads)

    # Convert to a 3D numpy array with shape (n_reads, seq_length, len(alphabet))
    reads_seq_binary = np.array(reads_seq_binary)
    reads_weights = np.ones(len(reads_seq_binary))
    qualities = np.where(qualities == 0, 2, qualities)

    return reads_seq_binary, reads_weights, qualities, read_mapping


def process_reads_unique(fname_fasta, fname_qualities, alphabet):
    with open(fname_qualities, "rb") as f:
        qualities = np.load(f, allow_pickle=True)

    alphabet_dict = {nucleotide: i for i, nucleotide in enumerate(alphabet)}
    unique_reads = {}
    read_mapping = {}

    for idx, seq in enumerate(SeqIO.parse(fname_fasta, "fasta")):
        seq_str = str(seq.seq)
        if seq_str not in unique_reads:
            binary_seq = np.eye(len(alphabet))[np.array([alphabet_dict.get(n, -1) for n in seq_str])]
            unique_reads[seq_str] = {
                'binary': binary_seq,
                'quality': qualities[idx],
                'weight': 1,
                'ids': [seq.id]
            }
            read_mapping[seq_str] = [(seq.id, 1)]
        else:
            unique_reads[seq_str]['quality'] = (unique_reads[seq_str]['quality'] * unique_reads[seq_str]['weight'] + qualities[idx]) / (unique_reads[seq_str]['weight'] + 1)
            unique_reads[seq_str]['weight'] += 1
            unique_reads[seq_str]['ids'].append(seq.id)
            read_mapping[seq_str].append((seq.id, 1))

    n_unique = len(unique_reads)
    seq_length = len(next(iter(unique_reads)))
    reads_seq_binary = np.empty((n_unique, seq_length, len(alphabet)), dtype=np.int8)
    reads_weights = np.empty(n_unique, dtype=int)
    output_qualities = np.empty((n_unique, len(qualities[0])))

    for i, (seq_str, read) in enumerate(unique_reads.items()):
        reads_seq_binary[i] = read['binary']
        reads_weights[i] = read['weight']
        output_qualities[i] = read['quality']

    output_qualities = np.where(output_qualities == 0, 2, output_qualities)
    return reads_seq_binary, reads_weights, output_qualities, read_mapping


###
### reads_log_error_proba
###

def compute_reads_log_error_matrix(
    theta, one_minus_theta, reads_seq_binary, size_alphabet
):
    log_theta = np.log(theta)
    log_theta = log_theta[:, :, np.newaxis]
    log_theta = np.tile(log_theta, (1, 1, size_alphabet))

    log_one_minus_theta = np.log(one_minus_theta)
    log_one_minus_theta = log_one_minus_theta[:, :, np.newaxis]
    log_one_minus_theta = np.tile(log_one_minus_theta, (1, 1, size_alphabet))

    final = np.einsum("NLB,NLB->NLB", log_theta, reads_seq_binary)
    final += np.einsum("NLB,NLB->NLB", log_one_minus_theta, 1 - reads_seq_binary)

    # if reads_list[n].seq_binary[l].sum(axis=0)=0 then "N" at position l then position l is ignored
    # there will be a zero in the row n,l
    # dimension: NxL
    all_N_pos = reads_seq_binary.sum(axis=2) > 0
    all_N_pos = all_N_pos[:, :, np.newaxis]
    all_N_pos = np.tile(all_N_pos, (1, 1, size_alphabet))
    # write zero where there is an "N"  in the position
    final[~all_N_pos] = 0

    return final  # dimension: NxLxB


def compute_reads_log_error_proba(qualities, reads_seq_binary, size_alphabet):
    """
    \log \theta_{n,l}^{r_{nl}^i}  \left( \frac{1-\theta_{n,l}}{B-1}\right)^{(1-r_{nl}^i)} \right)
    Do I need the weights? No this can be integrated later in update_eqs.py

    theta_{n,l} = probabiltiy that at positiion l in read n the base was called without error
    1 - theta_{n,l} = error probablity, e.g. probablity that at position l in read n the base was called erroneous.

    input-qualities:
    Q_{n,l} = confidence of the sequencer that base at position l in read n was called correctly.
    """
    theta = 1 - 10 ** (-qualities / 10)  # dimension: N X L
    one_minus_theta = (1 - theta) / (size_alphabet - 1)

    log_theta = np.log(theta)
    log_theta = log_theta[:, :, np.newaxis]
    log_theta = np.tile(log_theta, (1, 1, size_alphabet))

    log_one_minus_theta = np.log(one_minus_theta)
    log_one_minus_theta = log_one_minus_theta[:, :, np.newaxis]
    log_one_minus_theta = np.tile(log_one_minus_theta, (1, 1, size_alphabet))

    final = np.einsum("NLB,NLB->NLB", log_theta, reads_seq_binary)
    final += np.einsum("NLB,NLB->NLB", log_one_minus_theta, 1 - reads_seq_binary)

    # if reads_list[n].seq_binary[l].sum(axis=0)=0 then "N" at position l then position l is ignored
    # there will be a zero in the row n,l
    # dimension: NxL
    all_N_pos = reads_seq_binary.sum(axis=2) > 0
    all_N_pos = all_N_pos[:, :, np.newaxis]
    all_N_pos = np.tile(all_N_pos, (1, 1, size_alphabet))
    # write zero where there is an "N"  in the position
    final[~all_N_pos] = 0

    return final  # dimension: NxLxB

###
### reference
###

def load_reference_seq(reference_file, alphabet):
    for seq in SeqIO.parse(reference_file, "fasta"):
        return reference2binary(str(seq.seq).upper(), alphabet), seq.id


def reference2binary(reference_seq, alphabet):
    alphabet_dict = {nucleotide: i for i, nucleotide in enumerate(alphabet)}
    length_seq = len(reference_seq)
    reference_table = np.zeros((length_seq, len(alphabet)), dtype=np.int8)
    valid_indices = [i for i, base in enumerate(reference_seq) if base in alphabet_dict]
    reference_table[valid_indices, [alphabet_dict[reference_seq[i]] for i in valid_indices]] = 1
    return reference_table
