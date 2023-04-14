from .local_haplotype_inference.use_quality_scores.update_eqs import update_mean_haplo as uqs_update_mean_haplo
from .local_haplotype_inference.use_quality_scores.analyze_results import compute_unique_haplo_posterior
from .local_haplotype_inference.use_quality_scores import preparation as uqs_preparation
from .local_haplotype_inference.learn_error_params import preparation as lep_preparation
from .local_haplotype_inference.learn_error_params.update_eqs import update_mean_haplo as lep_update_mean_haplo

from Bio import SeqIO
import numpy as np
from typing import Optional, TextIO
import re
import math


alphabet = "ACGT-" # TODO hardcoded

def create_unique_haplo_var(support_handle):
    arr = []
    for i in SeqIO.parse(support_handle, "fasta"):
        arr.append(str(i.seq).upper())

    return arr

def ingest_sampler_output_to_calc_mean_cluster(haplo: list[str], corrected_handle: TextIO, N_s):

    mean_cluster = np.zeros((N_s, len(haplo)), dtype=int) # N_s x k

    for c_idx, c in enumerate(SeqIO.parse(corrected_handle, "fasta")):
        c_seq = str(c.seq).upper()
        for s_idx, s_seq in enumerate(haplo):
            if c_seq == s_seq:
                # probability that read i belongs to cluster j, only 0 and 1
                mean_cluster[c_idx][s_idx] = 1
                continue

    return mean_cluster

def parse_gamma_or_theta(datafile, var):
    for line in datafile:
        theta_or_gamma = re.search(f"^#{var}\\s=\\s.*", line)
        if theta_or_gamma != None:
            return float(theta_or_gamma.group(0).split("=")[1])

    raise ValueError("non theta and gamma found")

def ingest_sampler_results_gamma_theta(results_file_handle, inference_type):
    datafile = results_file_handle.readlines()
    if inference_type == "shorah":
        gamma = parse_gamma_or_theta(datafile, "gamma")
        theta = parse_gamma_or_theta(datafile, "theta")
        return (math.log(gamma), math.log(1-gamma)), (math.log(theta), math.log(1-theta))

    raise NotImplementedError() # TODO

def post_process_pooled_samples_mode(fref_in: str, freads_in: str,
                                     results_file_handle: TextIO,
                                     support_handle: TextIO,
                                     corrected_handle: TextIO,
                                     inference_type: str,
                                     fname_qualities: Optional[TextIO] = None) -> np.ndarray:

    assert inference_type == "use_quality_scores" or fname_qualities == None
    unique_modus = False # even off when on in sampler stage

    if inference_type == "use_quality_scores":
        reads_list, qualities = uqs_preparation.load_fasta_and_qualities(
            freads_in, # N_s - filter for sample
            fname_qualities,
            alphabet,
            unique_modus
        )
    else:
        reads_list = lep_preparation.load_fasta2reads_list(freads_in, alphabet, unique_modus)

    reads_seq_binary, reads_weights = uqs_preparation.reads_list_to_array(reads_list)

    assert sum(reads_weights) == len(reads_list) # TODO

    reference_binary = uqs_preparation.load_reference_seq(fref_in, alphabet)[0]
    unique_haplo_var = create_unique_haplo_var(support_handle)
    mean_z = ingest_sampler_output_to_calc_mean_cluster(unique_haplo_var, corrected_handle, len(reads_list))
    print("mean_z", mean_z)

    if inference_type == "use_quality_scores":
        reads_log_error_proba = uqs_preparation.compute_reads_log_error_proba(
            qualities, reads_seq_binary, len(alphabet)
        )
        mean_h = uqs_update_mean_haplo(
            reads_weights, # all 1, number of reads in window, numpy array?
            reference_binary, # from preparation -> load_reference_seq
            reads_log_error_proba, # L51 run dpm mfa
            mean_z,
            mean_log_gamma
        )
    else:
        mean_log_gamma, mean_log_theta = ingest_sampler_results_gamma_theta(
            results_file_handle,
            inference_type
        )

        mean_h = lep_update_mean_haplo(
            reads_seq_binary,
            reads_weights,
            reference_binary,
            mean_z,
            mean_log_theta,
            mean_log_gamma,
        )

    posterior = compute_unique_haplo_posterior( # dim k
        mean_h,
        unique_haplo_var, # dim k, each haplo is one string
        alphabet
    )

    return posterior, np.sum(mean_z, axis=0)

    # TODO write posterior in new support files
    # TODO if average reads or posterior == 0 do not list haplo in support file
    # TODO assert ave_reads == N_s - haplo that are 0