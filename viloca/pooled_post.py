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
import pickle

alphabet = "ACGT-" # TODO hardcoded

# TODO if average reads or posterior == 0 do not list haplo in support file
# TODO assert ave_reads == N_s - haplo that are 0

def _create_unique_haplo_var(support_handle):
    arr = []
    for i in SeqIO.parse(support_handle, "fasta"):
        arr.append(str(i.seq).upper())

    return arr

def _ingest_sampler_output_to_calc_mean_cluster(haplo: list[str], corrected_handle: TextIO):

    mean_cluster = np.empty((0, len(haplo)), dtype=int) # N_s x k

    for c in SeqIO.parse(corrected_handle, "fasta"):
        c_seq = str(c.seq).upper()
        row = np.zeros(len(haplo), dtype=int)
        for s_idx, s_seq in enumerate(haplo):
            if c_seq == s_seq:
                # probability that read i belongs to cluster j, only 0 and 1
                row[s_idx] = 1
                mean_cluster = np.vstack((mean_cluster, row))
                continue

    return mean_cluster

def _parse_gamma_or_theta(datafile, var):
    for line in datafile:
        theta_or_gamma = re.search(f"^#{var}\\s=\\s.*", line)
        if theta_or_gamma != None:
            return float(theta_or_gamma.group(0).split("=")[1])

    raise ValueError("non theta and gamma found")

def _ingest_sampler_results_gamma_theta(results_file_handle, inference_type):
    if inference_type == "shorah":
        datafile = results_file_handle.readlines()
        gamma = _parse_gamma_or_theta(datafile, "gamma")
        theta = _parse_gamma_or_theta(datafile, "theta")
        return (math.log(gamma), math.log(1-gamma)), (math.log(theta), math.log(1-theta))
    else:
        with results_file_handle as f:
            data = pickle.load(f)[0][1]
            if inference_type == "learn_error_params":
                return data["mean_log_gamma"], data["mean_log_theta"]
            if inference_type == "use_quality_scores":
                return data["mean_log_gamma"], None
            raise NotImplementedError()



def write_support_file_per_sample(support_handle: TextIO,
                                new_support_handle: TextIO,
                                posterior: np.ndarray,
                                average_reads:  np.ndarray):
    """Writes new support file with updated posterior and average reads

    Args:
        support_handle: File to be updated. Will not be overwritten.
        new_support_handle: Write destination.
        posterior: Dimension: 1 x number of haplotypes.
        average_reads: Dimension: 1 x number of haplotypes.
    """
    records = []
    for idx, s in enumerate(SeqIO.parse(support_handle, "fasta")):
        s.id = s.name = s.id.split("=")[0] + "=" + str(posterior[idx])
        s.description = (s.description.split("=")[0] + "=" + str(posterior[idx])
            + " ave_reads=" + str(average_reads[idx]))
        if average_reads[idx] != 0:
            records.append(s)

    SeqIO.write(records, new_support_handle, "fasta")


def recalculate_posterior_and_ave_reads(fref_in: str, freads_in: str,
                                     results_file_handle: TextIO,
                                     support_handle: TextIO,
                                     corrected_handle: TextIO,
                                     inference_type: str,
                                     fname_qualities: Optional[TextIO] = None
    ) -> tuple[np.ndarray, np.ndarray]:
    """Takes output of sampler to calculate new posterior value and average reads

    Args:
        fref_in: Reference in window.
        freads_in: Reads in window.
        results_file_handle: A sampler specific file. Depends on inference_type.
        support_handle:
        corrected_handle:
        inference_type:
        fname_qualities: File with quality scores. Only relevant for the
            inference_type use_quality_scores.

    Returns:
        tuple: updated posterior and average_reads
    """

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

    assert sum(reads_weights) == len(reads_list)

    reference_binary = uqs_preparation.load_reference_seq(fref_in, alphabet)[0]
    unique_haplo_var = _create_unique_haplo_var(support_handle)
    mean_z = _ingest_sampler_output_to_calc_mean_cluster(unique_haplo_var, corrected_handle)

    if inference_type == "use_quality_scores":
        mean_log_gamma, _ = _ingest_sampler_results_gamma_theta(
            results_file_handle,
            inference_type
        )
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
        mean_log_gamma, mean_log_theta = _ingest_sampler_results_gamma_theta(
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

def filter_fasta(fw, file_to_filter, sample_name):
    lock = False
    with open(file_to_filter) as f:
        for line in f:
            if line.startswith(f">__#{sample_name}#__") or lock == True:
                fw.write(line)
                lock = not lock
