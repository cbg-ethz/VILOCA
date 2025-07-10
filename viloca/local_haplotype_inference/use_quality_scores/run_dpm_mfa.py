#!/usr/bin/env python3

import sys
import os
import logging
import pickle
from numpy.random import default_rng

# my python-scripts
from . import preparation
from . import analyze_results
from . import cavi

logging.basicConfig(
    filename="viloca.log", encoding="utf-8", level=logging.INFO
)

def main(
    freads_in,
    fref_in,
    fname_qualities,
    output_dir,
    n_starts,
    K,
    alpha0,
    alphabet="ACGT-",
    unique_modus=True,
    convergence_threshold=1e-03,
    record_history=False,
    seed=27,
):

    window_id = freads_in.split("/")[-1][:-4]  # freads_in is absolute path
    output_name = output_dir + window_id + "-"

    os.makedirs(output_dir, exist_ok=True)

    reference_binary, ref_id = preparation.load_reference_seq(fref_in, alphabet)
    reads_seq_binary, reads_weights, qualities, read_mapping = preparation.load_and_process_reads(freads_in, fname_qualities, alphabet, unique_modus)
    reads_log_error_proba = preparation.compute_reads_log_error_proba(
        qualities, reads_seq_binary, len(alphabet)
    )

    if n_starts >1:
        result_list = cavi.multistart_cavi(
            K,
            alpha0,
            alphabet,
            reference_binary,
            reads_seq_binary,
            reads_weights,
            reads_log_error_proba,
            n_starts,
            output_name,
            convergence_threshold,
            record_history,
            seed
        )

    else:
        result_list = [
            cavi.run_cavi(
                K,
                alpha0,
                alphabet,
                reference_binary,
                reads_seq_binary,
                reads_weights,
                reads_log_error_proba,
                0,
                output_name,
                convergence_threshold,
                record_history,
                default_rng(seed=seed),
                seed
            )
        ]

    logging.info("reference " + fref_in)
    logging.info("reads " + freads_in)

    # Find best run
    sort_elbo = [
        (idx, state_run[1]["elbo"]) for idx, state_run in enumerate(result_list)
    ]
    # sort list of tuple by ELBO value
    sort_elbo.sort(key=lambda x: x[1], reverse=True)

    best_run_idx = sort_elbo[0][0]
    best_run_elbo = sort_elbo[0][1]
    logging.info("Maximal ELBO " + str(best_run_elbo) + "in run " + str(best_run_idx))

    sorted_results = [result_list[tuple_idx_elbo[0]] for tuple_idx_elbo in sort_elbo]
    exit_meassage = sorted_results[0][1]["exit_message"]
    logging.info("CAVI termination " + str(exit_meassage))

    if record_history:
        with open(output_name + "all_results.pkl", "wb") as f2:
            pickle.dump(sorted_results, f2)

        logging.info(
            "Results dicts of all runs written to " + output_name + "all_results.pkl"
        )

    state_curr_dict = result_list[best_run_idx][0]

    summary = analyze_results.summarize_results(
        state_curr_dict,
        alphabet,
        reads_seq_binary,
        reads_weights,
        read_mapping,
        reads_log_error_proba,
        reference_binary,
    )
    state_curr_dict.update(summary)

    # write output like in original shorah
    analyze_results.haplotypes_to_fasta(state_curr_dict, output_name + "support.fas")
    analyze_results.correct_reads(state_curr_dict, output_name + "cor.fas")

###
### example for testing
###
if __name__ == "__main__":
    main(
        freads_in="w-rep_a-p0-7210-8009.reads.fas",
        fref_in="w-rep_a-p0-7210-8009.ref.fas",
        fname_qualities="w-rep_a-p0-7210-8009.qualities.npy",
        output_dir='./',
        n_starts=1,
        K=100,
        alpha0=0.0001,
        alphabet="ACGT-",
        unique_modus=True,
        convergence_threshold=1e-03,
        record_history=False,
    )
