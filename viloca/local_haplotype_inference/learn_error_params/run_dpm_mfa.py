#!/usr/bin/env python3

import sys
import os
import pickle
import logging
from numpy.random import default_rng

# my python-scripts
from . import preparation
from . import analyze_results
from . import cavi

logging.basicConfig(
    filename="viloca.log", encoding="utf-8", level=logging.INFO
)


def gzip_file(f_name):
    """Gzip a file and return the name of the gzipped, removing the original"""
    import gzip

    f_in = open(f_name, "rb")
    f_out = gzip.open(f_name + ".gz", "wb")
    f_out.writelines(f_in)
    f_out.close()
    f_in.close()
    os.remove(f_in.name)

    return f_out.name


def main(freads_in, fref_in, output_dir, n_starts, K, alpha0, alphabet="ACGT-", unique_modus=True, record_history=False, seed=27):

    window_id = freads_in.split("/")[-1][:-4]  # freads_in is absolute path

    output_name = output_dir + window_id + "-"
    os.makedirs(output_dir, exist_ok=True)

    # Read in reads
    reference_seq, ref_id = preparation.load_reference_seq(fref_in)
    reference_binary = preparation.reference2binary(reference_seq, alphabet)
    reads_list = preparation.load_fasta2reads_list(freads_in, alphabet, unique_modus)
    reads_seq_binary, reads_weights = preparation.reads_list_to_array(reads_list)

    if n_starts >1:
        result_list = cavi.multistart_cavi(
            K,
            alpha0,
            alphabet,
            reference_binary,
            reference_seq,
            reads_list,
            reads_seq_binary,
            reads_weights,
            n_starts,
            output_name,
            seed
        )

    else:
        result_list = [cavi.run_cavi(
            K,
            alpha0,
            alphabet,
            reference_binary,
            reference_seq,
            reads_list,
            reads_seq_binary,
            reads_weights,
            0,
            output_name,
            default_rng(seed=seed),
            seed
        )]

    logging.info("len(result_list) " + str(len(result_list)))
    logging.info("reference " + fref_in)
    logging.info("reads " + freads_in)
    logging.info("lenght of sequences " + str(reads_list[0].seq_binary.shape[0]))
    logging.info("number of reads " + str(len(reads_list)))

    # Find best run
    sort_elbo = [
        (idx, state_run[1]["elbo"]) for idx, state_run in enumerate(result_list)
    ]
    sort_elbo.sort(key=lambda x: x[1], reverse=True)  # sort list of tuple by ELBO value

    max_idx = sort_elbo[0][0]
    max_elbo = sort_elbo[0][1]
    sort_results = [result_list[tuple_idx_elbo[0]] for tuple_idx_elbo in sort_elbo]

    logging.info("CAVI termination " + str(sort_results[0][2]["exit_message"]))

    if record_history:
        with open(output_name + "all_results.pkl", "wb") as f2:
            pickle.dump(sort_results, f2)
        logging.info(
            "Results dicts of all runs written to " + output_name + "all_results.pkl"
        )

    state_curr_dict = result_list[max_idx][1]
    logging.info("Maximal ELBO " + str(max_elbo) + "in run " + str(max_idx))

    # write output like in original shorah
    analyze_results.haplotypes_to_fasta(state_curr_dict, output_name + "support.fas")
    analyze_results.correct_reads(state_curr_dict, output_name + "cor.fas")


if __name__ == "__main__":
    main(
        sys.argv[1],
        sys.argv[2],
        sys.argv[3],
        int(sys.argv[4]),
        int(sys.argv[5]),
        float(sys.argv[6]),
        sys.argv[7],
    )
# freads_in, fref_in, output_dir, n_starts, K, alpha0, alphabet = 'ACGT-'
