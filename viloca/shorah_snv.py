#!/usr/bin/env python3

# Copyright 2007-2018
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
# Lukas Geyrhofer,
# Osvaldo Zagordi,
# Kerensa McElroy,
# ETH Zurich

# This file is part of ShoRAH.
# ShoRAH is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ShoRAH is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ShoRAH.  If not, see <http://www.gnu.org/licenses/>.


"""
    ------------
    Output:
    a file of raw snvs, parsed from the directory support,
    and a directory containing snvs resulting from strand
    bias tests with different sigma values

    Output Files:
    - raw_snv_collected.tsv: contains all SNVs from all windows. Each in a seperate row.
    - raw_snv.tsv: contains all SNVs
    - SNV.tsv: contains only SNVs that are covered by at least {min_windows_coverage}
               windows
    - SNVs_*.tsv, strand bias filter applied
    - raw_snvs_*.tsv: strand bias filter applied
    - SNVs_*_final.vcf: final file
    - SNVs_*_final.csv: file file
    - cooccurring_mutations.csv: contains all mutations linked to the respective haplotype they originate from
    ------------
"""

import glob
import os
import sys
import warnings
import pandas as pd
from collections import namedtuple
from dataclasses import dataclass
import logging
import numpy as np
from Bio import SeqIO
from re import search
from math import log10
import csv
import inspect
from datetime import date

import libshorah

SNP_id = namedtuple("SNP_id", [
    "pos", # a positon in (extended) reference space
    "var" # the variant
])

@dataclass
class SNV:
    chrom: str
    haplotype_id: str
    pos: int
    ref: str
    var: str
    freq: float = 0.0
    support: float = 0.0


standard_header_row = ["Chromosome", "Pos", "Ref", "Var", "Frq", "Pst"]


def _deletion_length(seq, char):
    """Determines the length of the deletion. Note that a sequence migth have
    more than one deletion
    seq: substring of the reconstructed haplotype
    char: character that is used to mark a deletion
    """
    count = 0
    for c in seq:
        if c == char:
            count += 1
        else:
            break
    return count

def _count_double_X(ref, seq, x):
    for i in reversed(range(x + 1)):
        if not ref[i] == seq[i] == "X":
            return x - i
    return 0

def _preprocess_seq_with_X(ref, seq):
    for idx, v in enumerate(seq):
        if v == "-" and ref[idx] == "X":
            seq = seq[:idx] + "X" + seq[idx+1:]
            assert len(seq) == len(ref)
    return seq

def _compare_ref_to_read(ref: str, seq: str, start, snp, av, post, chrom, haplotype_id):
    assert len(ref) == len(seq)

    pos = start
    tot_snv = 0
    aux_del = -1

    change_in_reference_space = 0

    for idx, v in enumerate(ref):  # iterate on the reference

        if v != seq[idx]:  # SNV detected, save it
            assert not (v != "X" and seq[idx] == "X")

            if seq[idx] == "-" or v == "X": # TODO what is if window starts like that?
                char = "-"
                relevant_seq = seq
                secondary_seq = ref
                if v == "X":
                    char = "X"
                    relevant_seq = ref
                    secondary_seq = seq
                # Avoid counting multiple times a long deletion in the same haplotype
                if idx > aux_del:
                    tot_snv += 1
                    # Check for gap characters and get the deletion
                    # length
                    del_len = _deletion_length(relevant_seq[idx:], char)
                    aux_del = idx + del_len
                    snp_id = SNP_id(
                        pos=pos, var=seq[idx:aux_del]
                    )

                    if snp_id in snp:
                        # Aggregate counts for long deletions which
                        # are observed in multiple haplotypes
                        snp[snp_id].freq += av
                        snp[snp_id].support += post * av
                    else:
                        # Report preceeding position as well
                        pos_prev = pos - 1
                        num_double_X = _count_double_X(ref, seq, pos_prev - start)
                        secondary_seq = secondary_seq[pos_prev - start - num_double_X] + secondary_seq[
                            (pos - start) : (pos_prev + del_len - start + 1)
                        ] # TODO pos_prev - 1 - beg might be out of range

                        snp[snp_id] = SNV(
                            chrom,
                            haplotype_id,
                            pos_prev - change_in_reference_space,
                            ref[pos_prev - start - num_double_X] if v =="X" else secondary_seq,
                            secondary_seq if v =="X" else relevant_seq[pos_prev - start - num_double_X],
                            av,
                            post * av,
                        )
            else:
                tot_snv += 1
                snp_id = SNP_id(pos=pos, var=seq[idx])
                if snp_id in snp:
                    snp[snp_id].freq += av
                    snp[snp_id].support += post * av
                else:
                    snp[snp_id] = SNV(
                        chrom,
                        haplotype_id,
                        pos - change_in_reference_space,
                        v,
                        seq[idx],
                        av,
                        post * av
                    )

        if v == "X":
            change_in_reference_space += 1

        pos += 1

    return tot_snv

def parseWindow(line, extended_window_mode, exclude_non_var_pos_threshold,
                working_dir, threshold=0.9):
    """SNVs from individual support files, getSNV will build
    the consensus SNVs
    It returns a dictionary called snp with the following structure
    key:   pos.allele (position on the reference file and mutated base)
    value: reference name, position, reference_base, mutated base,
           average number of reads, posterior times average n of reads
    """

    snp = {}
    reads = 0.0
    # winFile, chrom, beg, end, cov
    _, chrom, beg, end, _ = line.rstrip().split("\t")

    file_stem = "w-%s-%s-%s" % (chrom, beg, end)
    haplo_filename = os.path.join(working_dir, "haplotypes", file_stem + ".reads-support.fas")

    if extended_window_mode:
        ref_name = "extended-ref"
    elif exclude_non_var_pos_threshold > 0:
        ref_name = "envp-full-ref"
    else:
        ref_name = "ref"
    ref_filename = os.path.join("raw_reads", f"{file_stem}.{ref_name}.fas")

    start = int(beg)
    max_snv = -1

    try:
        with open(haplo_filename, "rt") as window, open(ref_filename, "rt") as ref:
            d = dict([[s.id, str(s.seq).upper()] for s in SeqIO.parse(ref, "fasta")])
            refSlice = d[chrom]

            for s in SeqIO.parse(window, "fasta"):
                seq = str(s.seq).upper()
                haplotype_id = str(s.id.split("|")[0]) + "-" + beg + "-" + end
                match_obj = search(r"posterior=(.*)\s*ave_reads=(.*)", s.description)
                post, av = float(match_obj.group(1)), float(match_obj.group(2))

                if post > 1.0:
                    warnings.warn("posterior = %4.3f > 1" % post)
                    logging.warning("posterior = %4.3f > 1" % post)

                # sequences in support file exceeding the posterior threshold
                if post < threshold:
                    continue

                reads += av

                tot_snv = _compare_ref_to_read(
                    refSlice, _preprocess_seq_with_X(refSlice, seq), start, snp, av, post, chrom, haplotype_id
                )

                if tot_snv > max_snv:
                    max_snv = tot_snv
    except FileNotFoundError as not_found:
        # Known bug
        logging.warning(f"{not_found.filename} was not found. All reads might start with N. Skipped.")

    logging.info("max number of snvs per sequence found: %d", max_snv)
    # normalize
    for _, v in snp.items():
        v.support /= v.freq
        v.freq /= reads

    return snp


def add_SNV_to_dict(all_dict, add_key, add_val):
    if add_key in all_dict.keys():
        all_dict[add_key].append(add_val)
    else:
        all_dict.update({add_key: [add_val]})

    return all_dict


def getSNV(extended_window_mode, exclude_non_var_pos_threshold, working_dir, window_thresh=0.9):
    """Parses SNV from all windows and output the dictionary with all the
    information.

    Writes file: collected_snv.tsv

    Note: SNVs that were discovered by several windows will occurr in several
    rows of the file (one row = one window).

    Returns a dict containing all SNVs, already grouped together per SNV_id.
    """

    all_snp = {}
    tmp = []

    # cycle over all windows reported in coverage.txt
    with open("coverage.txt") as cov_file, open(
        os.path.join(working_dir, "raw_snv_collected.tsv"), "w"
    ) as f_collect:
        f_collect.write("\t".join(standard_header_row) + "\n")
        for i in cov_file:
            snp = parseWindow(i, extended_window_mode, exclude_non_var_pos_threshold,
                              working_dir, window_thresh)
            winFile, chrom, beg, end, cov = i.rstrip().split("\t")
            # all_snp = join_snp_dict(all_snp, snp)
            for SNV_id, val in sorted(snp.items()):
                all_snp = add_SNV_to_dict(all_snp, SNV_id, val)
                f_collect.write(
                    "\t".join(
                        map(
                            str,
                            [
                                val.chrom,
                                val.pos,
                                val.ref,
                                val.var,
                                val.freq,
                                val.support,
                            ],
                        )
                    )
                    + "\n"
                )

    return all_snp

def writeRaw(all_snv, min_windows_coverage, working_dir):
    """Write the SNVs that were collected for each window into
    - raw_snv.tsv: contains all SNVs
    - SNV.tsv: contains only SNVs that are covered by at least {min_windows_coverage}
               windows
    """
    header_row = ["Chromosome", "Pos", "Ref", "Var"]

    max_number_window_covering_SNV = np.max(
        [len(val) for _, val in sorted(all_snv.items())]
    )

    header_row = header_row + [
        "Frq" + str(k + 1) for k in range(max_number_window_covering_SNV)
    ]
    header_row = header_row + [
        "Pst" + str(k + 1) for k in range(max_number_window_covering_SNV)
    ]

    with (open(os.path.join(working_dir, "raw_snv.tsv"), "w") as f_raw_snv,
          open(os.path.join(working_dir, "SNV.tsv"), "w") as f_SNV):
        f_raw_snv.write("\t".join(header_row) + "\n")
        f_SNV.write("\t".join(header_row) + "\n")

        for _, val in sorted(all_snv.items()):

            write_line = [val[0].chrom, val[0].pos, val[0].ref, val[0].var]
            freq_list = [single_val.freq for single_val in val]
            support_list = [single_val.support for single_val in val]

            number_window_covering_SNV = len(freq_list)

            if number_window_covering_SNV < max_number_window_covering_SNV:
                filler = (max_number_window_covering_SNV - len(freq_list)) * ["*"]
                freq_list += filler
                support_list += filler

            write_line = write_line + freq_list + support_list

            f_raw_snv.write("\t".join(map(str, write_line)) + "\n")

            if number_window_covering_SNV >= min_windows_coverage:
                f_SNV.write("\t".join(map(str, write_line)) + "\n")


def sb_filter(
    in_bam,
    file_to_append,
    out_file_prefix,
    sigma,
    amplimode="",
    drop_indels="",
    max_coverage=100000,
):  # TODO max_coverage is 10 times higher than in Cpp
    """run strand bias filter calling the external program 'fil'"""

    logging.debug("Running fil")
    logging.debug(f"{in_bam} {file_to_append} {out_file_prefix} {sigma} {max_coverage}")
    retcode = libshorah.fil(
        in_bam,
        file_to_append,
        out_file_prefix,
        sigma,
        max_coverage,
        False if amplimode == "" else True,
        False if drop_indels == "" else True,
    )
    return retcode


# is last column of final output file
def BH(p_vals, n):
    """performs Benjamini Hochberg procedure, returning q-vals'
    you can also see http://bit.ly/QkTflz
    """
    # p_vals contains the p-value and the index where it has been
    # found, necessary to assign the correct q-value
    q_vals_l = []
    prev_bh = 0
    for i, p in enumerate(p_vals):
        # Sometimes this correction can give values greater than 1,
        # so we set those values at 1
        bh = p[0] * n / (i + 1)
        bh = min(bh, 1)
        # To preserve monotonicity in the values, we take the
        # maximum of the previous value or this one, so that we
        # don't yield a value less than the previous.
        bh = max(bh, prev_bh)
        prev_bh = bh
        q_vals_l.append((bh, p[1]))
    return q_vals_l

## write cooccurring_mutations mutations file
def __mutations_in_haplotype(ref: str, seq: str, start, av, post, chrom, haplotype_id):
    """
    This function is adapted from _compare_ref_to_read.
    """

    assert len(ref) == len(seq)

    pos = start
    tot_snv = 0
    aux_del = -1

    tmp = [] # collect dicts of mutations on haplotype

    change_in_reference_space = 0

    for idx, v in enumerate(ref):  # iterate on the reference

        if v != seq[idx]:  # SNV detected, save it
            assert not (v != "X" and seq[idx] == "X")

            if seq[idx] == "-" or v == "X": # TODO what is if window starts like that?
                char = "-"
                relevant_seq = seq
                secondary_seq = ref
                if v == "X":
                    char = "X"
                    relevant_seq = ref
                    secondary_seq = seq
                # Avoid counting multiple times a long deletion in the same haplotype
                if idx > aux_del:
                    tot_snv += 1
                    # Check for gap characters and get the deletion
                    # length
                    del_len = _deletion_length(relevant_seq[idx:], char)
                    aux_del = idx + del_len

                    pos_prev = pos - 1
                    num_double_X = _count_double_X(ref, seq, pos_prev - start)
                    secondary_seq = secondary_seq[pos_prev - start - num_double_X] + secondary_seq[
                        (pos - start) : (pos_prev + del_len - start + 1)
                    ] # TODO pos_prev - 1 - beg might be out of range

                    # add deletion to list
                    snv_dict = {
                            "haplotype_id": haplotype_id,
                            "chrom": chrom,
                            "start": start,
                            "position": pos_prev - change_in_reference_space,
                            "ref": ref[pos_prev - start - num_double_X] if v =="X" else secondary_seq,
                            "var": secondary_seq if v =="X" else relevant_seq[pos_prev - start - num_double_X],
                            "reads": av,
                            "support": post,
                        }

                    tmp.append(snv_dict)
            else:
                tot_snv += 1

                snv_dict = {
                            "haplotype_id": haplotype_id,
                            "chrom": chrom,
                            "start": start,
                            "position": pos - change_in_reference_space,
                            "ref": v,
                            "var": seq[idx],
                            "reads": av,
                            "support": post,
                            #"support_normalize": post * av, # TODO Lara: not only post ??
                        }

                tmp.append(snv_dict)

        if v == "X":
            change_in_reference_space += 1

        pos += 1

    return tmp

def get_cooccuring_muts_haplo_df(haplo_filename, ref_filename, beg, end, chrom):

    reads = 0.0
    start = int(beg)
    max_snv = -1

    tmp_df = []

    with open(haplo_filename, "rt") as window, open(ref_filename, "rt") as ref:
        d = dict([[s.id, str(s.seq).upper()] for s in SeqIO.parse(ref, "fasta")])
        refSlice = d[chrom]

        for s in SeqIO.parse(window, "fasta"):
            seq = str(s.seq).upper()
            haplotype_id = str(s.id.split("|")[0]) + "-" + beg + "-" + end
            match_obj = search(r"posterior=(.*)\s*ave_reads=(.*)", s.description)
            post, av = float(match_obj.group(1)), float(match_obj.group(2))
            reads += av

            tot_snv = __mutations_in_haplotype(
                refSlice, seq, start, av, post, chrom, haplotype_id
            )
            if tot_snv == []:
                # haplotype is reference
                tot_snv = [{
                            "haplotype_id": "reference"+"-" + beg + "-" + end,
                            "chrom": chrom,
                            "start": start,
                            "reads": av,
                            "support": post,
                        }]

            tmp_df.append(pd.DataFrame(tot_snv))

    if len(tmp_df)>0:
        df = pd.concat(tmp_df)
        df['coverage']= reads
    else:
        df = pd.DataFrame(columns=['coverage'])

    return df

def write_cooccuring_muts_file():

    tmp_df =[]

    # iterate over haplotype files
    with open("coverage.txt") as cov_file:
        for line in cov_file:
            # winFile, chrom, beg, end, cov
            _, chrom, beg, end, _ = line.rstrip().split("\t")

            file_stem = "w-%s-%s-%s" % (chrom, beg, end)
            haplo_filename = os.path.join(working_dir, "haplotypes", file_stem + ".reads-support.fas")
            ref_name = os.path.join(working_dir, "haplotypes", file_stem + ".ref.fas")
            tmp_df.append(get_cooccuring_muts_df(fname_haplo, fname_ref, beg,end,chrom))

    pd.concat(tmp_df).to_csv("cooccurring_mutations.csv")

def main(args):
    """main code"""

    reference = args.f
    bam_file = args.b
    sigma = args.sigma
    increment = args.increment
    max_coverage = args.max_coverage
    ignore_indels = args.ignore_indels
    posterior_thresh = args.posterior_thresh
    path_insert_file = args.path_insert_file
    extended_window_mode = args.extended_window_mode
    min_windows_coverage = args.min_windows_coverage
    working_dir = args.working_dir
    exclude_non_var_pos_threshold = args.exclude_non_var_pos_threshold
    strand_bias_filter = args.strand_bias_filter

    assert os.path.isdir(args.working_dir) or args.working_dir == ""

    logging.info(str(inspect.getfullargspec(main)))

    # snpD_m is the file with the 'consensus' SNVs (from different windows)
    logging.debug("now parsing SNVs")
    all_SNVs = getSNV(extended_window_mode, exclude_non_var_pos_threshold,
                      working_dir, posterior_thresh)

    if path_insert_file is not None:
        min_windows_coverage = 1

    writeRaw(all_SNVs, min_windows_coverage, working_dir)

    with open(os.path.join(working_dir, "raw_snv.tsv")) as f_raw_snv:
        windows_header_row = f_raw_snv.readline().split("\t")
        windows_header_row[-1] = windows_header_row[-1].split("\n")[0]

    d = " -d" if ignore_indels else ""

    a = " -a" if increment == 1 else ""  # TODO when is increment == 1 (amplimode)

    # run strand bias filter
    retcode_n = sb_filter(
        bam_file,
        os.path.join(working_dir, "SNV.tsv"),
        os.path.join(working_dir, "SNVs_"),
        sigma,
        amplimode=a,
        drop_indels=d,
        max_coverage=max_coverage,
    )

    if retcode_n != 0:
        logging.error('sb_filter exited with error %d', retcode_n)
        sys.exit()

    snpFile = glob.glob(os.path.join(working_dir, "SNVs_*.tsv"))[0] # takes the first file only
    logging.debug(f"For BH - selected file: {snpFile}")

    write_list = []
    d = {}

    with open(snpFile) as f:
        for line_no, line in enumerate(f):
            parts = line.rstrip().split("\t")
            write_list.append(parts)
            idx = parts[0] + parts[1] + parts[2] + parts[3]
            if idx in d:
                d[idx][1].append(line_no)
            else:
                d[idx] = (float(parts[-1]), [line_no])

    p_vals = list(d.values())

    # sort p values, correct with Benjamini Hochberg and append to output
    p_vals.sort()
    q_vals = BH(p_vals, len(p_vals))

    for q, indices in q_vals:
        for i in indices:
            write_list[i].append(q)

    # Write ShoRAH csv output file
    if "csv" in args.format:
        csv_file = "{}_final.csv".format(os.path.splitext(snpFile)[0])
        header_row = windows_header_row + [
            "Fvar",
            "Rvar",
            "Ftot",
            "Rtot",
            "Pval",
            "Qval",
        ]
        with open(csv_file, "w") as cf:
            writer = csv.writer(cf)
            writer.writerow(header_row)
            # only print when q >= 5%
            for wl in write_list:
                if (wl[-1] >= 0.05) & (strand_bias_filter==True):
                    writer.writerow(wl)
                elif (strand_bias_filter==False):
                    writer.writerow(wl)

    max_number_window = int(windows_header_row[-1].split("Pst")[1])

    if "vcf" in args.format:
        VCF_file = f"{os.path.splitext(snpFile)[0]}_final.vcf"
        VCF_meta = [
            "##fileformat=VCFv4.2",
            f"##fileDate={date.today():%Y%m%d}",
            f"##source=VILOCA",
            f"##reference={args.f}",
        ]
        ref_m = dict([[s.id, str(s.seq).upper()] for s in SeqIO.parse(reference, "fasta")])
        for ref_name, ref_seq in ref_m.items(): # TODO can be removed? why?
            VCF_meta.append(
                f"##contig=<ID={ref_name},length={len(ref_seq)}>",
            )

        VCF_meta.extend(
            [
                '##INFO=<ID=Fvar,Number=1,Type=Integer,Description="Number of forward reads with variant">',
                '##INFO=<ID=Rvar,Number=1,Type=Integer,Description="Number of reverse reads with variant">',
                '##INFO=<ID=Ftot,Number=1,Type=Integer,Description="Total number of forward reads">',
                '##INFO=<ID=Rtot,Number=1,Type=Integer,Description="Total number of reverse reads">',
                '##INFO=<ID=Pval,Number=1,Type=Float,Description="P-value for strand bias">',
                '##INFO=<ID=Qval,Number=1,Type=Float,Description="Q-value for strand bias">',
            ]
        )
        VCF_meta.extend(
            [
                '##INFO=<ID=FreqX,Number=1,Type=Float,Description="Frequency of the variant in window X">',
                '##INFO=<ID=PostX,Number=1,Type=Float,Description="Posterior probability of the variant in window x">',
            ]
        )

        with open(VCF_file, "w") as vcf:
            vcf.write("\n".join(VCF_meta))
            # VCFv4.2 HEADER line
            vcf.write("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
            # Iterate over single SNV lines and write them to output file
            for wl in write_list:
                chrom = wl[0]
                pos = wl[1]
                base_ref = wl[2]
                base_var = wl[3]

                freqs_vcf = wl[4 : 4 + max_number_window]
                posts_vcf = wl[
                    4 + max_number_window : 4 + max_number_window + max_number_window
                ]

                # TODO: Why is the sampler returning mutations calls with freq==0 in all windows, is that a problem of the model?
                # TODO: not in csv
                sum_Freq = 0
                for freq in freqs_vcf:
                    try:
                        sum_Freq += float(freq)
                    except:
                        sum_Freq += 0

                # only print when q >= 5%
                q = wl[-1]
                if strand_bias_filter==True:
                    q = wl[-1]
                else:
                    q=1.0
                if (q >= 0.05) and (sum_Freq > 0):
                    info = (
                        f"Fvar={wl[-6]};Rvar={wl[-5]};Ftot={wl[-4]};"
                        f"Rtot={wl[-3]};Pval={wl[-2]};Qval={wl[-1]}"
                    )

                    freq_str = ";".join(
                        [f"Freq{i+1}={j}" for i, j in enumerate(freqs_vcf) if j != "*"]
                    )

                    post_str = ";".join(
                        [f"Post{i+1}={j}" for i, j in enumerate(posts_vcf) if j != "*"]
                    )

                    info = f"{freq_str};{post_str};{info}".replace("-", "0")

                    post_all = []
                    for freq, post in zip(freqs_vcf, posts_vcf):
                        if freq == "*":
                            pass
                        elif freq == "-":
                            post_all.append(0)
                        else:
                            post_all.append(min([1, float(post)]))
                        # Calculate posterior average
                        post_avg = sum(post_all) / len(post_all)
                    # Calculate a Phred quality score where the base calling
                    # error probabilities is set to (1 - posterior avg).
                    # Maximum is set to 100.
                    try:
                        qual_norm = -10 * log10(1 - post_avg)
                    except ValueError:
                        qual_norm = 100

                    vcf.write(
                        f"\n{wl[0]}\t{wl[1]}\t.\t{wl[2]}\t{wl[3]}"
                        f"\t{qual_norm}\tPASS\t{info}"
                    )
