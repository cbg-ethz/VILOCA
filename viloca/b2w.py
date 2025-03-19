import pysam
from typing import Optional, List, Tuple, Union
from viloca.tiling import TilingStrategy, EquispacedTilingStrategy
import numpy as np
import math
import logging
from multiprocessing import Pool, cpu_count
import os
import hashlib
import time
from collections import defaultdict
from pathlib import Path
#from viloca_rust import run_one_window_rust

def _write_to_file(lines, file_name):
    with open(file_name, "w") as f:
        f.writelines("%s\n" % l for l in lines)


def _calc_via_pileup(samfile, reference_name, maximum_reads):
    budget = dict()
    max_ins_at_pos = dict()
    indel_map = set() # TODO quick fix because pileup created duplicates; why?

    for pileupcolumn in samfile.pileup(
        reference_name,
        max_depth=1_000_000, # TODO big enough?
        stepper="nofilter",
        multiple_iterators=False, # TODO true makes faster?
        ignore_overlaps=False,
        adjust_capq_threshold=0,
        ignore_orphans=False,
        min_base_quality=0):

        budget[pileupcolumn.reference_pos] = min(
            pileupcolumn.nsegments,
            maximum_reads-1 # minus 1 because because maximum_reads is exclusive
        )

        max_at_this_pos = 0
        for pileupread in pileupcolumn.pileups:

            #assert pileupread.indel >= 0, "Pileup read indel is negative" TODO

            if pileupread.indel > 0 or pileupread.is_del:
                indel_map.add((
                    pileupread.alignment.query_name, # TODO is unique?
                    pileupread.alignment.reference_start, # TODO is unique?
                    hashlib.md5(pileupread.alignment.cigarstring.encode()).hexdigest(), # change to this hast to be compatible with rust
                    #hash(pileupread.alignment.cigarstring), # TODO is unique?
                    pileupcolumn.reference_pos,
                    pileupread.indel,
                    pileupread.is_del
                ))


            if pileupread.indel > max_at_this_pos:
                max_at_this_pos = pileupread.indel

        if max_at_this_pos > 0:
            max_ins_at_pos[pileupcolumn.reference_pos] = max_at_this_pos

    # ascending reference_pos are necessary for later steps
    indel_map = sorted(indel_map, key=lambda tup: tup[3])

    return budget, indel_map, max_ins_at_pos


def _build_one_full_read_no_extended_window_mode(
    full_read: List[str],
    full_qualities: Union[List[int], List[str], None],
    read_query_name: str | None,
    full_read_cigar_hash: str | None,
    first_aligned_pos: int,
    last_aligned_pos: int,
    indel_map: List[Tuple[str, int, str, int, int, int]],
    insertion_char: str
) -> Tuple[str, List[int]]:
    """
    Builds a full read sequence without extended window mode, handling indels.
    Args:
        full_read: List of characters representing the read sequence.
        full_qualities: Optional list of quality scores (ints or strs).
        read_query_name: Name of the read (or None).
        full_read_cigar_hash: Hash of the CIGAR string (or None).
        first_aligned_pos: First aligned position.
        last_aligned_pos: Last aligned position (unused).
        indel_map: List of tuples (name, start, cigar_hash, ref_pos, indel_len, is_del).
        insertion_char: Character to insert for deletions ("-" or "X").
    Returns:
        Tuple of (full read string, quality scores as list of ints).
    """
    # Import inside function to avoid circular dependency
    from viloca_rust import build_one_full_read_no_extended_window_mode_rust

    assert insertion_char in ["-", "X"], "Illegal char representation insertion"

    # Convert full_qualities to strings if present, None otherwise
    qualities_in = [str(q) for q in full_qualities] if full_qualities is not None else None

    # Call Rust function
    full_read_out, full_qualities_out = build_one_full_read_no_extended_window_mode_rust(
        first_aligned_pos,        # Position 1 (matches Rust)
        last_aligned_pos,         # Position 2 (matches Rust)
        full_read,                # Position 3 (matches Rust)
        indel_map,                # Position 4 (matches Rust)
        insertion_char,           # Position 5 (matches Rust)
        qualities_in,             # Position 6 (matches Rust)
        read_query_name,          # Position 7 (matches Rust)
        full_read_cigar_hash      # Position 8 (matches Rust)
    )

    # Convert qualities back to ints if present
    full_qualities_out = (
        [int(q) for q in full_qualities_out] if full_qualities_out is not None else None
    )

    return full_read_out, full_qualities_out


def _build_one_full_read_with_extended_window_mode(
    full_read: list[str],
    full_qualities: list[int] | list[str],
    read_query_name: str | None,
    full_read_cigar_hash: str | None,
    first_aligned_pos: int,
    last_aligned_pos: int,
    indel_map: list[tuple[str, int, str, int, int, int]],
    max_ins_at_pos: dict[int, int],
    insertion_char: str
) -> tuple[str, list[int]]:
    """Wrapper for Rust implementation of extended window mode read building."""
    # Import inside function to avoid circular dependency
    from viloca_rust import build_one_full_read_with_extended_window_mode_rust

    # Convert qualities to strings if present
    qualities_in = [str(q) for q in full_qualities] if full_qualities is not None else None

    # Call Rust function
    full_read_out, full_qualities_out = build_one_full_read_with_extended_window_mode_rust(
        first_aligned_pos,
        last_aligned_pos,
        full_read,
        indel_map,
        max_ins_at_pos,
        insertion_char,
        qualities_in,
        read_query_name,
        full_read_cigar_hash
    )

    # Convert qualities back to ints if present
    full_qualities_out = (
        [int(q) for q in full_qualities_out] if full_qualities_out is not None else None
    )

    return full_read_out, full_qualities_out

"""
def _run_one_window_rust(
    samfile,
    window_start,
    reference_name,
    window_length,
    control_window_length,
    minimum_overlap,
    permitted_reads_per_location,
    counter,
    exact_conformance_fix_0_1_basing_in_reads,
    indel_map,
    max_ins_at_pos,
    extended_window_mode,
    exclude_non_var_pos_threshold
):

    iter_reads = samfile.fetch(reference_name, window_start, window_start + window_length)

    # Extract read data into simple Python structures
    extracted_reads = []
    for read in iter_reads:
        if read.reference_start is None or read.reference_end is None:
            continue

        qname = read.query_name
        ref_start = read.reference_start
        ref_end = read.reference_end  # pysam reference_end is exclusive
        seq = read.query_sequence
        quals = list(read.query_qualities) if read.query_qualities else None
        cigarstring = read.cigarstring

        extracted_reads.append((qname, ref_start, ref_end, seq, quals, cigarstring))

    # Call the updated Rust function with extracted reads
    arr, arr_read_qualities_summary, arr_read_summary, pos_filter = run_one_window_rust(
        reads=extracted_reads,
        reference_name=reference_name,
        window_start=window_start,
        window_length=window_length,
        control_window_length=control_window_length,
        minimum_overlap=int(minimum_overlap),
        permitted_reads_per_location=dict(permitted_reads_per_location),
        exact_conformance_fix_0_1_basing_in_reads=exact_conformance_fix_0_1_basing_in_reads,
        indel_map=indel_map,
        max_ins_at_pos=max_ins_at_pos,
        extended_window_mode=extended_window_mode,
        exclude_non_var_pos_threshold=exclude_non_var_pos_threshold,
        counter=counter,
    )

    converted_arr = []
    for entry in arr:
        header, seq = entry.split("\n")
        name, pos = header[1:].split()
        converted_arr.append(f'>{name} {pos}\n{"".join(seq)}')

    # Convert qualities summary to numpy arrays if not None
    if arr_read_qualities_summary is not None:
        arr_read_qualities_summary = [np.array(q) for q in arr_read_qualities_summary]

    counter = window_start + window_length

    return converted_arr, arr_read_qualities_summary, arr_read_summary, pos_filter
"""

def update_tiling(tiling, extended_window_mode, max_ins_at_pos):
    """
    input tiling:
    window_start is 1-based
    max_ins_at_pos is 0-based

    return: tiling = [
            (window_start, original_window_length, control_window_length)
            for each window
            ]
    """
    updated_tiling = []

    for idx, (window_start, window_length) in enumerate(tiling):
        original_window_length = window_length
        if extended_window_mode:
            for pos, val in max_ins_at_pos.items():
                if window_start - 1 <= pos < window_start - 1 + original_window_length:
                    window_length += val
            updated_tiling.append((window_start,original_window_length, window_length))
        else:
            updated_tiling.append((window_start,original_window_length, window_length))

    return updated_tiling



def _run_one_window(samfile, window_start, reference_name, window_length,control_window_length,
        minimum_overlap, permitted_reads_per_location, counter,
        exact_conformance_fix_0_1_basing_in_reads, indel_map, max_ins_at_pos,
        extended_window_mode, exclude_non_var_pos_threshold):

    arr = []
    arr_read_summary = []
    arr_read_qualities_summary = []

    iter = samfile.fetch(
        reference_name,
        window_start, # 0 based
        window_start + window_length # arg exclusive as per pysam convention
    )

    # TODO: minimum_overlap
    # TODO: original_window_length
    # TODO: window_length
    original_window_length = window_length
    window_length = control_window_length
    original_minimum_overlap = minimum_overlap
    if extended_window_mode:
        # this is now done intilaly for all windows
        #for pos, val in max_ins_at_pos.items():
        #    if window_start <= pos < window_start + original_window_length:
        #        window_length += val
        minimum_overlap *= window_length/original_window_length

    if exclude_non_var_pos_threshold > 0:
        alphabet = "ACGT-NX"
        base_pair_distr_in_window = np.zeros((window_length, len(alphabet)), dtype=int)

    for read in iter:
        if (read.reference_start is None) or (read.reference_end is None):
            continue
        first_aligned_pos = read.reference_start # this is 0-based
        last_aligned_pos = read.reference_end - 1 #reference_end is exclusive


        if permitted_reads_per_location[first_aligned_pos] == 0:
            continue
        else:
            permitted_reads_per_location[first_aligned_pos] -= 1

        full_read = list(read.query_sequence)
        full_qualities = list(read.query_qualities) if read.query_qualities is not None else None

        for ct_idx, ct in enumerate(read.cigartuples):
            # 0 = BAM_CMATCH, 1 = BAM_CINS, 2 = BAM_CDEL, 7 = BAM_CEQUAL, 8 = BAM_CDIFF
            if ct[0] in [0,1,2,7,8]:
                pass
            elif ct[0] == 4: # 4 = BAM_CSOFT_CLIP
                for _ in range(ct[1]):
                    k = 0 if ct_idx == 0 else len(full_read)-1
                    full_read.pop(k)
                    if full_qualities is not None:
                        full_qualities.pop(k)

                if ct_idx != 0 and ct_idx != len(read.cigartuples)-1:
                    raise ValueError("Soft clipping only possible on the edges of a read.")
            elif ct[0] == 5: # 5 = BAM_CHARD_CLIP
                #logging.debug(f"[b2w] Hard clipping detected in {read.query_name}") # commented out because this happens too often
                pass
            else:
                raise NotImplementedError("CIGAR op code found that is not implemented:", ct[0])

        #full_read, full_qualities = _build_one_full_read(full_read, full_qualities,
        #    read.query_name, hashlib.sha1(read.cigarstring.encode()).hexdigest(), first_aligned_pos, last_aligned_pos,
        #    indel_map, max_ins_at_pos, extended_window_mode, "-")
        if extended_window_mode:
            full_read, full_qualities = _build_one_full_read_with_extended_window_mode(
                full_read, full_qualities, read.query_name,
                hashlib.md5(read.cigarstring.encode()).hexdigest(),
                first_aligned_pos, last_aligned_pos, indel_map, max_ins_at_pos, "-"
            )
        else:
            full_read, full_qualities = _build_one_full_read_no_extended_window_mode(
                full_read, full_qualities, read.query_name,
                hashlib.md5(read.cigarstring.encode()).hexdigest(),
                first_aligned_pos, last_aligned_pos, indel_map, "-"
            )

        if (first_aligned_pos + minimum_overlap < window_start + 1 + window_length
                and last_aligned_pos >= window_start + original_minimum_overlap - 2 # TODO justify 2
                #last_aligned_pos: is in the orginal reference genome space (not in the extended_window_mode-space)
                and len(full_read) >= minimum_overlap):

            num_inserts_right_of_read = 0
            num_inserts_left_of_read = 0
            num_inserts_left_of_window = 0
            if extended_window_mode:
                for pos, val in max_ins_at_pos.items():
                    if last_aligned_pos < pos < window_start + original_window_length:
                        num_inserts_right_of_read += val
                    if window_start <= pos < first_aligned_pos:
                        num_inserts_left_of_read += val
                    if first_aligned_pos <= pos < window_start:
                        num_inserts_left_of_window += val

            start_cut_out = window_start + num_inserts_left_of_window - first_aligned_pos
            end_cut_out = start_cut_out + window_length - num_inserts_left_of_read
            s = slice(max(0, start_cut_out), end_cut_out)

            cut_out_read = full_read[s]
            if full_qualities is None:
                #logging.warning("[b2w] No sequencing quality scores provided in alignment file. Run --sampler learn_error_params.")
                cut_out_qualities = None
            else:
                cut_out_qualities = full_qualities[s]

            k = (window_start + original_window_length - 1 - last_aligned_pos
                + num_inserts_right_of_read)

            if k > 0:
                cut_out_read = cut_out_read + k * "N"
                if full_qualities is not None:
                    cut_out_qualities = cut_out_qualities + k * [2]
                    # Phred scores have a minimal value of 2, where an “N” is inserted
                    # https://www.zymoresearch.com/blogs/blog/what-are-phred-scores
            if start_cut_out < 0:
                cut_out_read = (-start_cut_out + num_inserts_left_of_read) * "N" + cut_out_read
                if full_qualities is not None:
                    cut_out_qualities = (-start_cut_out + num_inserts_left_of_read) * [2] + cut_out_qualities


            assert len(cut_out_read) == window_length, (
                "read unequal window size",
                read.query_name, first_aligned_pos, cut_out_read, window_start, window_length, read.reference_end, len(cut_out_read)
            )
            if cut_out_qualities is not None:
                assert len(cut_out_qualities) == window_length, (
                    "quality unequal window size"
                )

            if exclude_non_var_pos_threshold > 0:
                for idx, letter in enumerate(cut_out_read):
                    base_pair_distr_in_window[idx][alphabet.index(letter)] += 1

            c = 0 if exact_conformance_fix_0_1_basing_in_reads == False else 1
            arr.append([read.query_name, first_aligned_pos + c, np.array(list(cut_out_read))])

            if cut_out_qualities is None:
                arr_read_qualities_summary = None
            else:
                arr_read_qualities_summary.append(np.asarray(cut_out_qualities))

        if read.reference_start >= counter and len(full_read) >= minimum_overlap:
            arr_read_summary.append(
                (read.query_name, read.reference_start + 1, read.reference_end, full_read)
            )

    pos_filter = None
    if exclude_non_var_pos_threshold > 0 and len(arr) > 0:
        pos_filter = (1 - np.amax(base_pair_distr_in_window, axis=1) / np.sum(base_pair_distr_in_window, axis=1)
            >= exclude_non_var_pos_threshold)
        logging.debug(f"Number of positions removed: {len(pos_filter) - pos_filter.sum()}")

        for idx, (_, _, rd) in enumerate(arr):
            arr[idx][2] = rd[pos_filter]
            arr_read_qualities_summary[idx] = arr_read_qualities_summary[idx][pos_filter] # TODO works?

    counter = window_start + window_length

    # TODO move out of this function
    convert_to_printed_fmt = lambda x: [f'>{k[0]} {k[1]}\n{"".join(k[2])}' for k in x]

#    return convert_to_printed_fmt(arr), arr_read_qualities_summary, arr_read_summary, counter, window_length, pos_filter
    return convert_to_printed_fmt(arr), arr_read_qualities_summary, arr_read_summary, pos_filter


def parallel_run_one_window(
    reference_filename,
    minimum_reads,
    tiling,
    region_end,
    idx,
    window_start,
    window_length,
    control_window_length,
    alignment_file,
    reference_name,
    win_min_ext,
    permitted_reads_per_location,
    counter,
    exact_conformance_fix_0_1_basing_in_reads,
    indel_map,
    max_ins_at_pos,
    extended_window_mode,
    exclude_non_var_pos_threshold,
):
    """
    build one window.
    """
    reffile = pysam.FastaFile(reference_filename)

    samfile = pysam.AlignmentFile(
        alignment_file,
        "r", # auto-detect bam/cram (rc)
        reference_filename=reference_filename,
        threads=1
    )

    reads = open(f"reads_{idx}.fas", "w")

    logging.info(f"Working on window (1-based) @ {window_start+1}")

#    (arr, arr_read_qualities_summary, arr_read_summary,
#    counter, control_window_length, pos_filter) = _run_one_window(
    (arr, arr_read_qualities_summary, arr_read_summary,
     pos_filter) = _run_one_window(
        samfile,
        window_start - 1, # make 0 based
        reference_name,
        window_length,
        control_window_length,
        math.floor(win_min_ext * window_length),
        dict(permitted_reads_per_location), # copys dict ("pass by value")
        counter,
        exact_conformance_fix_0_1_basing_in_reads,
        indel_map,
        max_ins_at_pos,
        extended_window_mode,
        exclude_non_var_pos_threshold
    )

    logging.debug(f"Window length: {control_window_length}")

    window_end = window_start + window_length - 1
    file_name = f'w-{reference_name}-{window_start}-{window_end}'

    # TODO solution for backward conformance
    if len(tiling) > 1:
        end_extended_by_a_window = region_end + (tiling[1][0]-tiling[0][0])*3
    else:
        end_extended_by_a_window = region_end + window_length*3

    for read in arr_read_summary:
        if idx == len(tiling) - 1 and read[1] > end_extended_by_a_window:
            continue
        # TODO reads.fas not FASTA conform, +-0/1 mixed
        # TODO global end does not really make sense, only for conformance
        # read name, global start, global end, read start, read end, read
        reads.write(
            f'{read[0]}\t{tiling[0][0]-1}\t{end_extended_by_a_window}\t{read[1]}\t{read[2]}\t{read[3]}\n'
        )

    reads.close()

    if (idx != len(tiling) - 1 # except last
        and len(arr) > 0) or len(tiling) == 1: # suppress output if window empty

        _write_to_file(arr, file_name + '.reads.fas')
        if arr_read_qualities_summary is not None:
            with open(file_name + '.qualities.npy', 'wb') as f:
                np.save(f, np.asarray(arr_read_qualities_summary, dtype=np.int64), allow_pickle=True)

        ref = reffile.fetch(reference=reference_name, start=window_start-1, end=window_end)

        if extended_window_mode:
            for file_name_comp, char in [("extended-ref", "X"), ("ref", "-")]:
                res_ref = _build_one_full_read(
                    list(ref), list(ref), None, None,
                    window_start-1, window_end-1,
                    indel_map, max_ins_at_pos, extended_window_mode, char
                )[0]

                k = max(0, control_window_length - len(res_ref))
                res_ref += k * "N"
                assert_condition = control_window_length == len(res_ref)

                if exclude_non_var_pos_threshold > 0 and file_name_comp == "ref":
                    _write_to_file([
                        f'>{reference_name} {window_start}\n' + res_ref
                    ], file_name + '.envp-full-ref.fas')

                    envp_ref = np.array(list(res_ref))
                    envp_ref[~pos_filter] = "="
                    _write_to_file([
                        f'>{reference_name} {window_start}\n' + "".join(envp_ref)
                    ], file_name + '.envp-ref.fas')

                    reduced_ref = np.array(list(res_ref))[pos_filter]
                    res_ref = "".join(reduced_ref)
                    assert_condition = (control_window_length ==
                                        len(reduced_ref) + len(pos_filter) - pos_filter.sum())

                _write_to_file([
                    f'>{reference_name} {window_start}\n' + res_ref
                ], f'{file_name}.{file_name_comp}.fas')

                assert assert_condition, (
                    f"""
                        Reference ({file_name_comp}) does not have same length as the window.
                        Location: {file_name}
                        Ref: {len(res_ref)}
                        Win: {control_window_length}
                    """
                )

        else:
            k = max(0, control_window_length - len(ref))
            ref += k * "N"

            if exclude_non_var_pos_threshold > 0:
                full_file_name = file_name + '.envp-full-ref.fas'
            else:
                full_file_name = file_name + '.ref.fas'

            _write_to_file([
                f'>{reference_name} {window_start}\n' + ref
            ], full_file_name)

            assert control_window_length == len(ref), (
                f"""
                    Reference does not have same length as the window.
                    Location: {file_name}
                    Ref: {len(ref)}
                    Win: {control_window_length}
                """
            )

            if exclude_non_var_pos_threshold > 0:
                envp_ref = np.array(list(ref))
                envp_ref[~pos_filter] = "="
                _write_to_file([
                    f'>{reference_name} {window_start}\n' + "".join(envp_ref)
                ], file_name + '.envp-ref.fas')
                reduced_ref = np.array(list(ref))[pos_filter]
                _write_to_file([
                    f'>{reference_name} {window_start}\n' + "".join(reduced_ref)
                ], file_name + '.ref.fas')

                assert (control_window_length == len(envp_ref) and
                        control_window_length == len(reduced_ref) + len(pos_filter) - pos_filter.sum()), (
                    f"""
                        Reference does not have same length as the window.
                        Location: {file_name}
                        Envp Ref: {len(envp_ref)}
                        Ref: {len(reduced_ref)}
                        Win: {control_window_length}
                    """
                )

        if len(arr) > minimum_reads:
            line = (
                f'{file_name}.reads.fas\t{reference_name}\t{window_start}\t'
                f'{window_end}\t{len(arr)}'
            )
            _write_to_file([line], f"coverage_{idx}.txt")


def run_window_wrapper(args):
    return parallel_run_one_window(*args)


def build_windows(alignment_file: str, tiling_strategy: TilingStrategy,
    win_min_ext: float, maximum_reads: int, minimum_reads: int,
    reference_filename: str,
    exact_conformance_fix_0_1_basing_in_reads: Optional[bool] = False,
    extended_window_mode: Optional[bool] = False,
    exclude_non_var_pos_threshold: Optional[float] = -1,
    maxthreads: Optional[int] = 1,
    reuse_files = False) -> None:
    """Summarizes reads aligned to reference into windows.
    Three products are created:
    #. Multiple FASTA files (one for each window position)
    #. A coverage file that lists all files in (1)
    #. A FASTA file that lists all reads used in (1)
        .. caution::
            ``reads.fas`` does not comply with the FASTA format.
    Args:
        alignment_file: Path to the alignment file in CRAM format.
        tiling_strategy: A strategy on how the genome is partitioned.
        win_min_ext: Minimum percentage of bases to overlap between reference
            and read to be considered in a window. The rest (i.e.
            non-overlapping part) will be filled with Ns.
        maximum_reads: Upper (exclusive) limit of reads allowed to start at the
            same position in the reference genome. Serves to reduce
            computational load.
        minimum_reads: Lower (exclusive) limit of reads allowed in a window.
            Serves to omit windows with low coverage.
        reference_filename: Path to a FASTA file of the reference sequence.
        exact_conformance_fix_0_1_basing_in_reads: Fixes an incorrect 0-basing
            of reads in the window file in the old C++ version. 1-basing is
            applied everywhere now. Set this flag to `False` only for exact
            conformance with the old version (in tests).
        extended_window_mode: Mode where inserts are not deleted but kept. The
            windows are instead extended with artificial/fake deletions.
        exclude_non_var_pos_threshold: Percentage value below which non-variable
            positions are excluded. Set to -1 if no position should be excluded
            (the default).
    """
    max_proc = min(max(cpu_count() - 1, 1), maxthreads)
    logging.info('CPU(s) count %u,  will run %u build_windows', cpu_count(), max_proc)


    assert 0 <= win_min_ext <= 1
    assert (0 <= exclude_non_var_pos_threshold <= 1 or
            exclude_non_var_pos_threshold == -1)

    pysam.index(alignment_file)
    samfile = pysam.AlignmentFile(
        alignment_file,
        "r", # auto-detect bam/cram (rc)
        reference_filename=reference_filename,
        threads=max_proc #1
    )
    # print this message since pysam print confusing error message with the functions above.
    print("Index file was created successfully.")
    #reffile = pysam.FastaFile(reference_filename) --> we need to read it in each child processo
    #counter = 0 #--> counter is now coputed initially for all windows
    reference_name = tiling_strategy.get_reference_name()
    tiling = tiling_strategy.get_window_tilings()
    region_end = tiling_strategy.get_region_end()

    permitted_reads_per_location, indel_map, max_ins_at_pos = _calc_via_pileup(
        samfile,
        reference_name,
        maximum_reads
    )

    tiling = update_tiling(tiling, extended_window_mode, max_ins_at_pos)
    #print("indel_map", indel_map)

    # generate counter for each window
    # counter = window_start - 1 + control_window_length, # make 0 based
    counter_list = [0] + [window_start - 1 + control_window_length for (window_start, window_length, control_window_length) in tiling]

    process_args = []
    for idx, (window_start, window_length, control_window_length) in enumerate(tiling):
        counter = counter_list[idx]
        if not (reuse_files and os.path.isfile(f"coverage_{idx}.txt")):
            args = (
                reference_filename,
                minimum_reads,
                tiling,
                region_end,
                idx,
                window_start,
                window_length,
                control_window_length,
                alignment_file,
                reference_name,
                win_min_ext,
                permitted_reads_per_location,
                counter,
                exact_conformance_fix_0_1_basing_in_reads,
                indel_map,
                max_ins_at_pos,
                extended_window_mode,
                exclude_non_var_pos_threshold
            )
            process_args.append(args)
        else:
            logging.info(f'[file already exists] Using window files generated on {time.ctime(Path(f"coverage_{idx}.txt").stat().st_mtime)}')


    max_proc = min(max(cpu_count() - 1, 1), maxthreads)
    logging.info('CPU(s) count %u, will run %u build_windows', cpu_count(), max_proc)

    with Pool(processes=max_proc) as pool:
        results = pool.map(run_window_wrapper, process_args)

    """
    for p in all_processes:
      p.start()

    for p in all_processes:
        p.join()
        if p.exitcode != 0:
            logging.debug("[b2w] A process was killed. Terminating the program.")
            exit(1)
    """
    logging.debug("[b2w] All processes completed successfully.")

    samfile.close()

    with open('reads.fas', 'w') as output_file:
        for file_path in [f"reads_{idx}.fas" for idx in range(len(tiling))]:
            with open(file_path, 'r') as input_file:
                for line in input_file:
                    output_file.write(line)
            os.remove(file_path)

    cov_arr = []
    with open("coverage.txt", 'w') as output_file:
        for file_path in [f"coverage_{idx}.txt" for idx in range(len(tiling))]:
            try:
                with open(file_path, 'r') as input_file:
                    for line in input_file:
                        output_file.write(line)
                os.remove(file_path)
            except Exception:
                pass

if __name__ == "__main__":
    import argparse

    # Naming as in original C++ version
    parser = argparse.ArgumentParser(description='b2w')
    parser.add_argument('-w', '--window_length', nargs=1, type=int,
        help='window length', required=True)
    parser.add_argument('-i', '--incr', nargs=1, type=int, help='increment',
        required=True)
    parser.add_argument('-m', nargs=1, type=float, help='minimum overlap in percent',
        required=True)
    parser.add_argument('-x', nargs=1, type=int,
        help='max reads starting at a position', required=True)
    parser.add_argument('-c', nargs=1, type=int,
        help='coverage threshold. Omit windows with low coverage.',
        required=True)

    parser.add_argument('-d', nargs='?',
        help='drop SNVs that are adjacent to insertions/deletions (alternate behaviour).',
        const=True)

    parser.add_argument('alignment_file', metavar='ALG', type=str)
    parser.add_argument('reference_filename', metavar='REF', nargs='?', default=None, type=str)
    parser.add_argument('region', metavar='REG', type=str)


    args = parser.parse_args()

    if args.d != None:
        raise NotImplementedError('This argument was deprecated.')

    eqsts = EquispacedTilingStrategy(args.region, args.window_length[0], args.incr[0])

    build_windows(
        alignment_file = args.alignment_file,
        tiling_strategy = eqsts,
        win_min_ext = args.m[0],
        maximum_reads = args.x[0], # 1e4 / window_length, TODO why divide?
        minimum_reads = args.c[0],
        reference_filename = args.reference_filename
    )
