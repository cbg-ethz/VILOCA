import pysam
from typing import Optional
from shorah.tiling import TilingStrategy, EquispacedTilingStrategy
import numpy as np
import math
import logging
from multiprocessing import Process, cpu_count
import os
import hashlib

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
                    hashlib.sha1(pileupread.alignment.cigarstring.encode()).hexdigest(),
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

def _build_one_full_read(full_read: list[str], full_qualities: list[int]|list[str],
    read_query_name: str|None, full_read_cigar_hash: str|None, first_aligned_pos,
    last_aligned_pos, indel_map, max_ins_at_pos,
    extended_window_mode: bool, insertion_char: str) -> tuple[str, list[int]]:

    assert insertion_char in ["-", "X"], "Illegal char represention insertion"

    all_inserts = dict()
    own_inserts = set()

    change_in_reference_space_ins = 0

    for name, start, cigar_hash, ref_pos, indel_len, is_del in indel_map:

        if name == read_query_name and start == first_aligned_pos and cigar_hash == full_read_cigar_hash:
            if indel_len > 0 and is_del == 1:
                logging.debug(f"[b2w] Del and ins at same position in {read_query_name} @ {ref_pos}")

            if indel_len > 0 and not extended_window_mode:
                for _ in range(indel_len):
                    full_read.pop(ref_pos + 1 - first_aligned_pos)
                    if full_qualities is not None:
                        full_qualities.pop(ref_pos + 1 - first_aligned_pos)

            if is_del == 1: # if del
                full_read.insert(ref_pos - first_aligned_pos + change_in_reference_space_ins, "-")
                if full_qualities is not None:
                    full_qualities.insert(ref_pos - first_aligned_pos + change_in_reference_space_ins, "2")

            if indel_len > 0 and extended_window_mode:
                own_inserts.add((ref_pos, indel_len))
                change_in_reference_space_ins += indel_len
                all_inserts[ref_pos] = max_ins_at_pos[ref_pos]



        if (extended_window_mode and
            (name != read_query_name or start != first_aligned_pos or cigar_hash != full_read_cigar_hash) and
            first_aligned_pos <= ref_pos <= last_aligned_pos and indel_len > 0):

            all_inserts[ref_pos] = max_ins_at_pos[ref_pos]

    if extended_window_mode:
        change_in_reference_space = 0
        own_inserts_pos = []
        own_inserts_len = []
        if len(own_inserts) != 0:
            [own_inserts_pos, own_inserts_len] = [list(t) for t in zip(*own_inserts)]

        for pos in sorted(all_inserts):
            n = all_inserts[pos] # TODO does all_inserts lead to the same behavior max_ins_at_pos
            if (pos, n) in own_inserts:
                change_in_reference_space += n
                continue

            L = max_ins_at_pos[pos]
            in_idx = pos + 1 - first_aligned_pos + change_in_reference_space
            if pos in own_inserts_pos:
                k = own_inserts_len[own_inserts_pos.index(pos)]
                L -= k
                in_idx += k
            for _ in range(L):
                full_read.insert(in_idx, insertion_char)
                if full_qualities is not None:
                    full_qualities.insert(in_idx, "2")

            change_in_reference_space += max_ins_at_pos[pos]

    full_read = ("".join(full_read))

    return full_read, full_qualities # TODO return same data type twice


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

        full_read, full_qualities = _build_one_full_read(full_read, full_qualities,
            read.query_name, hashlib.sha1(read.cigarstring.encode()).hexdigest(), first_aligned_pos, last_aligned_pos,
            indel_map, max_ins_at_pos, extended_window_mode, "-")


        if (first_aligned_pos + minimum_overlap < window_start + 1 + window_length
                and last_aligned_pos >= window_start + minimum_overlap - 2 # TODO justify 2
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
                logging.warning("[b2w] No sequencing quality scores provided in alignment file. Run --sampler learn_error_params.")
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
    exclude_non_var_pos_threshold
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



def update_tiling(tiling, extended_window_mode, max_ins_at_pos):
    """
    input tiling:

    return: tiling = [
            (window_start, original_window_length, control_window_length, counter)
            for each window
            ]
    """
    update_tiling = []

    for idx, (window_start, window_length) in enumerate(tiling):
        original_window_length = window_length
        if extended_window_mode:
            for pos, val in max_ins_at_pos.items():
                if window_start <= pos < window_start + original_window_length:
                    window_length += val
            update_tiling.append((window_start,original_window_length, window_length))
        else:
            update_tiling.append((window_start,original_window_length, window_length))

    return update_tiling


def build_windows(alignment_file: str, tiling_strategy: TilingStrategy,
    win_min_ext: float, maximum_reads: int, minimum_reads: int,
    reference_filename: str,
    exact_conformance_fix_0_1_basing_in_reads: Optional[bool] = False,
    extended_window_mode: Optional[bool] = False,
    exclude_non_var_pos_threshold: Optional[float] = -1,
    maxthreads: Optional[int] = 1) -> None:
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
    #reffile = pysam.FastaFile(reference_filename) --> we need to read it in each child process
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

    # generate counter for each window
    # counter = window_start - 1 + control_window_length, # make 0 based
    counter_list = [0] + [window_start- 1 + control_window_length for (window_start, window_length, control_window_length) in tiling]

    all_processes = []
    for idx, (window_start, window_length, control_window_length) in enumerate(tiling):
        counter = counter_list[idx]
        p = Process(
            target=parallel_run_one_window,
            args=(
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
                )
            )
        all_processes.append(p)

    for p in all_processes:
      p.start()

    for p in all_processes:
      p.join()

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
