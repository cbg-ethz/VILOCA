#!/usr/bin/env python3

# Copyright 2007-2018
# Niko Beerenwinkel,
# Nicholas Eriksson,
# Moritz Gerstung,
# Lukas Geyrhofer,
# Osvaldo Zagordi,
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

"""dec.py draws windows on the bam file, cuts them into separate fasta files,
    calls diri_sampler to correct them, and merges the correction into a
    file of corrected reads
"""
import os
import pipes
import sys
import logging
import re
import shutil
import numpy as np
from multiprocessing import Process, Pool, cpu_count
import glob
import time
import pysam
import random
import subprocess
from Bio import SeqIO
import gzip
from pathlib import Path
import tempfile

import libshorah

from . import shorah_snv
from . import b2w
from . import tiling
from . import pooled_pre
from . import pooled_post
from . import envp_post

# import local haplotype inference methods
from .local_haplotype_inference.use_quality_scores import run_dpm_mfa as use_quality_scores
from .local_haplotype_inference.learn_error_params import run_dpm_mfa as learn_error_params


#################################################
# a common user should not edit above this line #
#################################################
# parameters not controlled by command line options
fasta_length = 80   # controls line length in fasta files
#win_min_ext = 0.85  # if read covers at least win_min_ext fraction of
# the window, fill it with Ns
hist_fraction = 0.20  # fraction that goes into the history
min_quality = 0.9  # quality under which discard the correction
min_x_thresh = 10  # threshold of X in the correction
init_K = 20  # initial number of clusters in diri_sampler

#################################################
# a common user should not edit below this line #
#################################################

# Dictionary storing by read ID (key) all subsequences which are part of
# different windows
correction = {}
# Dictionary storing by read ID (key) the posterior for each of window
quality = {}

count = {
    'A': 0,
    'C': 0,
    'G': 0,
    'T': 0,
    'X': 0,
    '-': 0
}
clusters = [[]]
untouched = [[]]


def parse_aligned_reads(reads_file):
    """Parse reads from a file with aligned reads
    """
    out_reads = {}

    if not os.path.isfile(reads_file):
        logging.error('There should be a file here: %s', reads_file)
        sys.exit('There should be a file here: ' + reads_file)
    else:
        logging.info('Using file %s of aligned reads', reads_file)

    handle = open(reads_file)
    logging.debug('Parsing aligned reads')

    for h in handle:
        name, start, stop, mstart, mstop, this_m = h.rstrip().split('\t')
        if this_m == '':
            logging.warning('parsing empty read: %s', h)
        out_reads[name] = [None, None, None, None, []]
        # start of the region of interest (0-based indexing)
        out_reads[name][0] = int(start)
        # end of the region of interest (0-based inde
        out_reads[name][1] = int(stop)
        # start index of the aligned read w.r.t reference (1-based indexing)
        out_reads[name][2] = int(mstart)
        # end index of the aligned read w.r.t reference (1-based indexing)
        out_reads[name][3] = int(mstop)
        out_reads[name][4] = this_m

    return out_reads


def b2w_logging(run_settings):
    """run b2w to make windows from bam
    """
    bam, fasta, w, i, m, x, c, reg, ignore_indels = run_settings
    d = ' -d' if ignore_indels else ''
    my_arg = '-w %i -i %i -m %.4f -x %i -c %i%s %s %s %s' % \
        (w, i, m, x, c, d, bam, fasta, reg)
    logging.debug(f'To run standalone: python3 b2w.py {my_arg}')


def run_dpm(run_setting):
    """run the dirichlet process clustering
    """

    filein, j, a, seed, inference_type, n_max_haplotypes, n_mfa_starts, unique_modus, inference_convergence_threshold, record_history = run_setting

    ref_in = filein.split('.reads.')[0] + str('.ref.fas')
    fname_qualities = filein.split('.reads.')[0] + str('.qualities.npy')

    logging.debug('Running sampler')
    if inference_type == 'shorah': # run the original sampler of ShoRAH

        # dn = sys.path[0]
        #my_prog = shlex.quote(diri_exe)  # os.path.join(dn, 'diri_sampler')
        #my_arg = ' -i %s -j %i -t %i -a %f -K %d -R %d' % \
        #    (pipes.quote(filein), j, int(j * hist_fraction), a, init_K, seed)

        # TODO integration
        try:
            os.remove('./corrected.tmp')
            # os.remove('./assignment.tmp')
        except OSError:
            pass


        # runs the gibbs sampler for the dirichlet process mixture
        try:
            logging.debug(f"{filein} {j} {a} {int(j * hist_fraction)} {init_K} {seed}")
            retcode = libshorah.exec_dpm_sampler(
                pipes.quote(filein),
                j,
                a,
                int(j * hist_fraction),
                K_cluster_start=init_K,
                R_seed=seed
            )
            if retcode == 0:
                logging.debug(f'{filein} - Run finished successfully.')
            else:
                logging.error(f'{filein} - Run failed with return code %i.', retcode)
        except Exception as e:
            logging.error(f'{filein} - Run failed: {e}')

    elif inference_type == 'use_quality_scores':
        use_quality_scores.main(freads_in=filein,
                     fref_in=ref_in,
                     fname_qualities= fname_qualities,
                     output_dir='./',
                     n_starts=int(n_mfa_starts),
                     K=int(n_max_haplotypes),
                     alpha0=float(a),
                     alphabet = 'ACGT-',
                     unique_modus = unique_modus,
                     convergence_threshold = inference_convergence_threshold,
                     record_history = record_history,
                     seed=seed,
                     )

    elif inference_type == 'learn_error_params':
        learn_error_params.main(freads_in=filein,
                     fref_in=ref_in,
                     output_dir='./',
                     n_starts=int(n_mfa_starts),
                     K=int(n_max_haplotypes),
                     alpha0=float(a),
                     alphabet = 'ACGT-',
                     unique_modus = unique_modus,
                     #convergence_threshold = inference_convergence_threshold,
                     record_history = record_history,
                     seed=seed,
                     )
    logging.debug('Finished sampler')

    return


def correct_reads(chr_c, wstart, wend):
    """ Parses corrected reads (in fasta format) and correct the reads
    """
    # out_reads[match_rec.id][0] = qstart
    # out_reads[match_rec.id][1] = qstop
    # out_reads[match_rec.id][2] = mstart
    # out_reads[match_rec.id][3] = mstop
    # out_reads[match_rec.id][4] = Sequence...


    try:
        cor_file = 'w-%s-%s-%s.reads-cor.fas' % (chr_c, wstart, wend)
        if os.path.exists('corrected/' + cor_file + 'gz'):
            handle = gzip.open( # TODO to be removed
                cor_file + 'gz', 'rb' if sys.version_info < (3, 0) else 'rt')
        elif os.path.exists('corrected/' + cor_file):
            handle = open('corrected/' + cor_file, 'r')
        else:
            handle = open(cor_file, 'r')

        for seq_record in SeqIO.parse(handle, 'fasta'):
            read_id = seq_record.id
            try:
                correction[read_id][wstart] = list(str(seq_record.seq))
                quality[read_id][wstart] = \
                    float(seq_record.description.split('|')[1].split('=')[1])
                assert quality[read_id][wstart] <= 1.0, \
                    'try: quality must be < 1, %s' % cor_file
            except KeyError:
                correction[read_id] = {}
                quality[read_id] = {}
                correction[read_id][wstart] = list(str(seq_record.seq))
                quality[read_id][wstart] = \
                    float(seq_record.description.split('|')[1].split('=')[1])
                assert quality[read_id][wstart] <= 1.0, \
                    'except: quality must be < 1, %s' % cor_file
        handle.close()
        return
    except IOError:
        logging.warning('No reads in window %s?', wstart)
        return


def get_prop(filename):
    """fetch the number of proposed clusters from .dbg file
    """

    if os.path.exists(filename):
        h = open(filename)
    elif os.path.exists(filename + '.gz'):
        h = gzip.open(filename + '.gz', 'rb')
    elif os.path.exists('debug/' + filename):
        h = open('debug/' + filename)
    elif os.path.exists('debug/' + filename + '.gz'):
        h = gzip.open('debug/' + filename + '.gz', 'rb')
    else:
        return 'not found'

    prop = 'not found'
    for l in h:
        if l.startswith('#made'):
            prop = int(l.split()[1])
            break
    h.close()
    return prop


def base_break(baselist):
    """Break the tie if different corrections are found."""

    for c1 in count:
        count[c1] = 0
    for c in baselist:
        if c.upper() != 'N':
            count[c.upper()] += 1

    maxm = 0
    out = []
    for b in count:
        if count[b] >= maxm:
            maxm = count[b]
    for b in count:
        if count[b] == maxm:
            out.append(b)

    rc = random.choice(out)

    return rc


def win_to_run(alpha_w, seed, inference_type, n_max_haplotypes, n_mfa_starts, unique_modus, inference_convergence_threshold, record_history, reuse_files):
    """Return windows to run on diri_sampler."""

    rn_list = []
    try:
        file1 = open('coverage.txt')
    except IOError:
        sys.exit('Coverage file generated by b2w not found.')

    for f1 in file1:
        winFile, chr1, beg, end, cov = f1.rstrip().split('\t')
        j = min(300_000, int(cov) * 15)
        output_name = winFile.split(".fas")[0] + "-" + "cor.fas"
        if not (reuse_files and os.path.isfile(output_name)):
            rn_list.append((winFile, j, alpha_w, seed, inference_type, n_max_haplotypes, n_mfa_starts, unique_modus, inference_convergence_threshold, record_history))
        else:
            logging.info(f'[file already exits] Use {output_name} generated on {time.ctime(Path(output_name).stat().st_mtime)}')


    del end
    del(beg, chr1)
    return rn_list


def merge_corrected_reads(aligned_read):
    if aligned_read is None:
        print("empty window found", file=sys.stderr)
        return (None, [])

    ID = aligned_read[0]
    seq = aligned_read[1][4]
    corrected_read = correction.get(ID)
    merged_corrected_read = []

    if corrected_read is not None:
        # rlen: length of the orginal read
        rlen = len(seq)
        # rstart: start index of the original read with respect to reference
        rstart = aligned_read[1][2]
        # posterior: fraction of times a given read was assigned to current
        #            cluster among those iterations that were recorded
        posterior = quality.get(ID)

        # kcr: extract start index of the aligned reads in all windows
        # vcr: extract sequence of the aligned reads in all windows
        kcr = np.array(list(corrected_read.keys()), dtype=int)
        vcr = np.array([np.array(v) for v in corrected_read.values()], dtype=object) # FIXED dtype
        vcr_len = [v.size for v in vcr]


        for rpos in range(rlen):
            tp = rstart + rpos - kcr
            mask = np.logical_and(
                np.logical_and(np.greater_equal(tp, 0), np.less(tp, vcr_len)),
                np.greater(list(posterior.values()), min_quality)
            )
            if np.sum(mask > 0):
                # data structure needs this
                idx = np.argwhere(mask)
                this = [vcr[k][tp[k]]
                        for k in idx[0]]  # this is unlikely to be optimal
                if len(this) > 1:
                    corrected_base = base_break(this)
                else:
                    corrected_base = this[0]
            else:
                corrected_base = 'X'
            merged_corrected_read.append(corrected_base)

    return(ID, merged_corrected_read)

def move_files_into_dir(dir, files):
    for f in files:
        if os.stat(f).st_size > 0:
            try:
                os.remove(f"{dir}/{f}")
            except OSError:
                pass
            shutil.move(f, f"{dir}/")
        else:
            os.remove(f)



# def main(in_bam, in_fasta, win_length=201, win_shifts=3, region='',
#         max_coverage=10000, alpha=0.1, keep_files=True, seed=None):
def main(args):
    """
    Performs the error correction analysis, running diri_sampler
    and analyzing the result
    """

    in_bam = args.b
    in_fasta = args.f
    win_length = args.w # TODO remove this var
    win_shifts = args.win_shifts # TODO remove this var
    region = args.r
    max_coverage = args.max_coverage
    alpha = args.a
    cov_thrd = args.cov_thrd
    seed = args.seed
    ignore_indels = args.ignore_indels
    maxthreads = args.maxthreads
    path_insert_file = args.path_insert_file
    inference_type = args.inference_type
    n_max_haplotypes = args.n_max_haplotypes
    n_mfa_starts = args.n_mfa_starts
    unique_modus = args.unique_modus
    inference_convergence_threshold = args.conv_thres
    extended_window_mode = args.extended_window_mode
    exclude_non_var_pos_threshold = args.exclude_non_var_pos_threshold
    win_min_ext = args.win_min_ext
    reuse_files = args.reuse_files
    record_history = args.record_history

    logging.info(' '.join(sys.argv))

    if len(in_bam) == 1:
        in_bam = in_bam[0]
    else:
        in_bam = pooled_pre.pre_process_pooled(in_bam, in_fasta)

    # check options
    if win_length % win_shifts != 0:
        sys.exit('Window size must be divisible by win_shifts')
    if win_min_ext < 1 / win_shifts:
        logging.warning('Some bases might not be covered by any window')
    if max_coverage / win_length < 1:
        sys.exit('Please increase max_coverage')
    if not os.path.isfile(in_bam):
        sys.exit("File '%s' not found" % in_bam)
    if not os.path.isfile(in_fasta):
        sys.exit("File '%s' not found" % in_fasta)
    if seed is None:
        seed = np.random.randint(100, size=1)

    incr = win_length // win_shifts

    # run b2w

    logging.info('starting b2w')
    if not (reuse_files and os.path.isfile(f"coverage.txt")):
        try:
            if ignore_indels == True:
                raise NotImplementedError('This argument was deprecated.')
            b2w_logging((in_bam, in_fasta, win_length, incr, win_min_ext,
                max_coverage, cov_thrd, region, ignore_indels))

            if path_insert_file == None and region == "": # special case if no region defined
                samfile = pysam.AlignmentFile(
                    in_bam,
                    "r", # auto-detect bam/cram (rc)
                    reference_filename=in_fasta,
                    threads=1
                )
                if samfile.nreferences != 1:
                    raise NotImplementedError("There are multiple references in this alignment file.")
                strategy = tiling.EquispacedTilingStrategy(
                    f"{samfile.references[0]}:1-{samfile.lengths[0]}",
                    win_length,
                    incr,
                    False,
                    True
                )
            elif path_insert_file == None:
                strategy = tiling.EquispacedTilingStrategy(region, win_length, incr, True)
            else:
                strategy = tiling.PrimerTilingStrategy(path_insert_file)
                if region != "":
                    logging.warn(f"region is set to {region} but is not used with this tiling strategy")

            logging.info(f"Using tiling strategy: {type(strategy).__name__}")

            b2w.build_windows(
                in_bam,
                strategy,
                win_min_ext,
                max_coverage,
                cov_thrd,
                in_fasta,
                extended_window_mode=extended_window_mode,
                exclude_non_var_pos_threshold=exclude_non_var_pos_threshold,
                maxthreads=maxthreads,
                reuse_files=reuse_files
            )
            logging.info('finished b2w')

        except Exception as e:
            logging.debug(e)
            sys.exit('b2w run not successful')
    else:
        logging.info(f'[file already exits] Use coverage.txt generated on {time.ctime(Path("coverage.txt").stat().st_mtime)}')

    aligned_reads = parse_aligned_reads('reads.fas')
    if len(aligned_reads) == 0:
        msg = 'No reads found in the requested region %s' % region
        logging.debug(msg)
        print(msg, file=sys.stderr)
        logging.info('shotgun run ends with no processing')
        # technically it is a success: we did produce what we were asked for.
        # it just happens that we were asked to produce nothing and thus got nothing to output
        sys.exit(0)

    r = list(aligned_reads.keys())[0]
    gen_length = aligned_reads[r][1] - aligned_reads[r][0]

    if win_length > gen_length:
        sys.exit('The window size must be smaller than the genome region')

    logging.info('%s reads are being considered', len(aligned_reads))

    ############################################
    # Now the windows and the error correction #
    ############################################

    runlist = win_to_run(alpha, seed, inference_type, n_max_haplotypes, n_mfa_starts, unique_modus, inference_convergence_threshold, record_history, reuse_files)
    logging.info('will run on %d windows', len(runlist))
    # run diri_sampler on all available processors but one
    max_proc = max(cpu_count() - 1, 1)
    if maxthreads:
        max_proc = min(max_proc, maxthreads)
    logging.info('CPU(s) count %u, max thread limit %u, will run %u parallel dpm_sampler', cpu_count(), maxthreads, max_proc)

    all_processes = [Process(target=run_dpm, args=(run_set,)) for run_set in runlist]
    for p in all_processes:
        p.start()

    for p in all_processes:
        p.join()
        if p.exitcode != 0:
            logging.debug("[shotgun] A process was killed. Terminating program.")
            exit(1)

    logging.debug("[shotgun] All processes completed successfully.")

    # prepare directories
    for sd_name in ['debug', 'haplotypes', 'freq', 'sampling', 'work',
                    'corrected', 'raw_reads', 'inference']:
        try:
            os.mkdir(sd_name)
        except OSError:
            pass

    # parse corrected reads
    proposed = {}
    for i in runlist:
        winFile, j, a, s, inference_type, n_max_haplotypes, n_mfa_starts, unique_modus, inference_convergence_threshold, record_history = i
        del a  # in future alpha might be different on each window
        del s
        # greedy re match to handle situation where '.' or '-' appears in the
        # ID
        parts = re.match(
            r'^w-(?P<chrom>.*)-(?P<beg>\d+)-(?P<end>\d+).reads', winFile)
        chrom = parts.group('chrom')
        beg = int(parts.group('beg'))
        end = int(parts.group('end'))
        del parts
        logging.info('reading windows for start position %s', beg)
        # correct reads populates correction and quality, globally defined
        correct_reads(chrom, beg, end)
        stem = 'w-%s-%u-%u' % (chrom, beg, end)
        logging.info('this is window %s', stem)
        dbg_file = stem + '.dbg'
        # if os.path.exists(dbg_file):
        proposed[beg] = (get_prop(dbg_file), j)
        logging.info('there were %s proposed', str(proposed[beg][0]))

    move_files_into_dir("debug", glob.glob("./w*dbg"))
    move_files_into_dir("sampling", glob.glob("./w*smp"))
    move_files_into_dir("corrected", glob.glob("./w*reads-cor.fas"))
    move_files_into_dir("haplotypes", glob.glob("./w*reads-support.fas"))
    move_files_into_dir("freq", glob.glob("./w*reads-freq.csv"))
    move_files_into_dir("sampling", glob.glob("./w*smp"))
    raw_reads_files = glob.glob('./w*reads.fas') + glob.glob('./w*ref.fas') + glob.glob('./w*qualities.npy')
    move_files_into_dir("raw_reads", raw_reads_files)
    inference_files = glob.glob("./w*best_run.txt") + glob.glob("./w*history_run*.csv") + glob.glob("./w*results*.pkl")
    move_files_into_dir("inference", inference_files)

    ############################################
    ##      Print the corrected reads         ##
    ##
    ## correction[read_id][wstart] = sequence ##
    ## quality[read_id][wstart] = posterior   ##
    # ##########################################
    logging.info('Merging windows of corrected reads')
    # Multi-threaded version
    params = list(aligned_reads.items())
    pool = Pool(processes=max_proc)
    to_correct = pool.map(merge_corrected_reads, params)
    pool.close()
    pool.join()

    logging.info('All corrected reads have been merged')

    ccx = {}
    # handle case where bamfile has no dots in name
    cin_stem = re.sub(r'\.[^.]+$', r'', os.path.split(in_bam)[1])
    logging.debug('writing to file %s.cor.fas', cin_stem)
    with open('%s.cor.fas' % cin_stem, 'w') as fch:
        for ID, seq_list in to_correct:
            cor_read = ''.join(seq_list)
            init_x = len(cor_read.lstrip('-')) - len(cor_read.lstrip('-X'))
            fin_x = len(cor_read.rstrip('-')) - len(cor_read.rstrip('-X'))
            cx = seq_list.count('X') - init_x - fin_x
            ccx[cx] = ccx.get(cx, 0) + 1
            if cx <= min_x_thresh and cor_read.lstrip('-X') != '':
                fch.write('>%s %d\n' % (ID, aligned_reads[ID][2] + init_x
                                        - aligned_reads[ID][0]))
                cc = 0
                for c in cor_read.lstrip('-X'):
                    if c != 'X':
                        fch.write(str(c))
                        fch.flush()
                        cc = cc + 1
                        if cc % fasta_length == 0:
                            fch.write('\n')

                if cc % fasta_length != 0:
                    fch.write('\n')

    # write proposed_per_step to file
    ph = open('proposed.dat', 'w')
    ph.write('#base\tproposed_per_step\n')
    for kp in sorted(proposed):
        if proposed[kp] != 'not found' and 'not found' not in proposed[kp]:
            ph.write('%s\t%f\n' %
                     (kp, proposed[kp][0] / proposed[kp][1]))
    ph.close()

    logging.info('running snv.py')
    args.increment = win_length // win_shifts # TODO remove dependency on these vars

    # ENVP
    if exclude_non_var_pos_threshold > 0:
        with open("coverage.txt") as cov:
            for line in cov:
                window_file, _, _, _, _ = line.rstrip().split("\t")
                stem = window_file.split(".")[0]
                envp_post.post_process_for_envp(
                    open(f"raw_reads/{stem}.envp-full-ref.fas"),
                    open(f"raw_reads/{stem}.envp-ref.fas"),
                    f"haplotypes/{stem}.reads-support.fas",
                    f"haplotypes/{stem}.reads-support.fas" # overwrite
                )

    # Pooled
    b_list = args.b.copy()
    if len(b_list) > 1:
        for idx, i in enumerate(b_list):
            Path(f"sample{idx}/haplotypes").mkdir(parents=True, exist_ok=True)
            Path(f"sample{idx}/corrected").mkdir(parents=True, exist_ok=True)
            with open("coverage.txt") as cov:
                for line in cov:
                    window_file, _, _, _, _ = line.rstrip().split("\t")
                    stem = window_file.split(".")[0]
                    filtered_reads = tempfile.NamedTemporaryFile(mode="w", suffix=".fasta")
                    pooled_post.filter_fasta(filtered_reads, f"raw_reads/{stem}.reads.fas", f"sample{idx}")
                    filtered_reads.seek(0)
                    filtered_cor_reads_path = f"sample{idx}/corrected/{stem}.reads-cor.fas"
                    filtered_cor_reads = open(filtered_cor_reads_path, "w+")
                    pooled_post.filter_fasta(filtered_cor_reads, f"corrected/{stem}.reads-cor.fas", f"sample{idx}")
                    filtered_cor_reads.close()

                    posterior_and_avg = pooled_post.recalculate_posterior_and_ave_reads(
                        f"raw_reads/{stem}.ref.fas",
                        filtered_reads.name,
                        open(f"debug/{stem}.dbg") if inference_type == "shorah" else open(f"inference/{stem}.reads-all_results.pkl", "rb"),
                        open(f"haplotypes/{stem}.reads-support.fas"),
                        filtered_cor_reads_path,
                        inference_type,
                        None if inference_type != "use_quality_scores" else f"raw_reads/{stem}.qualities.npy" # TODO untested
                    )
                    filtered_reads.close()

                    pooled_post.write_support_file_per_sample(
                        open(f"haplotypes/{stem}.reads-support.fas"),
                        open(f"sample{idx}/haplotypes/{stem}.reads-support.fas", "w+"), # TODO
                        *posterior_and_avg
                    )

            args.b = i
            args.working_dir = f"sample{idx}"
            shorah_snv.main(args)
    else:
        args.b = b_list[0]
        args.working_dir = ""
        shorah_snv.main(args)

        # tidy snvs # TODO
        try:
            os.mkdir('snv')
        except OSError:
            os.rename('snv', 'snv_before_%d' % int(time.time()))
            os.mkdir('snv')

        for snv_file in glob.glob('./SNV*_final*')+ glob.glob('./cooccurring_mutations.csv'):
            shutil.copy(snv_file, 'snv/')
        for snv_file in glob.glob('./raw_snv*') + glob.glob('./SNV*.tsv'):
            shutil.move(snv_file, 'snv/')

    # now move all files that are not directly results into the debug directory
    shutil.move("inference", "work")
    shutil.move("raw_reads", "work")
    shutil.move("sampling", "work")
    shutil.move("debug", "work")
    shutil.move("freq", "work")
    shutil.move("corrected", "work")
    shutil.move("reads.fas", "work")
    shutil.move("proposed.dat", "work")
    shutil.move("snv", "work")
    shutil.move(glob.glob('*.cor.fas')[0], "work")

    logging.info('shotgun run ends')
    logging.info('VILOCA terminated')
    print("VILOCA terminated successfully.")
