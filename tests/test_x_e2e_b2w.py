import pytest
import filecmp
import os
import glob
from viloca import b2w, tiling
import math
import sys
import shutil

# Add necessary path for helper module
sys.path.append('/Users/lfuhrmann/Documents/Projects/ShoRAH2.0/Code_of_VILOCA/VILOCA/tests')
import helper_b2w_e2e

p = os.path.dirname(__file__)

def _collect_files(base_path):
    """Collect relevant output files for comparison."""
    spec_files = []
    for x in glob.glob(os.path.join(base_path, '*.reads.fas')):
        if os.path.getsize(x) > 0:  # Skip empty files
            spec_files.append(os.path.basename(x))
    spec_files.extend(['coverage.txt', 'reads.fas'])
    return spec_files

# Parameterize tests with both original cases and both threshold values
@pytest.mark.parametrize(
    "spec_dir,alignment_file,reference_file,region,window_length,overlap_factor,"
    "win_min_ext,maximum_reads,minimum_reads,extended_window_mode,exclude_non_var_pos_threshold",
    [
        # Original test cases (threshold=-1)
        ("data_1", "test_aln.cram", "test_ref.fasta", "HXB2:2469-2669",
         201, 3, 0.85, 497, 0, False, -1),  # Added False for extended_window_mode
        ("data_7", "reads.shotgun.bam", "reference.fasta", "MasterSequence:2000-2900",
         900, 3, 0.85, 20, 0, False, -1),

        # envp
        ("data_1", "test_aln.cram", "test_ref.fasta", "HXB2:2469-2669",
         201, 3, 0.85, 497, 0, False, 0.1),
        ("data_7", "reads.shotgun.bam", "reference.fasta", "MasterSequence:2000-2900",
         900, 3, 0.85, 20, 0, False, 0.1),  

         # Test cases ewm
         ("data_1", "test_aln.cram", "test_ref.fasta", "HXB2:2469-2669",
         201, 3, 0.85, 497, 0, True, -1),
        ("data_7", "reads.shotgun.bam", "reference.fasta", "MasterSequence:2000-2900",
         900, 3, 0.85, 20, 0, True, -1),

         # Test case ewm+envp
         ("data_1", "test_aln.cram", "test_ref.fasta", "HXB2:2469-2669",
         201, 3, 0.85, 497, 0, True, 0.1),
        ("data_7", "reads.shotgun.bam", "reference.fasta", "MasterSequence:2000-2900",
         900, 3, 0.85, 20, 0, True, 0.1)
    ],
    indirect=["spec_dir"]
)
def test_b2w(
    spec_dir, alignment_file, reference_file, region, window_length,
    overlap_factor, win_min_ext, maximum_reads, minimum_reads,
    extended_window_mode, exclude_non_var_pos_threshold
):
    # Validate window parameters
    assert window_length > 0 and window_length % overlap_factor == 0
    incr = window_length // overlap_factor
    minimum_overlap = math.floor(win_min_ext * window_length)
    extended_window_mode = False

    # Setup output directories
    helper_dir = os.path.join(p, spec_dir, "output_helper")
    b2w_dir = os.path.join(p, spec_dir, "output_b2w")
    os.makedirs(helper_dir, exist_ok=True)
    os.makedirs(b2w_dir, exist_ok=True)

    strategy = tiling.EquispacedTilingStrategy(region, window_length, incr, True)

    # Run helper implementation
    os.chdir(helper_dir)
    returncode = helper_b2w_e2e.build_windows(
        alignment_file=os.path.join(p, spec_dir, alignment_file),
        tiling_strategy=strategy,
        win_min_ext=win_min_ext,
        maximum_reads=maximum_reads,
        minimum_reads=minimum_reads,
        reference_filename=os.path.join(p, spec_dir, reference_file),
        exact_conformance_fix_0_1_basing_in_reads=False,
        extended_window_mode=extended_window_mode,
        exclude_non_var_pos_threshold=exclude_non_var_pos_threshold,
    )
    assert returncode is None

    # Run b2w implementation
    os.chdir(b2w_dir)
    b2w.build_windows(
        alignment_file=os.path.join(p, spec_dir, alignment_file),
        tiling_strategy=strategy,
        win_min_ext=win_min_ext,
        maximum_reads=maximum_reads,
        minimum_reads=minimum_reads,
        reference_filename=os.path.join(p, spec_dir, reference_file),
        exact_conformance_fix_0_1_basing_in_reads=False,
        extended_window_mode=extended_window_mode,
        exclude_non_var_pos_threshold=exclude_non_var_pos_threshold,
    )

    # Compare output files
    spec_files = _collect_files(helper_dir)
    match, mismatch, errors = filecmp.cmpfiles(
        helper_dir, b2w_dir, spec_files, shallow=False
    )
    assert len(mismatch) == 0 and len(errors) == 0

@pytest.fixture
def spec_dir(request):
    yield request.param

    # Cleanup after test
    helper_dir = os.path.join(p, request.param, "output_helper")
    b2w_dir = os.path.join(p, request.param, "output_b2w")
    if os.path.exists(helper_dir):
        shutil.rmtree(helper_dir)
    if os.path.exists(b2w_dir):
        shutil.rmtree(b2w_dir)
