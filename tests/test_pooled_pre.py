import pysam
import os
from viloca import pooled_pre

def test__annotate_alignment_file():
    out = "out.bam"
    pooled_pre._annotate_alignment_file(
        "data_1/test_aln.cram", "data_1/test_ref.fasta", out, "sample0")

    # TODO check closed status

    with pysam.AlignmentFile(
        out,
        "rb",
        threads=1
    ) as f:
        for i in f:
            assert i.query_name.startswith("__#sample0#__") == True

    # https://github.com/pysam-developers/pysam/issues/939
    # [E::idx_find_and_load] Could not retrieve index file for 'out.bam'
    # Warning may be ingnored

def test_pre_process_pooled():
    pooled_pre.pre_process_pooled([
        "data_1/test_aln.cram",
        "data_1/test_aln.cram"
    ], "data_1/test_ref.fasta")

    out = "merged.bam"

    a = [0,0,0]

    with pysam.AlignmentFile(
        out,
        "rb",
        threads=1
    ) as f:
        for i in f:
            if i.query_name.startswith("__#sample0#__"):
                a[0] += 1
            elif i.query_name.startswith("__#sample1#__"):
                a[1] += 1
            else:
                a[2] += 1

    os.remove(out)
    os.remove(out + ".bai")

    assert a[0] == a[1] != 0
    assert a[2] == 0
