import pysam
import os
from shorah import pooled_pre

def test_annotate_alignment_file():
    out = "out.bam"
    pooled_pre.annotate_alignment_file(
        "data_1/test_aln.cram", "data_1/test_ref.fasta", out, "sample0")

    # TODO check closed status

    with pysam.AlignmentFile(
        out,
        "rb",
        threads=1
    ) as f:
        for i in f:
            assert i.query_name.startswith("__#sample0#__") == True

    os.remove(out)
    os.remove(out + ".bai")