import pysam
import tempfile

def _annotate_alignment_file(
        alignment_filename, reference_filename, out_filename, sample_name):

    infile = pysam.AlignmentFile(
        alignment_filename,
        "r",
        reference_filename=reference_filename,
        threads=1
    )
    outfile = pysam.AlignmentFile(
        out_filename,
        "wb",
        reference_filename=reference_filename,
        template=infile,
        threads=1
    )
    for s in infile:
        s.query_name = f"__#{sample_name}#__{s.query_name}"
        outfile.write(s)

    outfile.close()
    infile.close()

def pre_process_pooled(alignment_filenames: list[str], reference_filename: str):
    arr = []
    for idx, i in enumerate(alignment_filenames):
        fp = tempfile.NamedTemporaryFile()
        sample_name = f"sample{idx}"
        _annotate_alignment_file(i, reference_filename, fp, sample_name)
        arr.append(fp)

    out_file = "merged.bam"

    pysam.merge("-f", "-o", out_file, *[i.name for i in arr])
    pysam.index(out_file)

    for i in arr:
        i.close()

    return out_file