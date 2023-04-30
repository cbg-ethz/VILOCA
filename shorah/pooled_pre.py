import pysam

def annotate_alignment_file(
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

    pysam.index(out_filename)