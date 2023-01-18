import pytest
from shorah.shorah_snv import _compare_ref_to_read, SNP_id, SNV


@pytest.mark.parametrize("ref, seq, spec", [
    ("AACTTA", "AGA--A", {
        SNP_id(pos=2, var='G'): SNV(chrom='HXB2', haplotype_id='some-id', pos=2, ref='A', var='G', freq=1, support=0.1),
        SNP_id(pos=3, var='A'): SNV(chrom='HXB2', haplotype_id='some-id', pos=3, ref='C', var='A', freq=1, support=0.1),
        SNP_id(pos=4, var='--'): SNV(chrom='HXB2', haplotype_id='some-id', pos=3, ref='CTT', var='A', freq=1, support=0.1)
    }),
    ("AA--G", "A---T", {
        SNP_id(pos=2, var='---'): SNV(chrom='HXB2', haplotype_id='some-id', pos=1, ref='AA--', var='A', freq=1, support=0.1),
        SNP_id(pos=5, var='T'): SNV(chrom='HXB2', haplotype_id='some-id', pos=5,ref='G', var='T', freq=1, support=0.1),
    }),
    ("AA--X-GC", "AA--X-TC", {
        SNP_id(pos=6, var='T'): SNV(chrom='HXB2', haplotype_id='some-id', pos=6, ref='G', var='T', freq=1, support=0.1),
    }),
    ("AA-X-GCXGG", "AA-G-TCX-G", {
        SNP_id(pos=4, var='G'): SNV(chrom='HXB2', haplotype_id='some-id', pos=3, ref='-', var='-G', freq=1, support=0.1),
        SNP_id(pos=6, var='T'): SNV(chrom='HXB2', haplotype_id='some-id', pos=5, ref='G', var='T', freq=1, support=0.1),
        SNP_id(pos=7, var='-'): SNV(chrom='HXB2', haplotype_id='some-id', pos=6, ref='CG', var='C', freq=1, support=0.1) # FIXME pos=8 ?
    }),
    ("AXXG", "AGGG", {
        SNP_id(pos=2, var='GG'): SNV(chrom='HXB2', haplotype_id='some-id', pos=1, ref='A', var='AGG', freq=1, support=0.1)
    }),
    ("ACXXXXGT", "ACGGXCCT", { # Not a realistic case
        SNP_id(pos=3, var='GGC'): SNV(chrom='HXB2', haplotype_id='some-id', pos=2, ref='C', var='CGGC', freq=1, support=0.1),
        SNP_id(pos=6, var='C'): SNV(chrom='HXB2', haplotype_id='some-id', pos=3, ref='G', var='C', freq=1, support=0.1),
    })
])
def test_compare_ref_to_read(ref, seq, spec):
    snp = {}

    # WARNING: Pos is 1-based!
    tot_snv = _compare_ref_to_read(ref, seq, 1, snp, 1, 0.1, "HXB2", "some-id")

    assert snp == spec

    assert tot_snv == len(snp)