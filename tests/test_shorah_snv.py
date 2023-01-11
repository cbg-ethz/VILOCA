import pytest
from shorah.shorah_snv import _compare_ref_to_read, SNP_id, SNV



def test_compare_ref_to_read():
    snp = {}
    tot_snv = _compare_ref_to_read("AACTTA", "AGA--A", 1, snp, 1, 0.1, "HXB2", "some-id")

    assert snp == {
        SNP_id(pos=2, var='G'): SNV(chrom='HXB2', haplotype_id='some-id', pos=2, ref='A', var='G', freq=1, support=0.1),
        SNP_id(pos=3, var='A'): SNV(chrom='HXB2', haplotype_id='some-id', pos=3, ref='C', var='A', freq=1, support=0.1),
        SNP_id(pos=4, var='--'): SNV(chrom='HXB2', haplotype_id='some-id', pos=3, ref='CTT', var='C', freq=1, support=0.1)
    }

    assert tot_snv == len(snp)