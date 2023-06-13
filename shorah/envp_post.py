import numpy as np
from Bio import SeqIO, Seq

def _post_process_for_envp_write_rec(
        full_ref_handle, ref_handle, support_handle):
    with full_ref_handle as template:
        full_ref = np.array(list(template.readlines()[1]))
        assert len(np.where(full_ref == "=")[0]) == 0
    with ref_handle as template:
        ref = np.array(list(template.readlines()[1]))
        excluded_pos = np.where(ref == "=")[0]

    records = []
    for s in SeqIO.parse(support_handle, "fasta"):
        seq = list(str(s.seq))
        for i in excluded_pos:
            seq.insert(i, full_ref[i])

        s.seq = Seq.Seq("".join(seq))
        records.append(s)

    return records

def post_process_for_envp(
        full_ref_handle, ref_handle, support_handle, new_support_handle):
    SeqIO.write(
        _post_process_for_envp_write_rec(
            full_ref_handle, ref_handle, support_handle),
        new_support_handle,
        "fasta"
    )
