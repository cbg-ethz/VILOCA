from unittest.mock import patch, mock_open
from shorah import envp_post

DEFAULT_MOCK_DATA = "default mock data"

def test__post_process_for_envp_write_rec():

    data_dict = {"file1": ">hap_0|posterior=1 ave_reads=18\nTGGCAACGACTCGTCACAATAAAAGTAGGGGGGCAACGGAAGC\n>hap_1|posterior=1 ave_reads=20.9658\nTGCAGCGACTCGTCACAATAAAGATAGGGGGGCAATTGGAAGC",
                "file2": ">ENVP-template\nT=GCAACGAC===TCGTCACAATAAAAGTAGGGGGGCAACT===GGAAGC",
                "file3": ">ENVP-template\nT*GCAACGAC***TCGTCACAATAAAAGTAGGGGGGCAACT***GGAAGC"}

    def open_side_effect(name):
        return mock_open(read_data=data_dict.get(name, DEFAULT_MOCK_DATA))()

    with patch(f"{__name__}.open", side_effect=open_side_effect):
        with open("file1") as sup, open("file2") as ref, open("file3") as full_ref:
            out = envp_post._post_process_for_envp_write_rec(full_ref, ref, sup)

    assert str(out[0].seq) == "T*GGCAACGA***CTCGTCACAATAAAAGTAGGGGGGCAAC***GGAAGC"
    assert str(out[1].seq) == "T*GCAGCGAC***TCGTCACAATAAAGATAGGGGGGCAATT***GGAAGC"
