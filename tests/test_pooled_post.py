from unittest.mock import patch, mock_open
from shorah import pooled_post
import numpy as np

DEFAULT_MOCK_DATA = "default mock data"

def test__ingest_sampler_output_to_calc_mean_cluster():

    data_dict = {"file1": ">hap_0|posterior=1 ave_reads=18\nTGGCAACGACCCCTCGTCACAATAAAAGTAGGGGGGCAACTAAAGGAAGC\n>hap_1|posterior=1 ave_reads=20.9658\nTGGCAGCGACCCCTCGTCACAATAAAGATAGGGGGGCAATTAAAGGAAGC",
                "file2": ">89.6-5766 2263|posterior=1\nTGGCAACGACCCCTCGTCACAATAAAAGTAGGGGGGCAACTAAAGGAAGC\n>89.6-8512 2268|posterior=1\nTGGCAACGACCCCTCGTCACAATAAAAGTAGGGGGGCAACTAAAGGAAGC\n>NL43-1404 2268|posterior=1\nTGGCAGCGACCCCTCGTCACAATAAAGATAGGGGGGCAATTAAAGGAAGC"}

    def open_side_effect(name):
        return mock_open(read_data=data_dict.get(name, DEFAULT_MOCK_DATA))()

    with patch(f"{__name__}.open", side_effect=open_side_effect):
        with open("file1") as sup, open("file2") as cor:
            out = pooled_post._ingest_sampler_output_to_calc_mean_cluster(pooled_post._create_unique_haplo_var(sup), cor, 3)


    np.testing.assert_array_equal(out, np.array([[1, 0],
                                                 [1, 0],
                                                 [0, 1]]))


def test__ingest_sampler_results_gamma_theta():
    with patch("builtins.open", mock_open(read_data="bla\n#gamma = 0.967662\n\nlolz\n#theta = 0.88\neof")):
        out = pooled_post._ingest_sampler_results_gamma_theta(open("debug/dbg.dbg"), "shorah")

        np.testing.assert_almost_equal(out[0][0], -0.032872426234558716)
        np.testing.assert_almost_equal(out[0][1], -3.431512269664515)

        np.testing.assert_almost_equal(out[1][0], -0.12783337150988489)
        np.testing.assert_almost_equal(out[1][1], -2.120263536200091)

def test_write_support_file_per_sample():
    data_dict = {"file1": ">hap_0|posterior=1 ave_reads=18\nTGGCAACGACCCCTCGTCACAATAAAAGTAGGGGGGCAACTAAAGGAAGC\n>hap_1|posterior=1 ave_reads=20.9658\nTGGCAGCGACCCCTCGTCACAATAAAGATAGGGGGGCAATTAAAGGAAGC",
                "file2": ""}

    def open_side_effect(name):
        return mock_open(read_data=data_dict.get(name, DEFAULT_MOCK_DATA))()

    with patch(f"{__name__}.open", side_effect=open_side_effect):
        with open("file1") as sup, open("file2") as new_sup:
            pooled_post.write_support_file_per_sample(sup, new_sup, [0.6, 0.8], [66, 88])
            # TODO

# pooled_post.post_process_pooled_samples_mode("raw_reads/w-HXB2-2938-3138.ref.fas", "raw_reads/w-HXB2-2938-3138.reads.fas",
#                                      open("debug/w-HXB2-2938-3138.dbg"),
#                                      open("support/w-HXB2-2938-3138.reads-support.fas"),
#                                      open("corrected/w-HXB2-2938-3138.reads-cor.fas"),
#                                      "shorah") # TODO