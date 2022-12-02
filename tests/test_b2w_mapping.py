from array import array
import pytest
from cigar import Cigar
from shorah import b2w

class MockAlignedSegment:
    def __init__(self, query_name: str, reference_start: int, query_sequence: str, cigartuples: str):
        self.query_name = query_name
        self.reference_start = reference_start
        self.query_sequence = query_sequence # 0 based
        self.reference_end = reference_start + len(query_sequence)
        self.query_qualities = array("B", [38] * len(query_sequence))
        self.cigartuples = []
        # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
        pysam_cigar = {
            "M": 0,
            "I": 1,
            "D": 2,
            "N": 3,
            "S": 4
            # etc
        }
        for i in list(Cigar(cigartuples).items()):
            self.cigartuples.append((pysam_cigar[i[1]], i[0]))

    def add_indels(self, indels_map):
        cnt = self.reference_start
        for i in self.cigartuples:
            if i[0] == 1: # insert TODO Justify -1
                indels_map.append((self.query_name, self.reference_start, cnt-1, i[1], 0)) # cnt-1
            elif i[0] == 2: # del
                indels_map.append((self.query_name, self.reference_start, cnt, 0, 1))
                cnt += i[1]
            else:
                cnt += i[1]

# window_start is 0 based
@pytest.mark.parametrize("mArr,spec,window_length,window_start,extended_window_mode", [
    ([
        MockAlignedSegment(
            "89.6-2108",
            2291,
            "AAGTAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGACATTGAGTTGCCAGGGAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGAGCAGATAGACATAGAAATCTGTGGACATAAAGCTAAAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCA",
            "3M1D61M1I6M1D175M"
        )
    ],[
        "CAGATGATACAGTATTAGAAGAATTGAG-TTGCCAGGGAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGAGCAGATAGACATAGAAATCTGTGGACATAAAGCTAAAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTT"
    ], 201, 2334, False),
    ([
        MockAlignedSegment(
            "Q1",
            2291,
            "AAGTAGGGGGGCAACTAAAG",
            "20M"
        ),
        MockAlignedSegment(
            "Q2",
            2293,
            "AAGTAGGGGGGCAACTAAAG",
            "20M"
        ),
        MockAlignedSegment(
            "Q3",
            2281,
            "AAGTAGGGGGGCAACTAAAG",
            "20M"
        )
    ],[
        "AAGTAGGGGGGCAAC",
        "NNAAGTAGGGGGGCA",
        "GCAACTAAAGNNNNN"
    ], 15, 2291, False),
    ([
        MockAlignedSegment(
            "Q1",
            2291,
            "AAAAACCCGGAAAAAAAAAA",
            "5M3I12M"
        ),
        MockAlignedSegment(
            "Q2",
            2291,
            "GGGCCAAAAAAAAAAAAAAAAA",
            "3M2I17M"
        )
    ],[
        "AAA--AACCCGGAAAAAAAA",
        "GGGCCAA---AAAAAAAAAA"
    ], 15, 2291, True),
    ([
        MockAlignedSegment(
            "Q1",
            2291,
            "AAAAACCCGGAAAAAAAAAA",
            "5M3I12M"
        ),
        MockAlignedSegment(
            "Q2",
            2291,
            "GGGCCAAAAAAAAAAAAAAAAA",
            "4M2I16M"
        )
    ],[
        "AAAA--ACCCGGAAAAAAA",
        "GGGCCAA---AAAAAAAAA"
    ], 14, 2291, True),
    ([
        MockAlignedSegment(
            "Q1",
            2291,
            "AAAAACCCGGAAAAAAAAAACCCCCCC",
            "5M3I3M3I13M"
        ),
        MockAlignedSegment(
            "Q2",
            2291,
            "GGGCCAAACCGTTAAAAAAAAAAAGGC",
            "8M5I14M"
        )
    ],[
        "AAAAACCCGGAAAA--AAAAAACCCCCCCNNNN",
        "GGGCC---AAACCGTTAAAAAAAAAAAGGCNNN"
    ], 25, 2291, True)
])
def test_some_func(mArr, spec, window_length, window_start, extended_window_mode, mocker):

    indels_map = []
    for m in mArr:
        m.add_indels(indels_map)
    indels_map.sort(key=lambda tup: tup[2])
    print(indels_map)

    max_indel_at_pos = [0] * (window_start + window_length)
    for el in indels_map:
        if el[4] == 0 and el[3] > 0:
            if max_indel_at_pos[el[2]] < el[3]:
                max_indel_at_pos[el[2]] = el[3]


    mock_samfile = mocker.MagicMock()
    mock_samfile.configure_mock(
        **{
            "fetch.return_value": mArr
        }
    )

    mock_dict = mocker.MagicMock()
    mock_dict.__getitem__.return_value = 42

    arr, arr_read_qualities_summary, arr_read_summary, counter = b2w._run_one_window(
        mock_samfile,
        window_start,
        "HXB2-does-not-matter",
        window_length,
        0,
        mock_dict,
        0,
        True,
        indels_map,
        max_indel_at_pos,
        extended_window_mode
    )
    print(arr)

    for idx, el in enumerate(arr):
        assert el.split("\n")[1] == spec[idx]