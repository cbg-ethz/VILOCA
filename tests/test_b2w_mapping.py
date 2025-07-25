from array import array
import pytest
from cigar import Cigar
from viloca import b2w
import hashlib

class MockAlignedSegment:
    def __init__(self, query_name: str, reference_start: int, query_sequence: str, cigarstring: str):
        self.query_name = query_name
        self.reference_start = reference_start
        self.query_sequence = query_sequence # 0 based
        self.query_qualities = array("B", [38] * len(query_sequence))
        self.cigarstring = cigarstring
        self.cigartuples = []
        # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
        num_ins = 0
        num_dels = 0
        pysam_cigar = {
            "M": 0,
            "I": 1,
            "D": 2,
            "N": 3,
            "S": 4
            # etc
        }
        for i in list(Cigar(cigarstring).items()):
            self.cigartuples.append((pysam_cigar[i[1]], i[0]))
            if i[1] == "I":
                num_ins += i[0]
            if i[1] == "D":
                num_dels += i[0]

        self.reference_end = reference_start + len(query_sequence) - num_ins + num_dels

    def add_indels(self, indels_map):
        cnt = self.reference_start
        for i in self.cigartuples:
            if i[0] == 1: # insert TODO Justify -1
                indels_map.append((self.query_name, self.reference_start, hashlib.md5(self.cigarstring.encode()).hexdigest(), cnt-1, i[1], 0)) # cnt-1
            elif i[0] == 2: # del
                for k in range(i[1]):
                    indels_map.append((self.query_name, self.reference_start, hashlib.md5(self.cigarstring.encode()).hexdigest(), cnt+k, 0, 1))
                cnt += i[1]
            else:
                cnt += i[1]

# window_start is 0 based
@pytest.mark.parametrize("mArr,spec,window_length,window_start,extended_window_mode,exclude_non_var_pos_threshold", [
    ([
        MockAlignedSegment(
            "89.6-2108",
            2291,
            "AAGTAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGACATTGAGTTGCCAGGGAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGAGCAGATAGACATAGAAATCTGTGGACATAAAGCTAAAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCA",
            "3M1D61M1I6M1D175M"
        )
    ],[
        "CAGATGATACAGTATTAGAAGAATTGAG-TTGCCAGGGAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGAGCAGATAGACATAGAAATCTGTGGACATAAAGCTAAAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTT"
    ], 201, 2334, False, -1),
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
    ], 15, 2291, False, -1),
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
    ], 15, 2291, True, -1),
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
    ], 14, 2291, True, -1),
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
    ], 25, 2291, True, -1),
    ([
        MockAlignedSegment(
            "Q1",
            2291,
            "ABCDEFGH",
            "8M"
        ),
        MockAlignedSegment(
            "Q2",
            2291,
            "ABGCDEFGH",
            "2M1I2D6M"
        )
    ],[
        "AB-CDEFGH",
        "ABG--CDEF"
    ], 8, 2291, True, -1),
    ([
        MockAlignedSegment(
            "Ref",
            2291,
            "AGGTACGTAC",
            "10M"
        ),
        MockAlignedSegment(
            "R1",
            2291,
            "ACAGGAATAC",
            "2M2D3M2I3M"
        ),
        MockAlignedSegment(
            "R2",
            2291,
            "ACAAGTAGGATAC",
            "2M2I5M1I3M"
        ),
        MockAlignedSegment(
            "R3",
            2291,
            "ACAATTGAC",
            "2M2I1D3M2D2M"
        )
    ],[
        "AG--GTACG--TAC",
        "AC----AGGAATAC",
        "ACAAGTAGGA-TAC",
        "ACAA-TTG----AC"
    ], 10, 2291, True, -1),
    ([
        MockAlignedSegment(
            "Ref",
            2291,
            "AGGTACGTAC",
            "10M"
        ),
        MockAlignedSegment(
            "R1",
            2291,
            "ACAAGTGAAATAC",
            "2M2I2M2D1M3I3M"
        ),
        MockAlignedSegment(
            "R2",
            2291,
            "AAGTAAACGTAC",
            "1M1D1I3M2I5M"
        )
    ],[
        "AG--GTA--CG---TAC",
        "ACAAGT----GAAATAC",
        "A-A-GTAAACG---TAC"
    ], 10, 2291, True, -1),
    ([
        MockAlignedSegment(
            "Ref",
            2291,
            "AGGTACGTAC",
            "10M"
        ),
        MockAlignedSegment(
            "R1",
            2291,
            "ACAAGTGAAATAC",
            "2M2I2M2D1M3I3M"
        ),
        MockAlignedSegment(
            "R2",
            2291,
            "AAGTAAACGTAC",
            "1M1D1I3M2I5M"
        )
    ],[
        "AG--GTA--CG---TACNN",
        "ACAAGT----GAAATACNN",
        "A-A-GTAAACG---TACNN"
    ], 12, 2291, True, -1),
    ([
        MockAlignedSegment(
            "Ref",
            91,
            "AGGTACGTAC",
            "10M"
        ),
        MockAlignedSegment(
            "R1",
            94,
            "TAACGTAACAGG",
            "2M1I4M1I4M"
        ),
        MockAlignedSegment(
            "R2",
            88,
            "TTAAATGGTACGT",
            "4M2I7M"
        )
    ],[
        "A--GGTA-CGTA-C",
        "NNNNNTAACGTAAC",
        "AATGGTA-CGTNNN"
    ], 10, 91, True, -1),
    ([
        MockAlignedSegment(
            "Q1",
            2291,
            "AAGAAGGGGGGCAACTAAAG",
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
        "AAGAGGGGGCAAC",
        "NNAGTAGGGGGCA",
        "GCACTAAANNNNN"
    ], 15, 2291, False, 0.1),
    ([
        MockAlignedSegment(
            "Ref",
            2291,
            "AGGTACGTAC",
            "10M"
        ),
        MockAlignedSegment(
            "R1",
            2291,
            "ACAAGTGAAATAC",
            "2M2I2M2D1M3I3M"
        ),
        MockAlignedSegment(
            "R2",
            2291,
            "AAGTAAACGTAC",
            "1M1D1I3M2I5M"
        )
    ],[
        "G--A--C---",
        "CAA----AAA",
        "-A-AAAC---"
    ], 12, 2291, True, 0.1),
])
def test_run_one_window(mArr, spec, window_length, window_start, extended_window_mode, exclude_non_var_pos_threshold, mocker):

    indels_map = []
    for m in mArr:
        m.add_indels(indels_map)
    indels_map.sort(key=lambda tup: tup[3])
    print(indels_map)

    max_indel_at_pos = dict()
    for el in indels_map:
        try:
            if el[5] == 0 and el[4] > 0:
                if max_indel_at_pos[el[3]] < el[4]:
                    max_indel_at_pos[el[3]] = el[4]
        except KeyError:
            max_indel_at_pos[el[3]] = el[4]

    mock_samfile = mocker.MagicMock()
    mock_samfile.configure_mock(
        **{
            "fetch.return_value": mArr
        }
    )

    mock_dict = mocker.MagicMock()
    mock_dict.__getitem__.return_value = 42

    # added by Lara
    original_window_length = window_length
    control_window_length = window_length

    if extended_window_mode:
        for pos, val in max_indel_at_pos.items():
            if window_start <= pos < window_start + original_window_length:
                control_window_length += val


    arr, _, _, _, = b2w._run_one_window(
        mock_samfile,
        window_start, # 0 based
        "HXB2-does-not-matter",
        window_length,
        control_window_length,
        0,
        mock_dict,
        0,
        True,
        indels_map,
        max_indel_at_pos,
        extended_window_mode,
        exclude_non_var_pos_threshold
    )
    print(arr)

    for idx, el in enumerate(arr):
        assert el.split("\n")[1] == spec[idx]
