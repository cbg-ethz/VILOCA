import pysam
import unittest
import tempfile
import hashlib

from viloca.b2w import _calc_via_pileup


class TestCalcViaPileup(unittest.TestCase):

    def setUp(self):
        # Create a temporary directory for the BAM file and index
        self.temp_dir = tempfile.TemporaryDirectory()
        self.bam_file_path = os.path.join(self.temp_dir.name, "test.bam")

        # Write a minimal BAM file with test data using pysam
        header = {
            'HD': {'VN': '1.0'},
            'SQ': [{'LN': 1000, 'SN': 'chr1'}]
        }
        with pysam.AlignmentFile(self.bam_file_path, "wb", header=header) as bam:
            # Add a single read with an adjusted query sequence and CIGAR string
            read = pysam.AlignedSegment()
            read.query_name = "read1"
            read.query_sequence = "ACGTACGTAACGT"  # Adjusted to match CIGAR string length (Option 1)
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 100
            read.mapping_quality = 60
            read.cigarstring = "4M1I4M1D4M"       # Matches adjusted query sequence length
            read.next_reference_id = -1
            read.next_reference_start = -1
            read.template_length = 0
            read.query_qualities = pysam.qualitystring_to_array("FFFFFFFFFFFFF")  # Adjusted length
            bam.write(read)

        # Create an index for the BAM file
        pysam.index(self.bam_file_path)

    def tearDown(self):
        # Clean up temporary directory and files
        self.temp_dir.cleanup()

    def test_calc_via_pileup(self):
        # Open the BAM file with pysam for testing
        with pysam.AlignmentFile(self.bam_file_path, "rb") as samfile:
            reference_name = "chr1"
            maximum_reads = 10

            # Call the function to be tested
            budget, indel_map, max_ins_at_pos = _calc_via_pileup(samfile, reference_name, maximum_reads)

            # Assertions to verify the output

            # Budget should have one entry for the position of the read's reference start
            self.assertIn(100, budget)
            self.assertEqual(budget[100], min(1, maximum_reads - 1))  # Only one read at position 100

            # Indel map should contain the insertion and deletion details
            self.assertEqual(len(indel_map), 2)  # One insertion and one deletion

            # Check insertion at position 103 (length 1)
            cigar_hash = hashlib.sha1("4M1I4M1D4M".encode()).hexdigest()
            self.assertIn(('read1', 100, cigar_hash, 103, 1, False), indel_map)

            # Check deletion at position 108 (length 0)
            self.assertIn(('read1', 100, cigar_hash, 108, 0, True), indel_map)

            # max_ins_at_pos should record the maximum insertion length at position 103
            self.assertIn(103, max_ins_at_pos)
            self.assertEqual(max_ins_at_pos[103], 1)  # Insertion of length 1

if __name__ == "__main__":
    unittest.main()
