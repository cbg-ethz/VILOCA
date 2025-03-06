import unittest
from unittest.mock import patch
import hashlib

from viloca.b2w import _build_one_full_read


class TestBuildOneFullRead(unittest.TestCase):

    def setUp(self):
        # Initial setup for all test cases
        self.full_read = list("ACGTACGT")
        self.full_qualities = list("FFFFFFFF")
        self.read_query_name = "test_read"
        self.full_read_cigar_hash = hashlib.sha1("4M1I3M".encode()).hexdigest()
        self.first_aligned_pos = 100
        self.last_aligned_pos = 107
        self.extended_window_mode = False

    def test_no_indels(self):
        """Test when no indels are present."""
        # Empty indel map and max insertions
        self.indel_map = set()
        self.max_ins_at_pos = {}

        result = _build_one_full_read(
            self.full_read,
            self.full_qualities,
            self.read_query_name,
            self.full_read_cigar_hash,
            self.first_aligned_pos,
            self.last_aligned_pos,
            self.indel_map,
            self.max_ins_at_pos,
            False,  # extended_window_mode disabled
            "-"
        )

        # Expected result should be unchanged
        expected_result = ("ACGTACGT", list("FFFFFFFF"))
        self.assertEqual(result, expected_result)
    def test_insertion_in_another_read(self):
        """Test when insertion is present in another read."""
        # Setup indel map and max insertions
        self.indel_map = {
            (
                self.read_query_name,        # Read name
                self.first_aligned_pos,     # Reference start position
                self.full_read_cigar_hash,  # CIGAR hash
                103,                        # Position of insertion
                1,                          # Length of insertion
                False                       # Not a deletion (is_del=False)
            )
        }
        self.max_ins_at_pos = {103: 1}  # Maximum insertion length at position 103

        result = _build_one_full_read(
            self.full_read,
            self.full_qualities,
            None,  # No specific read query name for reference sequence processing
            None,  # No cigar hash for reference sequence processing
            self.first_aligned_pos,
            self.last_aligned_pos,
            self.indel_map,
            self.max_ins_at_pos,
            True,  # extended_window_mode enabled
            "-"     # Insertion from another read
        )

        # Expected result after inserting '-' at position 103
        expected_result = ("ACGT-ACGT", ['F', 'F', 'F', 'F', '2', 'F', 'F', 'F', 'F'])
        self.assertEqual(result, expected_result)

    def test_insertion_in_another_read_no_extended_window_mode(self):
        """Test when insertion is present in another read."""
        # Setup indel map and max insertions
        self.indel_map = {
            (
                self.read_query_name,        # Read name
                self.first_aligned_pos,     # Reference start position
                self.full_read_cigar_hash,  # CIGAR hash
                103,                        # Position of insertion
                1,                          # Length of insertion
                False                       # Not a deletion (is_del=False)
            )
        }
        self.max_ins_at_pos = {103: 1}  # Maximum insertion length at position 103

        result = _build_one_full_read(
            self.full_read,
            self.full_qualities,
            None,  # No specific read query name for reference sequence processing
            None,  # No cigar hash for reference sequence processing
            self.first_aligned_pos,
            self.last_aligned_pos,
            self.indel_map,
            self.max_ins_at_pos,
            False,  # extended_window_mode disabled
            "-"     # Insertion from another read
        )

        # Expected result after inserting '-' at position 103
        expected_result = ("ACGTACGT", ['F', 'F', 'F', 'F',  'F', 'F', 'F', 'F'])
        self.assertEqual(result, expected_result)

    def test_insertion_in_this_read(self):
        """Test when insertion is present in this read."""
        # Setup indel map and max insertions
        self.indel_map = {(self.read_query_name, 100, self.full_read_cigar_hash, 103, 1, False)}
        self.max_ins_at_pos = {103: 1}

        result = _build_one_full_read(
            self.full_read,
            self.full_qualities,
            None,  # No specific read query name for reference sequence processing
            None,  # No cigar hash for reference sequence processing
            self.first_aligned_pos,
            self.last_aligned_pos,
            self.indel_map,
            self.max_ins_at_pos,
            True,  # extended_window_mode disabled
            "X"  # Insertion from this read
        )

        # Expected result after inserting 'X' at position 103
        expected_result = ("ACGTXACGT", list("FFFF2FFFF"))
        self.assertEqual(result, expected_result)

    def test_insertion_in_this_read_no_extended_window_mode(self):
        """Test when insertion is present in this read."""
        # Setup indel map and max insertions
        self.indel_map = {(self.read_query_name, 100, self.full_read_cigar_hash, 103, 1, False)}
        self.max_ins_at_pos = {103: 1}

        result = _build_one_full_read(
            self.full_read,
            self.full_qualities,
            None,  # No specific read query name for reference sequence processing
            None,  # No cigar hash for reference sequence processing
            self.first_aligned_pos,
            self.last_aligned_pos,
            self.indel_map,
            self.max_ins_at_pos,
            False,  # extended_window_mode disabled
            "X"  # Insertion from this read
        )

        # Expected result after inserting 'X' at position 103
        expected_result = ("ACGTACGT", list("FFFFFFFF"))
        self.assertEqual(result, expected_result)

    def test_extended_window_mode(self):
        """Test with extended window mode enabled."""
        # Enable extended window mode and use '-' as insertion character
        insertion_char = "-"

        result = _build_one_full_read(
            list("ACGTACGT"),  # Reference sequence as input for extended mode
            list("FFFFFFFF"),
            None,  # No specific read query name for reference sequence processing
            None,  # No cigar hash for reference sequence processing
            99,    # Start of extended window (adjust as needed)
            108,   # End of extended window (adjust as needed)
            set(),  # Empty indel map for reference sequence processing
            {},     # No max insertions for reference sequence processing
            True,   # extended_window_mode enabled
            insertion_char
        )

        # Example expectation; adjust based on actual behavior of extended window mode
        expected_result = ("ACGTACGT",list("FFFFFFFF"))

        # Verify that the result matches expectations
        self.assertEqual(result, expected_result)






if __name__ == '__main__':
    unittest.main()
