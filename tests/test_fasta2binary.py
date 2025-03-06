import pytest
import numpy as np
import sys
import os
from io import StringIO

from viloca.local_haplotype_inference.use_quality_scores.preparation  import load_and_process_reads
from viloca.local_haplotype_inference.use_quality_scores.preparation  import load_reference_seq, reference2binary


@pytest.fixture
def sample_fasta():
    return StringIO(">test_seq\nATCGNATCG\n")

def test_load_reference_seq(sample_fasta, tmp_path):
    # Write the sample fasta to a temporary file
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(sample_fasta.getvalue())

    alphabet = "ACGT"
    reference_binary, ref_id = load_reference_seq(str(fasta_file), alphabet)

    expected_binary = np.array([
        [1, 0, 0, 0],  # A
        [0, 0, 0, 1],  # T
        [0, 1, 0, 0],  # C
        [0, 0, 1, 0],  # G
        [0, 0, 0, 0],  # N
        [1, 0, 0, 0],  # A
        [0, 0, 0, 1],  # T
        [0, 1, 0, 0],  # C
        [0, 0, 1, 0],  # G
    ])

    np.testing.assert_array_equal(reference_binary, expected_binary)
    assert ref_id == "test_seq"
    print("Test Passed: Reference sequence loaded and converted to binary correctly.")

def test_reference2binary():
    alphabet = "ACGT"
    reference_seq = "ATCGNATCG"
    reference_binary = reference2binary(reference_seq, alphabet)

    expected_binary = np.array([
        [1, 0, 0, 0],  # A
        [0, 0, 0, 1],  # T
        [0, 1, 0, 0],  # C
        [0, 0, 1, 0],  # G
        [0, 0, 0, 0],  # N
        [1, 0, 0, 0],  # A
        [0, 0, 0, 1],  # T
        [0, 1, 0, 0],  # C
        [0, 0, 1, 0],  # G
    ])

    np.testing.assert_array_equal(reference_binary, expected_binary)
    print("Test Passed: Reference sequence converted to binary correctly.")


@pytest.fixture
def sample_fasta_file(tmp_path):
    fasta_content = (
        ">seq1\n"
        "ATCG\n"
        ">seq2\n"
        "GCTA\n"
        ">seq3\n"
        "GCTA\n"
    )
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)
    return str(fasta_file)

@pytest.fixture
def sample_qualities_file(tmp_path):
    qualities = np.array([[40, 40, 40, 40], [30, 30, 30, 30], [20, 20, 20, 20]])
    qualities_file = tmp_path / "test_qualities.npy"
    np.save(qualities_file, qualities)
    return str(qualities_file)

def test_fasta_to_reads_seq_binary(sample_fasta_file, sample_qualities_file):
    alphabet = "ACGT"

    print("unique_modus=True")
    # Process with old function
    unique_modus=True

    # Process with new function
    reads_seq_binary_new, reads_weights_new, qualities_new, _ = load_and_process_reads(sample_fasta_file, sample_qualities_file, alphabet, unique_modus)

    # Additional checks
    try:
        assert reads_seq_binary_new.shape == (2, 4, 4)  # 3 sequences, 4 nucleotides, 4 alphabet letters
        print("Test Passed: Binary sequence shape is correct.")
    except AssertionError:
        print("Test Failed: Binary sequence shape is incorrect.")

    try:
        assert reads_weights_new.shape == (2,)
        print("Test Passed: Read weights shape is correct.")
    except AssertionError:
        print("Test Failed: Read weights shape is incorrect.")

    try:
        assert qualities_new.shape == (2, 4)
        print("Test Passed: Qualities shape is correct.")
    except AssertionError:
        print("Test Failed: Qualities shape is incorrect.")

    # Check if binary representation is correct
    expected_binary = np.array([
        [[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]],  # ATCG
        [[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]]
        #[[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]],  # GCTA
        #[[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]   # TGCA
    ])
    try:
        np.testing.assert_array_equal(reads_seq_binary_new, expected_binary)
        print("Test Passed: Binary representation is correct.")
    except AssertionError:
        print("Test Failed: Binary representation is incorrect.")

    # Check if weights are all 1 (unique sequences)
    try:
        np.testing.assert_array_equal(reads_weights_new, np.array([1,2]))
        print("Test Passed: Weights are all ones for unique sequences.")
    except AssertionError:
        print("Test Failed: Weights are not all ones for unique sequences.")

    # Check if qualities are correct
    expected_qualities = np.array([
        [40, 40, 40, 40],
        [25, 25, 25, 25],
        #[20, 20, 20, 20]
        ])
    try:
        np.testing.assert_array_equal(qualities_new, expected_qualities)
        print("Test Passed: Qualities are correct.")
    except AssertionError:
        print("Test Failed: Qualities are incorrect.")

    # Process with old function
    print("unique_modus=False")
    unique_modus=False

    # Process with new function
    reads_seq_binary_new, reads_weights_new, qualities_new, _ = load_and_process_reads(sample_fasta_file, sample_qualities_file, alphabet, unique_modus)

    # Additional checks
    try:
        assert reads_seq_binary_new.shape == (3, 4, 4)  # 3 sequences, 4 nucleotides, 4 alphabet letters
        print("Test Passed: Binary sequence shape is correct.")
    except AssertionError:
        print("Test Failed: Binary sequence shape is incorrect.")

    try:
        assert reads_weights_new.shape == (3,)
        print("Test Passed: Read weights shape is correct.")
    except AssertionError:
        print("Test Failed: Read weights shape is incorrect.")

    try:
        assert qualities_new.shape == (3, 4)
        print("Test Passed: Qualities shape is correct.")
    except AssertionError:
        print("Test Failed: Qualities shape is incorrect.")

    # Check if binary representation is correct
    expected_binary = np.array([
        [[1, 0, 0, 0], [0, 0, 0, 1], [0, 1, 0, 0], [0, 0, 1, 0]],  # ATCG
        [[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]],  # GCTA
        [[0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0]]  # GCTA
        #[[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]   # TGCA
    ])
    try:
        np.testing.assert_array_equal(reads_seq_binary_new, expected_binary)
        print("Test Passed: Binary representation is correct.")
    except AssertionError:
        print("Test Failed: Binary representation is incorrect.")

    # Check if weights are all 1 (unique sequences)
    try:
        np.testing.assert_array_equal(reads_weights_new, np.ones(3))
        print("Test Passed: Weights are all ones for unique sequences.")
    except AssertionError:
        print("Test Failed: Weights are not all ones for unique sequences.")

    # Check if qualities are correct
    expected_qualities = np.array([[40, 40, 40, 40], [30, 30, 30, 30], [20, 20, 20, 20]])
    try:
        np.testing.assert_array_equal(qualities_new, expected_qualities)
        print("Test Passed: Qualities are correct.")
    except AssertionError:
        print("Test Failed: Qualities are incorrect.")
