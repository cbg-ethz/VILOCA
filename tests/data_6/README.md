### Sample files to test `VILOCA`

Use files in this directory to test shorah in shotgun mode. The reads data have been generated with V-pipe's benchmarking framework (simulated with parameters: ```seq_tech~illumina__seq_mode~shotgun__seq_mode_param~nan__read_length~90__genome_size~90__coverage~100__haplos~5@5@10@5@10@geom@0.75```)

The reads are from one single amplicon of length 90, meaning the reference is of length 90 and each read is of length 90bps.

To use the new inference method using the sequencing quality scores use:
```
viloca run -f reference.fasta -b reads.shotgun.bam -w 90 --mode use_quality_scores
```
To use the model that is learning the sequencing error parameter:
```
viloca run -f reference.fasta -b reads.shotgun.bam -w 90 --mode -learn_error_params
```
To run VILOCA with the insert file run:
```
viloca run -f reference.fasta -b reads.shotgun.bam -w 90 --mode use_quality_scores -z scheme.insert.bed
```
