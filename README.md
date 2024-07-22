VILOCA: VIral LOcal haplotype reconstruction and mutation CAlling for short and long read data
===============

VILOCA is an open source project for the analysis of next generation sequencing
data. It is designed to analyse genetically heterogeneous samples. Its tools
are written in different programming languages and provide error correction,
haplotype reconstruction and estimation of the frequency of the different
genetic variants present in a mixed sample.
VILOCA takes an alignment file as input, and subsequently generates mutation calls and local haplotypes.


---

### Installation
For installation miniconda is recommended: https://docs.conda.io/en/latest/miniconda.html.
We recommend to install VILOCA in a clean conda environment:
```
conda create --name env_viloca --channel conda-forge --channel bioconda libshorah
conda activate env_viloca
pip install git+https://github.com/cbg-ethz/VILOCA@master
```

### Example
To test your installation run VILOCA `tests/data_1`:
```
viloca run -b test_aln.cram -f test_ref.fasta -z scheme.insert.bed --mode use_quality_scores
```


If the sequencing amplicon strategy is known, we recommend using the amplicon-mode of the program, which takes as input the `<smth>.insert.bed` - file:
`viloca run -b test_aln.cram -f test_ref.fasta -z scheme.insert.bed --mode use_quality_scores`

If the sequencing quality scores are not trustable, the sequencing error parameters can also be learned:
`viloca run -b test_aln.cram -f test_ref.fasta -z scheme.insert.bed --mode learn_error_params`.

If there is no information on the sequencing amplicon strategy available, run:
`viloca run -b test_aln.cram -f test_ref.fasta --mode use_quality_scores`

### Parameters
There are several parameters available:  
`-b` [mandatory] sorted bam format alignment file  

`-f` [mandatory] reference genome in fasta format for mutation calling  

`-z` path to an (optional) insert file (primer tiling strategy), if available we highly recommend providing this file  

`--mode` mode to use:  
  - `learn_error_params`: model that is learning the error rate from the data  
  - `use_quality_scores`: model incorporating the sequencing quality scores that are passed through the alignment file  
  - `shorah`: use the tool ShoRAH (https://github.com/cbg-ethz/shorah)

`--extended_window_mode`: flag to call insertions (default: this flag is turned off)  

`--exclude_non_var_pos_threshold`: Percentage threshold for positions exclusion. Positions with base variations below this threshold will be excluded from the analysis, instead this position will be treated as if it only contains the reference base. This means that mutations of frequency < `exclude_non_var_pos_threshold` will not be called.

`--windowsize`: In case no insert file is provided, the genome is tiled into uniform local regions. `windowsize` determines the length of those local regions. It should be of roughly the length of the reads. This is also the length of the haplotypes that are produced.

### Output
`haplotypes` This directory contains the reconstructed local haplotypes as separate fasta files per local region.

`coverage.txt` List of each local region with start and end positions, and number of reads considered in the region.

`cooccurring_mutations.csv` All mutation calls including the information on which haplotype it occurred on.

## Development/CI with Docker
The following command in the root directory will let you interact with the project locally through Docker.
```bash
docker run --name viloca --rm -w="/usr/app" -it $(docker build -q .) bash
```
Run the following commands to copy the contents into the container and  install VILOCA inside Docker.
```bash
docker cp . viloca:/usr/app # run outside Docker
poetry install --only-root # run inside Docker
```

This is the same setup as used in the CI at [`.github/workflows/test.yaml`](.github/workflows/test.yaml).

## Profiling
```bash
poetry run python3 -m cProfile -m shorah shotgun ...
```
