VILOCA: VIral LOcal haplotype reconstruction and mutation CAlling for short and long read data
===============

VILOCA is an open source project for the analysis of next generation sequencing
data. It is designed to analyse genetically heterogeneous samples. Its tools
are written in different programming languages and provide error correction,
haplotype reconstruction and estimation of the frequency of the different
genetic variants present in a mixed sample.
VILOCA takes an alignment file as input, and subsequently generates mutation calls and local haplotypes.


The corresponding manuscript can be found at: https://academic.oup.com/nargab/article/6/4/lqae152/7912062

Fuhrmann L, Langer B, Topolsky I, Beerenwinkel N. VILOCA: sequencing quality-aware viral haplotype reconstruction and mutation calling for short-read and long-read data. NAR Genomics and Bioinformatics. 2024 Dec;6(4):lqae152.

---

### Installation
For installation miniconda is recommended: https://docs.conda.io/en/latest/miniconda.html.
We recommend to install VILOCA in a clean conda environment:
```
conda create --name env_viloca --channel conda-forge --channel bioconda viloca
conda activate env_viloca
```

If you want to install the `master` branch use:
```
conda create --name env_viloca --channel conda-forge --channel bioconda libshorah
conda activate env_viloca
pip install git+https://github.com/cbg-ethz/VILOCA@master
```

### Example
To test your installation run VILOCA `tests/data_1`:
```
viloca run -b test_aln.cram -f test_ref.fasta --mode use_quality_scores
```

Another example can be found in  `tests/data_6`:
If the sequencing amplicon strategy is known, we recommend using the amplicon-mode of the program, which takes as input the `<smth>.insert.bed` - file:
`viloca run -f reference.fasta -b reads.shotgun.bam -w 90 --mode use_quality_scores -z scheme.insert.bed`

If there is no information on the sequencing amplicon strategy available, run:
`viloca run -f reference.fasta -b reads.shotgun.bam -w 90 --mode use_quality_scores`

If the sequencing quality scores are not trustable, the sequencing error parameters can also be learned:
`viloca run -f reference.fasta -b reads.shotgun.bam -w 90 --mode learn_error_params`.


### Parameters
There are several parameters available:  
`-b` [mandatory] Input alignment file in sorted BAM format. The file must contain only primary alignments. Consider using a processing pipeline that filters out secondary and supplementary alignments.   

`-f` [mandatory] Reference genome in fasta format for mutation calling.  

`--mode` mode to use:  
  - `learn_error_params`: model learning the error rate from the data  
  - `use_quality_scores`: model incorporating the sequencing quality scores that are passed through the alignment file [recommend option]
  - `shorah`: use the tool ShoRAH (https://github.com/cbg-ethz/shorah)

`--windowsize` In case no insert file is provided, the genome is tiled into uniform local regions. `windowsize` determines the length of those local regions. It should be of roughly the length of the reads. This is also the length of the haplotypes that are produced. Any read that covers less than the minimum percentage of the local region—defined by the `--win_min_ext` parameter — will be excluded from the analysis.

`-z` Path to an file that defines the local regions used to segment the alignment (e.g.,`tests/data_1/scheme.insert.bed`). If this file is not specified, the alignment will be segmented into uniform regions of length defined by the `--windowsize` parameter. Any read that covers less than the minimum percentage of the local region—defined by the `--win_min_ext` parameter — will be excluded from the analysis.

`--win_min_ext` Minimum percentage of bases to overlap between reference and read to be considered in a window (default: 0.85). The rest (i.e. non-overlapping part) will be filled with Ns.

`-p`: Posterior threshold (default: 0.9) when calling variants from haplotypes.

`--extended_window_mode` Flag to call insertions (default: this flag is turned off)  

`--exclude_non_var_pos_threshold` Percentage threshold for positions exclusion. Positions with base variations below this threshold will be excluded from the analysis, instead this position will be treated as if it only contains the reference base. This means that mutations of frequency < `exclude_non_var_pos_threshold` will not be called.


### Output
`haplotypes` This directory contains the reconstructed local haplotypes as separate fasta files per local region.

`coverage.txt` List of each local region with start and end positions, and number of reads considered in the region.

`cooccurring_mutations.csv` The file contains one row per occurrence of a mutation in each haplotype, listing all haplotypes where a mutation is present. Note that the posterior threshold is not applied here.

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
poetry run python3 -m cProfile -m viloca run ...
```

### Applications

You can find several applications of VILOCA at:
- https://github.com/cbg-ethz/viloca_applications
- https://github.com/cbg-ethz/DCV-CrPV-cGAS-STING-pathway-data-analysis
