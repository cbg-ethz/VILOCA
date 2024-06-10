VILOCA: VIral LOcal haplotype reconstruction and mutation CAlling for short and long read data
===============

VILOCA is an open source project for the analysis of next generation sequencing
data. It is designed to analyse genetically heterogeneous samples. Its tools
are written in different programming languages and provide error correction,
haplotype reconstruction and estimation of the frequency of the different
genetic variants present in a mixed sample.

The corresponding manuscript can be found here: https://www.biorxiv.org/content/10.1101/2024.06.06.597712v1

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

### Applications

You can find several applications of VILOCA at https://github.com/cbg-ethz/viloca_applications.
