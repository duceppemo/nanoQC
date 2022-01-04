https://github.com/duceppemo/nanoQC

# nanoQC

TODO:
* Add pictures
* Create proper documentation
* Package properly
* more...

## Description
Quality control of Oxford Nanopore Technologies sequencing data basecalled qith Albacore. Works with fastq file(s) or the 'sequencing_summary.txt' file from Albacore.

## Installation
```
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels r

conda create -n nanoQC python=3

conda install numpy python-dateutil cython pandas matplotlib seaborn

git clone https://github.com/duceppemo/nanoQC.git

cd nanoQC

# Compile the cython module
python setup.py build_ext --inplace
```
```
# TODO -> create requirements.txt
```

## Usage
```
usage: nanoQC.py [-h] [-f /basecalled/folder/] [-s sequencing_summary.txt] -o
                 /qc/ [-t 48] [-p 4] [-i]

Create QC plots using nanopore sequencing or basecalling data

optional arguments:
  -h, --help            show this help message and exit
  -f /basecalled/folder/, --fastq /basecalled/folder/
                        Input folder with fastq file(s), gzipped or not
  -s sequencing_summary.txt, --summary sequencing_summary.txt
                        The "sequencing_summary.txt" file produced by
                        guppy_basecaller
  -o /qc/, --output /qc/
                        Output folder
  -t 48, --threads 48   Number of CPU. Default 48
  -p 4, --parallel 4    Number of samples to process in parallel.
  -i, --individual      Produce a nanoQC report for each files instead of
                        combining them in a single report

```
Note that nanoQC will look detect the keywords 'pass' and 'fail' from absolute path of the file(s). If found, graph will provide the quality information for each group.

With fastq file(s) as input:
```
nanoQC.py -f /fastq/folder/ -o /output/folder/
```
With 'sequencing_summary.txt' as input:
```
nanoQC.py -s sequencing_summary.txt -o /output/folder/
```
