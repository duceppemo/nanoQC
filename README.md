https://github.com/duceppemo/nanoQC

# nanoQC

## Description
Quality control of Oxford Nanopore Technologies sequencing data basecalled with Guppy. Input can be a folder containing fastq file(s) or the 'sequencing_summary.txt' file generated by Guppy. Using the summary file is faster, but will not produce any graph about %GC content.

## Installation
```bash
# Set conda channels for Bioconda
conda config --add channels conda-forge
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels r

# Create nanoQC virtual environment
mamba create -n nanoQC -c conda-forge -c bioconda -c plotly -y python=3.11.4 numpy=1.25.0 python-dateutil=2.8.2 \
    cython=0.29.36 pandas=2.0.3 matplotlib=3.7.1 seaborn=0.12.2 plotly=5.15.0 python-dateutil=2.8.2 scikit-learn=1.3.0 \
    fpdf=1.7.2

# Activate nanoQC virtual environment
conda activate nanoQC

## Clone nanoQC repo
git clone https://github.com/duceppemo/nanoQC.git

# Compile the cython module
cd nanoQC
python setup.py build_ext --inplace

# Test nanoQC
python nanoQC.py -h
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

## Examples
Note that nanoQC will look detect the keywords 'pass' and 'fail' from absolute path of the file(s). If found, graph will provide the quality information for each group.

With fastq file(s) as input:
```bash
nanoQC.py -f /fastq/folder/ -o /output/folder/
```
With 'sequencing_summary.txt' as input:
```bash
nanoQC.py -s sequencing_summary.txt -o /output/folder/
```
## Ouputs
* A bunch of `png` files.
* A `pdf` file containing all the images.