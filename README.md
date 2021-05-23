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

conda create -n nanoQC python=3.6
conda activate nanoQC
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

Note that nanoQC will look detect the key words 'pass' and 'fail' from absolute path of the file(s). If found, graph will provide the quality information for each group.

With fastq file(s) as input:
```
nanoQC_Guppy_v3.1.5.py -f /fastq/folder/ -o /output/folder/
```
With 'sequencing_summary.txt' as input:
```
nanoQC.py -s sequencing_summary.txt -o /output/folder/
```
