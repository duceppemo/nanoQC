# nanoQC

TODO:
* Add pictures
* Create proper documentation
* Package properly
* more...

## Description
Quality control of Oxford Nanopore Technologies sequencing data basecalled qith Albacore. Works with fastq file(s) or the 'sequencing_summary.txt' file from Albacore.

## Dependencies
This is a python3 program and depends on the following packages:
```
python3
numpy
pandas
matplotlib==3.1.0
seaborn==0.9.0
Cython
```
## Installation
Install dependecies:
```
# TODO -> create requirements.txt
```

Compile the cython module:
```
python setup.py build_ext --inplace
```

## Usage

Note that nanoQC will look detect the key words 'pass' and 'fail' from absolute path of the file(s). If found, graph will provide the quality information for each group.

With fastq file(s) as input:
```
nanoQC.py -f /fastq/folder/ -o /output/folder/
```
With 'sequencing_summary.txt' as input:
```
nanoQC.py -s sequencing_summary.txt -o /output/folder/
```
