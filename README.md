# nanoQC

TODO:
* Add pictures
* Create proper documentation
* Package properly
* more...

## Description
Quality control of Oxford Nanopore Technologies sequencing data basecalled qith Albacore. Works with fastq file(s) or the 'sequencing_summary.txt' file from Albacore.

## Installation
Install dependecies:
```
#install dependencies
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
nanoQC.py -s nanoQC_chunker_v2.py -o /output/folder/
```
