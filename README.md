# asqcan
## A combined pipeline for bacterial genome assembly, quality control and annotation

asqcan is a workflow pipeline for the automated assembly, quality control and annotation of bacterial genome sequences. Modern bacterial sequencing projects can involve a significant number of isolates and the process of assembling, running necessary QC and annotation can be time consuming. The asqcan pipeline seeks to automate this as much as possible. The current steps asqcan takes are:

1. Quality analysis of raw reads with FastQC
2. Genome assembly with spades
3. Quality analysis of assemblies with quast
4. Contamination and quality analysis of assemblies with blobtools
5. Annotation of assemblies using prokka

The asqcan pipleine runs these five steps on each .fastq or .fastq.gz reads file in the directory provided by the -i option. When asqcan completes, it generates a report on the success or failure of each step of the pipline (asqcan_rport.tsv). Successful steps will not be rerun on a subsequent execution, i.e. asqcan will detect successful steps and ignore them in future runs. 

###  Requirements

asqcan requires a linux-based system and the following:

- python (2.7)
- [GNU parallel](https://www.gnu.org/software/parallel/) (>=20170422)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (>=0.11.7)
- [spades](http://cab.spbu.ru/software/spades/) (>=3.11.1)
- [quast](http://quast.sourceforge.net/quast) (>=4.6.3)
- [blobtools](https://github.com/DRL/blobtools) (>=1.0)
- [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (>=2.7.1)
- [prokka](https://github.com/tseemann/prokka) (>=1.13)

### Installation

To download and install asqcan with all dependencies use conda:
```
conda install -c conda-forge -c bioconda asqcan
```

or pip (requires manual dependency installation): 
```
pip install git+https://github.com/bogemad/asqcan.git
```

or a manual install (again this requires you to manually install dependencies):
```
git clone https://github.com/bogemad/asqcan.git
cd asqcan
python setup.py install
```

### Usage

```
usage: asqcan [-h] -q READS_DIR -o OUTDIR -b DB [-t THREADS] [-m MEM] [-f]
              [--version] [-v]

required arguments:
  -q READS_DIR, --fastq-dir READS_DIR
                        Path to a directory with your interleaved fastq files.
  -o OUTDIR, --output-directory OUTDIR
                        Path to the output directory. A directory will be
                        created if one does not exist.
  -b DB, --blast_database DB
                        Path to the local nt blast database. This pipeline
                        requires you to download a copy of this database. See
                        https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_
                        TYPE=BlastDocs&DOC_TYPE=Download

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use for multiprocessing.
  -m MEM, --max_memory MEM
                        Maximum amount of RAM to assign to the pipeline in GB
                        (Just the number).
  -f, --force           Overwrite files in the output directories.
  --version             show program's version number and exit
  -v, --verbose         Increase verbosity on command line output (n.b.
                        verbose output is always saved to asqcan.log in the
                        output directory).
```
