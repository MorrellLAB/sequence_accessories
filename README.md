# sequence_accessories
#### Accessory scripts for [`sequence_handling`](https://github.com/MorrellLab/sequence_handling)
___

## What is `sequence_accessories`?

`sequence_accessories` is a set of scripts that complements `sequence_handling`. `sequence_handling` is designed to automate the processing of raw sequence data in FASTQ format. It also provides quality checks so users can diagnose potential errors in the data or processing. However, other tasks can be automated but fall outside the scope of `sequence_handling`. To keep `sequence_handling` focused, `sequence_accessories` was developed to handle those tasks.

## How does `sequence_accessories` work?

Much like `sequence_handling`, `sequence_accessories` is designed to process large amounts of data in parallel. However, `sequence_accessories` does **not** have an inherent dependence on the [Portable Batch System](http://www.pbsworks.com/). Instead of resource-heavy [handlers](https://github.com/MorrellLAB/sequence_handling#handlers), `sequence_accessories` uses lighter [accessories](https://github.com/MorrellLab/sequence_accessories#accessories) to do its processing. These accessories have fewer options than the handlers of `sequence_handling`. As such, `sequence_accessories` runs entirely on the command line without the help of a configuration file.

## Setting parameters

Unlike `sequence_handling`, `sequence_accessories` uses command-line arguments instead of parameters set in a configuration file. The format for arguments is `--parameter-name=value`. This differs from the common `-p value` syntax found in most programs; however, the equals-delimited syntax is easier to parse. In addition to arguments, some accessories also have flags. Flags follow the format `--flag`; to trigger a flag, simply pass it on the command line.

When reading documentation for `sequence_accessories`, required arguments are denoted by angle brackets around the value, while square brackets around the argument and value denote optional arguments. All flags are optional and are also denoted by square brackets.

## Accessories

Basic usage is done with the following command:

```bash
./sequence_accessories <accessory> [options]
```

Where `<accessory>` is one of the accessories listed below, and `[options]` are arguments for that accessory. A simple help message can be found by running:

```bash
./sequence_accessories
```

Detailed help for each accessory is available by running it without arguments.

### Available accessories

 - `ListGenerator` *(currently not ready; exits in dispatcher)*
 - `SummarizeStats`
 - `DumpFastq`
 - `SRADownloader`
 - `MergeBAM`
 - `SimpleCoverage`
 - `AddBAMLane`
 - `FreeBayes_SNP_Calls`
 - `SubsampleFastq`
 - `UG100_filter`
 - `PanDepthCoverage`

### SummarizeStats

The SummarizeStats accessory runs [SAMTools](https://github.com/samtools/samtools) idxstats on BAM files and creates a text file with the sequence length, number of mapped reads, and number of unmapped reads for every sample passed to it. If the BAM files are not indexed, SummarizeStats will generate CSI-format indexes for them.

#### Arguments:
 - `--sample-list=<sample_list>`: required list of BAM files to generate statistics for
 - `[--project=project]`: optional name for the output file, defaults to 'STATS'

#### Dependencies:
 - [SAMTools](http://www.htslib.org/) 1.3 or higher
 - [GNU Parallel](http://www.gnu.org/software/parallel/)

### DumpFastq

The DumpFastq accessory creates gzipped FASTQ files from SRA archives. This accessory can handle dumping to either single- or paired-end FASTQ files.

#### Arguments:
 - `--sample-list=<sample_list>`: required list of SRA archives to dump
 - `[--outdirectory=outdirectory]`: optional directory to dump the FASTQ files to
 - `[--paired]`: flag to dump to paired-end FASTQ files

#### Dependencies:
 - [fastq-dump](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=fastq-dump) from the [SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
 - [GNU Parallel](http://www.gnu.org/software/parallel/)

### SRADownloader

The SRADownloader accessory downloads SRA archives from the SRA FTP server. This accessory takes SR-/ER-/DR- numbers that correspond to an experiment, run, sample, or study.

#### Arguments:
 - `--sample-list=<sample_list>`: required list of SRA accession numbers to download
 - `--sample-type=<sample_type>`: required type of accession number given; can choose from: 'experiment', 'run', or 'study'
 - `[--outdirectory=outdirectory]`: optional directory to download SRA archives to
 - `[--validate]`: flag to use vdb-validate to run a checksum on the SRA archives within SRADownloader

#### Dependencies:
 - [lftp](http://lftp.tech/)
 - [GNU Parallel](http://www.gnu.org/software/parallel/)
 - [vdb-validate](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=vdb-validate) from the [SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc) if validating within SRADownloader

### MergeBAM

The MergeBAM accessory uses BAMtools to merge several BAM files into a single BAM file. This accessory can handle multiple merges at once.

#### Arguments
 - `--sample-list=<sample_list>`: required list of BAM files
 - `--name-table=<table>`: a table where the first column is the sample name for the merged BAM file, and the remaining columns are the names of BAM files that make up the merged BAM
 - `[--outdirectory=outdirectory]`: optional directory to place the merged BAM

#### Dependencies
 - [BAMtools](https://github.com/pezmaster31/bamtools)
 - [GNU Parallel](http://www.gnu.org/software/parallel/)

### SimpleCoverage

The SimpleCoverage accessory uses SAMTools to calculate coverage over BAM files. This accessory outputs a table summarizing the average depth across the entire sample.

#### Arguments
 - `--sample-list=<sample_list>`: required list of BAM files
 - `[--genome-size=genome_size]`: optional size of the genome in number of base pairs, will calculate automatically if not specified; if you have exome sequencing data, provide the exome size
 - `[--project=project]`: optional name for the output file, defaults to 'SimpleCoverage'
 - `[--outdirectory=outdirectory]`: optional directory to place output files

#### Dependencies
 - [awk](http://www.cs.princeton.edu/~bwk/btl.mirror/)
 - [SAMTools](http://www.htslib.org/)
 - [GNU Parallel](http://www.gnu.org/software/parallel/)

### AddBAMLane

Adds lane/read-group information to BAM files.

#### Dependencies
 - [SAMTools](http://www.htslib.org/)
 - [GNU Parallel](http://www.gnu.org/software/parallel/)

### FreeBayes_SNP_Calls

Runs a FreeBayes-based SNP calling workflow.

#### Notes
 - This accessory script currently contains internal, hard-coded paths and parameters.

#### Dependencies
 - [GNU Parallel](http://www.gnu.org/software/parallel/)
 - [FreeBayes](https://github.com/freebayes/freebayes)
 - [SAMTools](http://www.htslib.org/)
 - [BAMTools](https://github.com/pezmaster31/bamtools)
 - `ogap`
 - `bamleftalign`

### SubsampleFastq

Subsample FASTQ reads using `seqtk`.

#### Dependencies
 - [seqtk](https://github.com/lh3/seqtk)
 - [GNU Parallel](http://www.gnu.org/software/parallel/)

### UG100_filter

Runs quality filtering for the UG100 cohort variant calls.

#### Arguments
 - `INPUT_FILE`: input BCF/VCF file
 - `OUT_DIR`: output directory

#### Dependencies
 - [bcftools](https://samtools.github.io/bcftools/) 1.21 or higher

### PanDepthCoverage

Calculates coverage depth over BAM/CRAM files using [PanDepth](https://github.com/HuiyangYu/PanDepth). Supports whole-genome, gene-level (GFF/GTF), region-level (BED), and windowed coverage modes. The `--gff`, `--bed`, and `--window-size` options are mutually exclusive; if none is given, per-chromosome statistics are reported.

#### Arguments
 - `--sample-list=<sample_list>`: required list of BAM or CRAM files
 - `[--gff=gff_file]`: optional GFF/GTF file for gene-level coverage
 - `[--bed=bed_file]`: optional BED file for region-level coverage
 - `[--window-size=window_size]`: optional window size in bp for windowed coverage
 - `[--feature=feature_type]`: GFF/GTF feature to parse — `CDS` or `exon` (default: `CDS`)
 - `[--min-mapq=min_mapq]`: minimum mapping quality filter (default: `0`)
 - `[--threads=threads]`: PanDepth threads per sample (default: `3`)
 - `[--project=project]`: optional output file prefix (default: `PanDepthCoverage`)
 - `[--outdirectory=outdirectory]`: optional output directory

#### Dependencies
 - [PanDepth](https://github.com/HuiyangYu/PanDepth) v2.26 or higher
 - [GNU Parallel](http://www.gnu.org/software/parallel/)

## Future Accessories

<!--### RegionalSNPs

The RegionalSNPs accessory will call SNPs by region using Freebayes-->

### FastSanger

FastSanger functionality is currently provided by an external utility script:

- `/Users/pmorrell/Library/CloudStorage/Dropbox/Documents/Sandbox/PMorrell/Utilities/phd_to_fastq.py`

That script converts `.phd.1` files to FASTQ format for use with `sequence_handling`.
