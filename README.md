# sequence_accessories
#### Accessory scripts for [`sequence_handling`](https://github.com/MorrellLab/sequence_handling)
___

## What is `sequence_accessories`?

`sequence_accessories` is a series of scripts that complement the function of `sequence_handling`. `sequence_handling` is designed to automate the processing of raw sequence data in the FASTQ format. Furthermore, `sequence_handling` also provides quality checks to allow the user to diagnose potential errors with the data or processing. However, there are other tasks that can be automated but fall outside the scope of `sequence_handling`. In order to keep `sequence_handling` focused, `sequence_accessories` was developed to handle these tasks.

## How does `sequence_accessories` work?

Much like `sequence_handling`, `sequence_accessories` is designed to process large amounts of data in parallel. However, `sequence_accessories` does **not** have an inherent dependence on the [Portable Batch System](http://www.pbsworks.com/). Instead of resource-heavy [handlers](https://github.com/MorrellLAB/sequence_handling#handlers), `sequence_accessories` uses lighter [accessories](https://github.com/MorrellLab/sequence_accessories#accessories) to do its processing. These accessories have fewer options than the handlers of `sequence_handling`. As such, `sequence_accessories` runs entirely on the command line without the help of a configuration file.

## Setting parameters

Unlike `sequence_handling`, `sequence_accessories` uses command line based arguments instead of parameters set in a configuration file. The format for arguments is `--parameter-name=value`. This is different than the normal `-p value` syntax found in most programs; however, the equals-delimeted syntax is easier to write a parser for. In addition to traditional arguments, some accessories also have flags. Flags follow the format `--flag`; to trigger the flags, simply pass it on the command line.

When reading documentation for `sequence_accessories`, required arguments are denoted by angular brackets around the value while optional arguments are denoted by square braces around the argument and value. All flags are optional, and are also dentoed by the square braces.

## Accessories

Basic usage is done with the following command:

```bash
./sequence_accessories <accessory> [options]
```

Where `<accessory>` one of the accessories listed below and `[options]` are arguments for the accessory. A simple help message can be found by running:

```bash
./sequence_accessories
```

Detailed help for each accessory can be found by running the accessory without any arguments.

### SummarizeStats

The SummarizeStats accessory runs [SAMTools](https://github.com/samtools/samtools) idxstats on BAM files and creates a text file with the sequence length, number of mapped reads, and number of unmapped reads for every sample passed to it. If the BAM files are not indexed, then SummarizeStats will generate CSI-format indexes for the BAM files.

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

The SRADownloader accessory downloads SRA archives from the SRA's FTP server. This accessory is given SR-/ER-/DR- numbers that correspond to an experiment, run, sample, or study.

#### Arguments:
 - `--sample-list=<sample_list>`: required list of SRA accession numbers to downloads
 - `--sample-type=<sample_type>`: required type of accession number given; can choose from: 'experiment', 'run', 'sample', or 'study'
 - `[--outdirectory=outdirectory]`: optional directory to download SRA archives to
 - `[--validate]`: flag to use vdb-validate to run a checksum on the SRA archives within SRADownloader

#### Dependencies:
 - [lftp](http://lftp.tech/)
 - [GNU Parallel](http://www.gnu.org/software/parallel/)
 - [vdb-validate](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=vdb-validate) from the [SRA Toolkit](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc) if validating within SRADownloader

## Future Accessories

### ListGenerator

The ListGenerator accessory will create a sample list for use with `sequence_handling` and `sequence_accessories`

### RegionalSNPs

The RegionalSNPs accessory will call SNPs by region using Freebayes

### FastSanger

The FastSanger accessory will convert Sanger sequencing data to FastQ format to be used with `sequence_handling`
