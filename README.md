ExplorATE - *Explore Active Transposable Elements* - 
====


## Install
____________________
With git installed, simply type the following command to install `ExplorATE`

```sh
git clone https://github.com/FemeniasM/ExplorATE
```
#### Requirements
- ExplorATE uses Seqkit v2.1.0 for linux 64bit, if you use another version please download a [Seqkit](https://bioinf.shenwei.me/seqkit/download/) version for your system and replace it in the `/bin` folder
- We added the Python3 script from RepeatMasker RM2Bed.py to speed up the overlap resolution, and now python 3 and Pandas python package is required.

**In `mo` (model organism) and `nmo` (non-model organism) modes require previously installed:**
 
- bedtools (tested with v2.30.0)
- AWK (tested with v4.1.4)
- Salmon (tested with v1.4.0)

**In `nmo_in` and `nmo_all` additionally require installed:**

- HMMER (tested with v3.3.2)
- BLAST+ (tested with v2.2.29+)
- TransDecoder (tested with v5.5.0)
- RepeatMasker (tested with v4.1.0)

And protein databases:

- protein database (tested with Swiss-Prot, accessed March 22, 2019)
- peptides for protein domains database (tested with Pfam-A v.33.0, accessed March 22, 2019)


## Overview
____________________
ExplorATE uses the Selective Alignment strategy (SA) to filter co-transcribed transposons with genes, based on alignment scores. ExplorATE first identifies target TEs and decoy TE sequences, and second performs the quantification of the target TEs using the SA algorithm in Salmon.
ExplorATE allows the TE analysis in multiple organisms (with or without reference genome). If a reference genome is provided, the user could (1) use a set of target TEs from the intergenic regions of the genome, or (2) use a *de novo* transcriptome (and its RepeatMasker annotations) to define the target TEs. Further, ExplorATE uses the reference genome and a genome-derived RepeatMasker file to extract TEs and define decoy sequences. Users can resolve overlapping repeats before extract target TEs.
If the reference genome is not available, ExplorATE uses the repeats identified by RepeatMasker from a *de novo* transcriptome. This transcriptome-derived RepeatMasker file is processed to resolve overlapping repeats and define target TEs and decoy sequences.
Target TEs can be defined at the fragment level or at the transcript level. At the fragment level, all non-co-transcribed TEs are used as target. At the transcript level, a rule similar to Wicker's is used to establish the identity of the transcript. The algorithm uses this criterion to assign target transcripts based on the percentage of identity for a class/family of TEs (calculated as 100 - the percentage of divergence), the percentage for each TE class/family in the transcript (calculated as the ratio between the TE length and the transcript length), and the minimum transcript length. For example, ‘80-80-80’ Wicker-like rule is a selection criterion where a transcript is annotated as "target" f it contains a TE class/family with percentage of identity >80%, this TE class/family represents >80% of the transcript length, and the transcript is at least 80bp in length.
ExplorATE creates files for Salmon execution with Selective Alignment algorithm. The counts estimated by Salmon can be imported into R with specific functions from the R package [ExplorATE](https://github.com/FemeniasM/ExplorATEproject).

## Citation

Please cite [ExplorATE: a new pipeline to explore active transposable elements from RNA-seq data](https://doi.org/10.1093/bioinformatics/btac354), *Bioinformatics*, Volume 38, Issue 13, 2022, Pages 3361–3366.

### How to use
____________________
This shell script provides four execution modes: one (`mo`) for model organisms with a reference genome and three (`nmo`,`nmo_in`,`nmo_all`) for non-model organisms as described below. To run the program, the user must define the mode and the flags corresponding to each mode: 
```sh
bash ExplorATE [mode][flags]
```
#### Model organisms

The `mo` mode performs the analysis for model organisms. The allowed flags are:
```
Usage:  
ExplorATE mo [flags]

Flags:
   -p threads [N] (1 default)     Number of therads
   -k kmer [N] (31 default)       k-mer size
   -c chromosome alias file       Replace the name of the chromosomes in the RepeatMasker file using
                                  a tab separated file with the first column indicating the desired
                                  chromosome name (e.g. the name of the gtf file) and in the second
                                  column the name to replace in the RepeatMasker file. If the file 
                                  contains more columns they will be ignored.
   -b bedtools binary path        Path to bedtools binary file. It is assumed by default that the 
                                  program is in your $PATH
   -s salmon binary path          Path to bedtools binary file (salmon default)
   -f fasta genome                Path to fasta genome (mandatory)
   -g gtf file                    Path to the gtf file for the genome version used (mandatory)
   -r RepeatMasker .out file      Path to the RepeatMasker file for the genome version used (mandatory) 
   -e library format ['pe']['se'] Indicates the format of the libraries: 'pe' for paired-end and 'se'
                                  for single-end reads. Supported extensions are: <.fq> or <.fastq>
                                  or <.fq.gz> or <.fastq.gz>. Paired end file names should contain 
                                  _R1 _R2. Example: sample_R1.fq.gz, sample_R2.fq.gz (mandatory)
   -l folder with fastq files     Path to folder with fastq files (mandatory)
   -o output directory path       Path to output directory (mandatory)
   -t fasta transcriptome **      Path to de novo transcriptome. Only required when target TEst are 
                                  based on the de novo transcriptome.
   -u transcriptome-dervied       Path to transcriptome-derived RepeatMasker file. Only required  
      RepeatMasker file **        when target TEst are based on the de novo transcriptome.
   -v overlap resolution          Criteria for overlapping resolution. Supported arguments:
                                  ['higher_score']['longer_element']['lower_divergence'] 
   -a .align file                 Alignments file derived from RepeatMasker (for genome) if defined
                                  'lower_divergence' as overlap resolution 
   -h help                        Print help
** if the TEs targets are based on a de novo transcriptome
```
##### Runing `mo` mode:
The following command defines target TEs from intergenic regions:
```sh
bash ExplorATE mo -p 12 -f genome_hs.fa -g genemodel_hs.gtf -r repmask_hs.out -e pe -l reads -o out_hs -v 'higher_score'
```
The user can explore repeats from a *de novo* transcriptome to define target TEs with the following command:
```sh
bash ExplorATE mo -p 12  -t trme_hs.fa -u repmask_trme_hs.out -f genome_hs.fa -g genemodel_hs.gtf -r repmask_hs.out -e pe -l reads -o out_hs -v 'higher_score'
```
Note that the above command incorporates the `-t` and `-u` arguments for the *de novo* transcriptome and its corresponding RepeatMasker file.

#### Non-model organisms

When there is no reference genome, ExplorATE uses a RepeatMasker output file derived from the transcriptome, and uses TransDecoder gene models and BLAST annotations to identify overlapping TEs with genes. The user can run ExplorATE with files previously made with the `nmo` mode.  
```
Usage:  
ExplorATE nmo [flags]

Flags:
   -p threads [N] (1 default)     Number of therads
   -k kmer [N] (31 default)       k-mer size
   -b bedtools binary path        Path to bedtools binary file. It is assumed by default that the 
                                  program is in your $PATH
   -e library format ['pe']['se'] Indicates the format of the libraries: 'pe' for paired-end and 'se'
                                  for single-end reads. Supported extensions are: <.fq> or <.fastq>
                                  or <.fq.gz> or <.fastq.gz>. Paired end file names should contain 
                                  _R1 _R2. Example: sample_R1.fq.gz, sample_R2.fq.gz (mandatory)
   -l folder with fastq files     Path to folder with fastq files (mandatory)
   -o output directory path       Path to output directory (mandatory)
   -t de novo transcriptome       Path to de novo transcriptome file (mandatory) 
   -s salmon binary path          Path to bedtools binary file (salmon default)
   -n gene annotation file        Path to gene annotation file from BLAST in output format 6 
   -d TransDecoder gff3 file      Path to .gff3 gene models file from TransDecoder (mandatory)
   -w Wicker-like rule            Comma separated values indicating respectively
                                  -Percentage of identity: calculated as 100 minus the percentage 
                                  of divergence (from RepeatMasker file) for each TE class/family
                                  -Percentage of length: ratio between TE class/family length with 
                                  respect to total the transcript length
                                  -minimum length of the transcript: minimum transcript length
                                  ('0,0,0' default)
   -v overlap resolution          Criteria for overlapping resolution. Supported arguments:
                                  ['higher_score']['longer_element']['lower_divergence'] 
                                  ('higher_score' default)
   -x split repeats by            Indicates if the target TEs will be annotated by ['name']['family']
                                  ['subclass']
   -q annotate target TEs by      Indicates if the target TE sequences will be fragments or whole transcripts
                                  ['transcripts']['fragments'] ('transcripts' default)
   -a .align file                 Alignments file derived from RepeatMasker (for genome) if defined
                                  'lower_divergence' as overlap resolution 
   -h help                        Print help
```

The `nmo` mode is recommended when there is no reference genome, since the user can supervise the programs execution for create input files. Alternatively the user can generate the input files with ExplorATE. There are two modes for this purpose, (1) the user can create only the input files with the `nmo_in` mode, or (2) the user can create the input files and run the entire ExplorATE pipeline from the shell script with `nmo_all` mode. The `nmo_in` mode can be used for subsequent ExplorATE execution in R and the supported flags are:

```
Usage:  
ExplorATE nmo_in [flags]

Flags:
   -p threads [N] (1 default)     Number of therads
   -k kmer [N] (31 default)       k-mer size
   -r RepeatMasker binary path    Path to RepeatMasker binary file
   -o output directory path       Path to output directory (mandatory)
   -t de novo transcriptome       Path to de novo transcriptome file (mandatory) 
   -n blastp binary path          Path to blastp binary file. It is assumed by default that the 
                                  program is in your $PATH
   -m hmmer binary path           Path to hmmscan binary file. It is assumed by default that the 
                                  program is in your $PATH (require hmmer3)
   -d TransDecoder directory path Path to TransDecoder directory
   -u protein database            Path to protein database (e.g. SwissProt or Uniref90)
   -f Pfam database               Path to Pfam database
   -i RepeatMasker library        Path to custom library of repeats. Defines the argument to 
                                  RepeatMasker -lib. It is recommended to use a specific library for the 
                                  organism of interest.
   -j RepeatMasker species        If a custom library is not defined, a species closely related to the organism 
                                  of interest can be used. This flag defines the -species RepeatMasker argument 
                                  and supports the same specifications. See detailed information in the 
                                  RepeatMasker documentation.
   -h help                        Print help
```
To run the entire ExplorATE pipeline for non-model organisms users can use the `nmo_all` mode. This execution mode generates intermediate files and performs the quantification estimation with Salmon. The supported flags are:

```
Usage:  
ExplorATE nmo_in [flags]

Flags:
   -p threads [N] (1 default)     Number of therads
   -k kmer [N] (31 default)       k-mer size
   -e library format ['pe']['se'] Indicates the format of the libraries: 'pe' for paired-end and 'se'
                                  for single-end reads. Supported extensions are: <.fq> or <.fastq>
                                  or <.fq.gz> or <.fastq.gz>. Paired end file names should contain 
                                  _R1 _R2. Example: sample_R1.fq.gz, sample_R2.fq.gz (mandatory)
   -l folder with fastq files     Path to folder with fastq files (mandatory)
   -b bedtools binary path        Path to bedtools binary file. It is assumed by default that the 
                                  program is in your $PATH
   -r RepeatMasker binary path    Path to RepeatMasker binary file
   -o output directory path       Path to output directory (mandatory)
   -t de novo transcriptome       Path to de novo transcriptome file (mandatory) 
   -n blastp binary path          Path to blastp binary file. It is assumed by default that the 
                                  program is in your $PATH
   -m hmmer binary path           Path to hmmscan binary file. It is assumed by default that the 
                                  program is in your $PATH (require hmmer3)
   -d TransDecoder directory path Path to TransDecoder directory
   -u protein database            Path to protein database (e.g. SwissProt or Uniref90)
   -f Pfam database               Path to Pfam database
   -s salmon binary path          Path to bedtools binary file (salmon default)
   -i RepeatMasker library        Path to custom library of repeats. Defines the argument to 
                                  RepeatMasker -lib. It is recommended to use a specific library for the 
                                  organism of interest.
   -j RepeatMasker species        If a custom library is not defined, a species closely related to the organism 
                                  of interest can be used. This flag defines the -species RepeatMasker argument 
                                  and supports the same specifications. See detailed information in the 
                                  RepeatMasker documentation.
   -w Wicker-like rule            Comma separated values indicating respectively
                                  -Percentage of identity: calculated as 100 minus the percentage 
                                  of divergence (from RepeatMasker file) for each TE class/family
                                  -Percentage of length: ratio between TE class/family length with 
                                  respect to total the transcript length
                                  -minimum length of the transcript: minimum transcript length
                                  ('0,0,0' default)
   -v overlap resolution          Criteria for overlapping resolution. Supported arguments:
                                  ['higher_score']['longer_element']['lower_divergence'] 
                                  ('higher_score' default)
   -x split repeats by            Indicates if the target TEs will be annotated by ['name']['family']
                                  ['subclass']
   -q annotate target TEs by      Indicates if the target TE sequences will be fragments or whole transcripts
                                  ['transcripts']['fragments'] ('transcripts' default)
   -a .align file                 Alignments file derived from RepeatMasker (for genome) if defined
                                  'lower_divergence' as overlap resolution 
   -h help                        Print help
```

## Outputs
____________________

In the output directory a folder `quant_out` is created with the Salmon estimates and the reference file `references.csv`. These files are used to import estimates into R with functions from ExplorATE [package](https://github.com/FemeniasM/ExplorATEproject). Further, ExplorATE writes a tab separated file `repeats_in_UTRs.txt` with the transcripts that contain repeats in the UTR regions. The first column indicates the name of the transcript, the second and third columns indicate the position of the UTR region in the transcript (start and end), the fourth column indicates the type of feature (5'-UTR or 3'-UTR) , the fifth and sixth columns indicate the position of the repetition in the transcript (start and end), and finally the last column indicates the name, class and family of the repetition separated by colons (":", for example, "Plat_L3 : LINE: CR1 ").

## Examples
____________________
Users can refer to [vignette](https://femeniasm.github.io/ExplorATE_vignette/) and [user guide](https://femeniasm.github.io/ExplorATE_user_guide/) to run [test data](https://github.com/FemeniasM/ExplorATE_data_test). Intuitive examples for each execution modes of the shell script are given below:

Running `mo` mode
```sh
bash ExplorATE mo -p 12 -f genome.fa -g genemodel.gtf -r repmask.out -e pe \\
-l reads_folder -o out_folder -v 'higher_score'
```
Running `mo` mode with *de novo* transcriptome
```sh
bash ExplorATE mo -p 12 -f genome.fa -g genemodel.gtf -r repmask.out -e pe \\
-l reads_folder -o out_folder -v 'higher_score' -t trme -u repmask_trme.out
```
Running `nmo` mode
```sh
bash ExplorATE nmo -p 12 -b path/to/bedtools -s path/to/salmon -e pe \\
-l reads_folder -o out_folder -t trme.fa -r repmask_trme.out \\
-n blast_out.outfmt6 -d TransDecoder_out.gff3 -w 80,80,80 -v 'higher_score' \\
-x subclass -q transcripts
```
Running `nmo_in` mode
```sh
bash ExplorATE nmo_in -p 12 -n <blastp binary path> -m <hmmerscan binary path> \\
-r <RepeatMasker binary path> -d <TransDecoder directory path> \\
-u <SwissProt database> -f <Pfam database> -i <user-defined TE library> \\
-t trme.fa -o inputs_to_ExplorATE_pip
```
