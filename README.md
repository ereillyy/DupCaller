# DupCaller

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://github.com/YuheCheng62/DupCaller/blob/main/README.md) [![License](https://img.shields.io/badge/License-BSD%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause) [![Build Status](https://travis-ci.com/YuheCheng62/DupCaller.svg?branch=main)](https://app.travis-ci.com/YuheCheng62/DupCaller)[![Uptime Robot status](https://img.shields.io/uptimerobot/status/m795312784-02766a79f207f67626cef289)](https://stats.uptimerobot.com/jjqW4Ulymx)

DupCaller is a universal tool for calling somatic mutations and calculating somatic mutational burden from barcoded error-corrected next generation sequencing (ecNGS) data with matched normal (e.x. NanoSeq).

## Prerequisites
DupCaller requires python>=3.10 to run. Earlier versions may be sufficient to run DupCaller but have not been tested.
The complete DupCaller pipeline also requires the following tools for data preprocessing. The versions are used by the developer and other versions may or may not work.

- BWA version 0.7.17 (https://bio-bwa.sourceforge.net)
- GATK version 4.2.6 (https://github.com/broadinstitute/gatk/releases)
- Tabix for indexing compressed genomic files (recommended installation: `conda install bioconda::tabix`)

## Installation
The tool uses pip for installing scripts and prerequisites. To install DupCaller, simply clone this repository and install via pip:

```bash
git clone https://github.com/YuheCheng62/DupCaller.git
cd DupCaller
pip install .
```

## Pipeline
### Index reference genomes:
DupCaller uses a numpyrized reference genome to perform memory-efficient reference fetching and trinucleotide context indexing. To index the reference genome, run:
```bash
DupCaller.py index -f reference.fa
```
The command will generate three h5 files in the same folder: ref.h5, tn.h5 and hp.h5, which are numpyrized reference sequences, trinucleotide contexts and homopolymer lengths, respectively. Make sure that when running other DupCaller utilities, the three files are within the same folder as the reference genome. 


### Trim barcodes from reads:

DupCallerTrim.py is a script that can extract 5-prime barcodes from paired-end fastqs. The usage is as follows:

```bash
DupCaller.py trim -i read1.fq -i2 read2.fq -p barcode_pattern -o sample_name
```

where

'read1.fq' and 'read2.fq' are fastq files from read1 and read2 of the paired-end sequencing data, respectively. Both unzipped and gzip compressed files can be correctly processed.

'barcode_pattern' is the pattern of barcodes starting from the 5-prime end, with N representing a barcode base and X representing a skipped base. The notation is similar to what has been used in UMI-tools(https://github.com/CGATOxford/UMI-tools). For example, NanoSeq uses 3-base barcode followed by 4 constant bases, so the pattern should be NNNXXXX.

'sample_name' is the prefix of output paired fastqs. After the run completes, '{sample_name}\_1.fastq' and '{sample_name}\_2.fastq' will be generated. The barcodes will be recorded in each read name as {original_read_name}:{read1_barcode}+{read2_barcode}

If the matched normal is prepared in the same way as the sample, also apply the trimming with the same scheme to the matched normal fastqs. For traditional bulk normal, trimming is not needed.

### Align trimmed fastqs

Use a DNA NGS aligner, such as BWA-MEM, to align the trimmed fastqs of both sample and matched normal from the last step. Notice that GATK requires read group ID,SM and PL to be set, so adding those tags during bwa alignment is recommended. **Fastq tags must be kept. For bwa mem, this is using the -C option.** For example:

```bash
bwa mem -C -t {threads} -R "@RG\tID:{sample_name}\tSM:{sample_name}\tPL:ILLUMINA" reference.fa {sample_name}\_1.fastq {sample_name}\_2.fastq | samtools sort -@ {threads} > {sample_name}.bam
samtools index -@ {threads} {sample_name}.bam
```

where

'threads' is the number of cores used for aligning

'reference.fa' is the reference genome fasta file

'{sample_name}\_1.fastq' and '{sample_name}\_2.fastq' are trimmed fastq files from the last step.

### MarkDuplicates with optical duplicates tags and new read name configuration

Run GATK MarkDuplicates on sample and matched-normal bams. Notice that optical duplicates and PCR duplicates should be treated differently in ecNGS variant calling, so the "TAGGING_POLICY" of GATK MarkDuplicates should be set to OpticalOnly to differentiate optical duplicates from PCR duplicates. Additionally, the DUPLEX_UMI option should be set to true, and since the read name of trimmed fastq is modified, the READ_NAME_REGEX option should also be set to "(?:.*:)?([0-9]+)[^:]_:([0-9]+)[^:]_:([0-9]+)[^:]\_$". **Outdated versions of GATK do not have the --DUPLEX_UMI tag. Please update GATK to the latest version if this happens.** The MarkDuplicates command should look like this:

```bash
gatk MarkDuplicates -I sample.bam -O sample.mkdped.bam -M sample.mkdp_metrics.txt --READ_NAME_REGEX "(?:.*:)?([0-9]+)[^:]*:([0-9]+)[^:]*:([0-9]+)[^:]*$" --DUPLEX_UMI --TAGGING_POLICY OpticalOnly --BARCODE_TAG DB
```

### Variant Calling

After appropriate data preprocessing, DupCaller.py should be used to call somatic mutations. The usage depends on your experimental design, and here are some examples.

1. Whole genome/whole exome/reduced genome (e.g. NanoSeq) with a matched normal

```bash
DupCaller.py call -b ${sample}.bam -f reference.fa -o {output_predix} -p {threads} -n {normal.bam} -g germline.vcf.gz -m noise_mask.bed.gz
```

2. Mutagenesis panel without a matched normal

```bash
DupCaller.py call -b ${sample}.bam -f reference.fa -o {output_predix} -p {threads} -g germline.vcf.gz -m noise_mask.bed.gz -maf 0.1
```

Note that since DupCaller partitions jobs based on genomic regions, multithreading capability will be significantly compromised for small targeted panels. In this case, we suggest running at most one thread per distinct targeted region. 

**Input validation**: DupCaller now includes comprehensive input file validation that checks for the existence of all required files (BAM files, reference genome, h5 index files, and optional files) before starting analysis, providing clear error messages for missing files.

**Multi-threading improvements**: Recent updates include enhanced coverage file handling for multi-threaded execution, with automatic merging of adjacent region coverage files and improved memory efficiency. Coverage files from different threads are automatically combined post-processing with region-aware boundary detection.

Please see "Parameters" section for explanation of all parameters. See "Results" section for descriptions of all result files in the output folder

#### Parameters

**Required**

These options are required for each run.

| short option | long option | description |
| --- | --- | --- |
| -b | --bam | bam file of ecNGS data |
| -f | --reference | reference genome fasta file |
| -o | --output | prefix of the output files |

**Recommended**

These options should be understood by the user and customized accordingly. Some of them involve resources that should be used when available. All resources for GRCh38/hg38 and GRCm39/mm39 are provided at ... and should be used for the matching reference genome. 

| short option | long option | description | default |
| --- | --- | --- | --- |
| -r | --regions | contigs to consider for variant calling. The default is set for human. For any other species, please set the contigs accordingly. For example, for mouse, please set to "-r chr{1..19} chrX chrY" | default: chr{1..22} chrX chrY |
| -g | --germline | indexed germline vcf with AF field. | None |
| -p | --threads | number of threads | 1 |
| -n | --normalBam | a list of bam files of matched normals. When matched normal is not available, set the maximum allele frequency (-maf) to an appropriate value (e.g. 0.1) | None |
| -m | --noise | a list of bed interval files that mask noisy positions | None |
| -R | --regionfile | an inclusive bed file for specifying target regions | None |
| -maf | --maxAF | maximum allele fraction to call a somatic mutation. Must be set to an appropriate value when a matched normal (-n) is not available | 1 |
| -tt | --trimF | ignore mutation if it is less than n bps from the ends of template | 7 |
| -tr | --trimR | ignore mutation if it is less than n bps from the ends of read | 7 |
| -id | --indelBed | an indel enhanced panel of normal (ePoN) used for indel calling | None | 
| --naf | |maximum VAF in matched normal for a mutation to be called | 0.01 |
| --rescue | | output discarded variants with reason in the filter field | False |
| -nm | --nmflt | filter out any reads that have an edit distance larger than this value | 5 |
| -ax | --minMeanASXS | minimum mean AS-XS alignment score difference for a read group to be considered for calling (adjustable parameter for read quality filtering) | 50 |
| -gaf | --germlineAfCutoff | locations at which there is a germline mutation with population af larger than this threshold will be skipped | 0.001 |
| -d | --minNdepth | minimum coverage in normal for called variants | 10 |
**Advanced**

These are mostly variant calling model parameters and adjustment is unnecessary for general use.

| short option | long option | description | default |
| --- | --- | --- | --- |
| -AS | --amperrfile | pre-learned error profile for amplification SBS error | None | 
| -AI | --amperrfileindel | pre-learned error profile for amplification indel error | None | 
| -DS | --dmgerrfile | pre-learned error profile for SBS damage | None | 
| -DI | --dmgerrfileindel | pre-learned error profile for indel damage | None | 
| -mr | --mutRate | prior somatic mutation rate per base | 2.5e-7 |
| -ts | --thresholdSnv | log likelihood ratio threshold for SNV mutation calls | 0.5 |
| -ti | --thresholdIndel | log likelihood ratio threshold for indel mutation calls | 0.5 |
| -mq | --mapq | minimum mapq for an alignment to be considered | 40 |
| -w | --windowSize | genomic window size when calculating rough coverage and split bam files into equal regions. Adjust for smaller panel | 100000 |
| -bq | --minBq | bases with quality less than this number will be set to 6 | 18 |
| -aq | --minAltQual | minimum consensus quality of alt allele, if not 0, in a read group to be considered for training | 60 |
| --minRef | minimum consensus quality of ref allele, if not 0, in a read group to be considered for training | 2 |
| --minAlt | minimum consensus quality of alt allele, if not 0, in a read group to be considered for training | 2 |


### Mutation burden estimation
After mutation calling, mutational burden estimation can be performed within the folder:
```bash
DupCaller.py estimate -i sample -f reference.fa -r chr{1..22} chrX
```
Adjust the region according to the reference genome used. The estimation process now includes improved column ordering in output files and enhanced SBS (single base substitution) burden calculations with corrected trinucleotide context analysis.

#### Mean duplex coverage for each gene
For dNdScv coverage correction, DupCaller can output mean duplex depth for each gene, given a gene bed for option -gb:
```bash
DupCaller.py estimate -i sample -f reference.fa -r chr{1..22} chrX -gb {target}.bed
```
The gene bed should have the fourth column as {gene_name}_{exon_number} (e.g. tp53_1).

#### Re-estimate mutational burden of a given region:
To re-estimate trinucleotide-corrected mutational burden in certain regions without re-running the variant calling process, use the -rb option to provide the bed file for the desired region:
```bash
DupCaller.py estimate -i sample -f reference.fa -r chr{1..22} chrX -rb {re_estimate}.bed
```

#### Burden Estimation Parameters

**Required**

| short option | long option | description |
| --- | --- | --- |
| -i | --prefix | Input prefix of results from call command |

**Reference Options (choose one)**

| short option | long option | description | default |
| --- | --- | --- | --- |
| -f | --reference | Fasta file of reference genome | None |
| -ft | --refTrinuc | Precomputed trinucleotide composition of reference genome | None |

**Optional**

| short option | long option | description | default |
| --- | --- | --- | --- |
| -r | --regions | contigs to consider for trinucleotide calculation | chr{1..22} chrX |
| -ot | --outTrinuc | Output the computed trinucleotide composition file for future use | None |
| -c | --clonal | If True, mutations detected in more than one molecule will be considered as one mutation for calculating burdens | False |
| -d | --dilute | Set to true when sample and matched normal are from the same starting DNA material | False |
| -gb | --genebed | Gene bed file, if gene coverage needs to be calculated | None |
| -rb | --reestimatebed | Re-estimate bed file, if burden re-estimation is needed | None |


#### Results

**Core Output Files (from call command):**

**{sample_name}_snv.vcf**:
VCF of detected single nucleotide mutations in the sample. The VCF also includes multiple nucleotide mutations (MNVs).

**{sample_name}_indel.vcf**:
VCF of detected short insertion/deletion (indel) mutations in the sample.

**{sample_name}_coverage.bed.gz**:
Coverage file in BED format showing duplex coverage depths across genomic positions. For multi-threaded runs, coverage files from different threads are automatically merged with proper handling of region boundaries and tabix indexing.

**{sample_name}_trinuc_by_duplex_group.txt**:
Trinucleotide context counts grouped by duplex read numbers for mutation burden estimation.

**{sample_name}_duplex_group_stats.txt**:
Statistics for duplex groups including read counts and quality metrics.

**{sample_name}_stats.txt**:
General statistics file containing overall sequencing and analysis metrics.

**Error Profile Files (from call command):**

**{sample_name}.amp.tn.txt**:
Amplification error profile for trinucleotide contexts (SBS errors).

**{sample_name}.amp.id.txt**:
Amplification error profile for indel errors.

**{sample_name}.dmg.tn.txt**:
Damage error profile for trinucleotide contexts (SBS errors).

**{sample_name}.dmg.id.txt**:
Damage error profile for indel errors.

**Burden Estimation Files (from estimate command):**

**{sample_name}_sbs_burden.txt**:
Single base substitution (SBS) burden estimation including both uncorrected and corrected burden estimates with 95% confidence intervals.

**{sample_name}_indel_burden.txt**:
Indel burden estimation with 95% confidence intervals, including both masked and unmasked burden calculations.

**{sample_name}_sbs_96_corrected.txt**:
Corrected single base substitution counts organized by 96 trinucleotide contexts for mutational signature analysis.

**{sample_name}_sbs_burden_by_min_read_group_size.txt**:
Detailed burden estimates by minimum read group size with both uncorrected and corrected burden calculations and confidence intervals.

**{sample_name}_duplex_allele_counts.txt**:
Table showing duplex depths and allele counts for each unique mutation detected.

**Visualization Files (from estimate command):**

**{sample_name}_sbs_96_corrected.png**:
Plot showing the 96 trinucleotide context mutational signature.

**{sample_name}_sbs_burden_by_min_read_group_size.png**:
Plot showing burden estimates across different minimum read group sizes.

**Optional Files (conditional outputs):**

**{sample_name}_snv_flt.vcf** (when dilute mode is used):
Filtered SNV VCF file when dilute analysis is performed.

**{sample_name}_gene_coverage.txt** (when -gb option is used):
Mean duplex coverage for each gene, useful for dNdScv coverage correction.

**{sample_name}_sbs_burden_re_estimate.txt** (when -rb option is used):
Re-estimated burden for specific regions without re-running variant calling.

**{sample_name}_sbs_96_corrected_re_estimate.png** (when -rb option is used):
Mutational signature plot for re-estimated regions.

## Citation
Cheng, Y. et al. Improved Mutation Detection in Duplex Sequencing Data with Sample-Specific Error Profiles. bioRxiv (2025). https://doi.org/10.1101/2025.07.13.664565

## Copyright

Copyright (c) 2024, Yuhe Cheng [Alexandrov Lab] All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

## Contact Information
Yuhe Cheng (yuc211@ucsd.edu)
