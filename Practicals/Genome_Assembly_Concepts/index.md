# Genome Assembly Concepts
{:.no_toc}

* TOC
{:toc}

# Setup for Today

You will access a remote computer (a virtual machine or VM) using the command line.
The VM has been configured with all the tools and data required for today's hands-on practical.

## Accessing Your Virtual Machine

If you are using a **Mac or Linux** computer, the "Termial" is already available for you to use.
If you are using a **Windows** PC, please download either Putty ([https://www.putty.org]) or Mobaxterm ([https://mobaxterm.mobatek.net]) to use as your terminal.

We will circulate the login credentials for you to use for the hands-on.
This will include the **IP address** of the machine you are going to access as well as a **username** and **password** combination.

## Data

Attempting to do genome assembly with a real agriculturally relevant species isn't possible with the time we have.
Instead, we are going to look at using some very topical data, SARS-CoV-2 data, the virus causing COVID-19.

The advantage of this genome is that it's small (29,903 bp) and has a lot of different types of data available for it, including Illumina, Oxford Nanopore and PacBio.
We are using the following public data sets:

| Description         | Accession/URL                                                                                | Coverage |
|:--------------------|:--------------------------------------------------------------------------------------------:|---------:|
| RefSeq              | [NC_045512.2](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2?report=fasta)                 |          |
| **Illumina PE**     | [SRR11140748](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR11140748)                    |    6,353 |
| **Oxford Nanopore** | [SRR11140749](https://sra-download.ncbi.nlm.nih.gov/traces/sra53/SRR/010879/SRR11140749)     |   10,170 |
| **PacBio CCS**      | [SRR13144524](https://sra-download.ncbi.nlm.nih.gov/traces/sra19/SRR/012836/SRR13144524)     |  110,696 |

The amount of coverage provided by these data sets is *waaay* too much for our needs.
Therefore, we have provided subsamples of each for working on during this hands-on session.

## Software Environment

We have installed all the software you will need for this hands-on in a conda software environment called `agrigenomics`.
If you're interested to know about what tools and versions are installed, take a peek at the `/data/envs/agrigenomics.yml` file once you log into the VM.

To load the environment and make those tools available to you on the command line, you'll need to "activate" the conda environment:

```bash
conda activate agrigenomics
```

## Setting up Your Working Directory

Once you've logged into the VM with the credentials we've provided, you now need to setup a working directory with a copy of all the data you're going to work with.

```bash
# Setup project working directory
mkdir --parents ~/SARS-CoV-2/

# Get a copy of the data
cp --recursive /data/SARS-CoV-2/* ~/SARS-CoV-2/

# See what you have in you working directory
tree ~/SARS-CoV-2/
```

# Resequencing

## Illumina Reads

### 3' Quality Drop-Off

During the bridge-amplification stage, millions of clusters are created on the flowcell.
Each cluster comprises of 1,000-2,000 identical copies of the same template.
During the process of sequencing, the polymerases attached each of the thousands of copies in a cluster "advance" one base at a time.
At the start (5' end) all the polymerases are in perfect sync; generating a bright, clean, consistant light signal for detecting which base was incorporated.
However, as time progresses, some polymerases fall behind while some race infront.
The polymerases gradually get further and further out of phase with each other.
This leads to dimmer and less clear signals for detection, and thus lower quality base-calls.

The phasing problem seen with Illumina sequencing has improved dramatically over the years, resulting in read lengths increasing from 35 bp to 250 bp.

### Raw Read Quality Control

FastQC is a common tool used for assessing the quality of Illumina data.

```bash
cd ~/SARS-CoV-2/

fastqc ~/SARS-CoV-2/Illumina/SRR11140748_?_*x.fastq.gz
```

#### Questions

 * Take a look at the HTML files generated `~/SARS-CoV-2/Illumina/`
 * What sort of quality based trimming is required?
 * Are any adapter sequences detected?
 * What are the read lengths?

### Raw Read Quality and Adapter Trimming

We'll use `fastp` to quality and adapter trim our reads:

```bash
# Make the output directories ahead of time
mkdir --parents ~/SARS-CoV-2/qc_reads/Illumina

time fastp \
  --thread 2 \
  -i ~/SARS-CoV-2/Illumina/SRR11140748_1_50x.fastq.gz -I ~/SARS-CoV-2/Illumina/SRR11140748_2_50x.fastq.gz \
  -o ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_1_50x.fastq.gz --unpaired1 ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_1_50x.orphans.fastq.gz \
  -O ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_2_50x.fastq.gz --unpaired2 ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_2_50x.orphans.fastq.gz \
  --cut_right --cut_window_size 4 --cut_mean_quality 20 \
  --length_required 75 \
  --json ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_50x.json \
  --html ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_50x.html
```

### Post-Trimming FastQC

To assess how the trimming has performed, run FastQC across the paired reads output of fastp.

```bash
fastqc --threads 2 \
  ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_?_50x.fastq.gz
```

#### Questions

 - How does the raw and post-trimming QC reports differ?

### Read Alignment

There are many tools for performing Illumina PE read alignment (aka mapping) against a genome.
Today we will use a tool called BWA.
Conceptually, a read aligner aims to find the most likely place in the genome from which the read was generated.
Think of this as looking up a persons name in a phonebook.

#### Questions

Would it be easier to search for a name if:

 - *The names in the phonebook were sorted in alphabetical order or if they were in random order?*
 - *There was a tab index?*

![Tab index](images/phonebook-tabs.jpg)

#### Indexing the Reference Genome

Having the information stored in a certain way and having an index makes looking up names in the phonebook much faster.
In the context of aligning reads to a reference genome, having an "index" of the reference genome speeds up the process of finding and aligning reads to the genome.
The index of a genome is often tool-specific so we often need to create a different index for each tool we plan to use.

Like sorting and indexing a phonebook, this process is time consuming with larger genomes.
However, since it's only done once, and it provides significant speed-ups in read alignment, it is time well spent!

Index the reference genome:

```bash
# Index the reference genome for use with BWA
bwa index ~/SARS-CoV-2/reference/NC_045512.2.fasta.gz
```

##### Questions

 - *What 5 files were created by the `bwa index` command? Hint: look at the contents of `~/SARS-CoV-2/reference/`*

#### Align the Illumina Reads

Align the Illumina PE reads to the reference genome:

```bash
# Create an output directory ahead of time
mkdir --parents ~/SARS-CoV-2/alignments/

# Align reads
time bwa mem \
  -t 2 \
  ~/SARS-CoV-2/reference/NC_045512.2.fasta.gz \
  ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_1_50x.fastq.gz \
  ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_2_50x.fastq.gz \
| samtools sort \
  --threads 2 \
  -o ~/SARS-CoV-2/alignments/Illumina_PE_50x.bam
```

#### BAM File Visualisation

We are going to load the SARS-CoV-2 genome into a genome browser called IGV and then load the BAM file.
First, we need to create an index for the reference genome and an index file of the BAM file.
This ensures IGV can quickly load the reference sequence and corresponding read alignments for any coordinates we choose to navigate to.

```bash
# IGV-web requires an uncompressed reference sequence 
# and corresponding index file
pigz -dcp2 \
  < ~/SARS-CoV-2/reference/NC_045512.2.fasta.gz \
  > ~/SARS-CoV-2/reference/NC_045512.2.fasta

# Generate index of the FASTA-formatted genome file
samtools faidx ~/SARS-CoV-2/reference/NC_045512.2.fasta

# Index the BAM file
samtools index ~/SARS-CoV-2/alignments/Illumina_PE_50x.bam
```

Download the following files to your local computer:

 * `~/SARS-CoV-2/reference/NC_045512.2.fasta`
 * `~/SARS-CoV-2/reference/NC_045512.2.fasta.fai`
 * `~/SARS-CoV-2/alignments/Illumina_PE_50x.bam`
 * `~/SARS-CoV-2/alignments/Illumina_PE_50x.bam.bai`

Visit, [IGV-web](https://igv.org/app/) and load the genome from a `Local File ...` by selecting both the `NC_045512.2.fasta` and `NC_045512.2.fasta.fai` files.
Once the reference genome is loaded, load a "Track" from a `Local File ...` by selecting both the `Illumina_PE_50x.bam` and `Illumina_PE_50x.bam.bai` files.


## Oxford Nanopore Reads

Long reads have the potential for spanning "repeat" regions, untangling the mess that short reads make in these areas, provide evidence of structural variation etc.
Lets align the Nanopore reads to the reference genome.

```bash
# Index the reference genome for minimap2
minimap2 \
  -d ~/SARS-CoV-2/reference/NC_045512.2.fasta.gz.mmi \
  ~/SARS-CoV-2/reference/NC_045512.2.fasta.gz

# Align the Nanopore reads
time minimap2 \
  -ax map-ont \
  -t 2 \
  ~/SARS-CoV-2/reference/NC_045512.2.fasta.gz.mmi \
  ~/SARS-CoV-2/Nanopore/SRR11140749_1_50x.fastq.gz \
| samtools sort \
  --threads 2 \
  -o ~/SARS-CoV-2/alignments/Nanopore_50x.bam

# Index the BAM file
samtools index ~/SARS-CoV-2/alignments/Nanopore_50x.bam
```

Load the nanopore read alignment file into your IGV-web.

### Questions

 - *What do you make of the accuracy of the Nanopore reads compared to the Illumina reads?*
 - *Despite their poorer accuracy, do you still think Nanopore reads are useful?*
 - *What do you make of positions: `NC_045512.2:17373` and `NC_045512.2:20299`?*

## PacBio CCS Reads

PacBio CCS reads are longer than Illumina reads but typically shorter than Nanopore reads.
However, CCS reads have base accuracies far better than Nanopore and are approaching that of Illumina data.

```bash
# We already indexed the genome before, no need to do it again
#minimap2 \
#  -d ~/SARS-CoV-2/reference/NC_045512.2.fasta.gz.mmi \
#  ~/SARS-CoV-2/reference/NC_045512.2.fasta.gz

# Align the PacBio CCS reads
time minimap2 \
  -ax asm20 \
  -t 2 \
  ~/SARS-CoV-2/reference/NC_045512.2.fasta.gz.mmi \
  ~/SARS-CoV-2/PacBio/SRR13144524_1_50x.fastq.gz \
| samtools sort \
  --threads 2 \
  -o ~/SARS-CoV-2/alignments/PacBio_CCS_50x.bam

# Index the BAM file
samtools index ~/SARS-CoV-2/alignments/PacBio_CCS_50x.bam
```

# De Novo Genome Assembly

Wed 14th April 2021 marks the 18th year since the human genome was declared complete!
Sequencing technologies, computing power and genome assembly software have come an awful long way since then!
Up until about 5 years ago, Illumina sequencing was at the heart of most genome assembly projects, mainly because of its low cost and high throughput.
However, it has its problems.
Mainly due to the short reads it generates.
Repeats cannot be resolved/assembled if read lengths are shorter than the repeat length.
As a result Illumina-based assemblies tend to be highly fragmented, limiting their usefulness.

However, long-read technologies such as Oxford Nanopore and PacBio provide sequence reads orders of magnitude longer than Illumina.
Therefore, they are capable of spanning and resolving repeats in the genome.
However, they are more expensive.
The long reads are the result of observing a single read/polymerase at work, and as such have a high base-level error rate (15-20%).
PacBio's template has a circular topology so they polymerase can continue sequencing around and around the template, making multiple passes of the same insert and generating highly accurate CCS reads.

The utility of long reads has seen them used for genome assembly projects.
However, their high cost usually means they are paired with Illumina data in order to keep costs down.

## Hybrid Genome Assembly

We will use the tool SPAdes (spades.py) to perform a de novo assembly from: 1) Illumina plus Nanopore reads and 2) Illumina plus PabCio CCS reads. 

```bash
# Create the output directory ahead of time
mkdir --parents ~/SARS-CoV-2/de_novo_assembly/

# Illumina-only assembly
time spades.py \
  --threads 2 \
  -o ~/SARS-CoV-2/de_novo_assembly/illumina \
  -1 ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_1_50x.fastq.gz \
  -2 ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_2_50x.fastq.gz \
| tee ~/SARS-CoV-2/de_novo_assembly/illumina.log

# Illumina plus Nanopore assembly
time spades.py \
  --threads 2 \
  -o ~/SARS-CoV-2/de_novo_assembly/illumina_nanopore \
  -1 ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_1_50x.fastq.gz \
  -2 ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_2_50x.fastq.gz \
  --nanopore ~/SARS-CoV-2/Nanopore/SRR11140749_1_50x.fastq.gz \
| tee ~/SARS-CoV-2/de_novo_assembly/illumina_nanopore.log

# Illumina plus PacBio CCS assembly
time spades.py \
  --threads 2 \
  -o ~/SARS-CoV-2/de_novo_assembly/illumina_pacbio_ccs \
  -1 ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_1_50x.fastq.gz \
  -2 ~/SARS-CoV-2/qc_reads/Illumina/SRR11140748_2_50x.fastq.gz \
  -s ~/SARS-CoV-2/PacBio/SRR13144524_1_50x.fastq.gz \
| tee ~/SARS-CoV-2/de_novo_assembly/illumina_pacbio_ccs.log
```

## Comparing Assemblies to a Reference Genome

Compare the assembly to the SARS-CoV-2 RefSeq assembly using MUMmer’s nucmer tool:

```bash
# Create output directory ahead of time
mkdir --parents ~/SARS-CoV-2/comparisons

# Decompress the reference genome for MUMmer
pigz -dcp2 \
  < ~/SARS-CoV-2/reference/NC_045512.2.fasta.gz \
  > ~/SARS-CoV-2/reference/NC_045512.2.fasta

# Run nucmer from MUMmer package
nucmer \
  --maxmatch \
  --minmatch 100 \
  --mincluster 500 \
  --prefix ~/SARS-CoV-2/comparisons/illumina \
  ~/SARS-CoV-2/reference/NC_045512.2.fasta \
  ~/SARS-CoV-2/de_novo_assembly/illumina/contigs.fasta

nucmer \
  --maxmatch \
  --minmatch 100 \
  --mincluster 500 \
  --prefix ~/SARS-CoV-2/comparisons/illumina_nanopore \
  ~/SARS-CoV-2/reference/NC_045512.2.fasta \
  ~/SARS-CoV-2/de_novo_assembly/illumina_nanopore/contigs.fasta

nucmer \
  --maxmatch \
  --minmatch 100 \
  --mincluster 500 \
  --prefix ~/SARS-CoV-2/comparisons/illumina_pacbio_ccs \
  ~/SARS-CoV-2/reference/NC_045512.2.fasta \
  ~/SARS-CoV-2/de_novo_assembly/illumina_pacbio_ccs/contigs.fasta
```

### Visualisation of Delta Files

The online tool [Assemblytics](http://assemblytics.com) can be used to generate useful plots of the information contained within a .delta file.

