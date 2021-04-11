sudo apt update
sudo apt dist-upgrade -y

sudo apt install -y \
  python3 \
  tree \
  bmon \
  htop \
  screen

sudo su -

#####
# Change timezone
#####
dpkg-reconfigure tzdata

#####
# Setup conda
#####
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install to /opt/anaconda
bash ./Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

#####
# Add conda initialisation to /etc/skel/.bashrc
#####
cat >>/etc/skel/.bashrc <<'EOF'
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/opt/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/opt/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/opt/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/opt/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
EOF

#####
# Setup screenrc for new users
#####
cat >>/etc/skel/.screenrc <<'EOF'
startup_message off
caption string "%?%F%{= Bk}%? %C%A %D %d-%m-%Y %{= kB} %t%= %?%F%{= Bk}%:%{= wk}%? %n "
hardstatus alwayslastline
hardstatus string '%{= kG}[ %{G}%H %{g}][%= %{= kw}%?%-Lw%?%{r}(%{W}%n*%f%t%?(%u)%?%{r})%{w}%?%+Lw%?%?%= %{g}][%{R} %d/%m %{W}%c %{g}]'
altscreen on
# Let screen spawn its shell as the user's login shell. This will then also parse /etc/profile /etc/profile.d and ~/.profile
shell -$SHELL
EOF

#####
# Setup conda environment for all users
#####
mkdir -p /data/envs

cat >/data/envs/agrigenomics.yml <<'EOF'
name: agrigenomics
channels:
  - bioconda
  - conda-forge
  - defaults
dependencies:
  - python=3.8
  - fastp=0.20.1
  - fastqc=0.11.9
  - multiqc=1.10.1
  - pigz=2.6
  - sed=4.8
  - mawk=1.3.4
  - samtools=1.12
  - minimap2=2.18
  - spades=3.15.2
  - sra-tools=2.10.9
  - fastx_toolkit=0.0.14
  - filtlong=0.2.0
  - bwa=0.7.17
  - flye=2.8.3
  - mummer4=4.0.0rc1
EOF
#  
#  - mummer=3.23 
#  - r-base=4.0.3
#EOF

# Install environments
conda install mamba -n base -c conda-forge -y
mamba env create -f /data/envs/agrigenomics.yml
# If need to update the environment
#mamba env update --name agrigenomics --file /data/envs/agrigenomics.yml

#####
# Download from CloudStor
#####
curl https://cloudstor.aarnet.edu.au/plus/s/t5hn8HHErEiRlER/download \
  > agrigenomics.tar
tar -C / -xf agrigenomics.tar && rm agrigenomics.tar





#####
# Data Acquisition
#####
conda activate agrigenomics
mkdir -p /data/SARS-CoV-2/{reference,Illumina}

# SARS-CoV-2 Reference
# Download SARS-Cov-2 RefSeq
#####
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512&rettype=fasta&retmode=text" \
  | fasta_formatter \
  | bgzip \
  > /data/SARS-CoV-2/reference/NC_045512.2.fasta.gz
samtools faidx /data/SARS-CoV-2/reference/NC_045512.2.fasta.gz

# Info for Subsampling
#####
GENOME_SIZE=$(cut -f2 /data/SARS-CoV-2/reference/NC_045512.2.fasta.gz.fai)
ILLUMINA_TARGET_COVERAGES=(
  10
  50
)
NANOPORE_TARGET_COVERAGES=(
  10
  50
)
PACBIO_CCS_TARGET_COVERAGES=(
  10
  50
)

# SRA Accessions for Download
#####
ILLUMINA_PE_SRA_ACCESSIONS=(
  SRR11140748
)
NANOPORE_SRA_ACCESSIONS=(
  SRR11140749
)
PACBIO_CCS_ACCESSIONS=(
  SRR13144524
)

# Download SRA data
#####
for ACC in ${ILLUMINA_PE_SRA_ACCESSIONS[@]} ${NANOPORE_SRA_ACCESSIONS[@]} ${PACBIO_CCS_ACCESSIONS[@]}; do
  echo ${ACC}

  fastq-dump --split-files --split-e --gzip ${ACC}
done

#####
# Illumina PE Subsampling
#####
for ACC in ${ILLUMINA_PE_SRA_ACCESSIONS[@]}; do
  echo "Processing: ${ACC}"

  echo -n "  Calculating raw read coverage ... "

  if [ ! -f ${ACC}.tsv ]; then
    pigz -dcp1 ${ACC}_?.fastq.gz \
      | sed -n '2~4p' \
      | awk 'BEGIN{OFS="\t"}{n+=1; tot+=length($1)}END{print tot,n,tot/n}' \
      > ${ACC}.tsv
  fi
  INPUT_READS=$(cut -f2 ${ACC}.tsv | head -n1)
  INPUT_BP=$(cut -f1 ${ACC}.tsv | head -n1)
  INPUT_COVERAGE=$((INPUT_BP/GENOME_SIZE))

  printf "%'.f\n" "${INPUT_COVERAGE}"

  for TARGET_COVERAGE in ${ILLUMINA_TARGET_COVERAGES[@]}; do
    echo "  Target coverage: ${TARGET_COVERAGE}"
    
    if [ -f ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz ]; then
      echo "    SKIPPING: output file(s) already exist"
    else
      
      echo -n "    Subsampling ... "

      paste \
        <(pigz -dcp1 ${ACC}_1.fastq.gz | paste - - - -) \
        <(pigz -dcp1 ${ACC}_2.fastq.gz | paste - - - -) \
        | shuf -n $((INPUT_READS/2*TARGET_COVERAGE/INPUT_COVERAGE)) \
        | tee \
          >(cut -f 1-4 | tr "\t" "\n" | sed '3~4 s/^.\+$/+/' | pigz --best --processes 1 > ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz) \
          | cut -f 5-8 | tr "\t" "\n" | sed '3~4 s/^.\+$/+/' | pigz --best --processes 1 > ${ACC}_2_${TARGET_COVERAGE}x.fastq.gz
      echo "COMPLETE"
    fi
  done
done

#####
# Nanopore Subsampling
#####
for ACC in ${NANOPORE_SRA_ACCESSIONS[@]}; do
  echo "Processing: $ACC"

  echo -n "  Calculating raw read coverage ... "

  if [ ! -f ${ACC}.tsv ]; then
    pigz -dcp1 ${ACC}.fastq.gz \
      | sed -n '2~4p' \
      | awk 'BEGIN{OFS="\t"}{n+=1; tot+=length($1)}END{print tot,n,tot/n}' \
      > ${ACC}.tsv
  fi
  INPUT_READS=$(cut -f2 ${ACC}.tsv | head -n1)
  INPUT_BP=$(cut -f1 ${ACC}.tsv | head -n1)
  INPUT_COVERAGE=$((INPUT_BP/GENOME_SIZE))

  printf "%'.f\n" "${INPUT_COVERAGE}"

  for TARGET_COVERAGE in ${NANOPORE_TARGET_COVERAGES[@]}; do
    echo "  Target coverage: ${TARGET_COVERAGE}"
    
    if [ -f ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz ]; then
      echo "    SKIPPING: output file(s) already exist"
    else
      
      echo -n "    Subsampling ... "

      filtlong \
        --assembly /data/SARS-CoV-2/reference/NC_045512.2.fasta.gz \
        --target_bases $((GENOME_SIZE*TARGET_COVERAGE)) \
        ${ACC}.fastq.gz \
      | gzip \
      > ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz
      echo "COMPLETE"
    fi
  done
done

#####
# PacBio CCS Subsampling
#####
for ACC in ${PACBIO_CCS_ACCESSIONS[@]}; do
  echo "Processing: $ACC"

  echo -n "  Calculating raw read coverage ... "

  if [ ! -f ${ACC}.tsv ]; then
    pigz -dcp1 ${ACC}.fastq.gz \
      | sed -n '2~4p' \
      | awk 'BEGIN{OFS="\t"}{n+=1; tot+=length($1)}END{print tot,n,tot/n}' \
      > ${ACC}.tsv
  fi
  INPUT_READS=$(cut -f2 ${ACC}.tsv | head -n1)
  INPUT_BP=$(cut -f1 ${ACC}.tsv | head -n1)
  INPUT_COVERAGE=$((INPUT_BP/GENOME_SIZE))

  printf "%'.f\n" "${INPUT_COVERAGE}"

  for TARGET_COVERAGE in ${PACBIO_CCS_TARGET_COVERAGES[@]}; do
    echo "  Target coverage: ${TARGET_COVERAGE}"
    
    if [ -f ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz ]; then
      echo "    SKIPPING: output file(s) already exist"
    else
      
      echo -n "    Subsampling ... "

      pigz -dcp1 ${ACC}.fastq.gz \
        | paste - - - - \
        | shuf -n $((INPUT_READS*TARGET_COVERAGE/INPUT_COVERAGE)) \
        | tr "\t" "\n" | sed '3~4 s/^.\+$/+/' | pigz --best --processes 1 \
      > ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz
      
      #filtlong \
      #  --assembly /data/SARS-CoV-2/reference/NC_045512.2.fasta.gz \
      #  --target_bases $((GENOME_SIZE*TARGET_COVERAGE)) \
      #  ${ACC}.fastq.gz \
      #| gzip \
      #> ${ACC}_1_${TARGET_COVERAGE}x.fastq.gz
      echo "COMPLETE"
    fi
  done
done

#####
# Packaging prior to workshop
#####
mkdir -p /data/SARS-CoV-2/{Illumina,Nanopore,PacBio}/
mv SRR11140748* /data/SARS-CoV-2/Illumina/
mv SRR11140749* /data/SARS-CoV-2/Nanopore/
mv SRR13144524* /data/SARS-CoV-2/PacBio/

tar -c -f agrigenomics.tar /data/SARS-CoV-2

rm -f *.fastq.gz *.tsv


