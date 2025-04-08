# *In Silico* Guide Design and Screening for Viral Diagnostics with CRISPER-Cas13a

This pipeline was developed to design Cas13a guides for amplification free detection of RNA virses, as was done in [Fozouni, et al. Cell 2021]([10.1016/j.cell.2020.12.001](https://doi.org/10.1016/j.cell.2020.12.001)).  It can be adapted for Cas12 guide design as well. It is based on and uses elements from the following sources:

- Cas13a guide design (https://github.com/lareaulab/cas13a_guide_design/tree/main)
- In Silico screening (https://github.com/czbiohub-sf/sc2-guide-InSilicoSCR)
- Influenza detection (https://github.com/duopeng/Cas13a_guides_Influenza)


## Setup

#### clone repository
```shell
git clone https://github.com/qlangdon/screening_cas13a_guides
```

#### create conda environment
```shell
conda env create -f env.yml
```

## Step 1: Download genomes and metadata

[Nextstrain](https://nextstrain.org/) – pulls from NCBI and annotates clades
- Good for RSV and other common pathogens

[NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/find-data/virus) – useful if large dataset where subtype is commonly reported
- Good for Influenza and pathogens where you want to filter for genome completeness

## Step 2: Break genomes into every 20bp window and count frequency

1. Split initial fasta by subtype (or grouping of your choice) based on columns in the metadata file
2. Then split each accession into 20bp sliding windows and summarize strandedness and GC content

Example of **non-segmented** genome:
```shell
python /scripts/split_fasta_bySubtype_breakIntoWindows.py \
--metadata genome_metadata.tsv \
--fasta sequences.fasta \
--idColumn accession \
--subtypeColumn clade
```

Example of **segmented** genome where you want the count frequency to reflect both number of accessions and number of unique strains:
```shell
python /scripts/split_fasta_bySubtype_breakIntoWindows.py \
--metadata metadata.tsv \
--fasta sequences.fasta \
--idColumn Accession \
--subtypeColumn Genotype \
--strainColumn Isolate \
--segmented
```

3. Count every instance of that 20mer is found
- Filter for GC content (>25% and <75%)
- Count per group total and precent of group strains hit
- Rename alphanumerically to keep naming consistent

4. Sum over all subtypes (or continue by group if desired)
```shell
python /scripts/count_combine_guides_generic.py /
--metadataPrefix genome_metadata /
--pathPrefix windowed_genomes/ 
--output virus_name
```

5. Filter for top hit
- The previous step creates a list of every possible 20mer, even if it only exists in one genome
- Select just the targets that hit many strains, for example at least 50% of all strains or any group
```shell
Rscript /scripts/reformatFilterGuideHits.R virus_name_all 0.5
```

## Step 3: Check and filter for Cas (13a or 12) specifics 








