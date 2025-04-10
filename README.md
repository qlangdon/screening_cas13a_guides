# *In Silico* Guide Design and Screening for Viral Diagnostics with CRISPER-Cas13a

This pipeline was developed to design Cas13a guides for amplification free detection of RNA virses, as was done in [Fozouni, et al. Cell 2021](https://doi.org/10.1016/j.cell.2020.12.001). It is based on and uses elements from the following sources:

- Cas13a guide design https://github.com/lareaulab/cas13a_guide_design/tree/main
- In Silico screening https://github.com/czbiohub-sf/sc2-guide-InSilicoSCR
- Influenza detection https://github.com/duopeng/Cas13a_guides_Influenza

This workflow can be adapted for Cas12 guide design as well. For a working model of that see here:

- https://github.com/fletchlab-git/Cas-guide-pipeline

## Setup

#### clone repository
```shell
git clone https://github.com/qlangdon/screening_cas13a_guides
```

#### Other needed software

ViennaRNA RNAfold (https://www.tbi.univie.ac.at/RNA/RNAfold.1.html)

bowtie (https://bowtie-bio.sourceforge.net/index.shtml)

##### Most recent human transcriptome to check for cross-reactivity
```
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_rna.fna.gz

bowtie-build GRCh38_latest_rna.fna GRCh38_latest_rna
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

1. Check for secondary structure
- This step comes from https://github.com/lareaulab/cas13a_guide_design 

```
Rscript /scripts/score_RNAfold_crRNAs.R /
--enzyme Cas13a /
--rnafold /software/ViennaRNA-2.6.4/src/bin/RNAfold /
--out .
```

2. Filter to keep only those that will form a hairpin with the Cas stem sequence
- This step will use the output filtered from the total hits, so be sure to include the threshold used for that.
```
Rscript /scripts/filterGuidesFoldScore.R virus_name_al 0.5
```

3. Check for cross reactivity to human transcriptome
- This script comes from https://github.com/duopeng/Cas13a_guides_Influenza which was adapted from https://github.com/lareaulab/cas13a_guide_design
```
Rscript /scripts/align_bowtie.R --genome GRCh38_latest_rna --enzyme Cas13a --out .
```

4. Filter to remove any that map to the human transcriptome and further winnow down the top hits by within group or species inclusivity
- You can pick any cutoff theshold you would like, the script will report how many remain when given a percent cutoff
```
Rscript /scripts/filterGuidesHumanAlign_cutoffProp.R virus_name_all 0.9
```

## Step 4: User directed filtering

At this stage you have a working table  of potential crRNA guides (e.g. `virus_name_all_guides_windows_wHumanAlignZero_kept0.9.tsv`) that might still be quite long. This table will include the number of strains hit as well as proportion of strains hit for each group and in total. From this you can filter for guides that are specific for one group or will work globaly. The crRNA guides in this table pass the basic criteria for Cas13a or Cas12 guides and could be ordered and tested. 

Below is an additional, but more time and computational power intensive, step that could be taken to check more broadly for cross reactivity with other virses, pathogens, bacteria, or even more broadly. 

## Step 5: Check for cross-reactivity with a metagenomic classifier

1. This step uses KrakenUniq a metagenomics taxonomic classifier to assay if a crRNA guide target matches any other known sequence at any taxonomic level. The full workflow you can follow is found here:

https://github.com/czbiohub-sf/sc2-guide-InSilicoSCR

- Additonally you can add Kraken pre-built indexes to expand the range of databases tested. The "Standard" index is a good option. 
https://benlangmead.github.io/aws-indexes/k2

2. To classify if the guides are on-target or off-target you will need a list of the taxIDs that match your organism of interest. For this you will need the Taxonomy ID of  organism of interest to pull all daughter taxIDs from the Kraken taxDBs.
 
Example using Dengue, *Orthoflavivirus denguei* Taxonomy ID: [3052464](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=3052464&lvl=3&lin=f&keep=1&srchmode=1&unlock)
```
Rscript /scripts/getTaxIDsFromKrakenGeneric.R Dengue 3052464
```

3. This can then be combined with the classified output and the full table. Below is only an example script of a way to do this.
```
Rscript /scripts/parseSummariseKrakenOutputByClade.R /
virus_name /
virus_name_all_guides_windows_wHumanAlignZero_kept0.9.tsv /
targets.fa
```















