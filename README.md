# salmdetectpipeline

This repository holds the scripts to run the analyses used in [paper link].

This work was performed on Graham HPC provided by the Digital Research Alliance of Canada (formerly Compute Canada) and is formatted to run through SLURM job submissions. 

## Databases
Details on the creation or download of the Kraken 2 libraries are given in the publication's supplementary information.

## Creating or reconstituting libraries
Libraries are created from "master" files containing simulated shotgun sequencing reads. Details on the genomes or publicly available datasets used for these master files are available in the publication's supplementary information. They are not provided here due to file size.

### Create new libraries

To create libraries with background sequences from the "enterobac" mock microbiome and a particular number of Salmonella reads spiked in:
Edit information in 
`configs/mockcomm_step1_setup.yaml`
Give the path and names of the background community fastq files and the path and names of the Salmonella simulated reads fastq files.
Define the total library size and the number of Salmonella reads to include in that library.

Then, submit a job to create libraries at that spike-in level:
`SLURM_mockcomm_step1_setup.sh`,
which itself runs a snakefile:
`snakefile_ccmockcom_step1setup.py`

This must be repeated for each desired spike-in level.

### Save read names for future reconstitution
Fastq files of the generated libraries are large, and scratch directories on Compute Canada servers are routinely wiped, so external storage is necessary. To reduce storage space, lists of read names can be saved and the full fastq files "reconstitited" from this data using the below scripts and the master files.

Use the `printseqnames.py` script to get read names from fastq files. It will create a new text file of read names (1 per line) for each "R1_fq.gz" file in the directory where the script is called.

### Reconstitute library from read names
Using these read names and the master files described above, fastq files containing the library sequences can be re-made ("reconstituted"). If library fastqs are needed, edit and run:
`submit_slurm_recon.py`
(Be sure to edit the reps and spikes list in the python file to reflect the libraries you want).
This will submit one job per library via 
`SLURM_reconstitutelibs.sh`
which is a slurm submission script using bbmap's filterbyname.sh. There is no snakemake file for this step.

## Analysis steps

### Kr2-SSR pipeline

Once libraries are made, they are put through the analysis pipeline.

For most efficient use of HPC resources, the whole snakemake workflow including trimming, Kraken 2 classification, and comparison of putative Salmonella reads is given in a single SLURM submission per round + set of libraries. Use different SLURM submission files for the separate analysis rounds.

#### Round 1
 Round 1 of analysis aimed to find the best parameters and databases for Salmonella detection, as well as test the efficacy of SSR-confirmation for weeding out false positives. To do this, all libraries went through these combinations:
- Kraken2 confidence settings: 0, 0.25, 0.5, 0.75, 1
- Kraken2 databases: kr2bac, kr2std, kr2plrename

Round 1 relies on config file  `configs/mockcomm_Round1.yaml`
And snakemake file `snakefile_mockcomm_R1_step2analyze.py`

Submit a job for Round 1 using 
`SLURM_mockcomm_R1_analyze.sh` 
Be sure to adapt the time and memory allocations based on the number of libraries.


#### Round 2
 Round 2 of analysis aimed to find limits of detection using the best parameters found in Round 1:
- Kraken2 confidence settings: 0.25
- Kraken2 databases: kr2bac

Round 1 relies on config file  `configs/mockcomm_Round2.yaml` 
And snakemake file `snakefile_mockcomm_R2_step2analyze.py`
This workflow loads only a single database into temporary storage, whereas the workflow in round 1 loads 3 databases.

Submit a job for Round 2 using 
`SLURM_mockcomm_R2_analyze.sh` 
Be sure to adapt the time and memory allocations based on the number of libraries.

## Analysis results summary
After completing analysis on all the desired libraries, an additional SLURM job can be submitted to produce a summary table.

This analysis runs a snakemake workflow `snakefile_mockcomm_step3summarize.py`,
which relies on scripts
`scripts/count_checkhits.py`,
`scripts/summarize_hits.py`,
`scripts/summarize_reads.py`,
`scripts/summary_table.py`,
as well as the same config files used for analysis workflows.

To submit a summary table job for Round 1, use:
`SLURM_mockcomm_R1_step3summarize.sh`

To submit a summary table job for Round 2, use:
`SLURM_mockcomm_R2_step3summarize.sh`


### Metaphlan4 pipeline
First, download the necessary databases (on Compute Canada infrastructure, these very large files must be stored on scratch, which is subject to periodic wiping).
- There's a Compute-Canada-specific set of download instructions here, but it's incomplete/out of date: https://docs.alliancecan.ca/wiki/MetaPhlAn
- Downloading the few files from here http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/, as in instructions above, is insufficient. You need to also download the md5 file AND the bowtie-indexed databases here: http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/

```
parallel wget ::: http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.tar http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103_marker_info.txt.bz2 http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103_species.txt.bz 
```
```
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.md5
```
```
wget http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vJan21_CHOCOPhlAnSGB_202103_bt2.tar
```
```
http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/bowtie2_indexes/mpa_vJan21_CHOCOPhlAnSGB_202103_bt2.md5
```

Then run Metaphlan4 via
`SLURM_snakemake_metaphlan.sh`, which submits a job running
`snakefile_MCmetaphlan_step2.py`, which relies on configs/mphl4_configfile.yaml
Be sure to update the config before running! File names must be given in a list in the config file rather than being dynamically picked up from directory contents.


 

