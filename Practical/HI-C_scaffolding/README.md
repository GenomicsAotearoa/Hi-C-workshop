# Post-scaffolding using Hi-C reads – from quality control to visualisation

## Background
This README describes a workflow of how to do Hi-C data quality control, anchoring contigs/scaffolds onto chromosome-level (or post-scaffolding using Hi-C data) and chimeric scaffold correction based on visualization.

## Demonstrators
* Chen Wu (Plant and Food: chen.wu@plantandfood.co.nz)
* Shane Choi (Manaaki Whenua – Landcare Research: scho249@aucklanduni.ac.nz)

## 1. Logging into Jupyter
Click the link  https://jupyter.nesi.org.nz/hub/login and choose the options shown below:
* select project: Genomics Aotearoa Virtual Lab Training Access (nesi02659)
* Select walltime: 2 hs
* Select number of CPUs: 2
* Select memory size: 32 Gb
* Use GPU: No
[**NOTE**: modify them according to your datasets/projects]

## 2. Location of programs and scripts

The programs/scrips that are used in this tutorial include:
* Quality control tools: BWA, samblaster, samtools, matlock, hic_qc
*	Scaffolding tools: ALLHIC and Salsa2
*	Chimeric scaffold correction tool: Juicebox

Some of these programs/scripts have not been installed as modules on NeSI, rather they are installed in the directory under the nesi02659 project. Export this directory to your environment for easy use:

```
export software=/nesi/project/nesi02659/Hi-C/SW/
```

For example, you can use matlock after loading two modules like below:

```
ml load Python/3.8.2-gimkl-2020a
ml load GSL/2.4-GCC-7.4.0
$software/matlock/bin/matlock -h
```
Now, create and go to your working directory.

## 3.	Sample data description

For quick demonstration, we use the published corn leaf aphid Hi-C data to scaffold the genome contigs generated from PacBio long reads (article: https://academic.oup.com/gigascience/article/8/4/giz033/5429686). We will scaffold a PacBio assembly contains 689 contigs (326 Mb in size with N50 of 9Mb). Karyotype number of this species is four.

The Hi-C data and assembly contigs can be found from the directory below and let’s export the data directory to your environment:
```
export data=/scale_wlg_persistent/filesets/project/nesi02659/Hic_practical/
```
This directory also contains intermediate files generated from every one of the steps for you to compare with.

## 4.	Data quality control

### 4.1. Check sequencing quality (similar to normal Illumina libraries)

* Sequencing quality
*	Adapter contamination
*	Percentage of duplicates
```
ml load FastQC/0.11.9
fastqc $data/map/SRR7989292_1.fastq
fastqc $data/map/SRR7989292_2.fastq
```
Please check the fastqc reports and make an adjustment for whether some filtering/trimming tools should be applied. The raw reads will be used directly for the rest of steps in this tutorial.

### 4.2. Check the percentage of informative read pairs

The difference of Hi-C data quality control compared with normal Illumina library is to see whether there are enough read pairs in the Hi-C library that can help with scaffolding contigs. The two reads from a read pair aligned to the same contig are not useful for scaffolding, while those pairs with both reads aligned to different contigs are called informative reads which are efficient to scaffolding. PhaseGenomics recommends the rate of informative read pairs to be above 5%. However, data with lower threshold can be also tried to improve continuity of the assembly. The efficiency of scaffolding also depends on the quality of your contigs (assembly statistics, haplotigs, contaminations, etc...). 

PhaseGenomics has a workflow to calculate the percentage of informative read pairs in the HiC dataset and we are following every one of the steps from this link to generate data reports: https://phasegenomics.github.io/2019/09/19/hic-alignment-and-qc.html. You can also find detailed parameter description.

#### 4.2.1 Map reads to contigs

BWA is used for mapping with a specific parameter called “-5SP” designed for long-range read pairs mapping. We are piping BWA with samblaster (mark duplicates) and samtools (converting to compressed format; “-F 2316” means excluding reads that were not mapped, those with mate unmapped, not primary alignment and those that are supplementary alignment: https://broadinstitute.github.io/picard/explain-flags.html) 

[**NOTE**: if you have module loading problems such as incompatibility between modules, try to run ```ml purge``` to unload all modules first and then continue loading the necessary ones.]

```
ml load BWA/0.7.17-gimkl-2017a
ml load samblaster/0.1.26-GCC-9.2.0
ml load SAMtools/1.3.1-gimkl-2017a
bwa mem -5SP $data/map/contig.fasta $data/map/SRR7989292_1.fastq $data/map/SRR7989292_2.fastq | samblaster | samtools view -S -h -b -F 2316  hic.bam
```
Mapping statistics is shown below:
![Mapping statistics](https://github.com/GenomicsAotearoa/Hi-C-workshop/blob/main/Practical/HI-C_scaffolding/MappingStats.png)

#### 4.2.2 Filter alignment (generally not necessary, but we include here)

If your Hi-C report generated from the next step does not look good, you can try filtering your alignments including discarding those with poor mapping quality, etc. You can find parameter detailed from the PhaseGenomics web page mentioned above.
```
ml load Python/3.8.2-gimkl-2020a
ml load GSL/2.4-GCC-7.4.0
$software/matlock/bin/matlock bamfilt -i hic.bam -o filtered_hic.bam
```
The filtering statistics is shown below from the dataset we use in this tutorial:
![Filtering statistics](https://github.com/GenomicsAotearoa/Hi-C-workshop/blob/main/Practical/HI-C_scaffolding/filteringStats.png)

#### 4.2.3 Generate Hi-C data quality report

```hic_qc.py``` is written in Python 3. So, we need to load the appropriate module and run hic_qc.py on the mapped reads.
```
ml load Python/3.8.2-gimkl-2020a
pip install --upgrade scipy 
python $software/hic_qc/hic_qc.py -b filtered_hic.bam -r -o hic_qc -n 42114971
```
The detailed description of your data quality is reported in the file: ```hic_qc.html```


## 5.	Post-scaffolding using Hi-C data

We use the ```filtered_hic.bam``` generated above as the input for the following scaffolding steps. Salsa2 and ALLHIC have been used for several GA genome projects, however, there are other tools available for trying.

### 5.1. SALSA2 (https://github.com/marbl/SALSA)

For scaffolding human genome using 40x coverage data, SALSA requires 11.43 CPU hours with 21.43 GB for peak memory usage. CPU hours and memory can vary depending on number of CPUs and the quality of the assembly. You may want to make adjustment on your specific projects.

#### 5.1.1. Convert the bam file to bed file as Salsa2 doesn’t take bam format 
```
ml load BEDTools/2.29.2-GCC-9.2.0
bamToBed -i filtered_hic.bam  hic.bed
```
#### 5.1.2. Run Salsa2

Overview of the SALSA2 scaffolding algorithm:
![Salsa2 algorithm](https://github.com/GenomicsAotearoa/Hi-C-workshop/blob/main/Practical/HI-C_scaffolding/Salsa2Workflow.png)

Several parameters should be considered when you run Salsa:

*	```-e``` what restriction enzyme/enzymes you used for your Hi-C library
*	```-i``` how many scaffolding iterations do you want? Multiple iterations can be time-consuming.
*	```-m``` whether to correct mis-joint scaffolds after each iteration.

SALSA is written in Python 2, thus we need to load BEDtools and Python2 to run SALSA:
```
ml load Python/2.7.18-gimkl-2020a
```
Single RE example (used in this tutorial):
```
python $software/SALSA/run_pipeline.py -b hic.bed -a $data/map/contig.fasta_nt -l $data/map/contig.fasta.fai -o test_SALSA -e GATC
```
[**NOTE**: ```contig.fasta.fai``` is an index file, you can generate from ```samtools faidx contig.fasta```]

Multiple REs example:
```
python $software/SALSA/run_pipeline.py -b hic.bed -a $data/map/contig.fasta_nt -l $data/map/contig.fasta.fai -o test_SALSA -e GATC,ATTC
```

### 5.2 ALLHIC (https://github.com/tangerzhang/ALLHiC)

The unique feature of ALLHIC is that it can phase and then scaffold autopolyploid contigs into different chromosome-level haplotypes if you can provide allelic information/annotation (https://www.nature.com/articles/s41588-018-0237-2). It also has a workflow for scaffolding diploid assemblies, which we’ll use in this tutorial.

Several parameters should be considered when you run ALLHIC:
* Provide a group of cluster (karyotype) using ```-k```
*	Restriction enzyme using ```-e``` & ```--RE``` [NOTE: ALLHIC does not provide multi-RE setting at the moment]

**Step1**: ```partition``` separates all the contigs into clusters using average links:
```
$software/ALLHiC/bin/ALLHiC_partition -b filtered_hic.bam -r $data/map/contig.fasta -e GATC -k 4
```
**Step2**: ```extract``` does preprocessing: 1) extract inter-contig links; 2) extract intra-contig links and build a distribution; 3) count up the restriction sites to be used in normalization; 4) bundles the inter-contig links into pairs of contigs.
```
$software/ALLHiC/bin/allhic extract filtered_hic.bam $data/map/contig.fasta --RE GATC
```
**Step3**: ```optimize``` given a set of Hi-C contacts between contigs, as specified in the clm file, reconstruct the highest scoring ordering and orientations for these contigs. Using a bash loop here to run every one of the karyotype group:
```
for i in {1..4}; \
do "allhic optimize filtered_hic.counts_GATC.4g$i.txt \
filtered_hic.clm"; done
```
**Step4**: ```build``` generates fasta file of the assembly
```
$software/ALLHiC/bin/ALLHiC_build $data/map/contig.fasta
```
**Tips**: the Salsa2 assembly can be also as input to the ALLHIC program for post-scaffolding using Hi-C data.

## 6.	Visualization and mis-joint scaffold correction using JuiceBox

### 6.1	generate input files for JuiceBox

To visualize the contact map for your post-scaffolded assembly using Hi-C data, we need to map the Hi-C reads to the new assembly (i.e. Salsa2/ALLHIC output). Mapping is the same with quality control step described above.

**Step1**: generate genome assembly information files (taking Salsa2 output as an example)
```
python $software/juicebox_scripts/juicebox_scripts/makeAgpFromFasta.py $data/SALSA/scaffolds_FINAL.fasta out.agp
python $software/juicebox_scripts/juicebox_scripts/agp2assembly.py out.agp out.assembly
```
**Step2**: generate Hi-C contact file ```.hic``` [NOTE: the ```filtered_hic_new.bam``` file is the one you need to generate from mapping Hi-C data to your Salsa2 assembly using the quality control mapping step described above.]
```
$software/matlock/bin/matlock" bam2 juicer $/data/visualise/filtered_hic_new.bam out.links.txt
sort -k2,2 -k6,6 out.links.txt  out.sorted.links.txt ### This may take a while.
$software/3d-dna/visualize/run-assembly-visualizer.sh -p false out.assembly out.sorted.links.txt
```

### 6.2 JuiceBox (https://github.com/aidenlab/Juicebox) 

The online version JuiceBox has limited functions, it is recommended to download it on your desktop for running or running it on server.
```
java -jar $software/juicer/Juicebox_1.11.08.exe
```
It is a GUI tool, so you will need to select the .hic and .assembly file manually by clicking on FileOpen for .hic file and AssemblyImport Map assembly for .assembly file. If you are running on Jupyer, juicebox does not work. Please download the software to your desktop through this link: https://github.com/aidenlab/Juicebox/wiki/Download

When you open all the files, this is the contact map it will visualize:

![contact map](https://github.com/GenomicsAotearoa/Hi-C-workshop/blob/main/Practical/HI-C_scaffolding/Heatmap_JuiceBox.png)

There is a youtube video tutorial about how to use this tool specifically:
https://www.youtube.com/watch?v=Nj7RhQZHM18
