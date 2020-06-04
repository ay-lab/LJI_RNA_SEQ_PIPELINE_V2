**LJI RNASEQ QC PIPELINE V2**
======
* Niu Du (ndu [at] lji.org)
* AY and VIJAY labs
* La Jolla Institute for Immunology (LJI)
* La Jolla, CA USA
* Current version: 2.0 (05/27/2020)
------
#  Summary

The pipeline was developed for RNA-Seq read mapping and QC as part of the interactive work flow between sequence and bioinformatics team. The version 2 pipeline was implemented using Snakemake in Sun Grid Engine (SGE) environment, and it can be can be implemented in any cluster and cloud environments without significant modifications. Compare to version 1 pipeline, we added features that enables recording history of each QC iteration and integrated the 3 steps approach to 1 step operation. The workflow and pipeline structure are as follows:

**Workflow**

<img src='./img/workflow.png' width = 800>

**Running pipeline**

Prepare necessary files -> qsub pbs_submit.sh -> check report

**Pipeline structure**

<img src='./img/snakemake.png' width = 800>

# Pipeline setup



#### Step 1:
Make sure the following tools have been installed on the server. Please note that different versions of these tools may result in different results.

[fastp](https://github.com/OpenGene/fastp)

[STAR-seq mapping](https://github.com/alexdobin/STAR)

[Qualimap](http://qualimap.bioinfo.cipf.es/)

[Samtools](http://www.htslib.org/)

[Deeptools](https://deeptools.readthedocs.io/en/develop/index.html)



#### Step 2:
Make a clone of this git repository, and use the commend below to install [Snakemake](https://snakemake.readthedocs.io/en/v3.9.1/) and all dependencies for the pipeline. We recommend the  version 3.9.1 for stable execution.

```git clone https://github.com/ndu-UCSD/LJI_RNA_SEQ_PIPELINE_V2.git```

if your python version <= 3.6

```pip install snakemake==3.9.1```

if your python version >=3.7

```pip install snakemake-py37```

Then install all dependencies

```pip install -r requirements.txt```



# Data prepration
#### Samples and running sheet:

For running the pipeline you will need 1. all fastq files with proper names, 2. a sample run table with information on which samples to be used for the analysis, and 3. (optional) a meta table with meta information for all of the samples listed in the sample run table.  

Example of the fastq file structure for storage. The subfolder names here (original and reseq_1) are seq runs executed by the sequencing team.
```
1.Fastq_input
├── original
│   ├── S1_R1.fastq.gz
│   ├── S1_R2.fastq.gz
│   ├── S2_R1.fastq.gz
│   └── S2_R2.fastq.gz
└── reseq_1
    ├── S2_R1.fastq.gz
    └── S2_R2.fastq.gz
```
Example of the sample run csv file. The header designated the pipeline runs executed by the bioinfo team and inside the cells are fastq files to be used. In case of two items in one cell, e.g. S2 Run_3, fastq files from the "original" folder and the "reseq_1" folder will be merged and used in the pipeline.

* The sample name here must be consistent with names in the fastq files.


|sample_ID |Run_1 |Run_2 |Run_3 |
|----------|----------------|----------------|----------------|
|S1        |original        |original        |original        |   
|S2        |original        |reseq_1        |original,reseq_1         |

#### pipeline and config file:

Copy and paste the .sh and .json files from template folder to the working directory, and make proper changes to the configuration file.

* "app" - locations where the above listed apps were installed, and where the pipeline was downloaded to.

* "config" - directories and parameters for the next run; the example was configured for GRCH37.P13 as reference genome.    

* "QC_threshold" - The default QC parameters shown here are optimized for bulk RNA, Smart-Seq2 runs. Please consulate with sequence team before making changes. More details are available at [Change QC rules](#Change QC rules).


Example for LJI implementation.
```
{
    "app":{
        "fastp":"/mnt/BioHome/ndu/anaconda3/bin/fastp",
        "STAR":"/mnt/BioHome/ndu/anaconda3/bin/STAR",
        "samtools":"samtools",
        "bamCoverage":"/mnt/BioHome/ndu/anaconda3/bin/bamCoverage",
        "qualimap":"/mnt/BioHome/ndu/anaconda3/envs/snake/bin/qualimap",
        "pipeline":"/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/RNA_SEQ_Script/LJI_RNA_SEQ_PIPELINE_V2/pipeline"
    },

     "config": {
        "run_file":"/mnt/BioAdHoc/Groups/vd-vijay/VIJAY_LAB_RNA_SEQ/Project_SeAs/SeAs_CD4_Naive_Resting/sample_run.csv",
        "dirin": "/mnt/BioAdHoc/Groups/vd-vijay/VIJAY_LAB_RNA_SEQ/Project_SeAs/SeAs_CD4_Naive_Resting",
        "metadata_dir":"/mnt/BioAdHoc/Groups/vd-vijay/VIJAY_LAB_RNA_SEQ/Project_SeAs/SeAs_CD4_Naive_Resting/metadata_all.csv",
        "run_ID":"Run_1",
        "genome_version":"hg19",
        "bed_dir": "/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/hg19_GencodeCompV19.bed",
        "ref_dir": "/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/GRCH37.P13",
        "gtf_dir": "/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/GRCH37.P13/gencode.v19.annotation.gtf",
        "annotation_file": "/mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference/GRCh37_annotation.csv"
     },

     "QC_threshold": {
        "final_STAR_counts": 3000000.0,
        "uniquely_mapped_reads_perc": 60,
        "too_short_reads_perc": 50,
        "exonic_perc": 50,
        "STAR_counts_perc": 80,
        "Total_genes":[5000,9000],
        "t_rRNA_counts_perc": 10,
        "bias_5to3_prim":[0.95,1.25],
        "insert_median":[150,400],
        "gene_recovery_perc":0.9,
        "minimal_counts":"perc"
     }
}
```


# Check results

#### Step 1: overview
Once completed, you can locate QC reports and QC plots for each run in the ```4.Output``` folder within the working directory. The ```QC_report.html``` is the main index file containing links to all QC files required for checking and trouble shooting library qualities for the last run, which is recorded in the ```lastest_run.txt```. Records of QC report for the previous runs can be found in the  ``` QC_report_history``` folder.

```
4.Output/
├── check_QC_PCA.html
├── counts
│   ├── raw_counts.csv
│   ├── raw_counts_unfiltered.csv
│   └── TPM_counts.csv
├── lastest_run.txt
├── QC_plots
├── QC_report.csv
├── QC_report_history
│   ├── Run_1
│   ├── Run_2
│   └── Run_3
└── QC_report.html
```

To get an idea about the overall quality of RNASeq result, you can refer to the recommendation and outlier columns in ```QC_report.html```. In general, you should find most of your samples labeled with 1.Good and No in those two columns, respectively. If not, there might be  critical issues in sample prep or sequencing steps that need to be diagnosed separately. The links to **PCA plot**, **QC parameters** and **QC plots** are located on the top left corner of the report table, while the sample specific QC reports are located on the right side of their corresponding sample IDs.

<img src='img/QC_reports.png'>



#### Step 2: Examine PCA

The interactive PCA plot allows you to quickly inspect meta information of each sample. If no meta table was provided in the json file, you can still check the recommendation and Outlier information in respect to each sample. To select meta information, click on the dropdown menu and select the target item. You can make specific category invisible by click its symbol at the bottom left corner.

<img src='img/PCA_plot.png' width = 800>



#### Step 3: Examine QC table

The QC table includes all important QC parameters and should be your main source of information for manual QC. We used fairly stringent criteria for ranking the sample quality, therefore it is pretty safe to ignore a sample marked as a good sample and not a outlier (in most cases a good sample should not be a outlier). If a category 2 or 3 sample was not in the good list but close enough due to the stringent criteria, you might not need to re-sequence it. Besides these scenarios, the sample need to be resequenced or the library should be redone. If not possible then the sample should be eliminated. Below is the dictionary of QC parameters used.   

|Parameter                 |Note                                                                                                                   |
|--------------------------|-----------------------------------------------------------------------------------------------------------------------|
|total_reads               |All paired end reads; each pair is considered one read                                                                 |
|filtered_reads            |fastp filtered reads with low complexity and low quality reads removed                                                 |
|filtered_reads_perc       |percentage of remaining reads after fastq QC                                                                           |
|adaptor_trimm_perc        |percentage of reads that have adaptor content and have been trimmed                                                    |
|dup_rate                  |percentage of fastq sequence that have duplicated reads                                                                |
|uniquely_mapped_reads     |reads only mapped to one location of the genome model                                                                  |
|uniquely_mapped_reads_perc|percentage of uniquely mapped reads                                                                                    |
|spliced_reads             |total number of splicing events in each read                                                                           |
|anno_spliced_reads        |Splicing known in splice junction database                                                                             |
|too_short_reads           |the overlap between reads and genome is less than the minimal set level; in STAR the default is 60% of seqeuence length|
|too_short_reads_perc      |percentage of those too short reads                                                                                    |
|exonic_perc               |percentage of reads mapped to exonic regions                                                                           |
|intronic_perc             |percentage of reads mapped to intronic regions                                                                         |
|intergenic_perc           |percentage of reads mapped to intergenic regions                                                                       |
|bias_5_prim               |mean expression of 5' divided by mean expression of transcript                                                         |
|bias_3_prim               |mean expression of 3' divided by mean expression of transcript                                                         |
|bias_5to3_prim            |the ratio between both 5' and 3' biases                                                                                |
|STAR_counts               |STAR generated count table sum by samples                                                                              |
|STAR_counts_perc          |STAR counts in total uniquely mapped reads                                                                             |
|t_rRNA_counts_perc        |percentage of tRNA and rRNA in the total STAR counts                                                                   |
|protein_coding_perc       |percentage of protein coding gene counts in the total STAR counts                                                      |
|pseudogene_perc           |percentage of  pseudogene counts in the total STAR counts                                                              |
|long-noncoding_perc       |percentage of long-noncoding mRNA counts in the total STAR counts                                                      |
|short-noncoding_perc      |percentage of short-noncoding mRNA counts in the total STAR counts                                                     |
|final_STAR_counts         |STAR count counts sum by sample excluding tRNA and rRNA counts                                                         |
|insert_mean               | mean insert size of paired end reads                                                                                  |
|insert_median             | median insert size of paired end reads                                                                                |
|Total_genes                 | Total number of genes identified after mapping        |

#### Step 4: Examine QC plots
The QC plots folder contains three types of figures that can assist your diagnoses of RNASeq sample quality in step 3.
- The first type is a simple visualization of QC parameters with dash lines mark where the set thresholds are compare to the real data (QC parameter in file names);

<img src='img/type_1_example.png' width = 1000>

In these figures, Good sample points were condensed to the left side for allowing more space for the others that need to be examined.  
- The second type is the bam coverage plots and each file represents one quality level (file name start with bamcoverage);  
- The third type contains information that shows additional QC measures and here are the details:


  **STAR_minimal_counts_soft_threshold:** Total number of genes recovered was plotted as a function of STAR counts, and saturation function was used to fit the data. Once fitted, the algorithm first detected the minimal STAR counts required for yielding defined gene recovery percentage and compared it to the set level. The higher value of them was then used as the final minimal STAR counts threshold for separating samples that need to be resequenced (if less than the threshold and have no other issues). After that, normal distribution parameters were calculated for total gene numbers of all other samples that are above the minimal STAR counts threshold, and we considered any sample fell below 95% confidence interval need to be manually QCed. The minimal percentage of gene recovery can be specified in the json config file.

  scenario 1: soft threshold was used

  <img src='img/Saturation_case1.png' width = 800>

  scenario 2: soft threshold less than set level and the later one was used

  <img src='img/Saturation_case2.png' width = 800>


  **Spearman_correlation:** Pairwise spearman correlation of all samples was calculated for detecting outlier. We assume the mean spearman value of each sample should be within the normal distribution of all sample means with 95% confidence (one tail). If not then that sample will be marked as 'outlier'.   

  <img src='img/SP_corr.png' width = 800>  

# Reference genome


The general rule of making reference index should follow [STAR's tutorial](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), and any change/addition/deletion of the reference files must be properly documented.

Here are the source of implemented reference files where we downloaded reference genome from:
  1. Human genome version GRCh37: **Full** genome reference downloaded from [GENCODE Release 19 (GRCh37.p13)](https://www.gencodegenes.org/human/release_19.html), Bed file downloaded from [Rseqc reference hg19](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg19_GencodeCompV19.bed.gz/download).

  2. Human genome version GRCh38: **Primary assemble** genome reference downloaded from [GENCODE Release 32 (GRCh38.p13)](https://www.gencodegenes.org/human/release_32.html), Bed file downloaded from [Rseqc reference hg38](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/hg38_Gencode_V28.bed.gz/download).

  3. Mouse genome version GRCm38:  **Primary assemble** genome reference downloaded from [GENCODE Release 23](https://www.gencodegenes.org/mouse/release_M23.html), Bed file downloaded from [Rseqc reference mm10](https://sourceforge.net/projects/rseqc/files/BED/Mouse_Mus_musculus/mm10_Gencode_VM18.bed.gz/download).


  Importantly, STAR recommended to exclude patches and alternative haplotypes from the genome model (see tutorial section 2.2.1), therefore you should only use the primary assemble for mapping. In case only full genome model is available, you can use the following python script or your own code to trim the genome file.

  - create remove_patch.py and put in the following code


```
# /usr/bin/env python3

# remove PATCH and Alternative haplotypes from fasta file
# usage remove_patch.py in.fasta > output.fasta


from Bio import SeqIO
import sys

in_fasta = sys.argv[1]
ffile = SeqIO.parse(in_fasta, "fasta")
header_pattern = ['PATCH','HSCHR']
for seq_record in ffile:
      if not any([i in seq_record.description for i in header_pattern]):
          print(seq_record.format('fasta'))
```

  - For making reference index; please make proper changes according to your settings

```  
#!/bin/bash
#PBS -N STAR_gen_37
#PBS -o out_STAR_gen_37
#PBS -e err_STAR_gen_37
#PBS -q default
#PBS -l nodes=1:ppn=4
#PBS -l mem=40gb
#PBS -l walltime=20:00:00
cd  /mnt/BioAdHoc/Groups/vd-ay/RNASeq_Workflow/Reference

# trim off batches and alternative haplotypes if needed
./remove_patch.py GRCH37.P13/GRCh37.p13.genome.fa > GRCH37.P13/GRCh37.p13.genome.primary_assembly.fa

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ./GRCH37.P13 --genomeFastaFiles ./GRCH37.P13/GRCh37.p13.genome.primary_assembly.fa --sjdbGTFfile ./GRCH37.P13/gencode.v19.annotation.gtf --sjdbOverhang 100
```

  Separately, an annotation file should be made for counting reads by gene type (gene_type_4) and TPM calculation in the pipeline. A example human GRCh37 annotation file can be downloaded [here](./files/GRCh37_annotation.csv). To make the annotation table, you will need to execute the following steps:

  1. Getting gene name and type for each ensembl ID. You can export annotations table for [GRCh37](https://grch37.ensembl.org/biomart/martview/923bfa0c7a1727c3fe634eb8c422df78) and [GRCh38](http://uswest.ensembl.org/biomart/martview/e859cf18a85550949a12ba09c8ab117c) from ensembl biomart. Mouse gene annotations are also available from these links.
  2. Merge gene types so 4 categories. The dictionary for merging of gene types can be downloaded [here](files/dict_gene_type.csv) .
  3. Getting exon length for each gene; this is a common output from most counting tools such as featureCounts and HTSeq, or you can download from a confident source.



# Change QC rules

Depending on the sequencing method and sample type, rules optimized for ranking the quality of RNA seq results may vary. There are two separated parts of the QC ranking algorithm, the logic and thresholds. In general, the logic section should stay put while the thresholds can be more flexible. 

To change QC threshold:

In the configuration file QC_threshold section, the sole number for each item is the minimal cutoff (maximal for "too_short_reads_perc" and "t_rRNA_counts_perc") for that parameter, exceptions are 1. "Total_genes" - range of minimal total genes counts; example here defines that if calculated threshold more than 9000 then use 9000. 2. "bias_5to3_prim" - lower and upper limit. 3. "insert_median" - range of lower and upper limit. 4. "minimal_counts" - 'fixed' or 'perc'; if 'fixed' then use "final_STAR_counts" as cut off, if 'perc' then the threshold is based on saturation curve.

To change QC rule logic:

All QC rule functions are defined in scripts/make_report.py function RNA_QC; make sure proper version control when editing this section here.
