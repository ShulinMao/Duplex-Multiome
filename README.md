# Duplex-Multiome

Understanding the role of somatic mutations in cancer, aging, and healthy tissues requires high-accuracy methods that can characterize somatic mosaicism across diverse cell types while linking mutations to their functional consequences. Existing approaches for studying cell-type-specific somatic mutations within tissues are often low-throughput and prohibitively expensive. To address these challenges, we developed **Duplex-Multiome**, a novel technology that integrates **duplex consensus sequencing** with **single-nucleus ATAC-seq (snATAC-seq) and RNA-seq (snRNA-seq)** to simultaneously profile somatic single-nucleotide variants (sSNVs) and gene regulatory landscapes from the same nucleus.

![Duplex-Multiome Overview](images/Duplex_Multiome_schematic.png)

## Experimantal Protocol

We selected the 10X Single Cell Multiome ATAC + Gene Expression protocol as a foundation, a widely adopted method for profiling chromatin accessibility and gene expression across diverse tissues. However, a major challenge in applying this protocol to somatic mutation detection is the sequencing error rate, which ranges from 1 in $10^2$ to 1 in $10^3$ on standard sequencing platforms. This poses a great obstacle to somatic mutation calling in 10X Multiome data, especially in non-malignant tissues in which mutation rates are 1 in $10^7$ or even lower. Duplex consensus sequencing, in which products generated from both strands of an initial DNA molecule are sequenced independently, has emerged as the most effective means of reducing sequencing error. We therefore developed Duplex-Multiome, incorporating duplex consensus sequencing through a strand-tagging strategy to expand the capabilities of the 10X Multiome protocol.

Please see **our paper (need a link to the paper in the future)** for the detailed protocol for preparing Duplex-Multiome libraries.

![Duplex-Multiome Protocol](images/Protocol_schematic.png)

## Computational Pipeline  

Our pipeline processes and analyzes Duplex-Multiome data efficiently. The workflow includes:  

1. **Data Preprocessing**  
    - Read alignment  
    - Multiome data quality control (QC), normalization, batch effect correction, and integration  

2. **Variant Calling**  
    - Detection of somatic SNV candidates using duplex sequencing  
    - Error correction and quality filtering  
    - Recovery of somatic clonal SNVs from stringent filtering  

3. **Multiomic Integration**  
    - Multiome data cell-type annotation  
    - Linking mutations with snATAC-seq and snRNA-seq profiles  

4. **Mutational Pattern Analyses**  
    - Characterization of cell-type-specific mutational burden and spectra  
    - Exploration of clonal mutation landscapes and cell lineages  
    - Profiling the association between clonal mutations and gene expression changes

## Getting Started 

### Set up the computing environment
The analysis requires a SLURM computing cluster.

#### 1. Clone Our GitHub Repository  

To get started, clone this repository to your local machine.  
 
```bash
git clone https://github.com/ShulinMao/Duplex-Multiome.git
```

Only the `pipeline` and `scripts` folders will be used for following analyses. You could delete other
folders if want.

#### 2. Set Up the Environment with Mamba  
We use **Mamba** for environment management to ensure efficient and reproducible installation of dependencies.

If you don't have Mamba installed, you can install it via Conda or follow the instructions from :  

```bash
conda install -n base -c conda-forge mamba
```

Alternatively, you can install Micromamba, a lightweight version of Mamba:
```bash
curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
```

Once Mamba is installed, create a new environment for the Duplex-Multiome pipeline:
```bash
cd "path/to/your/repository"
mamba create env create -f ./pipeline/pipeline_preprocessing/environment.yml
mamba activate Duplex-Multiome
```

#### 3. Configure the R environment
We have included **R (v4.0.2)** in the Mamba environment, which is the R version used for all analyses. If you have completed the previous installation steps successfully, R should already be installed on your machine.  

To avoid conflicts between different R packages, we need to set up three separate libraries.  

##### 3.1 Install R packages required for most analyses

Create a folder for the packages in 3.1
```bash
cd "path/to/your/repository"
mkdir Rpackages
```
Start R and specify the package installation path to the newly created directory:
```R
.libPaths(c("path/to/your/repository/Rpackages/"))
```
Next, install the following R packages
- Seurat	v4.0.5
- SeuratDisk	v0.0.0.9019
- Signac	v1.5.0
- GenomicRanges	v1.40.0
- VariantAnnotation	v1.34.0
- future	v1.18.0
- harmony	v0.1.0
- EnsDb.Hsapiens.v86	v2.99.0
- BSgenome.Hsapiens.UCSC.hg38	v1.4.3
- truncnorm v1.0-9
- ggplot2	v3.3.2
- ggpubr	v0.4.0

##### 3.2 Install R packages required for mutational spectrum analysis

Create a folder for the packages in 3.2
```bash
cd "path/to/your/repository"
mkdir Rpackages_sig
```
Start R and specify the package installation path to the newly created directory:
```R
.libPaths(c("path/to/your/repository/Rpackages_sig/"))
```
Next, install the following packages in R
- MutationalPatterns	v3.0.1
- NMF	v0.25
- gridBase	v0.4-7

##### 3.3 Install R packages required for clonal mutation analysis

Create a folder for the packages in 3.3
```bash
cd "path/to/your/repository"
mkdir Rpackages_rstan
```
Start R and specify the package installation path to the newly created directory:
```R
.libPaths(c("path/to/your/repository/Rpackage_rstan/"))
```
Next, install the following packages in R
- rstan v2.32.6

#### 4. Download reference genome and other required materials
The analyses are performed using GRCh38 reference genome.
```bash
wget <link>
```

We use gnomAD v3.1.2 and dbSNP 147 (common SNPs) in data preprocessing and filtering. While we used these specific versions in our paper, we recommend using the latest versions for your own analyses, as they are continuously updated with improved variant annotations.
```bash
# We downloaded gnomAD data using ANNOVAR
# If you don't have ANNOVAR installed, follow the instructions on https://annovar.openbioinformatics.org/en/latest/user-guide/download/

<path to ANNOVAR>/annotate_variation.pl -downdb gnomad312_genome humandb -buildver hg38
# This will download a file named hg38_gnomad312_genome.txt
```

---
### Bulk WGS data processing

#### 1. WGS data alignment
We perform alignment using BWA v0.7.83 with the reference genome GRCh38. Data preprocessing follow the GATK Best Practices (https://gatk.broadinstitute.org/hc/en-us/articles/360035535912-Data-pre-processing-for-variant-discovery). 

Intermediate BAM files from the alignment are sorted using the SortSam function of Picard v2.26.10, and duplicated reads are marked using MarkDuplicates from Picard v2.26.10. Base quality score recalibration is then performed using Genome Analysis Toolkit (GATK) v4.3.0.0.

We provide a template script to conduct bulk WGS data processing under `pipeline/bulk_WGS_data_processing/alignment.sh`. You can explore and modify our scripts. If you need to run the pipeline on your own data, simply update the paths to match your source files.

Note: Running this script is not mandatory. You may use your own pipeline as long as it follows GATK Best Practices.

#### 2. Call high-confidence germline mutations
We use GATK HaplotypeCaller to identify germline mutation candidates from the corresponding WGS BAM file generated after the above step. To obtain high-confidence germline mutations, we:  
1. Extract heterozygous candidates.  
2. Retain candidates that are also present in dbSNP 147 (common SNPs).  

These variants serve as high-confidence germline mutations for a given sample. We will later use these variants to estimate mutation detection sensitivity.  


We also provide a template script to conduct bulk WGS data processing under `pipeline/bulk_WGS_data_processing/call_high_confidence_germline_mutations.sh`. You can explore and modify our script for your own analysis.

---
### Multiome data preprocessing



---
### Somatic SNV candidates calling

---
### Recovery of somatic clonal SNVs

---
### Mutational Pattern Analyses
