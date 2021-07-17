
## PAGEANT: Personal Access to Genome & Analysis of Natural Traits

<br/>


## #1. The "ACGTU" Philosophy & The "5Q" Design

> ### The hard core of PAGEANT is a suite of common bioinformatics tools and python modules, combined to manage and annotate user provided genetic data. 
> ### The main scripts include GUI.py and main.py and a suite of extra libraries. They work together to generate a user interface, manage bioinformatics processes and data flow, and eventually generate an easy-to-read report. 
> ### The "ACGTU" Philosophy:
> > - #### Academic quality standard, where state-of-the-art algorithms and millions of genetic variants could be used to calculate individual Polygenic Risk Score (PRS); 
> > - #### Confidential data run locally, without the need of sending genomic data to cloud servers; 
> > - #### Generalizable architecture and algorithm, where our “five-Q” design could scale the genetic report from dozens of traits to hundreds and thousands of traits; 
> > - #### Transparent source code for all underlying programming scripts; 
> > - #### User centric, where users have the full control to add/remove certain traits into/from a genetic report. 

> ### The "5Q" Design:
> > - #### (1) Quality assurance and control report of genetic data itself; 
> > - #### (2) Qualitative assessment for genetic characteristics with absolute or relatively high certainty; 
> > - #### (3) Quantitative assessment of health risk susceptibility based on PRS and population reference; 
> > - #### (4) Query of third-party variants databases such as ClinVAR and PharmGKB; 
> > - #### (5) QR code display of genetic variants of interest, such as those variants for genetic identification and for certain drug prescription.

<br/>


## #2. Install & Run
> ## Download:
> > - ### Download executables for [Windows](https://drive.google.com/file/d/147zOn5b9dqeojVbGJq_rbLKZw24NSKe_/view?usp=sharing), [Linux](https://drive.google.com/file/d/1_WUJwMuf7EAsAyW6Q4hfHeB8eE2LjLrH/view?usp=sharing), [Mac OS](https://drive.google.com/file/d/1njO2AKC8Z6PcwN1Zh6s6sVN9NUi32gfc/view?usp=sharing)
> > - ### Download [HapMap3 genotype](https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3), as population reference by default. LifeOver is not needed, because only SNP rsID is used by PAGEANT.
> > - ### Download [1000 genomes project (G1K) genotype](https://www.internationalgenome.org), to be used as population reference. This is only needed when the user genotype data is based on G1K imputation.
> ## Run:
> > - ### For Windows OS: the program could be run directly by double clicking "PAGEANT.exe".
> > - ### For Mac-OS:  run "[brew](https://brew.sh/) install zbar" first  to install necessary libraries, and then double click "PAGEANT".
> > - ### For Linux: the program could be run directly by typing "./PAGEANT".

<br/>


## #3. Example reports 

> ## QA/QC using principal components

![example report for QC](./images/Fig_PC.png)

> ## Quantitative traits

![example report for Quantitative traits](./images/Fig_Qt.png)

> ## QR code

![example report for QR code](./images/Fig_QR.png)

<br/>


## #4. Customize

> - ### The folder structures of PAGEANT is shown below. Advanced users could also follow this structure to customize the genetic report. 
> - ### For example, under "algorithm database" folder, there are 3 files for each trait folder: TRAIT.desc.txt for description text, TRAIT.jpg for a representative picture, TRAIT.snps.ref for a list of SNPs used and the relevant calculation rules. 
> - ### Users could follow this structure to add many new traits to be included in the genetic report.  

![Folder Structure](./images/Fig_folder.png)

<br/>


## Contact & Cite

> - ### [Jie Huang](jiehuang001@pku.edu.cn) MD MPH PhD, Department of Global Health, Peking University School of Public Health
> - ### Manuscript under review

