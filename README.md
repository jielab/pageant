
# PAGEANT: Personal Access to Genome and Analysis of Natural Traits

Contact: jiehuang001@pku.edu.cn (Jie Huang MD MPH PhD, Department of Global Health, Peking University School of Public Health)

<br />


# 1. Core functionalities (5Qs)

The hard core of PAGEANT is a suite of common bioinformatics software including VEP and PLINK to manage and annotate user provided genetic data. 
The main python script is used to generate user interface, manage the process and data flow, and eventually generate an easy-to-read report. 
<br />
<br />

# 2. Install and click to run

* Three zipped folders are downloadable from [Google Drive](https://drive.google.com/drive/folders/1utGpJNofmjqoV6TG8F9FqMv9iD-CKhwi?usp=sharing) for Linux, Mac, Windows, respectively. The hard core scripts include GUI.py and main.py and a suite of extra libraries.
* For Windows OS, the program could be run directly by double clicking "PAGEANT.exe".
* For MAC OS, please use brew (https://brew.sh/) to run "brew install zbar" to install necessary libraries first, and then double click "PAGEANT".
* For Linux, the program could be run directly by typing "./PAGEANT".
<br />
<br />

# 3. Example report of QA/QC

![example report for QC](./pictures/Fig_QC.png)
<br />

# 4. Example report of Quantitative traits

![example report for Quantitative traits](./pictures/figure6.png)
<br />

# 5. Customize PAGEANT

### 5.1 Download HapMap3 genotype, used as population reference by default

```
Open https://www.broadinstitute.org/medical-and-population-genetics/hapmap-3， 
click the 3 links under "A. SNP Genotype Data" in the section of "How To Download This Release".

Advanced users could use liftOver to convert the GRCh36 based .map file into GRCh37 based.
But this is not needed for PAGEANT, since only SNP rsID is used for querying and calculation.

```

## 5.2 download 1000 genomes project (1000G) genotype data, if needed

```

open https://www.internationalgenome.org, click "Data" menu on the top.
under "Available data" section, click "Phase 3" VCF files.

use wget to download files at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# rename the long file names to something short such as below:
mv ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz chr1.vcf.gz

# rename chromosome to add "chr" prefix, if needed
seq 1 22 | xargs -n1 -I % echo % chr% > chr_name_conv.txt
bcftools annotate chr1.vcf.gz --rename-chrs chr_name_conv.txt -Oz -o chr1.vcf.gz

# convert from build 37 positions to build 38 positions, if needed
gatk --java-options "-Xmx6g" LiftoverVcf -R Homo_sapiens_assembly38.fasta.gz -I chr1.vcf.gz -O chr1.b38.vcf.gz \
	-C hg19ToHg38.over.chain --REJECT rejected.vcf --DISABLE_SORT true

# create a small VCF subset, keeping only those samples of interest.
echo "NA20525" > sample.keep
plink2 --vcf chr1.vcf.gz --extract subset.snps --keep sample.keep --export vcf bgz id-paste=iid --out chr1.subset

```

## 5.3 Add or remove traits from the genetic report

Please follow the following design.

![Figure 3](./pictures/figure3b.png)
