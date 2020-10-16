
# PAGEANT (Personal Analysis of Genome and Annotation Toolkit ), is an open source platform that convert a personal genome into a health report.
PAGEANT differs from other similar tools, summarized by “A-C-G-T”: 
(1). Academic quality standard, where start-of-the-art algorithms and millions of genetic variants could be used to calculate PRS; 
(2). Confidential data run locally, no need to send genomic data to cloud servers; 
(3). Generalizable health assessment, easy to add/remove scope of health reports based on users' preference and comfort level. 
(4). Transparent source code for all underlying programming scripts..


Author: Zhisheng Liang MS, Jie Huang MD PhD, Department of Global Health, Peking University School of Public Health



# Demo:

# #1 下载参考人群的基因组数据
用户可以在[ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/)下载公开的千人基因组数据，然后作为报告分析做需要的参考人群基因组数据放在/Referenc目录下

# #2 进行分析
打开PAGEANT.exe，使用默认参数即可进行分析

# #3 查看报告
分析成功后，点击OK即可查看报告


# Steps:
# #1. Download 1000 genomes genotype data to use as references
start from 1000 genomes project main page https://www.internationalgenome.org. 
Then Click the "EBI FTP site" link under the "Alignments" section, and click "1000_genomes_project" link on the next page.
Users will directed to http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/. 
The "1000genomes.exome.GRCh38DH.alignment.index" file listed the FTP URL for 2,692 samples. 
The New York Genome Center (NYGC) released high-coverage (30x) data for a total of 3,202 samples. 
Users could download the aligned sequencing data for any set of samples from this link https://www.internationalgenome.org/data-portal/data-collection/30x-grch38, aligned to the GRCh38 reference genome. Once the CRAM file is downloaded, users could use samtools to extract certain regions of the genome to created a much smaller dataset, by using scripts such as below:

```

# create a 2genes.bed file with the following two rows, tab separated.
chr1    159204012       159206500       ACKR1
chr19   44905781        44909393        APOE

# run SAMTOOLS to extract the two genomic regions and create a new 2gene.bam file
samtools view -L 2genes.bed -O BAM -o 2genes.bam [raw-cram-file]

```
