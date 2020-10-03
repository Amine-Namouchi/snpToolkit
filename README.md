
Status: This tutorial is being updated, so keep tuned!

![img/snpToolkit_LOGO.png](img/snpToolkit_LOGO.png)

# About

**snpToolkit** is a computational framework written in Python 3. snpToolkit allow users to:

1. Visualize the content of their VCF files
2. Filter SNPs based on multiple criteria:
    - Coordinates of regions to exclude
    - Depth of coverage
    - Quality
    - The ratio corresponding to the number of reads that have the mutated allele / total number of reads at that particular position
3. Annotate SNPs using genome annotation data provided within a genbank file
4. Extract the distribution of all indels according to genome annotation
5. Visualize and explore the annotated SNPs for all analyzed files
6. Combine all snpToolkit output files generated using the annotate option and generate:
    - A table storing the distribution of all SNPs on each sample
    - A fasta file with all concatenated SNPs for each sample. such file can be used to build a phylogenetic tree.

snpToolkit detects automatically if the input vcf files were generated using samtools mpileup, gatk HaplotypeCaller or freebayes. Vcf files could be gzipped or not.

# How to install and update

Different python libraries need to be installed: Biopython, pysam,  pandas, plotly, dash,  tqdm and coloredlogs. The recommended way to install snptoolkit is:

```bash
pip install  git+git://github.com/Amine-Namouchi/snpToolkit.git
```

very soon it will be possible to install snptoolkit directly from pypi as follow: pip install snptoolkit



# snpToolkit menu

```bash
snptoolkit -h
usage: snptoolkit [-h] {explore,annotate,combine,viz}

    snpToolkit takes vcf files, as well as bam files (optional) as inputs. The vcf files could be generated using samtools/bcftools, 
		gatk HaplotypeCaller or freeBayes.
    Please visit https://github.com/Amine-Namouchi/snpToolkit for more information.

positional arguments:
  {explore,annotate,combine,viz}
                        commands
    explore             explore your vcf files before annotation
    annotate            Please provide one or multiple vcf files
    combine             combine snpToolkit output files in one alignment in fasta format
    viz                 visualize snptoolkit output files

optional arguments:
  -h, --help            show this help message and exit

GPLv3 licence | Amine Namouchi | amine.namouchi@gmail.com
```

For this tutorial we will have 10 vcf files to analyze and compare named sample1.vcf.gz to sample10.vcf.gz. In this example we will work with *Yersinia pestis* so we need the genbank file corresponding to the reference strain used in the alignment process to generate the vcf files. In this case we need the file GCF_000009065.1_ASM906v1_genomic.gbff that can be obtained from [NCBI.](https://www.ncbi.nlm.nih.gov/) 

# The explore command

```bash
snptoolkit explore -h
usage: snptoolkit explore [-h] -i IDENTIFIER

optional arguments:
  -h, --help     show this help message and exit

snpToolkit explore required options:
  -i IDENTIFIER  Provide the input vcf files
```

This command allows user to explore the SNPs on each of their vcf files. 

The option -i  allows to specify a common identifier in the vcf files names. If you want to explore all VCF files in a folder, you can use vcf  as identifier as it is present in all vcf file names (usually filename.vcf.gz). On the contrary, if you have added in the filenames of your vcf files, for example, the years of isolation of each sample, you can use the year you want as identifier. 

when you run the command:

```bash
$ snptoolkit explore -i vcf
[TIME][INFO] [snpToolkit is extracting your data and creating the different plots...]
progress: 100%|#########################################################################| 9/9 [00:00<00:00, 83.75it/s]
Dash is running on http://127.0.0.1:8050/

 * Serving Flask app "explore_snpToolkit" (lazy loading)
 * Environment: production
 * Running on http://127.0.0.1:8050/ (Press CTRL+C to quit)
```

snptoolkit will analyze all raw data on each VCF file in terms of SNPs and starts a web application that you access using the link mentioned above http://127.0.0.1:8050. For this example of 10 vcf files, it took less than a second to analyze all files. Figure 1 shows a screenshot of the generated dashboard to explore your data.

![img/Figure1.png](img/Figure1.png)

**Figure1**

For sample 5 for example, we can see that the total number of SNPs in the chromosome NC_003143.1 is 399 SNPs. This is the total raw number. Lets detail each column of the first table:

- If we apply just the depth filter (-d) when using the option annotate (see below), only 10 SNPs will be excluded as they have a coverage less than 3.
- If we consider 20 as a cutoff for the quality of each SNPs, the number drop to 356 SNPs
- If we only consider those that have a ratio (nb reads with mutated allele/total number of reads on that position) ≥ 0.9, the number of SNPs drops to 230.
- If all filters are used: depth ≥3, QUAL ≥20 and ratio≥0.9, the number of filtered SNPs will be equal to 215.

For the case of *Yersinia pestis,* there are 3 plasmids. For sample 5, there are SNPs on plasmid NC_003134.1 and NC_003131.1

The first plot in Figure1 shows the distribution of all SNPs based on Ratio (x axis) and Depth (y axis). The size of each circle is proportional to the quality of each SNP. The second plot complement the first plot as it give you an idea about the proportion of SNPs for the chromosome and each of the plasmids. For the chromosome NC_004143.1, we can see that there is a small proportion of SNPs located between 0.2 and 0.4, but most of the SNPs has a high ratio ≥ 0.9. 

To hide any of the data presented on each plot, you just need to select the name that you want. 

# The annotate command

```bash
snptoolkit annotate -h
usage: snptoolkit annotate [-h] -i IDENTIFIER -g GENBANK [-p PROCESSORS] [-f EXCLUDECLOSESNPS] [-q QUALITY] [-d DEPTH] [-r RATIO] [-e EXCLUDE]

optional arguments:
  -h, --help           show this help message and exit

snpToolkit annotate required options:
  -i IDENTIFIER        provide a specific identifier to recognize the file(s) to be analyzed
  -g GENBANK           Pleae provide a genbank file

snpToolkit annotate additional options:
  -p PROCESSORS        number of vcf files to be annotated in parallel default value [1]
  -f EXCLUDECLOSESNPS  exclude SNPs if the distance between them is lower then the specified window size in bp
  -q QUALITY           quality score to consider as a cutoff for variant calling. default value [20]
  -d DEPTH             minimum depth caverage. default value [3]
  -r RATIO             minimum ratio that correspond to the number of reads that has the mutated allele / total depth in that particular position. default
                       value [0]
  -e EXCLUDE           provide a tab file with genomic regions to exclude in this format: region start stop. region must correspond to the same name(s) of
                       chromsome and plasmids as in the genbank file
```

This command allows to filter and annotate all SNPs from each selected VCF files. Only two options are required:

1. option -i: The user need to specify a common identifier found on all VCF files he wants to analyze. If only one VCF file is to be analyzed, provide the file name. If all VCF files should be analyzed, the user  needs to provide e.g vcf as all vcf files will have at the end .vcf.gz of .vcf. 
2. -g: genbank file. The genbank file must include the fasta sequence for the chromosome and plasmids, if any. genbank files can be downloaded from NCBI

Several options are additional and are needed to filter SNPs:

- -f: To be able to exclude all SNPs that can be located in hotspot zones or short repeats, it is possible to specify an integer that will correspond to the minimum of distance between SNPs to be kept. if the distance between two SNPs is lower than the specified cutoff, both SNPs will be ignored.
- -q: Quality score to consider as a cutoff for variant calling. The default value is 20
- -d: Minimum depth coverage. The default value is 3.
- -r: $r = m/t$  where *m* is the umber of reads that carry the mutated allele and *t* is the total number of read on that position. If not specified all SNPs will be taken into account.
- -e: This is to specify a tab delimited file with the coordinates of the regions to be ignored when annotating SNPs. If we take the example of the genbank used in for this tutorial:

    ```bash
    $ grep 'LOCUS' /Users/amine/Documents/tutorials/snptoolkit/GCF_000009065.1_ASM906v1_genomic.gbff
    LOCUS       NC_003143            4653728 bp    DNA     circular CON 20-MAR-2020
    LOCUS       NC_003131              70305 bp    DNA     circular CON 20-MAR-2020
    LOCUS       NC_003134              96210 bp    DNA     circular CON 20-MAR-2020
    LOCUS       NC_003132               9612 bp    DNA     circular CON 20-MAR-2020
    ```

    as you can see there is one chromosome NC_003143 and three plasmids: NC_003131, NC_003134 and NC_003132. The tab delimited file should  look as follows:

    ```bash
    NC_003143.1	4016	4079
    NC_003143.1	7723	7758
    NC_003143.1	11562	19149
    NC_003143.1	25663	26698
    ```

    If there are regions on the plasmids sequences you can also add them in the same file.

Now time to run the annotate command

```bash
$ snptoolkit annotate -i vcf -g GCF_000009065.1_ASM906v1_genomic.gbff -d 5 -q 30 -r 0.9 -p 4

[15:37:30] [INFO] [4 CPUs requested out of 8 detected on this machine]
[15:37:30] [INFO] [snpToolkit is filtering and annotating your SNPs]
 89%|████████████████████████████████████████████████████████████████████████████████████| 9/9 [00:01<00:00,  2.67it/s]
[15:37:32] [INFO] [snpToolkit output files will be located in folders
				 snpToolkit_SNPs_output_... 
				 and snpToolkit_INDELS_output_...]
100%|████████████████████████████████████████████████████████████████████████████████████| 9/9 [00:02<00:00,  3.95it/s]
```

snptoolkit generates two folders with the  date and time stamp, one for SNPs and one for indels:

```bash
├── snpToolkit_INDELS_output_...
│   ├── sample3_snpToolkit_indels.txt
│   ├── sample9_snpToolkit_indels.txt
│   ├── sample10_snpToolkit_indels.txt
│   ├── sample1_snpToolkit_indels.txt
│   ├── sample2_snpToolkit_indels.txt
│   ├── sample4_snpToolkit_indels.txt
│   ├── sample5_snpToolkit_indels.txt
│   ├── sample6_snpToolkit_indels.txt
│   ├── sample7_snpToolkit_indels.txt
│   └── sample8_snpToolkit_indels.txt
├── snpToolkit_SNPs_output_...
    ├── sample3_snpToolkit_SNPs.txt
    ├── sample9_snpToolkit_SNPs.txt
    ├── sample10_snpToolkit_SNPs.txt
    ├── sample1_snpToolkit_SNPs.txt
    ├── sample2_snpToolkit_SNPs.txt
    ├── sample4_snpToolkit_SNPs.txt
    ├── sample5_snpToolkit_SNPs.txt
    ├── sample6_snpToolkit_SNPs.txt
    ├── sample7_snpToolkit_SNPs.txt
    └── sample8_snpToolkit_SNPs.txt
```

All generated output files are tab delimited. 

### Example of SNP output file:

```bash
##snpToolkit=version 2.2.1
##commandline= snptoolkit annotate -i vcf -g GCF_000009065.1_ASM906v1_genomic.gbff -d 5 -q 30 -r 0.9 -p 4
##VcfFile=sample5.vcf.gz
##Total number of SNPs before snpToolkit processing: 406
##The options -f and -e were not used
##Filtred SNPs. Among the 406 SNPs, the number of those with a quality score >= 30, a depth >= 5 and a ratio >= 0.9 is: 218
##After mapping, SNPs were located in:
##NC_003131.1: Yersinia pestis CO92 plasmid pCD1, complete sequence 70305 bp
##NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
##The mapped and annotated SNPs are distributed as follow:
##Location      Genes   RBS     tRNA    rRNA    ncRNA   Pseudogenes     intergenic      Synonymous      NonSynonumous
##SNPs in NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp 155     0       0       1       0       0       57      54      101
##SNPs in NC_003131.1: Yersinia pestis CO92 plasmid pCD1, complete sequence 70305 bp    2       0       0       0       0       0       3       1       1
##Syn=Synonymous NS=Non-Synonymous
##Coordinates   REF     SNP     Depth   Nb of reads REF Nb reads SNPs   Ratio   Quality Annotation      Product Orientation     Coordinates in gene     Ref codon       SNP codon       Ref AA  SNP AA  Coordinates protein     Effect  Location
82      C       A       36      0       34      1.0     138.0   intergenic      .       +       .       -       -       -       -       -       -       NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
130     G       C       28      0       27      1.0     144.0   intergenic      .       +       .       -       -       -       -       -       -       NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
855     G       A       69      0       62      1.0     228.0   YPO_RS01010|asnC        transcriptional regulator AsnC  -       411     ACC     AC[T]   T       T       137     Syn     NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
```

The first lines of the snptoolkit file for SNPs contain a summary and useful information. The SNPs annotation is organized in tab delimited table. The columns of this table are:

1. Coordinates: SNP coordinate 
2. REF: Reference allele
3. SNP: new allele in analyzed sample
4. Depth: total depth of coverage 
5. Nb of reads REF: number of reads with the reference allele
6. Nb reads SNPs: number of reads with the new allele
7. Ratio:  Nb reads SNPs/(Nb of reads REF+Nb reads SNPs)
8. Quality: quality score
9. Annotation: distribution within genes or intergenic
10. Product: functional product of the gene
11. Orientation: gene orientation
12. Coordinates in gene: coordinate of the SNP within the gene
13. Ref codon: reference codon, ACC in the example above
14. SNP codon: new codon, AC[T]
15. Ref AA: Amino Acid corresponding to reference codon 
16. SNP AA: Amino Acid corresponding to new codon
17. Coordinates protein: coordinate of the Amino Acid 
18. Effect: could be Synonymous (Syn) or Non-Synonymous (NS)
19. Location: ID of the chromosome and plasmids.