
The annotate command
====================

Done: Yes

.. code-block:: bash

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

This command allows to filter and annotate all SNPs from each selected VCF files. Only two options are required:


#. option -i: The user need to specify a common identifier found on all VCF files he wants to analyze. If only one VCF file is to be analyzed, provide the file name. If all VCF files should be analyzed, the user  needs to provide e.g vcf as all vcf files will have at the end .vcf.gz of .vcf. 
#. -g: genbank file. The genbank file must include the fasta sequence for the chromosome and plasmids, if any. genbank files can be downloaded from NCBI

Several options are additional and are needed to filter SNPs:


* -f: To be able to exclude all SNPs that can be located in hotspot zones or short repeats, it is possible to specify an integer that will correspond to the minimum of distance between SNPs to be kept. if the distance between two SNPs is lower than the specified cutoff, both SNPs will be ignored.
* -q: Quality score to consider as a cutoff for variant calling. The default value is 20
* -d: Minimum depth coverage. The default value is 3.
* -r: $r = m/t$  where *m* is the umber of reads that carry the mutated allele and *t* is the total number of read on that position. If not specified all SNPs will be taken into account.
* 
  -e: This is to specify a tab delimited file with the coordinates of the regions to be ignored when annotating SNPs. If we take the example of the genbank used in for this tutorial:

  .. code-block:: bash

       $ grep 'LOCUS' /Users/amine/Documents/tutorials/snptoolkit/GCF_000009065.1_ASM906v1_genomic.gbff
       LOCUS       NC_003143            4653728 bp    DNA     circular CON 20-MAR-2020
       LOCUS       NC_003131              70305 bp    DNA     circular CON 20-MAR-2020
       LOCUS       NC_003134              96210 bp    DNA     circular CON 20-MAR-2020
       LOCUS       NC_003132               9612 bp    DNA     circular CON 20-MAR-2020

    as you can see there is one chromosome NC_003143 and three plasmids: NC_003131, NC_003134 and NC_003132. The tab delimited file should  look as follows:

  .. code-block:: bash

       NC_003143.1 4016    4079
       NC_003143.1 7723    7758
       NC_003143.1 11562   19149
       NC_003143.1 25663   26698

    If there are regions on the plasmids sequences you can also add them in the same file.

Run
===

Now time to run the annotate command

.. code-block:: bash

   $ snptoolkit annotate -i vcf -g GCF_000009065.1_ASM906v1_genomic.gbff -d 5 -q 30 -r 0.9 -p 4

   [15:37:30] [INFO] [4 CPUs requested out of 8 detected on this machine]
   [15:37:30] [INFO] [snpToolkit is filtering and annotating your SNPs]
   100%|████████████████████████████████████████████████████████████████████████████████████| 9/9 [00:01<00:00,  2.67it/s]
   [15:37:32] [INFO] [snpToolkit output files will be located in folders
                    snpToolkit_SNPs_output_... 
                    and snpToolkit_INDELS_output_...]
   100%|████████████████████████████████████████████████████████████████████████████████████| 9/9 [00:02<00:00,  3.95it/s]

Outputs
=======

snptoolkit generates two folders with the  date and time stamp, one for SNPs and one for indels:

.. code-block:: bash

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

All generated output files are tab delimited. 

Example of SNP output file:
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

   ##snpToolkit=__version__
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

The first lines of the snptoolkit file for SNPs contain a summary and useful information. The SNPs annotation is organized in tab delimited table. The columns of this table are:


#. Coordinates: SNP coordinate 
#. REF: Reference allele
#. SNP: new allele in analyzed sample
#. Depth: total depth of coverage 
#. Nb of reads REF: number of reads with the reference allele
#. Nb reads SNPs: number of reads with the new allele
#. Ratio:  Nb reads SNPs/(Nb of reads REF+Nb reads SNPs)
#. Quality: quality score
#. Annotation: distribution within genes or intergenic
#. Product: functional product of the gene
#. Orientation: gene orientation
#. Coordinates in gene: coordinate of the SNP within the gene
#. Ref codon: reference codon, ACC in the example above
#. SNP codon: new codon, AC[T]
#. Ref AA: Amino Acid corresponding to reference codon 
#. SNP AA: Amino Acid corresponding to new codon
#. Coordinates protein: coordinate of the Amino Acid 
#. Effect: could be Synonymous (Syn) or Non-Synonymous (NS)
#. Location: ID of the chromosome and plasmids. 

IMPORTANT NOTE:  

In the example above, the total depth for the first SNP is 36, while the number of reads that carry the reference allele plus the number of reads that carry the new allele is equal to 34. The VCF file corresponding to that sample is generated using samtools mpileup. By default, samtools mpileup applies Phred-scaled probability of a read base being misaligned, known as BAQ. As indicated in samtools documentation, this greatly helps to reduce false SNPs caused by misalignments. 

The total depth shown by snpToolkit is the raw depth taking into account all reads (column 4). However, the columns 5 and 6 shows the number of reads with  Phred-scaled probability. The ratio in column 7 is based on only column 5 and 6.  I have made this decision to store as much information as possible from the original VCF file. 

If the VCF files where produced using samtools-mpileup with the option -B to skip Phred-scaled probability, you will not see such difference. 

Example of INDELS output file:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The indels output is in tab delimited format as follows:

.. code-block:: bash

   55732   CCGGGGCGGGGCGGGGCGG CCGGGGCGGGGCGG  62  0   20  1.0 228.0   intergenic  .   .   deletion    5   NC_003131.1: Yersinia pestis CO92 plasmid pCD1, complete sequence 70305 bp
   35188   T   TTC 41  0   32  1.0 228.0   intergenic  .   .   insertion   2   NC_003134.1: Yersinia pestis CO92 plasmid pMT1, complete sequence 96210 bp
   73418   GAA GA  72  0   68  1.0 228.0   intergenic  .   .   deletion    1   NC_003134.1: Yersinia pestis CO92 plasmid pMT1, complete sequence 96210 bp
   16  AC  A   13  0   13  1.0 149.0   intergenic  .   .   deletion    1   NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
   183029  CCAATAACAAT CCAATAACAATAACAAT   95  0   24  1.0 228.0   intergenic  .   .   insertion   6   NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
   266466  AGGGGGGGG   AGGGGGGGGG  40  1   25  0.96    66.0    CDS YPO_RS02340|YPO_RS02340 EscV/YscV/HrcV family type III secretion system export apparatus protein    insertion   1   NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
   552919  TGGGGGGG    TGGGGGGGG   93  0   71  1.0 122.0   CDS YPO_RS03585|tssM    type VI secretion system membrane subunit TssM  insertion   1   NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
   581519  GTTCAATTCAATTCAAT   GTTCAATTCAATTCAATTCAAT  31  0   9   1.0 228.0   intergenic  .   .   insertion   5   NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
   747924  AGGGGGGGG   AGGGGGGGGG  41  1   26  0.96    71.0    CDS YPO_RS04395|YPO_RS04395 pseudopilin insertion   1   NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp
   813977  GC  GCCTGGCCATC 54  0   10  1.0 228.0   CDS YPO_RS04755|YPO_RS04755 DASS family sodium-coupled anion symporter  insertion   9   NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp

for the case of the position 266466 for example

.. code-block:: bash

   266466  AGGGGGGGG   AGGGGGGGGG  40  1   25  0.96    66.0    CDS YPO_RS02340|YPO_RS02340 EscV/YscV/HrcV family type III secretion system export apparatus protein    insertion   1   NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp

The different columns are:


#. Coordinates (266466)
#. Reference (AGGGGGGGG)
#. Sample (AGGGGGGGGG)
#. Number of total reads (40)
#. Number of reads with reference sequence (1)
#. Number of reads with new sequence (25)
#. Ratio (0.96)
#. Quality score (66.0)
#. Location (CDS)
#. Gene or intergenic (YPO_RS02340|YPO_RS02340)
#. Gene product (EscV/YscV/HrcV family type III secretion system export apparatus protein)
#. Type of indel (insertion)
#. Number of nucleotide (1)
#. Sequence name (NC_003143.1: Yersinia pestis CO92, complete genome 4653728 bp)

While snpToolkit annotate indels, it is important to be careful and check any indels you are interested in before to elaborate any hypothesis and conclusion.
