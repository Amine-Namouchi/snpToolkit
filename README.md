
Please visit the full documentation at https://snptoolkit.readthedocs.io/en/latest/index.html 

![img/snpToolkit_LOGO.png](img/snpToolkit_LOGO.png)

**snpToolkit** is a computational framework written in Python 3. snpToolkit allows users to:


1. Visualize the content of their VCF files.
2. Filter SNPs based on multiple criteria:

   * Coordinates of regions to exclude
   * Depth of coverage
   * Quality
   * The ratio corresponding to the number of reads that have the mutated allele / total number of reads at that particular position.


3. Annotate SNPs using genome annotation data provided within a genbank file.
4. Extract the distribution of all indels according to genome annotation.
5. Visualize and explore the annotated SNPs for all analyzed files.
6. Combine all snpToolkit output files generated using the annotate option and generate:

   * A table storing the distribution of all SNPs on each sample
   * A fasta file with all concatenated SNPs for each sample. such file can be used to build a phylogenetic tree.

7. Analyse your data using two dimentionality reduction methods: PCA and UMAP.

snpToolkit detects automatically if the input vcf files were generated using samtools mpileup, gatk HaplotypeCaller or freebayes. Vcf files could be gzipped or not.
