.. image:: index/snpToolkit_LOGO.png
   :target: index/snpToolkit_LOGO.png
   :alt: index/snpToolkit_LOGO.png



**snpToolkit** is a computational framework written in Python 3. snpToolkit allows users to:


#. Visualize the content of their VCF files.
#. Filter SNPs based on multiple criteria:

   * Distance between SNPs
   * Coordinates of regions to exclude
   * Depth of coverage
   * Quality
   * The ratio corresponding to the number of reads that have the mutated allele / total number of reads at that particular position.


#. Annotate SNPs using genome annotation data provided within a genbank file.
#. Extract the distribution of all indels according to genome annotation.
#. Visualize and explore the annotated SNPs for all analyzed files.
#. Combine all snpToolkit output files generated using the annotate option and produce:

   * A table storing the distribution of all SNPs on each sample
   * A fasta file with all concatenated SNPs for each sample. such file can be used to build a phylogenetic tree.

#. Analyse your data using two dimentionality reduction methods: PCA and UMAP.

snpToolkit detects automatically if the input vcf files were generated using samtools mpileup, gatk HaplotypeCaller or freebayes. Vcf files can be in gzipped format or not.


********
Contents
********

.. toctree::
   :maxdepth: 2
   :numbered:
   
   1_How_to_Install
   2_snpToolkit_menu
   3_The_explore_command
   4_The_annotate_command
   5_The_viz_command
   6_The_combine_command
   7_The_analyse_command
