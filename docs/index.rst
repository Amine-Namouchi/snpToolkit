#############
snpToolkit
#############


.. image:: index/snpToolkit_LOGO.png
   :target: index/snpToolkit_LOGO.png
   :alt: index/snpToolkit_LOGO.png



**snpToolkit** is a computational framework written in Python 3. snpToolkit allow users to:


#. Visualize the content of their VCF files
#. Filter SNPs based on multiple criteria:

   * Coordinates of regions to exclude
   * Depth of coverage
   * Quality
   * The ratio corresponding to the number of reads that have the mutated allele / total number of reads at that particular position

#. Annotate SNPs using genome annotation data provided within a genbank file
#. Extract the distribution of all indels according to genome annotation
#. Visualize and explore the annotated SNPs for all analyzed files
#. Combine all snpToolkit output files generated using the annotate option and generate:

   * A table storing the distribution of all SNPs on each sample
   * A fasta file with all concatenated SNPs for each sample. such file can be used to build a phylogenetic tree.

#. Expand existent list of polymorphic sites when new SNPs output files are available

snpToolkit detects automatically if the input vcf files were generated using samtools mpileup, gatk HaplotypeCaller or freebayes. Vcf files could be gzipped or not.


********
Contents
********

.. toctree::
   :maxdepth: 2
   :numbered:
   
   1_ListTools
   2_ReadsFiltering
   3_Metagenomics_v2
   4_ReadsMapping_v2
   7_VariantsCall
   8_Filtering_SNPs
   9_DIY
   

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
