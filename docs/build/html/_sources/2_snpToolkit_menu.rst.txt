
snpToolkit menu
===============

.. code-block:: bash

   $ snptoolkit -h
   usage: snptoolkit [-h] {explore,annotate,combine,viz,expand} ...

       snpToolkit can takes vcf files, as well as bam files (optional) as inputs. The vcf files could be generated using samtools/bcftools, gatk HaplotypeCaller or freeBayes.
       Please visit https://snptoolkit.readthedocs.io/en/latest/index.html for more information.

   positional arguments:
     {explore,annotate,combine,viz,expand}
                           commands
       explore             Explore your vcf files before annotation
       annotate            Annotate one or multiple vcf files
       combine             Combine snpToolkit output files in one alignment in fasta format
       viz                 Visualize snptoolkit output files


   optional arguments:
     -h, --help            show this help message and exit
