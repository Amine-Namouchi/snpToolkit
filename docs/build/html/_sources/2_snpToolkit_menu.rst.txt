
snpToolkit menu
===============

.. code-block:: bash

   $ snptoolkit -h
   usage: snptoolkit [-h] {explore,annotate,combine,viz,expand} ...

       snpToolkit takes vcf files, as well as bam files (optional) as inputs. The vcf files could be generated using samtools/bcftools, gatk HaplotypeCaller or freeBayes.
       Please visit https://github.com/Amine-Namouchi/snpToolkit for more information.

   positional arguments:
     {explore,annotate,combine,viz,expand}
                           commands
       explore             explore your vcf files before annotation
       annotate            Annotate one or multiple vcf files
       combine             combine snpToolkit output files in one alignment in fasta format
       viz                 visualize snptoolkit output files
       expand              expand existent list of polymorphic sites when new SNP output files are availble

   optional arguments:
     -h, --help            show this help message and exit
