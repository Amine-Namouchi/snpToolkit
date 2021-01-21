
snpToolkit menu
===============

.. code-block:: bash

   $ snptoolkit -h
    usage: snptoolkit [-h] {explore,annotate,combine,viz,analyse} ...

        snpToolkit can takes vcf files, as well as bam files (optional) as inputs. The vcf files could be generated using samtools/bcftools, gatk HaplotypeCaller or freeBayes.
        Please visit https://snptoolkit.readthedocs.io/en/latest/index.html for more information.


    positional arguments:
      {explore,annotate,combine,viz,analyse}
                            commands
        explore             Explore your vcf files before annotation
        annotate            Annotate one or multiple vcf files
        combine             Identify polymorphic sites and create distribution table and alignment file in fasta format
        viz                 visualize snptoolkit output files
        analyse             analyse your SNPs data

    optional arguments:
      -h, --help            show this help message and exit
