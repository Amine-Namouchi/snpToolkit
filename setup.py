#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='snptoolkit',
      version='2.2.9',
      description='Explore, filter, annotate and analyse your SNPs/indels from VCF files',
      url='https://github.com/Amine-Namouchi/snpToolkit',
      author='Amine Namouchi',
      author_email='amine.namouchi@gmail.com',
      python_requires='>=3.0',
      long_description="""Please visit the full documentation at https://snptoolkit.readthedocs.io/en/latest/index.html
                        snpToolkit is a computational framework written in Python 3. snpToolkit allows users to:
                        - Visualize the content of their VCF files.
                        - Filter SNPs based on multiple criteria.
                        - Extract the distribution of all indels according to genome annotation.
                        - Visualize and explore the annotated SNPs for all analyzed files.
                        - Combine all snpToolkit output files generated using the annotate option.
                        - Analyse your data using two dimentionality reduction methods: PCA and UMAP.
                        snpToolkit detects automatically if the input vcf files were generated using samtools mpileup, gatk HaplotypeCaller or freebayes. Vcf files could be gzipped or not.""",
      license='GPLv3+',
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Programming Language :: Python :: 3',
          'Operating System :: OS Independent',
      ],
      keywords='SNPs filtering annotation variants calling',
      install_requires=['pysam', 'pandas', 'Biopython','dash','dash-daq','plotly','scikit-learn','umap-learn==0.4.6','llvmlite==0.34.0','numba==0.51.2','tqdm', 'coloredlogs'],
      scripts=['snptoolkit','./snpToolkit_modules/plot_polySites_output.py','./snpToolkit_modules/annotate_snpToolkit.py','./snpToolkit_modules/annotate_snpToolkit.py','./snpToolkit_modules/calls_snpToolkit.py', './snpToolkit_modules/argsLogger_snpToolkit.py','./snpToolkit_modules/combine_snpToolkit.py','./snpToolkit_modules/explore_snpToolkit.py','./snpToolkit_modules/plot_snptoolkit_output.py']
      )
