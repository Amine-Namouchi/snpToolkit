#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name='snpToolkit',
      version='1.0',
      description='SNPs filtering, annotation and combining',
      url='https://github.com/Amine-Namouchi/snpToolkit',
      author='Amine Namouchi',
      author_email='amine.namouchi@gmail.com',
      long_description='SNPs annotation, filtering and combining from vcf files',
      license='GPLv3+',
      classifiers=[
          'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
          'Programming Language :: Python :: 3',
          'Operating System :: OS Independent',
      ],
      keywords='SNPs filtering annotation variants calling',
      install_requires=['pysam', 'pandas', 'tqdm', 'coloredlogs', 'Biopython'],
      )
