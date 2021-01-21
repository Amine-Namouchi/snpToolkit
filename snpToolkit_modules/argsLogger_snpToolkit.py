#!/usr/bin/env python

# Copyright Amine Namouchi
# snpToolkit is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but without
# any warranty; without even the implied warranty of merchantability or fitness
# for a particular purpose. See the GNU General Public License for more details
#
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

__licence__ = 'GPLv3'
__author__ = 'Amine Namouchi'
__author_email__ = 'amine.namouchi@gmail.com'
__version__ = '2.2.8'

import argparse
import logging
import coloredlogs


def get_options():
    copyright = __licence__ + " licence | " + __author__ + " | " + __author_email__ + " | snpToolkit version: " + __version__

    description = """
    snpToolkit can takes vcf files, as well as bam files (optional) as inputs. The vcf files could be generated using samtools/bcftools, gatk HaplotypeCaller or freeBayes.
    Please visit https://snptoolkit.readthedocs.io/en/latest/index.html for more information.
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter, epilog=copyright)

    subparsers = parser.add_subparsers(help='commands', dest='command')

    explore = subparsers.add_parser(
        name='explore', help='Explore your vcf files before annotation', epilog=copyright)

    requiredOptions1 = explore.add_argument_group(
        'snpToolkit explore required options')

    requiredOptions1.add_argument('-i', required=True, type=str, dest='identifier',
                                  help='Provide the input vcf files')

    annotate = subparsers.add_parser(
        name='annotate', help='Annotate one or multiple vcf files', epilog=copyright)

    requiredOptions2 = annotate.add_argument_group(
        'snpToolkit annotate required options')

    requiredOptions2.add_argument('-i', required=True, type=str, dest='identifier',
                                  help='Provide a specific identifier to recognize the file(s) to be analyzed')

    requiredOptions2.add_argument('-g', required=True, type=str, dest='genbank',
                                  help='Provide a genbank file')

    AdditonalOptions2 = annotate.add_argument_group(
        'snpToolkit annotate additional options')
    
    AdditonalOptions2.add_argument('-p', required=False,  default=1, dest='processors', type=int,
                                   help='Number of vcf files to be annotated in parallel default value [1]')


    AdditonalOptions2.add_argument('-f', required=False, dest='excludeCloseSNPs', type=int,
                                   help='Exclude SNPs if the distance between them is lower then the specified window size in bp')

    AdditonalOptions2.add_argument('-q', required=False, default=20, type=int, dest='quality',
                                   help='Quality score to consider as a cutoff for variant calling. default value [20]')

    AdditonalOptions2.add_argument('-d', required=False, default=3, type=int,
                                   dest='depth',  help='Minimum depth caverage. default value [3]')

    AdditonalOptions2.add_argument('-r', required=False,
                                   dest='ratio', default=0.001, type=float, help='Minimum ratio that correspond to the number of reads that has the mutated allele / total depth in that particular position. default value [0]')

    AdditonalOptions2.add_argument('-e', required=False, type=str, dest='exclude',
                                   help='Please provide a tab file with genomic regions to exclude in this format: region	start stop. region must correspond to the same name(s) of chromsome or plasmids as indicated in the genbank file')

    combine = subparsers.add_parser(
        name='combine', help='Identify polymorphic sites and create distribution table and alignment file in fasta format', epilog=copyright)

    requiredOptions3 = combine.add_argument_group(
        'snpToolkit combine required options')

    requiredOptions3.add_argument('--loc', required=True, type=str, dest='location',
                                  help='Please provide for example the name of the chromosome or plasmid you want to create fasta alignemnt for')

    AdditonalOptions3 = combine.add_argument_group(
        'snpToolkit additional options')

    AdditonalOptions3.add_argument(
        '-r', required=False, type=float, dest='ratio', default=0.01, help='new versus reference allele ratio to filter SNPs from snpToolkit outputs. default [0]')

    AdditonalOptions3.add_argument('--bam', required=False, type=str,
                                   dest='bamFilter', nargs=3, help='provide the depth, ratio and the path to the folder containing the bam files. eg. 3 0.9 path')

    AdditonalOptions3.add_argument('--snps', required=False, type=str, dest='snps', default='all', choices=['ns', 's', 'all', 'inter'],
                                   help='Specify if you want to concatenate all SNPs or just synonymous (s), non-synonymous (ns) or intergenic (inter) SNPs. default [all]')

    AdditonalOptions3.add_argument('-e', required=False, type=str, dest='exclude',
                                   help='provide a yaml file with keywords and coordinates to be excluded') 


    viz = subparsers.add_parser(
        name='viz', help='visualize snptoolkit output files', epilog=copyright)

    requiredOptions4 = viz.add_argument_group('snpToolkit viz required options')
    requiredOptions4.add_argument('--dir', required=False, type=str, dest='directory',
                                  help='provide the path of the directory containing snptoolkit SNPs output files')


    analyze = subparsers.add_parser(
        name='analyse', help='analyse your SNPs data', epilog=copyright)
    requiredOptions5 = analyze.add_argument_group('snpToolkit analyze required options')

    requiredOptions5.add_argument('-p', required=True, type=str, dest='polymorphic_sites',
                                  help='provide the path of the polymorphic sites you want to analyze')
    requiredOptions5.add_argument('-c', required=False, type=str, dest='config',
                                  help='provide the path of the configuration file that contains the information to use for data visualization')


    return parser.parse_args()

def setupLogger():

    logger = logging.getLogger(__name__)
    coloredlogs.install(logger=logger, fmt='[%(asctime)s] [%(levelname)s] [%(message)s]', datefmt='%H:%M:%S', field_styles={
                        'levelname': {'color': 'white', 'bold': False, 'background': 'black'},  'asctime': {'color': 'green'}})
    filehandle = logging.FileHandler('snpToolkit.log')
    formatter = logging.Formatter(
        fmt='[%(asctime)s] [%(levelname)s] [%(message)s]', datefmt='%H:%M:%S')
    filehandle.setFormatter(formatter)
    logger.addHandler(filehandle)
    return logger
