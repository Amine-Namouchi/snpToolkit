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

# 6e1be7
__licence__ = 'GPLv3'
__author__ = 'Amine Namouchi'
__author_email__ = 'amine.namouchi@gmail.com'
__version__ = '2.1'

import argparse
import logging
import coloredlogs


def get_options():

    __licence__ = 'GPLv3'
    __author__ = 'Amine Namouchi'
    __author_email__ = 'amine.namouchi@gmail.com'
    __version__ = '2.1'

    copyright = __licence__ + " licence | " + __author__ + " | " + __author_email__

    description = """
    snpToolkit takes vcf files, as well as bam files (optional) as inputs. The vcf files could be generated using samtools/bcftools, gatk HaplotypeCaller or freeBayes.
    Please visit https://github.com/Amine-Namouchi/snpToolkit for more information.
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter, epilog=copyright)

    subparsers = parser.add_subparsers(help='commands', dest='command')


    explore = subparsers.add_parser(
        name='explore', help='explore your vcf files before annotation', epilog=copyright)

    requiredOptions1 = explore.add_argument_group(
        'snpToolkit explore required options')

    requiredOptions1.add_argument('-i', required=True, type=str, dest='identifier',
                                  help='Provide the input vcf files')

    annotate = subparsers.add_parser(
        name='annotate', help='Please provide one or multiple vcf files', epilog=copyright)

    requiredOptions2 = annotate.add_argument_group(
        'snpToolkit annotate required options')

    requiredOptions2.add_argument('-i', required=True, type=str, dest='identifier',
                                  help='provide a specific identifier to recognize the file(s) to be analyzed')

    requiredOptions2.add_argument('-g', required=True, type=str, dest='genbank',
                                  help='Pleae provide a genbank file')

    AdditonalOptions2 = annotate.add_argument_group(
        'snpToolkit annotate additional options')
    
    AdditonalOptions2.add_argument('-p', required=False,  default=1, dest='processors', type=int,
                                   help='number of vcf files to be annotated in parallel default value [1]')


    AdditonalOptions2.add_argument('-f', required=False, dest='excludeCloseSNPs', type=int,
                                   help='exclude SNPs if the distance between them is lower then the specified window size in bp')

    AdditonalOptions2.add_argument('-q', required=False, default=20, type=int, dest='quality',
                                   help='quality score to consider as a cutoff for variant calling. default value [20]')

    AdditonalOptions2.add_argument('-d', required=False, default=3, type=int,
                                   dest='depth',  help='minimum depth caverage. default value [3]')

    AdditonalOptions2.add_argument('-r', required=False,
                                   dest='ratio', default=0.001, type=float, help='minimum ratio that correspond to the number of reads that has the mutated allele / total depth in that particular position. default value [0]')

    AdditonalOptions2.add_argument('-e', required=False, type=str, dest='exclude',
                                   help='provide a tab file with genomic regions to exclude in this format: region	start stop. region must correspond to the same name(s) of chromsome and plasmids as in the genbank file')

    combine = subparsers.add_parser(
        name='combine', help='combine snpToolkit output files in one alignment in fasta format', epilog=copyright)

    requiredOptions3 = combine.add_argument_group(
        'snpToolkit combine required options')

    requiredOptions3.add_argument('--location', required=True, type=str, dest='location',
                                  help='provide for example the name of the chromosome or plamid you want to create fasta alignemnt for')

    AdditonalOptions3 = combine.add_argument_group(
        'snpToolkit additional options')

    AdditonalOptions3.add_argument(
        '-r', required=False, type=float, dest='ratio', default=0.01, help='SNP ratio')

    AdditonalOptions3.add_argument('-d', required=False, type=int, dest='depth',
                                   default=0, help='depth cutoff for cheking missing data')

    AdditonalOptions3.add_argument('--bam', required=False, type=str,
                                   dest='bamFolder', help='path to the folder containing the bam files')

    AdditonalOptions3.add_argument('--snps', required=False, type=str, dest='snps', default='all', choices=['ns', 's', 'all', 'inter'],
                                   help='Specify if you want to concatenate all SNPs or just synonymous (s), non-synonymous (ns) or intergenic (inter) SNPs. default [all]')

    AdditonalOptions3.add_argument('-e', required=False, type=str, dest='exclude',
                                   help='provide a yaml file with keywords and coordinates to be excluded') 



    viz = subparsers.add_parser(
        name='viz', help='visualize snptoolkit output files', epilog=copyright)

    requiredOptions4 = viz.add_argument_group('snpToolkit viz required options')
    requiredOptions4.add_argument('--dir', required=True, type=str, dest='directory',
                                  help='provide the path to the directory containing snptoolkit output files')


    expand = subparsers.add_parser(
        name='expand', help='expand curent phylogeny', epilog=copyright)

    requiredOptions4 = expand.add_argument_group('snpToolkit viz required options')
    requiredOptions4.add_argument('--dir', required=True, type=str, dest='directory',
                                  help='provide the path to the directory containing snptoolkit output files to be added')
    requiredOptions4.add_argument('-p', required=True, type=str, dest='polymorphic_sites',
                                help='provide the polymorphic sites file already generated')
    requiredOptions4.add_argument('-o', required=True, type=str, dest='output',
                                help='output')
    requiredOptions4.add_argument('--bam', required=True, type=str, dest='bamfiles_directory',
                                help='bam files directory of already analysed anciant DNA')

    requiredOptions4.add_argument('-l', required=True, type=str, dest='location',
                                help='name of chromosome, plasmid, contig')

    requiredOptions4.add_argument('-c', required=True, type=str, dest='cutoff',
                                help='cutoff  of coverage if bam file is provided ')
    requiredOptions4.add_argument('-e', required=True, type=str, dest='exclude',
                                help='exclude yaml ')

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
