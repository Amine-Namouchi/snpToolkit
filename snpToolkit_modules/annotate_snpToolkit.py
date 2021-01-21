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

from .version import version 
__licence__ = 'GPLv3'
__author__ = 'Amine Namouchi'
__author_email__ = 'amine.namouchi@gmail.com'
__version__ = '2.2.8'


import gzip
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from itertools import islice
from pathlib import Path
import glob
import ntpath
import sys
import logging
import coloredlogs
from tqdm import tqdm


standcode = {
    'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C', 'ttc': 'F', 'tcc': 'S',
    'tac': 'Y', 'tgc': 'C', 'tta': 'L', 'tca': 'S', 'taa': '#', 'tga': '#',
    'ttg': 'L', 'tcg': 'S', 'tag': '#', 'tgg': 'W', 'ctt': 'L', 'cct': 'P',
    'cat': 'H', 'cgt': 'R', 'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
    'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R', 'ctg': 'L', 'ccg': 'P',
    'cag': 'Q', 'cgg': 'R', 'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
    'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S', 'ata': 'I', 'aca': 'T',
    'aaa': 'K', 'aga': 'R', 'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
    'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G', 'gtc': 'V', 'gcc': 'A',
    'gac': 'D', 'ggc': 'G', 'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
    'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
}

Reverse_nuc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def parse_genbank_file(gb_filename):
    """
    parse_genbank_file takes as input file a genbank file and retreive for each entry (chromosome, plasmids)
    all data relative to each gene, rRNA, tRNA, ncRNA as follow:
    gene|synonymous  type  start   stop  orientation   product  nucleotide_sequence   protein_sequence
    """

    gbank = SeqIO.parse(open(gb_filename, "r"), "genbank")

    locus_order = []
    collectionDB = {}
    for genome in gbank:
        locus_order.append(genome.id)
        whole_nuc_seq=genome.seq
        collect_data = []
        for feature in genome.features:
            if(feature.type == "CDS"):
                if feature.location_operator == 'join':
                    #TODO: here it is important to pass the coordinates of the fragsment to join and the the other intromic' part. this is important as it is necessray to change the coordiantes of the SNPs according to the new generated sequence. 
                    #need to work more with the example sample1_n50_gatk.vcf
                    nucseq=''
                    for sub_location in feature.location.parts:
                        if sub_location.strand == 1:
                            nucseq = nucseq + whole_nuc_seq[sub_location.start:sub_location.end]
                else:
                    if 'gene' in feature.qualifiers.keys() and 'pseudo' not in feature.qualifiers.keys():
                        try:
                            locus = feature.qualifiers['locus_tag'][0]
                        except KeyError:
                            try:
                                locus = feature.qualifiers['gene_synonym'][0]
                            except KeyError:
                                locus = feature.qualifiers['gene'][0]
                        gene = feature.qualifiers['gene'][0]
        
                        start = format(feature.location.start.position + 1)
                        stop = format(feature.location.end.position)
                        strand = feature.location.strand
                        product = feature.qualifiers['product'][0]
                        nucseq = feature.extract(genome.seq)
                        if strand == 1:
                            collect_data.append(
                                [locus + '|' + gene, 'CDS', start, stop, '+', product, format(nucseq)])
                        else:
                            collect_data.append(
                                [locus + '|' + gene, 'CDS', start, stop, '-', product, format(nucseq)])
                    else:
                        if 'pseudo' not in feature.qualifiers.keys():
                            try:
                                locus = feature.qualifiers['locus_tag'][0]
                            except KeyError:
                                try:
                                    locus = feature.qualifiers['gene_synonym'][0]
                                except KeyError:
                                    locus = feature.qualifiers['gene'][0]
                            gene = feature.qualifiers['locus_tag'][0]
                            start = format(feature.location.start.position + 1)
                            stop = format(feature.location.end.position)
                            strand = feature.location.strand
                            product = feature.qualifiers['product'][0]
                            nucseq = feature.extract(genome.seq)
                            if strand == 1:
                                collect_data.append(
                                    [locus + '|' + gene, 'CDS', start, stop, '+', product, format(nucseq)])
                            else:
                                collect_data.append(
                                    [locus + '|' + gene, 'CDS', start, stop, '-', product, format(nucseq)])
            else:

                if feature.type in ['tRNA', 'rRNA', 'ncRNA']:
                    try:
                        ID1 = feature.qualifiers['locus_tag'][0]
                        ID2 = feature.qualifiers['locus_tag'][0]
                    except KeyError:
                        try:
                            ID1 = feature.qualifiers['gene'][0]
                            ID2 = feature.qualifiers['gene'][0]
                        except KeyError:
                            ID1=feature.type
                            ID2=feature.type
                    type = feature.type
                    start = format(feature.location.start.position + 1)
                    stop = format(feature.location.end.position)
                    strand = feature.location.strand
                    if 'product' in feature.qualifiers.keys():
                        product = feature.qualifiers['product'][0]
                    else:
                        if 'note' in feature.qualifiers.keys():
                            product = feature.qualifiers['note'][0]
                    nucseq = feature.extract(genome.seq)
                    if strand == 1:
                        collect_data.append(
                            [ID1 + '|' + ID2, type, start, stop, '+', product, format(nucseq), '*'])
                    else:
                        collect_data.append(
                            [ID1 + '|' + ID2, type, start, stop, '-', product, format(nucseq), '*'])

        collectionDB[genome.id + ': ' + genome.description +
                     ' {} bp'.format(len(genome.seq))] = collect_data
    collectionDB['locus_order'] = locus_order
    if collectionDB == {}:
        logging.error('Please check your genbank file')
        sys.exit(0)
    return collectionDB


def isVCF(filename):

    if filename.endswith(".vcf.gz"):
        filehandle = gzip.open(filename, 'rb')
        line = filehandle.readline()
        if b'##fileformat=VCF' in line:
            return True
        else:
            return False

    elif filename.endswith(".vcf"):
        filehandle = open(filename, 'r')
        line = filehandle.readline()
        if '##fileformat=VCF' in line:
            return True
        else:
            return False


class VCFtoolBox(object):
    """
    This class will handle vcf files to:
        1. validates if a given file is a valid vcf file
        2. distinguish how the vcf file was generated using mpileup or gatk
        3. extract all single nucleotides polymorphism
        4. extract all insertion and deletions
    """

    def __init__(self, filename):
        self.file_name = filename
        if self.file_name.endswith('.vcf'):
            self.fh = open(self.file_name, 'r')
            self.file_type = 'string'
        else:
            if self.file_name.endswith('.vcf.gz'):
                self.fh = gzip.open(self.file_name, 'rb')
                self.file_type = 'bytes'
            else:
                self.fh = zlib.decompress.open(self.file_name, 'rb')
                self.file_type = 'bytes'

    def vcf_generator(self):
        if self.file_type == 'string':
            for line in islice(self.fh, 10):
                if 'mpileup' in line:
                    return 'samtools'
                elif 'GATK' in line:
                    return 'GATK'
                elif 'freeBayes' in line:
                    return 'freeBayes'
                else:
                    pass
            return '%s will be skipped as it seems to be not a valid vcf file or generated with a program different than samtools, gatk or freeBayes!' % self.file_name

        else:
            for line in islice(self.fh, 10):
                if b'mpileup' in line:
                    return 'samtools'
                elif b'GATK' in line:

                    return 'GATK'
                elif b'freeBayes' in line:
                    return 'freeBayes'
                else:
                    pass
            return '%s will be skipped as it seems to be not a valid vcf file or generated with a program different than samtools, gatk or freeBayes!' % self.file_name

    def extract_all_snps(self):
        dna = ('A', 'C', 'T', 'G')
        if self.file_type == 'string':
            list_snps = [line.strip().split('\t') for line in self.fh if '##' not in line and 'REF	ALT' not in line and line.split(
                '\t')[3] in dna and line.split('\t')[4] in dna and line.split('\t')[3] != line.split('\t')[4]]
        else:
            list_snps = [str(line, 'utf-8').strip().split('\t') for line in self.fh if '##' not in str(line, 'utf-8') and 'REF	ALT' not in str(line, 'utf-8') and str(
                line, 'utf-8').split('\t')[3] in dna and str(line, 'utf-8').split('\t')[4] in dna and str(line, 'utf-8').split('\t')[3] != str(line, 'utf-8').split('\t')[4]]

        return list_snps

    def extract_all_indels(self):
        if self.file_type == 'string':
            list_indels = [line.strip().split('\t') for line in self.fh if '##' not in line and 'REF	ALT' not in line and 'INDELS' in line and line.split(
                '\t')[3] != line.split('\t')[4] and line.split('\t')[6] == 'PASS']
        else:
            list_indels = [str(line, 'utf-8').strip().split('\t') for line in self.fh if '##' not in str(line, 'utf-8') and 'REF	ALT' not in str(line, 'utf-8') and 'INDEL' in str(line, 'utf-8') and
                           str(line, 'utf-8').split('\t')[3] != str(line, 'utf-8').split('\t')[4] and str(line, 'utf-8').split('\t')[6] == 'PASS']

        return list_indels

    def extract_all_variation(self):
        if self.file_type == 'string':
            list_indels = [line.strip().split('\t') for line in self.fh if '##' not in line and 'REF	ALT' not in line and line.split(
                '\t')[3] != line.split('\t')[4]]
        else:
            list_indels = [str(line, 'utf-8').strip().split('\t') for line in self.fh if '##' not in str(line, 'utf-8') and 'REF	ALT' not in str(line, 'utf-8') and
                           str(line, 'utf-8').split('\t')[3] != str(line, 'utf-8').split('\t')[4]]

        return list_indels


class SNPfiltering (object):
    """
    This class allow:
        1. check the presence of specific SNP among all extracted SNPs
        2. Exclude SNPs if the distance between them is lower than a defined distance
        3. Exclude SNPs if they are located in regions specified by the user
        4. Exclude SNPs if the distance between them is lower than a defined distance and if they are located in regions specified by the user
    """

    def __init__(self, snps):
        self.snps = snps

    def check_presence(self, location, coordinate):
        for eachElem in self.snps:
            if location in eachElem[0] and coordinate == int(eachElem[1]):
                reference = eachElem[3]
                variant = eachElem[4]
                quality = eachElem[5]
                return '%s: %s%s%s' % (location, reference, str(coordinate), variant)
        return 'No SNP found in %s at posSNPition %d' % (location, coordinate)

    def FilterCloseSNPs(self, distance):
        """
        SNP_coordinates is set by defaut to 1 by reference to vcf files. this method will be used in case of tab input file too. user just set the column where the snp coordinates are located.
        distance refer to the distance betwen to adjacent SNPs
        """
        filteredSNPs = []

        i = 0
        while i < len(self.snps) - 1:
            if i == 0:
                if int(self.snps[i + 1][1]) - int(self.snps[i][1]) > distance:
                    filteredSNPs.append(self.snps[i])
            elif i == len(self.snps) - 1:
                if int(self.snps[i + 1][1]) - int(self.snps[i][1]) > distance:
                    filteredSNPs.append(self.snps[i])
                    filteredSNPs.append(self.snps[i + 1])
            else:
                if int(self.snps[i][1]) - int(self.snps[i-1][1]) > distance and int(self.snps[i + 1][1]) - int(self.snps[i][1]) > distance:
                    filteredSNPs.append(self.snps[i])

            i += 1
        return filteredSNPs

    def FilterSNPtoExclude(self, coordinates_to_exclude):
        filteredSNPs = []
        for content in self.snps:
            flag = False
            for eachRegion in coordinates_to_exclude:
                if content[0] == eachRegion[0] and int(eachRegion[2]) >= int(content[1]) >= int(eachRegion[1]):
                    flag = True
                    break
            if flag == False:
                filteredSNPs.append(content)
        return filteredSNPs
        
    def FilterCloseSNPsandInREgions(self, distance, coordinates_to_exclude):
        filteredSNPs1 = []
        for content in self.snps:
            flag = False
            for eachRegion in coordinates_to_exclude:
                if content[0] == eachRegion[0] and int(eachRegion[2]) >= int(content[1]) >= int(eachRegion[1]):
                    flag = True
                    break
            if flag == False:
                filteredSNPs1.append(content)

        filteredSNPs2 = []
        i = 0
        while i < len(filteredSNPs1) - 1:
            if i == 0:
                if int(filteredSNPs1[i + 1][1]) - int(filteredSNPs1[i][1]) > distance:
                    filteredSNPs2.append(filteredSNPs1[i])
            elif i == len(filteredSNPs1) - 1:
                if int(filteredSNPs1[i + 1][1]) - int(filteredSNPs1[i][1]) > distance:
                    filteredSNPs2.append(filteredSNPs1[i])
                    filteredSNPs2.append(filteredSNPs1[i + 1])
            else:
                if int(filteredSNPs1[i + 1][1]) - int(filteredSNPs1[i][1]) > distance and int(filteredSNPs1[i][1]) - int(filteredSNPs1[i - 1][1]) > distance:
                    filteredSNPs2.append(filteredSNPs1[i])
            i += 1

        return filteredSNPs1, filteredSNPs2


class SNPselect (object):

    """
    This class allow to select SNPs based on their quality and depth. It can handle data from mpileup, gatk or freebayes
    """

    def __init__(self, snps, QUAL, DEPTH, RATIO):
        self.snps = snps
        self.QUAL = QUAL
        self.DEPTH = DEPTH
        self.RATIO = RATIO

    def ExtractSNPinfo(self, origin):
        count_QUAL = 0
        count_DEPTH = 0
        count_RATIO = 0
        data_collection = {}
        units = list(set([elem[0] for elem in self.snps]))
        for eachUnit in units:
            extracted_data = []
            for content in self.snps:
                if eachUnit == content[0]:
                    coordiantes = int(content[1])
                    REF = content[3]

                    ALT = content[4]
                    quality = float(content[5])
                    if quality >= self.QUAL:
                        count_QUAL += 1
                        varInfo = content[7].split(';')
                        if origin == 'samtools':
                            for eachElem in varInfo:
                                if 'DP=' in eachElem:
                                    depth = int(eachElem.split('=')[1])
                                elif 'DP4=' in eachElem:
                                    DP4 = eachElem.split('=')[1].split(',')
                                    nbREADref = int(DP4[0]) + int(DP4[1])
                                    nbREADalt = int(DP4[2]) + int(DP4[3])

                            if depth >= self.DEPTH:
                                count_DEPTH += 1
                                if nbREADref == 0:
                                    ratio = 1.0
                                    extracted_data.append([coordiantes, REF, ALT, depth, nbREADref, nbREADalt, float(
                                        ('%.2f') % ratio), quality])
                                    count_RATIO += 1
                                else:
                                    ratio = float(
                                        nbREADalt*1.0/(nbREADref+nbREADalt))
                                    if ratio >= self.RATIO:
                                        extracted_data.append([coordiantes, REF, ALT, depth, nbREADref, nbREADalt, float(
                                            ('%.2f') % ratio), quality])
                                        count_RATIO += 1
                        elif origin == 'GATK':
                            for eachElem in varInfo:
                                if 'DP=' in eachElem:
                                    depth = int(eachElem.split('=')[1])
                            snpInfo = content[9].split(':')[1].split(',')
                            nbREADref = int(snpInfo[0])
                            nbREADalt = int(snpInfo[1])
                           
                            if depth >= self.DEPTH:
 
                                count_DEPTH += 1
                                if nbREADref == 0:
                                    ratio = 1.0
                                    extracted_data.append([coordiantes, REF, ALT, depth, nbREADref, nbREADalt, float(
                                        ('%.2f') % ratio), quality])
                                    count_RATIO += 1
                                else:
                                    ratio = float(
                                        nbREADalt*1.0/(nbREADref+nbREADalt))
                                    if ratio >= self.RATIO:
                                        extracted_data.append([coordiantes, REF, ALT, depth, nbREADref, nbREADalt, float(
                                            ('%.2f') % ratio), quality])
                                        count_RATIO += 1

                        else:  # freeBayes
                            for eachElem in varInfo:
                                if 'DP=' in eachElem:
                                    depth = int(eachElem.split('=')[1])
                                elif 'SRF=' in eachElem:
                                    nbREADref_F = int(eachElem.split('=')[1])
                                elif 'SRR=' in eachElem:
                                    nbREADref_R = int(eachElem.split('=')[1])
                                elif 'SAF=' in eachElem:
                                    nbREADalt_F = int(eachElem.split('=')[1])
                                elif 'SAR=' in eachElem:
                                    nbREADalt_R = int(eachElem.split('=')[1])

                            nbREADref = nbREADref_F + nbREADref_R
                            nbREADalt = nbREADalt_F + nbREADalt_R

                            if depth >= self.DEPTH:
                                count_DEPTH += 1
                                if nbREADref == 0:
                                    ratio = 1.0
                                    extracted_data.append([coordiantes, REF, ALT, depth, nbREADref, nbREADalt, float(
                                        ('%.2f') % ratio), quality])
                                    count_RATIO += 1
                                else:
                                    ratio = float(
                                        nbREADalt*1.0/(nbREADref+nbREADalt))
                                    if ratio >= self.RATIO:
                                        extracted_data.append([coordiantes, REF, ALT, depth, nbREADref, nbREADalt, float(
                                            ('%.2f') % ratio), quality])
                                        count_RATIO += 1
            data_collection[eachUnit] = extracted_data
        return data_collection, count_QUAL, count_DEPTH, count_RATIO


class SNPannotation (object):
    """
    This class will handle SNPs lists Objects to map them according to a provided genbank file by the user
    """

    def __init__(self, filteredSNPs):
        self.filteredSNPs = filteredSNPs  # list of filtered SNPs

    def mapSNPs(self, GenbankContent):
        """
        this method will allow to map all provided SNPs according to their distribution in Gene,
        RBS, tRNA, rRNA, ncRNA and intergenic. In addition SNPs will be grouped according to their
        location in chromosome or plasmids. This information is retreived from the LOCUS tag in genbank file. The output of this method is a dictionary with the foloowing structure:
        {LOCUS:{'Genes':[snps......],'RBS':[snps.......], etc}}
        """
        data_collector = {}
        gbank_accessions = GenbankContent.keys()
        for eachKEY in self.filteredSNPs.keys():
            if len(self.filteredSNPs[eachKEY]) > 0:
                collected_snp = self.filteredSNPs[eachKEY]
                snps_in_genes = []
                snps_in_RBS = []
                snps_in_tRNA = []
                snps_in_rRNA = []
                snps_in_ncRNA = []
                snps_in_pseudogene = []
                intergenic = []
                for eachAccession in gbank_accessions:
                    if eachAccession.split(':')[0] in eachKEY:
                        annotationDB = GenbankContent[eachAccession]
                        for eachSNP in collected_snp:
                            flag = False
                            posSNP = eachSNP[0]
                            for eachElem in annotationDB:
                                if posSNP <= int(eachElem[3]) and posSNP >= int(eachElem[2]):
                                    flag = True
                                    if eachElem[1] == 'Gene':
                                        snps_in_genes.append(eachSNP)
                                    elif eachElem[1] == 'RBS':
                                        snps_in_RBS.append(eachSNP)
                                    elif eachElem[1] == 'tRNA':
                                        snps_in_tRNA.append(eachSNP)
                                    elif eachElem[1] == 'rRNA':
                                        snps_in_rRNA.append(eachSNP)
                                    elif eachElem[1] == 'ncRNA':
                                        snps_in_ncRNA.append(eachSNP)
                                    elif eachElem[1] == 'pseudogene':
                                        snps_in_pseudogene.append(eachSNP)
                            if flag == False:
                                intergenic.append(eachSNP)
                        data_collector[eachAccession] = {'Gene': snps_in_genes, 'tRNA': snps_in_tRNA, 'RBS': snps_in_RBS,
                                                         'rRNA': snps_in_rRNA, 'ncRNA': snps_in_ncRNA, 'pseudogene': snps_in_pseudogene, 'intergenic': intergenic}

        return data_collector

    def mapAndannoINDELS(self, GenbankContent):
        data_collector = {}
        gbank_accessions = GenbankContent.keys()
        for eachKEY in self.filteredSNPs.keys():
            if len(self.filteredSNPs[eachKEY]) > 0:
                mappedINDELS = []

                for eachAccession in gbank_accessions:
                    if eachAccession.split(':')[0] in eachKEY or eachKEY in eachAccession:
                        annotationDB = GenbankContent[eachAccession]
                        for eachINDEL in self.filteredSNPs[eachKEY]:
                            flag = False
                            posINDEL = eachINDEL[0]
                            reference_allele = eachINDEL[1]
                            mutated_allele = eachINDEL[2]
                            total_depth = eachINDEL[3]
                            nb_reads_reference = eachINDEL[4]
                            nb_reads_mutated = eachINDEL[5]
                            ratio = eachINDEL[6]
                            quality = eachINDEL[7]
                            if len(reference_allele) > len(mutated_allele):
                                event = 'deletion'
                                size_event = len(
                                    reference_allele) - len(mutated_allele)
                            else:
                                event = 'insertion'
                                size_event = len(
                                    mutated_allele)-len(reference_allele)
                            for eachElem in annotationDB:
                                if posINDEL <= int(eachElem[3]) and posINDEL >= int(eachElem[2]):
                                    flag = True

                                    type = eachElem[1]
                                    gene_name = eachElem[0]
                                    gene_function = eachElem[5]

                                    collected_data = eachINDEL + [type, gene_name, gene_function,
                                                                  event, size_event]
                                    # collected_data = [posINDEL, reference_allele, mutated_allele, event, size_event, total_depth,
                                    #     nb_reads_reference, nb_reads_mutated, ratio, quality, type, gene_name, gene_function]
                                    mappedINDELS.append(collected_data)

                            if flag == False:
                                mappedINDELS.append(
                                    eachINDEL + ['intergenic', '.', '.', event, size_event])
                                cutoff_distance = 100
                                if posINDEL < int(annotationDB[0][2]):
                                    if annotationDB[0][4] == '+':
                                        distance_gene = int(
                                            annotationDB[0][2]) - posINDEL
                                        if distance_gene <= cutoff_distance:
                                            mappedINDELS.append(eachINDEL + ['intergenic', '.', str(
                                                distance_gene) + ' bp from ' + annotationDB[0][0], event, size_event])
                                        else:
                                            mappedINDELS.append(
                                                eachINDEL + ['intergenic', '.', '.', event, size_event])

                                    else:
                                        mappedINDELS.append(
                                            eachINDEL +
                                            ['intergenic', '.', '.', event, size_event])
                                else:
                                    if posINDEL > int(annotationDB[-1][3]):
                                        if annotationDB[-1][4] == '-':
                                            distance_gene = posINDEL - \
                                                int(annotationDB[0][3])
                                            if distance_gene <= cutoff_distance:
                                                mappedINDELS.append(eachINDEL + ['intergenic', '.', '+', str(
                                                    distance_gene) + ' bp from ' + annotationDB[-1][0], '.'])
                                            else:
                                                mappedINDELS.append(
                                                    eachINDEL +
                                                    ['intergenic', '.',
                                                     '+', '.', '.'])

                                        else:
                                            mappedINDELS.append(
                                                eachINDEL +
                                                ['intergenic', '.',
                                                 '+', '.', '.'])

                                k = 0
                                while k < len(annotationDB) - 1:
                                    if int(annotationDB[k][3]) < posINDEL < int(annotationDB[k + 1][2]):
                                        if annotationDB[k][4] == '+' and annotationDB[k + 1][4] == '+':
                                            distance_gene = int(
                                                annotationDB[k + 1][2]) - posINDEL
                                            if distance_gene <= cutoff_distance:
                                                mappedINDELS.append(eachINDEL + ['intergenic', '.', '+', str(
                                                    distance_gene) + ' bp from ' + annotationDB[k + 1][0], '.'])
                                            else:
                                                mappedINDELS.append(
                                                    eachINDEL +
                                                    ['intergenic', '.',
                                                     '+', '.', '.'])

                                        elif annotationDB[k][4] == '-' and annotationDB[k + 1][4] == '+':
                                            distance_gene1 = posINDEL - \
                                                int(annotationDB[k][3])
                                            distance_gene2 = int(
                                                annotationDB[k + 1][2]) - posINDEL
                                            if distance_gene1 <= cutoff_distance and distance_gene2 <= cutoff_distance:
                                                mappedINDELS.append(eachINDEL + ['intergenic', '.', '+', str(distance_gene1) + ' bp from ' + annotationDB[k][0] + ' and ' + str(
                                                    distance_gene2) + ' bp from ' + annotationDB[k + 1][0], '.'])
                                            elif distance_gene1 <= cutoff_distance and distance_gene2 > cutoff_distance:
                                                mappedINDELS.append(eachINDEL + ['intergenic', '.', '+', str(
                                                    distance_gene1) + ' bp from ' + annotationDB[k][0], '.'])
                                            elif distance_gene2 <= cutoff_distance and distance_gene1 > cutoff_distance:
                                                mappedINDELS.append(eachINDEL + ['intergenic', '.', '+', str(
                                                    distance_gene2) + ' bp from ' + annotationDB[k + 1][0], '.'])
                                            else:
                                                mappedINDELS.append(
                                                    eachINDEL +
                                                    ['intergenic', '.',
                                                     '+', '.', '.'])

                                        elif annotationDB[k][4] == '-' and annotationDB[k + 1][4] == '-':
                                            distance_gene = posINDEL - \
                                                int(annotationDB[k][3])
                                            if distance_gene <= cutoff_distance:
                                                mappedINDELS.append(eachINDEL + ['intergenic', '.', '+', str(
                                                    distance_gene) + ' bp from ' + annotationDB[k][0], '.'])
                                            else:
                                                mappedINDELS.append(
                                                    eachINDEL +
                                                    ['intergenic', '.',
                                                     '+', '.', '.'])

                                        else:

                                            mappedINDELS.append(
                                                eachINDEL +
                                                ['intergenic', '.',
                                                 '+', '.', '.'])

                                    k += 1

                        data_collector[eachAccession] = mappedINDELS

        return data_collector

    def mapAndannoSNPs(self, GenbankContent):
        data_collector = {}
        gbank_accessions = GenbankContent.keys()
        for eachKEY in self.filteredSNPs.keys():
            if len(self.filteredSNPs[eachKEY]) > 0:
                snps_in_genes = []
                snps_in_RBS = []
                snps_in_tRNA = []
                snps_in_rRNA = []
                snps_in_ncRNA = []
                snps_in_pseudogene = []
                intergenic = []
                for eachAccession in gbank_accessions:
                    if eachAccession.split(':')[0] in eachKEY or eachKEY in eachAccession:

                        annotationDB = GenbankContent[eachAccession]
                        for eachSNP in self.filteredSNPs[eachKEY]:
             
                            flag = False
                            posSNP = eachSNP[0]
                            total_depth = eachSNP[3]
                            nb_reads_reference = eachSNP[4]
                            nb_reads_mutated = eachSNP[5]
                            ratio = eachSNP[6]
                            quality = eachSNP[7]

                            for eachElem in annotationDB:
                                if posSNP <= int(eachElem[3]) and posSNP >= int(eachElem[2]):
                                    flag = True
                                    if eachElem[1] == 'CDS':
                                        codons_nuc_sequence = list(zip(range(1, len(eachElem[6]), 3), [
                                            eachElem[6][i:i+3] for i in range(0, len(eachElem[6]), 3)]))
                                        gene_name = eachElem[0]
                                        gene_fucntion = eachElem[5]
                                        gene_orientation = eachElem[4]
                                        if gene_orientation == '+':
                                            pos_snp_in_gene = (
                                                eachSNP[0] - int(eachElem[2])) + 1
                                            reference_allele = eachSNP[1]
                                            mutated_allele = eachSNP[2]
                                        #TODO:even if the gene is in the reverse strand i decided to keep the SNPs as according to genome. in the other hand, it will be better to change the codon as follorw to show the position of the SNP C(G)T or (C)GT or CG(T)
                                        else:
                                            pos_snp_in_gene = (
                                                int(eachElem[3]) - eachSNP[0]) + 1
                                            reference_allele = eachSNP[1]
                                            mutated_allele = eachSNP[2]

                                        for i in range(len(codons_nuc_sequence)):
                                            if codons_nuc_sequence[i][0] <= pos_snp_in_gene < codons_nuc_sequence[i][0] + 3:
                                                ref_codon = codons_nuc_sequence[i][1]
                                                position_nuc_ToBe_changed = (
                                                    pos_snp_in_gene - codons_nuc_sequence[i][0])
                                                if gene_orientation == '+':
                                                    new_codon = ref_codon[:position_nuc_ToBe_changed] + \
                                                        '['+mutated_allele + ']'+\
                                                        ref_codon[position_nuc_ToBe_changed+1:]
                                                    new_codon_x = ref_codon[:position_nuc_ToBe_changed] + \
                                                        mutated_allele +\
                                                        ref_codon[position_nuc_ToBe_changed+1:]
                                                    ref_AA = standcode[ref_codon.lower(
                                                    )]
                                                    new_AA = standcode[new_codon_x.lower(
                                                    )]
                                                else:
                                                    new_codon = ref_codon[:position_nuc_ToBe_changed] + \
                                                        '['+Reverse_nuc[mutated_allele] + ']'+\
                                                        ref_codon[position_nuc_ToBe_changed+1:]
                                                    ref_AA = standcode[ref_codon.lower(
                                                    )]
                                                    new_codon_x = ref_codon[:position_nuc_ToBe_changed] + \
                                                        Reverse_nuc[mutated_allele] + \
                                                        ref_codon[position_nuc_ToBe_changed+1:]                                                    
                                                    new_AA = standcode[new_codon_x.lower(
                                                    )]
                                                collected_data = [posSNP, reference_allele, mutated_allele, total_depth, nb_reads_reference, nb_reads_mutated, ratio,
                                                                  quality, gene_name, gene_fucntion, gene_orientation, pos_snp_in_gene, ref_codon, new_codon, ref_AA, new_AA, i+1]
                                                if ref_AA == new_AA:
                                                    collected_data.append('Syn')

                                                else:
                                                    collected_data.append('NS')
                                                snps_in_genes.append(
                                                    collected_data)

                                    elif eachElem[1] == 'RBS':
                                        snps_in_RBS.append(
                                            eachSNP + ['RBS'] + ['.'] * 9)

                                    elif eachElem[1] == 'tRNA':
                                        snps_in_tRNA.append(
                                            eachSNP + ['tRNA'] + ['.'] * 9)

                                    elif eachElem[1] == 'rRNA':
                                        snps_in_rRNA.append(
                                            eachSNP + ['rRNA'] + ['.'] * 9)

                                    elif eachElem[1] == 'ncRNA':
                                        snps_in_ncRNA.append(
                                            eachSNP + ['ncRNA'] + ['.'] * 9)

                                    elif eachElem[1] == 'pseudogene':
                                        codons_nuc_sequence = list(zip(range(1, len(eachElem[6]), 3), [
                                            eachElem[6][i:i+3] for i in range(0, len(eachElem[6]), 3)]))
                                        gene_name = eachElem[0]+'-pseudogene'
                                        gene_fucntion = eachElem[5]
                                        gene_orientation = eachElem[4]
                                        if gene_orientation == '+':
                                            pos_snp_in_gene = (
                                                eachSNP[0] - int(eachElem[2])) + 1
                                            reference_allele = eachSNP[1]
                                            mutated_allele = eachSNP[2]
                                        #TODO:even if the gene is in the reverse strand i decided to keep the SNPs as according to genome. in the other hand, it will be better to change the codon as follorw to show the position of the SNP C(G)T or (C)GT or CG(T)
                                        else:
                                            pos_snp_in_gene = (
                                                int(eachElem[3]) - eachSNP[0]) + 1
                                            reference_allele = eachSNP[1]
                                            mutated_allele = eachSNP[2]

                                        for i in range(len(codons_nuc_sequence)):
                                            if codons_nuc_sequence[i][0] <= pos_snp_in_gene < codons_nuc_sequence[i][0] + 3:
                                                ref_codon = codons_nuc_sequence[i][1]
                                                position_nuc_ToBe_changed = (
                                                    pos_snp_in_gene - codons_nuc_sequence[i][0])
                                                if gene_orientation == '+':
                                                    new_codon = ref_codon[:position_nuc_ToBe_changed] + \
                                                        '['+mutated_allele + ']'+\
                                                        ref_codon[position_nuc_ToBe_changed+1:]
                                                    new_codon_x = ref_codon[:position_nuc_ToBe_changed] + \
                                                        mutated_allele +\
                                                        ref_codon[position_nuc_ToBe_changed+1:]
                                                    ref_AA = standcode[ref_codon.lower(
                                                    )]
                                                    new_AA = standcode[new_codon_x.lower(
                                                    )]
                                                else:
                                                    new_codon = ref_codon[:position_nuc_ToBe_changed] + \
                                                        '['+Reverse_nuc[mutated_allele] + ']'+\
                                                        ref_codon[position_nuc_ToBe_changed+1:]
                                                    ref_AA = standcode[ref_codon.lower(
                                                    )]
                                                    new_codon_x = ref_codon[:position_nuc_ToBe_changed] + \
                                                        Reverse_nuc[mutated_allele] + \
                                                        ref_codon[position_nuc_ToBe_changed+1:]                                                    
                                                    new_AA = standcode[new_codon_x.lower(
                                                    )]
                                                collected_data = [posSNP, reference_allele, mutated_allele, total_depth, nb_reads_reference, nb_reads_mutated, ratio,
                                                                  quality, gene_name, gene_fucntion, gene_orientation, pos_snp_in_gene, ref_codon, new_codon, ref_AA, new_AA, i+1]
                                                if ref_AA == new_AA:
                                                    collected_data.append('Syn')

                                                else:
                                                    collected_data.append('NS')
                                                snps_in_genes.append(
                                                    collected_data)

                                        snps_in_pseudogene.append(eachSNP)
                            if flag == False:

                                cutoff_distance = 100
                                if posSNP < int(annotationDB[0][2]):
                                    if annotationDB[0][4] == '+':
                                        distance_gene = int(
                                            annotationDB[0][2]) - posSNP
                                        if distance_gene <= cutoff_distance:
                                            intergenic.append(eachSNP + ['intergenic', '.', '+', str(
                                                distance_gene) + ' bp from ' + annotationDB[0][0]] + ['.'] * 6)
                                        else:
                                            intergenic.append(
                                                eachSNP + ['intergenic', '.', '+', '.'] + ['-'] * 6)

                                    else:
                                        intergenic.append(
                                            eachSNP + ['intergenic', '.', '+', '.'] + ['-'] * 6)
                                else:
                                    if posSNP > int(annotationDB[-1][3]):
                                        if annotationDB[-1][4] == '-':
                                            distance_gene = posSNP -\
                                                int(annotationDB[0][3])
                                            if distance_gene <= cutoff_distance:
                                                intergenic.append(eachSNP + ['intergenic', '.', '+', str(
                                                    distance_gene) + ' bp from ' + annotationDB[-1][0]] + ['.'] * 6)
                                            else:
                                                intergenic.append(
                                                    eachSNP + ['intergenic', '.', '+', '.'] + ['-'] * 6)

                                        else:
                                            intergenic.append(
                                                eachSNP + ['intergenic', '.', '+', '.'] + ['-'] * 6)

                                k = 0
                                while k < len(annotationDB) - 1:
                                    if int(annotationDB[k][3]) < posSNP < int(annotationDB[k + 1][2]):
                                        if annotationDB[k][4] == '+' and annotationDB[k + 1][4] == '+':
                                            distance_gene = int(
                                                annotationDB[k + 1][2]) - posSNP
                                            if distance_gene <= cutoff_distance:
                                                intergenic.append(eachSNP + ['intergenic', '.', '+', str(
                                                    distance_gene) + ' bp from ' + annotationDB[k + 1][0]] + ['.'] * 6)
                                            else:
                                                intergenic.append(
                                                    eachSNP + ['intergenic', '.', '+', '.'] + ['.'] * 6)

                                        elif annotationDB[k][4] == '-' and annotationDB[k + 1][4] == '+':
                                            distance_gene1 = posSNP - \
                                                int(annotationDB[k][3])
                                            distance_gene2 = int(
                                                annotationDB[k + 1][2]) - posSNP
                                            if distance_gene1 <= cutoff_distance and distance_gene2 <= cutoff_distance:
                                                intergenic.append(eachSNP + ['intergenic', '.', '+', str(distance_gene1) + ' bp from ' + annotationDB[k][0] + ' and ' + str(
                                                    distance_gene2) + ' bp from ' + annotationDB[k + 1][0]] + ['.'] * 6)
                                            elif distance_gene1 <= cutoff_distance and distance_gene2 > cutoff_distance:
                                                intergenic.append(eachSNP + ['intergenic', '.', '+', str(
                                                    distance_gene1) + ' bp from ' + annotationDB[k][0]] + ['.'] * 6)
                                            elif distance_gene2 <= cutoff_distance and distance_gene1 > cutoff_distance:
                                                intergenic.append(eachSNP + ['intergenic', '.', '+', str(
                                                    distance_gene2) + ' bp from ' + annotationDB[k + 1][0]] + ['.'] * 6)
                                            else:
                                                intergenic.append(
                                                    eachSNP + ['intergenic', '.', '+', '.'] + ['-'] * 6)

                                        elif annotationDB[k][4] == '-' and annotationDB[k + 1][4] == '-':
                                            distance_gene = posSNP - \
                                                int(annotationDB[k][3])
                                            if distance_gene <= cutoff_distance:
                                                intergenic.append(eachSNP + ['intergenic', '.', '+', str(
                                                    distance_gene) + ' bp from ' + annotationDB[k][0]] + ['.'] * 6)
                                            else:
                                                intergenic.append(
                                                    eachSNP + ['intergenic', '.', '+', '.'] + ['.'] * 6)

                                        else:

                                            intergenic.append(
                                                eachSNP + ['intergenic', '.', '+', '.'] + ['.'] * 6)

                                    k += 1

                        data_collector[eachAccession] = {'Genes': snps_in_genes, 'tRNA': snps_in_tRNA, 'RBS': snps_in_RBS,
                                                         'rRNA': snps_in_rRNA, 'ncRNA': snps_in_ncRNA, 'Pseudogenes': snps_in_pseudogene, 'intergenic': intergenic}

        return data_collector
