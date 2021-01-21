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

import sys
from pathlib import Path
import pysam
from tqdm import tqdm




class CheckCoverage (object):

    def __init__(self, bamfile, location):
        self.bamfile = pysam.AlignmentFile(bamfile)
        self.location = location

    def location_exist(self):
        Listlocation = self.bamfile.references
        if self.location in Listlocation:
            return True
        else:
            return False

    def coverage(self, position):
        for pileup in self.bamfile.pileup(self.location, position, position + 1, truncate=True):
            return pileup.n


def checkBamindex(PathBamFiles):
    noindex = []
    FolderContent = Path(PathBamFiles)
    bamfiles = FolderContent.glob('*.bam')
    for eachFile in bamfiles:
        bam = pysam.AlignmentFile(str(eachFile))
        try:
            bam.check_index()
        except ValueError:
            noindex.append(str(eachFile))
    return [str(path) for path in noindex]


def checksnpToolkitBam(PathBamFiles):
    snpToolkitFolderContent = Path.cwd()
    bamFilterContent = Path(PathBamFiles)
    nameTObam = {}
    for name in snpToolkitFolderContent.glob('*_snpToolkit_SNPs.txt'):
        for bam in list(bamFilterContent.glob('*.bam')):
            if bam.stem in name.stem:
                nameTObam[name.name] = str(bam.parent) + '/' + bam.name
    return nameTObam


def excludeCloseSNPs(list_snps, regions_to_exclude):
    filteredSNPs = []
    Coordinates_to_exclude = []
    Keywords_to_exclude = []
    if 'COORDINATES' in regions_to_exclude.keys():
        Coordinates_to_exclude = [
            int(x) for x in regions_to_exclude['COORDINATES'].split(';')]
    if 'KEYWORDS' in regions_to_exclude.keys():
        Keywords_to_exclude = regions_to_exclude['KEYWORDS'].split(';')
    for snp in list_snps:
        if snp[0] not in Coordinates_to_exclude:
            if '|' in snp[3]:
                if snp[3].split('|')[0] not in Keywords_to_exclude and snp[3].split('|')[1] not in Keywords_to_exclude:
                    filteredSNPs.append(snp)
            else:
                if snp[3] not in Keywords_to_exclude:
                    filteredSNPs.append(snp)
    return filteredSNPs


def snp_distribution(list_snps, PathFiles, regions_to_exclude):

    if len(regions_to_exclude) > 0:
        filteredSNPs = excludeCloseSNPs(list_snps, regions_to_exclude)
    else:
        filteredSNPs = list_snps
    extract_files_content = {}
    for eachFile in PathFiles:
        with open(eachFile, 'r') as f:
            file_content = [[int(l.strip().split('\t')[0])] + l.strip().split(
                '\t')[1:3] + l.strip().split('\t')[8:-1]  for l in f if '##' not in l]
        extract_files_content[eachFile] = file_content

    for i in tqdm(range(len(filteredSNPs)), ascii=True, desc='progress'):
        eachSNP = filteredSNPs[i]
        searchSNP = eachSNP[:]
        for eachFile in extract_files_content.keys():

            if searchSNP in extract_files_content[eachFile]:
                eachSNP.append('1')
            else:
                eachSNP.append('0')

    return filteredSNPs


def snp_distribution_missing(location, ratio, snpType, list_snps, snpToolkitFiles, bamFilter, checkedBam, cutoff, checkRatio, regions_to_exclude):

    if len(regions_to_exclude) > 0:
        filteredSNPs = excludeCloseSNPs(list_snps, regions_to_exclude)
    else:
        filteredSNPs = list_snps

    Reverse_nuc = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    communFile = set(checkedBam.keys()).intersection(snpToolkitFiles)
    if len(communFile) > 0:
        FolderContent = Path(bamFilter)
        bamfiles = list(FolderContent.glob('*.bam'))
        for i in tqdm(range(len(filteredSNPs)), ascii=True, desc='progress'):
            eachSNP = filteredSNPs[i]
            searchSNP = eachSNP[:]
            if eachSNP[5] == '+':
                REF = eachSNP[1]
                SNP = eachSNP[2]
            else:
                REF = Reverse_nuc[eachSNP[1]]
                SNP = Reverse_nuc[eachSNP[2]]
            for eachFile in snpToolkitFiles:
                with open(eachFile, 'r') as f:
                    file_content = [[int(l.strip().split('\t')[0])] + l.strip().split(
                        '\t')[1:3] + l.strip().split('\t')[8:-1] for l in f if '##' not in l]
                if searchSNP in file_content:
                    eachSNP.append('1')
                else:
                    if eachFile in list(checkedBam.keys()):
                        bamFile = pysam.AlignmentFile(checkedBam[eachFile])
                        try:
                            depthATposition = 0
                            # bamFile.count_coverage() return four array.arrays of the same length in order A C G T
                            NucleotidesOrder = 'ACGT'
                            NumberOfReads = []
                            k = 0

                            for eachCov in bamFile.count_coverage(location, eachSNP[0] - 1, eachSNP[0]):
                                depthATposition = depthATposition + \
                                    list(eachCov)[0]
                                NumberOfReads.append(
                                    (list(eachCov)[0], NucleotidesOrder[k]))
                                k += 1
                            NumberOfReads.sort(reverse=True)
                            if depthATposition < cutoff:
                                eachSNP.append('?')
                            else:
                                if NumberOfReads[0][1] == SNP:
                                    if NumberOfReads[0][0] / depthATposition >= checkRatio:
                                        eachSNP.append('1')
                                    elif NumberOfReads[0][0] / depthATposition < checkRatio and NumberOfReads[0][0] / depthATposition >= 0.4 :
                                        eachSNP.append('?')
                                    else:
                                        eachSNP.append('0')
                                else:
                                    eachSNP.append('0')

                        except ValueError as e:
                            logging.error(location + ' was not detected in file ' +
                                          checkedBam[eachFile] + '. Please use samtools view -H on this bam file to check the exact name.')
                            sys.exit(0)
                    else:
                        eachSNP.append('0')
    else:
        list_snps_only_bamFiles = []
        for eachFile in list(communFile):
            with open(eachFile, 'r') as f:
                for l in f.readlines():
                    if l[:2] != '##':
                        content = l.strip().split('\t')
                        if float(content[6]) >= float(ratio):
                            snp_coordinate = int(content[0])
                            if snp_coordinate not in list_snps_only_bamFiles:
                                if location in l.split('\t')[-1] or location == l.split('\t')[-1]:
                                    if snpType == 'all':
                                        list_snps_only_bamFiles.append(
                                            snp_coordinate)
                                    elif snpType == 's':
                                        if l.split('\t')[17] == 'Syn':
                                            list_snps_only_bamFiles.append(
                                                snp_coordinate)
                                    elif snpType == 'ns':
                                        if l.split('\t')[17] == 'NS':
                                            list_snps_only_bamFiles.append(
                                                snp_coordinate)
                                    else:
                                        if l.split('\t')[8] == 'intergenic':
                                            list_snps_only_bamFiles.append(
                                                snp_coordinate)

        FolderContent = Path(bamFilter)
        bamfiles = list(FolderContent.glob('*.bam'))
        for i in tqdm(range(len(filteredSNPs)), ascii=True, desc='progress'):
            Flag = False
            if filteredSNPs[i][0] in list_snps_only_bamFiles:
                Flag = True
            eachSNP = filteredSNPs[i]
            searchSNP = eachSNP[:]
            if eachSNP[5] == '+':
                REF = eachSNP[1]
                SNP = eachSNP[2]
            else:
                REF = Reverse_nuc[eachSNP[1]]
                SNP = Reverse_nuc[eachSNP[2]]
            for eachFile in snpToolkitFiles:
                with open(eachFile, 'r') as f:
                    file_content = [[int(l.strip().split('\t')[0])] + l.strip().split(
                        '\t')[1:3] + l.strip().split('\t')[8:-1] for l in f if '##' not in l]
                if searchSNP in file_content:
                    eachSNP.append('1')
                else:
                    if Flag is True:
                        if eachFile in list(checkedBam.keys()):
                            bamFile = pysam.AlignmentFile(checkedBam[eachFile])
                            try:
                                depthATposition = 0
                                # bamFile.count_coverage() return four array.arrays of the same length in order A C G T
                                NucleotidesOrder = 'ACGT'
                                NumberOfReads = []
                                k = 0
                                # minimum quality of count_coverage is 15
                                for eachCov in bamFile.count_coverage(location, eachSNP[0] - 1, eachSNP[0]):
                                    depthATposition = depthATposition + \
                                        list(eachCov)[0]
                                    NumberOfReads.append(
                                        (list(eachCov)[0], NucleotidesOrder[k]))
                                    k += 1
                                NumberOfReads.sort(reverse=True)
                                if depthATposition <= cutoff:
                                    eachSNP.append('?')
                                else:
                                    if NumberOfReads[0][1] == SNP:
                                        if NumberOfReads[0][0] / depthATposition >= 0.8:
                                            eachSNP.append('1')
                                        elif 0.8 > NumberOfReads[0][0] / depthATposition >= 0.4:
                                            eachSNP.append('?')
                                        else:
                                            eachSNP.append('0')
                                    else:
                                        eachSNP.append('0')

                            except ValueError as e:
                                logging.error(location + ' was not detected in file ' +
                                              checkedBam[eachFile] + '. Please use samtools view -H on this bam file to check the exact name.')
                                sys.exit(0)
                        else:
                            eachSNP.append('0')
                    else:
                        eachSNP.append('0')

    return filteredSNPs


def reconstruct_fasta(filehandle, pad, names, list_snps):
    for i in tqdm(range(len(names)), ascii=True, desc='progress'):
        EachSample = names[i]
        sequence = ''
        for eachSNP in list_snps:
            if eachSNP[i + pad] == '1':
                sequence = sequence + eachSNP[2]
            elif eachSNP[i + pad] == '0':
                sequence = sequence + eachSNP[1]
            else:
                sequence = sequence + '?'
        filehandle.write('>' + EachSample + '\n' + sequence + '\n')
        i += 1
    filehandle.close()

