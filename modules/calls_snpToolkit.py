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
__version__ = '2.0.6'

from itertools import islice
from pathlib import Path
from datetime import datetime
from operator import itemgetter
from subprocess import Popen, PIPE, STDOUT
import glob
import ntpath
import sys
import os
import time
import logging
import coloredlogs
from tqdm import tqdm
from .argsLogger_snpToolkit import *
from .annotate_snpToolkit import *
from .combine_snpToolkit import *
# from .plot_snpToolkit import *
import pandas as pd
import pysam

logger = setupLogger()

def annotate(options):

    dna = ('A', 'C', 'T', 'G')
    allowed_format = {".vcf", ".vcf.gz", ".vcf.zip"}
    FilesToProcess = [FILE for FILE in glob.glob(
        '*' + options.identifier + '*') if FILE.endswith(tuple(allowed_format))]

    if len(FilesToProcess) == 0:
        logger.error(
            'No input file detected! Please check your(s) file(s) name(s)...')
        sys.exit(0)

    if options.excludeCloseSNPs == 0:
        logger.error('You must provide a value > 0 with the -f option')
        sys.exit(0)

    regions_to_exclude = []
    if options.exclude != None:
        if Path(options.exclude).exists():
            regions_to_exclude = pd.read_csv(
                options.exclude, sep='\t', header=None).values.tolist()
        else:
            logger.error('The file ' + options.exclude + ' does not exist')
            sys.exit(0)

    if options.genbank != None:
        if Path(options.genbank).exists():
            try:            
                annotationDB = parse_genbank_file(options.genbank)
            except BaseException as error:
                logger.error(
                    'Something went wrong when trying to extract data from your genbank file: {}'.format(error))
                sys.exit(0)
        else:
            logger.error('The file ' + options.genbank + ' does not exist')
            sys.exit(0)

    known_origins = {"samtools", "GATK", "freeBayes"}

    logger.info('snpToolkit is filtering and annotating your SNPs')

    warnings = []
    flag = False
    indelFile = False
    for i in tqdm(range(len(FilesToProcess)), ascii=True, desc='progress'):
        VcfFile = FilesToProcess[i]
        File_name = VcfFile.split('.vcf')[0]
        VcfFile_object = VCFtoolBox(VcfFile)
        vcf_source = VcfFile_object.vcf_generator()
        if vcf_source not in known_origins:
            warnings.append('Are you sure the vcf file ' + File_name +
                            ' was generated using samtools-mpileup, gatk-HaplotypeCaller or freeBayes!'.format(vcf_source))
        else:
            all_variations = VcfFile_object.extract_all_variation()
            indels_inVCF = []
            SNPs_in_VCF = []
            for EachVar in all_variations:
                if EachVar[3] in dna and EachVar[4] in dna:
                    SNPs_in_VCF.append(EachVar)
                else:
                    indels_inVCF.append(EachVar)
            raw_number_snps = len(SNPs_in_VCF)
           
            info1 = '##Total number of SNPs before snpToolkit processing: {}'.format(
                raw_number_snps)

            SNPs_to_process = SNPfiltering(SNPs_in_VCF)

            if options.excludeCloseSNPs != None and len(regions_to_exclude) == 0:
                Fitered_SNPs = SNPs_to_process.FilterCloseSNPs(
                    options.excludeCloseSNPs)
                info2 = '##By excluding SNPs that are closer than {} bp to each other, the number of remaining SNPs is: {}'.format(
                    options.excludeCloseSNPs, len(Fitered_SNPs))

            elif options.excludeCloseSNPs == None and len(regions_to_exclude) > 0:
                Fitered_SNPs = SNPs_to_process.FilterSNPtoExclude(
                    regions_to_exclude)
                info2 = '##By excluding SNPs that are located in the provided regions to exclude, the number of remaining SNPs is: {}'.format(
                    len(Fitered_SNPs))

            elif options.excludeCloseSNPs != None and len(regions_to_exclude) > 0:
                Fitered_SNPs = SNPs_to_process.FilterCloseSNPsandInREgions(
                    options.excludeCloseSNPs, regions_to_exclude)[1]
                info2 = '##By excluding SNPs that are closer than {} bp to each other, and also located in the provided regions to exclude, the number of remaining SNPs is: {}'.format(
                    options.excludeCloseSNPs, len(Fitered_SNPs))

            else:
                Fitered_SNPs = SNPs_in_VCF
                info2 = '##The options -f and -e were not used'


            PreProcessedRawSNPs = SNPselect(
                Fitered_SNPs, 0, 0, 0)
            extractedRawSNPs = PreProcessedRawSNPs.ExtractSNPinfo(vcf_source)

            PreProcessedSNPs = SNPselect(
                Fitered_SNPs, options.quality, options.depth, options.ratio)
            extractedSNPs = PreProcessedSNPs.ExtractSNPinfo(vcf_source)

            if options.ratio == 0.001:
                info3 = '##Filtred SNPs. Among the {} SNPs, the number of those with a quality score >= {}, a depth >= {} and a ratio >= {} is: {}'.format(
                    len(Fitered_SNPs), options.quality, options.depth, 0.0, sum([len(x) for x in extractedSNPs[0].values()]))
            else:
                info3 = '##Filtred SNPs. Among the {} SNPs, the number of those with a quality score >= {}, a depth >= {} and a ratio >= {} is: {}'.format(
                    len(Fitered_SNPs), options.quality, options.depth, options.ratio, sum([len(x) for x in extractedSNPs[0].values()]))

            genbank_accessions = annotationDB.keys()
            gbank_accessions_order = []
            for eachElem1 in annotationDB['locus_order']:
                for eachElem2 in genbank_accessions:
                    if eachElem1 in eachElem2:
                        gbank_accessions_order.append(eachElem2)
            SNPs_To_Map_and_Annotate = SNPannotation(extractedSNPs[0])
            Final_SNP_List = SNPs_To_Map_and_Annotate.mapAndannoSNPs(
                annotationDB)
            ##################################################################
            PreProcessedINDELS = SNPselect(
                indels_inVCF, options.quality, options.depth, options.ratio)
            extractedINDELS = PreProcessedINDELS.ExtractSNPinfo(vcf_source)
            INDELS_To_Map_and_Annotate = SNPannotation(extractedINDELS[0])
            Final_INDELS_List = INDELS_To_Map_and_Annotate.mapAndannoINDELS(
                annotationDB)
            indels_found = False
            for elem in Final_INDELS_List.keys():
                if len(Final_INDELS_List[elem]) > 1:
                    indels_found = True
                    break
            if indels_found is True:
                if indelFile is False:
                    directory_indels = 'snpToolkit_indels_output_' + \
                        datetime.now().strftime("%Y-%m-%d %H:%M:%S").replace(' ', '_at_').replace(':', '.')
                    indelFile = True

                    cmd = ['mkdir', directory_indels]
                    proc = Popen(cmd, stdout=PIPE, stderr=STDOUT)
                    proc.wait()
                    outputfile = open(directory_indels + '/' +
                                      File_name + '_snpToolkit_indels.txt', 'w')
                else:
                    cmd = ['mkdir', directory_indels]
                    proc = Popen(cmd, stdout=PIPE, stderr=STDOUT)
                    proc.wait()
                    outputfile = open(directory_indels + '/' +
                                      File_name + '_snpToolkit_indels.txt', 'w')

                for eachElem in Final_INDELS_List.keys():
                    for eachIndel in Final_INDELS_List[eachElem]:
                        outputfile.write(
                            '\t'.join([str(e) for e in eachIndel]) + '\n')
            ##################################################################

            info4 = '##After mapping, SNPs were located in: \n##' + \
                '\n##'.join(Final_SNP_List.keys())
            info5 = '##The mapped and annotated SNPs are distributed as follow:'

            summary_columns = [['Location', 'Genes', 'RBS', 'tRNA', 'rRNA', 'ncRNA',
                                'Pseudogenes', 'intergenic', 'Synonymous', 'NonSynonumous']]
            allSNPs = []
            all_data_collection = []
            for eachAccession in genbank_accessions:
                for location in Final_SNP_List.keys():
                    if eachAccession == location:
                        data_collection = []
                        snp_distribution = Final_SNP_List[location]
                        genomic_regions = snp_distribution.keys()
                        counts = [location]
                        Synonymous = 0
                        NonSynonumous = 0
                        for eachRegion1 in summary_columns[0][1:]:
                            for eachRegion2 in genomic_regions:
                                if eachRegion1 == eachRegion2:
                                    counts.append(
                                        format(len(Final_SNP_List[location][eachRegion1])))
                                    for eachSNP in Final_SNP_List[location][eachRegion1]:
                                        eachSNP.append(eachAccession)
                                        data_collection.append(eachSNP)
                                        allSNPs.append(eachSNP)
                                        if eachSNP[-2] == 'Syn':
                                            Synonymous += 1
                                        elif eachSNP[-2] == 'NS':
                                            NonSynonumous += 1
                        counts.append(format(Synonymous))
                        counts.append(format(NonSynonumous))
                        summary_columns.append(counts)
                        all_data_collection.append(data_collection)

            if len(all_data_collection) > 0:
                if flag is False:
                    directory = 'snpToolkit_SNPs_output_' + \
                        datetime.now().strftime("%Y-%m-%d %H:%M:%S").replace(' ', '_at_').replace(':', '.')
                    cmd = ['mkdir', directory]
                    proc = Popen(cmd, stdout=PIPE, stderr=STDOUT)
                    proc.wait()
                    flag = True
                    outputfile = open(directory + '/' +
                                      File_name + '_snpToolkit_SNPs.txt', 'w')
                    infos = ['##snpToolkit=version '+__version__, '##commandline= ' + ' '.join(
                        sys.argv[:]), '##VcfFile='+VcfFile, info1, info2, info3, info4, info5]
                    for info in infos:
                        outputfile.write(info + '\n')
                    outputfile.write('##'+'\t'.join(summary_columns[0]) + '\n')
                    for info in summary_columns[1:]:
                        outputfile.write('##SNPs in '+'\t'.join(info) + '\n')

                    outputfile.write(
                        '##Syn=Synonymous NS=Non-Synonymous'+'\n'+'##'+'\n')

                    header = ['Coordinates', 'REF', 'SNP', 'Depth', 'Nb of reads REF', 'Nb reads SNPs', 'Ratio', 'Quality', 'Location', 'Product',
                              'Orientation', 'Coordinates in gene', 'Ref codon', 'SNP codon', 'Ref AA', 'SNP AA', 'Coodinates protein', 'Effect', 'Location']

                    outputfile.write('\t'.join(header) + '\n')

                    for eachData in all_data_collection:
                        data_collection_sorted = sorted(
                            eachData, key=itemgetter(0))
                        for eachSNP in data_collection_sorted:
                            outputfile.write(
                                '\t'.join([format(eachElem) for eachElem in eachSNP]) + '\n')

                    outputfile.close()

                else:
                    outputfile = open(directory + '/' +
                                      File_name + '_snpToolkit_SNPs.txt', 'w')
                    infos = ['##snpToolkit=version '+__version__, '##commandline= ' + ' '.join(
                        sys.argv[:]), '##VcfFile='+VcfFile, info1, info2, info3, info4, info5]
                    for info in infos:
                        outputfile.write(info + '\n')
                    outputfile.write('##'+'\t'.join(summary_columns[0]) + '\n')
                    for info in summary_columns[1:]:
                        outputfile.write('##SNPs in '+'\t'.join(info) + '\n')

                    outputfile.write(
                        '##Syn=Synnonymous NS=Non-Synonymous'+'\n')

                    header = ['Coordinates', 'REF', 'SNP', 'Depth', 'Nb of reads REF', 'Nb reads SNPs', 'Ratio', 'Quality', 'Location', 'Product',
                                'Orientation', 'Coordinates in gene', 'Ref codon', 'SNP codon', 'Ref AA', 'SNP AA', 'Coodinates protein', 'Effect', 'Location']

                    outputfile.write('\t'.join(header) + '\n')

                    for eachData in all_data_collection:
                        data_collection_sorted = sorted(
                            eachData, key=itemgetter(0))
                        for eachSNP in data_collection_sorted:
                            outputfile.write(
                                '\t'.join([format(eachElem) for eachElem in eachSNP]) + '\n')

                    outputfile.close()

            else:
                warnings.append(
                    '{}: No SNPs were detected based on the specified filtering criteria'.format(File_name))

    if flag == True:
        if len(warnings) > 0:
            for warning in warnings:
                logger.warning(warning)
        logger.info(
            'snpToolkit output will be located in folder {}'.format(directory))
    else:
        if len(warnings) == True:
            for warning in warnings:
                logger.warning(warning)
        logger.info('snpToolkit output folder was not created')

 
       

def combine(options):
    regions_to_exclude = []
    if options.exclude is not None:
        import yaml
        if Path(options.exclude).exists():
            with open(options.exclude, 'r') as fh:
                try:
                    regions_to_exclude = yaml.load(fh)
                except yaml.YAMLError as exc:
                    logger.error(exc)
        else:
            logger.error('The file ' + options.exclude + ' does not exist')
            sys.exit(0)

    prefix = "_snpToolkit_SNPs.txt"
    choices = {'ns': 'Non-synonymous', 's': 'Sysnonymous',
               'inter': 'intergenic', 'all': 'SNPs'}

    checkOutput = Path(choices[options.snps] + '_polymorphic_sites.txt')
    if checkOutput.exists():
        logger.warning(choices[options.snps] + '_polymorphic_sites.txt exists already and was created on ' + time.ctime(
            checkOutput.stat().st_ctime) + '. This file will be overwritten. Press any key to continue or ctrl-c to exit!')
        try:
            input()
        except KeyboardInterrupt:
            sys.exit(0)

    if options.bamFolder != None:
        BamFilesToInclude = checksnpToolkitBam(options.bamFolder)
        if len(BamFilesToInclude) == 0:
            logger.error('None of the bam files names in the folder {} correspond to any of the snpToolkit outputs'.format(
                options.bamFolder))
            sys.exit(0)
        else:
            BamBai = checkBamindex(options.bamFolder)
            if len(BamBai) > 0:
                logger.error('index file missing for: \n' + '\n'.join(BamBai))
                sys.exit(0)

    All_snpToolkitFiles = glob.glob("*" + prefix)
    snpToolkitFiles = []
    for eachFile in All_snpToolkitFiles:
        if os.stat(eachFile).st_size == 0:
            logger.warning(eachFile+ ' is empty and will not be processed!')
        else:
            snpToolkitFiles.append (eachFile)

    if len(snpToolkitFiles) == 0:
        logger.error('No snpToolkit output files were detected')
        sys.exit(0)

    elif len(snpToolkitFiles) == 1:
        logger.error('At least two snpToolkit output files must be detected')
        sys.exit(0)

    else:
        logger.info('Searching for polymorphic sites...')
        SampleNames = []
        polymorphic_sites = []
        for eachFile in snpToolkitFiles:
            SampleNames.append(ntpath.basename(eachFile).split(prefix)[0])
            with open(eachFile, 'r') as f:
                for l in f.readlines():
                    if l[:2] != '##':
                        content = l.strip().split('\t')
                        if float(content[6]) >= float(options.ratio):
                            snp = [int(content[0])] + \
                                content[1:3] + content[8:-1]
                            if snp not in polymorphic_sites:
                                if options.location in l.split('\t')[-1] or options.location == l.split('\t')[-1]:
                                    if options.snps == 'all':
                                        polymorphic_sites.append(snp)
                                    elif options.snps == 's':
                                        if l.split('\t')[17] == 'Syn':
                                            polymorphic_sites.append(snp)
                                    elif options.snps == 'ns':
                                        if l.split('\t')[17] == 'NS':
                                            polymorphic_sites.append(snp)
                                    else:
                                        if l.split('\t')[8] == 'intergenic':
                                            polymorphic_sites.append(snp)

        if len(polymorphic_sites) == 0:
            logger.error(
                'No polymorphic sites found! are you sure you provided the right name when using the --location option?')
            sys.exit(0)

        polymorphic_sites.sort()
        header = '##ID Coordinates REF SNP Location Product Orientation NucPosition REF-codon NEW-codon REF-AA NEW-AA ProPostion Type'.split(
            ' ')
        infos = ['##snpToolkit=version '+__version__, '##commandline= ' + ' '.join(
            sys.argv[:]), '##location='+options.location, '##Number of polymorphic sites= ' + str(len(polymorphic_sites)), '\t'.join(header + SampleNames)]

        logger.info(choices[options.snps] +
                    ' polymorphic sites distribution. Please wait...')

        if options.bamFolder == None:
            distribution_result = snp_distribution(
                polymorphic_sites, snpToolkitFiles, regions_to_exclude)

        else:
            distribution_result = snp_distribution_missing(
                options.location, options.ratio, options.snps, polymorphic_sites, snpToolkitFiles, options.bamFolder, BamFilesToInclude, options.depth, regions_to_exclude)

        logger.info(
            'Creating ' + choices[options.snps] + '_polymorphic_sites.txt')
        outputfile1 = open(choices[options.snps] +
                           '_polymorphic_sites.txt', 'w')
        
        
        


        
        for info in infos:
            outputfile1.write(info + '\n')

        reference = ''
        i = 1
        for eachSNP in distribution_result:
            reference = reference + eachSNP[1]
            outputfile1.write('snp{}'.format(i) + '\t' +
                              '\t'.join([str(x) for x in eachSNP]) + '\n')
            i += 1
        outputfile1.close()

        logger.info('Creating ' + choices[options.snps] + '_alignment.fasta')
        outputfile2 = open(choices[options.snps] + '_alignment.fasta', 'w')
        outputfile2.write('>' + options.location + '\n' + reference + '\n')
        reconstruct_fasta(outputfile2, 13, SampleNames, distribution_result)



def expand(options):
    FilesToAdd = [FILE for FILE in glob.glob(
        options.directory + '/*_snpToolkit_SNPs.txt')]

    with open(options.polymorphic_sites, 'r') as input:
        PolyMorphicSites = [l.strip().split('\t') for l in input]

    firstSampleSet = PolyMorphicSites[4][14:]
    aDNAcoordiantes = []
    for s in aDNA:
        i = 0
        while i < len(firstSampleSet):
            if s == firstSampleSet[i]:
                aDNAcoordiantes.append(i)
            i += 1

    NewaDNAcoordiantes = []
    for s in NewaDNA:
        i = 0
        while i < len(firstSampleSet):
            if s == firstSampleSet[i]:
                NewaDNAcoordiantes.append(i)
            i += 1

    with open(options.output+'-combined.txt', 'w') as input:
        for elem in PolyMorphicSites:
            input.write('\t'.join(elem) + '\n')

    isHEADER = False
    for i in tqdm(range(len(FilesToAdd)), ascii=True, desc='progress'):
        with open(options.output+'-combined.txt', 'r') as input:
            PolyMorphicSites = [l.strip().split('\t') for l in input]
        if isHEADER is False:
            header = PolyMorphicSites[4]
            nbSamples = len(PolyMorphicSites[4][14:])
            NbPoly = len(PolyMorphicSites[5:])
            SNPcontent = PolyMorphicSites[5:]
        else:
            header = PolyMorphicSites[0]
            nbSamples = len(PolyMorphicSites[0][14:])
            NbPoly = len(PolyMorphicSites[1:])
            SNPcontent = PolyMorphicSites[1:]

        sample_name = FilesToAdd[i].split('/')[-1].split('_snpToolkit')[0]
        with open(FilesToAdd[i], 'r') as f:
            SNPnewSample = [l.strip().split('\t')
                            for l in f if '##' not in l and options.location in l]

        # TODO: check for the new aDNA the position of polymorphic sites
        SNPpos = [x[0] for x in SNPnewSample]
        header.append(sample_name)
        for polySNP in SNPcontent:
            if polySNP[1] in SNPpos:
                polySNP.append('1')
                SNPpos.remove(polySNP[1])
            else:
                if sample_name in NewaDNA:
                    bamFile = pysam.AlignmentFile(
                        FilesToAdd[i].split('_snpToolkit')[0]+'.bam')
                    try:
                        depthATposition = 0
                        # bamFile.count_coverage() return four array.arrays of the same length in order A C G T
                        NucleotidesOrder = 'ACGT'
                        NumberOfReads = []
                        k = 0
                        for eachCov in bamFile.count_coverage(options.location, int(polySNP[1]) - 1, int(polySNP[1])):
                            depthATposition = depthATposition + \
                                list(eachCov)[0]
                            NumberOfReads.append(
                                (list(eachCov)[0], NucleotidesOrder[k]))
                            k += 1
                        NumberOfReads.sort(reverse=True)
                        if depthATposition < int(options.cutoff):
                            polySNP.append('?')
                        else:
                            if NumberOfReads[0][1] == polySNP[3]:
                                if NumberOfReads[0][0] / depthATposition >= 0.6:
                                    polySNP.append('1')
                                elif NumberOfReads[0][0] / depthATposition >= 0.4:
                                    polySNP.append('?')
                                else:
                                    polySNP.append('0')
                            else:
                                polySNP.append('0')
                    except ValueError as e:
                        logging.error(
                            'Please use samtools view -H on this bam file to check the exact name.')
                else:
                    polySNP.append('0')

        y = 0
        for eachSNP1 in SNPpos:
            for eachSNP2 in SNPnewSample:
                if eachSNP1 == eachSNP2[0]:
                    transfert = ['snp'+str(NbPoly+y)]

                    SNPinfo = eachSNP2[:3]+eachSNP2[8:-1]
                    SNPdist = []
                    k = 0
                    while k < nbSamples:
                        if str(k) in aDNAcoordiantes:
                            SNPdist.append('?')
                        else:
                            SNPdist.append('0')
                        k += 1
                    SNPdist.append('1')
                    combined = transfert+SNPinfo+SNPdist
                    SNPcontent.append(combined)
                    y += 1
        with open(options.output+'-combined.txt', 'w') as f:
            isHEADER = True
            f.write('\t'.join(header)+'\n')
            for eachSNP in SNPcontent:
                f.write('\t'.join(eachSNP)+'\n')

    with open(options.output+'-combined.txt', 'r') as f:
        AllPolyMorphicSites = [l.strip().split('\t') for l in f]
    samples = AllPolyMorphicSites[0]

    with open(options.output+'-reconstracted.fasta', 'w') as fastaoutput:

        i = 14
        while i < len(samples):

            sequence = ''
            for SNP in AllPolyMorphicSites[1:]:
                if SNP[i] == '0':
                    sequence = sequence + SNP[2]
                elif SNP[i] == '1':
                    sequence = sequence + SNP[3]
                else:
                    sequence = sequence + '?'

            fastaoutput.write('>'+samples[i]+'\n'+sequence+'\n')
            i = i+1



if __name__ == '__main__':
    app.run_server()