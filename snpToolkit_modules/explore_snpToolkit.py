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


import glob
import argparse
import logging
import coloredlogs
from tqdm import tqdm
import pandas as pd
import random
import plotly.graph_objs as go
import plotly.offline as pyo
import plotly.figure_factory as ff
import numpy as np 
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies  import Input, Output
import dash_table as dt
from operator import itemgetter
from subprocess import Popen, PIPE, STDOUT
from annotate_snpToolkit import *
from argsLogger_snpToolkit import *
from plotly.colors import n_colors


color_choice = "aliceblue, antiquewhite, aqua, aquamarine, azure, beige, bisque, black, blanchedalmond, blue, blueviolet, brown, burlywood, cadetblue, chartreuse, chocolate, coral, cornflowerblue,cornsilk, crimson, cyan, darkblue, darkcyan,darkgoldenrod, darkgray, darkgrey, darkgreen,darkkhaki, darkmagenta, darkolivegreen, darkorange,darkorchid, darkred, darksalmon, darkseagreen,darkslateblue, darkslategray, darkslategrey, darkturquoise, darkviolet, deeppink, deepskyblue,dimgray, dimgrey, dodgerblue, firebrick,floralwhite, forestgreen, fuchsia, gainsboro, ghostwhite, gold, goldenrod, gray, grey, green,greenyellow, honeydew, hotpink, indianred, indigo,ivory, khaki, lavender, lavenderblush, lawngreen,lemonchiffon, lightblue, lightcoral, lightcyan,lightgoldenrodyellow, lightgray, lightgrey,lightgreen, lightpink, lightsalmon, lightseagreen,lightskyblue, lightslategray, lightslategrey,lightsteelblue, lightyellow, lime, limegreen,linen, magenta, maroon, mediumaquamarine, mediumblue, mediumorchid, mediumpurple,mediumseagreen, mediumslateblue, mediumspringgreen,mediumturquoise, mediumvioletred, midnightblue, mintcream, mistyrose, moccasin, navajowhite, navy,oldlace, olive, olivedrab, orange, orangered,orchid, palegoldenrod, palegreen, paleturquoise,palevioletred, papayawhip, peachpuff, peru, pink,plum, powderblue, purple, red, rosybrown,royalblue, saddlebrown, salmon, sandybrown,seagreen, seashell, sienna, silver, skyblue,slateblue, slategray, slategrey, snow, springgreen,steelblue, tan, teal, thistle, tomato, turquoise,violet, wheat, white, whitesmoke, yellow,yellowgreen".split(',')


logger = setupLogger()

parser = argparse.ArgumentParser()
parser.add_argument('-i', dest='identifier', required=True, help='provide input file')


args= parser.parse_args()
dna = ('A', 'C', 'T', 'G')
allowed_format = {".vcf", ".vcf.gz", ".vcf.zip"}
FilesToProcess = [FILE for FILE in glob.glob(
    '*' + args.identifier + '*') if FILE.endswith(tuple(allowed_format))]

if len(FilesToProcess) == 0:
    logger.error(
        'No input file detected! Please check your(s) file(s) name(s)...')
    sys.exit(0)
known_origins = {"samtools", "GATK", "freeBayes"}
logger.info('snpToolkit is extracting your data and creating the different plots...')
warnings = []
flag = False
indelFile = False
vcf_data_collection = {}
isolates_collection = []
for i in tqdm(range(len(FilesToProcess)), ascii=True, desc='progress'):
    VcfFile = FilesToProcess[i]
    File_name = VcfFile.split('.vcf')[0]
    isolates_collection.append ({'label':File_name,'value':File_name})
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
        PreProcessedRawSNPs = SNPselect(SNPs_in_VCF, 0, 0, 0)
        extractedRawSNPs = PreProcessedRawSNPs.ExtractSNPinfo(vcf_source)
        vcf_data_collection [File_name] = extractedRawSNPs

diffContigs = []
for each_vcf_File in vcf_data_collection.keys():
    for eachElem in vcf_data_collection[each_vcf_File][0].keys():
        if eachElem not in diffContigs:
            diffContigs.append(eachElem)
diffContigs.sort()
color_picker={}
randNum = list(set(range(1, len(color_choice))))
random.shuffle(randNum)
uniqueRandNum = list(set(randNum))
for x in range (len (diffContigs)):
    color_picker[diffContigs[x]]= uniqueRandNum[x]

data ={}
data2 = {}
all_snps_stats = {}
for each_vcf_File in vcf_data_collection.keys():
    vcf_content = vcf_data_collection[each_vcf_File][0]
    allSNP = []
    allSNP2 = []
    statsSNPs =[]
    if len (vcf_content.keys())>1:
        colors = n_colors('rgb(5, 200, 200)', 'rgb(200, 10, 10)', len(vcf_content.keys()), colortype='rgb')
    else:
        colors = ['rgb(62, 155, 184)']
    snpsSTATS=[]
    for eachElem,color in zip(vcf_content.keys(),colors):
        df_raw = pd.DataFrame(vcf_content[eachElem],[x[0] for x in vcf_content[eachElem]],['Position','REF','SNP','Depth','Depth reference','Depth SNPs','Ratio','Quality'])
        snpsSTATS.append ([eachElem,len(df_raw),len(df_raw[df_raw.Depth >= 3]),len(df_raw[df_raw.Quality >= 20]),len(df_raw[df_raw.Ratio >= 0.9]),len(df_raw[(df_raw.Ratio >= 0.9) & (df_raw.Depth >= 3) & (df_raw.Quality >= 20)])])
        bubble_color = color_picker[eachElem]
        Qrange = [df_raw['Quality'].between(1, 20), df_raw['Quality'].between(21, 50), df_raw['Quality'].between(51, 100), df_raw['Quality'].between(101, 1000000)] 
        values = [5, 10, 15, 20]
        df_raw['Qvalues'] = np.select(Qrange, values)
        allSNP.append(go.Scatter(x=df_raw["Ratio"],y=df_raw["Depth"],mode="markers",name=eachElem,marker=dict(size=df_raw['Qvalues'],color=color,opacity=0.3,showscale=False),showlegend=True))#color=df_raw['Quality']
        allSNP2.append(go.Violin(x=df_raw["Ratio"],name=eachElem,marker=dict(opacity=0.5,color=color),points=False,orientation='h', side='positive',showlegend=True))

    all_snps_stats[each_vcf_File]=pd.DataFrame(columns=['Location', 'Total SNPs','Nb SNP if depth >= 3','Nb SNP if QUAL >= 20','Nb SNP if Ratio >= 0.9','All filters on'],data=snpsSTATS)

    data[each_vcf_File]=allSNP
    data2[each_vcf_File]=allSNP2



app = dash.Dash()

colors={'background':'#111111','text':'#ffffff','textHeader':'#111111'}
app.layout = html.Div(children=[
    html.H2('snpToolkit plots -  VCFs SNPs content',style={'color':colors['textHeader'],'text-align':'center'}),
    dcc.Dropdown(id='isolate',options=isolates_collection,value=isolates_collection[0]['value']),
    html.H3('Number of identified SNPs',style={'color':colors['textHeader']}),
    html.Div(dt.DataTable(columns = [{"name": i, "id": i} for i in ['Location', 'Total SNPs','Nb SNP if depth >= 3','Nb SNP if QUAL >= 20','Nb SNP if Ratio >= 0.9','All filters on']],id='snpTable', style_table={ 'align':'center'},style_cell={'textAlign': 'center', 'border': '1px solid grey'},    style_header={
        'color':'white',
        'backgroundColor': 'rgb(62, 155, 184)',
        'border': '1px solid black'
    })),
    dcc.Graph(id='graph'),
    html.Div(children=[
    dcc.Graph(id='graph2')
    ])])


@app.callback (Output('snpTable','data'),[Input('isolate','value')])
def update_table(selected_isolate):
    return all_snps_stats[selected_isolate].to_dict('rows')

@app.callback (Output('graph','figure'),[Input('isolate','value')])
def update_graph(selected_isolate):
    layout = go.Layout(title='Depth vs Ratio',xaxis={'title':'Ratio'},yaxis={'title':'Depth'},hovermode='closest')
    return {'data':data[selected_isolate],'layout':layout}

@app.callback (Output('graph2','figure'),[Input('isolate','value')])
def update_graph2(selected_isolate):
    layout2 = go.Layout(title='Frequency distribution',xaxis={'title':'Ratio','range':[0,1]}, yaxis={'showticklabels':False}, hovermode=False,autosize=False,height=300,xaxis_showgrid=False, xaxis_zeroline=False)
    return {'data':data2[selected_isolate],'layout':layout2}



if __name__ == '__main__':
    app.run_server()
    

