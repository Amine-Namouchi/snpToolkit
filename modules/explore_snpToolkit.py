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
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies  import Input, Output
from operator import itemgetter
from subprocess import Popen, PIPE, STDOUT
from annotate_snpToolkit import *
from argsLogger_snpToolkit import *
from argsLogger_snpToolkit import *


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
for each_vcf_File in vcf_data_collection.keys():
    vcf_content = vcf_data_collection[each_vcf_File][0]
    allSNP = []
    for eachElem in vcf_content.keys():
        df_raw = pd.DataFrame(vcf_content[eachElem],[x[0] for x in vcf_content[eachElem]],['Position','REF','SNP','Depth','Depth reference','Depth SNPs','Ratio','Quality'])
        allSNP.append(go.Scatter(x=df_raw["Ratio"],y=df_raw["Depth"],mode="markers",name=eachElem,marker=dict(color=color_picker[eachElem],size=df_raw['Quality'],opacity=0.5,sizeref=5, line_width=1,showscale=False),showlegend=True))#color=df_raw['Quality']
        #allSNP.append(go.Histogram(x=df_raw["Ratio"],name=eachElem,marker=dict(color='black',opacity=0.2),showlegend=True))
    data[each_vcf_File]=allSNP


app = dash.Dash()

colors={'background':'#111111','text':'#ffffff','textHeader':'#111111'}

app.layout = html.Div(children=[
    html.H1('snpToolkit plots',style={'color':colors['textHeader'],'text-align':'center'}),
    dcc.Graph(id='graph'),
    dcc.Dropdown(id='isolate',options=isolates_collection,value=File_name)])
    


@app.callback (Output('graph','figure'),[Input('isolate','value')])

def update_graph(selected_isolate):
    layout = go.Layout(title='Depth vs Ratio',xaxis={'title':'Ratio'},yaxis={'title':'Depth'},hovermode='closest')
    return {'data':data[selected_isolate],'layout':layout}


if __name__ == '__main__':
    app.run_server()
    

