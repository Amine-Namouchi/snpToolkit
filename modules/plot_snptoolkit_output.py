#!/usr/bin/env python
import glob
import logging
import coloredlogs
import argparse
import plotly.graph_objs as go
import plotly.offline as pyo
import plotly.figure_factory as ff
from plotly.colors import n_colors
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies  import Input, Output
import dash_table as dt
from itertools import islice
import pandas as pd 




#TODO: adjust the size of the input text. color dots in scatter according to SNPS effect. density plot and table filter using coordinates not working!


parser = argparse.ArgumentParser()
parser.add_argument('--dir', dest='directory', required=True, help='provide the path to the directory containing snptoolkit output files')
args= parser.parse_args()

list_files = glob.glob(args.directory+'/*_snpToolkit_SNPs.txt')
snptoolkitInputs =[]
for f in list_files:
    snptoolkitInputs.append ({'label':f.split('/')[-1].split('_snpToolkit_')[0],'value':f.split('/')[-1].split('_snpToolkit_')[0]})


colors={'background':'#262626','text':'#ffffff','textHeader':'#111111'}
app = dash.Dash()

# snpToolkit=version 2.0.6
# commandline= ../snpToolkit annotate -i vcf.gz -g /Users/amine/Documents/from_google_drive/HAARLEM/reference/GCF_000195955.2_ASM19595v2_genomic.gbff -r 0.9 -f 50
# VcfFile=TUN228.vcf.gz
# Total number of SNPs before snpToolkit processing: 1248
# By excluding SNPs that are closer than 50 bp to each other, the number of remaining SNPs is: 915
# Filtred SNPs. Among the 915 SNPs, the number of those with a quality score >= 20, a depth >= 3 and a ratio >= 0.9 is: 874
# After mapping, SNPs were located in: 
# NC_000962.3: Mycobacterium tuberculosis H37Rv, complete genome 4411532 bp
# The mapped and annotated SNPs are distributed as follow:
# Location	Genes	RBS	tRNA	rRNA	ncRNA	Pseudogenes	intergenic	Synonymous	NonSynonumous
# SNPs in NC_000962.3: Mycobacterium tuberculosis H37Rv, complete genome 4411532 bp	763	0	0	0	2	0	111	304	459
# Syn=Synonymous NS=Non-Synonymous

# Coordinates	Ref	SNP	Depth	Nb of reads Ref	Nb reads SNPs	Ratio	Quality	Location	Product	Orientation	Coordinates annotation	Ref codon	SNP codon	Ref AA	SNP AA	Coodinates Protein	Effect	Distribution
header = ['Coordinates', 'REF', 'SNP', 'Depth', 'Nb of reads REF', 'Nb reads SNPs', 'Ratio', 'Quality', 'Annotation', 'Product',
            'Orientation', 'Coordinates in gene', 'Ref codon', 'SNP codon', 'Ref AA', 'SNP AA', 'Coodinates protein', 'Effect', 'Location']

app.layout = html.Div(
    [
        html.Div([
        dcc.Dropdown(
            id='name-dropdown',
            options=snptoolkitInputs,
            value = snptoolkitInputs[0]['value']
            ),
            ],style={'width': '50%', 'display': 'inline-block'}),
        html.Div([
        dcc.Dropdown(
            id='opt-dropdown',
            ),
            ],style={'width': '50%', 'display': 'inline-block'}
        ),
        html.Hr(),
        dcc.RadioItems(id='selected-type',options=[
            {'label': 'All SNPs', 'value': 'ALL'},
            {'label': 'Synonymous SNPs', 'value': 'SYN'},
            {'label': 'Non-synonymous SNPs', 'value': 'NS'},
            ],value='ALL',labelStyle={'display': 'inline-block'}
        ),
        dcc.Input(placeholder='Enter a window size to compute SNPs density...',type='number',id="windowsize",value=''),
        dcc.Graph(id='graph'),
        html.Hr(),
        html.Div(dt.DataTable(columns = [{"name": i, "id": i} for i in header],
        id='snpTable', 
        filter_action="native",
        page_size=20,
        style_table={ 'align':'center','overflowY': 'scroll'},
        style_cell={'textAlign': 'center', 'border': '1px solid grey'},
        style_header={
        'color':'white',
        'backgroundColor': 'rgb(62, 155, 184)',
        'border': '1px solid black'
        },
        style_data_conditional=[
            {
                'if': {
                'column_id': 'Effect',
                'filter_query': '{Effect} eq "NS"'
            },
            'backgroundColor': '#FBBF4C',
            'color': 'white',
        },
        {
                'if': {
                'column_id': 'Effect',
                'filter_query': '{Effect} eq "Syn"'
            },
            'backgroundColor': '#51A8C7',
            'color': 'white',
        },
        {
                'if': {
                'column_id': 'Effect',
                'filter_query': '{Effect} eq "."'
            },
            'backgroundColor': '#0e2f2e',
            'color': '#0e2f2e',
        }
        ])),
    ]
)
@app.callback(
    Output('opt-dropdown', 'options'),
    [Input('name-dropdown', 'value')]
)
def update_dropdown(name):
    SNPsLocation =[]
    for eachFile in list_files:
        if name in eachFile:
            with open (eachFile, 'r') as f:
                headerFile = islice(f, 12)
                for c in headerFile:
                    if '##SNPs in ' in c:     
                        SNPsLocation.append ({'label': c.replace('##SNPs in ','').split('\t')[0], 'value': c.replace('##SNPs in ','').split('/t')[0]})
    return SNPsLocation


@app.callback(
    Output('opt-dropdown', 'value'),
    [Input('opt-dropdown', 'options')]
)
def default_selection(option2):
    if len (option2) > 1:
        return ''
    else:
        return option2[0]['value']


@app.callback (Output('graph','figure'),[Input('name-dropdown','value'),Input('opt-dropdown','value'),Input('selected-type','value'),Input('windowsize','value')])
def update_graph(sample,location,snpType,window):
    for eachFile in list_files:
        if sample in eachFile:
            df = pd.read_csv(eachFile,sep='\t',skiprows=12)
    print (df["Depth"].rolling(10000).mean())
    df1 = df.loc[df['Effect'] == 'NS']
    df2 = df.loc[df['Effect'] == 'Syn']
    layout = go.Layout(title='Depth vs Ratio',xaxis={'title':'Coordinates'},yaxis={'title':'Depth'},hovermode='closest')
    if snpType == 'ALL':
        data = [go.Scatter(name="Others",x=df["Coordinates"],y=df["Depth"],mode="markers",opacity=0.5,marker={"color":"grey"}),go.Scatter(name='NS',x=df1["Coordinates"],y=df1["Depth"],mode="markers",opacity=0.8,marker={"color":"#FBBF4C"}),go.Scatter(name='Syn',x=df2["Coordinates"],y=df2["Depth"],mode="markers",opacity=0.8,marker={"color":"#51A8C7"})]
    elif snpType == 'NS':
        
        data = [go.Scatter(x=df1["Coordinates"],y=df1["Depth"],mode="markers",opacity=0.8,marker={"color":"#FBBF4C"})]
    else 
        data = [go.Scatter(x=df2["Coordinates"],y=df2["Depth"],mode="markers",opacity=0.8,marker={"color":"#51A8C7"})]
    return {'data':data,'layout':layout}


@app.callback (Output('snpTable','data'),[Input('name-dropdown','value'),Input('opt-dropdown','value')])
def update_table(sample,location):
    for eachFile in list_files:
        if sample in eachFile:
            df = pd.read_csv(eachFile,sep='\t',skiprows=12)
    return df.to_dict('rows')
        
            
 


if __name__ == '__main__':
    app.run_server()
