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




#TODO: adjust the size of the input text. color dots in scatter according to SNPS effect.


parser = argparse.ArgumentParser()
parser.add_argument('--dir', dest='directory', required=True, help='provide the path to the directory containing snptoolkit output files')
args= parser.parse_args()

list_files = glob.glob(args.directory+'/*_snpToolkit_SNPs.txt')
snptoolkitInputs =[]
for f in list_files:
    snptoolkitInputs.append ({'label':f.split('/')[-1].split('_snpToolkit_')[0],'value':f.split('/')[-1].split('_snpToolkit_')[0]})


colors={'background':'#262626','text':'#ffffff','textHeader':'#111111'}
app = dash.Dash()


header = ['##Coordinates', 'REF', 'SNP', 'Ref codon', 'SNP codon', 'Ref AA', 'SNP AA','Effect','Annotation', 'Product',
            'Orientation', 'Coordinates in gene',  'Coordinates protein','Depth', 'Nb of reads REF', 'Nb reads SNPs', 'Ratio', 'Quality']

app.layout = html.Div(
    [
        html.Div([
        html.H1("snpToolkit",style={'textAlign': 'center'}),
        html.H3("Visualization of snpToolkit annotated output file",style={'textAlign': 'center', 'color': '#5e5e5e'}),
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
        #TODO: dcc.Input(placeholder='Enter a window size to compute SNPs density...',type='number',id="windowsize",value=''),
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
                headerFile = islice(f, 30)
                for c in headerFile:
                    if '##SNPs in ' in c:  
                        SNPsLocation.append ({'label': c.strip().replace('##SNPs in ','').split('\t')[0], 'value': c.strip().replace('##SNPs in ','').split('\t')[0]})
    
    return SNPsLocation


@app.callback(
    Output('opt-dropdown', 'label'),
    [Input('opt-dropdown', 'options')]
)
def default_selection(option2):
    if len (option2) > 1:
        return ''
    else:
        return option2[0]['value']


@app.callback (Output('graph','figure'),[Input('name-dropdown','value'),Input('opt-dropdown','value'),Input('selected-type','value')])##TODO: Input('windowsize','value')
def update_graph(sample,location,snpType):
    for eachFile in list_files:
        if sample in eachFile:
            fh = open (eachFile,'r')
            line = fh.readline()
            x = 0
            while x !='':
                if line.startswith('##Coordinates')==True:
                    break
                x+=1
                line = fh.readline()
            df = pd.read_csv(eachFile,sep='\t',skiprows=x)

    df1 = df.loc[(df['Effect'] == 'NS') & (df['Location'] == location)][header]
    df2 = df.loc[(df['Effect'] == 'Syn') & (df['Location'] == location)][header]

    layout = go.Layout(title='Depth vs Ratio',xaxis={'title':'Coordinates'},yaxis={'title':'Depth'},hovermode='closest')
    if snpType == 'ALL':
        data = [go.Scatter(name="Others",x=df[(df["Location"]==location)]["##Coordinates"],y=df[(df["Location"]==location)]["Depth"],mode="markers",opacity=0.5,marker={"color":"grey"}),go.Scatter(name='NS',x=df1["##Coordinates"],y=df1["Depth"],mode="markers",opacity=0.8,marker={"color":"#FBBF4C"}),go.Scatter(name='Syn',x=df2["##Coordinates"],y=df2["Depth"],mode="markers",opacity=0.8,marker={"color":"#51A8C7"})]
    elif snpType == 'NS':
        
        data = [go.Scatter(x=df1["##Coordinates"],y=df1["Depth"],mode="markers",opacity=0.8,marker={"color":"#FBBF4C"})]
    else: 
        data = [go.Scatter(x=df2["##Coordinates"],y=df2["Depth"],mode="markers",opacity=0.8,marker={"color":"#51A8C7"})]
    return {'data':data,'layout':layout}


@app.callback (Output('snpTable','data'),[Input('name-dropdown','value'),Input('opt-dropdown','value')])
def update_table(sample,location):
    for eachFile in list_files:
        if sample in eachFile:
            fh = open (eachFile,'r')
            line = fh.readline()
            x = 0
            while x !='':
                if line.startswith('##Coordinates')==True:
                    break
                x+=1
                line = fh.readline()
            df = pd.read_csv(eachFile,sep='\t',skiprows=x)
    return df[(df['Location'] == location)][header].to_dict('rows')
        

if __name__ == '__main__':
    app.run_server()
