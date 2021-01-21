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
import dash_daq as daq
from dash.dependencies  import Input, Output
import dash_table as dt
from itertools import islice
import pandas as pd 
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import umap
import plotly.express as px




parser = argparse.ArgumentParser()
parser.add_argument('-p', required=False, dest='polymorphic_sites',
                                help='provide the path of the polymorphic sites you want to analyze')
parser.add_argument('-conf', required=False, dest='config',
                                help='provide the path of the configuration file that contains the information to use for data visualization')
args= parser.parse_args()

col_to_delete = 'Coordinates	REF	SNP	Location	Product	Orientation	NucPosition	REF-codon	NEW-codon	REF-AA	NEW-AA	ProPostion	Type'.split('\t')
df1 = pd.read_csv(args.polymorphic_sites, skiprows=4,sep='\t',index_col=0)
for col in col_to_delete:
    del df1[col]

features =[]
if args.config != None:
    metadata = pd.read_csv(args.config,sep='\t',index_col=0)
    for f in list(metadata.columns):
        features.append ({'label':f,'value':f})
else:
    metadata = pd.DataFrame([],[])
    features.append ({'label':'','value':''})

df_transpose = df1.T

#################### PCA ###################
scaler = StandardScaler()
scaler.fit(df_transpose)
scaled_data = scaler.transform(df_transpose)
pca = PCA(n_components=2)
pca.fit(scaled_data)
projections = pca.transform(scaled_data)
pca_df = pd.DataFrame({'X':projections[:,0],'Y':projections[:,1]},df_transpose.index)
if args.config != None:
    dfmerge_pca=pd.merge(pca_df,metadata,left_index=True, right_index=True)
    dfmerge_pca['samples'] = list(dfmerge_pca.index)
else:
    dfmerge_pca=pca_df
    dfmerge_pca['samples'] = list(dfmerge_pca.index)

#################### PCA ###################


#################### UMAP ###################
umap_2d = umap.UMAP(n_components=2,random_state=0)
umap_obj = umap_2d.fit(df_transpose)
projections = umap_obj.transform(df_transpose)
umap_df = pd.DataFrame({'X':projections[:,0],'Y':projections[:,1]},df_transpose.index)
if args.config != None:
    dfmerge_umap=pd.merge(umap_df,metadata,left_index=True, right_index=True)
    dfmerge_umap['samples'] = list(dfmerge_umap.index)
else:
    dfmerge_umap=umap_df
    dfmerge_umap['samples'] = list(dfmerge_umap.index)
#################### UMAP ###################



app = dash.Dash()
debug=True
if args.config != None:
    app.layout = html.Div(
        [
            html.Div([
            html.H1("snpToolkit",style={'textAlign': 'center'}),
            html.H3("Dimension reduction using PCA and UMAP",style={'textAlign': 'center', 'color': '#5e5e5e'}),
            html.Hr(),
            dcc.RadioItems(id='selected-type',options=[
                {'label': 'PCA', 'value': 'pca'},
                {'label': 'UMAP', 'value': 'umap'},
                ],value='pca',labelStyle={'display': 'inline-block'}, 
                style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}
            ),
            html.H4("Select a feature from the dropdown menu for color mapping:",style={'color': '#5e5e5e'}),
            dcc.Dropdown(
                id='name-dropdown',
                options=features),
                ],style={'width': '100%', 'display': 'inline-block'}),
        
            html.Div([
            html.Br(),
            dcc.Graph(id='graph'),
            ]),        
            daq.ToggleSwitch(
            id='mode-switch',
            labelPosition='left',
            label='Dark mode',
            value=False,
            style={'display': 'flex', 'align-items': 'left', 'justify-content': 'left'},
            vertical=False
            ),
        ])

else:
    app.layout = html.Div(
        [
            html.Div([
            html.H1("snpToolkit",style={'textAlign': 'center'}),
            html.H3("Dimension reduction using PCA and UMAP",style={'textAlign': 'center', 'color': '#5e5e5e'}),
            html.Hr(),
            dcc.RadioItems(id='selected-type',options=[
                {'label': 'PCA', 'value': 'pca'},
                {'label': 'UMAP', 'value': 'umap'},
                ],value='pca',labelStyle={'display': 'inline-block'}, 
                style={'display': 'flex', 'align-items': 'center', 'justify-content': 'center'}
            ),
            html.H4("No metadata was provided for color mapping",style={'color': '#5e5e5e'}),
            dcc.Dropdown(
                id='name-dropdown',
                options=features),
                ],style={'width': '100%', 'display': 'inline-block'}),
        
            html.Div([
            html.Br(),
            dcc.Graph(id='graph'),
            ]),        
            daq.ToggleSwitch(
            id='mode-switch',
            labelPosition='left',
            label='Dark mode',
            value=False,
            style={'display': 'flex', 'align-items': 'left', 'justify-content': 'left'},
            vertical=False
            ),
        ])

@app.callback(
    Output('opt-dropdown', 'options'),
    [Input('name-dropdown', 'value')]
)

@app.callback (Output('graph','figure'),[Input('selected-type','value'),Input('name-dropdown','value'),Input('mode-switch','value')])
def update_graph(analysis,feature,mode):
    themes = ["plotly", "plotly_white", "plotly_dark", "ggplot2", "seaborn", "simple_white", "none"]
    if mode == False:
        display_mode="simple_white"
    else:
        display_mode="plotly_dark"

    if analysis == 'pca':
        fig = px.scatter(dfmerge_pca,x="X",
                        y="Y",
                        opacity=0.7,
                        color=feature,
                        hover_name="samples",
                        height=600,
                        template=display_mode,
                        )
        if display_mode == 'plotly_dark':
            fig.update_xaxes(showgrid=False,showline=True,zeroline=False)
            fig.update_yaxes(showgrid=False,showline=True,zeroline=False)


    elif analysis == 'umap':
        fig = px.scatter(dfmerge_umap,x="X",
                y="Y",
                opacity=0.7,
                color=feature,
                hover_name="samples",
                height=600,
                template=display_mode,
                )
        if display_mode == 'plotly_dark':
            fig.update_xaxes(showgrid=False,showline=True,zeroline=False)
            fig.update_yaxes(showgrid=False,showline=True,zeroline=False)
    return fig



if __name__ == '__main__':
    app.run_server()
