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
__version__ = '2.2.5'


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




parser = argparse.ArgumentParser()

parser.add_argument('-p', required=False, type=str, dest='polymorphic_sites',
                                help='provide the path of the polymorphic sites you want to analyze')
parser.add_argument('-conf', required=False, type=str, dest='config',
                                help='provide the path of the configuration file that contains the information to use for data visualization')
args= parser.parse_args()

df = pd.read_csv(args.polymorphic_sites, skiprows=4,sep='\t')
print (df.index)
print (list(df.columns))


# app.layout = html.Div(
#     [
#         html.Div([
#         dcc.Dropdown(
#             id='name-dropdown',
#             options=snptoolkitInputs,
#             value = snptoolkitInputs[0]['value']
#             ),
#             ],style={'width': '50%', 'display': 'inline-block'}),
#         html.Div([
#         dcc.Dropdown(
#             id='opt-dropdown',
#             ),
#             ],style={'width': '50%', 'display': 'inline-block'}
#         ),
#         html.Hr(),
#         dcc.RadioItems(id='selected-type',options=[
#             {'label': 'All SNPs', 'value': 'ALL'},
#             {'label': 'Synonymous SNPs', 'value': 'SYN'},
#             {'label': 'Non-synonymous SNPs', 'value': 'NS'},
#             ],value='ALL',labelStyle={'display': 'inline-block'}
#         ),
#         #TODO: dcc.Input(placeholder='Enter a window size to compute SNPs density...',type='number',id="windowsize",value=''),
#         dcc.Graph(id='graph'),
#         html.Hr(),
#         html.Div(dt.DataTable(columns = [{"name": i, "id": i} for i in header],
#         id='snpTable', 
#         filter_action="native",
#         page_size=20,
#         style_table={ 'align':'center','overflowY': 'scroll'},
#         style_cell={'textAlign': 'center', 'border': '1px solid grey'},
#         style_header={
#         'color':'white',
#         'backgroundColor': 'rgb(62, 155, 184)',
#         'border': '1px solid black'
#         },
#         style_data_conditional=[
#             {
#                 'if': {
#                 'column_id': 'Effect',
#                 'filter_query': '{Effect} eq "NS"'
#             },
#             'backgroundColor': '#FBBF4C',
#             'color': 'white',
#         },
#         {
#                 'if': {
#                 'column_id': 'Effect',
#                 'filter_query': '{Effect} eq "Syn"'
#             },
#             'backgroundColor': '#51A8C7',
#             'color': 'white',
#         },
#         {
#                 'if': {
#                 'column_id': 'Effect',
#                 'filter_query': '{Effect} eq "."'
#             },
#             'backgroundColor': '#0e2f2e',
#             'color': '#0e2f2e',
#         }
#         ])),
#     ]
# )


# if __name__ == '__main__':
#     app.run_server()