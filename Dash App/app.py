# Run the entire script and the console will give a link in which the App can be viewed

# Import packages
import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import plotly.express as px
import pandas as pd

# Load data
df = pd.read_csv("Combined_Dataset_thesis_final.csv")
df.rename(columns={'Network.size': 'Network'}, inplace = True)


# Create app
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

app.layout = html.Div(children=[
        html.H1(children='SARS-CoV-2 Transmission Dashboard', 
            style={'textAlign': 'center'}
            ),

    html.Div(children="Set sliders to make a selection of interest:", 
             style={'textAlign': 'center'}
             ),
    
       html.Div([
           dcc.Graph(id='fig1'),
       ]) ,
       html.Div([
           html.H6('R0'),
           dcc.Slider(
               id='slider-R',
               min=df['R'].min(),
               max=df['R'].max(),
               step=.5,
               value=df['R'].min(),
               marks={1: '1',
                      1.5: '1.5',
                      2: '2',
                      2.5: '2.5',
                      3: '3'}
           ),
           html.H6('Network size'),
           dcc.Slider(
               id='slider-Network',
               min=df['Network'].min(),
               max=df['Network'].max(),
               step=None,
               value=df['Network'].min(),
               marks={str(network): str(network) for network in df['Network'].unique()})
       ])
])

# Create color scheme
color_scheme = ["#544949", "#89ae90", "#efbc74", "#ef3222", "#b74242"]

@app.callback(
    Output('fig1', 'figure'),
    [Input('slider-R', 'value'), Input('slider-Network', 'value')])
def update_figure(selected_R, selected_Network):
    filtered_df = df[df.R == selected_R]
    filtered_df = filtered_df[filtered_df.Network == selected_Network]



    fig = px.box(filtered_df, x="Bubble.size", y="Percentage.infected",
                     color="PCR.frequency", template = "seaborn",
             color_discrete_sequence= color_scheme,
             category_orders={"Bubble.size": ["5", "10", "25", "50", "none"],
                              "PCR.frequency": ["workday", "semi-weekly",
                                                "weekly", "monthly", "none"]},
             labels={
                     "Bubble.size": "Bubble size",
                     "Percentage.infected": "% infected",
                     "PCR.frequency": "PCR frequency"
                 })

    fig.update_layout(transition_duration=500)

    return fig


if __name__ == '__main__':
    app.run_server(debug=True)