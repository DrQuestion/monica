import pickle

import pandas as pd
import numpy as np
import plotly.offline as py
import plotly.graph_objs as go
import matplotlib.pyplot as mplt

from monica.genomes.aligner import ALIGNMENT_PICKLE, normalizer, alignment_to_data_frame


ALIGNMENT = normalizer(pickle.load(open(ALIGNMENT_PICKLE, 'rb'))) #if wish to change to not normalized have to switch here, but wouldn't make sense
DATAFRAME = alignment_to_data_frame(ALIGNMENT)
OPACITY = 0.75


# TODO: add host species '_(host)' suffix when host present by using an host parameter default None
# TODO: add optionality by interactive button to switch to CVD-friendly color palette (purples?)
# TODO: represent on another barplot aside the unmapped material in bases sequenced colored in grey
# TODO: add hover table showing information on the selected taxunit, like score distribution among strains and barcode
# TODO: (fix by clicking + highlighting the selected bar?)


def dataframe_by_taxunit(dataframe=DATAFRAME):
    by_taxunit=dict()
    for index, series in dataframe.iterrows():
        taxunit=index[0]
        if taxunit not in by_taxunit:
            by_taxunit[taxunit]=series.fillna(value=0).values
        else:
            by_taxunit[taxunit]+=series.fillna(value=0).values
    return by_taxunit


def color_generator(n_elements, palette):
    cmap = mplt.get_cmap(name=palette)
    colors = cmap(np.linspace(0, 1, n_elements))
    lines = list(map(tuple, colors))
    lines = list(map(lambda x: 'rgba'+str(x), lines))
    spaces = colors
    for i in range(n_elements):
        spaces[i][3] = OPACITY
    spaces = list(map(tuple, spaces))
    spaces = list(map(lambda x: 'rgba'+str(x), spaces))
    return spaces, lines


def plotter(alignment, palette='jet'):
    x = list(ALIGNMENT.keys())
    bars = []
    colors_spaces, colors_lines = color_generator(len(alignment), palette)
    for taxunit, color_space, color_line in zip(alignment.keys(), colors_spaces, colors_lines):
        bars.append(go.Bar(
            x=x,
            y=list(alignment[taxunit]),
            name=taxunit,
            marker=dict(
                color=color_space,
                line=dict(
                    color=color_line,
                    width=1.5
                )
            ),
            opacity=0.75,
            hovertext=taxunit,
            hoverinfo='y+text'
        ))

    layout = go.Layout(
        barmode='stack',
        hovermode='closest',
        legend=dict(
            orientation='h'
        )
    )

    fig = go.Figure(data=bars, layout=layout)
    py.plot(fig, filename='barplot')


if __name__ == '__main__':

    df=dataframe_by_taxunit()
    plotter(df)
