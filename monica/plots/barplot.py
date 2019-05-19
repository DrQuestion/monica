import os
import numpy as np
import plotly.offline as py
import plotly.graph_objs as go
import matplotlib.pyplot as mplt


OPACITY = 0.75


# TODO: add host species '_(host)' suffix when host present by using an host parameter default None
# TODO: add optionality by interactive button to switch to CVD-friendly color palette (purples?)
# TODO: represent on another barplot aside the unmapped material in bases sequenced colored in grey
# TODO: add hover table showing information on the selected taxunit, like score distribution among strains and barcode
# TODO: (fix by clicking + highlighting the selected bar?)


def _by_taxunit(dataframe):
    by_taxunit = dict()
    for index, series in dataframe.iterrows():
        taxunit = index[0]
        if taxunit not in by_taxunit:
            by_taxunit[taxunit] = series.fillna(value=0).values
        else:
            by_taxunit[taxunit] += series.fillna(value=0).values
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


def plotter(alignment_df, output_folder=None, palette='jet',
            host=None, guests=None, mode=None, show_legend=True, auto_open=True):
    x = [col for col in alignment_df.columns]
    alignment_by_taxunit = _by_taxunit(alignment_df)
    bars = []
    if host:
        host = '_'.join(host.split(sep=' '))
        if guests:
            print(guests)
            # guests = [g for g in map(lambda guest: '_'.join(guest.split(sep=' ')), guests)]
            # guests = [g for g in map(lambda guest: '_'.join(guest.split(sep=' ')), list(guests))]
            guests = [g for g in guests]
            title_text = 'Guests: {}; host: {}; analysis mode: {}'.format(guests, host, mode)
        else:
            # Unlikely? Overnight with host given?
            title_text = 'Host: {}; analysis mode: {}'.format(host, mode)
    elif guests:
        print(guests)
        guests = [g for g in guests]
        # guests = [g for g in map(lambda guest: '_'.join(guest.split(sep=' ')), list(guests))]
        print(guests)
        title_text = 'Guests: {}; analysis mode: {}'.format(guests, mode)
    else:
        # Overnight no host given?
        title_text = 'Analysis mode: {}'.format(mode)
    colors_spaces, colors_lines = color_generator(len(alignment_by_taxunit), palette)
    for taxunit, color_space, color_line in zip(alignment_by_taxunit.keys(), colors_spaces, colors_lines):
        if taxunit == host:
            name = taxunit + '_(host)'
        else:
            name = taxunit
        bars.append(go.Bar(
            x=x,
            y=list(alignment_by_taxunit[taxunit]),
            name=name,
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
        title=go.layout.Title(
            font=dict(
                size=24,
                color='#436eee',
                family='droid serif, bold'
            ),
            text=title_text,
            xref='paper',
            x=0
        ),
        barmode='stack',
        hovermode='closest',
        showlegend=show_legend,
        legend=dict(
            orientation='h'
        )
    )

    fig = go.Figure(data=bars, layout=layout)
    py.plot(fig, filename=os.path.join(output_folder, 'monica.barplot.html'), auto_open=auto_open)
