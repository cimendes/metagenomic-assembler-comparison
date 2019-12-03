#!/usr/bin/env python3
"""

This script takes csv tables by species and plots them as series of scatter plots

Edited by the King of Plots through dark arts and the use of some kind of putrid concoctions.
"""

import sys
import os

import pandas as pd
from plotly import subplots
import plotly.graph_objs as go
from plotly.offline import plot
import plotly.figure_factory as ff


def get_file_names(files):
    result = []
    for file in files:
        result.append(os.path.basename(file).split('.')[0].replace('_', ' '))
    return result


def main():

    try:
        csv_tables = sys.argv[1:]

    # If no arguments were given, just print the header line.
    except IndexError:
        print('No files provided. Exiting...')
        sys.exit(0)

    num_cols = 3
    num_rows = len(csv_tables)/num_cols
    num_rows = int(num_rows) if len(csv_tables)%num_cols == 0 else int(num_rows)+1

    # define number and organization of subplots
    meta_subplots = subplots.make_subplots(rows=num_rows,
                                           cols=num_cols,
                                           shared_yaxes=True,
                                           subplot_titles=get_file_names(csv_tables),
                                           horizontal_spacing=0.04,
                                           vertical_spacing=0.08)

    r = 1
    c = 1
    # colors for each subplot
    colors = ['#abd9e9', '#80cdc1', '#d6604d', '#4575b4', '#006837', '#fee090',
              '#bababa', '#b2abd2', '#66bd63', '#a50026', '#d73027', '#fdae61']

    # call Cthulhu and beg him to make this work
    for file in csv_tables:
        filename = os.path.basename(file).split('.')[0].replace('_', ' ')

        print('Processing... {0}'.format(filename))

        # import table with data
        data = pd.read_csv(file)
        coverage = list(data['Breadth of coverage'])
        contigs = list(data['Contigs'])
        assemblers = list(data['Assembler'])

        # create a tracer for each assembler data point
        # simpler to manage colors and legends
        tracers = []
        legend = True if r == 1 and c == 1 else False
        for a in range(len(assemblers)):
            tracer = go.Scatter(x=[contigs[a]],
                                y=[coverage[a]],
                                name=assemblers[a],
                                showlegend=legend,
                                legendgroup='group{0}'.format(a),  # group legends
                                opacity=1,
                                mode='markers',
                                marker=dict(color=colors[a],
                                            size=7,
                                            line=dict(width=1,color='black'))
                                )

            tracers.append(tracer)

        # add assemblers tracers to subplot
        for t in tracers:
            meta_subplots.add_trace(t, r, c)

        # define xaxes attributes
        meta_subplots.update_xaxes(showgrid=False,
                                   showline=True,
                                   linecolor='black',
                                   linewidth=2,
                                   row=r,
                                   col=c,
                                   ticks='outside',
                                   tickcolor='black',
                                   tickwidth=2,
                                   tickfont=dict(color='black',
                                                 size=12),
                                   range=[0, 1200],
                                   tickvals=[0, 300, 600, 900, 1200])

        # add axes title to middle column/bottom row subplot
        if c == 2 and r == num_rows:
            meta_subplots.update_xaxes(title_text='Number of Contigs',
                                       titlefont=dict(size=18, color='black'),
                                       row=r,
                                       col=c)

        # define yaxes attributes
        meta_subplots.update_yaxes(showgrid=False,
                                   showline=True,
                                   linecolor='black',
                                   linewidth=2,
                                   row=r,
                                   col=c,
                                   tickfont=dict(color='black',
                                                 size=12),
                                   range=[0, 1.05],
                                   tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1])

        # only add ticks to first column subplots
        if c == 1:
            meta_subplots.update_yaxes(ticks='outside',
                                       tickwidth=2,
                                       tickcolor='black',
                                       row=r,
                                       col=c)
            # add axes title to first column/middle row subplot
            if r == 2:
                meta_subplots.update_yaxes(title_text='Breadth of Coverage (%)',
                                           titlefont=dict(size=18, color='black'),
                                           row=r,
                                           col=c)

        # define the next subplot that will be populated with tracer data
        c += 1
        if c > num_cols:
            r += 1
            c = 1

    # change background color for all subplots
    meta_subplots['layout']['plot_bgcolor'] = 'rgb(255,255,255)'

    # change legend attributes
    meta_subplots['layout']['legend']['font'] = dict(color='black', size=12)
    #meta_subplots['layout']['legend']['orientation'] = 'v'
    #meta_subplots['layout']['legend']['x'] = 0.80
    #meta_subplots['layout']['legend']['y'] = 0
    meta_subplots['layout']['legend']['borderwidth'] = 2

    # change annotations attributes
    for i in meta_subplots['layout']['annotations']:
        i['font']['color'] = 'black'
        i['font']['size'] = 16

    # it worked and the end result looks surprisingly ok...?
    # Praise Lord Cthulhu
    plot(meta_subplots, filename='scandalous_plots.html', auto_open=True)


if __name__ == '__main__':

    # do the beep boop dance!
    main()
