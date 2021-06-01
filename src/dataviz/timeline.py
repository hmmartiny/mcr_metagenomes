#!/usr/bin/env python3
import matplotlib.pyplot as plt
import pandas as pd

def make_timeline(ax, dates, labels):
    """Create a timeline on an axis

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The ax to plot the timeline on.
    dates : list
        List of dates marking the point in time.
    labels : list
        Labels to plotted for the time points.
    """
    
    levels = np.repeat(1, len(dates))
    
    markerline, _, _ = ax.stem(
        dates, levels,
        #linefmt="C3-",
        basefmt='k-',
        use_line_collection=True
    )
    
    markerline.set_ydata(np.zeros(len(dates)))

    vert = np.array(['top', 'bottom'])[(levels >0).astype(int)]
    for d, l, r, va in zip(dates, levels, labels, vert):
        ax.annotate(
            r, xy=(d, l), xytext=(-3, np.sign(l)*2),
            textcoords="offset points",
            va=va,
            ha='center'
        )

    # format axis
    ax.get_yaxis().set_visible(False)
    for spine in ["left", "top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.margins(y=0.1)

def make_bars(ax, data, x, ys, width, colors, alpha=.5):
    """Make side by side bars on an axis

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to plot the bars on.
    data : pandas.DataFrame
        DataFrame with data for x and y
    x : str
        Name of column with values on x-axis
    ys : list of strings
        List of names of column containg values to be used as heights for the bars
    width : float
        Value for how much to shift each bar away from the previous one.
    colors : List of strings
        List of colors for each type of bar
    alpha : float, optional
        Opacity of bars, by default .5
    """

    shift=0
    for y, c in zip(ys, colors):
        ax.bar(
            data[x] + shift,
            data[y],
            width = width,
            color=c,
            label=y,
            alpha=alpha
        )
        
        shift += width

def make_stacked_bars(ax, data, labels, index, width):
    """Make stacked barplot

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axis to create the stacked bars on
    data : pd.DataFrame
        DataFrame with data to be plotted
    labels : list
        [description]
    index : list
        [description]
    width : float
        Width of bar
    """
    
    ax.bar(index, data[labels[0]], width, label=labels[0])
    bottom = np.array(data[labels[0]])
    
    i = 1
    while i < len(labels):
        ax.bar(index, data[labels[i]], width, label=labels[i], bottom=bottom)
        bottom += np.array(data[labels[i]])
        i += 1

def plot_timeline(cm, cm_closed, gene_timeline, timesamples, samples_ys, figsize=(12, 7), legend_kwargs={}, ylabel_kwargs={}):
    
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    gs = fig.add_gridspec(ncols=1, nrows=3, height_ratios=[.15, 1, .5])

    ax_timeline = fig.add_subplot(gs[0])
    ax_alr = fig.add_subplot(gs[1], sharex=ax_timeline)
    ax_samples = ax_alr.twinx()
    ax_bar = fig.add_subplot(gs[2], sharex=ax_timeline)
    
    # plot the timeline according to the litterature
    timeline = gene_timeline.groupby('Year discovered')['Gene'].apply(
        lambda x: "\n".join([gene_formatter(g) for g in x])
    )
    years = timeline.index.astype('int').tolist()
    labels = timeline.values.tolist()

    make_timeline(
        ax = ax_timeline,
        dates = years,
        labels = labels
    )

    # plot the count of samples with and without gene hits    
    make_bars(
        ax = ax_samples,
        data = timesamples,
        x = 'year',
        ys = samples_ys,
        width = .2,
        colors=['dimgray', 'lightgray']
    )
    
    # plot ALR data
    for gene in cm:
        ax_alr.plot(
            cm.index,
            cm[gene],
            '-o',
            label=gene_formatter(gene),
            lw = 1.5,
            alpha = .8
        )
        
    # plot relative abundance(s) in the composition
    make_stacked_bars(
        ax = ax_bar,
        data = cm_closed,
        labels = cm_closed.columns,
        index = cm_closed.index,
        width = .5
    )
    
    # show all xticks
    # TODO: If there are gaps in the year, the tick is not shown. Should be fixed.
    xticks = sorted(list(set(cm_closed.index.tolist() + years)))
    for ax in [ax_timeline, ax_alr, ax_bar]:
        ax.set_xticks(xticks)
    
    # hide axis on middle plot
    ax_alr.get_xaxis().set_visible(False)
    
    # create figure legends
    ax_alr.legend(**legend_kwargs)
    ax_samples.legend(loc='lower right', **legend_kwargs)
    
    # create labels for the y-axes
    ax_alr.set_ylabel('ALR', **ylabel_kwargs)
    ax_bar.set_ylabel('Relative abundance', **ylabel_kwargs)
    ax_samples.set_ylabel(r"$\log_{10}$(Count) (bars)", **ylabel_kwargs)
    
    # add title for discovery
    ax_timeline.text(xticks[0], .1, "Timeline of reporting", ylabel_kwargs)

    plt.close(fig)
    return fig