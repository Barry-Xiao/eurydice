import numpy as np
import pandas as pd
import skbio

from matplotlib import rcParams
import matplotlib.colors as mpc
import matplotlib.pyplot as plt
import seaborn as sn

rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42

hide_axes = dict(left=False, 
                 right=False, 
                 labelleft=False, 
                 labelright=False, 
                 length=0, labelsize=0)

axis_in = dict(left=True, 
               right=True, 
               labelleft=False,
               labelright=False,
               tickdir='in',
               labelsize=0)


def build_pcoa_plot(figpc: plt.Figure, 
                    pc0:str, pc1:str, hue:str, 
                    data: pd.DataFrame,
                    means: pd.DataFrame = None,
                    palette: dict=None, 
                    hue_order: list=None, 
                    pc_ratio: int=4, 
                    perc_exp: pd.Series=None, 
                    ax_lab: list=[1, 2], 
                    margin: str='box', 
                    margin_rolling: float=0.05, 
                    scatter_kws: dict=None, 
                    mean_kws: dict=None, 
                    box_kws: dict=None, 
                    trace_kws: dict=None, 
                    ylabel_kws: dict=None):
    """
    Builds Ordination plot with marginal axes

    Parameters
    ----------
    figpc : Figure
        Matplotlib figure; assumed to be square
    pc0, pc1: str
        Columns in data containing the x and y coordinates, respectively. the
        columsn are expected to be numeric, the names shoudl be strings.
    hue: str, optional
        Columns in data that should be used to select the colors hue for the 
        plot. This can be continous, or categorical
    data: pd.DataFrame
        A dataframe containing the information for all samples include the 
        PC coordinates (given in `pc0` and `pc`) and the `hue`.
        This is likely the output a merge between a metadata file and 
        the samples information from an skbio.OrdinationResults object.
    means: DataFrame, optional
        A collapsed version of data that groups in the data into key centroid
        points, however centroids are calculated. Should contain `pc0`, `pc1`,
        and `hue`, if applicable
    palette: str, dict, optional
        A matplotlib compatiable color palette whcih should be used to color 
        the data.
    hue_order : ndarray like, optional
        The order for the `hue` groups
    pc_ratio: int, optional
        The ratio between the main PC axis and the marginal axis.
    perc_exp : pd.Series, optional
        The percent explained by each eaxis. Likely the proportion_explained 
        parameter from a  skbio.OrdinationResults object. If not supplied, 
        the axes will not be labeled with the percent expalined. 
    box_kws, scatter_kws, trace_kws, ylabel_kws: dict, optional
        Formatting arguments the bxoplots and scatter plots, respectively.
        Arguments here will not be passed to the mean scatter plot. See
        sn.scatterplot and sn.boxplot for more details nad options.
    """
    # ### AXES ###
    # Sets up the gridspec and axes
    gspc = figpc.add_gridspec(2, 2, 
                              height_ratios=[1, pc_ratio], 
                              width_ratios=[pc_ratio, 1])
    axm = figpc.add_subplot(gspc[1, 0])
    ax0 = figpc.add_subplot(gspc[0, 0], sharex=axm)
    ax1 = figpc.add_subplot(gspc[1, 1], sharey=axm)

    # ### SCATTER ####
    # Plots the main scatter plot
    _build_pcoa_scatter(axm, data, pc0, pc1, hue, palette, hue_order, 
                        scatter_kws)
    
    # If desired, plots the mean or highlight scatter
    if means is not None:
        _build_pcoa_scatter(axm, means, pc0, pc1, hue, palette, hue_order,
                            mean_kws)
    # Cleans up the axes
    _tidy_main_ax(axm, pc0, pc1, perc_exp, ax_lab, ylabel_kws)

    # ### MARGINAL AXES ###
    if margin == 'box':
        _build_pcoa_marginal_boxes(ax0, ax1, data, pc0, pc1, hue, 
                                   palette, hue_order, box_kws)
    elif margin == 'trace':
        _build_pcoa_marginal_traces(ax0, ax1, data, pc0, pc1, hue, 
                                    margin_rolling, trace_kws)
    # Cleans up the marginal axes
    _tidy_marginal(ax0, ax1)
    
    

def _build_pcoa_scatter(ax_m, data, pc0='pc0', pc1='pc1', hue=None,
                        palette=None, hue_order=None, scatter_kws=None):
    """
    Builds the scatter plot portion of a marginal ordination plot

    Parameters
    ----------
    ax_m: plt.Axes
        Axis where the data should be plotted
    data: pd.DataFrame
        A dataframe containing the information for all samples include the 
        PC coordinates (given in `pc0` and `pc`), the hue, if applicable, and 
        the style, if applciable. This is likely the output a merge between
        a metadata file and the samples information from an 
        skbio.OrdinationResults object. 
       pc0, pc1: str
        Columns in data containing the x and y coordinates, respectively. the
        columsn are expected to be numeric, the names shoudl be strings.
    hue: str, optional
        Columns in data that should be used to selec the colors (hue) and 
        marker style (style) for the plot. If omitted, thse are arguments 
        will be ignored and everything will be one color or one style.
    palette: str, dict, optional
        A matplotlib compatiable color palette whcih should be used to color 
        the data.
    hue_order : ndarray like, optional
        The order for the `hue` groups
    
    Additional kwargs passed to plt.scatter_plot
    """
    
    # Sets up the default for the scatter parameters
    scatter_kwargs =  dict(s=2, 
                           edgecolor="None", 
                           alpha=0.5)
    if scatter_kws is not None:
        scatter_kwargs.update(scatter_kws)
    scatter_kwargs['legend'] = False

    # Builds the scatter plot
    sn.scatterplot(ax=ax_m, x=pc0, y=pc1, 
                   data=data,
                   palette=palette,
                   hue=hue,
                   hue_order=hue_order,
                   **scatter_kwargs)


def _build_pcoa_marginal_boxes(ax0, ax1, data, pc0, pc1, hue, 
                               palette=None, hue_order=None, box_kws=None):
    """
    Plots the marginal boxplots for categorical data
    """
    # Sets up the boxplot parameters
    box_kwargs = dict(linewidth=1, fliersize=1, width=0.75)
    if box_kws is not None:
        box_kwargs.update(box_kws)
    box_kwargs['order'] = hue_order
    box_kwargs['palette'] = palette

    sn.boxplot(x=pc0, y=hue, ax=ax0, data=data, orient='h', **box_kwargs)
    sn.boxplot(y=pc1, x=hue, ax=ax1, data=data, orient='v', **box_kwargs)


def _build_pcoa_marginal_traces(ax0, ax1, data, pc0, pc1, hue, rolling=0.05, 
                                trace_kws=None):
    """
    Continous trace of values along the axis of an ordination plot
    """
    # Sets up the default line kwargs
    line_kwargs = dict(linewidth=1, color='k')
    if trace_kws is not None:
        line_kwargs.update(trace_kws)
    
    # Calculates the smoothing factor
    if rolling < 1:
        smooth = int(np.round(rolling * len(data), 0))
    else:
        smooth = rolling

    data[hue] = data[hue].astype(float)
    
    # Calculates the rolling mean for the continous value along the data
    x0 = data.sort_values(pc0, ascending=True)[pc0]
    y0 = data.sort_values(pc0, ascending=True)[hue].rolling(smooth).mean()

    y1 = data.sort_values(pc1, ascending=True)[pc1]
    x1 = data.sort_values(pc1, ascending=True)[hue].rolling(smooth).mean()

    ax0.plot(x0, y0, **line_kwargs)
    ax1.plot(x1, y1, **line_kwargs)
    
def _tidy_marginal(ax0, ax1):
    """
    Cleans up the PC axes
    """
    # Clean sup the marginal axes
    marg_lim = np.array([
        np.min(np.hstack([ax0.get_ylim(), ax1.get_xlim()])),
        np.max(np.hstack([ax0.get_ylim(), ax1.get_xlim()]))
        ])
    ax0.set_ylim(marg_lim)
    ax1.set_xlim(marg_lim)

    for ax_ in [ax0.xaxis, ax0.yaxis, ax1.xaxis, ax1.yaxis]:
        ax_.set_tick_params(**hide_axes)
        ax_.get_label().set_visible(False)
    for ax in [ax0, ax1]:
        # Hides a legend, if there is one
        if ax.get_legend() is not None:
            ax.get_legend().set_visible(False)
        # Hides the axis facecolor
        ax.set_facecolor('None')
        # Despines the axis
        sn.despine(ax=ax, left=True, right=True, top=True, bottom=True)


def _tidy_main_ax(axm, pc0, pc1, perc_exp=None, pc_labels=[1, 2], kwargs=None):
    """
    Tidies the main axis for the PCoA plot
    """
    # Sets up defaults
    text_kwargs = dict(size=10)
    if kwargs is not None:
        text_kwargs.update(**kwargs)
    
    # Formats the axes
    axm.xaxis.set_tick_params(**axis_in)
    axm.yaxis.set_tick_params(**axis_in)

    # Formats the xlabel
    if perc_exp is None:
        x_lab = 'PC {0:1.0f}'.format(*pc_labels)
        y_lab = 'PC {1:1.0f}'.format(*pc_labels)
        perc_exp = dict()
    else:
        x_lab = f'PC {{0:1.0f}} ({{{{{pc0}:>4.1%}}}})'.format(*pc_labels)
        y_lab = f'PC {{1:1.0f}} ({{{{{pc1}:>4.1%}}}})'.format(*pc_labels)

    axm.set_xlabel(x_lab.format(**perc_exp), **text_kwargs)
    axm.set_ylabel(y_lab.format(**perc_exp), **text_kwargs)


def _check_kwargs(*args):
    """
    Checks the kwargs are not None, and provides an empty dict if they are
    """
    arg_new = [{None: dict}.get(v, v) for v in args]
    
    return arg_new


def build_categorical_legend(figl, order, palette, n_rows=5,
                              label_lookup=None,
                              scatter_kws=None, text_kws=None):
    """
    Makes a categorical legend figure    
    """
    scatter_kwargs = dict(marker='o', s=50, edgecolor='None')
    if scatter_kws is not None:
        scatter_kwargs.update(scatter_kws)

    text_kwargs = dict(ha='left', va='center', size=8)
    if text_kws is not None:
        text_kwargs.update(text_kws)

    if label_lookup is not None:
        labels = [label_lookup.get(o, o) for o in order]
    else:
        labels = order
    
    if isinstance(palette, str):
        palette =  sn.color_palette(palette, n_colors=len(order))
        colors = {v: to_hex(c) for v, c in zip(*(order, palette))}
    elif isinstance(palette, (list, tuple, np.ndarray)):
        colors = {v: c for v, c in zip(*(order, palette))}
    elif isinstance(palette, dict):
        colors = palette
    else:
        raise ValueError('The color palette must be a dictionary, list, or '
                         'color palette string.')
    colors = pd.Series(colors)[order]
    p = np.arange(len(order))
    x = np.floor(p / n_rows)
    y = p - x * n_rows
    
    axl = figl.add_subplot(1,1,1)
    axl.scatter(x=x, y=y, c=colors.values, **scatter_kwargs)

    for xx, yy, ll in zip(*(x, y, labels)):
        axl.text(x=xx + 0.1, y=yy, s=ll, **text_kwargs)
    axl.set_ylim(y.max() + 0.5, -0.7)
    axl.set_xlim(-0.1, x.max() + 1.05)
    axl.xaxis.set_tick_params(bottom=False, labelbottom=False, top=False, 
                              labeltop=False, length=0, labelsize=0)
    axl.yaxis.set_tick_params(bottom=False, labelbottom=False, top=False, 
                              labeltop=False, length=0, labelsize=0)

def _plot_pcoa_continous(fig, x, y, hue, data, palette='viridis', 
                         perc_exp=None, hue_scale=None,
                         palette_min=None, palette_max=None,
                         smooth=15, scatter_kws=None, line_kws=None):
    """
    Plots a PCoA with a continous marginal axis showing the distribution
    of data
    """
    # Sets up the colormap
    data[hue] = data[hue].astype(float)
    palette2 = _build_continous_cmap(hue, palette, data, hue_scale, 
                                     min_=palette_min, max_=palette_max)
    grouper = data.index.names 
    data['color'] = pd.concat(
         axis=0, 
         objs=[(data[hue] >= thresh) * 1 for thresh in palette2.keys()]
    ).groupby(grouper).sum() - 1
    data['color'] = data['color'] * hue_scale + min(palette2.keys())

    # Sets up the default kwargs
    scatter_kwargs = dict(legend=False, s=5, edgecolor='None')
    if scatter_kws is not None:
        scatter_kwargs.update(scatter_kws)
    scatter_kwargs['palette'] = palette2
    scatter_kwargs['legend'] = False

    line_kwargs = dict(linewidth=1, color='k')
    if line_kws is not None:
        line_kwargs.update(line_kws)
                                     
    # Sets up the axes for plotting
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 4], width_ratios=[4, 1])
    axm = fig.add_subplot(gs[1, 0])
    ax0 = fig.add_subplot(gs[0, 0], sharex=axm, facecolor='None')
    ax1 = fig.add_subplot(gs[1, 1], sharey=axm, facecolor='None')

    # Plots the scatter and marginal traces
    sn.scatterplot(ax=axm, x=x, y=y, data=data, 
                   hue='color', 
                   **scatter_kwargs
                  )
    ax0.plot(data.sort_values([x])[x],
             data.sort_values([x])[hue].rolling(smooth).mean(),
             **line_kwargs) 
    ax1.plot(data.sort_values([y])[hue].rolling(smooth).mean(),
             data.sort_values([y])[y],
             **line_kwargs)

    # Cleans up the axes
    for ax_ in [ax0.xaxis, ax0.yaxis, ax1.xaxis, ax1.yaxis]:
        ax_.set_tick_params(bottom=False, labelbottom=False, 
                            top=False, labeltop=False, 
                            length=0, 
                            labelsize=0)
        ax_.get_label().set_visible(False)
    ax1.set_xlim(ax0.get_ylim())
    axm.xaxis.set_tick_params(left=True, right=True, labelleft=False, 
                              labelright=False, tickdir='in')
    axm.yaxis.set_tick_params(left=True, right=True, labelleft=False, 
                              labelright=False, tickdir='in')
    ax1.set_xlim(np.array(ax1.get_xlim()) + 
                 np.diff(ax1.set_xlim()) * np.array([-0.2, 0.2]))
    ax0.set_ylim(ax1.get_xlim())
    sn.despine(ax=ax0, left=True, right=True, top=True, bottom=True)
    sn.despine(ax=ax1, left=True, right=True, top=True, bottom=True)

    # Adds axis label     
    if perc_exp is not None:
        axm.set_xlabel('PC 1 ({0:>4.1%})'.format(*perc_exp), size=10)
        axm.set_ylabel('PC 2 ({1:>4.1%})'.format(*perc_exp), size=10)
    else:
        axm.set_xlabel('PC 1')
        axm.set_ylabel('PC 2')

    data.drop(columns=['color'], inplace=True)

def build_continous_legend(figl, hue, data, palette='viridis', hue_scale=None, 
                           num_groups=10,  palette_min=None, 
                           palette_max=None, num_labels=5, 
                           tickformat='{x:1.0f}', text_kws=None):
    """
    Makes a legend figure for a PCoA with continous data
    """
    # Sets up the text arguments
    text_kwargs = dict(ha='center', va='bottom', size=8)
    if text_kws is not None:
        text_kwargs.update(text_kws)
    
    # Gets the color palette
    palette2 = \
        _build_continous_cmap(hue, palette, data, hue_scale, num_groups,
                              palette_min, palette_max)
    palette2 = pd.Series(palette2)
    
    # Determines the linespace
    num_groups =  len(palette2)
    bound = 1 / (num_groups)
    left = np.arange(0, 1, bound)
    center = left + bound / 2

    # Plots the parts
    ax = figl.add_subplot(1,1,1)
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.05, 1.05)

    ax.barh(y=np.ones(left.shape) * 0.375,
             left=left,
             height=0.25,
             width=np.ones(left.shape) * bound,
             color=palette2.values
            )

    # Set sup the labels
    position = np.linspace(center.min(), center.max(), num_labels)
    values = np.linspace(palette2.index.min(), palette2.index.max(), 
                         num_labels)
    for x, l in zip(*(position, values)):
        ax.text(x, 0.625, tickformat.format(x=l), **text_kwargs)

    # Clanes up the axis
    ax.yaxis.set_tick_params(left=False, right=False, labelleft=False, 
                             labelright=False, labelsize=0, length=0)
    ax.xaxis.set_tick_params(left=False, right=False, labelleft=False, 
                             labelright=False, labelsize=0, length=0)
