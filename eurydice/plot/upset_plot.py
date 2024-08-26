import itertools as it

import numpy as np
import pandas as pd

from matplotlib import colormaps
import matplotlib.colors as mpc
import matplotlib.pyplot as plt
import seaborn as sn

hide_axes = dict(left=False, right=False, labelleft=False, labelright=False, 
                 length=0, labelsize=0)

def plot_upset_part_only(figup, num_part, num_inter, poly_id, 
                         label_group_size=True, link_colors=None, 
                         spec_labels=None, labelsize=7,
                         unit='Participants\n($\\times 10^{3}$)',
                         off_text_part=1000, off_text_inter=500,
                         grouped_colors=None, group_size_kws=dict(),
                         overlap_scatter_kws=None, overlap_link_kws=None,
                         part_scale=1000):
    """
    Builds an upset plot based on the input matrix
    """
    # Set supt he axes for plotting
    # Sets up the gridspec and axes
    gs_up = figup.add_gridspec(2, 3,
                           height_ratios=[1, 1],
                           width_ratios=[1.5, 1, 2.5])

    # Sets up the overlaps
    ax_ov = figup.add_subplot(gs_up[-1, -1], facecolor='None')
    ax_ss = figup.add_subplot(gs_up[-1,  0], sharey=ax_ov, facecolor='None')
    ax_sx = figup.add_subplot(gs_up[-2, 0], sharex=ax_ss, facecolor='None')
    ax_sl = figup.add_subplot(gs_up[-1, 1], sharey=ax_ov, facecolor='None')
    ax_is = figup.add_subplot(gs_up[0,  -1], sharex=ax_ov, facecolor='None')
    ax_iy = figup.add_subplot(gs_up[0, -2], sharey=ax_is, facecolor='None')

    # Cleans up the axes and formats them. We'll re-format later, but we'll
    # start ehre
    for ax_ in figup.axes:
        ax_.xaxis.set_tick_params(**hide_axes)
        ax_.yaxis.set_tick_params(**hide_axes)
    sn.despine(figup, left=True, right=True, top=True, bottom=True)

    # Overlap matrix
    _plot_poly_id(ax_ov, poly_id, link_colors, overlap_scatter_kws, 
                  overlap_link_kws)

    # Number of participants bar chart
    _plot_upset_bar_h(ax_ss, num_part, num_part_colors=grouped_colors, 
                      color_axis=0, text=label_group_size, 
                      label_size=labelsize, text_off=off_text_part, 
                      **group_size_kws)
    # Numer of intersections bar chart
    _plot_upset_bar_v(ax_is, num_inter, num_inter_colors=grouped_colors, 
                      color_axis=0, text=label_group_size, 
                      label_size=labelsize, text_off=off_text_inter, 
                      **group_size_kws)

    # Handles label axes
    ss_max = \
        np.round(num_part.sum(axis=1).max() / part_scale) * part_scale + 1
    ss_xticks = \
        np.arange(0, ss_max, part_scale)
    ax_ss.set_xticks(ss_xticks)
    for x in ss_xticks:
        ax_sx.text(x=x, y=0.05, s=f'{x / part_scale:1.0f}', 
                   ha='center', va='bottom', size=8) 
    ax_sx.text(x=np.diff(ax_sx.get_xlim()) / 2 + ax_sx.get_xlim()[0],
               y=0.25, s=unit,
               size=8, ha='center', va='bottom',)
    sn.despine(ax=ax_sx, bottom=False, left=True, right=True, top=True, 
               trim=True)
    ax_sx.xaxis.set_tick_params(bottom=True, tickdir='in', length=3)

    iyt = \
        np.arange(0, 
                  np.round(num_inter.sum(axis=1).max() / part_scale) * part_scale + 1, 
                  part_scale)
    ax_iy.set_yticks(iyt)
    sn.despine(ax=ax_iy, left=True, right=False, bottom=True, top=True, trim=True)
    ax_iy.yaxis.set_tick_params(right=True, left=False, labelright=False, 
                                labelleft=False, tickdir='in', length=4)
    for y_ in iyt:
        ax_iy.text(0.85, y_, f'{y_ / part_scale:1.0f}', size=8, ha='right', 
                   va='center')
    ax_iy.text(x=0.4, 
               y=ax_iy.get_ylim()[0] + (np.diff(ax_iy.get_ylim()) / 2), 
               s=unit,  size=8, ha='center', va='center', rotation=90)

    # ### Overlap labels ####
    if spec_labels is None:
        spec_labels = dict()
    for y, (l, v) in enumerate(num_part.iterrows()):
         ax_sl.text(0.5, y, spec_labels.get(l, l), 
                    ha='center', va='center', size=8)
    ax_ov.set_ylim(ax_ov.get_ylim()[::-1])
    

def build_upset_matrix(overlapped: pd.DataFrame, 
                       intersect_cols: list=None,
                       group_order: list=None,
                       inter_order: (list, str)='counts',
                       count_group: str=None,
                       ):
    """
    Builds the linkage matrices for an upset plot 

    Parameters
    ----------
    overlapped: DataFrame
        A binary dataframe where each row represents an observation unit and 
        each column indicates an exposure group to overlap. 
    intersect_cols: array-like, optional
        The columns which should be used as the overlap. If no list is 
        supposed, all columns will be used.
    group_order: array-like, optional
        If the groups should be presented in some specific order, that 
        should be passed here. Otherwise, they will be ordered based on 
        group size
    inter_order:  {'counts', 'group_size', 'magic', list-like},
        A list of the exponentiated group labels to show in the intersection
        plot.
            - "counts": the number of observations in the overlapped group
            - "group_size": the number of groups in the overlap
            - "magic": hopefully pretty, need a better description
            - list_like: custom order
    
    Returns
    -------
    
    
    Raises
    ------
    ValueError
        inter_order is not an acceptable value 
    """
    # Sets hte default for the intersection
    if intersect_cols is None:
        intersect_cols = overlapped.columns
    if count_group is not None:
        intersect_check = \
            (overlapped.set_index(count_group)[intersect_cols].copy() > 0) * 1
    else:
        intersect_check = (overlapped[intersect_cols].copy() > 0) * 1

    # Calculates and sorts the number of groups overall
    num_group = intersect_check.sum(axis=0)
    if group_order is None:
        num_group.sort_values(ascending=False, inplace=True)
        group_order = list(num_group.index)
    else:
        num_group = num_group[group_order]

    # Calculates and sorts the intersections
    group_ids = np.arange(0, len(num_group)) + 1
    num_inter = \
        (intersect_check[group_order] * np.power(2, group_ids)).sum(axis=1)
    num_inter = num_inter.value_counts()

    # Calculates the poly_id binary matrix
    poly_id = _build_poly_id(group_order)
    poly_id = poly_id.loc[group_order, num_inter.index]

    # Gets the ordering for the binary matrix
    inter_order = \
        _build_ordering(poly_id, group_order, num_inter, inter_order)
    poly_id = poly_id.loc[group_order, inter_order].copy()
    num_inter = num_inter.loc[inter_order]
    
    # Groups the count plot
    if count_group is not None:
        num_group = intersect_check.groupby(count_group)[group_order].sum().T
        num_group = num_group.loc[group_order]
        num_inter = \
            (intersect_check * np.power(2, group_ids)).sum(axis=1).to_frame()
        num_inter['counter'] = 1
        num_inter.rename(columns={0: 'intersect'}, inplace=True)
        num_inter = \
            num_inter.groupby([count_group, 'intersect'])['counter'].sum()
        num_inter = \
            num_inter.unstack(0).fillna(0).loc[inter_order].astype(int)
        num_inter.index.set_names(None, inplace=True)
    else:
        num_inter = num_inter.to_frame().loc[inter_order]
        num_group = num_group.to_frame().loc[group_order]
    poly_id = poly_id.loc[group_order, inter_order]

    return num_group, num_inter, poly_id
    


def _build_ordering(poly_id, group_order, num_inter, inter_order):
    """
    Builds a list for the order of the data
    """
    ordering = poly_id.loc[group_order, num_inter.index].T.copy()
    ordering['group_size'] = ordering.sum(axis=1)
    ordering['counts'] = num_inter.values
    ordering['label'] = ordering.index
    
    if isinstance(inter_order, list):
        pass
    elif inter_order == 'counts':
        ordering.sort_values(['counts', 'label'], 
                             ascending=[False, True], inplace=True)
        inter_order = list(ordering.index.values)
    elif inter_order == 'group_size':
        ordering.sort_values(['group_size', 'label'],
                             ascending=[False, True], 
                             inplace=True)
        inter_order = list(ordering.index.values)
        
    elif inter_order == 'magic':
        # Gets the group label, which is built off the y index, for each of
        # of the columns and then sorts based on that grouping.
        order1 = pd.DataFrame.from_dict(
            orient='columns', 
            data={c: ordering[c] * ordering['label'] for c in group_order})
        order1.replace({0: order1.max().max() + 1}, inplace=True)
        order1.sort_values(group_order, ascending=True, inplace=True)
        inter_order = list(order1.index)
    else:
        raise ValueError('The inter_order must either be "counts", '
                         '"group_size", "magic", or a custom list')

    return inter_order


def _build_poly_id(group_order):
    """
    Expands data for a poly id matrix
    """
    group_ids =  np.arange(len(group_order)) + 1
    poly_id = {0: pd.Series([], dtype=int)}
    for n_j in np.arange(len(group_ids)) + 1:
        for i in it.combinations(group_ids, int(n_j)):
            poly_id[np.power(2, np.array(i)).sum()] = \
                pd.Series([True] * int(n_j),
                          index=[group_order[v - 1] for v in i])
    poly_id = pd.DataFrame(poly_id).fillna(False) * 1
    poly_id = poly_id.loc[group_order].copy()
    
    return poly_id


def _build_upset_coords(poly_id):
    """
    Puts together a 2D matrix for the upset plot and linkages  

    Parameters
    ----------
    poly_id: DataFrame

    Returns
    -------
    """
    
    poly_x = pd.DataFrame(
        data=(poly_id.values * 0 + 1) * np.arange(0, len(poly_id.columns)),
        index=poly_id.index,
        columns=poly_id.columns)
    poly_y = pd.DataFrame(
        data=(poly_id.values * 0 + 1) * \
            np.atleast_2d(np.arange(0, len(poly_x))).T,
        index=poly_id.index,
        columns=poly_id.columns)
    linkage_y = np.vstack([poly_y.mask(poly_id == 0).min(axis=0),
                           poly_y.mask(poly_id == 0).max(axis=0)])
    linkage_x = np.vstack([poly_x.iloc[0].values] * 2)

    return poly_x, poly_y, linkage_x, linkage_y


def _plot_poly_id(ax_ov, poly_id, link_colors=None, scatter_kws=None, 
                  link_kws=None):
    """
    Plots the overlap linkage matrix for an upset plot
    """
    # Gets the linkage colors
    if isinstance(link_colors, (list, tuple, np.ndarray)):
        pc1 = link_colors[1]
        pc0 = link_colors[0]
    else:
        pc1 = '#525252'
        pc0 = '#bfbfbf'
        
    # Sets up the keyword arguments for the scatter plot
    linkage_kwargs = dict(linewidth=4, zorder=0, color=pc1) 
    scatter_kwargs = dict(s=100, zorder=1, marker='o')
    if isinstance(scatter_kws, dict):
        scatter_kwargs.update(scatter_kws)
    if isinstance(link_kws, dict):
        linkage_kwargs.update(link_kws)

    # Gets the linkage values
    poly_x, poly_y, linkage_x, linkage_y = _build_upset_coords(poly_id)
    
    # #Makes the overlap scatter
    ax_ov.plot(linkage_x, linkage_y, **linkage_kwargs)
    ax_ov.scatter(x=np.hstack(poly_x.values), 
                  y=np.hstack(poly_y.values), 
                  c=np.hstack(poly_id.replace({0: pc0, 1: pc1}).values),
                  **scatter_kwargs
                  )

def _check_bar_colors(bar_group, color_group=None, color_dim=1):
    """
    Matches the dimension of the colors to the bar chart
    """
   
    if color_group is None:
        # If no value is supplied, we're going to use a gray color.
        colors = pd.DataFrame(data=np.zeros(bar_group.shape), 
                              index=bar_group.index,
                              columns=bar_group.columns)
        colors.replace({0: '#bdbdbd'}, inplace=True)
    elif isinstance(color_group, str) and mpc.is_color_like(color_group):
        # If we only have a color string, just build a matrix with that color
        colors = pd.DataFrame(data=np.zeros(bar_group.shape), 
                              index=bar_group.index,
                              columns=bar_group.columns)
        colors.replace({0: color_group}, inplace=True)
    elif (isinstance(color_group, str) and (color_group in list(colormaps)) 
              and (color_dim in {0, '0', 'columns'})):
        # If a colormap was passed and it should be displayed along columns,
        # we get the number of columns and build a stacked hex colormap. Note
        # that bar_dims gives the number of rows first. 
        bar_dim = bar_group.shape
        num_colors = bar_dim[1]
        color_list = np.array([
            mpc.to_hex(c) 
            for c in sn.color_palette(color_group, n_colors=num_colors)
            ])
        colors = pd.DataFrame(data=np.vstack([color_list] * bar_dim[0]),
                              index=bar_group.index,
                              columns=bar_group.columns)
    elif (isinstance(color_group, str) and (color_group in list(colormaps)) 
              and (color_dim in {1, '1', 'index'})):
        # If a colormap was passed and it should be displayed along columns,
        # we get the number of columns and build a stacked row colormap. Note
        # that bar_dims gives the number of rows first.
        bar_dim = bar_group.shape
        num_colors = bar_dim[0]
        color_list = np.array([
            mpc.to_hex(c) 
            for c in sn.color_palette(color_group, n_colors=num_colors)
            ])
        colors = pd.DataFrame(data=np.vstack([color_list] * bar_dim[1]).T,
                              index=bar_group.index,
                              columns=bar_group.columns)
    elif (isinstance(color_group, str) and (color_group in list(colormaps)) 
            and (color_dim in {-1, '-1', 'both'})):
        # If a colormap was passed and it should be displayed along both axes,
        # we get the number of columns and reshape into a column
        bar_dim = bar_group.shape
        num_colors = np.product(bar_dim)
        color_list = np.array([
            mpc.to_hex(c) 
            for c in sn.color_palette(color_group, n_colors=num_colors)
            ])
        colors = pd.DataFrame(data=color_list.reshape(bar_dim),
                              index=bar_group.index,
                              columns=bar_group.columns)
    elif isinstance(color_group, dict) and (color_dim in {0, '0', 'columns'}):
        bar_dim = bar_group.shape
        color_list = bar_group.columns.to_frame().replace(color_group).values
        color_list = color_list.flatten()
        colors = pd.DataFrame(data=np.vstack([color_list] * bar_dim[0]),
                              index=bar_group.index,
                              columns=bar_group.columns)
    elif isinstance(color_group, dict) and (color_dim in {1, '1', 'index'}):
        bar_dim = bar_group.shape
        color_list = bar_group.index.to_frame().replace(color_group).values
        color_list = color_list.flatten()
        colors = pd.DataFrame(data=np.vstack([color_list] * bar_dim[1]).T,
                              index=bar_group.index,
                              columns=bar_group.columns)
    elif (isinstance(color_group, pd.DataFrame) & 
              (color_group.shape == bar_group.shape)):
        # If the color group is already a matching palette, use that
        colors = color_group
    else:
        raise ValueError('The provided colors or color dimension are not '
                         'copesectic.\nPlease re-evaluate your life choices '
                         'and the function arguments')

    return colors


def _plot_upset_bar_h(ax_ss, num_part, num_part_colors=None, color_axis=0, 
                      text=True,  text_colors='k', text_off=50,
                      label_size=8, num_group_kws=dict()):
    """
    Plots the horizontal stacked barchart
    """
    # Gets the colors
    bar_colors = _check_bar_colors(num_part, num_part_colors, color_axis)
    text_colors = _check_bar_colors(num_part, text_colors, color_axis)

    # Plots the barplot
    y = np.arange(0, len(num_part))
    left_ = num_part.cumsum(axis=1) - num_part
    for c, v in num_part.items():
        ax_ss.barh(y=y, width=v, color=bar_colors[c], left=left_[c],
                   **num_group_kws)
    ax_ss.set_xlim(ax_ss.get_xlim()[::-1])
    ax_ss.set_yticks(y)

    # Adds group labels
    if text:
        spacer = 0.025 * np.absolute(np.diff(ax_ss.get_xlim()))
        for col_, vals in num_part.items():
            for y_, (i, v_) in zip(*(y, vals.items())):
                if v_ > 0:
                    x_ = (v_ + spacer) * (v_ < text_off) + \
                         (v_ > text_off) * (v_ / 2) + \
                         left_.loc[i, col_]
                    ha_ = {True: 'right', False: 'center'}[v_ < text_off]
                    ax_ss.text(x=x_, y=y_, s=v_, 
                               color=text_colors.loc[i, col_],
                               va='center',
                               ha=ha_,
                               size=label_size,
                               )


def _plot_upset_bar_v(ax_is, num_inter, num_inter_colors=None, color_axis=0,
                      text=True, text_colors='k', label_size=8, text_off=1000,
                      num_inter_kws=dict()):
    """
    Plots stacked barplot for interactions
    """

    # Gets the colors
    bar_colors = _check_bar_colors(num_inter, num_inter_colors, color_axis)
    text_colors = _check_bar_colors(num_inter, text_colors, color_axis)

    # Plots the barplot
    x = np.arange(0, len(num_inter))
    bottom_ = num_inter.cumsum(axis=1) - num_inter

    for c, v in num_inter.items():
        ax_is.bar(x=x, bottom=bottom_[c], height=v, color=bar_colors[c], 
                  **num_inter_kws)

    if text: 
        spacer = 0.05 * np.absolute(np.diff(ax_is.get_ylim()))
        for col_, vals in num_inter.items():
            for x_, (i, v_) in zip(*(x, vals.items())):
                if v_ > 0:
                    y_ = (v_ + spacer) * (v_ < text_off) + \
                         (v_ > text_off) * (v_ / 2) + \
                         bottom_.loc[i, col_]
                    ha_ = {True: 'bottom', False: 'center'}[v_ < text_off]
                    ax_is.text(x=x_, y=y_, s=v_, 
                               color=text_colors.loc[i, col_],
                               ha='center',
                               va=ha_,
                               rotation=90,
                               size=label_size,
                               )