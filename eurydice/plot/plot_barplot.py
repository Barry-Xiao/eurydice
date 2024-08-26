import biom
import pandas as pd 
import numpy as np

import scripts.barplot_hardcode as bhc

hide_axes = dict(left=False, 
                 right=False, 
                 labelleft=False, 
                 labelright=False, 
                 length=0, 
                 labelsize=0)

axis_in = dict(left=True, 
               right=True, 
               labelleft=False,
               labelright=False,
               tickdir='in',
               labelsize=0)

def collapse_taxa_2lvl(otu_q2, taxa_q2, keep_big=None, keep_small=None, 
                       big_group=1, small_group=4):
    """
    Collapses an OTU table to a taxonomic barplot using two taxonomic levels

    This function will collapse an OTU table into a nested taxonomic bar plot
    where taxa within the larger group/level will be collapsed into 
    level-specific taxonomic groups.

    Parameters
    ----------
    otu_q2 : QIIME 2 Artifact of FeatureTable[Frequency] semantic type
        The count x sample table for the samples of interest. If the data
        needs to be averaged, the data should be collapsed before its 
        passed into this function. Sample order does not matter.  
        (See the q2-feature-table plugin group function for collaping 
        by a metadata category)
    taxa_q2 : QIIME 2 Artifact of FeatureData[Taxonomy] semantic type
        An Artifact providing the taxonomic assignment for the features in 
        `otu_q2`
    keep_big: 1d ndarray-like
        The list of larger taxonomic level groups to keep. By default, this
        will default ot the hard coded phylum level assignments for ECHO 
        extant gut samples annotated with Silva 138.1.
    keep_small: 1d ndarray-like
        The list of the finer taxonomic level to keep unique. Anything in this
        list will be retained as a unique bar; other members of `keep_big`
        will be collapsed into a single group. All other taxa will be
        combined into an "other" category.
        By default, this will default ot the hard coded family level 
        assignments for ECHO extant gut samples annotated with Silva 138.1.
    large_group, small_group: int
        An integer designation for the taxonomic level for the nested 
        collapse: 1=phylum; 2=class; 3=order; 4=family; 5=genus. This should
        correspond to the taxonomy list supplied in keep_big and keep_small,
        respectively
    
    Returns
    -------
    pd.DataFrame
        A dataframe of the collapsed data with the samples as rows and 
        collapsed taxonomy as the columns.

    Also See
    --------
    collapse_taxa_1lvl
    """
    # Sets up the filtering function 
    if keep_big is None:
        keep_big = bhc.feces_keep_phyla
    if keep_small is None:
        keep_small = bhc.fecal_tax_order_hard
    # Pulls in level data
    little_str, big_str, other_str = _build_strings(big_group, small_group)
    big_lvl = f'tax_{big_group:1.0f}'

    # Defines the filter function
    def filter_fun(id_, md):
        fam_label = little_str.format(**md)
        if fam_label in keep_small:
            return fam_label
        elif md[big_lvl] in keep_big:
            return big_str.format(**md)
        else:
            return other_str

    collapse =  _collapse_otu_table(otu_q2, taxa_q2, filter_fun, keep_small)

    return collapse


def collapse_taxa_1lvl(otu_q2, taxa_q2, keep=None, group=5):
    """
    Collapses taxonomy based on a single group level

    This function will collapse an OTU table to a single taxonomic level,
    grouping anything not in a pre-designated set of taxonomic levels into
    an "other" cateogry at the same taxonomic level. 

    Parameters
    ----------
    otu_q2 : QIIME 2 Artifact of FeatureTable[Frequency] semantic type
        The count x sample table for the samples of interest. If the data
        needs to be averaged, the data should be collapsed before its 
        passed into this function. Sample order does not matter.  
        (See the q2-feature-table plugin group function for collaping 
        by a metadata category)
    taxa_q2 : QIIME 2 Artifact of FeatureData[Taxonomy] semantic type
        An Artifact providing the taxonomic assignment for the features in 
        `otu_q2`
    keep: 1d ndarray-like
        The list of larger taxonomic level groups to keep. By default, this
        will use the hard coded map for nasal samples in ECHO with Silva
        138.1 taxonomy.
    group: int
        An integer designation for the taxonomic level for the nested 
        collapse: 1=phylum; 2=class; 3=order; 4=family; 5=genus. This should
        correspond to the taxonomy list supplied in keep_big and keep_small,
        respectively. 
        By default, this is set to match the hard coded map for nasal 
        samples in ECHO with Silva 138.1 taxonomy.
    
    Returns
    -------
    pd.DataFrame
        A dataframe of the collapsed data with the samples as rows and 
        collapsed taxonomy as the columns.

    Also See
    --------
    collapse_taxa_2lvl
    """
    # Sets up the filtering function 
    if keep is None:
        keep = bhc.nasal_order
    # Pulls in level data
    little_str, _, other_str = _build_strings(1, group)

    # Defines the collapse function
    def filter_fun(id_, md):
        fam_label = little_str.format(**md)
        if fam_label in keep:
            return fam_label
        else:
            return other_str

    collapse =  _collapse_otu_table(otu_q2, taxa_q2, filter_fun, keep)

    return collapse


def construct_bar_position(meta, collapsed, sort_columns=None, 
                           space_columns=None, levels=None, 
                           step_scale=5):
    """
    Constructs a list of barplot positions

    Parameters
    ----------
    meta: pd.DataFrame
        Sample information, indexed by SampleID which matches the unit in 
        collapsed
    collapsed: pd.DataFrame
        Taxonomic dataframe with taxa as columns and samples as rows
    columns: list
        The columns in meta that should be used to order the data. We assume
        an asscending order for sorting, so code appropriately.
    levels: list
        The taxonomic levels which should be collapsed for sorting withi
        a group
    step_scale: numeric, dict
        The padding width to add between each column group. This can be 
        supplied as a constant (number) or as a dictionary if the size should
        vary. 

    Returns
    -------
    pd.DataFrame
        The sample order for the bar plot
    """
    # If no level is provided, we'll sort based on the first column
    if levels is None:
        levels = [collapsed.columns[0]]
    if sort_columns is None:
        sort_columns = ['ref']
    if space_columns is None:
        space_columns = ['ref']

    columns = np.unique(np.hstack([sort_columns, space_columns]))
    meta['ref'] = 1

    # Pulls out the taxonomic order
    tax_order = collapsed.index.to_frame()
    tax_order.rename({0: 'SampleID'}, inplace=True)

    # Updates the data with the metdata 
    tax_order[columns] = meta[columns]
    meta.drop(columns=['ref'], inplace=True)

    # Adds in the "best" order
    tax_order['best'] = collapsed[levels].sum(axis=1)

    # Gets the list of the columns to combine
    all_sort = list(np.hstack([columns, 'best']))
    tax_order.sort_values(all_sort, ascending=True, inplace=True)
    
    # Sorts by the order and gets the duplicated data    
    if isinstance(step_scale, (int, float)):
        step_scale = {col_: i * step_scale 
                       for i, col_ in enumerate(space_columns[::-1])}
    elif isinstance(step_scale, dict):
        step_scale['ref'] = 0
    else:
        raise ValueError('the step_scale must be an integer or dictionary '
                         'of column names')
        
    # Calculates the x position and intervals assuming each bar is 1 unit wide
    tax_order['x'] = pd.concat(
        axis=1,
        objs=[(~tax_order[col_].duplicated(keep='first')) * step_scale[col_]
               for col_ in space_columns]
    ).sum(axis=1) + 1
    tax_order['x'] = tax_order['x'].cumsum()
    tax_order['x'] = tax_order['x'] - tax_order['x'].min()

    return tax_order


def _collapse_otu_table(otu_q2, taxa_q2, filter_fun, keep):
    """
    Helper function to collapse OTU table to a taxonomic level
    """
    # Formats the OTU table for collapse
    otu_table = otu_q2.view(biom.Table)
    otu_table = otu_table.norm(axis='sample')
    
    # Adds the taxonomy to the OTU table
    taxa = \
        taxa_q2.view(pd.Series).loc[otu_table.ids(axis='observation')].copy()
    taxa = taxa.apply(lambda x: pd.Series(x.split('; '))).add_prefix('tax_')
    otu_table.add_metadata(taxa.to_dict(orient='index'), 
                           axis='observation')
    
    # Collapses the data
    collapsed = otu_table.collapse(filter_fun, axis='observation', norm=False)
    collapsed = collapsed.to_dataframe().sparse.to_dense()

    # Gets an updated order argument because there may be some levels missing
    keep_order = pd.Series(keep)
    keep_order = keep_order[keep_order.isin(collapsed.index)].values

    collapsed = collapsed.loc[keep_order]
    collapsed.columns.set_names('SampleID', inplace=True)

    return collapsed.T


def _build_strings(large_level, small_level):
    """
    Sets up a formatting string for taxonomy
    """
    lvl = ['d', 'p', 'c', 'o', 'f', 'g', 's']
    level_strs = [f'{{tax_{i+1:1.0f}}}' for i in np.arange(small_level)]
    other_strs = [f'{{{i + 1:1.0f}}}__Other'.format(*lvl) 
                  for i in np.arange(small_level)]
    small_str = ';'.join(level_strs)
    other_str = ';'.join(other_strs)
    big_str = np.hstack([level_strs[:large_level], 
                         other_strs[large_level:]])
    big_str = ';'.join(big_str)

    return small_str, big_str, other_str


def _plot_barplot(ax, x, collapse, group_order, palette, **bar_kws):
    """
    Builds the barplot
    """
    collapse = collapse.loc[x.index, group_order]
    bottom = collapse.cumsum(axis=1) - collapse

    bar_kwargs = dict(width=1, edgecolor='None')
    bar_kwargs.update(**bar_kws)
    
    for o, height in collapse.items():
        ax.bar(x=x['x'], 
               height=height, 
               bottom=bottom[o], 
               color=palette[o],
               **bar_kwargs
              )
               
def _fecal_barplot_legend_wide(ax, level=0, wraps=[],
                               tax_order=bhc.fecal_tax_order_hard,
                               palette=bhc.feces_palette_hard):
    """
    Generates a progrmatic legend for stacked barplots
    """
    taxa_levels = pd.DataFrame([v.split(';') for v in tax_order],
                               index=tax_order)
    taxa_levels['x'] = (~taxa_levels[level].duplicated(keep='first') & 
                        ~taxa_levels[level].isin(wraps)) * 1
    taxa_levels['x'] = taxa_levels['x'].cumsum() - 1
    taxa_levels['y'] = 1
    taxa_levels['y'] =  taxa_levels.groupby('x')['y'].cumsum() + \
        taxa_levels[level].isin(wraps) * 1 - 1
    taxa_levels['c'] = pd.Series(palette)


    # Color block scatter
    ax.scatter(x='x', y='y', c='c', marker='s', data=taxa_levels,)
    ax.set_xlim(-0.05, taxa_levels['x'].max() + 1.05)
    ax.set_ylim(taxa_levels['y'].max() + 1, -2)
    ax.xaxis.set_tick_params(**hide_axes)
    ax.yaxis.set_tick_params(**hide_axes)
    # Adds the text
    for _, [p_, c_, o_, f_, x_, y_, c_] in taxa_levels.iterrows():
        ax.text(x=x_ + 0.10, y=y_, s=f_.replace("__", ". ").replace("_", ' '),
                ha='left', 
                va='center',
                size=7)
    for _, [p_, x_, y_] in taxa_levels[[level, 'x', 'y']].drop_duplicates(level).iterrows():
        ax.text(x=x_ + 0.5, y=y_-1, s=p_.split('__')[-1],
                ha='center',
                va='center',
                size=9,
                weight='bold',
                )

