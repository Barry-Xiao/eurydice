import pandas as pd
import numpy as np
import scipy


from .suppress_counts import suppress_counts


def build_table1(ref_col, meta, meta_dd, col_order, count_col, 
                 var_label=None, col_labels=None, type_remap=None, test=True,
                 categorical_kws=None, 
                 mean_kws=None, median_kws=None, range_kws=None):
    """
    Builds a table 1 based on the desired columns

    Paramaters
    ----------
    ref_col: str
        The column that should be used to group the data (i.e. cohort, 
        outcome, reference)
    meta: DataFrame
        The dataframe with the per-participant characteristics containing
        the ref_col, test_col, and count_col
    meta_dd: DataFrame
        The data dictionary associated with the microbiome data. Each row
        (index) should be a column name and the data type is assumed to be 
        in a column called `dtype` and information for the categorical 
        variables coded in `values`
    col_order: list
        The columns that should be included in the table
    count_col: str
        The name of a reference column which should be used for counting the
        samples. This could be a dummy column, but it cannot have an NaN
        value
    var_label: dict, optional
        A nested dictiomary of dictionaries of categorical columns where
        the column dictionary keys the numeric code to the description
    col_labels: dict, optional
        A dictionary to rename column with a tidier or more intrepretable
        label
    type_remap: dict, optional
        A dictionary to remap the column type based on column name. Format
        options are "integer", "float" and "categorical".
    categorical_kws: dict, optional
        Keyword arguments surrounding categorical representations of data.
        see build_categorical_suppression for accepted arguments
    mean_kws, median_kws, range_kws: dict, optional
        Keywords argument to format the mean, median, and range displays
        for continous data, respectively. See build_continous_suppression
        for options

    Return
    ------
    DataFrame
        A tidied data table 1 in a data frame
    """
    
    # Sets up the arguments for each data type
    cat_kwargs, mean_kwargs, median_kwargs, range_kwargs = \
        _build_kwargs(count_col, categorical_kws, mean_kws, median_kws, 
                      range_kws)

    # Sets up the data kwargs
    data_kwargs = dict(meta=meta, ref_col=ref_col)

    # Checks the column order and classifies the column type
    col_order = [c for c in col_order if c in meta.columns]
    col_types = _classify_col_type(col_order, meta_dd, type_remap)

    # Sets up the var label data
    if var_label is None:
        def split_x(x):
            if (',' in x) and ('|' in x):
                x = x.replace(' | ', '|')
                z = [y.split(',', 1) for y in x.split('|')]
                z = {k.replace(" ", ''): v for k, v in z}
                return z
            else:
                return np.nan
        var_label = meta_dd.loc[col_types['categorical'], 'values'].dropna()
        var_label = var_label.apply(split_x).dropna()

    # ### Summarizes the data ###
    # Categorical
    combo = [
        build_categorical_suppression(col_, **cat_kwargs, **data_kwargs, 
                                     var_names=var_label.get(col_, None))
        for col_ in col_types['categorical']
        ]
    # Float
    combo.extend([
        build_continous_suppression(col_, **cont_kwargs, **data_kwargs,
                                    int_=False)
        for col_ in col_types['float']
        for cont_kwargs in [mean_kwargs, median_kwargs, range_kwargs]
    ])
    # Integers
    combo.extend([
        build_continous_suppression(col_, **cont_kwargs, **data_kwargs,
                                    int_=True)
        for col_ in col_types['integer']
        for cont_kwargs in [mean_kwargs, median_kwargs, range_kwargs]
    ])
    combo = pd.concat(axis=0, objs=combo)
    combo.reset_index(inplace=True)
    combo['codes'] = combo['codes'].astype(int)
    combo.set_index(['column', 'codes', 'labels'], inplace=True)
    combo.sort_index(ascending=True, inplace=True)
    combo = combo.loc[col_order]

    # Cleans up labeling for the p-value column
    combo.columns.set_names([ref_col, 'group_size', 'stat'], inplace=True)
    if test:
        combo['99'] = combo['99'].fillna(method='ffill', axis=1)
        combo.drop(columns=[('99', np.nan, 'c')], inplace=True)
        combo[('99', '', 'c')] = combo[('99', '', 'c')].apply(_tidy_p)
        combo.fillna('', inplace=True)
    else:
        combo.drop(columns=['99'], level=ref_col, inplace=True)
    
    # Cleans up the header labeling
    group_size = {k[1]: '(n={1:1.0f})'.format(*k) for k in combo.columns 
                  if isinstance(k[1], (int, float))}

    combo.rename(columns=group_size, level='group_size', inplace=True)
    combo.rename(columns={'-1': 'All', '99': 'p-value'}, 
                 level=ref_col, 
                 inplace=True)
    
    # Cleans up column names
    combo.rename(index={89: -9, 85: -5, 96: -6}, 
                 level='codes',
                 inplace=True)
    if col_labels is not None:
        combo.rename(index=col_labels,
                     level='column',
                     inplace=True)
    combo.fillna('', inplace=True)
    
    return combo



def build_categorical_suppression(test_col, ref_col, count_col, meta, 
                                  test=True, var_names=None, p_col='p1',
                                  suppress_kws=None, no_l5_perc=False):
    """
    Builds a per-variable entry for a supression table 

    Parameters
    ----------
    ref_col: str
        The column that should be used to group the data (i.e. cohort, 
        outcome, reference)
    test_col: str
        The column in the data that's being described (age, sex, etc)
    meta: DataFrame
        The dataframe with the per-participant characteristics containing
        the ref_col, test_col, and count_col
    count_col: str, optional
        Column to be counted. Should be included int he dataframe and not be 
        either the ref_col or the test_col. An identifier column is ideal.
    test: bool, optional
        Should a statistical test (chi2-test), ignroing any missing values,
        be preformed.
    var_names: dict, optioanl
        A dictionary of columns and the code value substitutions
    p_col: {'p1', 'p2'}
        How the percentage value should be calculated. P1 will include all
        values and calculate a percentage including the missing data. p2 will
        excluding missing from the total and percentage calculation. 
    supress_kws: dict, optional
        Key word arguments to pass to `supress_counts`
    no_l5_perc: bool, optional
        Suppresses percentages for values less than 5

    Returns
    -------
    DataFrame
        A wide dataframe where each row is a group in the test_col and the 
        rows are a categorical level. If a test was supposed to be included, 
        an empty header row with the p-value is provided.

    """
    # Sets up the suprpession arguments
    supress_kwargs = dict(total=True, 
                          round_up=True, 
                          squish_cols=['-6', '-9'])
    if suppress_kws is not None:
        supress_kwargs.update(suppress_kws)
    supress_kwargs['total'] = True
    
    # Assumes we have a dictionary for the label, if we dont then just hides it
    if var_names is None:
        var_names = dict()

    missing_resp ={'-9': '89', '-8': '89', '-7': '89', 
                   '-6': '86', '-5': '85'}
    missing_labels = {'89': 'missing', 
                       '86': 'not applicable',
                       '85': 'inconsistent',
                       '-1': '',
                       }

    # Constructs the count table
    counts = meta.fillna('89').replace(missing_resp)
    counts = counts.groupby([ref_col, test_col])[count_col].count()
    counts = counts.unstack().fillna(0)

    # Calculates counts without missing data
    drop_cols = [c for c in missing_labels.keys() if c in counts.columns]
    counts_nm = counts.copy().drop(columns=drop_cols)
    total_nm = counts_nm.sum(axis=1)
    total_nm['total'] = total_nm.sum()

    # Peforms a chi-square test, if desired and possible
    if test:
        counts_t = counts_nm.loc[counts_nm.sum(axis=1) > 0]
        x2, p, dof, _ = scipy.stats.chi2_contingency(counts_t)
    else:
         p = np.nan
    test_df = pd.DataFrame([['-1', '99', p, False]], 
                           columns=['codes', ref_col, 'c', 'm'])


    # Supresses the counts and converts to a long form matrix for clean up.
    # We also calculat the percentage based on the outcome group or grouping
    # variable. 
    count_m, mask_m, _ = suppress_counts(counts, **supress_kwargs)
    mask_dm = pd.DataFrame({'c': count_m.unstack(), 
                            'm': mask_m.unstack()})
    totals = mask_dm.copy().xs('total', level=test_col)['c']
    mask_dm.drop('total', level=test_col, inplace=True)
    mask_dm.rename({'-9': '89'}, level=test_col, inplace=True)
    mask_dm.reset_index(inplace=True)
    
    mask_dm['group_size'] = mask_dm[ref_col].replace(totals)
    mask_dm['group_nm'] = mask_dm[ref_col].replace(total_nm).mask(
        mask_dm[test_col].isin(missing_labels.keys()))
    mask_dm.rename(columns={test_col: 'codes'}, inplace=True)
    
    # Calculates the percent values
    mask_dm['p1'] = mask_dm['c'] / mask_dm['group_size']
    mask_dm['p1'] = mask_dm['p1'].mask( mask_dm['p1'] > 1, 1)
    mask_dm['p2'] = mask_dm['c'] / mask_dm['group_nm']
    mask_dm['p2'] = mask_dm['p2'].mask( mask_dm['p2'] > 1, 1)

    # Masks the values so we know what was rounded
    def _perc_str(x, str_='({x:>6.1%})'):
        if pd.isnull(x): 
            return ''
        else:
            return str_.format(x=x)
    
    mask_dm['p'] = mask_dm[p_col].apply(_perc_str).mask(
        mask_dm['m'],  
        mask_dm[p_col].apply(_perc_str, str_='(<{x:>5.1%})'))
    mask_dm['p'] = \
        mask_dm['p'].mask(mask_dm['m'] & (mask_dm['c'] == 5) & no_l5_perc, '')
    mask_dm['c'] = mask_dm['c'].apply(_perc_str, str_='{x:1.0f}').mask(
        mask_dm['m'], 
        mask_dm['c'].apply(_perc_str, str_='<{x:1.0f}'))

    # If there is a test result, let's display that
    mask_dm = pd.concat(axis=0, objs=[mask_dm, test_df])
    mask_dm.reset_index(drop=True, inplace=True)

    # Cleans up the labeling
    mask_dm['column'] = test_col
    mask_dm['labels'] = mask_dm['codes'].replace(var_names)
    mask_dm.replace({'labels':missing_labels}, inplace=True)
    mask_dm.replace({ref_col: {'total': '-1'}}, inplace=True)
    
    mask_w = mask_dm.pivot(index=['column', 'codes', 'labels'],
                           columns=[ref_col, 'group_size'],
                           values=['c', 'p'])
    mask_w.drop(columns=[('p', '99', np.nan)], inplace=True)
    mask_w.rename(columns={np.nan: ''}, level='group_size', inplace=True)
    mask_w.columns.set_names(['stat', ref_col, 'group_size'], inplace=True)
    mask_w.columns = \
        mask_w.columns.reorder_levels([ref_col, 'group_size', 'stat'])
    mask_w.sort_index(axis='columns', ascending=True, inplace=True)

    return mask_w


def build_continous_suppression(test_col, ref_col, meta, test=None,
                                c_col='mean', p_col='std', int_=True,
                                pad_=1, ctype=0, hide0miss=True, label=None,
                                drop_header=False):
    """
    Summarizes continous data and suppresses missing values

    parameters
    ----------
    test_col: str
        The column in `meta` which contains the continous data we want to 
        describe by group (i.e. age, bmi, etc) 
    ref_col: str
        The column that should be used to group the data (i.e. cohort, 
        outcome, reference)
    meta: DataFrame
        The dataframe with the per-participant characteristics containing
        the ref_col, test_col, and count_col
    test: function, optional
        This should be a function that take a variable length list of values 
        and returns a test statistic and a p-value. If no test is provided,
        no test will be performed.
    c_col, p_col: str
        The column produced by the decribe function that should be used to 
        describe the data. Options include {'mean', 'std', 'min', 'max', 
        '25%', 'median', '75%', 'iqr', 'range'}
    pad_: int, optional
        The number of additional decimal places hwich should be added to the
        standard deviation measurement.
    ctype: int, optional
        The code value for the description. Useful if multiple descriptive
        values are being produced (i.e. mean (std) and median [IQR])
    hide0miss: bool, optional
        If there are no missing values, this will drop the display row for 
        missing values
    label, optioanl
        How the described data should be labeled. By default, this value will
        be "c_col (p_col)".

    Returns
    -------
    DataFrame
        A wide dataframe where each column is a group in the test_col and the 
        rows are a sumary of the continous variable level. If a test was 
        supposed to be included, an empty header row with the p-value is 
        provided.
    """
    # Cleans the data for summarization
    meta2 = meta.set_index([ref_col], append=True)[[test_col]]
    meta2['count_col'] = '-1'
    meta2.replace({'-6': np.nan, '-7': np.nan, '-8': np.nan, '-9': np.nan, 
                   '1799': np.nan, -6: np.nan, -7: np.nan, -8: np.nan, 
                   -8: np.nan, -9: np.nan, 1799: np.nan},
                  inplace=True)

    # Determines the missing cells 
    miss_ = (meta2.isna() * 89).astype(str).groupby([ref_col, test_col])
    miss_ = miss_['count_col'].count().unstack().fillna(0)
    drop_miss = len(miss_.columns) == 1
    if drop_miss:
        miss_['89'] = 0
    miss_c, miss_m, _ = suppress_counts(miss_, total=True)

    meta2 = meta2.astype(float)

    # Tests, if appropriate
    if test is not None:
        f, p, = test(
            *meta2[test_col].dropna().groupby(ref_col).apply(
                lambda x: x.values).values
        )
    else:
        p = np.nan
    test_df = pd.DataFrame(
        [['99', '99',  p, np.nan, '-1']],
        columns=[ref_col, 'group_size', 'c', 'p', 'code']
    )

    # Summarizes the data and cleans it up for rounding
    describe_cat = pd.concat([meta2.groupby(ref_col)[test_col].describe(),
                              meta2.groupby('count_col')[test_col].describe()
                             ]).drop(columns=['count'])
    
    describe_cat.rename({-1: '-1'}, inplace=True)
    round_pad =_round_data(describe_cat['std'], pad=pad_)
    scaler = np.power(10, round_pad)
    s_fac = np.absolute(np.min([0, -round_pad])).astype(int)
    if int_ & (s_fac > 0): 
        round_pad =_round_data(describe_cat['std'], pad=pad_)
        scaler = np.power(10, round_pad)
        int_cols = ['mean', 'min', 'max', '25%', '50%', '75%']
        describe_cat[int_cols] = describe_cat[int_cols].dropna().round(0).astype(int)
        describe_cat2 = describe_cat[int_cols].applymap(
            lambda x: f'{x:>5.0f}' if not pd.isnull(x) else x)
        describe_cat['std'] = \
            np.round(describe_cat['std'] * scaler, 0) / scaler
        describe_cat2['std'] = describe_cat['std'].apply(
            lambda x: f'{{0:>5.{s_fac:1d}f}}'.format(x) if not pd.isnull(x) 
            else x)
    else:
        describe_cat = np.round(describe_cat * scaler, 0) / scaler
        # Cleans up the text version fo the columns so we have a nice
        # version for working
        
        describe_cat2 = describe_cat.copy().applymap(
            lambda x: f'{{0:>5.{s_fac:1d}f}}'.format(x) if not pd.isnull(x) 
                else x)

    # Formats clean text strings
    describe_cat2['iqr'] = \
        describe_cat2.dropna(subset=['mean']).apply(
            lambda x: '[{0}, {1}]'.format(*x[['25%', '75%']].values), axis=1)
    describe_cat2['range'] = \
        describe_cat2.dropna(subset=['mean']).apply(
            lambda x:'[{0}, {1}]'.format(*x[['min', 'max']].values), axis=1)
    describe_cat2.rename(columns={'50%': 'median'} , inplace=True)
    
    # Adds in additional columns
    describe_cat2['group_size'] = miss_c['total']
    describe_cat2.loc['-1', 'group_size'] = miss_c.loc['total', 'total']
    miss_str =  miss_c['89'].apply(lambda x: f'{x:1.0f}').mask(
        miss_m['89'], miss_c['89'].apply(lambda x: f'<{x:1.0f}'))
    describe_cat2['missing'] = miss_str
    describe_cat2.loc['-1', 'missing'] = miss_str['total']
    describe_cat2.fillna('--', inplace=True)

    # Puts together the long form matrix and adds codes
    demo_wide = pd.concat(axis=0, objs=[
        describe_cat2[['group_size', c_col, p_col]].rename(
            columns={c_col: 'c', p_col: 'p'}),
        describe_cat2[['group_size', 'missing']].rename(
            columns={'missing': 'c'}),
        ])
    demo_wide.index.set_names(ref_col, inplace=True)
    demo_wide.reset_index(inplace=True)

    # Pulls out hte coding
    demo_wide['code'] = \
        (demo_wide['p'].isna() * 89 + demo_wide['p'].notna() * ctype)
    demo_wide['code'] = demo_wide['code'].astype(int).astype(str)

    demo2 = pd.concat(axis=0, objs=[demo_wide, test_df])
    demo2['code'] = demo2['code'].astype(str)
    if label is None:
        label = f'{c_col} ({p_col})'
    demo2['label'] = demo2['code'].replace({'-1': '', 
                                            '89': 'missing',
                                            str(ctype): label})
    demo2['column'] = test_col
    
    demo2 = demo2.pivot(index=['column', 'code', 'label'],
                        columns=[ref_col, 'group_size'],
                        values=['c', 'p'])
    demo2.columns.set_names(['stat', ref_col, 'group_size'], inplace=True)
    demo2.columns = \
        demo2.columns.reorder_levels([ref_col, 'group_size', 'stat'])
    demo2.sort_index(axis='columns', ascending=True, inplace=True)
    demo2.drop(columns=[('99', '99', 'p')], inplace=True)
    demo2.rename(columns={'99': ''}, level='group_size', inplace=True)

    if drop_header:
        demo2.drop(['89', '-1'], level='code', inplace=True)
    elif hide0miss & drop_miss:
        demo2.drop(['89'], level='code', inplace=True)
    return demo2


def _build_kwargs(count_col, cat_kws=None, mean_kws=None, median_kws=None, 
                  range_kws=None, test=True):
    """
    Builds the list of default for column types
    """
    # Cleans up mean arguments, assuming that the test will be on medians 
    mean_kwargs = dict(test=None, 
                       pad_=1, 
                       hide0miss=False,
                       drop_header=True)
    if mean_kws is not None:
        mean_kwargs.update(mean_kws)

    # Range key words
    range_kwargs = dict(test=None, 
                        ctype=2, 
                        pad_=1, 
                        label='min, max', 
                        hide0miss=False,
                        drop_header=True)
    if range_kws is not None:
        range_kwargs.update(range_kws)
    range_kwargs['c_col'] = 'min' 
    range_kwargs['p_col'] = 'max' 

    # Median arguments
    if test:
        median_test = scipy.stats.kruskal
    else:
        median_test = None
    median_kwargs = dict(test= scipy.stats.kruskal,
                         ctype=1, 
                         pad_=1, 
                         label='median [IQR]',
                         hide0miss=True,
                         )
    if median_kws is not None:
        median_kwargs.update(median_kws)
    if not test:
        median_kwargs['test'] = None
    median_kwargs['c_col'] = 'median'
    median_kwargs['p_col'] = 'iqr'

    # Count arguments
    cat_kwargs = dict(p_col='p2', 
                      test=True)
    if cat_kws is not None:
        cat_kwargs.update(cat_kws)
    cat_kwargs['count_col'] = count_col
    cat_kwargs['test'] = test

    return cat_kwargs, mean_kwargs, median_kwargs, range_kwargs


def _check_cols(data, cols):
    """
    Checks columns exist within the reference database
    """
    if cols is None:
        cols = []
    cols2 = [c for c in cols if c in data.columns]

    return cols2


def _classify_col_type(test_cols, data_dict, remap_=None):
    """
    Classififes the samples based on the column type
    """
    col_types = data_dict.loc[test_cols, 'dtype'].copy().to_dict()    
    if remap_ is not None:
        col_types.update(remap_)
    col_types = pd.Series(col_types, name='dtype')
    col_types = col_types.reset_index()

    col_dict = {'categorical': [], 'float': [], 'integer': []}
    col_type2 = col_types.groupby('dtype')['index'].apply(lambda x: x.values)

    col_dict.update(col_type2.to_dict())
    
    return col_dict


def _round_data(vals, pad=1):
    """
    Identifies the rounding degree
    """
    logged = \
        -np.floor(np.log10(np.absolute(vals.replace({0: np.nan}).dropna())))
    rounded = logged + pad

    return rounded.max().max()


def _tidy_p(p):
    """
    Tidies p-values
    """
    if pd.isnull(p):
        return ''
    elif p >= 0.01:
        return f'{p:1.2f}'
    elif p >= 0.001:
        return f'{p:1.3f}'
    elif p < 1e-16:
        return '<1e-16'
    else:
        return f'{p:1.0e}'
    