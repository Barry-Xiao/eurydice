import pandas as pd
import numpy as np
from scipy.stats import norm as z

from ._rubins import (_pool_est, _calc_fmi_rvi, _pool_dof)



def pool_adonis(adonis, 
                ci_alpha=0.05, 
                p_method='mean', 
                name_level='name', 
                index_levels=None, 
                model_level=None,
                impute_col='imputation', 
                r2_col='r2', 
                p_col='pval',
                df_col='df'
                ):
    """
    Pools Adonis imputed MICE results using Harel's R2 pooling

    Parameters
    ----------
    adonis_res : pd.DataFrame
        A dataframe with a multi level index including the adonis results.
        The frame should be indexed by a multiple level index whcih is 
        based on the `index_col` plus the `impute_level`, where the 
        variable name is in the `name_level` and the imputation identifier 
        is in the `impute_level` level.
        The dataframe should also include columns summarizing the 
        correlation coeffecient (`r2_col`), p-value (`p_col`), and 
        degrees of freedom (`df_col`).
    ci_alpha: float, (0, 1)
        The critical value to calculate the confidence interval for R2
        values
    p_method: {'mean', 'median', 'max'}
        The method for pooling the p-value across imputations. Options
        are mean, median, and max.
    name_level: str
        The level name for the index portion identifying each unique
        variable. This must be part of the index_levels list.
    model_level: str
        And identifier for a model or equation where the number of
        observations is contingent on the model. Relevant only when
        `adonis` contains observations from multiple models with different
        sample sizes. If this is included, it must be part of the 
        index_levels list.

    impute_col: str
        The name of a column containing a dummy variable identifying 
        each unique imputation. 
    df_col, r2_col, p_col: str
        Column names in adonis_res file that describes the variance 
        explained (r2_col), permutative p-value (p_col), and the degrees
        of freedom (df_col)

    Returns
    -------
    pd.DataFrame
    

    References
    ----------
    [1] Harel (2009). "The estimation of R2 and adjusted R2 in incomplete 
        data sets using multiple imputation". Journal of Applied Statitics. 
        36:1109. doi: 10.1080/02664760802553000
    """
    # Checks the p method
    p_opts = {'mean': np.mean, 'max': np.max, 'median': np.median}
    p_fun = p_opts.get(p_method, 'fail')
    if p_fun in {'fail'}:
        raise ValueError(f'The p_method is {p_method}.\nAcceptable '
                         'values are "mean", "median", and "max".\n'
                         'Please supply a correct option for p_method.')
    # Checks the adonis data frame index
    _check_adonis_cols(adonis, r2_col, p_col, df_col, impute_col)
    index_levels = _check_adonis_idx(adonis, name_level, model_level, 
                                     index_levels)
    adonis = adonis.copy().set_index(impute_col, append=True)

    # Determines the number of observations per model
    num_obs = _get_num_obs(adonis, name_level, model_level, df_col)
    num_obs.drop(['Total'], level=name_level, inplace=True)

    # Extracts the degrees of freedom for each model
    dof = adonis.groupby(index_levels)[df_col].first()
                    
    # Uses the fisher z transform and casts the R2 value to a wide
    # matrix to make life better
    r2 = adonis[r2_col].drop(['Total'], level=name_level)

    # Converts the num_obs and r2 to fisher z transformed data and
    # casts to wide so each column is an imputation
    q_long, v_long = _fisher_z(r2, num_obs)
    q = q_long.reset_index().pivot(index=index_levels, 
                                   columns=impute_col, 
                                   values=r2_col)
    v = v_long.reset_index().pivot(index=index_levels,
                                   columns=impute_col,
                                   values=v_long.name)
    num_impute = len(q.columns)

    # Pools the q values
    q_pool = _pool_est(q)
                    
    # Calculates the pooled variance
    var_w, var_b, var_f, var_t = _pool_variance_simp(q, v, num_impute)
    bse_pool = np.sqrt(var_t)

    # Calculates the CI on the q-value
    q_width = z.ppf(1 - ci_alpha / 2) * bse_pool
    q_ci_lo = q_pool - q_width
    q_ci_hi = q_pool + q_width

    # Gets the fmi and rvi
    fmi, rvi = _calc_fmi_rvi(var_b, var_f, var_t, var_w, num_impute)
    
    # Calculates the adjusted dof based on the rvi and fmi
    dof_br = _pool_dof(rvi, fmi, num_impute, dof)

    # Calculates a pooled p-value, based on the specified mode
    p_pool = adonis.groupby(index_levels)[p_col].apply(p_fun)

    # Combines the results into a summary
    summary = pd.DataFrame({'dof_obs': dof,
                            'dof_br': dof_br,
                            'r2': _release_z(q_pool),
                            'r2_ci_lo': _release_z(q_ci_lo),
                            'r2_ci_hi': _release_z(q_ci_hi),
                            f'p_{p_method}': p_pool,
                            'fmi': fmi,
                            'rvi': rvi,
                            })
    return summary


def _oxford_comma(x, final_join=' and ', spacer=', '):
    if len(x) < 3:
        str_ = final_join.join(x)
    else:
        penultimate = spacer.join(x[:-1])
        last = x[-1]
        str_ = f'{penultimate}{spacer}{final_join}{last}'
    str_ = str_.replace('  ', ' ')
    return str_


def _check_adonis_cols(adonis_df, r2_col, p_col, df_col, impute_col):
    """
    Checks the supplied adonis dataframe is compatible with the function
    """
    # Checks the required columns are present and raises an error they
    # they arent
    col_des = {'p_col': p_col, 
               'r2_col': r2_col, 
               'df_col': df_col,
               'impute_col': impute_col}
    cols = adonis_df.columns
    col_pass = np.all([c in cols for c in col_des.values()])
    if not col_pass:
        missed = [cat for cat, col in col_des.items() if col not in cols]
        miss_cols = _oxford_comma([col_des[c] for c in missed])
        miss_cats = _oxford_comma(missed)

        raise ValueError(f'We expected to find columns named {miss_cols} '
                         'in the adonis dataframe, but they were absent.'
                         '\nPlease make sure the column names supplied '
                         f'for the {miss_cats} are included in the '
                         'dataframe')


def _check_adonis_idx(adonis_df, name_level, model_level, index_levels):
    """
    Checks the index levels in the adonis dataframe
    """
    adonis_idx = adonis_df.index.to_frame().copy().reset_index(drop=True)
    idx_cols = adonis_idx.columns

    # If the column levels are not supplied, builds a relevant variable
    if (index_levels is None) and (model_level is not None):
        index_levels = [model_level, name_level]
    elif (index_levels is None):
        index_levels = [name_level]
    
        
    # Checks the model_level is in the index levels, if needed
    if ((model_level is not None) and (model_level not in index_levels)):
        raise ValueError((f'A value was given for the model level '
                          f'({model_level}), but this variable was not '
                          'included in the the listed index levels.\n'
                          'Please make sure your model_level in in your '
                          'index_levels list.'
                          ))
    # Checks the name varaible is in the index_levels list
    if name_level not in index_levels:
        raise ValueError(f'The name_level ({name_level}) must '
                         'be one of the index_levels\n'
                         'Please check your name_level and index_levels')

    # Checks the index_levels are in the index
    missed_cols = [c for c in index_levels if c not in idx_cols]
    if len(missed_cols) > 0:
        missed_list = _oxford_comma(missed_cols)
        raise ValueError(f'The index_list includes {missed_list}, but '
                         'these are not levels in the adonis dataframe '
                         'index. Please check the index_list columns '
                         'and the adonis dataframe index.')

    return index_levels


def _get_num_obs(adonis, name_level='name', model_level=None, df_col='df'):
    """
    Gets the number of observations based on the total degrees of freedom
    """
    dof = \
        adonis.xs('Total', level='name', axis='index')[df_col].copy() + 1
    if model_level is None:
        total_dof = {'ref_': dof.drop_duplicates().values}
    else:
        total_dof = dof.groupby(model_level).apply(
            lambda x: x.unique()).to_dict()

    # Checks there is only one dof measure for each mdoel
    if np.any([len(v) > 1 for v in total_dof.values()]):
        raise ValueError('Each model should represent only one set of '
                         'input data.\nPlease review the values in '
                         'model_level and correct the label so the '
                         'number of observations are consistent')
    else:
        total_dof = {k: v[0] for k, v in total_dof.items()}

    # Gets the number of observations for each model
    num_obs = adonis.index.to_frame()
    if model_level is None:
        num_obs['num_obs'] = total_dof['ref_']
    else:
        num_obs['num_obs'] = \
            num_obs[model_level].replace(total_dof)

    return num_obs['num_obs']
    


def _pool_r2(adonis, num_obs, alpha=0.05, name_col='name', 
             impute_col='imputation',
             index_cols=['name', 'var_order'], r2_col='r2'):
    """
    """

    # Casts the data into a correlation matrix and gets the variance
    # in that parameter, based on a fisher Z transform
    q_long, v_long = \
        _fisher_z(adonis[r2_col].drop(['Total'], level=name_col), 
                  num_obs.drop(['Total'], level=name_col))
                 
    # Pivots the variable to longer so we have each column as an imputed
    # transformed value
    q_wide = q_long.reset_index()
    q_wide = q_wide.pivot(index=index_cols,
                          columns=impute_col,
                          values=r2_col)
    v_wide = v_long.reset_index()
    v_wide = v_wide.pivot(index=index_cols,
                          columns=impute_col,
                          values=v_long.name)
    m = len(q_wide.columns)
    
    # # Calculates the pooled q estimate
    q_pool = _pool_est(q_wide)

    # # Calculates the pooled variance estimate
    var_w, var_b, var_f, var_t = _pool_variance_simp(q_wide, v_wide, m)
    bse_t = np.sqrt(var_t)

    # Gets the fmi and rvi
    fmi, rvi = _calc_fmi_rvi(var_b, var_f, var_t, var_w, m)

    # Calculates the CI on the q-value
    q_width = scipy.stats.norm.ppf(1 - alpha / 2) * bse_t
    q_ci_lo = q_pool - q_width
    q_ci_hi = q_pool + q_width

    # Summarizes the results
    r2_summary = pd.DataFrame({
        'r2': _release_z(q_pool),
        'r2_ci_lo':  _release_z(q_ci_lo),
        'r2_ci_hi': _release_z(q_ci_hi),
        'fmi': fmi,
        'rvi': rvi,
    })

    return r2_summary


def _fisher_z(r2, num_obs):
    """
    Performs a fisher's Z transform for R2 correlation data

    Parameters
    ----------
    r2: pd.Series
        The correlation coeffecients for the analysis
    num_obs: pd.Series
        The number of observations in the orginal data set

    Returns
    pd.Series
        The z-transformed correlation coeffecient
    pd.Series
        The standard deviation in the z-transformed correlation coeffecient

    References
    ----------
    [1] Harel (2009). "The estimation of R2 and adjusted R2 in incomplete 
        data sets using multiple imputation". Journal of Applied Statitics. 
        36:1109. doi: 10.1080/02664760802553000
    [2] Fisher Z transform. (2 June 2024) In Wikipedia.
        https://en.wikipedia.org/wiki/Fisher_transformation 

    Also See
    --------
    _throwback_z
    """
    ra = np.sqrt(r2)
    q = 0.5 * np.log((1 + ra) / (1 - ra))
    v = np.square(1 / np.sqrt(num_obs - 3))

    return q, v


def _release_z(q):
    """
    Undoes the Fisher's Z transformation

    Parameters
    ----------
    q: pd.Series
        The fisher z-transforemd correlation coeffecient
    
    Returns
    -------
    pd.Series
        The untransfromed correlation coeffecient values    

    References
    ----------
    [1] Harel (2009). "The estimation of R2 and adjusted R2 in incomplete 
        data sets using multiple imputation". Journal of Applied Statitics. 
        36:1109. doi: 10.1080/02664760802553000
    [2] Fisher Z transform. (2 June 2024) In Wikipedia.
        https://en.wikipedia.org/wiki/Fisher_transformation 

    Also See
    --------
    _fisher_z
    """
    ra = np.tanh(q)
    r2 = np.square(ra)

    return r2


def _pool_variance_simp(params, variance, m):
    """
    Pools the variance estimates, not variance/covariance
    """
    var_w = variance.mean(axis=1)

    coef_pool = _pool_est(params).to_frame().values
    vcov_b = np.dot((params - coef_pool), (params - coef_pool).T) / (m-1)
    var_b = pd.Series(np.diag(vcov_b), index=var_w.index)
    var_f = var_b / m
    var_t = var_b + var_w + var_f

    return var_w, var_b, var_f, var_t
    
    
