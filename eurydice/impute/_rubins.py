import pandas as pd
import numpy as np
from scipy.stats import t as t_dist


def pool_fits(fits_, alpha=0.05):
    """
    Pools fit parameter estimates after multiple imputation
    
    Parameters
    ----------
    fits: list of statsmodels fits
        A list of statsmodels fits from imputed data that need to be pooled
        for an effect estimate
    alpha: float
        The critical value for the two-tailed threshhold
    
    Returns
    -------
    pd.DataFrame
        A dataframe summarizing the fit result including:
            `param`: the pooled parameter estimates
            `bse`: the pooled standard error estimate
            `t`: the t-value for the fit data
            `p`: the p-value based on the calculated t-value for a t-
                distribution with Barnard and Rubin (1999) adjusted degrees of 
                freedom.
            `width`: the width of the 95% confidence interval, based on a t
                distribution with Barnard and Rubin (1999) adjusted degrees
                of freedom
            `fmi`: The Fraction of Missing data estimate (uncorrected)
            `rvi`: The relative variance increase due ot missingness
            `dof_obs`: Imputation-corrected residual degrees of freedom
            `dof_r`: The Rubuin degree of freedom for the pooled estimate
    
    References
    ----------
    [1] Enders, CK. (2022) "Chapter 7. Multiple Imputation." Applied 
        Missing Data Analysis. 2nd ed. New York, USA: Guilford Press. Online
        version. ISBN 9781462549863.

    [2] Hyemens, MW and Eekhout, I (2019). Applied Missing Data Analysis with
        SPSS and (R)studio. 1st ed. Amsterdam, Neatherlands. Online.
        https://bookdown.org/mwheymans/bookmi/

    [3] Diasbato, David. (2021) "Analyzing and Pooling Results From Multiply 
        Imputed Data". David Diasbato, blogspot. 13 February 2021.
        http://www.daviddisabato.com/blog/2021/2/13/analyzing-and-pooling-
        results-from-multiply-imputed-data

    """
    # Pulls out the parameters and variace-covariance matrix for the fits
    params = pd.DataFrame({i: fit.params for i, fit in enumerate(fits_)})
    vcovar = pd.DataFrame({i: fit.cov_params().unstack() 
                           for i, fit in enumerate(fits_)})
    # Gets hte number of imputations and number of remaining degrees of 
    # freedom from the original modeling
    m = len(fits_)
    dof_res = fits_[0].df_resid
    
    # Gets the parameter estimate
    param_pool = _pool_est(params)
    
    # Gets the variance estimates
    var_w, var_b, var_f, var_t = _pool_variance(params, vcovar, m)
    bse_pool = np.sqrt(var_t)
    
    # Fraction of missing data and relative increase in variance
    fmi, rvi = _calc_fmi_rvi(var_b, var_f, var_t, var_w, m)
    
    # Gets the degrees of freedom
    dof_br = _pool_dof(rvi, fmi, m, dof_res)
    dof_br.fillna(dof_res, inplace=True)
    
    # Gets the t-statistics and p-value
    t_pool = param_pool / bse_pool
    t_crit = t_dist.ppf(1 - alpha / 2, dof_br)
    ci_width = t_crit * bse_pool
    p_val = t_dist.pdf(t_pool, dof_br)
    
    summary = pd.DataFrame({
        'param': param_pool,
        'bse': bse_pool,
        't': t_pool,
        'p': p_val,
        'width': ci_width,
        'fmi': fmi,
        'rvi': rvi,
        'dof_br': dof_br,
        # 'dof_res': dof_res,
        })
    
    return summary


def _pool_est(params):
    """
    Pools the parameter estimate
    
    References
    ----------
    eq. 7.12 ref [1]_
    pg 9 ref [3]_
    """
    pooled = params.mean(axis=1)

    return pooled


def _pool_variance(params, vcovar, m):
    """
    Calculates the pooled variance for imputed data
    
    References
    ----------
    eq.7.13 - 
    """
    # Within variance - mean of the variance (eq. 7.13; ...)
    vcov_w = vcovar.mean(axis=1).unstack()
    
    # Between imputation variance (eq 7.14)
    coef_pool = _pool_est(params).to_frame().values
    vcov_b = pd.DataFrame(
        data=np.dot((params - coef_pool), (params - coef_pool).T) / (m - 1),
        index=params.index,
        columns=params.index
        )
    
    # Finite covaraince (eq 7.15, infer)
    vcov_f = vcov_b / m
    
    # Total variance is the within, between, and finite eq 7.15
    vcov_t = vcov_w + vcov_b + vcov_f
    
    var_t = pd.Series(np.diag(vcov_t), index=vcov_t.index)
    var_w = pd.Series(np.diag(vcov_w), index=vcov_w.index)
    var_b = pd.Series(np.diag(vcov_b), index=vcov_b.index)
    var_f = pd.Series(np.diag(vcov_f), index=vcov_f.index)
    
    return var_w, var_b, var_f, var_t

    
def _calc_fmi_rvi(var_b, var_f, var_t, var_w, m):
    """
    Calculates the fraction of missing data and relative increase in variance
    
    References
    ----------
    [1]
    [2]
    """
    fmi = (var_b + var_f) / var_t
    df_r = (m - 1) * (1 + 1 / np.square(fmi))
    rvi = ((var_b + var_f) / var_t) * (df_r + 1) / (df_r + 3) + 2 / (df_r + 3)
    
    return fmi, rvi


def _pool_dof(riv, fmi, m, dof_res):
    """
    Calculates the pooled degrees of freedom
    
    References
    ----------
    [1] https://www.daviddisabato.com/blog/2021/2/13/analyzing-and-pooling-results-from-multiply-imputed-data
    """
    v = (m - 1) / np.square(fmi)
    vhat = (1 - fmi) * ((dof_res + 1) / (dof_res + 3)) * dof_res
    dof_br = np.power(1 / v + 1 / vhat, -1)


    return dof_br
    # # Gets the Rubins dof (eq 7.19)
    # dof_r = (m - 1) * (1 + 1 / np.square(fmi))
    # # Gets the observed df (eq 7.21)
    # dof_obs = dof_res * (dof_res + 1)/(dof_res + 3)*(1 - riv)
    # dof_br = (dof_r * dof_obs) / (dof_r + dof_obs)
    
    # return dof_r, dof_obs, dof_br
