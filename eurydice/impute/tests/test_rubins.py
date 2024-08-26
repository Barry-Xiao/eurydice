from unittest import TestCase, main

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import statsmodels.formula.api as smf

from test_data import (m1, m2, m3, m4, m5)
from eurydice.impute.rubins import (pool_fits,
                                    _pool_variance,
                                    _pool_est,
                                    _calc_fmi_rvi,
                                    _pool_dof,
                                    _fisher_z,
                                    _throwback_z,
                                    )



class SuppressionTest(TestCase):
    """
    Examples based on 
    https://www.daviddisabato.com/blog/2021/2/13/analyzing-and-pooling-results-from-multiply-imputed-data
    """
    def setUp(self):
        # Sets up the data sets
        imputed_cols = ['rating', 'complaints', 'privileges', 'learning', 
                        'raises', 'critical', 'advance']
        self.imputed = dict(
            m1=pd.DataFrame(data=m1.copy(), columns=imputed_cols),
            m2=pd.DataFrame(data=m2.copy(), columns=imputed_cols),
            m3=pd.DataFrame(data=m3.copy(), columns=imputed_cols),
            m4=pd.DataFrame(data=m4.copy(), columns=imputed_cols),
            m5=pd.DataFrame(data=m5.copy(), columns=imputed_cols),
            )
        # Fits the models
        formula = ' + '.join(['rating ~ complaints + privileges + learning', 
                               'raises + critical + advance'])
        self.fits = {k: smf.ols(formula, data).fit() 
                      for k, data in self.imputed.items()}
        # Extracts the model parameters
        self.params = pd.DataFrame({k: fit.params 
                                     for k, fit in self.fits.items()})
        self.vcov = pd.DataFrame({k: fit.cov_params().unstack() 
                                  for k, fit in self.fits.items()})
        self.m = 5

        # Sets up the known values
        self.param_est = pd.Series({'Intercept': 8.122396499761342,
                                     'complaints': 0.6583405579558438,
                                     'privileges': 0.0015854276635162035,
                                     'learning': 0.2412253160040742,
                                     'raises': 0.013925587536421258,
                                     'critical': 0.07394788777761074,
                                     'advance': -0.17583697065400106},
                                     )

        self.vcov_w = pd.Series({'Intercept': 144.08223120420462,
                                 'complaints': 0.026950430125796544,
                                 'privileges': 0.016184925976195738,
                                 'learning': 0.02757608817783755,
                                 'raises': 0.056244394778998244,
                                 'critical': 0.020180979881808937,
                                 'advance': 0.030029162539505865})
        self.vcov_b = pd.Series({'Intercept': 28.03950702646562,
                                 'complaints': 0.0029856570606628507,
                                 'privileges': 0.008135089372141656,
                                 'learning': 0.008392196632026528,
                                 'raises': 0.007886985587655163,
                                 'critical': 0.0024516494329062466,
                                 'advance': 0.002533061518116835})
        self.vcov_f = pd.Series({'Intercept': 5.607901405293124,
                                 'complaints': 0.0005971314121325701,
                                 'privileges': 0.001627017874428331,
                                 'learning': 0.0016784393264053054,
                                 'raises': 0.0015773971175310326,
                                 'critical': 0.0004903298865812493,
                                 'advance': 0.000506612303623367})
        self.vcov_t = pd.Series({'Intercept': 177.72963963596337,
                                 'complaints': 0.03053321859859197,
                                 'privileges': 0.025947033222765724,
                                 'learning': 0.03764672413626938,
                                 'raises': 0.06570877748418444,
                                 'critical': 0.023122959201296433,
                                 'advance': 0.033068836361246065})

        self.p_rvi = pd.Series({'Intercept': 0.20298843284248158,
                                'complaints': 0.12327429626133453,
                                'privileges': 0.41161472617265077,
                                'learning': 0.2907950029690884,
                                'raises': 0.1526032151670191,
                                'critical': 0.13410151981023477,
                                'advance': 0.0956999723924996})

        self.p_fmi = pd.Series({'Intercept': 0.18931793538026304,
                                'complaints': 0.11734067475482719,
                                'privileges': 0.3762321172813234,
                                'learning': 0.2675036457881243,
                                'raises': 0.14403528824538236,
                                'critical': 0.12723195564530287,
                                'advance': 0.09191958823511696})

        self.r_fmi = np.array([0.18931794, 0.11734067, 0.37623212, 0.26750365,
                               0.14403529, 0.12723196, 0.09191959])

        self.dof_br = pd.Series({'Intercept': 14.911722612376684,
                                 'complaints': 17.603984531158236,
                                 'privileges': 9.017227507721728,
                                 'learning': 12.166606096731558,
                                 'raises': 16.60747013972903,
                                 'critical': 17.236956539688986,
                                 'advance': 18.524847158713516}) 
        self.dof_br_r = np.array([14.911723, 17.603985,  9.017228, 12.166606,
                                  16.607470, 17.236957, 18.524847])

        # Known summary
        self.summary = pd.DataFrame(
            data=np.array([[ 8.12239650e+00,  1.33315280e+01,  6.09262230e-01, 3.22609280e-01,  2.84301378e+01,  1.89317935e-01, 2.02988433e-01,  1.49117226e+01],
                           [ 6.58340558e-01,  1.74737571e-01,  3.76759592e+00, 1.60694543e-03,  3.67702869e-01,  1.17340675e-01, 1.23274296e-01,  1.76039845e+01],
                           [ 1.58542766e-03,  1.61080828e-01,  9.84243551e-03, 3.88034498e-01,  3.64284056e-01,  3.76232117e-01, 4.11614726e-01,  9.01722751e+00],
                           [ 2.41225316e-01,  1.94027638e-01,  1.24325235e+00, 1.77851432e-01,  4.22108742e-01,  2.67503646e-01, 2.90795003e-01,  1.21666061e+01],
                           [ 1.39255875e-02,  2.56337234e-01,  5.43252626e-02, 3.92371073e-01,  5.41799674e-01,  1.44035288e-01, 1.52603215e-01,  1.66074701e+01],
                           [ 7.39478878e-02,  1.52062353e-01,  4.86299773e-01, 3.47257968e-01,  3.20487885e-01,  1.27231956e-01, 1.34101520e-01,  1.72369565e+01],
                           [-1.75836971e-01,  1.81848388e-01, -9.66942694e-01, 2.43382469e-01,  3.81274886e-01,  9.19195882e-02, 9.56999724e-02,  1.85248472e+01]]),
            columns=['param', 'bse', 't', 'p', 'width', 'fmi', 'rvi', 'dof_br'],
            index=self.dof_br.index)

        # Know values vs R
        self.r_mice = pd.DataFrame(
            data=np.array([[ 8.1223965  , 13.331528,  0.60926223 , 14.911723, 0.551521122],
                           [ 0.658340558, 0.1747376,  3.767595918, 17.603985, 0.001456847],
                           [ 0.001585428, 0.1610808,  0.009842436,  9.017228, 0.992361317],
                           [ 0.241225316, 0.1940276,  1.243252348, 12.166606, 0.237199557],
                           [ 0.013925588, 0.2563372,  0.054325263, 16.60747 , 0.95732419 ],
                           [ 0.073947888, 0.1520624,  0.486299773, 17.236957, 0.632879679],
                           [-0.175836971, 0.1818484, -0.966942694, 18.524847, 0.34602789 ]],),
            index=self.param_est.index,
            columns=['param', 'bse', 't', 'df', 'p-value']
            )




    def test_pool_fits(self):
        summary = pool_fits(list(self.fits.values()))
        pdt.assert_frame_equal(summary, self.summary)
        # Chekcs against the pooled mice calculation
        npt.assert_array_almost_equal(summary['param'].values,
                                     self.r_mice['param'].values)
        npt.assert_array_almost_equal(summary['bse'].values,
                                     self.r_mice['bse'].values)
        npt.assert_array_almost_equal(summary['t'].values,
                                     self.r_mice['t'].values)
        npt.assert_array_almost_equal(summary['dof_br'].values,
                                     self.r_mice['df'].values)

    def test_pool_est(self):
        test = _pool_est(self.params)
        pdt.assert_series_equal(test, self.param_est)
        # We allow 6 decimal places of percision in the array vs the R values
        npt.assert_array_almost_equal(test.values, 
                                      self.r_mice['param'].values)

    def test_pool_variance(self):
        var_w, var_b, var_f, var_t = \
            _pool_variance(self.params, self.vcov, self.m)
        pdt.assert_series_equal(var_w, self.vcov_w)
        pdt.assert_series_equal(var_b, self.vcov_b)
        pdt.assert_series_equal(var_f, self.vcov_f)
        pdt.assert_series_equal(var_t, self.vcov_t)

    def test_calc_fmi_rvi(self):
        fmi, rvi = _calc_fmi_rvi(self.vcov_b, self.vcov_f, self.vcov_t, 
                                 self.vcov_w, self.m)
        pdt.assert_series_equal(rvi, self.p_rvi)
        pdt.assert_series_equal(fmi, self.p_fmi)
        npt.assert_array_almost_equal(fmi.values, self.r_fmi)

    def test_calc_dof(self):
        dof_br = _pool_dof(self.p_rvi, self.p_fmi, self.m, 
                                           self.fits['m1'].df_resid)
        pdt.assert_series_equal(dof_br, self.dof_br)
        npt.assert_array_almost_equal(dof_br.values, self.dof_br_r)



    


if __name__ == '__main__':
    main()

