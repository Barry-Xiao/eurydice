from unittest import TestCase, main

import os
import warnings

import numpy as np
import numpy.testing as npt
import pandas as pd
import pandas.testing as pdt
import statsmodels.formula.api as smf


from eurydice.impute.harel_adonis import (pool_adonis,
                                          _fisher_z,
                                          _release_z,
                                          _pool_variance_simp,
                                          _pool_r2,
                                          _get_num_obs,
                                          _oxford_comma,
                                          _check_adonis_cols,
                                          _check_adonis_idx,
                                          )

warnings.simplefilter(action='ignore', 
                      category=pd.errors.PerformanceWarning)

class AdonisPoolTest(TestCase):
    def setUp(self):
        # The correlation testing and z-testing
        self.fisher_r2 = pd.Series([0.25])
        self.num_obs = pd.Series([12])
        self.fisher_q =  pd.Series(np.array([0.5493061443340549]))
        self.fisher_var = pd.Series(np.array([0.11111111111111]))

        imputed = pd.read_csv(
             os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                'data/test_res.tsv'),
             sep='\t', 
             dtype=str
             )
        imputed.set_index(['name', 'var_order'], inplace=True)
        self.impute_res = imputed.astype(float)

    def test_pool_adonis_perr(self):
        known = ('The p_method is Jesper.\nAcceptable values are "mean", '
                 '"median", and "max".\nPlease supply a correct option '
                 'for p_method.')
        with self.assertRaises(ValueError) as e_:
            pool_adonis(self.impute_res, p_method='Jesper')
        self.assertEqual(known, str(e_.exception))

    def test_pool_adonis(self):
        known = pd.DataFrame(
            data=np.array([[2.39000000e+02, 2.36471917e+02, 9.22849967e-01,
                            9.02063747e-01, 9.39371134e-01,         np.nan,
                            2.20509005e-03, 2.20616820e-03],
                           [1.00000000e+00, 4.99810567e-01, 7.84331842e-03,
                            1.29042888e-03, 4.42434149e-02, 2.80000000e-03,
                            3.78857946e-04, 3.78889830e-04],
                           [3.00000000e+00, 1.99937671e+00, 1.19279284e-02,
                            2.27267347e-04, 5.29836414e-02, 3.65400000e-01,
                            3.11624567e-04, 3.11646141e-04],
                           [3.00000000e+00, 1.99892407e+00, 1.91005582e-02,
                            2.05947215e-04, 6.65030355e-02, 1.30000000e-03,
                            5.37899838e-04, 5.37964100e-04],
                           [2.00000000e+00, 1.19852576e+00, 1.25530568e-02,
                            1.50620749e-04, 5.42624016e-02, 3.60000000e-03,
                            1.22833636e-03, 1.22867124e-03],
                           [1.00000000e+00, 4.99635004e-01, 2.57162720e-02,
                            1.36805651e-03, 7.77982566e-02, 1.00000000e-03,
                            7.29961663e-04, 7.30079986e-04]]),
            columns=['dof_obs', 'dof_br', 'r2', 'r2_ci_lo', 'r2_ci_hi', 'p_mean', 'fmi', 'rvi'],
            index=pd.MultiIndex.from_arrays(
                np.array([['Residual', 'age', 'alcohol_use', 'birth_country', 'bmi_class', 'ethnicity'],
                          ['6', '1', '4', '5', '3', '2']]),
                names=['name', 'var_order'])
            )
        test = pool_adonis(self.impute_res, 
                           index_levels=['name', 'var_order'])
        pdt.assert_frame_equal(known, test)

    
    def test_oxford_comma_2(self):
        test = _oxford_comma(['cats', 'dogs'])
        known = 'cats and dogs'
        self.assertEqual(test, known)

    def test_oxford_comma_3(self):
        test = _oxford_comma(['cats', 'dogs', 'kids'])
        known = 'cats, dogs, and kids'
        self.assertEqual(test, known)

    def test_check_adonis_cols_pass(self):
        _check_adonis_cols(self.impute_res, 'r2', 'pval', 'df', 
                           'imputation')

    def test_check_adonis_cols_fail(self):
        known = ('We expected to find columns named R2 and DOF in the '
                 'adonis dataframe, but they were absent.\nPlease make '
                 'sure the column names supplied for the r2_col and '
                 'df_col are included in the dataframe')
        with self.assertRaises(ValueError) as e_:
            _check_adonis_cols(self.impute_res, 'R2', 'pval', 'DOF', 
                               'imputation')
        self.assertEqual(known, str(e_.exception))

    def test_check_adonis_index_pass_no_model(self):
        test = _check_adonis_idx(self.impute_res, 'name', None, None)
        known = ['name']
        self.assertTrue(isinstance(test, list))
        npt.assert_array_equal(known, test)

    def test_check_adonis_index_miss_model(self):
        test = _check_adonis_idx(self.impute_res, 'name', 'var_order', None)
        known = ['var_order', 'name']
        self.assertTrue(isinstance(test, list))
        npt.assert_array_equal(known, test)

    def test_check_adons_index_missing_model(self):
        known = (f'A value was given for the model level (model), but '
                 'this variable was not included in the the listed index '
                 'levels.\nPlease make sure your model_level in in your '
                 'index_levels list.')
        with self.assertRaises(ValueError) as e_:
            _check_adonis_idx(self.impute_res, 'name', 'model', 
                              ['name', 'var_order'])
        self.assertEqual(known, str(e_.exception))

    def test_check_adons_index_missing_name(self):
        known = (f'The name_level (model) must be one of the '
                 'index_levels\nPlease check your name_level and '
                 'index_levels')
        with self.assertRaises(ValueError) as e_:
            _check_adonis_idx(self.impute_res, 'model',  'name',
                              ['name', 'var_order'])
        self.assertEqual(known, str(e_.exception))

    def test_check_adonis_miss_levels(self):
        known = ('The index_list includes model, but '
                 'these are not levels in the adonis dataframe '
                 'index. Please check the index_list columns '
                 'and the adonis dataframe index.')
        with self.assertRaises(ValueError) as e_:
            _check_adonis_idx(self.impute_res, 'name', None, 
                              ['name', 'model'])
        self.assertEqual(known, str(e_.exception))
        
    def test_get_num_obs_no_model_pass(self):
        data = self.impute_res.loc[['age', 'Total']].copy()
        data.set_index('imputation', inplace=True, append=True)
        known = pd.Series(np.array([250.] * len(data)), 
                          index=data.index, 
                          name='num_obs')
        test = _get_num_obs(data)
        pdt.assert_series_equal(test, known)

    def test_get_num_obs_model_pass(self):
        # Sets up data with two differetn model levels
        data_a = self.impute_res.loc[['age', 'Total']].copy()
        data_b = self.impute_res.loc[['age', 'Total']].copy()
        data_a['model'] = 'A'
        data_b['model'] = 'B'
        data_b.replace({"df": {249: 499}}, inplace=True)
        data = pd.concat(axis=0, objs=[data_a, data_b])
        data.set_index(['model', 'imputation'], append=True, inplace=True)
        
        known = pd.Series(
            data=np.hstack([[250.] * len(data_a), [500.] * len(data_b)]),
            index=data.index,
            name='num_obs')

        test = _get_num_obs(data,  model_level='model')

        pdt.assert_series_equal(test, known)

    def test_get_num_obs_model_fail(self):
        # Sets up data with two differetn model levels
        data_a = self.impute_res.loc[['age', 'Total']].copy()
        data_b = self.impute_res.loc[['age', 'Total']].copy()
        data_a['model'] = 'A'
        data_b['model'] = 'B'
        data_b.replace({"df": {249: 499}}, inplace=True)
        data = pd.concat(axis=0, objs=[data_a, data_b])
        data.set_index(['model', 'imputation'], append=True, inplace=True)

        err_str = ('Each model should represent only one set of '
                   'input data.\nPlease review the values in '
                   'model_level and correct the label so the '
                   'number of observations are consistent')
        
        # If there isnt an appropriate model level, we should get an error
        with self.assertRaises(ValueError) as e_:
            _get_num_obs(data,  model_level=None)
        self.assertEqual(str(e_.exception), err_str)
    
    def test_fisher_z(self):
        test_q, test_v = _fisher_z(self.fisher_r2, self.num_obs)

        pdt.assert_series_equal(self.fisher_q, test_q)
        pdt.assert_series_equal(self.fisher_var, test_v)

    def test_release_z(self):
        test_r2 = _release_z(self.fisher_q)
        pdt.assert_series_equal(self.fisher_r2, test_r2)

    def test_pool_variance_simp(self):
        var = pd.DataFrame(
            data=np.array([[0.00404858] * 10]),
            index=['A']
        )
        param = pd.DataFrame(
            data=np.array([[0.08973622, 0.08796234, 0.08785097, 0.09045728, 
                            0.09035581, 0.08890988, 0.08730415, 0.08727498, 
                            0.08937271, 0.08872715]]),
            index=['A']
        )
        k_var_w = pd.Series({'A': 0.00404858})
        k_var_b = pd.Series({'A': 1.39492858e-06})
        k_var_f = pd.Series({'A': 1.39492858e-07})
        k_var_t = pd.Series({'A': 0.00405012})

        var_w, var_b, var_f, var_t = _pool_variance_simp(param, var, 10)
        
        pdt.assert_series_equal(k_var_w, var_w)
        pdt.assert_series_equal(k_var_b, var_b)
        pdt.assert_series_equal(k_var_f, var_f)
        pdt.assert_series_equal(k_var_t, var_t)
    

if __name__ == '__main__':
    main()
