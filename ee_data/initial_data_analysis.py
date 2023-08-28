import pandas as pd
import numpy as np
import scipy.stats as stats


commercial = [
    '6_1_1_1',
    '1_1_1_1',
    'aa_1',
    '5_1_1_2',
    '190_1_1_2',
    '6_1_1_21',
    '3_1_1_21',
    '1_1_1_7',
    '1_1_1_3',
    '1_1_1_22',
    '1_1_1_26',
    '1_1_1_24',
    '11_1_1_30',
    '1_1_1_30',
    '1_1_1_27',
    '6_1_1_14',
    '1_1_1_29',
    '1_1_1_15',
    '1_1_1_20',
    '1_1_1_18',
    'aa_18'
]

nns = ['30_3_1_11',
        '48_2_1_10',
        '124_2_1_11',
        '77_1_2_11',
        '189_1_2_11']

training_set_succesful = [
    '1_4_1_2',
    '3_1_2_18',
    '14_1_2_14',
    '16_3_1_9',
    '22_4_4_28',
    '56_2_1_1',
    '73_3_1_29',
    '154_1_2_15',
    '171_2_2_17',
    '185_1_2_2',
    '185_2_1_10',
    '187_1_4_2',
    '187_1_4_30',
    '200_3_1_21',
    '225_1_1_2',
    '225_1_1_13',
    '249_4_4_3',
    '250_1_3_12',
    '252_1_1_8',
    '254_2_2_11',
]

joc_ligands = [
    '90_1_1_17',
    '7_1_2_2',
    '7_2_1_2',
    '14_1_1_13',
    '281_4_4_2',
    '227_2_2_2'
]

if __name__ == '__main__':
    # read in data
    all_data1 = pd.read_csv('data/processed_data_initial_screening_1a-1e.csv', header = 0, index_col= 0)
    all_data1 = all_data1.loc[:, ['ee']]

    # separate out commercial and training sets
    only_ts = all_data1.loc[[i for i in all_data1.index if '_'.join(i.split('_')[:-1]) in training_set_succesful + joc_ligands], :]
    only_com = all_data1.loc[[i for i in all_data1.index if ('_'.join(i.split('_')[:-1]) in commercial) ], :]

    b1a = only_ts.loc[ [i for i in only_ts.index if i[-2:] == '1A'], :]
    b1b = only_ts.loc[ [i for i in only_ts.index if i[-2:] == '1B'], :]
    b1c = only_ts.loc[ [i for i in only_ts.index if i[-2:] == '1C'], :]
    b1f = only_ts.loc[ [i for i in only_ts.index if i[-2:] == '1D'], :]
    b1g = only_ts.loc[ [i for i in only_ts.index if i[-2:] == '1E'], :]

    ts_data = [b1a.to_numpy().ravel(), b1b.to_numpy().ravel(), b1c.to_numpy().ravel(), b1f.to_numpy().ravel(), b1g.to_numpy().ravel()]

    c1a = only_com.loc[ [i for i in only_com.index if i[-2:] == '1A'], :]
    c1b = only_com.loc[ [i for i in only_com.index if i[-2:] == '1B'], :]
    c1c = only_com.loc[ [i for i in only_com.index if i[-2:] == '1C'], :]
    c1f = only_com.loc[ [i for i in only_com.index if i[-2:] == '1D'], :]
    c1g = only_com.loc[ [i for i in only_com.index if i[-2:] == '1E'], :]

    com_data = [c1a.to_numpy().ravel(), c1b.to_numpy().ravel(), c1c.to_numpy().ravel(), c1f.to_numpy().ravel(), c1g.to_numpy().ravel()]

    print('********* ANDERSON-DARLING TEST FOR UTS DATA ********************')
    for product in zip(ts_data, ['1A', '1B', '1C', '1D', '1E']):
        # print(product[0])
        test_stat, critical_value, significance_levels = stats.anderson(product[0], dist = 'norm')
        print(f'Product {product[1]} UTS A value = {test_stat}')
        print(f'Product {product[1]} UTS critical value = {critical_value}')
        print(f'significance levels: {significance_levels}')
        print()

    print('********* ANDERSON-DARLING TEST FOR COMMERCIAL DATA ********************')
    for product in zip(com_data, ['1A', '1B', '1C', '1D', '1E']):
        test_stat, critical_value, significance_levels = stats.anderson(product[0], dist = 'norm')
        print(f'Product {product[1]} commerical A value = {test_stat}')
        print(f'Product {product[1]} commercial critical value = {critical_value}')
        print(f'significance levels: {significance_levels}')
        print()

    print('********* LEVENES TEST UTS VS. COMMERCIAL ********************')
    for i, product in enumerate(['1A', '1B', '1C', '1D', '1E']):
        test_stat, p_value = stats.levene(ts_data[i], com_data[i])
        print(f'Product {product} Levene value for UTS vs commercial = {test_stat}')
        print(f'Product {product} p value = {p_value}')
        # print(significance_levels)
        print()

    print('********* ONE-TAIL F TEST UTS VS. COMMERCIAL ********************')
    for i, product in enumerate(['1A', '1B', '1C', '1D', '1E']):
        f = np.var(ts_data[i], ddof=1)/np.var(com_data[i], ddof=1)
        nun = ts_data[i].shape[0] - 1
        dun = com_data[i].shape[0] - 1
        p_value = 1- stats.f.cdf(f, nun, dun)
        print(f'F test for {product}')
        print(f'critical value {stats.f.ppf(q=1-.05, dfn=nun, dfd=dun)}')
        print(f'f: {f}')
        print(f'TS more varied than commercial? {f > stats.f.ppf(q=1-.10, dfn=nun, dfd=dun)}')
        print(f'p value: {p_value}')

