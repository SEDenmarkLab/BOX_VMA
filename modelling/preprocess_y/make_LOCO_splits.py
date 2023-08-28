import pandas as pd
from sklearn.preprocessing import PowerTransformer
import pickle
import os
import numpy as np

ee_data_sheet = '../../ee_data/data/processed_data_initial_screening_1a-1e.csv'
i = 'y_data_ee.csv'
o_train = 'LOCO_splits/train/'
o_test = 'LOCO_splits/test/'


# read in the y data
data = pd.read_csv(ee_data_sheet, header = 0, index_col=0)
# get rid of the dr column
data.drop(['dr'], axis = 1, inplace = True)

# calculate er
data['er'] = (1.0 + (data['ee']/100.0)) / (1.0 - (data['ee']/100.0))
# just storing this temporarily
data['pph of desired enantiomer'] = 100.0 / (data['er'] + 1)
data['pph of desired enantiomer'] = data['er']*data['pph of desired enantiomer']

# convert to ddg, -20C = 253 K, R = 0.0019872 kcal/molK
data['delta delta G (kcal/mol)'] = np.log(data['er']) * 253.15 * 0.0019872

# save summary file
with open('y_data_summary.csv', 'wb') as o:
    data.to_csv(o)

# get just our ddg values
all_data = data.loc[:, ['delta delta G (kcal/mol)']]

# transform with Yeo-Johnson (sklearn default)
transformer = PowerTransformer()
all_data = pd.DataFrame(transformer.fit_transform(all_data.to_numpy()), index = all_data.index)

# save the transformer for later
with open('transformer.pkl', 'wb') as o:
    pickle.dump(transformer, o)

# all of the catlyst labels
cats = list(set(['_'.join(label.split('_')[:-1]) for label in all_data.index]))


# # make LOCO train/test splits
for cat in cats:
    i = cat

    test_set = all_data.loc[[label for label in all_data.index if label[:-3] == cat], :]
    train_set = all_data.loc[[label for label in all_data.index if label[:-3] != cat], :]

    # save out
    if not os.path.exists(o_train):
        os.makedirs(o_train)
        os.makedirs(o_test)

    with open(o_train + 'train_' + str(i) + '.csv', 'wt', newline = '') as otr:
        train_set.to_csv(otr, header = None)

    with open(o_test + 'test_' + str(i) + '.csv', 'wt', newline = '') as ote:
        test_set.to_csv(ote, header = None)