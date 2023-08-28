import pandas as pd
import numpy as np
import os

from sklearn.svm import SVR, NuSVR

import pickle

subs_one_hot = {'1A': [1, 0 , 0, 0,0],
        '1B': [0, 1, 0, 0,0],
        '1C': [0,0,1,0,0],
        '1D': [0,0,0,1,0],
        '1E' : [0,0,0,0,1]
}

nns = ['30_3_1_11',
        '48_2_1_10',
        '124_2_1_11',
        '77_1_2_11',
        '189_1_2_11']
nns_bridge = [
        '77_1_2_12',
        '189_1_2_12'] 

prediction = ['188_2_1_12']


def append_one_hot_subs(X):

    onehot = pd.DataFrame(subs_one_hot).transpose()
    onehot.columns = [-1, -2, -3, -4, -5]

    onehot = pd.DataFrame( onehot.loc[ [ i.split('_')[-1] for i in X.index ], : ].to_numpy(), index = X.index, columns = onehot.columns )
  
    rtn = pd.concat([X, onehot  ], axis = 1)
    return rtn

def append_subs_labels_copy_desc(X):
    rtn_index = []
    rtn = []
    for row in X.index:
        
        label_1a = row + '_1A'
        label_1b = row + '_1B'
        label_1c = row + '_1C'
        label_1f = row + '_1D'
        label_1g = row + '_1E'
        rtn_index = rtn_index + [ label_1a, label_1b, label_1c, label_1f, label_1g]
        
        feature_1a = list(X.loc[row, :])
        rtn.append(feature_1a)
        feature_1b = list(X.loc[row, :])
        rtn.append(feature_1b)
        feature_1c = list(X.loc[row, :])
        rtn.append(feature_1c)
        feature_1f = list(X.loc[row, :])
        rtn.append(feature_1f)
        feature_1g = list(X.loc[row, :])
        rtn.append(feature_1g)

    return pd.DataFrame(rtn, index = rtn_index, columns = X.columns)


def preprocess_y(in_data):
    
    
    # read in the y data
    data = pd.read_csv(in_data, header = 0, index_col=0)
    # get rid of the dr column
    data.drop(['dr'], axis = 1, inplace = True)

    # calculate er
    data['er'] = (1.0 + (data['ee']/100.0)) / (1.0 - (data['ee']/100.0))
    # just storing this temporarily
    data['pph of desired enantiomer'] = 100.0 / (data['er'] + 1)
    data['pph of desired enantiomer'] = data['er']*data['pph of desired enantiomer']

    # convert to ddg, -20C = 253 K, R = 0.0019872 kcal/molK
    data['delta delta G (kcal/mol)'] = np.log(data['er']) * 253.15 * 0.0019872

    
    return data

if __name__ == '__main__':
    


    ''' Here's the section of directories, models, and feature redeuction methods you need to specify.'''
    with open('preprocess_y/transformer.pkl', 'rb') as o:
        Ytransform = pickle.load(o)
    all_in_silico_X_file = 'aggregated_features/all/all_top20_feats_26D.csv'

    os.environ['MKL_NUM_THREADS'] = "4"
    os.environ['OMP_NUM_THREADS'] = "4"


    model = None
    with open('out_test_aggregated_features/SVR/modelling_top20_feats_26D.csv/SVR.pkl', 'rb') as o:
        model = pickle.load(o)

    X_data = pd.read_csv(all_in_silico_X_file, index_col=0, header=0)
    print(X_data)

    X_data = X_data.loc[nns + nns_bridge + prediction, :]
    print(X_data)

    X_data = append_subs_labels_copy_desc(X_data)
    X_data = append_one_hot_subs(X_data)

    print(X_data)


    yhat = model.predict(X_data.to_numpy())
    real_test_y = pd.DataFrame( Ytransform.inverse_transform(yhat.reshape(len(yhat), 1)), index = X_data.index)


    real_y = pd.read_csv('../../ee_data/data/processed_data_fa_and_nn_1a-1e.csv', index_col = 0, header = 0)
    real_y = real_y.loc[real_test_y.index, ['ee']]
    print(real_y)

    real_y_processed = preprocess_y('../../ee_data/data/processed_data_fa_and_nn_1a-1e.csv')
    
    real_y_processed = real_y_processed.loc[real_test_y.index, 'delta delta G (kcal/mol)']
    print(real_y_processed)

    out_df = pd.concat([real_y_processed, real_test_y], axis = 1)
    out_df.columns = ['Y', 'Yhat']
    out_df['absolute error'] = np.abs(out_df['Y'] - out_df['Yhat'])

    out_df['er_pred'] = np.exp((-1 * out_df['Yhat']) / (253.15 * 0.0019872))
    out_df['ee_pred'] = -1* (100*out_df['er_pred'] - 100) / (1 + out_df['er_pred'])

    out_df['er'] = np.exp((-1 * out_df['Y']) / (253.15 * 0.0019872))
    out_df['ee'] = -1* (100*out_df['er'] - 100) / (1 + out_df['er'])

    out_df['absolute error ee'] = np.abs(out_df['ee'] - out_df['ee_pred'])

    print(out_df)

    print(pd.concat([real_y['ee'], out_df['ee']], axis = 1))
    if not os.path.exists('predictions/'):
        os.makedirs('predictions')

    with open('predictions/nn_and_prediction.csv', 'w') as o:
        out_df.to_csv(o )



    
