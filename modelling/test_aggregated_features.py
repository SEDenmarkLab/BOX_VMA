import pandas as pd
import numpy as np
import os
from csv import reader, writer
from sklearn.preprocessing import MinMaxScaler, PowerTransformer
from sklearn.feature_selection import VarianceThreshold, SelectFromModel, RFE
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from random import randint
import random
from sklearn.linear_model import LassoCV, RidgeCV
from sklearn.svm import SVR, NuSVR
from sklearn.model_selection import GridSearchCV, RandomizedSearchCV
from scipy.stats import loguniform, uniform
import shutil
import subprocess
from sklearn.base import clone
import pickle
subs_one_hot = {'1A': [1, 0 , 0, 0,0],
        '1B': [ 0, 1, 0, 0,0],
        '1C': [0,0,1,0,0],
        '1D': [0,0,0,1,0],
        '1E' : [0,0,0,0,1]
}


def my_read_csv(file):
    rtn = pd.read_csv(file, header = None)
    rtn.set_index(0, inplace=True)
    rtn.columns = pd.RangeIndex(stop = len(rtn.columns))
    return(rtn)



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

# FOR ONE-HOT CONTROL
def onehot_features(X):
    rtn = []
    for i in range(len(X.index)):
        lis = [0] * len(X.index)
        lis[i] = 1
        rtn.append(lis)
    return pd.DataFrame(rtn, index = X.index)

# FOR RANDOM CONTROL
def randomize_features(X):
    minmax = MinMaxScaler().fit(X.to_numpy())
    min_vector = minmax.data_min_
    max_vector = minmax.data_max_
    random_feats = np.zeros(X.shape)

    # sample uniform random between min max of each column
    for coli in range(random_feats.shape[1]):
        random_feats[:, coli] = np.random.uniform(min_vector[coli], max_vector[coli], random_feats.shape[0])

    return pd.DataFrame(random_feats, index = X.index)



if __name__ == '__main__':
    


    ''' Here's the section of directories, models, and feature redeuction methods you need to specify.'''

    splits_dir = 'preprocess_y/LOCO_splits/' # a directory containing train/ and test/ subdirectories with transformed labelled Y data for splits. Data csv files are formatted as 'train_<ID>' where ID is split index as integer
    X_data_dir = 'aggregated_features/for_modelling/' # all x features for each catalyst. Only catalyst descriptors. NOTE: MODIFY CODE TO ADD SUBSTRATE DESCRIPTORS. SEE BELOW.
    # this is just to get the inverse tranformer (from full data set, not train/test splits) for later. The train/test splits are already transformed (I do it with yeo-johnson).
    with open('preprocess_y/transformer.pkl', 'rb') as t:
        Ytransform = pickle.load(t)
    out_dir = 'out_test_aggregated_features/' # directory where you want the models written. This directory will be made if it does not exist.

    # some potentially useful hyperparameter domains
    params1 = {
        'learning_rate': [0.01, 0.1, 1.0],
        'n_estimators' : [10, 100, 1000],
        'subsample': [0.9, 1.0]
        }
    distributions = {
        'C' : loguniform(1e-4, 10.0),
        'epsilon' : loguniform(1e-4, 10.0),
        'kernel' : ('linear', 'rbf'),
        'gamma' : ('auto', 'scale')
        }
    distributions2 = {
        'C' : loguniform(1e-3, 10.0),
        'nu' : uniform(1e-4, 1.0),
        'kernel' : ('linear', 'rbf'),
        'gamma' : ('auto', 'scale')
        }
    lassoparams = {
        'alpha' : uniform(loc = 0, scale = 5)
    }
    random_grid = {'n_estimators': [int(x) for x in np.linspace(start = 200, stop = 2000, num = 20)],
               'max_features': ['auto', 'sqrt'],
               'max_depth': [int(x) for x in np.linspace(10, 110, num = 11)] + [None],
               'min_samples_split': [2, 5, 10],
               'min_samples_leaf': [1, 2, 4],
               'bootstrap': [True, False]}

    
    # models contains all of the models that we want to test.
    models = {
                 'SVR' : RandomizedSearchCV(estimator = SVR(), param_distributions = distributions, random_state = 1, n_iter = 100 ),
                'nuSVR' : RandomizedSearchCV(estimator = NuSVR(),  param_distributions = distributions2, random_state = 1, n_iter = 100),
                # 'LassoCV' : LassoCV(max_iter = 100000, random_state = 1),
                # 'PLS1' : PLSRegression(1),
                'PLS2' : PLSRegression(2),
                'PLS3' : PLSRegression(3),
                # 'ElasticNet' : ElasticNet(random_state=1),
                'PLS4' : PLSRegression(4),
                # 'PLS5' : PLSRegression(5),
                # 'PLS5' : PLSRegression(7),
                # 'PLS5' : PLSRegression(9),
                # 'DecisionTree' : DecisionTreeRegressor(random_state = 1),
                'LassoCV': LassoCV(max_iter = 1000, random_state = 1),
                'RidgeCV': RidgeCV(),
                # 'GaussianProcess': GaussianProcessRegressor(random_state=1),
                'RandomForest': RandomForestRegressor(n_estimators = 100, random_state = 1),
                'GBR': GradientBoostingRegressor(n_estimators = 100, random_state = 1)

              
        }

    ''' end main user specified bock. There are some things you may need to change below depending on you modelling scheme '''
    
    # f are csv files of different ligand feature spaces
    for f in os.listdir(X_data_dir):
        # start iterating through train/test splits
        for train_file in os.listdir(splits_dir + 'train/'):

                # the train/test split files are formatted as 'train_x.csv' where x is splitID
            splitID = '_'.join(train_file.split('_')[1:]).split('.')[0]

            Y_data_train = my_read_csv(splits_dir +'train/' + train_file)
            Y_test = my_read_csv(splits_dir +'test/test_' + str(splitID) + '.csv')

            X_data = pd.read_csv(X_data_dir + f, header = 0, index_col = 0)

            # CHANGE THIS PART FOR RANDOM/ONE HOT CONTROLS
            X_data = append_subs_labels_copy_desc(X_data)
            X_data = append_one_hot_subs(X_data)

            X_data_train = X_data.loc[[i for i in Y_data_train.index], :]
            X_test = X_data.loc[[i for i in Y_test.index], :]

                
            # interate through models
            for model in models.keys():
                
                # make dfs into numpy
                X = X_data_train.to_numpy()
                Y = Y_data_train.to_numpy().ravel()
                Xt = X_test.to_numpy()
                Yt = Y_test.to_numpy().ravel()
              
                # fit model, predict
                modelc = clone(models[model])
                modelc.fit(X,Y)

                y_hat_train = modelc.predict(X)
                y_hat_test = modelc.predict(Xt)
                # make our y's back into delta delta G's. comment this out if you don't care/ don't have an inverse_transform function.
                real_train_y = Ytransform.inverse_transform(Y.reshape(len(Y), 1))
                real_train_y_hat = Ytransform.inverse_transform(y_hat_train.reshape(len(y_hat_train), 1))
                real_test_y = Ytransform.inverse_transform(Yt.reshape(len(Yt), 1))
                real_test_y_hat = Ytransform.inverse_transform(y_hat_test.reshape(len(y_hat_test), 1))

                # add the predictions to the out data
                train_y_df_out = pd.DataFrame(np.append(real_train_y, real_train_y_hat, axis=1), index = Y_data_train.index) 
                test_y_df_out = pd.DataFrame(np.append(real_test_y, real_test_y_hat, axis = 1), index = Y_test.index ) 

                #make directories if we need them.
                if not os.path.exists(out_dir + model + '/' + str(f) +'/train'):
                    try:
                        os.makedirs(out_dir + model + '/' + str(f) +'/train')
                    except FileExistsError:
                        pass
                if not os.path.exists(out_dir + model + '/' + str(f) +'/test'):
                    try:
                        os.makedirs(out_dir + model + '/' + str(f) +'/test')
                    except FileExistsError:
                        pass

                # write train/test splits out
                with open(out_dir + model + '/' + str(f) +'/train/train_' + str(splitID) +'.csv', 'wt', newline = '') as out_train:
                    train_y_df_out.to_csv(out_train, header = [ 'Y', 'Y hat'])
                with open(out_dir + model + '/' + str(f) +'/test/test_' + str(splitID) +'.csv', 'wt', newline = '') as out_test:
                    test_y_df_out.to_csv(out_test, header = ['Y', 'Y hat'])
                
                # write out model
                with open(f'{out_dir}{model}/{f}/{model}.pkl', 'wb') as o:
                    pickle.dump(modelc, o)


    
