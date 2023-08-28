import pandas as pd
import numpy as np
import os
import pickle
from sklearn.feature_selection import RFE
from sklearn.cross_decomposition import PLSRegression
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import LassoCV, RidgeCV
from sklearn.svm import SVR, NuSVR
from sklearn.model_selection import RandomizedSearchCV
from scipy.stats import loguniform, uniform
from multiprocessing import Pool

subs_one_hot = {'1A': [1, 0 , 0, 0,0],
        '1B': [ 0, 1, 0, 0,0],
        '1C': [0,0,1,0,0],
        '1D': [0,0,0,1,0],
        '1E' : [0,0,0,0,1]
}

# a stupid function from back when I was clunky with pandas
def my_read_csv(file):
    rtn = pd.read_csv(file, header = None)
    rtn.set_index(0, inplace=True)
    rtn.columns = pd.RangeIndex(stop = len(rtn.columns))
    return(rtn)


# appends one hot descriptors to a ligand_substrate labelled df
def append_one_hot_subs(X):
    onehot = pd.DataFrame(subs_one_hot).transpose()
    onehot.columns = [-1, -2, -3, -4, -5]
    onehot = pd.DataFrame( onehot.loc[ [ i.split('_')[-1] for i in X.index ], : ].to_numpy(), index = X.index, columns = onehot.columns )
    rtn = pd.concat([X, onehot  ], axis = 1)
    return rtn



# pool handler does the hard work
def pool_handler(args):
    
    # the train/test split files are formatted as 'train_x.csv' where x is splitID
    splits_dir = args[0]
    # where we save files
    out_dir = args[1]
    # the training file
    train_file = args[2]
    # the directory where our descriptors are
    x_in_dir = args[3]
    # dictionary of models to test
    models = args[4]
    # number of feats to select
    feat = args[5]
    # the Ytransform, for inverse transform later
    Ytransform = args[6]
    

    # the train/test split files are formatted as 'train_x.csv' where x is splitID
    splitID = '_'.join(train_file.split('_')[1:]).split('.')[0]

    
    # read in the data
    Y_data_train = my_read_csv(splits_dir +'train/' + train_file)
    X_data_train = pd.read_csv(x_in_dir + 'train/train_' + str(splitID) + '.csv', header = 0, index_col= 0)
    Y_test = my_read_csv(splits_dir +'test/test_' + str(splitID) + '.csv')
    X_test = pd.read_csv(x_in_dir + 'test/test_' + str(splitID) + '.csv', header = 0, index_col= 0)

    # add in the one hot descriptors
    X_data_train = append_one_hot_subs(X_data_train)
    X_test = append_one_hot_subs(X_test)
    
    # eliminate features
    eliminator = RFE(RidgeCV(), n_features_to_select = feat)
    eliminator.fit(X_data_train.to_numpy(), Y_data_train.to_numpy().ravel())
    kept_features = eliminator.get_support(indices = True)
    # get the kept feature names
    kept_features = [X_data_train.columns[i] for i in kept_features] 
    
    # save the feature names we kept for later
    kept_featuresmsg = 'kept features for split ' +str(splitID) + ',' + ','.join([str(i) for i in kept_features]) +'\n'  
    if not os.path.exists(out_dir + str(feat) + '_feats/'):
        try:
            os.makedirs(out_dir + str(feat) + '_feats/')
        except FileExistsError:
            pass
    
    # this is opened for appending so all of the pool workers can add to it.
    with open(out_dir + str(feat) + '_feats/kept_feats_log.csv', 'a') as logfile:
        logfile.write(kept_featuresmsg)
    
    # feature selection
    X_train_red = X_data_train.loc[:, [i for i in X_data_train.columns if i in kept_features ] ]
    X_test_red = X_test.loc[:, [i for i in X_test.columns if  i in kept_features ]  ]

    # iterate through models
    for model in models.keys():

        train_y_df_out = pd.DataFrame()
        test_y_df_out = pd.DataFrame()


        Xtr_df = X_train_red
        Ytr_df = Y_data_train
        Xte_df = X_test_red
        Yte_df = Y_test

        
        # fit model, predict
        models[model].fit(Xtr_df.to_numpy(),Ytr_df.to_numpy().ravel())
        y_hat_train = models[model].predict(Xtr_df.to_numpy())
        y_hat_test = models[model].predict(Xte_df.to_numpy())
        # make our y's back into delta delta G's. comment this out if you don't care/ don't have an inverse_transform function.
        real_train_y = Ytransform.inverse_transform(Ytr_df.to_numpy().reshape(len(Ytr_df.to_numpy()), 1))
        real_train_y_hat = Ytransform.inverse_transform(y_hat_train.reshape(len(y_hat_train), 1))
        real_test_y = Ytransform.inverse_transform(Yte_df.to_numpy().reshape(len(Yte_df.to_numpy()), 1))
        real_test_y_hat = Ytransform.inverse_transform(y_hat_test.reshape(len(y_hat_test), 1))

        # add the predictions from this cluseter to the out data
        train_y_df_out = pd.DataFrame(np.append(real_train_y, real_train_y_hat, axis=1), index = Ytr_df.index)
        test_y_df_out = pd.DataFrame(np.append(real_test_y, real_test_y_hat, axis = 1), index = Yte_df.index )


        #make directories if we need them.
        if not os.path.exists(out_dir + str(feat) + '_feats/' + model +  '/train'):
            try:
                os.makedirs(out_dir + str(feat) + '_feats/' + model +  '/train')
            except FileExistsError:
                pass
        if not os.path.exists(out_dir + str(feat) + '_feats/' + model +  '/test'):
            try:
                os.makedirs(out_dir + str(feat) + '_feats/' + model +  '/test')
            except FileExistsError:
                pass
           

        # write train/test results out
        with open(out_dir + str(feat) + '_feats/' + model +  '/train/train_' + str(splitID) +'.csv', 'wt', newline = '') as out_train:
            train_y_df_out.to_csv(out_train, header = [ 'Y', 'Y hat'])
        with open(out_dir + str(feat) + '_feats/' + model +  '/test/test_' + str(splitID) +'.csv', 'wt', newline = '') as out_test:
            test_y_df_out.to_csv(out_test, header = ['Y', 'Y hat'])




if __name__ == '__main__':
    


    #t he section of directories, models, and feature redeuction methods you need to specify.
    splits_dir = 'preprocess_y/LOCO_splits/'
    x_in_dir = 'RFERidge_reduced_feats/10_feats/' # a directory containing train/ and test/ subdirectories with transformed labelled Y data for splits. Data csv files are formatted as 'train_<ID>' where ID is split index as integer
    # this is just to get the inverse tranformer (from full data set, not train/test splits) for later. The train/test splits are already transformed (I do it with yeo-johnson).
    with open('preprocess_y/transformer.pkl', 'rb') as t:
        Ytransform = pickle.load(t)
    # directory where you want the models written. This directory will be made if it does not exist.
    out_dir = 'RFERidge_reduced_feats_subs_dependent/' 

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
               'max_depth': [int(x) for x in np.linspace(10, 110, num = 11)].append(None),
               'min_samples_split': [2, 5, 10],
               'min_samples_leaf': [1, 2, 4],
               'bootstrap': [True, False]}

    # specify iterable of integers for number of reduced features to test for feature reduction/elimination. Models will be made for each reduced space.
    feature_reduction_range = [i for i in range (4,40)] 

    # define the hyperparameter domains for our models
    # models contains all of the models that we want to test.
    models = {
                'SVR' : RandomizedSearchCV(estimator = SVR(), param_distributions = distributions, random_state = 1, n_iter = 1000 ),
                'nuSVR' : RandomizedSearchCV(estimator = NuSVR(),  param_distributions = distributions2, random_state = 1, n_iter = 1000),
                'PLS2' : PLSRegression(2),
                'PLS3' : PLSRegression(3),
                'PLS4' : PLSRegression(4),
                'LassoCV': LassoCV(max_iter = 100000, random_state = 1),
                'RidgeCV': RidgeCV(),
                'RandomForest': RandomForestRegressor(n_estimators = 1000, random_state = 1),
                'GBR': GradientBoostingRegressor(n_estimators = 1000, random_state = 1)
              
        }

    os.environ['MKL_NUM_THREADS'] = "4"
    os.environ['OMP_NUM_THREADS'] = "4"

    # iterate through our feature reduction range and process splits
    for feat in feature_reduction_range:
        args = [ (splits_dir, out_dir, f, x_in_dir, models, feat, Ytransform) for f in os.listdir(splits_dir + 'train/') ]
        with Pool(processes=32) as p:
            p.map(pool_handler, args)
    
    print('Sucess!')
