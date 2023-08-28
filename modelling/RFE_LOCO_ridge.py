import pandas as pd
import os
from sklearn.feature_selection import RFE

from sklearn.linear_model import  RidgeCV

from sklearn.base import clone
from multiprocessing import Pool
import pickle

def my_read_csv(file):
    rtn = pd.read_csv(file, header = None)
    rtn.set_index(0, inplace=True)
    rtn.columns = pd.RangeIndex(stop = len(rtn.columns))
    return(rtn)

# this function expands the X labels to include the substrates, without adding any additional columns
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


# pool handler does hard work - see the comments in substrate-dependent RFE for more explanation.
def pool_handler(args):
    

    splits_dir = args[0]
    out_dir = args[1]
    train_file = args[2]
    X_data_file = args[3]
    feature_reduction_range = args[4]
    eliminator = args[5]
    Yinverse = args[6]

    
    splitID = '_'.join(train_file.split('_')[1:]).split('.')[0]
    

    X_data = pd.read_csv(X_data_file, header = 0, index_col = 0)

    # we will copy the descriptors, but not add any substrate descriptors yet.
    X_data = append_subs_labels_copy_desc( X_data )



    Y_data_train = my_read_csv(splits_dir +'train/' + train_file)
    X_data_train = X_data.loc[[l for l in Y_data_train.index], :]


    X_data_train = X_data.loc[[l for l in Y_data_train.index], :]
    Y_test = my_read_csv(splits_dir +'test/test_' + str(splitID) + '.csv')
    X_test = X_data.loc[ [l for l in Y_test.index], : ]

    
    # separate the data into individual product clusters
    Xtr_1A = X_data_train.loc[[i for i in Y_data_train.index if i.split('_')[-1] == '1A' ], :]
    Xtr_1B = X_data_train.loc[[i for i in Y_data_train.index if i.split('_')[-1] == '1B' ], :]
    Xtr_1C = X_data_train.loc[[i for i in Y_data_train.index if i.split('_')[-1] == '1C' ], :]
    Xtr_1F = X_data_train.loc[[i for i in Y_data_train.index if i.split('_')[-1] == '1D' ], :]
    Xtr_1G = X_data_train.loc[[i for i in Y_data_train.index if i.split('_')[-1] == '1E' ], :]
    
    Ytr_1A = Y_data_train.loc[[i for i in Y_data_train.index if i.split('_')[-1] == '1A' ], :]
    Ytr_1B = Y_data_train.loc[[i for i in Y_data_train.index if i.split('_')[-1] == '1B' ], :]
    Ytr_1C = Y_data_train.loc[[i for i in Y_data_train.index if i.split('_')[-1] == '1C' ], :]
    Ytr_1F = Y_data_train.loc[[i for i in Y_data_train.index if i.split('_')[-1] == '1D' ], :]
    Ytr_1G = Y_data_train.loc[[i for i in Y_data_train.index if i.split('_')[-1] == '1E' ], :]

    Xte_1A = X_test.loc[[i for i in X_test.index if i.split('_')[-1] == '1A' ], :]
    Xte_1B = X_test.loc[[i for i in X_test.index if i.split('_')[-1] == '1B' ], :]
    Xte_1C = X_test.loc[[i for i in X_test.index if i.split('_')[-1] == '1C' ], :]
    Xte_1F = X_test.loc[[i for i in X_test.index if i.split('_')[-1] == '1D' ], :]
    Xte_1G = X_test.loc[[i for i in X_test.index if i.split('_')[-1] == '1E' ], :]
    
    Yte_1A = Y_test.loc[[i for i in Y_test.index if i.split('_')[-1] == '1A' ], :]
    Yte_1B = Y_test.loc[[i for i in Y_test.index if i.split('_')[-1] == '1B' ], :]
    Yte_1C = Y_test.loc[[i for i in Y_test.index if i.split('_')[-1] == '1C' ], :]
    Yte_1F = Y_test.loc[[i for i in Y_test.index if i.split('_')[-1] == '1D' ], :]
    Yte_1G = Y_test.loc[[i for i in Y_test.index if i.split('_')[-1] == '1E' ], :]

    clusters = [(Xtr_1A, Ytr_1A, '1A', Xte_1A, Yte_1A), (Xtr_1B, Ytr_1B, '1B', Xte_1B, Yte_1B),
                 (Xtr_1C, Ytr_1C, '1C', Xte_1C, Yte_1C), (Xtr_1F, Ytr_1F, '1D', Xte_1F, Yte_1F), 
                 (Xtr_1G, Ytr_1G, '1E', Xte_1G, Yte_1G)]

    
    # process features.
    for f in feature_reduction_range:
        if not os.path.exists(out_dir + str(f) + '_feats/train/' ):
            os.makedirs(out_dir + str(f) + '_feats/train/')
            os.makedirs(out_dir + str(f) + '_feats/test/')
        
        # we will save out some useful files
        with open(out_dir + str(f) + '_feats/featslog' +'.csv', 'a')  as featfile:
            with open(out_dir + str(f) + '_feats/performancelog_total' +'.csv', 'a') as perflog:
                all_kept_features = []
                featfile.write(f'for split: {splitID}\n')

                # process each cluster individually
                for cluster in clusters:
                    subs_eliminator = clone(eliminator)
                    subs_eliminator.set_params(n_features_to_select = f)
                    subs_eliminator.fit(cluster[0].to_numpy(), cluster[1].to_numpy().ravel())

                    # save out test statistics
                    with open(out_dir + str(f) + f'_feats/performancelog_{cluster[2]}' +'.csv', 'a') as perflogsubs:
                        if len(cluster[3]) > 0:
                            yhat = subs_eliminator.predict(cluster[3].to_numpy())
                            real_test_y = Yinverse.inverse_transform(cluster[4].to_numpy().reshape(len(cluster[4].to_numpy()), 1))
                            real_test_y_hat = Yinverse.inverse_transform(yhat.reshape(len(yhat), 1))
                            perflog.write(f'{splitID}_{cluster[2]}, {real_test_y[0][0]}, {real_test_y_hat[0][0]},\n'  )
                            perflogsubs.write(f'{splitID}_{cluster[2]}, {real_test_y[0][0]}, {real_test_y_hat[0][0]},\n')


                    # write out the features we kept for future reference
                    kept_features = subs_eliminator.get_support(indices = True)
                    featfile.write(str(cluster[2]) + '/' + str(f) + 'feats (kept features indices),' + str(list(kept_features))[1:-1] +'\n'  )
                    all_kept_features.extend(list(kept_features))
                
                all_kept_features = list(set(all_kept_features))
                featfile.write('total features from all substrates: ' + str(len(all_kept_features)) + '\n')
                featfile.write('\n')

                # reduce features
                X_train_red = X_data_train.iloc[:, all_kept_features]
                X_test_red = X_test.iloc[:, all_kept_features]

                # save
                with open(out_dir + str(f) + '_feats/train/train_' + str(splitID) + '.csv', 'wt', newline = '') as o:
                    X_train_red.to_csv(o)
                with open(out_dir + str(f) + '_feats/test/test_' + str(splitID) + '.csv', 'wt', newline = '') as o:
                    X_test_red.to_csv(o)


if __name__ == '__main__':
    


    ''' Here's the section of directories, models, and feature redeuction methods you need to specify.'''

    splits_dir = 'preprocess_y/LOCO_splits/' # a directory containing train/ and test/ subdirectories with transformed labelled Y data for splits. Data csv files are formatted as 'train_<ID>' where ID is split index as integer
    X_data_file = 'preprocess_x/aso_aeif_for_modelling_nocorr95.csv' # all x features for each catalyst. Only catalyst descriptors. NOTE: MODIFY CODE TO ADD SUBSTRATE DESCRIPTORS. SEE BELOW.
    # this is just to get the inverse tranformer (from full data set, not train/test splits) for later. The train/test splits are already transformed (I do it with yeo-johnson).
    with open('preprocess_y/transformer.pkl', 'rb') as t:
        Ytransform = pickle.load(t)
    out_dir = 'RFERidge_reduced_feats/' # directory where you want the models written. This directory will be made if it does not exist.
    
    feature_reduction_range = [10]
  
    feature_eliminator = RFE(RidgeCV())

    

    args = [ (splits_dir, out_dir, f, X_data_file, feature_reduction_range, clone(feature_eliminator), Ytransform) for f in os.listdir(splits_dir + 'train/') ]

    # =======================================================================================
    #       THIS IS VERY IMPORTANT
    #       By default sklearn relies on environment variables to do shit in parallel
    #       Those conflict with a lot of quantum chemistry installations
    #       such as XTB, CREST and like
    #       THIS OVERRIDES RANDOM BULLSHIT VARIABLES
    #       MAKE SURE THAT MKL_NUM_THREADS * POOL_WORKERS <= 128
    # =======================================================================================
    os.environ['MKL_NUM_THREADS'] = "4"
    os.environ['OMP_NUM_THREADS'] = "4"

    with Pool(processes=32) as p:
        p.map(pool_handler, args)

    print('Success!')


    

    
