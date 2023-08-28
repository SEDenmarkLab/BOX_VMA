import pandas as pd
import re
import os

in_dir = 'RFERidge_reduced_feats_subs_dependent/'
out_dir = 'aggregated_features/'
# feature_range = [i for i in range(4,40)]
feature_range = [26]


if __name__ == '__main__':

    # feature spaces
    modfeatures = pd.read_csv('preprocess_x/aso_aeif_for_modelling_nocorr95.csv', header = 0, index_col = 0)
    # allfeatures = pd.read_csv('preprocess_x/all_aso_aeif_nocorr95.csv', header = 0, index_col = 0)

    for f in feature_range:
        filename = in_dir + str(f) + '_feats/kept_feats_log.csv'

        data = pd.read_csv(filename, index_col = 0, header = None).astype(str)

        # tally the votes
        counter = {}
        for vote in data.values.ravel():
            # skip over substrate one hot descriptors
            if ''.join(re.findall("[a-zA-Z]+", vote)) == '':
                continue
            # vote
            if vote in counter:
                counter[vote] = counter[vote] + 1
            else: 
                counter[vote] = 1


        # sort by the most votes
        counter = pd.DataFrame(counter, index = [str(f) + 'D feature space']).transpose()
        counter.sort_values([str(f) + 'D feature space'], axis = 0, inplace = True, ascending = False)
        
  
        if not os.path.exists(out_dir):
            os.makedirs(f'{out_dir}/for_modelling/')
            os.makedirs(f'{out_dir}/all/')
        
        # do modelling features
        for i in range(4, 21, 2):
            selected_features = modfeatures.loc[:, counter.index[0:i]]
            print(selected_features)
            
            with open (out_dir + f'for_modelling/modelling_top{i}_feats_' + str(f) + 'D.csv', 'w') as out:
                selected_features.to_csv(out)


        # do whole in silico library features
        
        # selected_features = allfeatures.loc[:, counter.index[0:20]]
        # print(selected_features)
        
        # with open (out_dir + 'all/all_top20_feats_' + str(f) + 'D.csv', 'w') as out:
        #     selected_features.to_csv(out)
        
