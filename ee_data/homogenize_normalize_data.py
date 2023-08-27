from csv import reader, writer
import os
from string import whitespace
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
from sklearn.preprocessing import PowerTransformer
import re
import os
os.environ['OMP_NUM_THREADS'] = '8'


label_key = 'catalyst_label_key.csv'

product_key = {
    'PDT1' : '2B',
    'PDT6' : '1A',
    'PDT4' : '1B',
    'PDT5' : '1C',
    '1E' : '1E',
    '1D' : '1D'
}

# a function to make it look nice for the paper
def pretty_process_initial_data():
    inp = 'data/unprocessed_data_initial_screening.csv'
    output = 'data/pretty_unprocessed_data_initial_screening_1a-1e.csv'
    # read in the label key. It has three columns, Merck's label, our labs label, and the chemistry at the 4-position of the BOX ligand.
    label_key_data = pd.read_csv(label_key, header=0, index_col=0)
    print(label_key_data)
     # read in the unprocessed data.
    data = pd.read_csv(inp, header=0)
    print(data)
    # exit()

    # first make our ee column into strings
    data['ee'] = data['ee'].values.astype('str')
    # use some regex to throw out any characters that are not digits 1-9, decimal points, or minus signs
    data['ee'] = [re.sub('[^0-9.-]', '', i) for i in data['ee']]
    # replace our empty rows with nan  
    data['ee'].replace('', np.nan, inplace= True)
    # drop nan rows, now there are 566 rows
    data.dropna(axis=0, inplace=True)
    data['ee'] = data['ee'].values.astype('float')
    # do same for dr, but we won't drop the rows with nan
    data['dr'] = data['dr'].values.astype('str')
    data['dr'] = [re.sub('[^0-9.-]', '', i) for i in data['dr']]
    data['dr'].replace('', np.nan, inplace= True)
    data['dr'] = data['dr'].values.astype('float')
    print(data)
    # exit()

    # now we make the merck labels column into strings
    data['Merck label'] = data['Merck label'].values.astype('str')
    # remove those annoying eol characters - these are unnecessary anyway. also strip trailing whitespace
    data['Merck label'] = [i.split('\n')[0].strip() for i in data['Merck label']]
    print(data)

    # now assign all of the product labels correctly
    data['pdt'] = [product_key[i] for i in data['pdt']]
    print(data)

    # note the stereochemistry
    data['stereo'] = [label_key_data.loc[i, ['stereochemistry at 4 position']].values[0] for i in data['Merck label']]
    print(data)
    # assign the Denmark labels correctly
    data['Denmark label'] = [label_key_data.loc[i, ['Denmark label']].values[0] for i in data['Merck label']]
    print(data)
    # # apply the correction to the S,S catalysts - get rows with S stereochem
    # s_data = data.where(data['stereo'] == 'S')
    # print(s_data)
    # # invert the ee value
    # s_data['ee'] = s_data['ee'] * -1
    # print(s_data)
    # # keep the rows that are R unchanged, else, replace with the value in s data
    # print(data)
    # data = data.where(data['stereo'] == 'R', s_data)
    # print(data)

    # lets go ahead and drop the 2B data. These are mostly racemic, derived from non-methylated furan.
    print(data)
    data = data.where(data['pdt'] != '2B').dropna(axis = 0)
    print(data)

    # assign unique reaction handles for each one
    data['reaction handle'] = data['Denmark label'] + '_' + data['pdt']
    print(data)
    # exit()
        
    # now we average
    data_avg = data.groupby('reaction handle').mean()
    # pd.options.display.max_rows = 500
    print(data_avg)
    data_avg['cat'] =['_'.join(i.split('_')[:-1]) for i in data_avg.index]
    data_avg['prod'] =['_'.join(i.split('_')[-1:]) for i in data_avg.index]
    data_avg['prod'] =[i.lower() for i in data_avg['prod']]
    data_avg['prod'] =['3' + i[-1] for i in data_avg['prod']]
    # print(data)
    data_avg = data_avg.sort_values(by = 'cat', ascending=False)
    data_avg = data_avg.sort_values(by = 'prod', ascending=True)
    print(data_avg)
    # save the homogenized, corrected data
    with open(output, 'w') as o:
        data_avg.to_csv(o)
    
# a function to make it look nice for the paper
def pretty_process_nn_and_fa():
    inp = 'data/unprocessed_data_fa_and_nn.csv'
    output = 'data/pretty_unprocessed_data_fa_and_nn_1a-1e.csv'
    # read in the label key. It has three columns, Merck's label, our labs label, and the chemistry at the 4-position of the BOX ligand.
    label_key_data = pd.read_csv(label_key, header=0, index_col=0)
    print(label_key_data)
     # read in the unprocessed data.
    data = pd.read_csv(inp, header=0)
    print(data)
    # exit()

    # first make our ee column into strings
    data['ee'] = data['ee'].values.astype('str')
    # use some regex to throw out any characters that are not digits 1-9, decimal points, or minus signs
    data['ee'] = [re.sub('[^0-9.-]', '', i) for i in data['ee']]
    # replace our empty rows with nan  
    data['ee'].replace('', np.nan, inplace= True)
    # drop nan rows, now there are 566 rows
    data.dropna(axis=0, inplace=True)
    data['ee'] = data['ee'].values.astype('float')
    # do same for dr, but we won't drop the rows with nan
    data['dr'] = data['dr'].values.astype('str')
    data['dr'] = [re.sub('[^0-9.-]', '', i) for i in data['dr']]
    data['dr'].replace('', np.nan, inplace= True)
    data['dr'] = data['dr'].values.astype('float')
    print(data)
    # exit()

    # now we make the merck labels column into strings
    data['Merck label'] = data['Merck label'].values.astype('str')
    # remove those annoying eol characters - these are unnecessary anyway. also strip trailing whitespace
    data['Merck label'] = [i.split('\n')[0].strip() for i in data['Merck label']]
    print(data)

    # now assign all of the product labels correctly
    data['pdt'] = [product_key[i] for i in data['pdt']]
    print(data)

    # note the stereochemistry
    data['stereo'] = [label_key_data.loc[i, ['stereochemistry at 4 position']].values[0] for i in data['Merck label']]
    print(data)
    # assign the Denmark labels correctly
    data['Denmark label'] = [label_key_data.loc[i, ['Denmark label']].values[0] for i in data['Merck label']]
    print(data)
    # # apply the correction to the S,S catalysts - get rows with S stereochem
    # s_data = data.where(data['stereo'] == 'S')
    # print(s_data)
    # # invert the ee value
    # s_data['ee'] = s_data['ee'] * -1
    # print(s_data)
    # # keep the rows that are R unchanged, else, replace with the value in s data
    # print(data)
    # data = data.where(data['stereo'] == 'R', s_data)
    # print(data)

    # lets go ahead and drop the 2B data. These are mostly racemic, derived from non-methylated furan.
    print(data)
    data = data.where(data['pdt'] != '2B').dropna(axis = 0)
    print(data)

    # assign unique reaction handles for each one
    data['reaction handle'] = data['Denmark label'] + '_' + data['pdt']
    print(data)
    # exit()
        
    # now we average
    data_avg = data.groupby('reaction handle').mean()
    # pd.options.display.max_rows = 500
    print(data_avg)
    data_avg['cat'] =['_'.join(i.split('_')[:-1]) for i in data_avg.index]
    data_avg['prod'] =['_'.join(i.split('_')[-1:]) for i in data_avg.index]
    data_avg['prod'] =[i.lower() for i in data_avg['prod']]
    data_avg['prod'] =['3' + i[-1] for i in data_avg['prod']]
    # print(data)
    data_avg = data_avg.sort_values(by = 'cat', ascending=False)
    data_avg = data_avg.sort_values(by = 'prod', ascending=True)
    print(data_avg)
    # save the homogenized, corrected data
    with open(output, 'w') as o:
        data_avg.to_csv(o)

def process_initial_data():
    inp = 'data/unprocessed_data_initial_screening.csv'
    output = 'data/processed_data_initial_screening_1a-1e.csv'
    # read in the label key. It has three columns, Merck's label, our labs label, and the chemistry at the 4-position of the BOX ligand.
    label_key_data = pd.read_csv(label_key, header=0, index_col=0)
    print(label_key_data)

    # read in the unprocessed data.
    data = pd.read_csv(inp, header=0)
    print(data)
    # exit()

    # first make our ee column into strings
    data['ee'] = data['ee'].values.astype('str')
    # use some regex to throw out any characters that are not digits 1-9, decimal points, or minus signs
    data['ee'] = [re.sub('[^0-9.-]', '', i) for i in data['ee']]
    # replace our empty rows with nan  
    data['ee'].replace('', np.nan, inplace= True)
    # drop nan rows, now there are 566 rows
    data.dropna(axis=0, inplace=True)
    data['ee'] = data['ee'].values.astype('float')
    # do same for dr, but we won't drop the rows with nan
    data['dr'] = data['dr'].values.astype('str')
    data['dr'] = [re.sub('[^0-9.-]', '', i) for i in data['dr']]
    data['dr'].replace('', np.nan, inplace= True)
    data['dr'] = data['dr'].values.astype('float')
    print(data)
    # exit()

    # now we make the merck labels column into strings
    data['Merck label'] = data['Merck label'].values.astype('str')
    # remove those annoying eol characters - these are unnecessary anyway. also strip trailing whitespace
    data['Merck label'] = [i.split('\n')[0].strip() for i in data['Merck label']]
    print(data)

    # now assign all of the product labels correctly
    data['pdt'] = [product_key[i] for i in data['pdt']]
    print(data)

    # note the stereochemistry
    data['stereo'] = [label_key_data.loc[i, ['stereochemistry at 4 position']].values[0] for i in data['Merck label']]
    print(data)
    # assign the Denmark labels correctly
    data['Denmark label'] = [label_key_data.loc[i, ['Denmark label']].values[0] for i in data['Merck label']]
    print(data)
    # apply the correction to the S,S catalysts - get rows with S stereochem
    s_data = data.where(data['stereo'] == 'S')
    print(s_data)
    # invert the ee value
    s_data['ee'] = s_data['ee'] * -1
    print(s_data)
    # keep the rows that are R unchanged, else, replace with the value in s data
    print(data)
    data = data.where(data['stereo'] == 'R', s_data)
    print(data)

    # lets go ahead and drop the 2B data. These are mostly racemic, derived from non-methylated furan.
    print(data)
    data = data.where(data['pdt'] != '2B').dropna(axis = 0)
    print(data)

    # assign unique reaction handles for each one
    data['reaction handle'] = data['Denmark label'] + '_' + data['pdt']
    print(data)
        
    # now we average
    data = data.groupby('reaction handle').mean()
    pd.options.display.max_rows = 500

    data = data.sort_values(by = 'ee', ascending=False)
    # there are only 234 data points now... the one that was thrown out was 73_1_3_29_1A, due to overlapping peaks in the chromatogram
    print(data)
    # exit()

    # save the homogenized, corrected data
    with open(output, 'w') as o:
        data.to_csv(o)

def process_nn_and_fa():
    inp = 'data/unprocessed_data_fa_and_nn.csv'
    output = 'data/processed_data_fa_and_nn_1a-1e.csv'
    # read in the label key. It has three columns, Merck's label, our labs label, and the chemistry at the 4-position of the BOX ligand.
    label_key_data = pd.read_csv(label_key, header=0, index_col=0)
    print(label_key_data)

    # read in the unprocessed data.
    data = pd.read_csv(inp, header=0)
    print(data)
    # exit()

    # first make our ee column into strings
    data['ee'] = data['ee'].values.astype('str')
    # use some regex to throw out any characters that are not digits 1-9, decimal points, or minus signs
    data['ee'] = [re.sub('[^0-9.-]', '', i) for i in data['ee']]
    # replace our empty rows with nan  
    data['ee'].replace('', np.nan, inplace= True)
    # drop nan rows, now there are 566 rows
    data.dropna(axis=0, inplace=True)
    data['ee'] = data['ee'].values.astype('float')
    # do same for dr, but we won't drop the rows with nan
    data['dr'] = data['dr'].values.astype('str')
    data['dr'] = [re.sub('[^0-9.-]', '', i) for i in data['dr']]
    data['dr'].replace('', np.nan, inplace= True)
    data['dr'] = data['dr'].values.astype('float')
    print(data)
    # exit()

    # now we make the merck labels column into strings
    data['Merck label'] = data['Merck label'].values.astype('str')
    # remove those annoying eol characters - these are unnecessary anyway. also strip trailing whitespace
    data['Merck label'] = [i.split('\n')[0].strip() for i in data['Merck label']]
    print(data)

    # now assign all of the product labels correctly
    data['pdt'] = [product_key[i] for i in data['pdt']]
    print(data)

    # note the stereochemistry
    data['stereo'] = [label_key_data.loc[i, ['stereochemistry at 4 position']].values[0] for i in data['Merck label']]
    print(data)
    # assign the Denmark labels correctly
    data['Denmark label'] = [label_key_data.loc[i, ['Denmark label']].values[0] for i in data['Merck label']]
    print(data)
    # apply the correction to the S,S catalysts - get rows with S stereochem
    s_data = data.where(data['stereo'] == 'S')
    print(s_data)
    # invert the ee value
    s_data['ee'] = s_data['ee'] * -1
    print(s_data)
    # keep the rows that are R unchanged, else, replace with the value in s data
    print(data)
    data = data.where(data['stereo'] == 'R', s_data)
    print(data)

    # lets go ahead and drop the 2B data. These are mostly racemic, derived from non-methylated furan.
    print(data)
    data = data.where(data['pdt'] != '2B').dropna(axis = 0)
    print(data)

    # assign unique reaction handles for each one
    data['reaction handle'] = data['Denmark label'] + '_' + data['pdt']
    print(data)
        
    # now we average
    data = data.groupby('reaction handle').mean()
    pd.options.display.max_rows = 500

    data = data.sort_values(by = 'ee', ascending=False)
    # there are only 234 data points now... the one that was thrown out was 73_1_3_29_1A, due to overlapping peaks in the chromatogram
    print(data)
    # exit()

    # save the homogenized, corrected data
    with open(output, 'w') as o:
        data.to_csv(o)



if __name__ == '__main__':

    process_initial_data()
    process_nn_and_fa()
    pretty_process_initial_data()
    pretty_process_nn_and_fa()
