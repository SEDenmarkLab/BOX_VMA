
import os
from multiprocessing import Pool
import pandas as pd

# starting ligands before conformer generation
in_dir = '../main_library/in_silico_library_uff_min_checked1/'

# where the separated conformers are written out
out_dir1 = 'out1/out_conformers1_separated/'
# mxyz files
out_dir_mxyz1 = 'out1/out_conformers1_mxyz/'
# log file from conformer generation
log_file = 'confgen1.out'

# make a nice dataframe from log file
log_df = []
with open(log_file, 'r') as o:
  for line in o:
    spl = line.split(' ')
    # if spl[0] == 'processing':
    #   continue
    if spl[0] == 'Success':
      dfrow = [spl[2], 'pass', spl[4], spl[7]]
    if spl[0] == 'Failure':
      dfrow = [spl[2], 'fail', spl[4], '']
    log_df.append(dfrow)
# label columns
log_df = pd.DataFrame(log_df, columns = ['file', 'pass?', 'num rotatable bonds', 'num confs'])
# save log file cleaned up
with open('out1_log.csv', 'w') as o:
  log_df.to_csv(o)
print(log_df)

failures = log_df.loc[log_df['pass?'] == 'fail']
print(failures)

successes = log_df.loc[log_df['pass?'] == 'pass']
print(successes)

print(len(os.listdir(out_dir_mxyz1)))
outset = set([i.split('.')[0] for i in os.listdir(out_dir_mxyz1)])
# print(outset)

good_files = set([i.split('.')[0] for i in successes.loc[:, 'file']])


print(len(outset.difference(good_files)))
print(outset.difference(good_files))







