
# from rdkit import Chem
# from rdkit.Chem import AllChem
# from rdkit.Chem.AllChem import AlignMol, EmbedMolecule, EmbedMultipleConfs
# from rdkit.Chem.rdForceFieldHelpers import UFFGetMoleculeForceField
# from rdkit.Chem import Draw,PyMol
# from rdkit.Chem.Draw import IPythonConsole
import random
import os
from tqdm import tqdm
from multiprocessing import Pool
import shutil
import pandas as pd


# starting ligands before conformer generation
in_dir = 'in_silico_fa_library_uff_min_checked/'

# where the separated conformers are written out
out_dir1 = 'out_fa_confs/out_conformers_fa_separated/'
# mxyz files
out_dir_mxyz1 = 'out_fa_confs/out_conformers_fa_mxyz/'
# log file from conformer generation
log_file = 'confgen_fa.out'

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







