import os
from rdkit import Chem
import shutil
from multiprocessing import Pool

# change these lines to your input/output directories
input_dir = 'out1/out_conformers1_separated/'
checked_output = 'out1/out_conformers1_separated_checked/'
failed_labelling = 'out1/out_conformers1_separated_failed/'
pool_size = 64

# start 1973453

def pool_handler(inp):
    # read in the file. You have to do it like this cause RDKit is awesome. https://github.com/rdkit/rdkit/issues/5357
    try:
        mol = Chem.rdmolfiles.MolFromMolFile(input_dir + inp, sanitize = False)
        Chem.SanitizeMol(mol)
    except Exception:
        print(f'{inp} failed file reading!')
        return 0
    try:
        Chem.rdCIPLabeler.AssignCIPLabels(mol)
    except Exception:
        # s = [a.GetProp('_CIPCode') for a in chiral_centers]
        print(f'{inp} failed cip labelling!')
        return 0

    # get the chiral center atoms
    chiral_centers = [atom for atom in mol.GetAtoms() if atom.HasProp('_CIPCode')]
    s = [a.GetProp('_CIPCode') for a in chiral_centers]
    
    # break it up into 4,5 and any other stereocenters
    pos4 = []
    pos5 = []
    others = []
    for atom in chiral_centers:
        if 'N' in [a.GetSymbol() for a in atom.GetNeighbors()]:
            pos4.append(atom.GetProp('_CIPCode'))
        elif 'O' in [a.GetSymbol() for a in atom.GetNeighbors()]:
            pos5.append(atom.GetProp('_CIPCode'))
        else:
            others.append(atom.GetProp('_CIPCode'))
    
    # # Make sure we have partitioned sets correctly
    try:
        # assert len(others) % 2 == 0 # this will fail for the admantyl ones
        assert len(pos5) % 2 == 0
        assert len(pos4) % 2 == 0
    except AssertionError:
        print(f'{inp} failed chiral center partition! {s}')
        return 0

    ## groups with sulfur , 33 = 2-thiophene, 34 = 2-furanyl, 55 = sulfide, 273 = ester, are R config by CIP rules
    special_R4s = ['33', '34', '55', '273']
    pos4_label = inp.split('.')[0].split('_')[0]

    # all the catalysts should be (R) configured at both 4,4' positions, if there's any S we toss it.
    if ('S' in pos4) and (pos4_label not in special_R4s):
        return -1
    # make sure 4pos are the same configuration
    elif len(set(pos4)) != 1:
        return -1
    # 5-position should have same CIP configuration if it is C2 symmetric
    elif (len(pos5) !=0) and (len(set(pos5)) != 1):
        return -1
    elif (pos4_label == 'aa') and ('R' in pos5):
        return -1
    # if it passes all of this it is chiral, C2 symmetric
    else: 
        shutil.move(input_dir + inp, checked_output + inp)
        return 1


if __name__ == '__main__':
    args = os.listdir(input_dir)
    print(len(args))

    if not os.path.exists(checked_output):
        os.makedirs(checked_output)
    if not os.path.exists(failed_labelling):
        os.makedirs(failed_labelling)

    with Pool(pool_size) as p:
        rtn = p.map(pool_handler, args)
        failed = [args[i] for i, rtncode in enumerate(rtn) if rtncode == 0]
        print(f'Failed for stupid reason:\n{failed}\n')
        for i in failed:
            shutil.move(input_dir + i, failed_labelling + i)
    
    print('Success!')