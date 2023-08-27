import molli as ml
from openbabel import openbabel as ob
import os
from tqdm import tqdm
import numpy as np

indir1 = 'out_fa_confs/out_conformers_fa_separated_checked/'

indir2 = '../main_conformers/out1/out_conformers1_separated_checked/'

# from Blake
def load_obmol(fname, input_ext: str = 'xyz') -> ob.OBMol:
    '''
    This function takes any file and creates an openbabel style mol format
    '''
    conv = ob.OBConversion()
    obmol = ob.OBMol()
    conv.SetInFormat(input_ext)
    conv.ReadFile(obmol, fname)

    obmol.ConnectTheDots()
    obmol.PerceiveBondOrders()

    return obmol

# from Blake
def obmol_to_mol2(obmol: ob.OBMol):
    '''
    This returns a basic mol2 block
    '''

    conv = ob.OBConversion()
    conv.SetOutFormat('mol2')
    return conv.WriteString(obmol)

def make_fa_collection():
    
        print('made fa only library!')
        return

def append_library():
    return None

if __name__ == '__main__':

    
    

    lib = ml.ConformerLibrary.new("out_fa_confs/conformers_no_linker_fa_only.mlib")

    with lib:
        conformer_dict = {}
                ## first pair up the conformer files with the core structure in a dictionary.
        for f in os.listdir(indir1):
            core = '_'.join(f.split('_')[:-1])
            if core in conformer_dict:
                conformer_dict[core].append(f)
            else:
                conformer_dict[core] = [f]

        ## sort the conformers so that they are in order of increasing energy
        for core in conformer_dict:
            conformer_dict[core].sort(key = lambda x: int(x.split('.')[0].split('_')[-1]))

        print(f'processing {len(conformer_dict)} core files!')
        # iterate through cores
        for core in tqdm(conformer_dict):
            print(f'processing {core}...')
            # list of conformer files
            conformer_files = conformer_dict[core]
            # load mol file with openbabel
            obmol2s = [load_obmol(indir1 + f, input_ext = 'mol') for f in conformer_files]
            # convert mol into mol2 block
            mol2s = [obmol_to_mol2(mol) for mol in obmol2s]
            # read mol2 block into molli Molecule with the core name assigned
            mollimols = [ml.Molecule.loads_mol2(mol, name = core) for mol in mol2s]

            # delete the two bromines and the carbon linker in all of the structures
            for mol in mollimols:

                # get the atoms we need to remove
                bromines = [i for i in mol.yield_atoms_by_element('Br')]
                assert len(bromines) == 2
                carbon = next(mol.connected_atoms(bromines[0]))

                # remove the atoms. I OVERRODE DEL_ATOM IN MOLECULE CLASS TO DEAL WITH ATOMIC_CHARGES
                for a in bromines:
                    mol.del_atom(a)
                mol.del_atom(carbon) 

            # make the ensemble
            ensemble = ml.ConformerEnsemble(
                mollimols[0],
                n_conformers=len(mollimols),
                n_atoms=mollimols[0].n_atoms,
                name=mollimols[0].name,
                atomic_charges=mollimols[0].atomic_charges
            )

            

            # add the coordinates of each conformer to the ensemble
            for i, m in enumerate(mollimols):
                ensemble._coords[i] = m.coords

            # ensembles.append(ensemble)
            # add to FA only ensemble
            lib.append(ensemble.name, ensemble)
            # add this ensemble to the larger library
            # r4pos = ensemble.name.split('_')[0]

    all_fa_cats = [i.name for i in lib]
    assert len(set(all_fa_cats)) == len(all_fa_cats)
    print(f'# fas: {len(set(all_fa_cats))}')

    lib2 = ml.ConformerLibrary(path = "out_fa_confs/conformers_no_linker_main_with_fa.mlib", readonly=False)
    # lib2 = ml.ConformerLibrary(path = "../main_conformers/out1/conformers_no_linker.mlib", readonly=True)
    full_lib_cats = [i.name for i in lib2]
    assert len(set(full_lib_cats)) == len(full_lib_cats)
    print(f'# main: {len(set(full_lib_cats))}')

    overlap = set(full_lib_cats).intersection(set(all_fa_cats))
    print(f'# redundant = {len(overlap)}')
    print(overlap)

    to_add = set(all_fa_cats).difference(set(full_lib_cats))
    print(f'# to add = {len(to_add)}')

    # exit()
    with lib2:
        for e in lib:
            if e.name in to_add:
                lib2.append(e.name, e)

    full_lib_cats = [i.name for i in lib2]
    assert len(set(full_lib_cats)) == len(full_lib_cats)
    print(f'# main after adding: {len(set(full_lib_cats))}')



    print('Success!')


