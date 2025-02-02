import molli as ml
from openbabel import openbabel as ob
import os
from tqdm import tqdm

indir = 'out1/out_conformers1_separated_checked/'

# from Blake Ocampo
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

# from Blake Ocampo
def obmol_to_mol2(obmol: ob.OBMol):
    '''
    This returns a basic mol2 block
    '''

    conv = ob.OBConversion()
    conv.SetOutFormat('mol2')
    return conv.WriteString(obmol)


if __name__ == '__main__':

    conformer_dict = {}

    # this will be our conformer library
    lib = ml.ConformerLibrary.new("out1/conformers_no_linker.mlib",)

    ## first pair up the conformer files with the core structure in a dictionary.
    for f in os.listdir(indir):
        core = '_'.join(f.split('_')[:-1])
        if core in conformer_dict:
            conformer_dict[core].append(f)
        else:
            conformer_dict[core] = [f]

    ## sort the conformers so that they are in order of increasing energy
    for core in conformer_dict:
        conformer_dict[core].sort(key = lambda x: int(x.split('.')[0].split('_')[-1]))

    print(f'processing {len(conformer_dict)} core files!')
    
    with lib:
        # iterate through cores
        for core in tqdm(conformer_dict):
            print(f'processing {core}...')
            # list of conformer files
            conformer_files = conformer_dict[core]
            # load mol file with openbabel
            obmol2s = [load_obmol(indir + f, input_ext = 'mol') for f in conformer_files]
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
                atomic_charges=mollimols[0].atomic_charges,
            )

            # add the coordinates of each conformer to the ensemble
            for i, m in enumerate(mollimols):
                ensemble._coords[i] = m.coords

            # add this ensemble to the larger library
            lib.append(ensemble.name, ensemble)

    print('Success!')


