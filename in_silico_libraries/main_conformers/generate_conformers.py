
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from multiprocessing import Pool

in_dir = '../main_library/in_silico_library_uff_min_checked1/'

out_dir = 'out1/out_conformers1_separated/'
out_dir_mxyz = 'out1/out_conformers1_mxyz/'
reference_structure_file = 'alignment_core.mol'


def get_oxazaline_alignment_atoms(mol):
    """ Function that takes in a free bisoxazoline ligand and identifies and returns the atom indices of the nitrogen in one ring, C2 of that ring,
    the bridging methylene, the C2 of the other ring, and the nitrogen of the other ring, in that order. 

    argument mol: RDKit molecule

    returns: list of atom indices of the nitrogen in one ring, C2 of that ring,
    the bridging methylene, the C2 of the other ring, and the nitrogen of the other ring, in that order. 

    """
    # find all the nitrogens, we will use different cases of these to identify the structure
    nitrogens = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
    # print(len(nitrogens))
    first_N, second_N, c2, bridge, other_c2 = None, None, None, None, None
    oxazoline_Ns = []

    # the standard case, when the only two nitrogens in a molecule are the oxazolines
    if len(nitrogens) == 2:
        oxazoline_Ns = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'N']
        
    # the case when the bridging group is nitrile
    elif len(nitrogens) == 3:
        # print('here')
        # the nitrile N only has one neighbor
        for atom in nitrogens:
            if len(mol.GetAtomWithIdx(atom).GetNeighbors()) != 1:
                oxazoline_Ns.append(atom)
        

    # the case where we have nitrogens in other parts of the molecule, e.g. pyridine rings
    elif len(nitrogens) == 4:
        oxazoline_Ns = []
        # iterate through nitrogens
        for atom in nitrogens:
            # print('here1')
            neighbors = mol.GetAtomWithIdx(atom).GetNeighbors()
            # iterate through nitrogen neighbors
            for neighbor in neighbors:
                # oxazoline C2 should be SP2
                if str(neighbor.GetHybridization()) == 'SP2':
                    # oxazoline C2 should have exactly 3 neighbors, C bridge, N, and O
                    if sorted([i.GetSymbol() for i in neighbor.GetNeighbors()]) == sorted(['C', 'N', 'O']):
                        # print('here2')
                        oxazoline_Ns.append(atom)
                        break
    
    # the cases where we have oxazoline nitrogens, nitrogens at 4-position, and a nitrile bridging group, and the case
    # with oxazoline nitrogens, nitrogens at 4-position, and bridging group 29 with 2 nitrogens can be treated similarly (5 and 6 total nitrogens)
    elif (len(nitrogens) == 5) or (len(nitrogens) == 6):
        oxazoline_Ns = []
        # iterate through nitrogens
        for atom in nitrogens:
            neighbors = mol.GetAtomWithIdx(atom).GetNeighbors()
            # this is the nitrile, skip it
            if len(neighbors) == 1:
                continue
            # iterate through nitrogen neighbors
            for neighbor in neighbors:
                # oxazoline C2 should be SP2
                if str(neighbor.GetHybridization()) == 'SP2':
                    # oxazoline C2 should have exactly 3 neighbors, C bridge, N, and O
                    if sorted([i.GetSymbol() for i in neighbor.GetNeighbors()]) == sorted(['C', 'N', 'O']):
                        # print('here2')
                        oxazoline_Ns.append(atom)
                        break

    else:
        raise Exception( f'{mol.name} failed! nitrogen count: {len(nitrogens)}' )


    assert len(oxazoline_Ns) == 2
    first_N, second_N = oxazoline_Ns
    # get c2 next to first N
    
    # this logic still works for the cbr2 linkers (which are CSp3)
    c2 = [atom.GetIdx() for atom in mol.GetAtomWithIdx(first_N).GetNeighbors() if str(atom.GetHybridization()) == 'SP2'][0]
    other_c2 = [atom.GetIdx() for atom in mol.GetAtomWithIdx(second_N).GetNeighbors() if str(atom.GetHybridization()) == 'SP2'][0]
    # print([i.GetSymbol() for i in mol.GetAtomWithIdx(c2).GetNeighbors()])
    # get bridge
    bridge = set([atom.GetIdx() for atom in mol.GetAtomWithIdx(c2).GetNeighbors()]).intersection(set([atom.GetIdx() for atom in mol.GetAtomWithIdx(other_c2).GetNeighbors()]))
    assert len(bridge) == 1
    bridge = bridge.pop()

    assert second_N in [i.GetIdx() for i in mol.GetAtomWithIdx(other_c2).GetNeighbors()]
    
    return first_N, c2, bridge, other_c2, second_N

""" 
The pool handler does the heavy lifting
"""
## see https://www.rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
def pool_handler(args):
    # separate our arguments
    reference = args[0] # reference mol object
    arg = args[1] # file for reading
    ref_alignment_atoms = args[2] # alignment atoms from ref
    threshold = 15.0 # energy cutoff, change if you want different

    try: # try to read in the file, calc num rotatable bonds
        # read in the file, set chiral tags
        mol = AllChem.MolFromMol2File(in_dir + arg, sanitize = False)
        mol.name = arg.split('.')[0]
        Chem.SanitizeMol(mol)
        AllChem.AssignCIPLabels(mol)
        # get the number of rotatable bonds, this will determine how many confs to make
        num_rotatable_bonds = AllChem.CalcNumRotatableBonds(mol)
        n_confs = 0
        if num_rotatable_bonds <= 7:
            n_confs = 50
        if num_rotatable_bonds >=8:
            n_confs = 150
    except Exception as exp:
        print(f'Failure for {arg} on file reading/rotatable bond calculation: {exp!s}')
        return

    try: # try to embed conformers
        # set some embedding parameters
        ps = AllChem.ETKDGv2()
        ps.maxIterations = 10000
        ps.randomSeed = 42
        ps.enforceChirality = True
        ps.useRandomCoords = True
        # add Hs
        mol = Chem.AddHs(mol, addCoords = True)
        # generate n_conformers
        conf_ids = Chem.rdDistGeom.EmbedMultipleConfs(mol, numConfs = n_confs, params = ps)
    except Exception as exp:
        print(f'Failure for {arg} on conformer embedding: {exp!s}')
        return

    try: # try to optimize conformer distribution
        # optimize with MMFF
        AllChem.MMFFSanitizeMolecule(mol)
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads = 1, maxIters = 200)
    except Exception as exp:
        print(f'Failure for {arg} on conformer geometry optimization: {exp!s}')
        return

    try: # try to align conformers with reference
        # align the conformers
        mol_alignment_atoms = list(get_oxazaline_alignment_atoms(mol))
        mol_map = [i for i in zip(mol_alignment_atoms, ref_alignment_atoms)]
        AllChem.AlignMol(mol, reference, atomMap = mol_map) # this line aligns the first conf to the reference mol
        AllChem.AlignMolConformers(mol ,atomIds = mol_alignment_atoms) # this line aligns mol confs to first conf
    except Exception as exp:
        print(f'Failure for {arg} on alignment: {exp!s}')
        return


    try: # try to calculate energies for conformers
        # sort the conformer list by increasing energy
        energy_list = []
        mol_props = AllChem.MMFFGetMoleculeProperties(mol)
        for cid in list(conf_ids):
            # print(f'minimming conf {cid}')
            ff = AllChem.MMFFGetMoleculeForceField(mol, mol_props, confId = cid) 
            energy_list.append(ff.CalcEnergy())
    except Exception as exp:
        print(f'Failure for {arg} on energy calculation: {exp!s}')
        return

    # sorted conformer energy list
    sorted_conformers = sorted(zip(list(conf_ids), energy_list), key = lambda x: x[1], reverse = False)
    min_energy = min([i[1] for i in sorted_conformers])
    # filter conformers based on energy threshold
    filtered_conformers = [i for i in sorted_conformers if (i[1] - min_energy < threshold)]

    try:
        # write out multi xyz files
        with open(out_dir_mxyz + arg.split('.')[0] + '.xyz', 'w') as out:
            for cid, energy in filtered_conformers:
                out.write(Chem.MolToXYZBlock(mol, confId = cid))
        # write out separated mol files
        for i, tuple in enumerate(filtered_conformers):
            Chem.MolToMolFile(mol, out_dir + arg.split('.')[0] + '_' + str(i) + '.mol',  confId = tuple[0])
    except Exception as exp:
        print(f'Failure for {arg} on file writing: ' + exp.ToString())
    else:
        print(f'Success for {arg} with {num_rotatable_bonds} rotatable bonds! {len(filtered_conformers)} conformers succesfully generated!')
    return

    

if __name__ == '__main__':

    # make directories
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if not os.path.exists(out_dir_mxyz):
        os.makedirs(out_dir_mxyz)

    reference_mol = AllChem.MolFromMolFile(reference_structure_file)
    files = os.listdir(in_dir)
    # it will save us some time to call this only once here
    alignment_atoms = list(get_oxazaline_alignment_atoms(reference_mol))

    args = [i for i in zip([reference_mol for i in range(len(files))], files, [alignment_atoms for i in range(len(files))])]

    
    with Pool(100) as p:
        p.map(pool_handler, args)
    print('Success!')
    # job id: 1082045


