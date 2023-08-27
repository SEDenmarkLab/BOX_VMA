import molli as ml
from tqdm import tqdm
import os
import subprocess
from multiprocessing import Pool

# file of labelled cores, with bridging group in place (2d), attachment points labelled, A1: R4, A2: R5 syn, A3: R2 anti.
cf_cores_with_bridge = ml.ftypes.CDXMLFile('chemdraws/cores.cdxml')
# 4 position fragments
cf_4position = ml.ftypes.CDXMLFile('chemdraws/focused_library_r4_2.cdxml')

# this is where we save the files
out_mol2_dir = 'in_silico_fa_library/'
out_dir_uff_min = 'in_silico_fa_library_uff_min/'

def make_combinatorial_library():

    print('Making combinatorial library...\n')

    for item in tqdm(cf_cores_with_bridge.keys()):

        core = cf_cores_with_bridge[item]

        for subs in cf_4position.keys():

            pos4 = cf_4position[subs]

            # print(pos4.dumps_mol2())
            # print(core.dumps_mol2())
            # exit()

            # make core structure - A1 are the attachment points for 4 position
            m1 = ml.Molecule.join(core, pos4, 'A1', 'AP0', optimize_rotation=True)
            core_with_bridge = ml.Molecule.join(m1, pos4, 'A1', 'AP0', optimize_rotation=True)
            core_with_bridge.add_implicit_hydrogens()
            with open(f'{out_mol2_dir}{pos4.name}_{core.name}.mol2', 'w') as o:
                core_with_bridge.dump_mol2(o)
    return 1


# we will use pool to make this faster
def minimize_obabel(pool_size = 64):

    print('Minimizing library...')
    args = os.listdir(out_mol2_dir)
        
    with Pool(pool_size) as p:
        p.map(pool_handler_minimize_obabel, args)
    
    return 1

# the pool handler uses obabel command line for speed
def pool_handler_minimize_obabel(arg):
    if arg in os.listdir(out_dir_uff_min):
        return
    else:
        print(f'processing file {arg}...')
        # usinf uff force field
        with open(out_dir_uff_min + arg, 'w') as o:
            subprocess.run(['obminimize', '-ff', 'uff', '-o', 'mol2', f'{out_mol2_dir}{arg}'], stdout= o) 
        return 
    


if __name__ == '__main__':

    # make directories if we need them
    if not os.path.exists(out_mol2_dir):
        os.makedirs(out_mol2_dir)
    if not os.path.exists(out_dir_uff_min):
        os.makedirs(out_dir_uff_min)


    if make_combinatorial_library() == 1:
        print("Successfully made FA combinatorial library!")

    if minimize_obabel(pool_size=128) == 1:
        print("Successfully minimized combinatorial library!")



    print('Successfully made and minimized in silico library!')
            