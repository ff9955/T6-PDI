from src.pyAOM_utils import *


# STO parameters for GTO-STO projection step, refers to exponential coefficients of STO basis functions
STO_proj_dict={
    'H':{'STOs':1  ,'1s':1.0000},
    'C':{'STOs':1+3,'2s':1.6083,'2p':1.442657},
    'N':{'STOs':1+3,'2s':1.9237,'2p':1.646703},
    'O':{'STOs':1+3,'2s':2.2458,'2p':1.858823},
    'F':{'STOs':1+3,'2s':2.5638,'2p':2.136394},
    'S':{'STOs':1+3,'3s':2.1223,'3p':1.651749},
}

# general configuration
config={
    'label': '6T',
#    'single_mol_xyz': 'single_molecules/6T_opt_geom_EnGn.xyz',
    'single_mol_xyz': 'single_molecules/6T_opt_geom_EnGn.xyz',
    'basis': 'DZVP-GTH',
    'potential': 'GTH-PBE',
    'MO': 73,
    'cp2k_template': 'templates/sp_GGA_template.inp',
    'cp2k_basis_file': 'cp2k_files/GTH_BASIS_SETS',
}
# CP2K DFT parameters
CP2K_QS_GGA_parameters={
    'PROJECT_NAME': '6T',
    'XC_FUNCTIONAL': 'PBE',
    'BASIS_SET_FILE_NAME': 'GTH_BASIS_SETS',
    'POTENTIAL_FILE_NAME': 'POTENTIAL',
    'MGRID_CUTOFF': 450.0,
    'MGRID_REL_CUTOFF': 75.0,
    'SCF_MAX_SCF': 200,
    'MO_CUBES_NHOMO': 6,
    'MO_CUBES_NLUMO': 6,
    'MO_CUBES_STRIDE': 5,
}
# default cp2k:qs output path
cp2k_out=f'cp2k_output/{CP2K_QS_GGA_parameters["PROJECT_NAME"]}.log'

# read a single molecule xyz file
mymol=single_molecule(config['single_mol_xyz'])
# prepare CP2K:QS input file for single molecule
mymol.prep_cp2k_single(CP2K_QS_GGA_parameters,config['label'],config['cp2k_template'],config['basis'],config['potential'])
print(f'Note: make sure CP2K:QS output gets in "{cp2k_out}"')

# provided we carried out the DFT calculation and placed the CP2K output in the appropriate directory,
# read CP2K:QS output and retrieve MO info
mymol.get_cp2k_info(config['MO'],cp2k_out,config['cp2k_basis_file'],config['basis'])

# STO init and overlap matrix calculation
mymol.initialize_STOs(STO_proj_dict)
# carry out GTO->STO projection
mymol.project()

# store AOM coefficients
mymol.save_AOM(config['label'])

# projection completeness
mymol.orb_compl_dict

dict_anumber = {
    "C": 6,
    "H": 6,
    "S": 16,
    "N": 7,
    "O": 8,
}

names = mymol.species
coeff = mymol.AOM_pi_coeffs

for index,n in enumerate(mymol.species):
    print(names[index], dict_anumber[names[index]], "1", "0.00000", coeff[index])

# define STO parameters for AOM calculations
AOM_dict={
    'H':{'STOs':1  ,'1s':1.0000},
    'C':{'STOs':1+3,'2s':1.6083,'2p':1.385600},
    'N':{'STOs':1+3,'2s':1.9237,'2p':1.617102},
    'O':{'STOs':1+3,'2s':2.2458,'2p':1.505135},
    'F':{'STOs':1+3,'2s':2.5638,'2p':1.665190},
    'S':{'STOs':1+3,'3s':2.1223,'3p':1.641119},
}

# locate dimers: define the directory with dimers
# suggested file structure: create a directory with the molecule name in 'dimers'
# and use the following naming for each dimer: <name>_<No>_<comment>.xyz
dimer_dir='dimers/OT6_crystal_pairs/'
# define the AOM files for both dimer fragments
frag1_AOM_file='output/6T/AOM_COEFF.dat'
frag2_AOM_file='output/6T/AOM_COEFF.dat'
# define MOs
frag1_MO=73
frag2_MO=73
#--------------------------------------------------------------------------------



#Sab_array=[]
#for file in sorted([i for i in os.listdir(dimer_dir) if i.endswith('.xyz')==True]):
#    print(file)
#    dimer_xyz_file=f'{dimer_dir}/{file}'
#    val=Sab(dimer_xyz_file,frag1_AOM_file,frag2_AOM_file,frag1_MO,frag2_MO,AOM_dict)
#    Sab_array.append(val)
#    print(f'{file} Sab: {val}')

print(os.listdir(dimer_dir))
name=["/P_pair.xyz",  "/Px_pair.xyz",  "/T1_pair.xyz"]
#name=["/P.xyz",  "/T1_close.xyz",  "/T2_close.xyz"]

Sab_array=[]
for file in name:
    #print(file)
    dimer_xyz_file = dimer_dir + file
    print(dimer_xyz_file)

    val=Sab(dimer_xyz_file,frag1_AOM_file,frag2_AOM_file,frag1_MO,frag2_MO,AOM_dict)
    Sab_array.append(val)
    print(f'{file} Sab: {val}')

# define STO parameters for AOM calculations
AOM_dict={
    'H':{'STOs':1  ,'1s':1.0000},
    'C':{'STOs':1+3,'2s':1.6083,'2p':1.385600},
    'N':{'STOs':1+3,'2s':1.9237,'2p':1.617102},
    'O':{'STOs':1+3,'2s':2.2458,'2p':1.505135},
    'F':{'STOs':1+3,'2s':2.5638,'2p':1.665190},
    'S':{'STOs':1+3,'3s':2.1223,'3p':1.641119},
}

# locate dimers: define the directory with dimers
# suggested file structure: create a directory with the molecule name in 'dimers'
# and use the following naming for each dimer: <name>_<No>_<comment>.xyz
dimer_dir='dimers/'
# define the AOM files for both dimer fragments
frag1_AOM_file='output/6T/AOM_COEFF.dat'
frag2_AOM_file='output/PDI/HOMO_AOM_COEFF.dat'
# define MOs
frag1_MO=73
frag2_MO=70
#--------------------------------------------------------------------------------



#Sab_array=[]
#for file in sorted([i for i in os.listdir(dimer_dir) if i.endswith('.xyz')==True]):
#    print(file)
#    dimer_xyz_file=f'{dimer_dir}/{file}'
#    val=Sab(dimer_xyz_file,frag1_AOM_file,frag2_AOM_file,frag1_MO,frag2_MO,AOM_dict)
#    Sab_array.append(val)
#    print(f'{file} Sab: {val}')

print(os.listdir(dimer_dir))
name=["/interface_dimer.xyz", "interface_problem_stuc283-305.xyz", "interface_problem_stuc283-305-step1.xyz", "interface_problem_stuc283-329.xyz", "interface_problem_stuc283-305-allameno2.xyz" ]
#name=["/P.xyz",  "/T1_close.xyz",  "/T2_close.xyz",]

Sab_array=[]
for file in name:
    #print(file)
    dimer_xyz_file = dimer_dir + file
    print(dimer_xyz_file)

    val=Sab(dimer_xyz_file,frag1_AOM_file,frag2_AOM_file,frag1_MO,frag2_MO,AOM_dict)
    Sab_array.append(val)
    print(f'{file} Sab: {val}')
