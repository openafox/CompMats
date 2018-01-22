"""Dimer: Diffusion along rows"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import pdb
import subprocess
import shlex
from math import sqrt, pi

# from ase import Atoms, Atom
from ase.io import Trajectory
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.calculators.emt import EMT
from ase.dimer import DimerControl, MinModeAtoms, MinModeTranslate
from ase.lattice.surface import fcc111
from ase.lattice.surface import fcc100
from ase.lattice.surface import add_adsorbate

from ase import *
from ase.visualize import view
from ase.io import *
from ase.io.trajectory import *
import sys
import time
import os
from ase.md import VelocityVerlet
from ase.calculators.lammpsrun import *  # LAMMPS
from ase.constraints import FixAtoms
from ase.lattice.surface import surface
from ase.lattice.cubic import FaceCenteredCubic
# ################remember to set up ASE LAMMPS############
# add to .profile or .bashrc or .bash_profile file:
# export LAMMPS_COMMAND=/bin/lmp
# close Terminal
# and reopen
#######################

def histplot():

    data = np.loadtxt("./data/reps1/reps.dat",delimiter="\t",skiprows=1)
    i = data[:,0]
    PE = data[:,1]
    dist = data[:,2]
    plt.hist(PE, 50, normed=1, facecolor='green', alpha=0.75)
    plt.xlabel('Energy [eV]')
    plt.ylabel('Probability [number/600]')
    #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
    #plt.ylim([-1,10])
    #  plt.axis([40, 160, 0, 0.03])
    plt.grid(True)
    plt.show()


############ dimer ##############################
path = os.path.dirname(os.path.abspath(__file__))
print(path)
starttime = time.strftime("%y%m%d%H%M%S")
a=3.9242
# Setting up the initial image:
# initial = fcc111('Pt', size=(4, 4, 2), a=3.9242, vacuum=10.0)
# add_adsorbate(initial, 'Pt', 1.611, 'hollow')
"""Pt_atoms = FaceCenteredCubic(directions=[[0,1,0], [0,0,1,], [1,0,0]],
                            size=(4,4,1),latticeconstant=a,
                            symbol='Pt', pbc=(1,1,1))
initial = surface(Pt_atoms, (0,0,1), 3)
# [1,-1,0], [1,1,-2], [1,1,1]], # for 111 surface"""
initial = fcc100('Pt', size=(4, 4, 5), a=a)
add_adsorbate(initial, 'Pt', 1.611, 'hollow')  #add adatom on surface
at = initial[len(initial)-6]
at_pos = at.position
at_ind = at.index
print(at)
print(at_pos)
print(at_ind)
# del initial[len(initial)-7]  # delete first atom ie create vacancy
initial.center(vacuum=10.0, axis=2)
# Save it
# write('Molecule-%s-%s-start.xyz' % (str(11), str(1)), initial)
N = len(initial)  # number of atoms

# Freeze the "slab"
# Make a mask of zeros and ones that select fixed atoms - the two
# bottom layers:
# mask = initial.positions[:, 2] - min(initial.positions[:, 2]) < 1 * a 
# mask = [atom.index < 4*4 for atom in initial]
mask = [atom.tag > 4 for atom in initial]
initial.set_constraint(FixAtoms(mask=mask))
# initial.set_pbc((True,True,True))
print(initial.get_pbc)
# view(initial)

# Set up LAMMPS
# calc = LAMMPS('initial', tmp_dir='./', keep_tmp_files=True,
# no_data_file=False) #files = ['../CH.airebo'] \n min_style fire \n  unfix
# fix_nve','run': '1', , 'minimize': '10e-12 10e-14 100000 100000'
calc = LAMMPS('initial',
              # tmp_dir='./data/lmp/',
              # keep_tmp_files=False,
              parameters={'pair_style': 'eam',
                          'pair_coeff': ['* * ' + path + '/data/Pt_u3.eam']})
# Calculate using LAMMPS:

initial.set_calculator(calc)

# Relax the initial state:
QuasiNewton(initial).run(fmax=0.05)  # I don't get why but with this
# minimiztion it does not work correctly
e0 = initial.get_potential_energy()
print(e0, type(e0))

# ##### Run multiple repitionions of the dimer while creating energy bins ####
  
adatom = initial[len(initial)-1]
adatom_pos_in = adatom.position
print("adatom pos=", adatom_pos_in)
reps = 10000
rep = 0

fpath = './data/reps1/'
with open(fpath + "reps.dat","w") as f:
    f.write("Rep\tPE Diff\tAdatom displacment\n")
while True:
    runner = initial.copy()
    runner.set_calculator(calc)   # for some reason this needs reset
    fname = '%s%03.0f_dimer' %(fpath,rep)
    trajname = fname + '.traj'  
    print(trajname)

    traj = Trajectory(trajname, 'w', runner)
    traj.write()

    # Making dimer mask list:
    # d_mask = [atom.index == 69 for atom in initial]
    d_mask = [False] * (N - 1) + [True]  # [0, 0, 0, 0, 1]
    # d_mask = initial.positions[:, 2] - min(initial.positions[:, 2]) < 1 * a 
    # Set up the dimer (gaus is for ramdom rotation std sets up the level):
    d_control = DimerControl(
            initial_eigenmode_method='gauss',
            # initial_eigenmode_method='displacement',
            displacement_method='gauss',  # 'vector',
            # displacement_method='vector',
            gauss_std=0.1,
            trial_angle=pi / 12.0,
            logfile=fname + '.log',
            mask=d_mask)
    """ All the avaiable parameters are:
    parameters = {'eigenmode_method': 'dimer',
                  'f_rot_min': 0.1,
                  'f_rot_max': 1.00,
                  'max_num_rot': 1,
                  'trial_angle': pi / 4.0,
                  'trial_trans_step': 0.001,
                  'maximum_translation': 0.1,
                  'cg_translation': True,
                  'use_central_forces': True,
                  'dimer_separation': 0.0001,
                  'initial_eigenmode_method': 'gauss''displacement',
                  'extrapolate_forces': False,
                  'displacement_method': 'gauss''vector,
                  'gauss_std': 0.1,
                  'order': 1,
                  'mask': None, # NB mask should not be a "parameter"
                  'displacement_center': None,
                  'displacement_radius': None,
                  'number_of_displacement_atoms': None}"""
    d_atoms = MinModeAtoms(runner, d_control)
    """# Enable with 'displacement' - 'vector'
    # Displacement settings:
    # initial vector:
    displacement_vector = np.zeros((N, 3))
    # Strength of displacement along y axis = along row:
    displacement_vector[-1, 0] = 0.1
    # The direction of the displacement needed if usint displacement vector method 
    # displacement_vector[-1, a], where a can be 0 for x, 1 for y and 2 for z.
    d_atoms.displace(displacement_vector=displacement_vector)"""

    # Converge to a saddle point:
    dim_rlx = MinModeTranslate(d_atoms,
                               trajectory=traj)  # , logfile=None)
    dim_rlx.run(fmax=0.001)

    diff = runner.get_potential_energy() - e0
    """print(('The energy barrier is %f eV.' % diff))
    args = shlex.split('ase-gui ' + trajname)
    print(args)
    subprocess.call(args)"""

    adatom = runner[len(initial)-1]
    adatom_pos_f = adatom.position
    adatom_dis = np.sqrt(np.square(adatom_pos_in - adatom_pos_f).sum())
    with open(fpath + "reps.dat","a+") as f:
        f.write(str(rep) + "\t" + str(diff) + "\t" + str(adatom_dis) + "\n")
    rep += 1
    if rep == reps:
        break
histplot()
