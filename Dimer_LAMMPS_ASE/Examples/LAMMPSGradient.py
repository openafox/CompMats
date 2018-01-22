#!/usr/bin/python
#IMPORT PYTHON PACKAGES
#import paramiko
import numpy as np
import math 
from ase import *
from ase.visualize import view
from ase.io import *
from ase.io.trajectory import *
import sys
import time
from ase import Atoms, units
from ase.optimize import QuasiNewton
from ase.md import VelocityVerlet
#from ase.optimize OR ase.md import Langevin
from ase.calculators.lammpsrun import * #LAMMPS
from ase.constraints import FixAtoms
#IMPORT NEB SPECIFIC PACKAGES

#Cis
def Gradient(Elements, Image, ImageNo, GuessNo, output, dummy_variable): #output is for labeling the process


#Creates a directory for each image, at each stage

	omitoutput = os.system( 'mkdir -p  ./Molecule-%s/Image%s' % (str(GuessNo), str(ImageNo)) )
        os.chdir( 'Molecule-%s/Image%s' % (str(GuessNo), str(ImageNo)) )

#If a gradient has already been generated for a given image it will skip it
        SaveGradient = 'Gradient-%s-%s' % (str(GuessNo),str(ImageNo))


        os.system('ls --ignore \'INPUT*\' > CheckGradTemp')
        if SaveGradient in open('CheckGradTemp').read():
                os.system('rm CheckGradTemp')

                f1 = open('Gradient-%s-%s' % (str(GuessNo),str(ImageNo)), 'r')
		Ftemp = f1.readlines()
		F = []
		for i in range(0, len(Ftemp)):
			force = filter(None, Ftemp[i][2:-2].split(' '))
			Newline = []
			for k in range(0, len(force)):
				Newline.append(float(force[k]))
			F.append(Newline)		
 
                print 'Gradient obtained for Guess No %s and Image No %s' % ( str(GuessNo),str(ImageNo))
	else:		
		os.system('rm CheckGradTemp')	

# Set up the molecule by specifying atomic positions and atomic number

		Molecule = Atoms(numbers = Elements, positions=Image, pbc=(0,0,0))
                write('Molecule-%s-%s-start.xyz' % (str(GuessNo),str(ImageNo)), Molecule)
# Center the molecule in the cell with some vacuum around
		Molecule.center(vacuum=6.0)
		##FIX:
		const = FixAtoms(mask = Molecule.positions[22,:]) #blah
		Molecule.set_constraint(const)  #HERE!!! 

# Run the relaxation for each energy shift, and print out the
# corresponding total energy, bond length and angle
		starttime = time.time()
		calc = LAMMPS('Molecule', parameters = { 'pair_style' : 'airebo 3.0 1 1', 'pair_coeff' : 
                    ['* * ../../CH.airebo C H'], 'mass' : ['1 12.0107','2 1.00794']}, tmp_dir='./', keep_tmp_files=True, no_data_file=False)
                #files = ['../CH.airebo'] \n min_style fire \n  unfix fix_nve','run': '1', , 'minimize': '10e-12 10e-14 100000 100000'
	                 
		Molecule.set_calculator(calc)

		F = Molecule.get_forces()
		#print calculation_required(Molecule, 'energy')
		E = Molecule.get_potential_energy()
		

                #print 'F-%s-%s' % (str(GuessNo), str(ImageNo)), 'is', F, 'type F is; ',type(F), np.shape(F)

		f2 = open(SaveGradient, 'w').write(str(F))
		f3 = open('Energy', 'w').write(str(E))

                AngleCNNC = Molecule.get_dihedral([3, 1, 0, 2])
                AngleCCNN = Molecule.get_dihedral([7, 3, 1, 0])
                AngleNNCC = Molecule.get_dihedral([1, 0, 2, 6])

                AngleCNN = Molecule.get_angle([3, 1, 0])
                AngleNNC = Molecule.get_angle([1, 0, 2])
                AngleCCN = Molecule.get_angle([7, 3, 1])

                f4 = open('CNNCandCCNNandNNCC', 'w').write(str(AngleCNNC)+ '\t' + str(AngleCCNN) + '\t' + str(AngleNNCC))
                f5 = open('CNNandNNCandCCN', 'w').write(str(AngleCNN)+ '\t' + str(AngleNNC)+ '\t' + str(AngleCCN))

#endtime = time.time()
	#walltime = endtime - starttime
	#print 'Wall time: %.5f' % walltime
	#print                                # Make the output more readable 	
	#client.close()
	output.put((ImageNo,F))
	#return F
