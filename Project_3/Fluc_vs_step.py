import numpy as np
import matplotlib.pyplot as plt
import subprocess
import os.path      #for checking if a file exist

def run_lmp(t_stp):
#could write some short code to check that the dtvar is commented out but not going to
    run_str = "lmp -var dt " + str(t_stp) + " < in.Si-VACF"
    print run_str
    subprocess.call(run_str, shell=True)
    #p = subprocess.Popen(run_str, ) 

def make_var(Name):
    print Name
    subprocess.call(["./Python-Tools-For-LAMMPS/log2txt.py log.lammps " + Name +" Time TotEng"], shell=True)      #conver log file to a text file
    """with open ('hessian.d', "r") as myfile:
        for line in myfile:
            h.append(list(np.fromstring(line, dtype=float, sep=' ')))"""
def get_var(Name):
    m = np.loadtxt(Name)         #much more efficient then python open and line read above
    #print m
    #print m[:,1]
    DE = np.std(m[:,1])
    print "DE: ", DE
    return DE


def fulc_plot():
    F = 0.0005
    Fstep = 0.0001
    Fstop = 0.016
    F_its = []
    Vars = []
    while True:
        print F
        Name = "./data/time_v_Etot_dt-" + str(F)[2:] + ".d"         
        if not os.path.exists(Name):
            print "testing"
            run_lmp(F)
            make_var(Name)

        F_its.append(F) 
        Vars.append(get_var(Name))
        print "#####################################################"
        print "#####################################################"
        print "#####################################################" 
        F += Fstep
        if F > Fstop:
            break
    
    # Plot FvsVars
    F_its = np.array(F_its)
    Vars = np.array(Vars)

    print "F_its: ", F_its
    print "Vars: ", Vars
    Vars = Vars/(3*216*(1./1605)*300)           #deg freedom * numb atoms /kb *temp
    print "Vars: ", Vars

    ax = plt.subplot2grid((1,1),(0,0))
    ax.plot(np.log10(F_its), np.log10(Vars),'o')
    ax.set_xlabel('log(time step[ps])')
    ax.set_ylabel('log(temp fluctuation)')
    #ax.set_xscale('symlog')
    #ax.set_yscale('symlog')
    plt.show()


fulc_plot()
