import os, datetime
from ase.io import read
from ase.calculators.vasp import Vasp
from ase import Atoms
from multiprocessing import Pool

ceria = read('POSCAR_CeO2')
aual2 = read('POSCAR_AuAl2')
calc = Vasp(xc='pbe',
            istart=0,
            icharg=2,
            ispin=2,
            encut=400.,
            prec='Normal',
            ediff=1.e-4,
            gga='PE',
            ismear=1,
            sigma=0.2,
            lreal=False,
            isym=1,
            nsw=0,
            ibrion=-1,
            potim=0.5,
            isif=2,
            ediffg=-0.03)

ceria.set_calculator(calc)
aual2.set_calculator(calc)

def worker(atoms):
    now = datetime.datetime.now().strftime("%Y_%m_%d-%H%M%S_%f")
    os.system('mkdir {0}'.format(now))
    os.chdir('{0}'.format(now))
    os.system('cp -rf ../run_slurm_vasp.sh .')
    energy = atoms.get_potential_energy()
    #os.chdir('..')
    return energy

p = Pool(2)
energies = p.map(worker, [ceria, aual2])

print "calculations done!"
print energies
