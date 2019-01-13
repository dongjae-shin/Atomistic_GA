import os
from random import random
from ase.io import read, write
from ase.db import connect
from ase.visualize import view
from ase import Atoms
from ase.calculators.vasp import Vasp
#from ase.calculators.emt import EMT

from ase.ga.utilities import closest_distances_generator
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.data import DataConnection

from utilities_alloy import *

# Genetic algorithm parameters
population_size = 20
mutation_probability = 0.3
num_generations = 200
num_parents = 8
num_offspring = 8
ediff = 10e-5
nlimit = 30

pop_file_name_i = 'initial_pop.traj'
pop_file_name_f = 'final_pop.traj'
db_file_name = 'initial_pop.db'

f = open('EMIN', 'a+')
print "{0:20}:{1:10.5f}".format('population_size', population_size)
print "{0:20}:{1:10.5f}".format('mutation_probability', mutation_probability)
print "{0:20}:{1:10.5f}".format('num_generations', num_generations)
print "{0:20}:{1:10.5f}".format('num_parents',     num_parents)
print "{0:20}:{1:10.5f}".format('num_offspring',   num_offspring)
print "{0:20}:{1:10.5f}".format('ediff',           ediff)
print "{0:20}:{1:10.5f}".format('nlimit',          nlimit)

f.write("{0:20}:{1:10.5f}\n".format('population_size', population_size))
f.write("{0:20}:{1:10.5f}\n".format('mutation_probability',
                                   mutation_probability))
f.write("{0:20}:{1:10.5f}\n".format('num_generations', num_generations))
f.write("{0:20}:{1:10.5f}\n".format('num_parents',     num_parents))
f.write("{0:20}:{1:10.5f}\n".format('num_offspring',   num_offspring))
f.write("{0:20}:{1:10.5f}\n".format('ediff',           ediff))
f.write("{0:20}:{1:10.5f}\n".format('nlimit',          nlimit))


atoms = read('POSCAR')
atom_numbers = atoms.numbers
all_atom_types = get_atom_types(atoms)

# Setting calculator
queue = True
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

#calc = EMT()

# Setting operation classes
sg = StartGenerator(atoms, atom_numbers)
pairing = CoordinatePairing(atoms)
mutations = CoordinateSwapMutation(all_atom_types, atoms)

# For filtering unphysical candidate
#blmin = closest_distances_generator(all_atom_types,
#                                    ratio_of_covalent_radii=0.7)

# For diversity of the population pool
comp = InteratomicDistanceComparator(n_top=0,
                                     pair_cor_cum_diff=0.015,
                                     pair_cor_max=0.7,
                                     dE=0.02,
                                     mic=True)

# Continue from the final population
if os.path.exists(pop_file_name_f):
    population = [read(pop_file_name_f, index=i) for i in
            range(population_size)]
    print "population read from {0}.".format(pop_file_name_f)
# Creating initial population
else:
    population = [sg.get_new_candidate()]

    # Garanteeing the diversity of the pool
    size = len(population)
    count = 0
    while (size < population_size):

        new_candidate = sg.get_new_candidate()
        print "ordering: ", new_candidate.info['ordering']
        print "population size: ", size

        isNew = True
        for i in range(len(population)):
            cum_diff, max_diff = comp.__compare_structure__(new_candidate, 
                                                            population[i])
            if (cum_diff < comp.pair_cor_cum_diff
                    and max_diff < comp.pair_cor_max):
                # compare only the structure, not energy here
                isNew = False
                print "already exists"
                break

        if isNew:
            print "got new one!"
            population = population + [new_candidate]
            size = len(population)

        print count
        count += 1
    print "diverse initial population obtained!"

    print "initial population is written as {0}.".format(pop_file_name_i)
    write(pop_file_name_i, population)

#=============================================================================
#    for a in population:
#        print a.info['ordering']
#        db.write(a, diverse=True)

# Main loop
print "Main loop starts..."
Emin = []
count = 0
for generation in range(num_generations):
    print "generation: {0}".format(generation)

    # Parents Selection
    ss = SurvivorSelector(population, calc, queue=queue)
    if generation == 0:
        energies = ss._get_energies()
        print "energies:\n", energies
        Emin = Emin + [min(energies)]
        print "{0} E_min: {1:10.5f}".format(generation, Emin[generation])
        f.write("{0} E_min: {1:10.5f}\n".format(generation, Emin[generation]))
    parents = ss.get_parents(num_parents)

    # Mating parents to create a new candidate
    offspring = []
    size = len(offspring)
    while size < num_offspring:
        # Crossover
        new_candidate = pairing.cross(parents[size%num_parents],
                                      parents[(size+1)%num_parents])
        # Further mutation
        if random() < mutation_probability:
            new_candidate = mutations.mutate(new_candidate)

        # Check if there's already same structure in the population or the
        # generated offspring
        isNew = True
        for i in range(len(population)):
            cum_diff, max_diff = comp.__compare_structure__(new_candidate, 
                                                            population[i])
            if (cum_diff < comp.pair_cor_cum_diff
                    and max_diff < comp.pair_cor_max):
                # compare only the structure, not energy here
                isNew = False
                print "already exists in the population"
                break

        if isNew and size != 0:
            for j in range(len(offspring)):
                cum_diff, max_diff = comp.__compare_structure__(new_candidate, 
                                                                offspring[j])
                if (cum_diff < comp.pair_cor_cum_diff
                        and max_diff < comp.pair_cor_max):
                    # compare only the structure, not energy here
                    isNew = False
                    print "already exists in the generated offspring"
                    break

        if isNew:
            print "got new one!"
            offspring += [new_candidate]
            size = len(offspring)

#    offspring = [pairing.cross(parents[i%num_parents], 
#                               parents[(i+1)%num_parents]) 
#                 for i in range(num_offspring)]
#
#    # Further mutations
#    for a in offspring:
#        if random() < mutation_probability:
#            a = mutations.mutate(a)
            
    lowest = ss.where_are_lowest(num_offspring)
    print "lowest indices: ", lowest

    for i in range(num_offspring):
        population[lowest[i]] = offspring[i]
        print "candidate {0} was replaced.".format(lowest[i])

    # Renew populatoin
    ss = SurvivorSelector(population, calc, queue=queue)
    energies = ss._get_energies()
    print "energies:\n", energies
    Emin = Emin + [min(energies)]
    print "{0} E_min: {1:10.5f}".format(generation+1, Emin[generation+1])
    f.write("{0} E_min: {1:10.5f}\n".format(generation+1, Emin[generation+1]))

    write(pop_file_name_f, population)

    # Additional termination criterion
    if abs(Emin[generation+1] - Emin[generation]) < ediff:
        print "E_min difference is smaller than {0:10.5f}.".format(ediff)
        count += 1
    else:
        count = 0
        continue

    if count == nlimit:
        print "E_min didn't change for {0} generations.".format(count)
        break

print "latest population is written as {0}.".format(pop_file_name_f)
write(pop_file_name_f, population)
