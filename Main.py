import time, sys
from Population import *
from Genetics import *
from Plant import *
from math import *
import numpy as np
import random, cPickle
import Generation_Summary



if len(sys.argv) != 2:
	print('Plant Population Simulation')
	print('Usage:', argv[0], 'init')
	print('init = name of initialization file')
else:
	init_file = open(sys.argv[1],'r')
	title = init_file.readline().strip()  #title string
	maiting_system = int(init_file.readline()) #int-this is the key for the SI_functions dictionary 
	sphere = eval(init_file.readline()) #bool to determine if using grid or sphere, no longer relavent can remove
	n_pix = int(init_file.readline()) #int number of cells in the grid
	init_num_plants = int(init_file.readline()) #int initial number of plants 
	generations = int(init_file.readline()) #int number of generations to run the simulation
	data_rate = int(init_file.readline()) #int to determine how often data will be collected
	n_markers = int(init_file.readline()) #int for number of markers in the chromosomes ex: n_markers=3=[S,M1,M2,M3]
	n_pollen = int(init_file.readline()) #int for how many pollen will be produced by each plant
	n_ovules = int(init_file.readline()) #int for how many ovules each plant will have
	mutation_p = eval(init_file.readline()) #float mutation rate for S and M loci
	del_mut_p = eval(init_file.readline()) #float mutation rate for a deleterious mutation on chromosome 2
	p_mean = eval(init_file.readline()) #float for parameter of pollen dispersal
	s_mean = eval(init_file.readline()) #float for parameter of seed dispersal
	init_file.close()	



#establish width and height of grid
for i in xrange(int(ceil(sqrt(n_pix)+1)),0,-1):  
	if n_pix%i == 0:
		x = i
		y = n_pix/i
		break
x_max = x #int
y_max = y #int

def grid_disperse(mean, x, y):
	angle = random.uniform(0,2*pi) #float, chooses a random angle to represent the direction
	r = np.random.exponential(mean) #float choose a distance based on an exponential distribution and p_mean or s_mean (floats)
	x = int(round(x+r*cos(angle))) #int convert polar coordinates to rectangular and add distance to current position
	y = int(round(y+r*sin(angle))) #int 
	if x >= 0 and x < x_max and y >= 0 and y < y_max:
		return x * y_max + y
	else: return -1


disperse = grid_disperse #no longer relevant since there is only one disperse function

def recombination(plant, s):
	c = [plant.genotype_ch1()[s], plant.genotype_ch1()[not s]] #list of lists [[int,int,int...],[int,int,int...]]
	c2 = random.choice(plant.genotype_ch2()) #list with one int item [int]
	pix = plant.position() #int
	if len(c[0]) == 1:
		return [c[0],plant.parents()[s], c2, pix]  #[int,int,int...],[int,int,int...],[int], int] 
	else:
		r = random.random() #float
		if r > 0.5:
			return [c[0],plant.parents()[s], c2, pix] #[int,int,int...],[int,int,int...],[int], int] 
		else:
			n_markers = len(c[0]) #int
			i = int((n_markers) + ceil(log(r,2))) #int
			if i < 1: i = 1
			haplotype = c[0][0:i] + c[1][i:] #list [int,int,int...]
			ancestry = plant.parents()[s][0:i] + plant.parents()[not s][i:] #list [int,int,int...]
			return [haplotype, ancestry, c2, pix] #[int,int,int...],[int,int,int...],[int], int] 

def n_si(maternal_plant,paternal_plant,pollen_haplotype,dominance_rank): #(Plant object,Plant object, int(0/1), list of floats)
	'''no self incompatibility'''
	return True
	
def p_si(maternal_plant,paternal_plant, pollen_haplotype, dominance_rank):#(Plant object,Plant object, int(0/1), list of floats)
	'''physical: prevented from mating with itself'''
	m_plant = maternal_plant.position() #int
	p_plant = paternal_plant.position() #int
	if m_plant == p_plant: return False
	return True
	
def g_si(maternal_plant,paternal_plant, pollen_haplotype, dominance_rank): #(Plant object,Plant object, int(0/1), list of floats)
	'''Gametophytic self-incompatibility determined by pollen haplotype'''
	m_plant = maternal_plant.s_genotype() #tuple (int,int) method from Plant class
	pollen = paternal_plant.s_genotype()[pollen_haplotype] #int
	if m_plant[0] == pollen: return False
	if m_plant[1] == pollen: return False
	return True
	
def s_si(maternal_plant,paternal_plant, pollen_haplotype, dominance_rank):#(Plant object,Plant object, int(0/1), list of floats)
	'''Codominant Sporophytic self-incompatibility determined by diploid genotype '''
	m_plant = maternal_plant.s_genotype() #tuple (int,int)
	pollen = paternal_plant.s_genotype() #tuple (int,int)
	if m_plant[0] == pollen[0]: return False
	if m_plant[0] == pollen[1]: return False
	if m_plant[1] == pollen[0]: return False
	if m_plant[1] == pollen[1]: return False
	return True
	
def b_si(maternal_plant,paternal_plant, pollen_haplotype, dominance_rank): #(Plant object,Plant object, int(0/1), list of floats)
	'''Brassic Dom-codom Sporophytic self-incompatibility determined by domiant phenotype of pollen'''
	m_plant = maternal_plant.s_genotype() #tuple (int,int)
	pollen = paternal_plant.s_genotype() #tuple (int,int)
	if dominance_rank[pollen[0]] > dominance_rank[pollen[1]]:
		pollen = pollen[0]
	else:
		pollen = pollen[1]
	if m_plant[0] == pollen: return False
	if m_plant[1] == pollen: return False
	return True 
		
SI_functions = {0:n_si,1:p_si,2:g_si,3:s_si,4:b_si}	
compatibility = SI_functions[maiting_system] #sets which SI function will be used to establish compatibility

while n_pix < init_num_plants:
	print'Error: The starting number of plants: ' , init_num_plants , ' is greater than the number of pixels: ', n_pix
	init_num_plants = int(input("Enter a number < or = to the number of pixels: "))
	
def main(title, maiting_system, n_pix, sphere, init_num_plants, generations, n_markers, n_ovules, n_pollen, mutation_p, del_mut_p, p_mean, s_mean, data_rate):
	
	#----------------------Set Parameters-------------------------------------
	#set average time between mutations
	n_loci = n_markers + 1 #int
	s_mut_p = mutation_p #float
	if n_markers:
		m_mut_p = s_mut_p/n_markers #float

	
	#keep count of the total number of alleles
	total_s_alleles = 2*init_num_plants  #int
	total_m_alleles = n_markers * 2 * init_num_plants #int
	
	#start first mutation countdown
	s_mut_count = new_mut_count(s_mut_p) #int
	m_mut_count = new_mut_count(m_mut_p) #int
	d_mut_count = new_mut_count(del_mut_p) #int 
	
	#set functions to decrement the mutation countdown
	s_decrement_p = decrement_count(1,n_pollen)
	m_decrement_p = decrement_count(n_markers, n_pollen)
	d_decrement_p = decrement_count(1, n_pollen)
	s_decrement_o = decrement_count(1,n_ovules)
	m_decrement_o = decrement_count(n_markers, n_ovules)
	d_decrement_o = decrement_count(1, n_ovules)

	dominance_rank = [random.random() for i in xrange(total_s_alleles)] #list of floats

	#set file information
	pop_file_name = title+'.pkl' #str
	#pop_file = open(pop_file_name, 'w')
	outfile = open(pop_file_name, 'wb') #file
	
	paternal = 0
	maternal = 1
	#----------------------Initialize Population-------------------------------
	
   	pop = Grid(n_pix, n_ovules, title, x_max, y_max) #Population Grid Object
	init_plant_population(pop,init_num_plants,n_loci)
	gen = 0
	#-----------------------------(Disperse Pollen)----------------------------
	print('-'*50)
	print('Gen # \t Time')
	print('-'*50)
	gen_count = 1
	
	

	while generations:
		start_time = time.clock()
		if gen_count % data_rate == 0: collect_data = True
		else: collect_data = False
		
		for pix in xrange(n_pix):
			if pop.cell_weights()[gen][pix]: #list of floats, zero if there is no plant in location
				paternal_plant = pop.population()[gen][pix] #Plant object at location in pop
				x = pix//y_max #int
				y = pix%y_max #int
				check_S_mutation = False
				check_M_mutation = False
				check_D_mutation = False
				#Check if any pollen from this plant will have a mutation and decrement counts
				#Note: this only allows one mutation to occur per plant per loci
				if s_mut_count - n_pollen < 0:
					#determine index of pollen that will have S allele mutation
					s_mut_index = s_mut_count #int
					#reset count adding the remainder
					s_mut_count = new_mut_count(s_mut_p)+ n_pollen-s_mut_index-1 #int
					check_S_mutation = True
				else: 
					s_mut_count = s_decrement_p(s_mut_count) #int
				if m_mut_count - n_pollen*n_markers < 0:
					#determine index of pollen that will have a marker mutation
					m_mut_index = m_mut_count//n_markers #int
					marker_index = m_mut_count % n_markers #int
					#reset count adding the remainder
					m_mut_count = new_mut_count(m_mut_p) + n_pollen*n_markers - m_mut_count -1 #int
					check_M_mutation = True
				else:
					m_mut_count = m_decrement_p(m_mut_count) #int
				if d_mut_count - n_pollen < 0:
					#determine index of pollen that will have d allele mutation
					d_mut_index = d_mut_count #int
					#reset count adding the remainder
					d_mut_count = new_mut_count(del_mut_p) + n_pollen-d_mut_index-1 #int
					check_D_mutation = True
				else: 
					d_mut_count = d_decrement_p(d_mut_count) #int

				pollen_positions = [disperse(p_mean,x,y) for pollen in xrange(n_pollen)] #list of disperal positions-int. -1 if out of range
				#new list to determine if there is a plant in that position 
				occupied_positions = [pop.cell_weights()[gen][position] if position >= 0 else 0 for position in pollen_positions] #list of ints
				for i,occupied in enumerate(occupied_positions):
					if occupied: 
						new_pix = pollen_positions[i] #int
						maternal_plant = pop.population()[gen][new_pix] #plant object
						pollen_haplotype = random.randint(0,1)
						compatible = compatibility(maternal_plant,paternal_plant, pollen_haplotype, dominance_rank)
						if maiting_system == 2 and check_S_mutation and i == s_mut_index:
							compatible = True
						if compatible:
							o = random.randint(0,n_ovules-1)
							ovule_weight = pop.ovule_weights()[new_pix][o] #float
							pollen_weight = random.random() #float
							if ovule_weight < pollen_weight:
								pop.ovule_weights()[new_pix][o] = pollen_weight #float
								pollen = recombination(paternal_plant, pollen_haplotype) #[int,int,int...],[int,int,int...],[int], int]
								pop.ovules()[new_pix][o][paternal] = pollen
								if check_S_mutation:
									if i == s_mut_index:
										print pollen
										mutate(pollen[0],total_s_alleles) #([int,int,int...], int)
										print pollen
										total_s_alleles += 1
										dominance_rank.append(random.random())
								if check_M_mutation:
									if i == m_mut_index:
										mutate(pollen[0],total_m_alleles,marker_index,1) #([int,int,int...],int,int,int)
										total_m_alleles += 1
								if check_D_mutation:
									if i == d_mut_index:
										pollen[2]=[1]
										
		

	
		for pix in xrange(n_pix):
			if pop.cell_weights()[gen][pix]:
				x = pix//y_max  #int
				y = pix%y_max	#int
				maternal_plant = pop.population()[gen][pix] #plant object
				check_S_mutation = False
				check_M_mutation = False
				check_D_mutation = False
				#Check if any ovules from this plant will have a mutation and decrement counts
				#Note: this only allows one mutation to occur per plant per loci
				if s_mut_count - n_ovules < 0:
					#determine index of ovule that will have S allele mutation
					s_mut_index = s_mut_count #int
					#reset count adding the remainder
					s_mut_count = new_mut_count(s_mut_p)+ n_ovules-s_mut_index-1 #int
					check_S_mutation = True
				else: 
					s_mut_count = s_decrement_o(s_mut_count) #int
				if m_mut_count - n_ovules*n_markers < 0:
					#determine index of ovule that will have a marker mutation
					m_mut_index = m_mut_count//n_markers #int
					marker_index = m_mut_count % n_markers #int
					#reset count adding the remainder
					m_mut_count = new_mut_count(m_mut_p) + n_ovules*n_markers - m_mut_count -1 #int
					check_M_mutation = True
				else: 
					m_mut_count = m_decrement_o(m_mut_count)#int
				if d_mut_count - n_pollen < 0:
					#determine index of ovule that will have d allele mutation
					d_mut_index = d_mut_count #int
					#reset count adding the remainder
					d_mut_count = new_mut_count(del_mut_p) + n_ovules-d_mut_index-1 #int
					check_D_mutation = True
				else: 
					d_mut_count = d_decrement_o(d_mut_count) #int

				occupied_ovules = pop.ovule_weights()[pix] #list of floats [float,float,float...n_ovules]
				seed_disperse = [disperse(s_mean,x,y) if ovule else -2 for ovule in occupied_ovules] #list of ints
				for o,new_pix in enumerate(seed_disperse):
					if new_pix >= 0:
						seed_weight = random.random() #float
						cell_weight = pop.cell_weights()[not gen][new_pix] #float
						if seed_weight > cell_weight:
							ovule_haplotype = random.choice((0,1)) #int
							ovule = recombination(maternal_plant,ovule_haplotype) #[int,int,int...],[int,int,int...],[int], int]
							pop.ovules()[pix][o][maternal] = ovule
							if check_S_mutation:
								if o == s_mut_index:
									mutate(ovule[0],total_s_alleles)
									total_s_alleles += 1
									dominance_rank.append(random.random())
							if check_M_mutation:
								if o == m_mut_index:
									mutate(ovule[0],total_m_alleles,marker_index,1)
									total_m_alleles += 1
							if check_D_mutation:
								if o == d_mut_index:
									ovule[2]=[1]
							if not lethal_homozygous(pop.ovules()[pix][o][paternal][2],pop.ovules()[pix][o][maternal][2]):
								plant = Plant(new_pix) #plant object
								plant.set_genotype_ch1([pop.ovules()[pix][o][paternal][0],pop.ovules()[pix][o][maternal][0]]) #list of lists [[int,int,int...],[int,int,int...]]
								plant.set_genotype_ch2([pop.ovules()[pix][o][paternal][2],pop.ovules()[pix][o][maternal][2]]) #list of lists [[int],[int]]
								plant.set_parents(pop.ovules()[pix][o][paternal][3], pop.ovules()[pix][o][maternal][3], n_loci)
								plant.set_grandparents(([pop.ovules()[pix][o][paternal][1], pop.ovules()[pix][o][maternal][1]]))
								pop.cell_weights()[not gen][new_pix] = seed_weight
								pop.population()[not gen][new_pix] = plant
		

		print gen_count,'\t',time.clock()-start_time
		
		weights = [int(ceil(i)) if i else i for i in pop.cell_weights()[not gen]] #list of ints
		#pop_file.write(str(weights)+"\n")

		if collect_data:
			current_pop = Generation_Summary.current_plants(pop.population()[not gen], weights)
			cPickle.dump(current_pop,outfile)
		

		#reset for next generation
		pop.reset_ovule_weights()
		pop.reset_cell_weights(gen)
		gen = not gen
		generations -= 1						
		gen_count += 1
	
	#pop_file.close()	
	outfile.close()

		

main(title, maiting_system, n_pix, sphere, init_num_plants, generations,n_markers,n_ovules,n_pollen,mutation_p, del_mut_p, p_mean,s_mean,data_rate)