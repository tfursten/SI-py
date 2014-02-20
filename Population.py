import Plant
#import Grid_Graph
import random
import numpy as np
from math import *





class Population(object):
	def __init__(self, n_pix, n_ovules, title):
		self.__n_pix = n_pix
		self.__n_ovules = n_ovules
		self.__population = [[None for x in xrange(n_pix)], [None for x in xrange(n_pix)]]
		self.__cell_weights = [[0 for i in xrange(n_pix)], [0 for i in xrange(n_pix)]]
		self.__ovule_weights = [[0 for i in xrange(n_ovules)] for j in xrange(n_pix)]
		self.__ovules = [[[None,None] for x in xrange(n_ovules)] for y in xrange(n_pix)]  
		self.__title = title

	def __str__(self):
		return str(self.__population)
		
	def reset_cell_weights(self,gen):
		self.__cell_weights[gen] = [0 for i in xrange(self.__n_pix)]
	
	def reset_ovule_weights(self):
		self.__ovule_weights = [[0 for i in xrange(self.__n_ovules)] for j in xrange(self.__n_pix)]

	def population(self):
		return self.__population
	
	def cell_weights(self):
		return self.__cell_weights
		
	def ovule_weights(self):
		return self.__ovule_weights
		
	def ovules(self):
		return self.__ovules

	def title(self):
		return self.__title

	def n_pix(self):
		return self.__n_pix
		

class Grid(Population):
	def __init__(self, n_pix, n_ovules, title, x_max, y_max):
		Population.__init__(self, n_pix, n_ovules, title)
		self.__x_max = x_max
		self.__y_max = y_max

	#def projection(self, file_name, title):
		#Grid_Graph.grid(self.__x_max, self.__y_max, file_name, title)

	def x_max(self):
		return self.__x_max

	def y_max(self):
		return self.__y_max





		

		
#------------------------------------------------------------------------------------------------------		


	
def init_genotype(n, n_loci, init_num_plants):
	nx2_plants = 2*init_num_plants
	c1 = [n+i for i in range(0,nx2_plants*n_loci, nx2_plants)]
	c2 = [n+i+1 for i in range(0,nx2_plants*n_loci, nx2_plants)]
	return [c1,c2]


			
def init_plant_population(population,init_num_plants, n_loci):
	"""a function to initialize the first generation of plants"""
	n_pix = population.n_pix()
	pop = population.population()[0]
	weight = population.cell_weights()[0]
	location = random.sample(xrange(n_pix), init_num_plants)
	#location = random.sample(xrange(4*(n_pix/12),8*(n_pix/12)), init_num_plants)
	i = 0
	for pix in location:
		plant = Plant.Plant(pix)
		plant.set_genotype_ch1(init_genotype(i, n_loci, init_num_plants))
		plant.set_genotype_ch2([[0],[0]])
		plant.set_parents(pix,pix,n_loci)
		pop[pix] = plant
		weight[pix] = random.random()
		i += 2
	population.population()[0] = pop
	population.cell_weights()[0] = weight





