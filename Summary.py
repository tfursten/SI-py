import sys, cPickle
from math import *
from collections import Counter
import pygame as pg
import numpy as np
import matplotlib.pyplot as plt


if len(sys.argv) > 3:
	print('Plant Population Simulation')
	print('Usage:', argv[0], 'init')
	print('init = name of initialization file')
else:
	init_file = open(sys.argv[1],'r')
	title = init_file.readline().strip()
	maiting_system = int(init_file.readline())
	sphere = eval(init_file.readline())
	n_pix = int(init_file.readline())
	init_num_plants = int(init_file.readline())
	generations = int(init_file.readline())
	data_rate = int(init_file.readline())
	n_markers = int(init_file.readline())
	n_pollen = int(init_file.readline())
	n_ovules = int(init_file.readline())
	mutation_p = eval(init_file.readline())
	del_mut_p = eval(init_file.readline())
	p_mean = eval(init_file.readline())
	s_mean = eval(init_file.readline())
	init_file.close()
	run_animation = bool(eval(sys.argv[2]))
print run_animation

for i in xrange(int(ceil(sqrt(n_pix)+1)),0,-1):  
	if n_pix%i == 0:
		x = i
		y = n_pix/i
		break
x_max = x
y_max = y

n_loci = n_markers + 1

def grid(x_max,y_max,grid, name):
	
	#Define some colors
	black = (   0,   0,   0)
	white = ( 255, 255, 255)
	green = (   0, 255,   0)
	red   = ( 255,   0,   0)
	
	#Set the width and height of each grid location
	width = 10
	height = 10
	#set the margin
	margin = 1
	
	
	#Initialize pygame
	pg.init()
	
	#Set the height and width of the screen
	size = [(x_max*(height+margin)+1),(y_max*(width+margin)+1)]
	screen = pg.display.set_mode(size)
	
	#Set title of the screen
	gen = 1
	pg.display.set_caption(name + str(gen))
	
	#Loop until the user clicks the close button.
	done = False
	
	#Manage how fast the screen updates
	clock = pg.time.Clock()
	
	game_array = [np.array(i).reshape((x_max,y_max)).tolist() for i in grid]
	#-----------Main Program Loop---------------

	
	'''while done == False:
		for event in pg.event.get(): #User did something
			if event.type == pg.QUIT: #If user clicked close
				done = True #Flag that we are done so we exit this loop'''
	while gen <= len(game_array):
							
	#------------------------DRAW----------------------------
			
		for i in xrange(len(game_array)):	
			screen.fill(black)
			pg.display.set_caption(name + str(gen))	
			for x in xrange(x_max):
				for y in xrange(y_max):
					color = white
					if game_array[i][x][y] == 1:
						color = green
					if game_array[i][x][y] == 2:
						color = red
					pg.draw.rect(screen, color, [(margin+width)* x + margin,(margin+height)* y + margin, width,height])
			
	
	#update screen with what has been drawn
			
			pg.display.flip()	
			pg.image.save(screen, name+str(gen)+'.png')
			gen += 1
			pg.time.wait(300)
		
		
	pg.quit()

#---------------------------------------------------------------------------------------------------
def effective_alleles(allele_dict):
	A_e_list = []
	for i in xrange(n_loci):
		counts = allele_dict[i].values()
		total_alleles = sum(counts)
		freq = [i/float(total_alleles) for i in counts]
		A_e = 1/sum([i**2 for i in freq])
		A_e_list.append(A_e)
	return A_e_list

def allele_counter(allele_dict,allele):
	if allele in allele_dict:
		allele_dict[allele] += 1
	else:
		allele_dict[allele] = 1


t='\t'
infile_name = title+'.pkl'
infile = open(infile_name, 'rb')
outfile_1_name = title+'.txt'
outfile_2_name = title+'alleles.pkl'
outfile1 = open(outfile_1_name, 'w')
outfile2 = open(outfile_2_name, 'wb')
if run_animation:
	pop_grids = []
	ch2_grids = []
outstr = '#Generation \t n_Pop \t n_S_alleles \t'
for i in xrange(n_markers):
	outstr +='Marker_'+str(i+1)+t
outstr += 'S_Ae \t'
for i in xrange(n_markers):
	outstr += 'M'+str(i+1)+'_Ae'+t
outstr += 'pIDB_S \t'
for i in xrange(n_markers):
	outstr += 'pIDB_'+str(i+1)+t
outstr += '\n'
outfile1.write(outstr)
for gen in xrange(generations/data_rate):
	plant_list = cPickle.load(infile)
	pop_n = len(plant_list)
	if run_animation:
		pop_grid = [0 for i in xrange(n_pix)]
		ch2_genotype_grid = [0 for i in xrange(n_pix)]
	loci_dict = {loci:{} for loci in xrange(n_loci)}
	IBD = [0 for i in xrange(n_loci)]
	for plant in plant_list:
		ch2_genotype = plant.genotype_ch2()
		ch1_genotype = plant.genotype_ch1()
		grandparents = plant.grandparents()
		if run_animation:
			position = plant.position()
			pop_grid[position] = 1
			if ch2_genotype == [[0],[0]]:
				ch2_genotype_grid[position] = 1
			elif ch2_genotype == [[1],[0]]:
				ch2_genotype_grid[position] = 2
			elif ch2_genotype == [[0],[1]]:				
				ch2_genotype_grid[position] = 2
		for chromosome in ch1_genotype:
			for loci in xrange(n_loci):
				allele_counter(loci_dict[loci],chromosome[loci])
		for i in xrange(n_loci):
			if ch1_genotype[0][i] == ch1_genotype[1][i]:
				if grandparents[0][i] == grandparents[1][i]:
					IBD[i] += 1

	p_IBD = [i/float(pop_n) for i in IBD]

	if run_animation:
		#pop_grids.append(pop_grid)
		ch2_grids.append(ch2_genotype_grid)
	A_e = effective_alleles(loci_dict)
	outdata = str((gen+1)*data_rate)+t+str(pop_n)+t
	for i in xrange(n_loci):
		outdata += str(len(loci_dict[i]))+t
	for i in xrange(n_loci):
		outdata += str(A_e[i]) + t
	for i in xrange(n_loci):
		outdata += str(p_IBD[i]) + t
	outdata += '\n'
	outfile1.write(outdata)

	cPickle.dump(loci_dict, outfile2)
if run_animation:
	grid(x_max,y_max, pop_grids, title+'pop')
	grid(x_max,y_max, ch2_grids, title+'ch2')
infile.close()
outfile1.close()
outfile2.close()






