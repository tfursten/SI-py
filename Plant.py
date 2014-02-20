import random

class Plant(object):
	def __init__(self,pixel):
		"""initialize the plants and their positions. """
		self.__pixel = pixel
		self.__genotype_ch1 = [[],[]]
		self.__genotype_ch2 = [[],[]]
		self.__parents = []
		self.__grandparents = [[],[]]
		
		
	def __str__(self):
		return str(self.__genotype_ch1)
		
	def set_genotype_ch1(self,chromosomes):
		self.__genotype_ch1 = chromosomes

	def set_genotype_ch2(self,chromosomes):
		self.__genotype_ch2 = chromosomes

	def set_parents(self,pix_dad, pix_mom, n_loci):
		dad = [pix_dad for i in xrange(n_loci)]
		mom = [pix_mom for i in xrange(n_loci)]
		self.__parents = [dad,mom]

	def set_grandparents(self, grandparents):
		self.__grandparents = grandparents

	def genotype_ch1(self):
		return self.__genotype_ch1

	def genotype_ch2(self):
		return self.__genotype_ch2

	def s_genotype(self):
		return self.__genotype_ch1[0][0],self.__genotype_ch1[1][0]
	
	def parents(self):
		return self.__parents

	def grandparents(self):
		return self.__grandparents

	def position(self):
		return self.__pixel


		
	