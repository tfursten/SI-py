import numpy

def new_mut_count(p):
	"""A function which sets the mutation count down"""
	return numpy.random.geometric(p)
	
def time_to_mutate(mut_count,threshold=1):
	if mut_count < threshold:
		return True
	else:
		return False 
		
def mutate(ch1,total_alleles,locus=0,shift=0):
	ch1[locus+shift] = total_alleles
	return ch1
	
		
def decrement_count(n,m):return lambda x: (x - n*m)

def lethal_homozygous(c2_paternal,c2_maternal):
	if c2_paternal==[1] and c2_maternal == [1]: return True
	else: return False


	