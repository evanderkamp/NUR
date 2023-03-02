import numpy as np
import matplotlib.pyplot as plt
import timeit

#exercise 1

def Poisson(l, k): 
	"""A function that returns the Poisson probability for an integer k and a mean position l"""
	#convert l and k to 32 bits for low memory usage
	l = l.astype(np.float32)
	k = k.astype(np.int32)
	
	Pois = np.zeros(len(l))
	
	for i in range(len(l)):
		#if k is too big (>20 we get an overflow error), go to logspace 
		#use ln(k!) ~ ln(1) + ln(2) + ... + ln(k)
		
		if k[i] > 20:
			logk_fact = np.sum(np.log(np.linspace(1,k[i],k[i])))
			logP = np.float32((np.log(l[i])*k[i]  -l[i]) - logk_fact)
			Pois[i] = np.exp(logP)
		else:
			Pois[i] = np.float32(l[i]**k[i] * np.exp(-l[i]) / np.math.factorial(k[i]))
		
	return Pois


ls = np.array([1,5,3,2.6,101])
ks = np.array([0,10,21,40,200])


Poissout = Poisson(ls, ks)

#save the output as a txt file
np.savetxt("NUR1_Poisson.txt", np.transpose(ls,ks,Poissout))



