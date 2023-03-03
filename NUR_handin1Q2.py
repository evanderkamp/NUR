import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import timeit

#exercise 2

#functions for row swapping, adding, and scaling
#I made the functions together with my sister Liz van der Kamp (s2135752) during the working classes

def rowswap(M, i, j):
	"""swap row i and j of matrix M"""
	B = np.copy(M).astype("float64")

	row = B[i, :].copy()
	B[i,:] = B[j, :]
	B[j, :] = row
	return B

def rowswapvec(M, i, j):
	"""swap indices i and j of a vector M"""
	B = np.copy(M).astype("float64")

	row = B[i].copy()
	B[i] = B[j]
	B[j] = row
	return B
    

def rowadd(M, i, j, scale):
	"""add row i to row j of matrix M scale times"""
	B = np.copy(M).astype("float64")

	row = B[j, :].copy()
	B[j, :] = row+ scale*B[i,:]
	return B


def rowscale(M, i, scale):
	"""scale row i of matrix M scale times"""
	B = np.copy(M).astype("float64")

	row = B[i, :].copy()
	B[i,:] = row*scale
	return B

def LUdecomp(M, x):
	"""LU decomposition given a matrix M and x, returns x_end satisfying M*x_end = x"""
	#make copies so we dont mess with any global variables accidentally
	M = M.copy()
	x = x.copy()
	x_new = x.copy()
	rows = M.shape[0]
	columns = M.shape[1]

	pivot_ind = np.zeros(columns).astype("int")
	index = np.linspace(0,rows-1, rows).astype("int")
	parit = 1

	#check if the matrix is singular
	for j in range(columns):
		if np.all(M[:,j] == 0):
			print("Matrix is singular")
			return M, np.zeros(rows), index

    	#get pivots by choosing the maximum and put the pivots on the diagonal
	for i in range(columns):
		maxim = np.abs(M[i,i])
		piv = i
		for k in range(i+1, rows):
			xij = np.abs(M[k,i])
			if xij > maxim:
				maxim = xij
				piv = k
                
		pivot_ind[i] = piv
		if piv != i:
			M = rowswap(M,i,piv)
			#keep track of what we swap
			index = rowswapvec(index,i,piv).astype("int")
			parit *= -1
            
		xii = M[i,i]
		#get the LU matrix
		for n in range(i+1, rows):
			LUni = M[n,i] / xii
			M[n,i] = LUni
			for m in range(i+1, rows):
				M[n,m] -= LUni * M[i,m]
        
        
	#get the solution
	x_end = np.zeros(rows)
	y_end = np.zeros(rows)

	#forward substitution
	for p in range(rows):
		ay = 0
		for l in range(0,p):
			ay += M[p,l]*y_end[l]
		y_end[p] = x_new[index[p]] -ay

	#backward substitution
	for w in range(rows):
		back = rows-(w+1)
		ax = 0
		for q in range(back,rows):
			ax += M[back,q]*x_end[q]
		x_end[back] = 1/M[back,back] * (y_end[back] -ax)
        
	
	#return the index vector so we know what we swapped
	return M, x_end, index




#get the data
data=np.genfromtxt(os.path.join(sys.path[0],"Vandermonde.txt"),comments='#',dtype=np.float64)
x=data[:,0]
y=data[:,1]

xint=np.linspace(x[0],x[-1],1001) #x values to interpolate at

#vandermonde matrix:
VM = np.empty((len(x),len(x)))
for i in range(len(x)):
	VM[:,i] = x**i


#solve for c with LU decomp

LUm, c, indexs = LUdecomp(VM, y)
print("Solution for c: ", c)

#save the output as a txt file
np.savetxt("NUR1_c_sol1.txt", c)

#get the y values for the interpolated x

xint_power = [xint**k for k in range(len(x))]
x_power = [x**k for k in range(len(x))]

yint = np.sum(c[:,None]*xint_power, axis=0)
ysol = np.sum(c[:,None]*x_power, axis=0)



#plot the results
fig,(ax1,ax2) = plt.subplots(2, figsize=(10,5), sharex=True)
ax1.scatter(x,y, label="data points")
ax1.plot(xint,yint, label="interpolated polynomial (LU)", color='crimson')
ax1.set_ylabel("y")
ax1.set_ylim(-400,400)
ax1.legend()

ydiff = np.abs(y - ysol)

ax2.scatter(x,ydiff)

ax2.set_xlabel("x")
ax2.set_ylabel(r"|y$_i$ - y(x)|")
plt.savefig("NUR1_Q2_plot1.pdf")
plt.close()

#2b

def polinterp(x, x_i, y_i):
	"""Interpolate data x_i and y_i with a polynomial of order length(x_i) using Neville's algorithm"""

	#bisection
	ys = np.zeros(len(x))
	length = len(x_i)
	M= length
	half = int(length/2)
	for j in range(len(x)):
		xj = x[j]
		start = 0
		half2 = int(length/2)
		i =2
		while half2 != start:
            
			if xj < x_i[half2]:
				half2 = int(half2-(half2-start)/2)
				i *=2
			else:
				start = half2
				half2 = half2+int(length/i)

			if half2 >= length-1:
				half2 = length-2
        
        	#make sure it goes well for odd and even M
		begin = half2-int((M-1)/2)
		end = half2+int((M-1)/2)+(1-M%2)

            	#if the closest x is at the start/end of the data set, use the first/last M data pts to interpolate with
		if end >= length-1:
			end = length-1
			begin = length-M

		if begin <= 0:
			begin= 0 
			end= M-1

        	#interpolate
		y_p = y_i.copy()
		for k in range(1,M):
			for l in range(M-k):   
				y_p[l] = ((x_i[l+k] - xj)*y_p[l] + (xj - x_i[l])*y_p[l+1])/(x_i[k+l] - x_i[l])
		ys[j] = y_p[0]
        
	return ys


y_Neville = polinterp(xint, x, y)
#also interpolate so we have a set at the data points of the same size
y_comp = polinterp(x, x, y)

fig,(ax1,ax2) = plt.subplots(2, figsize=(10,5), sharex=True)
ax1.scatter(x,y, label="data points")
ax1.plot(xint,yint, label="interpolated polynomial (LU)", color='crimson')
ax1.plot(xint,y_Neville, label="interpolated polynomial (Neville)", color='black')
ax1.set_ylabel("y")
ax1.set_ylim(-400,400)
ax1.legend()


ydiff2 = np.abs(y_comp - y)

ax2.scatter(x,ydiff2, label="difference data pts and Neville polynomial")
ax2.legend()

ax2.set_xlabel("x")
ax2.set_ylabel(r"|y$_i$ - y(x)|")
plt.savefig("NUR1_Q2_plot2.pdf")
plt.close()


#The first bit of Nevilles algorithm and the LU decomposition do not match (on the range x~0-20) 
#In LU decomposition the error with y_i is larger than in Neville's algorithm because in Neville we already iterate multiple times over solutions to reduce the error, unlike in LU decomposition


#2c



def LUiteration(M, LU, x, b, index, iters):
	"""Improves a given LU decomposition of matrix M, M*x = b, by iterating iters times"""
	#make copies so we dont mess with any global variables accidentally
	M = M.copy()
	LU = LU.copy()
	x_iter = x.copy()
	rows = M.shape[0]
	
	
	for i in range(iters):
		#sum with axis=1 gives matrix multiplication
		del_b = np.sum(M*x_iter, axis=1) - b

	#get the solution for x_end (deltax, the error in the solution)
		x_end = np.zeros(rows)
		y_end = np.zeros(rows)

	#forward substitution
		for p in range(rows):
			ay = 0
			for l in range(0,p):
				ay += LU[p,l]*y_end[l]
			y_end[p] = del_b[index[p]] -ay
	#backward substitution
		for w in range(rows):
			back = rows-(w+1)
			ax = 0
			for q in range(back,rows):
				ax += LU[back,q]*x_end[q]
			x_end[back] = (1/LU[back,back]) * (y_end[back] -ax)
	#substract the error
		x_iter -= x_end
	
	#return the index vector so we know what we swapped
	return x_iter


c_improv = LUiteration(VM, LUm, c, y, indexs, iters=10)

print("10x iterated solution for c:", c_improv)
#save the output as a txt file
np.savetxt("NUR1_c_sol10.txt", c_improv)



yint10 = np.sum(c_improv[:,None]*xint_power, axis=0)
ysol10 = np.sum(c_improv[:,None]*x_power, axis=0)


fig,(ax1,ax2) = plt.subplots(2, figsize=(10,5), sharex=True)
ax1.scatter(x,y, label="data points")
ax1.plot(xint,yint, label="interpolated polynomial (LU1)", color='crimson')
ax1.plot(xint,yint10, label="interpolated polynomial (LU10)", color='black')
ax1.set_ylabel("y")
ax1.set_ylim(-400,400)
ax1.legend()


ydiff3 = np.abs(ysol10 - y)

ax2.scatter(x,ydiff3, label="difference data pts and LU (10 iters)")
ax2.legend()

ax2.set_xlabel("x")
ax2.set_ylabel(r"|y$_i$ - y(x)|")
plt.savefig("NUR1_Q2_plot3.pdf")
plt.close()


#2d

#time how long it takes to solve everything 10x
number = 20


starttime_a = timeit.default_timer()
for t in range(number):
	LUm, c, indexs = LUdecomp(VM, y)
aver_ta = (timeit.default_timer() - starttime_a)/number

print("The average time taken to do LU decomposition (2a) :", aver_ta)


starttime_b = timeit.default_timer()
for t in range(number):
	y_Neville = polinterp(xint, x, y)
aver_tb = (timeit.default_timer() - starttime_b)/number

print("The average time taken to do Neville's algorithm (2b) :", aver_tb)


starttime_c = timeit.default_timer()
for t in range(number):
	c_improv = LUiteration(VM, LUm, c, y, indexs, iters=10)
aver_tc = (timeit.default_timer() - starttime_c)/number

print("The average time taken to do improve LU solution 10x (2c) :", aver_tc)


#The LU decomposition is and improving it 10x take about the same time, which is both more than 10x faster than Neville's algorithm, so doing the LU decomposition + improvement is most efficient because for Neville you have to do bisection for every point and then interpolate and iterate many times, which is slow if you want to do it for many points
#LU decomposition just uses 1 matrix you have to calculate once and thus improving the solution is efficient because you can use the same matrix every time
#From the plots, Neville's algorithm yields the most accurate solution because it already iterates a lot of times to more precisely calculate the polynomial

np.savetxt("NUR1_Q2_avertimes.txt", [aver_ta, aver_tb, aver_tc])


