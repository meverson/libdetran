"""




"""

from __future__ import division

def get_dlop(coarse_structure):
	""" 

	coarse_structure	

        dlp = {}

	rhos[0][0]


    DLP group numbers are local to coarse groups

	"""


	dlp = {}
	for n,num in enumerate(coarse_structure):		
		for k in range(num):
			dlp[(k,0,n)] = 1.0
			dlp[(k,1,n)] = (num-1.0 - 2.0*k) / (num-1.0)
		for m in range(1,num-1):
			for k in range(num):
				c0=m*(num-1+m+1)
				c1=(2*m+1)*(num-1-2*k)
				c2=(m+1)*(num-1-m)
				dlp[(k,m+1,n)]=(c1*dlp[(k,m,n)]-c0*dlp[(k,m-1,n)])/c2

	rhos = []
	for i,num in enumerate(coarse_structure):
		rhos.append([])
		rhos[-1].append(num)
		for m in range(num):
			if m == 0: continue
			prod = 1.0
			for mm in range(1,m+1):
				prod *= (num+mm)/(num-mm)
			rho_m = num/(2.0*m+1.0)*prod
			rhos[-1].append(rho_m) 	


	return dlp,rhos

