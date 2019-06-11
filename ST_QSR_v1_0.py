'''
This code calculates all RCC8 spatial relations on two sample maps
ST_QSR_v1_0 (version 1.0)
'''
from skimage.measure import regionprops
import numpy as np
import json
import pickle
from scipy.spatial.distance import cdist
#import matplotlib.pyplot as plt

#ST_Rel = {}

#-----------------------------------------------------------------------
def overlay(x,y, preIm, postIm, archPath):
	
	xshape_props=regionprops(x+1)
    
	
	
	ST_Rel = {"PreImage":preIm,"PostImage":postIm}
	#print(ST_Rel)
	[rows, cols] = x.shape
	#print rows
	#print cols

	min_x = np.min(x)
	min_y = np.min(y)

	max_x = np.max(x)
	max_y = np.max(y)

	#print max_x
	overlay_mat = np.zeros(shape=(max_x+1,max_y+1), dtype='int')
	#overlay_mat = [[0 for i in range(max_x+1)] for i in range(max_y+1)]
	#print(overlay_mat.shape)

	for row in range(0,rows):
		for col in range(0,cols):
			r = x[row][col]
			c = y[row][col]
			overlay_mat[r][c]+=1

	#print(np.matrix(overlap_mat))
	#print np.sum(overlap_mat)
	#print "a b -> overlay (Area)"
	#print "--------------------"

   
	for row in range(min_x, max_x+1):
		for col in range(min_y, max_y+1):
			#print "O(a_%d, b_%d):%d"%(row,col,overlap_mat[row][col])
			
			s_area = xshape_props[row].area #get the area of a source region
			
			if overlay_mat[row][col] == 0:
				[SR, count] = refineDR(x,y,row,col)
				if row in ST_Rel.keys():
					ST_Rel[row][SR].append((col,count))
				else:
					ST_Rel[row] = {"ST_DC":[],"ST_EC":[],"ST_PO":[],"ST_NTPP":[],"ST_TPP":[],"ST_NTPPi":[],"ST_TPPi":[],"ST_EQ":[]}
					ST_Rel[row][SR].append((col,count))
				#print("%s %d %d %d"%(SR, row,col,count))   
			else:
				SR = refineOR(x,y,row,col)#row-> source region at t-1 and col-> target region at t
				#print("%s %d %d %d"%(SR,row,col,overlay_mat[row][col]))
				overlay_fraction = float(overlay_mat[row][col])/s_area
				if row in ST_Rel.keys():
					ST_Rel[row][SR].append((col,overlay_fraction))
				else:
					ST_Rel[row] = {"ST_DC":[],"ST_EC":[],"ST_PO":[],"ST_NTPP":[],"ST_TPP":[],"ST_NTPPi":[],"ST_TPPi":[],"ST_EQ":[]}
					ST_Rel[row][SR].append((col,overlay_fraction))
          
	#print(ST_Rel[0])
	with open(archPath+'\\'+preIm+'_'+postIm+'_ST_Rel.pkl', 'wb') as fhandle:
		pickle.dump(ST_Rel, fhandle)
	
	fhandle.close()
	
	#np.save(preIm+'_'+postIm+'_ST_Rel.npy', ST_Rel) 
	'''
	with open('ST_Rel.txt', 'rb') as handle:
		b = pickle.loads(handle.read())
	print(b[0]['ST_EC'])
	print(b[1]['ST_PO'])
	print(b[2]['ST_NTPPi'])
	print(b[3]['ST_DC'])
	'''
#-----------------------------------------------------------------------
def getSetClosure(s):
	ts = set(s)
	
	for x in s:
		xl = (x[0]-1, x[1])
		xr = (x[0]+1, x[1])
		xt = (x[0], x[1]-1) 
		xb = (x[0], x[1]+1)
		xtl = (x[0]-1, x[1]-1)
		xtr = (x[0]+1, x[1]-1)
		xbl = (x[0]-1, x[1]+1)
		xbr = (x[0]+1, x[1]+1)
		
		if xl not in s:
			ts.add(xl)
		if xr not in ts:
			ts.add(xr)
		if xt not in ts:
			ts.add(xt)
		if xb not in ts:
			ts.add(xb)
		if xtl not in ts:
			ts.add(xtl)
		if xtr not in ts:
			ts.add(xtr)
		if xbl not in ts:
			ts.add(xbl)
		if xbr not in ts:
			ts.add(xbr)
	
	return ts



def getPixCount_EC(s,t):
	return len(s.intersection(t))
	
#-----------------------------------------------------------------------
#routines for overlap and proper-parts
def isPart(s,t):
    if s.issubset(t) and len(s)!=0:
       return True
    else:
       return False
       
       
      
def isProperPart(s,t):
    if isPart(s,t) and not(s==t):
       #print "s is a proper part of t"
       return True
    else:
       return False
	   

def isTangentialProperPart(s,t):
    if isProperPart(s,t) and len(getSetClosure(s).difference(t))!=0:
       #print "s is a proper part of t"
       return True
    else:
       return False


def isTangentialProperPartIn(s,t):
    if isTangentialProperPart(t,s):
       #print "s is a proper part of t"
       return True
    else:
       return False
	   

def isNonTangentialProperPart(s,t):
    if isProperPart(s,t) and len(getSetClosure(s).difference(t))==0:
       #print "s is a proper part of t"
       return True
    else:
       return False

def isNonTangentialProperPartIn(s,t):
    if isNonTangentialProperPart(t,s):
       #print "s is a proper part of t"
       return True
    else:
       return False


def isProperPartI(s,t):
    if isPart(t,s) and not(s==t):
       #print "s is a proper part of t"
       return True
    else:
       return False


def isOverlap(s,t):
    if len(s.intersection(t))!=0:
       #print "s overlaps t"
       return True
    else:
       return False


def isPartialOverlap(s,t):
    if isOverlap(s,t) and not isPart(s,t) and not isPart(t,s):
       #print "s partially overlaps t"
       return True
    else:
       return False



def isEqual(s,t):
    if isPart(s,t) and isPart(t,s):
       #print "s equal t"
       return True
    else:
       return False

	   
	   

	

#-----------------------------------------------------------------------
def refineOR(x,y,s,t):
	xshape_props=regionprops(x+1)
	yshape_props=regionprops(y+1)
    
	s_coords = xshape_props[s].coords
	t_coords = yshape_props[t].coords
    #print s_coords
    #print t_coords
	s_set = set(tuple(r) for r in s_coords)
	t_set = set(tuple(r) for r in t_coords)
    
    #print s_set
    #print t_set
    
    #t_set = set(t_coords)   #print( s_set.intersection(t_set) )
    #print(s_set.issubset(t_set))
    #print(t_set.issubset(s_set))
    #print(t_set == s_set)
	
	if isPartialOverlap(s_set, t_set):
		return "ST_PO"
	elif isTangentialProperPart(s_set, t_set):
		return "ST_TPP"
	elif isNonTangentialProperPart(s_set, t_set):
		return "ST_NTPP"
	elif isTangentialProperPartIn(s_set, t_set):
		return "ST_TPPi"
	elif isNonTangentialProperPartIn(s_set, t_set):
		return "ST_NTPPi"
	elif isEqual(s_set, t_set):
		return "ST_EQ"
#-----------------------------------------------------------------------

def refineDR(x,y,s,t):
	xshape_props=regionprops(x+1)
	yshape_props=regionprops(y+1)
    
	s_coords = xshape_props[s].coords
	t_coords = yshape_props[t].coords
    #print s_coords
    #print t_coords
	s_set = set(tuple(r) for r in s_coords)
	t_set = set(tuple(r) for r in t_coords)

	
	s_set_cl = getSetClosure(s_set)
	
	if len(s_set_cl.intersection(t_set))==0:
		dist = cdist(s_coords, t_coords, metric='sqeuclidean') 
		min_dist = dist.min()**0.5
		return ["ST_DC", min_dist]
	elif len(s_set.intersection(t_set))== 0 and len(s_set_cl.intersection(t_set))>0:
		c = getPixCount_EC(s_set_cl, t_set)
		return ["ST_EC", c]
	
#-----------------------------------------------------------------------

'''

#main routine starts
#Example labelled arrays a and b

a = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0],
              [0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0],
              [0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0],
              [0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0],
              [0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0],
              [0, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 0],
              [0, 0, 3, 3, 3, 3, 2, 2, 2, 2, 2, 0],
              [0, 0, 3, 3, 3, 3, 2, 2, 2, 2, 2, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

b = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, 1, 1, 3, 3, 3, 3, 3, 3, 0],
              [0, 0, 1, 1, 3, 3, 3, 3, 3, 3, 3, 0],
              [0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 0, 0],
              [0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 0, 0],
              [0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 0, 0],
              [0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
			  
c = np.array([[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 1, 1, 1, 3, 3, 3, 3, 3, 3, 0],
              [0, 0, 1, 1, 3, 3, 3, 3, 3, 3, 3, 0],
              [0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 0, 0],
              [0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 0, 0],
              [0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 0, 0],
              [0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
              [4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0],
              [4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0],
              [4, 4, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0]])


#overlay(a,c, 'a', 'c','D:\\27_Programming\\OBIA_python\\Modules\\test')
overlay(a,b, 'a', 'b','D:\\27_Programming\\OBIA_python\\Modules\\test')

'''











