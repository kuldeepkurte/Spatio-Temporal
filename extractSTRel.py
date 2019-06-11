'''
This script will identify the pre-flood image and post-flood images and calls ST_QSR_v1_0 python script, which will extract all the QSR (RCC8)

Basically the following relations are established
0->1
0->2
0->3
0->4 
'''
from ST_QSR_v1_0 import *
import numpy as np
import sys
import csv
import time

program_name = sys.argv[0]
arguments = sys.argv[1:]
count = len(arguments)
#print(arguments[0])
dirPath = arguments[0]
fname = dirPath+"\ST_meta.csv"

id
imageid=''
segPath=''
ST_im_dict = {}
with open(fname) as f:
	content = csv.reader(f)
	ST_im_dict = {int(row[0]):row[1] for row in content}

#print(len(ST_im_dict))
ST_im_count = len(ST_im_dict)
preIm = ST_im_dict.get(0)
postIm =''

pre_map = np.loadtxt(dirPath+"\Segments\\"+preIm+'\map.txt', dtype='int', delimiter=' ')
post_map = np.zeros(shape=np.shape(pre_map), dtype='int')
#print(post_map.shape)

for stimid in range(1,ST_im_count):
	#print(ST_im_dict.get(stimid))
	postIm = ST_im_dict.get(stimid)
	post_map = np.loadtxt(dirPath+"\Segments\\"+postIm+'\map.txt', dtype='int', delimiter=' ')
	tic = time.time()
	overlay(pre_map,post_map, preIm, postIm, dirPath)
	toc = time.time()
	print("%s-%s %.2fs"%(preIm,postIm,(toc-tic)))

