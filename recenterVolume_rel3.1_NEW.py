#!/usr/bin/env python

import optparse
import itertools
import math
import os,sys
import glob
#import heapq
import numpy as np

def setupParserOptions():
	parser = optparse.OptionParser()
	parser.set_usage("%prog -f run_data.star")
	parser.add_option("-f",dest="f1",type="string",metavar="FILE",
		help=".star file")
	parser.add_option("-s", dest="shft", type="int", metavar="INT",
		help="multiple of 8nm to shift")

	options,args = parser.parse_args()

	if len(args) > 1:
		parser.error("Unknown commandline options: " +str(args))

	if len(sys.argv) < 2:
		parser.print_help()
		sys.exit()

	params={}

	for i in parser.option_list:
		if isinstance(i.dest,str):
			params[i.dest] = getattr(options,i.dest)
	return params
	
def getMat(rot, tilt, psi):
	[rot, tilt, psi]=map(np.radians,[rot, tilt, psi])

	M1=np.array([ [np.cos(rot), -np.sin(rot), 0],\
		[np.sin(rot), np.cos(rot), 0],\
		[0, 0, 1]])

	M2=np.array([ [np.cos(tilt), 0, -np.sin(tilt)],\
		[0, 1, 0],\
		[np.sin(tilt), 0,  np.cos(tilt)]])

	M3=np.array([ [np.cos(psi), -np.sin(psi), 0],\
		[np.sin(psi), np.cos(psi), 0],\
		[0, 0, 1]])

	return np.matmul(M3,np.matmul(M2,M1))

def recenter(rot,tilt,psi,offX,offY,d=[0,0,0]):
	M=getMat(-rot,-tilt,-psi)
	d1= np.matmul(M,d)
	return (offX+d1[0],offY+d1[1])

def mainloop(params):
	f1 = file(params['f1'])
	l1 = f1.readlines()

        shft = params['shft']
        d=[0,0,shft*82.5]

	fout = file("%s_shft%dnm.star"%(params['f1'][:-5],shft*8),"w")

	# first write out the header
	#for i in range(100):
	#	t1 = l1[i].split()
	#	if len(t1) < 13:
	#		fout.write("%s"%l1[i])

	for line in l1:
		l=line.split()
		if len(l)>1:
			if l[1][0]=='#':
				if l[0]=='_rlnAngleRot': rotN=int(l[1].split('#')[1])-1
				elif l[0]=='_rlnAngleTilt': tiltN=int(l[1].split('#')[1])-1
				elif l[0]=='_rlnAnglePsi': psiN=int(l[1].split('#')[1])-1
				elif l[0]=='_rlnOriginXAngst': oxN=int(l[1].split('#')[1])-1
				elif l[0]=='_rlnOriginYAngst': oyN=int(l[1].split('#')[1])-1
			else:
				if len(l)>15:
					(l[oxN],l[oyN])=recenter(float(l[rotN]),float(l[tiltN]),float(l[psiN]),float(l[oxN]),float(l[oyN]),d)
					
		for item in l:
			fout.write(str(item)+" ")
		fout.write("\n")

	f1.close()
	fout.close()



if __name__ == "__main__":
	params = setupParserOptions()
	mainloop(params)
