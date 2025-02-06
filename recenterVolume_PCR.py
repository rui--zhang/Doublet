#!/usr/bin/env python

import optparse
import itertools
import math
import os,sys
import glob
#import heapq
import numpy as np
from EMAN2 import *

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
	
	
def rotate(psi,theta,phi,matrix):
        t1 = Transform({"type":"spider","psi":psi,"theta":theta,"phi":phi})
        #v = Transform([0.99738439,0.00085389,0.07227484,-21.43914131,-0.00921683,0.99326980,0.11545630,-25.72880227,-0.07168982,-0.11582045,0.99067966,207.54504114])
        v = Transform(matrix)
        t2 = t1*v.inverse()
        d = t2.get_params("spider")
        return d["psi"],d["theta"],d["phi"]

def mainloop(params):
	f1 = file(params['f1'])
	l1 = f1.readlines()

        shft = params['shft']
	# change 160 to 138
        d=[0,0,shft*138]

	fout = file("%s_tmp.star"%(params['f1'][:-5]),"w")

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
				if len(l)>10:
					(l[oxN],l[oyN])=recenter(float(l[rotN]),float(l[tiltN]),float(l[psiN]),float(l[oxN]),float(l[oyN]),d)
					
		for item in l:
			fout.write(str(item)+" ")
		fout.write("\n")

	f1.close()
	fout.close()
	
	
def fre2relion(params):
	shft = params['shft']
        matrix1 = [0.99738439,0.00085389,0.07227484,-21.43914131,-0.00921683,0.99326980,0.11545630,-25.72880227,-0.07168982,-0.11582045,0.99067966,207.54504114]
        matrix2 = [0.99737468,-0.00925357,-0.07182003,36.06412169,0.00087058,0.99326215,-0.11588593,49.61470832,0.07240848,0.11551917,0.99066257,-201.14123756]

	f2 = file("%s_tmp.star"%(params['f1'][:-5]))
	ll2 = f2.readlines()
	l2 = [x for x in ll2 if len(x.split())>11]

	fout2 = file("%s_PCR_v3_%dnm.star"%(params['f1'][:-5],shft*16),"w")

	n2 = len(l2)

	# first write out the header
	for i in range(100):
		t2 = ll2[i].split()
		if len(t2) < 15:
			fout2.write("%s"%ll2[i])
	
	#if (num_df1 != 10) or (num_df2 != 11) or (num_angast != 12) or (num_amp != 17):
	#	print "Error, columns don't match!"

	#apix = params['apix']

#_rlnImageName #1
#_rlnMicrographName #2
#_rlnCoordinateX #3
#_rlnCoordinateY #4
#_rlnAngleRot #5
#_rlnAngleTilt #6
#_rlnAnglePsi #7
#_rlnOriginXAngst #8
#_rlnOriginYAngst #9
#_rlnDefocusU #10
#_rlnDefocusV #11
#_rlnDefocusAngle #12
#_rlnPhaseShift #13
#_rlnCtfBfactor #14
#_rlnOpticsGroup #15
#_rlnRandomSubset #16
#_rlnClassNumber #17
#_rlnHelicalTubeID #18

	for i in range(n2):
		t2 = l2[i].split()
		psi = float(t2[6])
		theta = float(t2[5])
		phi = float(t2[4])
		# this is the key, should be <
		if shft < 0:
			psi_,theta_,phi_ = rotate(psi,theta,phi,matrix1)
		else:
			psi_,theta_,phi_ = rotate(psi,theta,phi,matrix2)
#		fout.write('%s\t%s\t%s\t%s\t%s\t %.6f\t%.6f\t%.6f\t %s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n'%(t2[0],t2[1],t2[2],t2[3],t2[4],df1,df2,angast,t2[8],t2[9],t2[10],t2[11],t2[12],t2[13],t2[14],t2[15],t2[16],1,phi,theta,psi,-1*shx/apix,-1*shy/apix))
		#fout.write('%s\t%s\t%s\t%s\t%s\t %s\t%s\t%s\t%s\t%s\t %s\t%s\t%s\t%d\t%.6f\t %.6f\t%.6f\t%.6f\t%.6f\n'%(t2[0],t2[1],t2[2],t2[3],t2[4],t2[5],t2[6],t2[7],t2[8],t2[9],t2[10],t2[11],t2[12],1,phi,theta,psi,-1*shx/apix,-1*shy/apix))
		fout2.write('%s\t%s\t%s\t%s\t%.6f\t %.6f\t%.6f\t%s\t%s\t%s\t %s\t%s\t%s\t%s\t%s\t %s\t%s\n'%(t2[0],t2[1],t2[2],t2[3],phi_,theta_,psi_,t2[7],t2[8],t2[9],t2[10],t2[11],t2[12],t2[13],t2[14],t2[15],t2[16]))

	f2.close()
	fout2.close()



if __name__ == "__main__":
	params = setupParserOptions()
	mainloop(params)
	fre2relion(params)
