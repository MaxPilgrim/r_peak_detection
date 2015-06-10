#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import shutil
import os
import codecs
import sys
import matplotlib.pyplot as plt

DATA_PATH = 'data/data_1.in'
FILTER_PATH = 'filter/FIR_kernel.in'

n = 6000

def readECG():
	lines = open(DATA_PATH,'r').readlines()
	ecg = map(lambda x: float(x) , lines) #* 255
	baseline = 0
	c = 0
	for e in ecg:
		if e != -10000000 :
			baseline += e
			c += 1
	baseline /= c
	ecg = map(lambda x: x if x != -10000000 else baseline, ecg)
	'''	
	a = min(ecg)
	b = max(ecg)
	d = 127.0
	c = -127.0
	
	ecg = map(lambda x: (x - a) / (b - a) * (d - c) + c, ecg) 

	f = open('data.out','w')
	for e in ecg:
		f.write(str(e) + "\n")
	f.close()

	ecg = map(lambda x: x * 255, ecg)
	'''	
	#for testing
	return ecg[0:n]

def filterBandPassFIR(ecg):
	ker = open(FILTER_PATH,'r').readlines()
	ker = map(str.strip, ker)
	ker = map(float, ker)
	
	newEcg = []
	for i in range(len(ker), len(ecg)):
		v = 0.0
		for j in range(0,len(ker)):
			v += ker[j] * ecg[i - j]
		newEcg.append(v)

	return newEcg

def getQRS(ecg):
	d1 = []
	#count squared double differences
	for i in range(len(ecg) - 1):
		d1.append(ecg[i + 1] - ecg[i])
		
	d = []
	for j in range(len(ecg) - 2):
		x = d1[j + 1] - d1[j]
		d.append(x * x)
	dMax = max(d)
	#print "dmax = ", dMax
	#print "before filter d.len = ", len(d)
	d = zip(d, range(len(d)))
	d = filter(lambda x: x[0] > 0.09 * dMax, d)
	#print "after filter d.len = ", len(d)
	d = sorted(d, key = lambda x: x[0], reverse = True)

	qrs = [d[0]]
	for peak in d:
		#check if current peak is in qrs region [-75ms; 75ms]
		flag = True
		for qrsPeak in qrs:
			if abs(qrsPeak[1]-peak[1]) * 2 < 75 :
				flag = False
				break;
		if flag :
			qrs.append(peak)

	#qrs = sorted(qrs, key = lambda x: x[1])

	# print "qrs regions  = ", qrs
	if qrs[0][1] < 10 :
		return qrs[1:]
	else :
		return qrs

def getRPeaks(ecg, qrs) :
	rPeaks = []
	for qrsReg in qrs :
		start = max(0, qrsReg[1] - 38)
		end = min(len(ecg), qrsReg[1] + 38)
		ecg_min = 1000000000
		ecg_max = -1000000000	
		#count mean 
		for i in range(start, end) :
			ecg_min = min(ecg_min, ecg[i])
			ecg_max = max(ecg_max, ecg[i])
		ecg_mean = (ecg_max + ecg_min) / 2.0
		#find local maximum of relative magnitudes
		ecg_max = -100000000
		ecg_max_ind = 0
		for i in range(start, end) :
			if ecg[i] - ecg_mean > ecg_max :
				ecg_max = ecg[i] - ecg_mean
				ecg_max_ind = i
		if ecg_max_ind > 0 :
			rPeaks.append(ecg_max_ind)



	#print rPeaks
	return rPeaks



def filterRPeaks(ecg, rPeaks):
	#Any peaks detected within 200 ms of the first is considered as noise peak and eliminated

	newRPeaks = [rPeaks[0]]
	for i in range(1,len(rPeaks)):
		flag = True
		for peak in newRPeaks :
			if abs(rPeaks[i] - peak) < 100 :
				flag = False
				break
		if flag :
			newRPeaks.append(rPeaks[i])
	rPeaks = sorted(newRPeaks)
	return rPeaks

def main():
	plotFlag = True

	tm = time.time()

	ecg = readECG()
	
	t = [0.002 * x for x in range(len(ecg))]
	#print ecg 
	#before filter
	#if plotFlag :
	#	plt.figure(1)
	#	plt.plot(t[:n], ecg[:n])


	# print "read data in ", time.time() - tm
	#need to filter signal
	ecgFIR = filterBandPassFIR(ecg)
	# print "data filtered in ", time.time() - tm
	n = len(ecg)
	m = len(ecgFIR)

	if plotFlag :
		plt.figure(2)
		plt.plot(t[60:n], ecg[60:n],t[30:m+30], ecgFIR[:m], 'r-')


	#plt.show()

	#quit()

	m = len(ecgFIR)
	
	qrs = getQRS(ecgFIR)
	# print "found qrs in ", time.time() - tm
	
	rPeaks = getRPeaks(ecg, qrs) 	
	
	# print "found R-peaks in ", time.time() - tm
	
	
	rPeaks = filterRPeaks(ecg, rPeaks)

	# print "fitered R-peaks in ", time.time() - tm
	
	print "R-peaks = ", rPeaks
	
	if plotFlag :
		plt.figure(3)
		#plt.plot(t[:1000], ecg[:1000])
		tt = map(lambda x: x * 0.002, rPeaks);
		d = []
		for item in rPeaks :
			d.append(ecg[item])
		plt.plot(t[:m], ecg[:m],'b-', tt, d, 'ro')
	plt.show()
	
	return

main()