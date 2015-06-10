#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import shutil
import os
import codecs
import sys
import math
import matplotlib.pyplot as plt


DATA_PATH = 'data/data_2.in'
FILTER_PATH = 'filter/FIR_kernel_27.in'


n = 20000
lambda_k = 0.99
c_k = 4
lambda_D = 0.99
lambda_Th = 0.99

filterDelay = 14
filterN = 27


def sign(x):
    if (x < 0):
        return -1
    else:
        return 1

def readECG():
    lines = open(DATA_PATH,'r').readlines()
    ecg = map(lambda x: float(x) , lines) #* 255
    '''
    baseline = 0
    c = 0
    for e in ecg:
        if e != -10000000 :
            baseline += e
            c += 1
    baseline /= c
    ecg = map(lambda x: x if x != -10000000 else baseline, ecg)
    ''' 
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
    '''

    #ecg = map(lambda x: x * 255, ecg)
        
    #for testing
    return ecg[0:n]

def filterBandPassFIR(ecg):
    ker = open(FILTER_PATH,'r').readlines()
    ker = map(str.strip, ker)
    ker = map(float, ker)
    global filterN 
    filterN = len(ker)
    global filterDelay
    filterDelay = filterN / 2

    newEcg = []
    for i in range(len(ker), len(ecg)):
        v = 0.0
        for j in range(0,len(ker)):
            v += ker[j] * ecg[i - j]
        newEcg.append(v)

    return newEcg

def nonLinearFilter(ecg):
    return map(lambda x: sign(x) * x * x, ecg)

def addHFS(input):
    z = []
    k_prev = 0.0
    for i in range(len(input)):
        k = lambda_k * k_prev + (1 - lambda_k) * abs(input[i]) * c_k
        z.append(input[i] + pow(-1, i) * k)
        k_prev = k
    return z

def computeFeature(z):
    d = []
    d_prev = 0
    for i in range(1,len(z)):
        dd = abs(sign(z[i]) - sign(z[i - 1])) / 2
        new_d = lambda_D * d_prev + (1 - lambda_D) * dd
        d.append(new_d) 
        d_prev = new_d
    return d
def computeTheta(d):
    Th = []
    th_prev = 0.0
    for i in range(len(d)):
        new_Th = lambda_Th * th_prev + (1 - lambda_Th) * d[i]
        Th.append(new_Th) 
        th_prev = new_Th
    return Th

def getEvents(D, Th):
    events = [] #each event is a tuple: (start, end)
    start = 0
    inEvent = False
    needToCombine = False
    lastEvent = (-10000, -1000)
    for i in range(len(D)):
        if D[i] < Th[i] and not inEvent :
            #new event detected
            start = i
            inEvent = True
            if (i - lastEvent[1]) < 42 :
                start = events[-1][0]
                needToCombine = True   
            else :
                needToCombine = False 
            continue
            #need to check distance from last event
            
        if D[i] > Th[i] and inEvent :
            #event ended
            lastEvent = (start, i)
            if needToCombine :
                events[-1] = lastEvent
            else :
                events.append(lastEvent)
            needToCombine = False
            inEvent = False
            start = 0
    if inEvent :
        events.append((start, len(D)))
    return events

def getRpeaks(events, y):
    rPeaks = []
    for event in events:
        start = event[0]
        end = event[1]
        y_max = -10000000000
        y_max_ind = 0
        for i in range(start, end) :
            if y[i]  > y_max :
                y_max = y[i]
                y_max_ind = i
        if y_max_ind > 0 :
            rPeaks.append(y_max_ind)
    return rPeaks

def main():
    plotFlag = True

    tm = time.time()

    ecg = readECG()
    
    t = [0.002 * x for x in range(len(ecg))]


    print "read data in ", time.time() - tm
    #need to filter signal
    
    ecgFIR = filterBandPassFIR(ecg)
    y = nonLinearFilter(ecgFIR)


    print "data filtered in ", time.time() - tm
    n = len(ecg)
    m = len(ecgFIR)

    print 'n = ', n
    print 'm = ', m
    print 'filter N = ', filterN
    print 'filterDelay = ', filterDelay
    if plotFlag :
        plt.figure(1)
        plt.plot(t[0:n], ecg[0:n], t[filterDelay:n - (filterN - filterDelay)], ecgFIR[0:m], 'r-')


    #adding high-frequency seq
    z = addHFS(y)
    print "high-frequency seq added in ", time.time() - tm

    D = computeFeature(z)
    print "d computed in ", time.time() - tm

    Th = computeTheta(D)

    events = getEvents(D, Th)

    print "events = ", events
    

    rPeaks = getRpeaks(events, y)

    print "R peaks = ", rPeaks



    
    d = (n - m) / 2 
    print d

    
    if plotFlag :


        plt.figure(2)
        #plt.plot(t[0:n], ecg[0:n],t[d : n - d - 1], z[0:m], 'r-')
        #plt.plot(t[0:m], y[0:m],t[0 : m - 1], D[0:m - 1], 'r-')
        plt.plot(t[0:len(D)],Th[0:len(D)], t[0 : len(D)], D, 'r-')
        #plt.plot(t[0:len(D)],events,'r-')
        
        plt.figure(3)
        #plt.plot(t[:1000], ecg[:1000])
        tt = map(lambda x: x * 0.002, rPeaks);
        d = []
        for item in rPeaks :
            d.append(y[item])
        plt.plot(t[:m], y[:m],'b-', tt, d, 'ro')


    plt.show()
    quit()

    
    
    return

main()  