#!/usr/bin/env python
# -*- coding: utf-8 -*-
import time
import shutil
import os
import codecs
import sys
import math
import matplotlib.pyplot as plt

DATA_PATH = 'data/data_1.in'

n = 10000

filterDelay = 14
filterN = 27



def sign(x):
    if (x < 0):
        return -1
    else:
        return 1

def readECG():
    lines = open(DATA_PATH,'r').readlines()
    
    #lines = map(lambda x: str.split(str(x))[2], lines)

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
    '''
    
    #ecg = map(lambda x: x * 255, ecg)
        
    #for testing
    return ecg[0:n]

def filterBandPassFIR(ecg, filterPATH):
    ker = open(filterPATH,'r').readlines()
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


def getFeature(w, indx):
    feature = []
    p_prev = 0.0
    for ind in indx :
        p_prev += abs(w[ind][0])
    for i in range(1, len(w[1])):
        # if (i % 32 == 0) :
        p = 0.0
        for ind in indx :
            p += abs(w[ind][i])
        c = p
        p = (p_prev + p) / 2.0
        p_prev = c
        feature.append(p)

    return feature


def getMean(x,indx):
    mean = 0
    for ind in indx :
        mean += x[ind]
    mean /= len(indx)
    return mean



def processLevel1(p):
    peaks = []
    for i in range(1, len(p) - 1) :
        if p[i - 1] <= p[i] and p[i] > p[i + 1] :
            #adding index of inflexion point
            peaks.append(i) 
    return peaks





def main():
    plotFlag = True

    tm = time.time()

    ecg = readECG()
    
    t = [0.002 * x for x in range(len(ecg))]


    # print "read data in ", time.time() - tm
    #need to filter signal
    
    ecgFIR_1 = filterBandPassFIR(ecg, 'filter/FP_FIR_1.in')
    ecgFIR_2 = filterBandPassFIR(ecg, 'filter/FP_FIR_2.in')
    ecgFIR_3 = filterBandPassFIR(ecg, 'filter/FP_FIR_3.in')
    ecgFIR_4 = filterBandPassFIR(ecg, 'filter/FP_FIR_4.in')
    ecgFIR = [[], ecgFIR_1, ecgFIR_2, ecgFIR_3, ecgFIR_4]

    # print "data filtered in ", time.time() - tm
    n = len(ecg)
    m = len(ecgFIR_1)

    # print 'n = ', n
    # print 'm = ', m
    # print 'filter N = ', filterN
    # print 'filterDelay = ', filterDelay
    if plotFlag :
        plt.figure(1)
        plt.plot(t[0:n], ecg[0:n], t[filterDelay:n - (filterN - filterDelay)], ecgFIR_4[0:m], 'r-')

    #counting features
    P_1 = getFeature(ecgFIR, [1, 2, 3])
    P_2 = getFeature(ecgFIR, [1, 2, 3, 4])
    P_3 = getFeature(ecgFIR, [2, 3, 4])

    # print "len P_1 == ", len(P_1)

    # print "P_1 = ", P_1
    if plotFlag :
        plt.figure(2)
        plt.plot(t[0:len(P_3)], P_3, 'r-')

    # processing peaks
    sPeaksChan_1 = [] #here we store indexes
    nPeaksChan_1 = []
    sPeaksChan_2 = []
    nPeaksChan_2 = []
    sPeaks = []
    nPeaks = []
    SLChan_1 = 0
    NLChan_1 = 0
    SLChan_2 = 0
    NLChan_2 = 0
    SL = 0
    NL = 0

    T_1 = 0.08
    T_2 = 0.7
    T_4 = 0.3

    peaks = processLevel1(P_1)

    # print 'events = ', peaks

    d_1 = []
    d_2 = []
    d = []
    rPeaks = []

    for peak in peaks :
        #level 2 starts
        DChan_1 = 0.0
        
        if (P_2[peak] > SLChan_1) : 
            DChan_1 = 1.0
        else :
            if (P_2[peak] < NLChan_1) :
                DChan_1 = 0.0
            else :
                DChan_1 = (P_2[peak] - NLChan_1) / (SLChan_1 - NLChan_1)
        isBeat_1 = False
        
        if (DChan_1 > T_1) :
            #peak is true for channel 1
            isBeat_1 = True
            sPeaksChan_1.append(peak)
            SLChan_1 = getMean(P_2, sPeaksChan_1)    
        else :
            isBeat_1 = False
            nPeaksChan_1.append(peak)
            NLChan_1 = getMean(P_2, nPeaksChan_1)
        #channel 1 processed
        # print "DChan_1 = ", DChan_1, "isBeat = ",isBeat_1


        DChan_2 = 0.0
        # print "peak = ", peak, "; P_2 = ", P_2[peak], "; SLChan_2 = ", SLChan_2, "; NLChan_2 = ", NLChan_2
        if (P_2[peak] > SLChan_2) : 
            DChan_2 = 1.0
        else :
            if (P_2[peak] < NLChan_2) :
                DChan_2 = 0.0
            else :
                DChan_2 = (P_2[peak] - NLChan_2) / (SLChan_2 - NLChan_2)
        isBeat_2 = False
        if (DChan_2 > T_2) :
            #peak is true for channel 2
            isBeat_2 = True
            sPeaksChan_2.append(peak)
            SLChan_2 = getMean(P_2, sPeaksChan_2)    
        else :
            isBeat_2 = False
            nPeaksChan_2.append(peak)
            NLChan_2 = getMean(P_2, nPeaksChan_2)
        #channel 2 processed
        #level 2 done
        # print "DChan_2 = ", DChan_2, "isBeat = ",isBeat_2

        #test
        d_1.append(DChan_1)
        d_2.append(DChan_2)


        #level 3 starts
        isBeat = False
        if isBeat_2 :
            isBeat = True
        else :
            if isBeat_1 :
                #need to comare deltas
                delta_1 = (DChan_1 - T_1) / (1 - T_1)
                delta_2 = (T_2 - DChan_2) / T_2
                isBeat = delta_1 > delta_2
        #level 3 done

        #level 4 starts
        isActualRPeak = False
        if isBeat :
            isActualRPeak = True
            sPeaks.append(peak)
            SL = getMean(P_3, sPeaks)
        else :
            D = 0.0
            if (P_3[peak] > SL) : 
                D = 1.0
            else :
                if (P_3[peak] < NL) :
                    D = 0.0
                else :
                    D = (P_3[peak] - NL) / (SL - NL)
            if (D > T_4) :
                #peak is true for channel 2
                sPeaks.append(peak)
                isActualRPeak = True
                SL = getMean(P_3, sPeaks)    
            else :
                nPeaks.append(peak)
                NL = getMean(P_3, nPeaks)
            d.append(D)
        #level 4 done
        if isActualRPeak :
            realPeakInd = peak + filterDelay + 1
            if len(rPeaks) == 0:
                rPeaks.append(realPeakInd)
                continue
            if abs(realPeakInd  - rPeaks[-1]) < 125:
                if ecg[realPeakInd] > ecg[rPeaks[-1]] :
                    rPeaks[-1] = realPeakInd
            else :
                rPeaks.append(realPeakInd)


    # if plotFlag :
    #     plt.figure(3)
    #     plt.plot(t[0:len(d_1)], d_1,t[0:len(d_2)], d_2,t[0:len(d)], d, 'r-')

     

    # print "signal peaks channel 1 = ", sPeaksChan_1
    # print "SL channel 1 = ", SLChan_1
    # print "NL channel 1 = ", NLChan_1

    # print "###################################"
    # print "signal peaks channel 2 = ", sPeaksChan_2
    # print "SL channel 2 = ", SLChan_2
    # print "NL channel 2 = ", NLChan_2



    # print "###################################"
    # print "peaks = ", sPeaks   
    # print "###################################"
    print "R peaks = ", rPeaks



    
    d = (n - m) / 2 
    print d

    
    if plotFlag :

        plt.figure(2)
        #plt.plot(t[:1000], ecg[:1000])
        tt = map(lambda x: x * 0.002, rPeaks);
        d = []
        for item in rPeaks :
            d.append(ecg[item])
        plt.plot(t[:n], ecg[:n],'b-', tt, d, 'ro')


    plt.show()
    quit()

    
    
    return

main()  