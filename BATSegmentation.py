#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
main function for Brown Adipose Tissue segmentation

Created on Fri Sep  8 08:26:03 2017

@author: louis
"""

import numpy as np
import dicom
import os
import string
from scipy import ndimage

def readTxtIntoList(filename):
   flist = []
   with open(filename) as f:
         flist = f.read().splitlines()
   return flist

def WriteListtoFile(filelist, filename):
    with open(filename, 'w') as f:
        for i in filelist:
            f.write(i+'\n')
    return 1

# preprocessing
filepathFatFraction = '', 
filepathT2Star = '', 
filepathFat = '', 
ilepathWater = '', 

listFatFraction = readTxtIntoList(filepathFatFraction)
listT2Star = readTxtIntoList(filepathT2Star)
listFat = readTxtIntoList(filepathFat)
listWater = readTxtIntoList(ilepathWater)


ifLogOutput = True

if len(listFatFraction)== len(filepathT2Star)== len(filepathT2Star)== len(filepathT2Star):
    N = len(listFatFraction)
else:
    raise ValueError('the length of the files should be same')
for i in xrange(N):
    WFprocessing = BATpreprocessingWF()
    WFprocessing.inputImage(listFatFraction[i], listT2Star[i], listFat[i], listWater[i])
    WFprocessing.computeWFmap(ifLogOutput)
    WFprocessing.removeBackgroundmap(threshold = 28, ifLogOutput)
    WFprocessing.removeSkinVoxels(FFthreshold = 30, T2Sthreshold = 0.125, iternum = 5,\
                                  ClossingSize = 5, ErosionSize = 7, ifLogOutput)

# multi atlas registartion
atlasImageListdir = '' 
atlasLabelListdir = '' 
atlasNum = 10
segmentation = MultiAtlasSegmentation(atlasImageListdir, atlasLabelListdir, atlasNum)
segmentation.readAtlasImagetoList()
segmentation.readAtlasLabeltoList()

affinePara={}
movingimagedir = ''
imageOutputDir = ''
labelOutputDir = ''
movingimageList = readTxtIntoList(movingimagedir)
for movIm in movingimageList:
    segmentation.AffineRegistartion(movIm, imageOutputDir, labelOutputDir, affinePara)
    
BsplinePara={}
movingimagedir = ''
imageOutputDir = ''
labelOutputDir = '' 
movingimageList = readTxtIntoList(movingimagedir)
for movIm in movingimageList:
    segmentation.BsplineRegistartion(self, movIm, imageOutputDir, labelOutputDir, BsplinePara)

for i in xrange(len(movingimageList)):
    labelListdir = ''
    fusedLabeldir = ''
    threshold = 0.5
    segmentation.FusionofSegmentation(labelListdir,threshold,fusedLabeldir)

# PostProcessing

filepathLabel = '' 
filepathFatFraction = '' 
filepathT2Star = ''
filepathFat = ''
filepathWater = ''

listLabel = readTxtIntoList(filepathLabel)
listFatFraction = readTxtIntoList(filepathFatFraction)
listT2Star = readTxtIntoList(filepathT2Star)
listFat = readTxtIntoList(filepathFat)
listWater = readTxtIntoList(ilepathWater)

ifLogOutput = True

if len(listLabel)==len(listFatFraction)== len(filepathT2Star)== len(filepathT2Star)== len(filepathT2Star):
    N = len(listFatFraction)
else:
    raise ValueError('the length of the files should be same')
dilationPara = {}
ClossingPara = {}
ErosionPara = {}
for i in xrange(len(N)):
    postProcessing = postRefinement()
    postProcessing.inputImage(listLabel[i], listFatFraction[i], listT2Star[i], listFat[i], listWater[i])
    postProcessing.fineAdjustSegmentation(thresholdFF =, erosionIterNum =, T2Sthreshold =, ifLogOutput)
    postProcessing.removeBoneadnAir(threshold =, ifLogOutput)
    postProcessing.finalRefine(dilationPara, ClossingPara, ErosionPara, ifLogOutput)

# calculate dice score
