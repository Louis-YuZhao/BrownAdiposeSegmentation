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
def main():
    # preprocessing
    inputroot = '/media/louis/Volume/ResearchData/BATSegmentationData'
    outputroot = '/media/louis/Volume/ProgramWorkResult/BATSegresult'
    filepathFatFraction = inputroot + '/ImageSplitFF/FileList.txt' 
    filepathT2Star = inputroot + '/ImageSplitT2S/FileList.txt' 
    filepathFat = inputroot + '/ImageSplitF/FileList.txt'
    ilepathWater = inputroot + '/ImageSplitW/FileList.txt'

    listFatFraction = readTxtIntoList(filepathFatFraction)
    listT2Star = readTxtIntoList(filepathT2Star)
    listFat = readTxtIntoList(filepathFat)
    listWater = readTxtIntoList(ilepathWater)

    ifLogOutput = True

    if len(listFatFraction)== len(filepathT2Star)== len(filepathT2Star)== len(filepathT2Star):
        N = len(listFatFraction)
    else:
        raise ValueError('the length of the files should be same')

    outputdir = outputroot + '/WFPreProcessing'
    WFprocessing = BATpreprocessingWF(outputdir)
    for i in xrange(N):
        WFprocessing.inputImage(listFatFraction[i], listT2Star[i], listFat[i], listWater[i])
        WFprocessing.computeWFmap(ifLogOutput)
        WFprocessing.removeBackgroundmap(threshold = 28, ifLogOutput)
        WFprocessing.removeSkinVoxels(FFthreshold = 30, T2Sthreshold = 0.125, iternum = 5,\
                                    ClossingSize = 5, ErosionSize = 7, ifLogOutput)
    # recovering the preprocessed result
    preprocessedFile = outputroot + '/WFPreProcessing'
    ImageList = []
    for root, dirnames, filenames in os.walk(preprocessedFile):
        for name in enumerate(filenames):
            if name.endswith('_WF_RSkin.dcm'):
                ImageDir = os.path.join(root, name)
                ImageList.append(ImageDir)
    WriteListtoFile(ImageList, preprocessedFile + '/preProcessedWFImage.txt')

    # multi atlas registartion
    atlasNum = 10

    readTxtIntoList(preprocessedFile + '/preProcessedWFImage.txt')
    WriteListtoFile(ImageList[: atlasNum], preprocessedFile + '/atlasImageListdir.txt')
    WriteListtoFile(ImageList[atlasNum :], preprocessedFile + '/testImageListdir.txt')

    readTxtIntoList(inputroot + '/ImageSplitLeft/FileList.txt')
    WriteListtoFile(ImageList[: atlasNum], inputroot + '/ImageSplitLeft/atlasLabelListdir.txt')
    WriteListtoFile(ImageList[atlasNum :], preprocessedFile + '/testLabelListdir.txt')

    atlasImageListdir = preprocessedFile + '/atlasImageListdir.txt' 
    atlasLabelListdir = inputroot + '/ImageSplitLeft/atlasLabelListdir.txt' 
    
    segmentation = MultiAtlasSegmentation(atlasImageListdir, atlasLabelListdir, atlasNum)
    segmentation.readAtlasImagetoList()
    segmentation.readAtlasLabeltoList()

    affinePara={}
    affinePara['Optimizer']='AdaptiveStochasticGradientDescent'
    affinePara['NumberOfResolutions']='4'
    affinePara['MaximumNumberOfIterations']='3000'
    affinePara['FinalBSplineInterpolationOrder']='1'

    BsplinePara={}
    BsplinePara['Optimizer']='AdaptiveStochasticGradientDescent'
    BsplinePara['NumberOfResolutions']='4'
    BsplinePara['MaximumNumberOfIterations']='7000'
    BsplinePara['FinalBSplineInterpolationOrder']='1'
    BsplinePara['FinalGridSpacingInPhysicalUnits']='20.0 20.0 20.0'

    fixedimagedir = preprocessedFile + '/testImageListdir.txt'
    fixedimageList = readTxtIntoList(fixedimagedir)
    fusedLabeldir = outputroot + '/fusedLabels/'
    if not os.path.exists(affineImageOut):
        subprocess.call('mkdir ' + '-p ' + affineImageOut, shell=True)
    threshold = 0.5

    for fixIm in fixedimageList:
        name, ext = os.path.splitext(fixIm)
        baseName = os.path.basename(name)
        
        # affine segmentation
        affineImageOut = outputroot + '/MASAffine/alignedImage/' + baseName
        affineLabelOut = outputroot + '/MASAffine/alignedLabel/' + baseName
        if not os.path.exists(affineImageOut):
            subprocess.call('mkdir ' + '-p ' + affineImageOut, shell=True)
        if not os.path.exists(affineLabelOut):
            subprocess.call('mkdir ' + '-p ' + affineLabelOut, shell=True)
        
        segmentation.AffineRegistartion(fixIm, affineImageOut, affineLabelOut, affinePara)
        
        # Bspline segmentation
        BsplineimageOut = outputroot + '/MASBspline/alignedImage/' + baseName
        BsplinelabelOut = outputroot + '/MASBspline/alignedLabel/' + baseName 
        if not os.path.exists(BsplineimageOut):
            subprocess.call('mkdir ' + '-p ' + BsplineimageOut, shell=True)
        if not os.path.exists(BsplinelabelOut):
            subprocess.call('mkdir ' + '-p ' + BsplinelabelOut, shell=True)

        segmentation.BsplineRegistartion(fixIm, BsplineimageOut, BsplinelabelOut, BsplinePara)

        # Label Fussion
        labelListdir = outputroot + '/MASBspline/alignedLabel/' + baseName + 'FileList.txt'
        detailfusedLabeldir = fusedLabeldir + baseName +'_fusedLabel.nrrd'
        segmentation.FusionofSegmentation(labelListdir,threshold, detailfusedLabeldir)

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

if __name__ == "__main__":
    main()
