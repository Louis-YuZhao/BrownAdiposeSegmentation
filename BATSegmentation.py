#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
main function for Brown Adipose Tissue segmentation

Created on Fri Sep  8 08:26:03 2017

@author: louis
"""

import numpy as np
import os
import string
import subprocess
import preprocessingV2 as pp
import MultiAtlasSegmentation as MAS
import postRefinement as pr

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
    '''
    whole pipline for multi atlas segmentation
    '''    
    # preprocessing
    inputroot = '/media/louis/Volume/ResearchData/BATSegmentationData'
    outputroot = '/media/louis/Volume/ProgramWorkResult/BATSegresult'
    filepathFatFraction = inputroot + '/ImageSplitFF/FileList.txt' 
    filepathT2Star = inputroot + '/ImageSplitT2S/FileList.txt' 
    filepathFat = inputroot + '/ImageSplitF/FileList.txt'
    ilepathWater = inputroot + '/ImageSplitW/FileList.txt'
    
    listFF = readTxtIntoList(filepathFatFraction)
    listT2Star = readTxtIntoList(filepathT2Star)
    listFat = readTxtIntoList(filepathFat)
    listWater = readTxtIntoList(ilepathWater)
    
    IFLogOutput = True
    
    if len(listFF) == len(listT2Star) == len(listFat) == len(listWater):
        N = len(listFF)
    else:
        raise ValueError('the length of the files should be same')
    
    # PreProcessing WF image
    outputdir = outputroot + '/WFPreProcessing'
    WFprocessing = pp.BATpreprocessingWF(outputdir)
    for i in xrange(N):
        WFprocessing.inputImage(listFF[i], listT2Star[i], listFat[i], listWater[i])
        WFprocessing.computeWFmap(ifLogOutput = IFLogOutput)
        WFprocessing.removeBackgroundmap(threshold = 300, ifLogOutput = IFLogOutput)
        WFprocessing.removeSkinVoxels(FFthreshold = 60, T2Sthreshold = 70, iternum = 2,\
                                      filterSize = (5,5), ClossingSize = np.ones((5,5)),\
                                      ErosionSize = np.ones((3,3)), ifLogOutput = IFLogOutput)
        print "processing num: %d" %(i)
    # recovering the preprocessed result
    ImageList = []
    for root, dirnames, filenames in os.walk(outputdir):
        for name in filenames:
            if name.endswith('_WF_RSkin.nrrd'):
                ImageDir = os.path.join(root, name)
                ImageList.append(ImageDir)
    ImageList.sort()
    WriteListtoFile(ImageList, outputdir + '/preProcessedWFImage.txt')
    
    # PreProcessing T2_Star image
    outputdir = outputroot + '/T2SPreProcessing'
    T2Sprocessing = pp.BATpreprocessingT2S(outputdir)
    for i in xrange(N):
        T2Sprocessing.inputImage(listFF[i], listT2Star[i], listFat[i], listWater[i])
        T2Sprocessing.removeBackgroundmap(threshold = 300, ifLogOutput = IFLogOutput)
        T2Sprocessing.removeSkinVoxels(FFthreshold = 60, T2Sthreshold = 70, iternum = 2,\
                                       filterSize = (5,5), ClossingSize = np.ones((5,5)),\
                                       ErosionSize = np.ones((3,3)), ifLogOutput = IFLogOutput)
        T2Sprocessing.reduceNoise(filterSize = (5,5), ifLogOutput = IFLogOutput)
        print "processing num: %d" %(i)
    # recovering the preprocessed result
    ImageList = []
    for root, dirnames, filenames in os.walk(outputdir):
        for name in filenames:
            if name.endswith('_T2S_RSkin.nrrd'):
                ImageDir = os.path.join(root, name)
                ImageList.append(ImageDir)
    ImageList.sort()
    WriteListtoFile(ImageList, outputdir + '/preProcessedR2SImage.txt')
    
    #%%
    # multi atlas registartion
    atlasNum = 10
    
    preprocessedFile = '/media/louis/Volume/ProgramWorkResult/BATSegresult/WFPreProcessing'
    ImageList = readTxtIntoList(preprocessedFile + '/preProcessedWFImage.txt')
    WriteListtoFile(ImageList[: atlasNum], preprocessedFile + '/atlasImageListdir.txt')
    WriteListtoFile(ImageList[atlasNum :], preprocessedFile + '/testImageListdir.txt')
    
    ImageList = readTxtIntoList(inputroot + '/ImageSplitLeft/FileList.txt')
    WriteListtoFile(ImageList[: atlasNum], inputroot + '/ImageSplitLeft/atlasLabelListdir.txt')
    WriteListtoFile(ImageList[atlasNum :], preprocessedFile + '/testLabelListdir.txt')
    
    atlasImageListdir = preprocessedFile + '/atlasImageListdir.txt' 
    atlasLabelListdir = inputroot + '/ImageSplitLeft/atlasLabelListdir.txt' 
    
    segmentation = MAS.MultiAtlasSegmentation(atlasImageListdir, atlasLabelListdir, atlasNum)
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
    fusedLabeldir = outputroot + '/fusedLabels'
    if not os.path.exists(fusedLabeldir):
        subprocess.call('mkdir ' + '-p ' + fusedLabeldir, shell=True)
    
    Fusionthreshold = 0.5
    fusedLabelList = []
    
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
        labelListdir = outputroot + '/MASBspline/alignedLabel/' + baseName + '/FileList.txt'
        detailfusedLabeldir = fusedLabeldir + '/fusedLabel_' + baseName +'.nrrd'
        fusedLabelList.append(detailfusedLabeldir)
        segmentation.FusionofSegmentation(labelListdir, Fusionthreshold, detailfusedLabeldir)
        
    WriteListtoFile(fusedLabelList, fusedLabeldir + '/FileList.txt')
    
    #%%
    # PostProcessing
    filepathLabel = outputroot + '/fusedLabels/FileList.txt' 
    refListLabel = readTxtIntoList(filepathLabel)
    
    refListFF = []
    refListT2Star = []
    refListFat = []
    refListWater = []
    
    ext = '.dcm'
    for i in refListLabel:    
        name, exttemp = os.path.splitext(i)
        BaseName = os.path.basename(name)
        imageBaseName = string.join(BaseName.split("_")[1:-2], "_")
    
        refListFF.append(inputroot + '/ImageSplitFF/' + imageBaseName +'_FF' + ext)
        refListT2Star.append(inputroot + '/ImageSplitT2S/' + imageBaseName +'_T2_STAR'+ ext)
        refListFat.append(inputroot + '/ImageSplitF/' + imageBaseName +'_F'+ ext)
        refListWater.append(inputroot + '/ImageSplitW/' + imageBaseName +'_W'+ ext)
    
    IFLogOutput = True
    
    if len(refListLabel)==len(refListFF)==len(refListT2Star)==len(refListFat)==len(refListWater):
        N = len(refListFF)
    else:
        raise ValueError('the length of the files should be same')
    
    dilationPara = {}
    dilationPara['IfDilation'] = True
    dilationPara['structure'] = np.ones((3,3))
    dilationPara['iterations'] = 1
    ClossingPara = {}
    ClossingPara['IfClossing'] = True
    ClossingPara['structure'] = np.ones((3,3))
    ClossingPara['iterations'] = 1
    ErosionPara = {}
    ErosionPara['IfErosion'] = True
    ErosionPara['structure'] = np.ones((3,3))
    ErosionPara['iterations'] = 1
    
    refinedLabeldir = '/media/louis/Volume/ProgramWorkResult/BATSegresult/refinedLabels'
    if not os.path.exists(refinedLabeldir):
        subprocess.call('mkdir ' + '-p ' + refinedLabeldir, shell=True)
    
    for i in xrange(N):
        postProcessing = pr.postRefinement(refinedLabeldir)
        postProcessing.inputImage(refListLabel[i], refListFF[i], refListT2Star[i],\
                                  refListFat[i], refListWater[i])
        postProcessing.fineAdjustSegmentation(thresholdFF = 60, erosionIterNum = 1,\
                                              T2Sthreshold = 70, ifLogOutput = IFLogOutput)
        postProcessing.removeBoneadnAir(threshold = 300, ifLogOutput = IFLogOutput)
        postProcessing.finalRefine(dilationPara, ClossingPara, ErosionPara,\
                                   ifLogOutput = IFLogOutput)
        print "processing num: %d" %(i)

# calculate dice score
if __name__ == "__main__":
    main()