#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
class for Brown Adipose Tissue 

Created on Fri Sep  8 08:26:03 2017

@author: louis
"""

import numpy as np
import os
import string
import subprocess
import SimpleITK as sitk
from scipy import ndimage

class postRefinement(object):
    """ class for fine adjustment and refinement of segmentation """
    def __init__(self, outputdir):
        self.outputdir = outputdir
        if not os.path.exists(self.outputdir):
            subprocess.call('mkdir ' + '-p ' + self.outputdir, shell=True)
    
    def inputImage(self, filepathLabel, filepathFatFraction, filepathT2Star, filepathFat, filepathWater):
        
        # label
        self.filepathLab = filepathLabel
        self.imageLab = sitk.ReadImage(self.filepathLab)
        self.arrayLab = sitk.GetArrayFromImage(self.imageLab)
        self.shapeLab = np.shape(self.arrayLab)

        # fat fraction(FF)
        self.filepathFF = filepathFatFraction
        self.imageFF = sitk.ReadImage(self.filepathFF)
        self.arrayFF = sitk.GetArrayFromImage(self.imageFF)
        self.shapeFF = np.shape(self.arrayFF)
        print "maximum value of FatFraction: %f"%(np.amax(np.amax(np.amax(self.arrayFF))))
        print "minimum value of FatFraction: %f"%(np.amin(np.amin(np.amin(self.arrayFF))))
        
        # T2 star(T2S)
        self.filepathT2S = filepathT2Star
        self.imageT2S = sitk.ReadImage(self.filepathT2S)
        self.arrayT2S = sitk.GetArrayFromImage(self.imageT2S)
        self.shapeT2S = np.shape(self.arrayT2S)
        print "maximum value of T2Star: %f"%(np.amax(np.amax(np.amax(self.arrayT2S))))
        print "minimum value of T2Star: %f"%(np.amin(np.amin(np.amin(self.arrayT2S))))
        
        # Fat(F)
        self.filepathF = filepathFat
        self.imageF = sitk.ReadImage(self.filepathF)
        self.arrayF = sitk.GetArrayFromImage(self.imageF)
        self.shapeF = np.shape(self.arrayF)
        print "maximum value of Fat: %f"%(np.amax(np.amax(np.amax(self.arrayF))))
        print "minimum value of Fat: %f"%(np.amin(np.amin(np.amin(self.arrayF))))
        
        # Water(W)
        self.filepathW = filepathWater
        self.imageW = sitk.ReadImage(self.filepathW)
        self.arrayW = sitk.GetArrayFromImage(self.imageW)
        self.shapeW = np.shape(self.arrayW)
        print "maximum value of Water: %f"%(np.amax(np.amax(np.amax(self.arrayW))))
        print "minimum value of Water: %f"%(np.amin(np.amin(np.amin(self.arrayW))))

    def fineAdjustSegmentation(self, thresholdFF, erosionIterNum, T2Sthreshold, ifLogOutput):

        # Automated fine adjustment of segmentation
        # FFthreshold >= 40%, R2_star <= 50 s-1 (T2_star >= 0.2S)
        idxFF = (self.arrayFF >= thresholdFF)
        self.arrayLab[~idxFF]= np.uint8(0)
       
        ErosionSkinArea = np.zeros_like(self.arrayLab)
        x,y,z = np.shape(self.arrayLab)

        for i in xrange(x):
            ErosionSkinArea[i,:,:]= ndimage.binary_erosion(self.arrayLab[i,:,:],\
            iterations = erosionIterNum).astype(self.arrayLab.dtype)    
        
        idxR2S = (self.arrayT2S >= T2Sthreshold)
        ErosionSkinArea[~idxR2S]= np.uint8(0)

        self.arrayLab = ErosionSkinArea        
        self.imageLab = sitk.GetImageFromArray(self.arrayLab)
        self.imageLab.SetOrigin(self.imageLab.GetOrigin())                               
        self.imageLab.SetSpacing(self.imageLab.GetSpacing())                                
        self.imageLab.SetDirection(self.imageLab.GetDirection())
        
        if ifLogOutput != False:
            dirname = self.outputdir + '/'             
            name, ext = os.path.splitext(self.filepathLab)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[1:-2], "_")
            outputdirection = dirname + outBaseName + '_Lab_fineAdjust' + ext             
            sitk.WriteImage(self.imageLab, outputdirection)
        else:
            return self.arrayLab

    def removeBoneadnAir(self, threshold, ifLogOutput):
        """ remove the bone and the air of Lab map 
        """
        waterAndfatArray = self.arrayW + self.arrayF
        idx = (waterAndfatArray <= threshold) # idx of the bone and the air area. 
        self.arrayLab[idx] = 0
 
        self.imageLab = sitk.GetImageFromArray(self.arrayLab)
        self.imageLab.SetOrigin(self.imageLab.GetOrigin())                               
        self.imageLab.SetSpacing(self.imageLab.GetSpacing())                                
        self.imageLab.SetDirection(self.imageLab.GetDirection())

        if ifLogOutput != False:
            dirname = self.outputdir + '/'             
            name, ext = os.path.splitext(self.filepathLab)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[1:-2], "_")
            outputdirection = dirname + outBaseName + '_Lab_RBandA' + ext             
            sitk.WriteImage(self.imageLab, outputdirection)
        else:
            return self.arrayLab       

    def finalRefine(self, dilationPara, ClossingPara, ErosionPara, ifLogOutput):
        """ remove the skin voxel
        """
        firstdilationArea = np.zeros_like(self.arrayLab)
        secondClossingArea = np.zeros_like(self.arrayLab)
        thirdErosionArea = np.zeros_like(self.arrayLab)
        x,y,z = np.shape(self.arrayLab)

        if dilationPara['IfDilation'] == True:
            DilationSize = dilationPara['structure']
            iternum = dilationPara['iterations']
            for i in xrange(x):
                firstdilationArea[i,:,:]= ndimage.binary_dilation(self.arrayLab[i,:,:],\
                structure = DilationSize, iterations = iternum).astype(self.arrayLab.dtype)
            self.arrayLab = firstdilationArea
        if ClossingPara['IfClossing'] == True:
            ClossingSize = ClossingPara['structure']
            iternum = ClossingPara['iterations']
            for i in xrange(x):
                secondClossingArea[i,:,:] = ndimage.binary_closing(self.arrayLab[i,:,:],\
                structure = ClossingSize, iterations = iternum).astype(self.arrayLab.dtype)
            self.arrayLab = secondClossingArea

        if ErosionPara['IfErosion'] == True:
            ErosionSize = ErosionPara['structure']
            iternum = ErosionPara['iterations']
            for i in xrange(x):
                thirdErosionArea[i,:,:] = ndimage.binary_erosion(self.arrayLab[i,:,:],\
                structure = ErosionSize, iterations = iternum).astype(self.arrayLab.dtype)
            self.arrayLab = thirdErosionArea

        self.imageLab = sitk.GetImageFromArray(self.arrayLab)
        self.imageLab.SetOrigin(self.imageLab.GetOrigin())                               
        self.imageLab.SetSpacing(self.imageLab.GetSpacing())                                
        self.imageLab.SetDirection(self.imageLab.GetDirection())

        if ifLogOutput != False:
            dirname = self.outputdir + '/'             
            name, ext = os.path.splitext(self.filepathLab)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[1:-2], "_")
            outputdirection = dirname + outBaseName + '_Lab_finalRefine' + ext             
            sitk.WriteImage(self.imageLab, outputdirection)
        else:
            return self.arrayLab       