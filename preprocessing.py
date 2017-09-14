#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
class for Brown Adipose Tissue preprocessing

Created on Fri Sep  8 08:26:03 2017

@author: louis
"""

import numpy as np
import dicom
import os
import string
import subprocess
from scipy import ndimage


class BATpreprocessingWF(object):
    """ class for Brown Adipose Tissue preprocessing """
    def __init__(self, outputdir):
        self.outputdir = outputdir
        if not os.path.exists(self.outputdir):
            subprocess.call('mkdir ' + '-p ' + self.outputdir, shell=True)

    def inputImage(self, filepathFatFraction, filepathT2Star, filepathFat, filepathWater):
        
        # Fat fraction(FF)
        self.filepathFF = filepathFatFraction
        self.imageFF = dicom.read_file(self.filepathFF)
        self.arrayFF = self.imageFF.pixel_array
        self.shapeFF = np.shape(self.arrayFF)
        print "maximum value of FatFraction: %f"%(np.amax(np.amax(np.amax(self.arrayFF))))
        print "minimum value of FatFraction: %f"%(np.amin(np.amin(np.amin(self.arrayFF))))
        
        # T2 star(T2S)
        self.filepathT2S = filepathT2Star
        self.imageT2S = dicom.read_file(self.filepathT2S)
        self.arrayT2S = self.imageT2S.pixel_array
        self.shapeT2S = np.shape(self.arrayT2S)
        print "maximum value of T2Star: %f"%(np.amax(np.amax(np.amax(self.arrayT2S))))
        print "minimum value of T2Star: %f"%(np.amin(np.amin(np.amin(self.arrayT2S))))
        
        # Fat(F)
        self.filepathF = filepathFat
        self.imageF = dicom.read_file(self.filepathF)
        self.arrayF = self.imageF.pixel_array
        self.shapeF = np.shape(self.arrayF)
        print "maximum value of Fat: %f"%(np.amax(np.amax(np.amax(self.arrayF))))
        print "minimum value of Fat: %f"%(np.amin(np.amin(np.amin(self.arrayF))))
        
        # Water(W)
        self.filepathW = filepathWater
        self.imageW = dicom.read_file(self.filepathW)
        self.arrayW = self.imageW.pixel_array
        self.shapeW = np.shape(self.arrayW)
        print "maximum value of Water: %f"%(np.amax(np.amax(np.amax(self.arrayW))))
        print "minimum value of Water: %f"%(np.amin(np.amin(np.amin(self.arrayW))))        

    def computeWFmap(self, ifLogOutput):
        """ canculate the WF map """
        maxitemFF = (np.amax(np.amax(np.amax(self.arrayFF))))
        self.arrayWF = maxitemFF - self.arrayFF
        self.imageWF = dicom.read_file(self.filepathFF)
        self.imageWF.pixel_array = self.arrayWF
        self.imageWF.Pixelimage = self.imageWF.pixel_array.tostring()
        self.shapeWF = np.shape(self.arrayWF)
                
        if ifLogOutput != False:
            dirname = self.outputdir             
            name, ext = os.path.splitext(self.filepathFF)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
            outputdir = dirname + outBaseName + '_WF' + ext             
            self.imageWF.save_as(outputdir)
        else:
            return self.imageWF.pixel_array
    
    def removeBackgroundmap(self, threshold, ifLogOutput):
        """ remove background of WF map 
            remove the bone and the air
        """
        waterAndfatArray = self.arrayW + self.arrayF
        idx = (waterAndfatArray <= threshold) # idx of the bone and the air area. 
        self.arrayWF[idx] = 0
 
        self.imageWF.pixel_array = self.arrayWF
        self.imageWF.Pixelimage = self.imageWF.pixel_array.tostring()

        if ifLogOutput != False:
            dirname = self.outputdir             
            name, ext = os.path.splitext(self.filepathFF)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
            outputdir = dirname + outBaseName + '_WF_RBG' + ext             
            self.imageWF.save_as(outputdir)
        else:
            return self.imageWF.pixel_array       
    
    def removeSkinVoxels(self, FFthreshold, T2Sthreshold, iternum, ClossingSize, ErosionSize, ifLogOutput):
        """ remove the skin voxel
        """
        # skin FFthreshold <= 30%, R2_star >= 80 s-1 (T2_star <= 0.0125S)
        idxFF = (self.arrayFF <= FFthreshold)
        idxT2S = (self.arrayT2S <= T2Sthreshold)
        suspectSkinArea = idxFF & idxT2S
        firstErosionSkinArea = np.zeros_like(suspectSkinArea)
        secondClossingArea = np.zeros_like(suspectSkinArea)
        thirdErosionArea = np.zeros_like(suspectSkinArea)
        x,y,z = self.shapeFF

        for i in xrange(x):
            firstErosionSkinArea[i,:,:]= ndimage.binary_erosion(suspectSkinArea[i,:,:],\
            iterations = iternum).astype(suspectSkinArea.dtype)
        
        for i in xrange(x):
            secondClossingArea[i,:,:] = ndimage.binary_closing(firstErosionSkinArea[i,:,:],\
            structure = ClossingSize).astype(suspectSkinArea.dtype)
        
        for i in xrange(x):
            thirdErosionArea[i,:,:] = ndimage.binary_erosion(secondClossingArea[i,:,:],\
            structure = ErosionSize).astype(suspectSkinArea.dtype)

        afterRemoveSkin = np.zeros(np.shape(self.arrayWF))
        afterRemoveSkin[~suspectSkinArea] = self.arrayWF[~suspectSkinArea]
        idxRemoveSkin = (thirdErosionArea>0)
        afterRemoveSkin[idxRemoveSkin] = self.arrayWF[~idxRemoveSkin]

        self.arrayWF = afterRemoveSkin
        self.imageWF.pixel_array = self.arrayWF
        self.imageWF.Pixelimage = self.imageWF.pixel_array.tostring()


        if ifLogOutput != False:
            dirname = self.outputdir             
            name, ext = os.path.splitext(self.filepathFF)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
            outputdir = dirname + outBaseName + '_WF_RSkin' + ext             
            self.imageWF.save_as(outputdir)
        else:
            return self.imageWF.pixel_array  


class BATpreprocessingT2S(object):
    """ class for Brown Adipose Tissue """
    def __init__(self):
        pass
    
    def inputImage(self, filepathFatFraction, filepathT2Star, filepathFat, filepathWater):
        
        # fat fraction(FF)
        self.filepathFF = filepathFatFraction
        self.imageFF = dicom.read_file(self.filepathFF)
        self.arrayFF = self.imageFF.pixel_array
        self.shapeFF = np.shape(self.arrayFF)
        print "maximum value of FatFraction: %f"%(np.amax(np.amax(np.amax(self.arrayFF))))
        print "minimum value of FatFraction: %f"%(np.amin(np.amin(np.amin(self.arrayFF))))
        
        # T2 star(T2S)
        self.filepathT2S = filepathT2Star
        self.imageT2S = dicom.read_file(self.filepathT2S)
        self.arrayT2S = self.imageT2S.pixel_array
        self.shapeT2S = np.shape(self.arrayT2S)
        print "maximum value of T2Star: %f"%(np.amax(np.amax(np.amax(self.arrayT2S))))
        print "minimum value of T2Star: %f"%(np.amin(np.amin(np.amin(self.arrayT2S))))
        
        # Fat(F)
        self.filepathF = filepathFat
        self.imageF = dicom.read_file(self.filepathF)
        self.arrayF = self.imageF.pixel_array
        self.shapeF = np.shape(self.arrayF)
        print "maximum value of Fat: %f"%(np.amax(np.amax(np.amax(self.arrayF))))
        print "minimum value of Fat: %f"%(np.amin(np.amin(np.amin(self.arrayF))))
        
        # Water(W)
        self.filepathW = filepathWater
        self.imageW = dicom.read_file(self.filepathW)
        self.arrayW = self.imageW.pixel_array
        self.shapeW = np.shape(self.arrayW)
        print "maximum value of Water: %f"%(np.amax(np.amax(np.amax(self.arrayW))))
        print "minimum value of Water: %f"%(np.amin(np.amin(np.amin(self.arrayW))))   

        self.imageT2SProcessed = dicom.read_file(self.filepathT2S)
        self.arrayT2SProcessed = self.imageT2SProcessed.pixel_array  
    
    def removeBackgroundmap(self, threshold, ifLogOutput):
        """ remove background of WF map 
            remove the bone and the air
        """
        waterAndfatArray = self.arrayW + self.arrayF
        idx = (waterAndfatArray <= threshold) # idx of the bone and the air area. 
        self.arrayT2SProcessed[idx] = 0
 
        self.imageT2SProcessed.pixel_array = self.arrayT2SProcessed
        self.imageT2SProcessed.Pixelimage = self.imageT2SProcessed.pixel_array.tostring()

        if ifLogOutput != False:
            dirname = os.path.dirname(self.filepathFF)             
            name, ext = os.path.splitext(self.filepathFF)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
            outputdir = dirname + outBaseName + '_T2S_RBG' + ext             
            self.imageT2SProcessed.save_as(outputdir)
        else:
            return self.imageT2SProcessed.pixel_array       
    
    def removeSkinVoxels(self, FFthreshold, T2Sthreshold, iternum, ClossingSize, ErosionSize, ifLogOutput):
        """ remove the skin voxel
        """
        # skin FFthreshold <= 30%, R2_star >= 80 s-1 (T2_star <= 0.0125S)
        idxFF = (self.arrayFF <= FFthreshold)
        idxT2S = (self.arrayT2S <= T2Sthreshold)
        suspectSkinArea = idxFF & idxT2S
        firstErosionSkinArea = np.zeros_like(suspectSkinArea)
        secondClossingArea = np.zeros_like(suspectSkinArea)
        thirdErosionArea = np.zeros_like(suspectSkinArea)
        x,y,z = self.shapeFF

        for i in xrange(x):
            firstErosionSkinArea[i,:,:]= ndimage.binary_erosion(suspectSkinArea[i,:,:],\
            iterations = iternum).astype(suspectSkinArea.dtype)
        
        for i in xrange(x):
            secondClossingArea[i,:,:] = ndimage.binary_closing(firstErosionSkinArea[i,:,:],\
            structure = ClossingSize).astype(suspectSkinArea.dtype)
        
        for i in xrange(x):
            thirdErosionArea[i,:,:] = ndimage.binary_erosion(secondClossingArea[i,:,:],\
            structure = ErosionSize).astype(suspectSkinArea.dtype)

        afterRemoveSkin = np.zeros(np.shape(self.arrayT2SProcessed))
        afterRemoveSkin[~suspectSkinArea] = self.arrayT2SProcessed[~suspectSkinArea]
        idxRemoveSkin = (thirdErosionArea>0)
        afterRemoveSkin[idxRemoveSkin] = self.arrayT2SProcessed[~idxRemoveSkin]

        self.arrayT2SProcessed = afterRemoveSkin
        self.imageT2SProcessed.pixel_array = self.arrayT2SProcessed
        self.imageT2SProcessed.Pixelimage = self.imageT2SProcessed.pixel_array.tostring()


        if ifLogOutput != False:
            dirname = os.path.dirname(self.filepathFF)             
            name, ext = os.path.splitext(self.filepathFF)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-3], "_")
            outputdir = dirname + outBaseName + '_T2S_RSkin' + ext             
            self.imageT2SProcessed.save_as(outputdir)
        else:
            return self.imageT2SProcessed.pixel_array  

    def reduceNoise(self, filterSize, ifLogOutput):
        """ reduce noise by using median-filter
        """
        input = slef.arrayT2SProcessed
        slef.arrayT2SProcessed = ndimage.filters.median_filter(input, size=filterSize,\
         footprint=None, output=None, mode='reflect', cval=0.0, origin=0)

        self.imageT2SProcessed.pixel_array = self.arrayT2SProcessed
        self.imageT2SProcessed.Pixelimage = self.imageT2SProcessed.pixel_array.tostring()

        if ifLogOutput != False:
            dirname = os.path.dirname(self.filepathFF)             
            name, ext = os.path.splitext(self.filepathFF)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-3], "_")
            outputdir = dirname + outBaseName + '_T2S_Mfilter' + ext             
            self.imageT2SProcessed.save_as(outputdir)
        else:
            return self.imageT2SProcessed.pixel_array  