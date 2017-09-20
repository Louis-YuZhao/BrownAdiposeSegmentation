#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
class for Brown Adipose Tissue preprocessing

Created on Fri Sep  8 08:26:03 2017

@author: louis
"""

import numpy as np
import SimpleITK as sitk
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

    def computeWFmap(self, ifLogOutput):
        """ canculate the WF map """
        self.arrayFF[self.arrayFF < -10] = -10
        self.arrayFF[self.arrayFF > 110] = 110
        maxitemFF = (np.amax(np.amax(np.amax(self.arrayFF))))
        self.arrayWF = (maxitemFF - 10) - self.arrayFF
        self.shapeWF = np.shape(self.arrayWF)
        self.imageWF = sitk.GetImageFromArray(self.arrayWF)
        self.imageWF.SetOrigin(self.imageFF.GetOrigin())                               
        self.imageWF.SetSpacing(self.imageFF.GetSpacing())                                
        self.imageWF.SetDirection(self.imageFF.GetDirection())   
                        
        if ifLogOutput != False:
            dirname = self.outputdir + '/'            
            name, ext = os.path.splitext(self.filepathFF)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-1], "_")
            outputdirection = dirname + outBaseName + '_WF' + '.nrrd'             
            sitk.WriteImage(self.imageWF, outputdirection)
        else:
            return self.arrayWF
    
    def removeBackgroundmap(self, threshold, ifLogOutput):
        """ remove background of WF map 
            remove the bone and the air
        """
        waterAndfatArray = self.arrayW + self.arrayF
        idx = (waterAndfatArray <= threshold) # idx of the bone and the air area. 
        self.arrayWF[idx] = 0
        self.imageWF = sitk.GetImageFromArray(self.arrayWF)
        self.imageWF.SetOrigin(self.imageFF.GetOrigin())                               
        self.imageWF.SetSpacing(self.imageFF.GetSpacing())                                
        self.imageWF.SetDirection(self.imageFF.GetDirection())   

        if ifLogOutput != False:
            dirname = self.outputdir + '/'            
            name, ext = os.path.splitext(self.filepathFF)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-1], "_")
            outputdirection = dirname + outBaseName + '_WF_RBG' + '.nrrd'             
            sitk.WriteImage(self.imageWF, outputdirection)
        else:
            return self.arrayWF       
    
    def removeSkinVoxels(self, FFthreshold, T2Sthreshold, iternum, filterSize,\
                         ClossingSize, ErosionSize, ifLogOutput):
        """ remove the skin voxel
        """
        # skin FFthreshold <= 30%, R2_star >= 80 s-1 (T2_star <= 0.0125S)
        idxFF = (self.arrayFF <= FFthreshold)
        idxT2S = (self.arrayT2S >= T2Sthreshold)
        safeSkinArea = ~(idxFF & idxT2S)
        firstErosionSkinArea = np.zeros_like(safeSkinArea)
        secondClossingArea = np.zeros_like(safeSkinArea)
        thirdErosionArea = np.zeros_like(safeSkinArea)
        x,y,z = self.shapeFF

        for i in xrange(x):
            firstErosionSkinArea[i,:,:]= ndimage.binary_erosion(safeSkinArea[i,:,:],\
            iterations = iternum).astype(safeSkinArea.dtype)
        
        for i in xrange(x):
            secondClossingArea[i,:,:] = ndimage.binary_closing(firstErosionSkinArea[i,:,:],\
            structure = ClossingSize).astype(safeSkinArea.dtype)
        
#        secondClossingArea[i,:,:] = ndimage.filters.median_filter(secondClossingArea[i,:,:], size=filterSize,\
#         footprint=None, output=None, mode='reflect', cval=0.0, origin=0)
        
        for i in xrange(x):
            thirdErosionArea[i,:,:] = ndimage.binary_erosion(secondClossingArea[i,:,:],\
            structure = ErosionSize).astype(safeSkinArea.dtype)

        afterRemoveSkin = np.zeros_like(self.arrayWF)
        finalIndex = (safeSkinArea > 0)
        afterRemoveSkin[finalIndex] = self.arrayWF[finalIndex]
        
        self.arrayWF = afterRemoveSkin
        self.imageWF = sitk.GetImageFromArray(self.arrayWF)
        self.imageWF.SetOrigin(self.imageFF.GetOrigin())                               
        self.imageWF.SetSpacing(self.imageFF.GetSpacing())                                
        self.imageWF.SetDirection(self.imageFF.GetDirection())   

        if ifLogOutput != False:
            dirname = self.outputdir + '/'            
            name, ext = os.path.splitext(self.filepathFF)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-1], "_")
            outputdirection = dirname + outBaseName + '_WF_RSkin' + '.nrrd'             
            sitk.WriteImage(self.imageWF, outputdirection)
        else:
            return self.arrayWF              
#----------------------------------------------------------------------
#        # afterRemoveSkin
#        self.arrayTest = np.float32(afterRemoveSkin)
#        self.imageTest = sitk.GetImageFromArray(self.arrayTest)
#        self.imageTest.SetOrigin(self.imageFF.GetOrigin())                               
#        self.imageTest.SetSpacing(self.imageFF.GetSpacing())                                
#        self.imageTest.SetDirection(self.imageFF.GetDirection())   
#
#        if ifLogOutput != False:
#            dirname = self.outputdir + '/'            
#            name, ext = os.path.splitext(self.filepathFF)
#            inBaseName = os.path.basename(name)
#            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
#            outputdirection = dirname + outBaseName + '_WF_Test_afterRemoveSkin' + '.nrrd'             
#            sitk.WriteImage(self.imageTest, outputdirection)
#        else:
#            return self.arrayTest
#----------------------------------------------------------------------
#------------------------------------------------------------------------------        
#        # idxFF
#        self.arrayTest = np.float32(idxFF)
#        self.imageTest = sitk.GetImageFromArray(self.arrayTest)
#        self.imageTest.SetOrigin(self.imageFF.GetOrigin())                               
#        self.imageTest.SetSpacing(self.imageFF.GetSpacing())                                
#        self.imageTest.SetDirection(self.imageFF.GetDirection())   
#
#        if ifLogOutput != False:
#            dirname = self.outputdir + '/'            
#            name, ext = os.path.splitext(self.filepathFF)
#            inBaseName = os.path.basename(name)
#            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
#            outputdirection = dirname + outBaseName + '_WF_Test_idxFF' + '.nrrd'             
#            sitk.WriteImage(self.imageTest, outputdirection)
#        else:
#            return self.arrayTest
#        
#        # idxT2S
#        self.arrayTest = np.float32(idxT2S)
#        self.imageTest = sitk.GetImageFromArray(self.arrayTest)
#        self.imageTest.SetOrigin(self.imageFF.GetOrigin())                               
#        self.imageTest.SetSpacing(self.imageFF.GetSpacing())                                
#        self.imageTest.SetDirection(self.imageFF.GetDirection())   
#
#        if ifLogOutput != False:
#            dirname = self.outputdir + '/'            
#            name, ext = os.path.splitext(self.filepathFF)
#            inBaseName = os.path.basename(name)
#            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
#            outputdirection = dirname + outBaseName + '_WF_Test_idxT2S' + '.nrrd'             
#            sitk.WriteImage(self.imageTest, outputdirection)
#        else:
#            return self.arrayTest
#        
#        # suspectSkinArea        
#        self.arrayTest = np.float32(safeSkinArea)
#        self.imageTest = sitk.GetImageFromArray(self.arrayTest)
#        self.imageTest.SetOrigin(self.imageFF.GetOrigin())                               
#        self.imageTest.SetSpacing(self.imageFF.GetSpacing())                                
#        self.imageTest.SetDirection(self.imageFF.GetDirection())   
#
#        if ifLogOutput != False:
#            dirname = self.outputdir + '/'            
#            name, ext = os.path.splitext(self.filepathFF)
#            inBaseName = os.path.basename(name)
#            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
#            outputdirection = dirname + outBaseName + '_WF_Test_suspectSkinArea' + '.nrrd'             
#            sitk.WriteImage(self.imageTest, outputdirection)
#        else:
#            return self.arrayTest
#        
#        # firstErosionSkinArea        
#        self.arrayTest = np.float32(firstErosionSkinArea)
#        self.imageTest = sitk.GetImageFromArray(self.arrayTest)
#        self.imageTest.SetOrigin(self.imageFF.GetOrigin())                               
#        self.imageTest.SetSpacing(self.imageFF.GetSpacing())                                
#        self.imageTest.SetDirection(self.imageFF.GetDirection())   
#
#        if ifLogOutput != False:
#            dirname = self.outputdir + '/'            
#            name, ext = os.path.splitext(self.filepathFF)
#            inBaseName = os.path.basename(name)
#            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
#            outputdirection = dirname + outBaseName + '_WF_Test_firstErosionSkinArea' + '.nrrd'             
#            sitk.WriteImage(self.imageTest, outputdirection)
#        else:
#            return self.arrayTest
#        
#        # secondClossingArea        
#        self.arrayTest = np.float32(secondClossingArea)
#        self.imageTest = sitk.GetImageFromArray(self.arrayTest)
#        self.imageTest.SetOrigin(self.imageFF.GetOrigin())                               
#        self.imageTest.SetSpacing(self.imageFF.GetSpacing())                                
#        self.imageTest.SetDirection(self.imageFF.GetDirection())   
#
#        if ifLogOutput != False:
#            dirname = self.outputdir + '/'            
#            name, ext = os.path.splitext(self.filepathFF)
#            inBaseName = os.path.basename(name)
#            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
#            outputdirection = dirname + outBaseName + '_WF_Test_secondClossingArea' + '.nrrd'             
#            sitk.WriteImage(self.imageTest, outputdirection)
#        else:
#            return self.arrayTest
#        
#        # thirdErosionArea        
#        self.arrayTest = np.float32(thirdErosionArea)
#        self.imageTest = sitk.GetImageFromArray(self.arrayTest)
#        self.imageTest.SetOrigin(self.imageFF.GetOrigin())                               
#        self.imageTest.SetSpacing(self.imageFF.GetSpacing())                                
#        self.imageTest.SetDirection(self.imageFF.GetDirection())   
#
#        if ifLogOutput != False:
#            dirname = self.outputdir + '/'            
#            name, ext = os.path.splitext(self.filepathFF)
#            inBaseName = os.path.basename(name)
#            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
#            outputdirection = dirname + outBaseName + '_WF_Test_thirdErosionArea' + '.nrrd'             
#            sitk.WriteImage(self.imageTest, outputdirection)
#        else:
#            return self.arrayTest

#------------------------------------------------------------------------------
class BATpreprocessingT2S(object):
    """ class for Brown Adipose Tissue """
    def __init__(self, outputdir):
        self.outputdir = outputdir
        if not os.path.exists(self.outputdir):
            subprocess.call('mkdir ' + '-p ' + self.outputdir, shell=True)
       
    def inputImage(self, filepathFatFraction, filepathT2Star, filepathFat, filepathWater):
        
        # Fat fraction(FF)
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
    
    def removeBackgroundmap(self, threshold, ifLogOutput):
        """ remove background of WF map 
            remove the bone and the air
        """
        waterAndfatArray = self.arrayW + self.arrayF
        idx = (waterAndfatArray <= threshold) # idx of the bone and the air area. 
        self.arrayT2SProcessed = sitk.GetArrayFromImage(self.imageT2S)
        self.arrayT2SProcessed[idx] = 0

        self.imageT2SProcessed = sitk.GetImageFromArray(self.arrayT2SProcessed)
        self.imageT2SProcessed.SetOrigin(self.imageT2S.GetOrigin())                               
        self.imageT2SProcessed.SetSpacing(self.imageT2S.GetSpacing())                                
        self.imageT2SProcessed.SetDirection(self.imageT2S.GetDirection())   

        if ifLogOutput != False:
            dirname = self.outputdir + '/'             
            name, ext = os.path.splitext(self.filepathT2S)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
            outputdirection = dirname + outBaseName + '_T2S_RBG' + '.nrrd'             
            sitk.WriteImage(self.imageT2SProcessed, outputdirection)
        else:
            return self.imageT2SProcessed
  
    def removeSkinVoxels(self, FFthreshold, T2Sthreshold, iternum, filterSize,\
                         ClossingSize, ErosionSize, ifLogOutput):
        """ remove the skin voxel
        """
        # skin FFthreshold <= 30%, R2_star >= 80 s-1 (T2_star <= 0.0125S)
        idxFF = (self.arrayFF <= FFthreshold)
        idxT2S = (self.arrayT2S >= T2Sthreshold)
        safeSkinArea = ~(idxFF & idxT2S)
        firstErosionSkinArea = np.zeros_like(safeSkinArea)
        secondClossingArea = np.zeros_like(safeSkinArea)
        thirdErosionArea = np.zeros_like(safeSkinArea)
        x,y,z = np.shape(self.arrayT2S)

        for i in xrange(x):
            firstErosionSkinArea[i,:,:]= ndimage.binary_erosion(safeSkinArea[i,:,:],\
            iterations = iternum).astype(safeSkinArea.dtype)
        
        for i in xrange(x):
            secondClossingArea[i,:,:] = ndimage.binary_closing(firstErosionSkinArea[i,:,:],\
            structure = ClossingSize).astype(safeSkinArea.dtype)
        
#        secondClossingArea[i,:,:] = ndimage.filters.median_filter(secondClossingArea[i,:,:], size=filterSize,\
#         footprint=None, output=None, mode='reflect', cval=0.0, origin=0)
        
        for i in xrange(x):
            thirdErosionArea[i,:,:] = ndimage.binary_erosion(secondClossingArea[i,:,:],\
            structure = ErosionSize).astype(safeSkinArea.dtype)

        afterRemoveSkin = np.zeros_like(self.arrayT2S)
        finalIndex = (safeSkinArea > 0)
        afterRemoveSkin[finalIndex] = self.arrayT2S[finalIndex]
        
        self.arrayT2SProcessed = afterRemoveSkin
        self.imageT2SProcessed = sitk.GetImageFromArray(self.arrayT2SProcessed)
        self.imageT2SProcessed.SetOrigin(self.imageT2S.GetOrigin())                               
        self.imageT2SProcessed.SetSpacing(self.imageT2S.GetSpacing())                                
        self.imageT2SProcessed.SetDirection(self.imageT2S.GetDirection())   

        if ifLogOutput != False:
            dirname = self.outputdir + '/'                
            name, ext = os.path.splitext(self.filepathT2S)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
            outputdirection = dirname + outBaseName + '_T2S_RSkin' + '.nrrd'            
            sitk.WriteImage(self.imageT2SProcessed, outputdirection)
        else:
            return self.imageT2SProcessed   

    def reduceNoise(self, filterSize, ifLogOutput):
        """ reduce noise by using median-filter
        """        
        x,y,z = np.shape(self.arrayT2S)
        inputarray = self.arrayT2SProcessed
        # 2D median-filter
        for i in xrange(x):
            self.arrayT2SProcessed[i,:,:] = ndimage.filters.median_filter(inputarray[i,:,:], size=filterSize,\
                                  footprint=None, output=None, mode='reflect', cval=0.0, origin=0)
        # 3D median-filter
#        self.arrayT2SProcessed = ndimage.filters.median_filter(inputarray, size=filterSize,\
#         footprint=None, output=None, mode='reflect', cval=0.0, origin=0)

        self.imageT2SProcessed = sitk.GetImageFromArray(self.arrayT2SProcessed)
        self.imageT2SProcessed.SetOrigin(self.imageT2S.GetOrigin())                               
        self.imageT2SProcessed.SetSpacing(self.imageT2S.GetSpacing())                                
        self.imageT2SProcessed.SetDirection(self.imageT2S.GetDirection())   

        if ifLogOutput != False:
            dirname = self.outputdir + '/'             
            name, ext = os.path.splitext(self.filepathT2S)
            inBaseName = os.path.basename(name)
            outBaseName = string.join(inBaseName.split("_")[0:-2], "_")
            outputdirection = dirname + outBaseName + '_T2S_Mfilter' + '.nrrd'             
            sitk.WriteImage(self.imageT2SProcessed, outputdirection)
        else:
            return self.imageT2SProcessed   