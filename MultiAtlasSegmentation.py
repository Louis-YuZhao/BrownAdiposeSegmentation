#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
class for multi atlas segmentation base on elastix 

Created on Fri Sep  8 08:26:03 2017

@author: louis
"""

import numpy as np
import SimpleITK as sitk 
import os
import subprocess
import time
import string
from scipy import ndimage

class MultiAtlasSegmentation(object):
    """ class for multi atlas segmentation base on elastix """
    
    def __init__(self, atlasImageListdir, atlasLabelListdir, atlasNum):
        self.atlasImageListdir = atlasImageListdir
        self.atlasLabelListdir = atlasLabelListdir
        self.atlasNum = atlasNum

    def readAtlasImagetoList(self):
        """
        read atlas list from given atlas dir to list
        """
        self.atlasImageList = []
        with open(self.atlasImageListdir) as f:
            self.atlasImageList = f.read().splitlines()
        return self.atlasImageList
    
    def readAtlasLabeltoList(self):
        """
        read atlas list from given atlas dir to list
        """
        self.atlasLabelList = []
        with open(self.atlasLabelListdir) as f:
            self.atlasLabelList = f.read().splitlines()
        return self.atlasLabelList
   
    def __readFiletoList(self, filename):
        outputList = []
        with open(filename) as f:
            outputList = f.read().splitlines()
        return outputList

    def __WriteListtoFile(self, filelist, filename):
        with open(filename, 'w') as f:
            for i in filelist:
                f.write(i+'\n')
        return 1

    def AffineRegistartion(self, movingimagedir, imageOutputDir, labelOutputDir, affinePara):
        """
        multi atlas segmentation
        affine registration
        """
        s = time.time()
        
        logFileDir = imageOutputDir+'/RegTrans' + '_RUN_' + '.log'
        logFile = open(logFileDir, 'w')
        
        # calculating the registation parameters
        movingImage = sitk.ReadImage(movingimagedir)      
        selx = sitk.ElastixImageFilter()
        selx.SetMovingImage(movingImage)
        
        paraTrans = sitk.GetDefaultParameterMap("translation")
        paraAffine = sitk.GetDefaultParameterMap("affine")
        paraTrans['Metric'] = ['AdvancedNormalizedCorrelation']
        paraAffine['Metric'] = ['AdvancedNormalizedCorrelation']
        for key, value in affinePara.items():
            paraAffine[key] = [value]
        
        RegImOutputList = []
        for i in range(len(self.atlasImageList)):
            fixedIm = self.atlasImageList[i]
            name, ext = os.path.splitext(fixedIm)
            baseName = os.path.basename(name)
            transformPara = imageOutputDir + '/ImAffineReg_' + baseName
            if not os.path.exists(transformPara):
                subprocess.call('mkdir ' + '-p ' + transformPara, shell=True)
        
            # Run registration                      
            fixedImage = sitk.ReadImage(fixedIm)
            selx.SetFixedImage(fixedImage)
            selx.SetLogFileName(logFileDir)
            selx.SetOutputDirectory(transformPara)                        
            selx.SetParameterMap(paraTrans)
            selx.AddParameterMap(paraAffine)
            selx.Execute()
            
            # write the result image to the file 
            RegImOutput= imageOutputDir+'/ImAffineReg_' + baseName + '.nrrd'
            RegImOutputList.append(RegImOutput)
            sitk.WriteImage(selx.GetResultImage(), RegImOutput)
            print ("the %d th registration is finished" % (i))
            
        self.__WriteListtoFile(RegImOutputList, imageOutputDir+"/FileList.txt")  
        
        # aligning the labels
        stran = sitk.TransformixImageFilter()

        transformParaMaplist = []    
        for i in self.atlasLabelList:
            name, ext = os.path.splitext(i)
            labelBaseName = os.path.basename(name)
            imageBaseName = string.join(labelBaseName.split("_")[0:-2], "_")
            transformParaMaplist.append(imageOutputDir+'/ImAffineReg_' + imageBaseName) 

        RegLabOutputList = []
        
        for i in range(len(self.atlasLabelList)):              
            
            movingIm = self.atlasLabelList[i]
            movingImage = sitk.ReadImage(movingIm)
                    
            transformParadir_Trans = transformParaMaplist[i] + '/TransformParameters.0.txt'       
            paraTrans = sitk.ReadParameterFile(transformParadir_Trans)
            paraTrans['FinalBSplineInterpolationOrder'] = ['0']
            
            transformParadir_Affine = transformParaMaplist[i] + '/TransformParameters.1.txt' 
            paraAffine = sitk.ReadParameterFile(transformParadir_Affine)
            paraAffine['FinalBSplineInterpolationOrder'] = ['0']

            stran.SetMovingImage(movingImage)
            stran.SetTransformParameterMap(paraTrans)
            stran.AddTransformParameterMap(paraAffine)
            stran.SetLogFileName(logFileDir)
            stran.SetOutputDirectory(labelOutputDir)
            stran.Execute()        

            # write the result image to the file 
            name, ext = os.path.splitext(movingIm)
            labelBaseName = os.path.basename(name)
            LabImOutput = labelOutputDir + '/LabelAffineReg_' + labelBaseName + '.nrrd'
            RegLabOutputList.append(LabImOutput)
            sitk.WriteImage(stran.GetResultImage(), LabImOutput)
            print ("the %d th registration is finished" % (i))  

        self.__WriteListtoFile(RegLabOutputList, labelOutputDir+"/FileList.txt")    

        e = time.time()
        l = e - s

        print 'affine registration is finished/n'
        print 'Total running time: %f mins'%(l/60.0)

    def BsplineRegistartion(self, movingimagedir, imageOutputDir, labelOutputDir, BsplinePara):
        """
        multi atlas segmentation
        Bspline registration
        """
        s = time.time()
        
        logFileDir = imageOutputDir+'/RegTrans' + '_RUN_' + '.log'
        logFile = open(logFileDir, 'w')
        
        # calculating the registation parameters
        movingImage = sitk.ReadImage(movingimagedir)      
        selx = sitk.ElastixImageFilter()
        selx.SetMovingImage(movingImage)
        
        paraBspline = sitk.GetDefaultParameterMap("Bspline")
        paraBspline['Metric'] = ['AdvancedNormalizedCorrelation']
        for key, value in BsplinePara.items():
            paraBspline[key] = [value]
        
        RegImOutputList = []
        for i in range(len(self.atlasImageList)):
            fixedIm = self.atlasImageList[i]
            name, ext = os.path.splitext(fixedIm)
            baseName = os.path.basename(name)
            transformPara = imageOutputDir + '/ImBsplineReg_' + baseName
            if not os.path.exists(transformPara):
                subprocess.call('mkdir ' + '-p ' + transformPara, shell=True)
        
            # Run registration                      
            fixedImage = sitk.ReadImage(fixedIm)
            selx.SetFixedImage(fixedImage)
            selx.SetLogFileName(logFileDir)
            selx.SetOutputDirectory(transformPara)                        
            selx.AddParameterMap(paraBspline)
            selx.Execute()
            
            # write the result image to the file 
            RegImOutput= imageOutputDir+'/ImBsplineReg_' + baseName + '.nrrd'
            RegImOutputList.append(RegImOutput)
            sitk.WriteImage(selx.GetResultImage(), RegImOutput)
            print ("the %d th registration is finished" % (i))
            
        self.__WriteListtoFile(RegImOutputList, imageOutputDir+"/FileList.txt")  
        
        # aligning the labels
        stran = sitk.TransformixImageFilter()

        transformParaMaplist = []    
        for i in self.atlasLabelList:
            name, ext = os.path.splitext(i)
            labelBaseName = os.path.basename(name)
            imageBaseName = string.join(labelBaseName.split("_")[0:-2], "_")
            transformParaMaplist.append(imageOutputDir+'/ImBsplineReg_' + imageBaseName) 

        RegLabOutputList = []
        
        for i in range(len(self.atlasLabelList)):              
            
            movingIm = self.atlasLabelList[i]
            movingImage = sitk.ReadImage(movingIm)
                   
            transformParadir_Bspline = transformParaMaplist[i] + '/TransformParameters.0.txt' 
            paraBspline = sitk.ReadParameterFile(transformParadir_Bspline)
            paraBspline['FinalBSplineInterpolationOrder'] = ['0']

            stran.SetMovingImage(movingImage)
            stran.AddTransformParameterMap(paraBspline)
            stran.SetLogFileName(logFileDir)
            stran.SetOutputDirectory(labelOutputDir)
            stran.Execute()        

            # write the result image to the file 
            name, ext = os.path.splitext(movingIm)
            labelBaseName = os.path.basename(name)
            LabImOutput = labelOutputDir + '/LabelBsplineReg_' + labelBaseName + '.nrrd'
            RegLabOutputList.append(LabImOutput)
            sitk.WriteImage(stran.GetResultImage(), LabImOutput)
            print ("the %d th registration is finished" % (i))  

        self.__WriteListtoFile(RegLabOutputList, labelOutputDir+"/FileList.txt")    

        e = time.time()
        l = e - s
        print 'Bspline registration is finished/n'
        print 'Total running time: %f mins'%(l/60.0)

    def FusionofSegmentation(self,labelListdir,threshold,fusedLabeldir):
        labelList = self.__readFiletoList(labelListdir)
        image = sitk.ReadImage(labelList[0])
        posibilityMap = np.zeros_like(sitk.GetArrayFromImage(image))
        for label in labelList:
            posibilityMap += sitk.GetArrayFromImage(sitk.ReadImage(label))
        posibilityMap = posibilityMap/np.float(len(labelList))

        image_array = posibilityMap 
        image_array[image_array > threshold] = np.uint8(1)
        image_array[image_array <= threshold] = np.uint8(0)

        img = sitk.GetImageFromArray(image_array)
        img.SetOrigin(image.GetOrigin())
        img.SetSpacing(image.GetSpacing())
        img.SetDirection(image.GetDirection())
        sitk.WriteImage(img,fusedLabeldir)