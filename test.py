# -*- coding: utf-8 -*-
import numpy as np
import SimpleITK as sitk

#filedir = "/media/louis/Volume/ResearchData/BATSegmentationData/ImageSplit/BAT_001_IM_0007_FF.dcm"
#ds = dicom.read_file(filedir)
#
#print "maximum value of FatFraction: %f"%(np.amax(np.amax(np.amax(ds.pixel_array))))
#print "minimum value of FatFraction: %f"%(np.amin(np.amin(np.amin(ds.pixel_array))))
#
##image = sitk.ReadImage(filedir)

#reader = sitk.ImageFileReader()
#reader.SetFileName(filedir)
#image1 = reader.Execute()
#imagearray = sitk.GetArrayFromImage(image1)
#print imagearray.shape
#print "maximum value of FatFraction: %f"%(np.amax(np.amax(np.amax(imagearray))))
#print "minimum value of FatFraction: %f"%(np.amin(np.amin(np.amin(imagearray))))


#from medpy.io import load
#image_data, image_header = load(filedir)
#image_data.shape
#image_data.dtype


filedir = "/media/louis/Volume/ResearchData/BATSegmentationData/ImageSplit/BAT_001_IM_0007_FF.dcm"
outdir = '/media/louis/Volume/ProgramWorkResult/BATSegresult/WFPreProcessing'
image = sitk.ReadImage(filedir)
array = sitk.GetArrayFromImage(image)
newArray = array + 1
newImage = sitk.GetImageFromArray(newArray)
newImage.SetOrigin(image.GetOrigin())
newImage.SetDirection(image.GetDirection())
newImage.SetSpacing(image.GetSpacing())
sitk.WriteImage(newImage, outdir + '/output.nrrd')

print "maximum value of FatFraction: %f"%(np.amax(np.amax(np.amax(newArray))))
print "minimum value of FatFraction: %f"%(np.amin(np.amin(np.amin(newArray))))


