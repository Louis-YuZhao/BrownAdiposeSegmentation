# -*- coding: utf-8 -*-
import dicom
ds = dicom.read_file("/media/louis/Volume/ResearchData/BATSegmentationData/True/IM_0007.dcm")
ds.pixel_array

