# -*- coding: utf-8 -*-
#import dicom
#ds = dicom.read_file("/media/louis/Volume/ResearchData/BATSegmentationData/True/IM_0007.dcm")
#ds.pixel_array

import os
for root, dirs, files in os.walk("/media/louis/Volume/PersonalData/Deutsche/Netzwerk1.2", topdown=False):
    for name in files:
        print(os.path.join(root, name))
    for name in dirs:
        print(os.path.join(root, name))