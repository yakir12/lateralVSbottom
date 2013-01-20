#!/bin/usr/env python
# fix the threshold images in 'new'
import subprocess
import glob
import re
home = "/home/yg32/Documents/PostDoc/Holodeck/Code/"
images = glob.glob(home+"new"+"/*.png")
for image in images:
    cmd = "convert "+image+" -morphology Close Disk:10 -morphology Smooth Disk:10 -threshold 50% -morphology EdgeIn Disk:1 "+image.replace("new","new2") 
    subprocess.call(cmd,shell=True)

