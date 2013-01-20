#!/usr/bin/env python
# Write table with all the image names and the correct fields 
import csv
import glob
import re
path = "/home/yg32/Documents/PostDoc/Holodeck/Code/"
fldrs = glob.glob(path+"2processedata/*")
table = 'table.csv'
f1 = open(path+table, 'w')
tab = csv.writer(f1)
titles = ('file','cuttlefish','background','year','month','day','hour','minute','second','screen','view distance','replicate','folder')
tab.writerow(titles)
for fldr in fldrs:
    images = glob.glob(fldr+"/*.png")
    for image in images:
        namext = os.path.basename(image)
        name = os.path.splitext(namext)[0]
        name = name.replace('Checkerboard_','Checkerboard')
        parts = re.split('_',name)
        if len(parts) == 10: parts.insert(9,'nan')
        parts.insert(0,namext)
        parts.append(fldr)
        tab.writerow(parts)
f1.close()
