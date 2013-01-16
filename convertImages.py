#!/bin/usr/env python
# Convert all the images 
import subprocess
import shutil
import glob
import re
import os

home = "/home/yg32/Documents/PostDoc/Holodeck/Code/"
path1 = 'data'
path2 = 'processedata'
table = 'table.txt'

shutil.rmtree(home+path2)
os.makedirs(home+path2)

mask = home+"mask.png"

fldrs = glob.glob(home+path1+"/*")

fldrs = fldrs[2:3]

bashCommand = "convert -size 200x200 xc:black -fill white -draw 'circle 100,100,100,5' +antialias -quality 100% -colors 2 -colorspace Gray """+mask
subprocess.call(bashCommand,shell=True)
    
for fldr in fldrs:
    fld = fldr.split('/')[-1]
    os.makedirs(home+path2+"/"+fld)
    images = glob.glob(fldr+"/*.jpg")  
    bashCommand = """identify -format "%[fx:w]" """+images[0]
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    w =process.communicate()[0]
    w = w.replace("\n", " ")
    w = w.replace('"', "")
    w = w.replace(' ', "")
    ww = int(w)
    bashCommand = """identify -format "%[fx:h]" """+images[0]
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    h =process.communicate()[0]
    h = h.replace("\n", " ")
    h = h.replace('"', "")
    w = w.replace(' ', "")
    hh = int(h)
    # bashCommand = "convert -size "+w+"x"+h+""" xc:black -fill white -draw 'circle """+str(ww/2)+","+str(hh/2)+","+str(min([ww,hh])/2)+""",10' +antialias -quality 100% -colors 2 -colorspace Gray """+mask
    # subprocess.call(bashCommand,shell=True)   

    for image in images:
        image2 = image
        image2 = image2.replace(path1,path2)
        image2 = image2.replace(".jpg",".png")
        
        # bashCommand = "convert "+image+" -strip -set colorspace RGB "+mask+" -compose Bumpmap -composite -colorspace LAB -auto-level -format png -set colorspace RGB -sample 200x200 "+image2 
        # bashCommand = "convert -define jpeg:size=400x400 "+image+"[200x200] -strip -set colorspace RGB -colorspace LAB \( -channel R +clone -blur 0x20 -compose minus_dst -composite \) "+mask+" -compose Bumpmap -composite -auto-level -format png -colors 8 -set colorspace RGB "+image2
        #        bashCommand = "convert -define jpeg:size=400x400 "+image+"[200x200] -strip -set colorspace RGB  -contrast-stretch 15% "+image2+"; "+home+"denoise "+image2+" "+image2+"; convert "+image2+" -colorspace LAB \( -channel R +clone -blur 0x20 -compose minus_dst -composite \) "+mask+" -compose Bumpmap -composite -format png -set colorspace RGB "+image2 
        #        kaka = subprocess.Popen(bashCommand,shell=True)
        bashCommand = "convert -define jpeg:size=100x100 "+image+"[200x200] -set colorspace RGB -strip -colorspace Gray -contrast-stretch 15% \( +clone -blur 0x20 \) -compose Divide_Src -composite "+mask+" -compose Bumpmap -composite -format png "+image2

        with open(os.devnull, "w") as fnull:
            result = subprocess.call(bashCommand, stdout = fnull, stderr = fnull, shell=True)
        #        kaka = process.communicate()[0]
        # bashCommand = "convert -define jpeg:size=200x200 "+image2+"[100x100] "+image2
        # subprocess.call(bashCommand,shell=True)

        
        # execfile(home+'python1.py')

# Write table with all the image names and the correct fields 

fldrs = glob.glob(home+path2+"/*")

os.remove(home+table)
f1 = open(home+table,'w+')
titles = ('file','cuttlefish','background','year','month','day','hour','minute','second','screen','view distance','replicate')
txt =  ("%s \t"*(len(titles)-1)+"%s") % tuple(titles)
print >>f1, txt
for fldr in fldrs:
    images = glob.glob(fldr+"/*.png")
    for image in images:
        namext = os.home.basename(image)
        name = os.home.splitext(namext)[0]
        name = name.replace('Checkerboard_','Checkerboard')

        parts = re.split('_',name)
        if len(parts) == 10: parts.insert(9,'nan')
        parts.insert(0,namext)
        txt =  ("%s \t"*(len(parts)-1)+"%s") % tuple(parts)
        print >>f1, txt

f1.close()

subprocess.Popen('Rscript '+home+'try1.R',shell=True)

# images = glob.glob(home+"new"+"/*.png")  
# for image in images:
#     bashCommand = home+"ptilethresh -p 0.2 "+image+" "+image #[1:-4]+"z.png"
#     with open(os.devnull, "w") as fnull:
#             result = subprocess.call(bashCommand, stdout = fnull, stderr = fnull, shell=True)

