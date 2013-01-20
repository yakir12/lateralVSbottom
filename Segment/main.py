import glob
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import subprocess
import os

from numpy.linalg import eig, inv

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])

def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))

def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])


home = "/home/yg32/Documents/PostDoc/Holodeck/Code/"
path0 = 'data'
path1 = '2processedata'
fldrs = glob.glob(home+path1+"/*")
for fldr in fldrs:
    files = glob.glob(fldr+"/*.png")
    #    files = files[0:1]
    n = len(files)
    axes = np.zeros((n,2))
    center = np.zeros((n,2))
    phi = np.zeros((n,1))
    
    for i in range(n):
        im = files[i]
        I = plt.imread(im)
        sz1 = I.shape
        p = np.where(I)
        x = p[1]
        y = p[0]
        
        a = fitEllipse(x,y)
        center[i,] = ellipse_center(a)
        phi[i] = ellipse_angle_of_rotation(a)
        axes[i,] = ellipse_axis_length(a)
        # if axes[i,0] < axes[i,1]:
        #     phi2[i] = phi[i]+np.pi/2
        # else:
        #     phi2[i] = phi[i]
statind = np.indices((n,2))
idx = np.argsort(axes, axis = 1)
axes = axes[statind[0],idx]

mdat = np.ma.masked_array(axes,np.isnan(axes))
mm = np.mean(mdat,axis=0)
axes =  mm.filled(np.nan)

phi = np.concatenate((phi,phi+np.pi/2),axis=1)
phi = phi[statind[0],idx]
phi = np.delete(phi,1,1)


            
            # xx = x-center[0]
            # yy = y-center[1]

            # s = np.sqrt(xx**2+yy**2)

            # xx = xx/s
            # yy = yy/s

            # R = np.arctan2(yy,xx)

            # a, b = axes
            # xx = center[0] + a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi)
            # yy = center[1] + a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi)

            # plt.clf()
            # plt.figure
            # plt.imshow(I)
 
            # plt.plot(xx,yy,'.')
            # plt.show()
kill = np.isnan(np.sum(center,axis = 1))
for fldr in fldrs:
    files = glob.glob(fldr+"/*.png")
    #    files = files[0:1]
    n = len(files)    
    for i in range(n):
        if not kill[i]:
            im = files[i]      
            namext = os.path.basename(im)
            name = os.path.splitext(namext)[0]
            name0 = fldr.replace('2processedata','data')+"/"+name+".jpg"
            name2 = fldr.replace('2processedata','32processedata')+"/"+name+".png"
            sz0 = plt.imread(name0).shape
            # cmd = "convert "+name0+" \\( -size "+str(sz0[0])+"x"+str(sz0[1])+""" xc:black -fill white -draw "push graphic-context translate """+str(sz0[0]/sz1[0]*center[0])+","+str(sz0[1]/sz1[1]*center[1])+" rotate "+str(180*phi/np.pi)+" fill white  ellipse 0,0 "+str(sz0[0]/sz1[0]*axes[0])+","+str(sz0[1]/sz1[1]*axes[1])+""" 0,360 pop graphic-context" -colors 2 \\) -compose Bumpmap -composite -format png """ +home+"a.png"
            cmd = "convert "+name0+""" -resize 300x300 -set colorspace RGB -auto-level \\( -size 300x300 xc:black -fill white -draw "push graphic-context translate """+str(np.asscalar(center[i,0]))+","+str(np.asscalar(center[i,1]))+" rotate "+str(np.asscalar(180*phi[i]/np.pi))+" fill white  ellipse 0,0 "+str(axes[0])+","+str(axes[1])+""" 0,360 pop graphic-context" -colors 2 \\) -compose Bumpmap -composite -format png -fuzz 1% -trim -background black -rotate """+str(np.asscalar(-180*phi[i]/np.pi))+" -trim "+name2
            #cmd = "convert "+name0+""" -resize 300x300 -set colorspace RGB -auto-level """+im+""" -fill white -draw 'color """+str(center[0])+","+str(center[1])+""" floodfill' -compose Bumpmap -composite -format png -fuzz 1% -trim -background black -rotate """+str(-180*phi/np.pi)+" -trim "+name2
            subprocess.Popen(cmd,shell=True)
# sstot = sum((y-mean(y))**2)
# sserr = sum((y-yy)**2)
# r2 = 1-sserr/sstot
            # r2 = sum(sqrt((yy-y)**2+(xx-x)**2))

            # #table.writerow([im,phi,center[0],center[1],axes[0],axes[1]])
            # if r2 < 100:
            #         cmd = "convert "+im+
                        
            #         subprocess.Popen('Rscript '+home+'try1.R',shell=True)
