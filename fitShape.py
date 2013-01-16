import glob
import numpy as np
import scipy
from numpy.linalg import eig, inv
import matplotlib.pyplot as plt
import math
from pylab import *


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
path1 = 'processedata'

fldrs = glob.glob(home+path1+"/*")

fldrs = fldrs[2:3]


# import csv
# cvsfile = open(home+'table.csv', 'wb')
# table = csv.writer(cvsfile, delimiter=',')

for fldr in fldrs:
    images = glob.glob(fldr+"/*.png")
    for im in images:

            I = plt.imread(im)
            p = np.where(I)
            x = p[0]
            y = p[1]

            a = fitEllipse(x,y)
            center = ellipse_center(a)
            phi = ellipse_angle_of_rotation(a)
            axes = ellipse_axis_length(a)
            
            xx = x-center[0]
            yy = y-center[1]

            s = sqrt(xx**2+yy**2)

            xx = xx/s
            yy = yy/s

            R = np.arctan2(yy,xx)
#R.sort()

            a, b = axes
            xx = center[0] + a*np.cos(R)*np.cos(phi) - b*np.sin(R)*np.sin(phi)
            yy = center[1] + a*np.cos(R)*np.sin(phi) + b*np.sin(R)*np.cos(phi)

# sstot = sum((y-mean(y))**2)
# sserr = sum((y-yy)**2)
# r2 = 1-sserr/sstot
            r2 = sum(sqrt((yy-y)**2+(xx-x)**2))

            #table.writerow([im,phi,center[0],center[1],axes[0],axes[1]])
            if r2 < 100:
                    cmd = "convert "+im+
                        
                    subprocess.Popen('Rscript '+home+'try1.R',shell=True)
   

cvsfile.close()

