import pyfits
from numpy import log10,std,ceil,floor,histogram2d,meshgrid,array
from pylab import *

hdus = pyfits.open("/usr/local/surveysim/obs/spire_fls_dr2.fits")
hdus.info()
table = hdus[1].data
print hdus[1].data.names

f1 = table["F250"]
f2 = table["F350"]
f3 = table["F500"]

colors = array([log10(f1/f3),log10(f2/f3),log10(f1/f2)])
fluxes = array([log10(f1),log10(f2),log10(f3)])

y_range = array([min(colors.flatten()), max(colors.flatten())])
x_range = array([min(fluxes.flatten()),max(fluxes.flatten())])
histrange=array([x_range,y_range])
print histrange

num=1
for color in colors:
    for flux in fluxes:
        y=color
        x=flux
        dx=3.49*std(x)/pow(len(x),0.3333333)
        dy=3.49*std(y)/pow(len(y),0.3333333)
        nx=floor((x_range[1]-x_range[0])/dx)
        ny=floor((y_range[1]-y_range[0])/dy)
        H,xedges,yedges = histogram2d(x,y,range=histrange,bins=[nx,ny])
        clf()
        imshow(H, interpolation='nearest', origin='low',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        savefig("obs"+str(num)+".eps")
        num = num+1

hdus = pyfits.open("output.fits")
table = hdus[3].data

f1 = table["F1"]
f2 = table["F2"]
f3 = table["F3"]

colors = array([log10(f1/f3),log10(f2/f3),log10(f1/f2)])
fluxes = array([log10(f1),log10(f2),log10(f3)])

num=1
for color in colors:
    for flux in fluxes:
        y=color
        x=flux
        dx=3.49*std(x)/pow(len(x),0.3333333)
        dy=3.49*std(y)/pow(len(y),0.3333333)
        nx=floor((x_range[1]-x_range[0])/dx)
        ny=floor((y_range[1]-y_range[0])/dy)
        H,xedges,yedges = histogram2d(x,y,range=histrange,bins=[nx,ny])
        clf()
        imshow(H, interpolation='nearest', origin='low',
               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        savefig("model"+str(num)+".eps")
        num = num+1
