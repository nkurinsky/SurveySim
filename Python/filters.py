#!/usr/bin/env python

import os

def read_filters():
    codedir=os.getcwd()+'/trunk/';
#read-in filter transmission curves/ use the FSPS filter set
    with open (codedir+'filters/FILTER_LIST','r') as f: 
        flines=f.readlines()

    filter_id,filter_choices=[],[]

    for fline in flines:
        if fline.strip():
            ncol=len(fline.split())
            if 3<ncol<6:
                tmp1,tmp2,tmp3,tmp4=fline.split()[0:4]
                if (tmp4 != '(from') and (tmp2 != '2MASS') and (tmp2 != 'Steidel') and (tmp2 != 'Stromgren') and (tmp2 != 'Idealized') and (tmp2 != 'SCUBA') and (tmp2 != 'JWST') :
                    filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3+'_'+tmp4)
                else:
                    filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3)
            if ncol >= 6:
                tmp1,tmp2,tmp3,tmp4,tmp5=fline.split()[0:5]
                if (tmp4 != '(from') and (tmp3 != 'IRAC') and (tmp3 != 'WFC3') and (tmp2 != '2MASS') and (tmp2 != 'Steidel') and (tmp2 != 'Stromgren') and (tmp2 != 'Idealized') and (tmp2 != 'SCUBA') and (tmp2 != 'JWST') :
                    filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3+'_'+tmp4)
                else:
                    if (tmp3 != 'IRAC') and (tmp3 != 'WFC3'):
                        filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3)
                if (tmp3 == 'IRAC') or (tmp3 == 'WFC3'):
                    filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3+'_'+tmp5)
            if ncol == 3:
                tmp1,tmp2,tmp3=fline.split()[0:3]
                filter_id.append(int(tmp1));filter_choices.append(tmp2+'_'+tmp3)

    return filter_id,filter_choices   

def fill_filters(f_id):
    codedir=os.getcwd()+'/trunk/';
#read-in filter transmission curves for selected bands
    with open (codedir+'filters/allfilters.dat','r') as f: 
        flines=f.readlines()

    fcount=-1
    lam1,lam2,lam3,trans1,trans2,trans3=[],[],[],[],[],[]
    for fline in flines:
        if fline.strip():
            ncol=len(fline.split())
            if ncol > 2:
                fcount=fcount+1
            if ncol == 2:
                if fcount == f_id[0]: 
                    tmp1,tmp2=fline.split()[0:2]
                    lam1.append(float(tmp1)),trans1.append(float(tmp2))
                if fcount == f_id[1]:
                    tmp1,tmp2=fline.split()[0:2]
                    lam2.append(float(tmp1)),trans2.append(float(tmp2))
                if fcount == f_id[2]:
                    tmp1,tmp2=fline.split()[0:2]
                    lam3.append(float(tmp1)),trans3.append(float(tmp2))

    return lam1,lam2,lam3,trans1,trans2,trans3
