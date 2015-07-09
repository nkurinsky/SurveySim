import os

def filterDir():
    ssdir=os.getenv("SURVEYSIMPATH")
    if(ssdir != None):
        if(os.path.exists(ssdir)):
            return ssdir+"/filters/"
    instdir='/usr/local/surveysim/filters/'
    instdir2=os.getenv("HOME")+"/local/surveysim/filters/"
    if(os.path.exists(instdir)):
        return instdir
    if(os.path.exists(instdir2)):
        return instdir2
    else:
        wd=os.getcwd()
        topdir="SurveySim"
        index=wd.find(topdir)+len(topdir)
        return wd[0:index]+'/trunk/filters/'

def read_filters():
    #read-in filter transmission curves/ use the FSPS filter set
    with open (filterDir()+'FILTER_LIST','r') as f:
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
    #read-in filter transmission curves for selected bands
    with open (filterDir()+'allfilters.dat','r') as f: 
        flines=f.readlines()

    nfilt=len(f_id)
    fcount=0
    if(nfilt ==3):
        lam1,lam2,lam3,trans1,trans2,trans3=[],[],[],[],[],[]
    if(nfilt ==4):
        lam1,lam2,lam3,lam4,trans1,trans2,trans3,trans4=[],[],[],[],[],[],[],[]

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
                if(nfilt == 4):
                    if fcount == f_id[3]:
                        tmp1,tmp2=fline.split()[0:2]
                        lam4.append(float(tmp1)),trans4.append(float(tmp2))

 
    if nfilt == 1:
        return lam1,trans1
    if nfilt == 2:
        return lam1,lam2,trans1,trans2
    if nfilt == 3:
        return lam1,lam2,lam3,trans1,trans2,trans3
    if nfilt == 4:
        return lam1,lam2,lam3,lam4,trans1,trans2,trans3,trans4

def getFilterID(name):
    ids,names=read_filters()
    for i in range(0,len(ids)):
        if(name in names[i]):
            return ids[i],names[i]
    raise LookupError("Filter \""+name+"\" not found in library")

def getFilterName(fid):
    ids,names=read_filters()
    for i in range(0,len(ids)):
        if(fid == ids[i]):
            return names[i]
    raise LookupError("Filter ID \""+str(fid)+"\" not found in library")

def getFilterIDs(instrument):
    ids=[]
    names=[]
    fid,fnames=read_filters()
    for i in range(0,len(fid)):
        if(instrument in fnames[i]):
            ids.append(fid[i])
            names.append(fnames[i])
    return ids,names
