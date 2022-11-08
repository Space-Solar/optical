import time
from numpy import *
from os import path, getcwd
from os import name as osname

class data:
    def __init__(self,lam=0,T=0,R=0):
        self.lam=lam
        self.T=T
        self.R=R

class Multilayer:#
    def __init__(self,filename="Untitled",name="Untitled",saved=True,lmin=200,lmax=1100):
        self.filename=filename
        self.name=name
        self.saved=saved
        self.lmin=lmin
        self.lmax=lmax
        #self.index_dir=index_dir
        #self.struture_dir=structure_dir
        self.layer=[]
    def add(self,l):
        self.layer.append(l)
        self.recalculate()
    def remove(self,l):
        self.layer.remove(l)
        self.recalculate()
    def insert(self,l):
        i=self.layer.index(l)+1
        emptylayer=layer()
        self.layer.insert(i,emptylayer)
        self.recalculate()
    def clear(self):
        del self.layer[:]
        self.filename="Untitled"
        self.name="Untitled"
        self.saved=True
    def recalculate(self):
        EMA(self)

class Layer:
    def __init__(self,name="Unnamed",file1="Not_selected",file2="Not_selected",file3="Not_selected",fr1="100",fr2="0",fr3="0",thickness="0",incoherent="0",roughness="0",ln=[],lk=[]):
        self.name=name
        self.file1=file1 # path to in3 file for first component
        self.file2=file2 # path to in3 file for second component
        self.file3=file3 # path to in3 file for third component
        self.fr1=fr1 # % fraction of first EMA
        self.fr2=fr2 # % fraction of second EMA
        self.fr3=fr3 # % fraction of third EMA
        self.thickness=thickness # layer thickness in nm
        self.incoherent=incoherent#True if the layer has to be treated as inchoerent
        self.roughness=roughness #top roughness in A -- should this be in nm? GUI says nm #TODO
        self.ln=ln
        self.lk=lk

def EMA(multilayer):
    #computes n and k for a mixed material using
    #effective medium approximation (EMA)
    for layer in multilayer.layer:
        fr1,fr2,fr3=int(layer.fr1),int(layer.fr2),int(layer.fr3)
        if fr1==100:#just one medium, no EMA
            layer.ln,layer.lk=ReadIn3(layer.file1)
        else:#EMA 2 or 3 mediums
            n,k=[],[]
            nt,kt=ReadIn3(layer.file1)
            n.append(nt)
            k.append(kt)
            nt,kt=ReadIn3(layer.file2)
            n.append(nt)
            k.append(kt)
            n.append(0)
            k.append(0)
            if fr3 == 0:#just 2 mediums
                n[2],k[2]=nt,kt
            else:#O.K. 3 mediums
                n[2],k[2]=ReadIn3(layer.file3)

            ln,lk=[],[]
            for j in range(3):
                lnt,lkt=[],[]
                for h in range(len(n[j])):
                    lnt.append(n[j][h][0])
                for h in range(len(k[j])):
                    lkt.append(k[j][h][0])
                ln.append(lnt)
                lk.append(lkt)

            lmin=max(min(ln[0]),min(ln[1]),min(ln[2]),min(lk[0]),min(lk[1]),min(lk[2]))
            lmax=min(max(ln[0]),max(ln[1]),max(ln[2]),max(lk[0]),max(lk[1]),max(lk[2]))

            ln=ln[0]+ln[1]+ln[2]+[lmin,lmax]
            lk=lk[0]+lk[1]+lk[2]+[lmin,lmax]

            ln=sortRemoveDupes(ln)
            lk=sortRemoveDupes(lk)

            lnt=ln[:]
            ln=[]
            for j in lnt:
                if lmin <= j <= lmax:
                    ln.append(j)

            lkt=lk[:]
            lk=[]
            for j in lkt:
                if lmin <= j <= lmax:
                    lk.append(j)

            nn,nk=[],[]
            fr=[fr1,fr2,fr3]
            for lam in ln:
                ne=[]
                for h in range(3):
                    nr=Interpol(n[h],lam)
                    ni=Interpol(k[h],lam)
                    ne.append(nr+ni*(1j))
                nn.append([lam, EMA3(ne, fr).real])

            for lam in lk:
                ne=[]
                for h in range(3):
                    nr=Interpol(n[h],lam)
                    ni=Interpol(k[h],lam)
                    ne.append(nr+ni*(1j))
                nk.append([lam, EMA3(ne, fr).imag])

            layer.ln=nn
            layer.lk=nk

def ReadIn3(filename):
    #load n and k from an index file *.in3
    n=[]
    k=[]
    filename=filename.replace("\\","/")#for win compatibility
    in_file=open(path.normpath(filename),"r")
    line=in_file.readline()
    line=in_file.readline()
    line=in_file.readline()
    nn=int(CleanRecord(line))
    for i in range(nn):
        line=in_file.readline()
        answer=CleanRecord(line).split()
        n.append([float(answer[0])/10.0, float(answer[1])])

    line=in_file.readline()
    nk=int(CleanRecord(line))
    for i in range(nk):
        line=in_file.readline()
        answer=CleanRecord(line).split()
        k.append([float(answer[0])/10.0, float(answer[1])])
    in_file.close()

    return n,k

def CleanRecord(a):
    if a[-2:]=="\r\n":#remove from string Ms format end of record
        a=a[:-2]
    elif a[-1:]=="\n":#remove from string Unix format end of record
        a=a[:-1]
    a=a.replace(","," ")
    return a

def CleanRecord2(a):
    if a[-2:]=="\r\n":#remove from string Ms format end of record
        a=a[:-2]
    elif a[-1:]=="\n":#remove from string Unix format end of record
        a=a[:-1]
    return a

def Interpol(a,l):
    #linear interpolation
    for i in range(len(a)):
        if a[i][0]==l:
            n=a[i][1]
            break
        if a[i][0]>l:
            n = a[i-1][1] + 1.0*(l - a[i-1][0]) / (a[i][0] - a[i-1][0]) * \
            (a[i][1] - a[i-1][1])
            break
    return n

def PrepareList(multilayer,l):
    #prepare list for globscatmatr
    b=[]
    for layer in multilayer.layer:
        n=[]
        for lam in l:
            nr=Interpol(layer.ln,lam)
            ni=Interpol(layer.lk,lam)
            n.append(nr+ni*(1j))
        #b.append([int(float(layer.thickness)), n, int(layer.incoherent), int(layer.roughness)])
        b.append([float(layer.thickness), n, int(layer.incoherent), float(layer.roughness)])
    return b

def sortRemoveDupes(lst):
    """Sort the list, and remove duplicate symbols.
    """
    if len(lst) == 0:
        return lst
    lst.sort()
    lst = [lst[0]] + [lst[i] for i in range(1, len(lst))
            if lst[i] != lst[i - 1]]
    return lst

def EMA3(n, fr):
    #effective medium approximation 3 medium
    e1,e2,e3=n[0]**2,n[1]**2,n[2]**2
    f1,f2,f3=fr[0]/100.0,fr[1]/100.0,fr[2]/100.0
    nguess=f1*n[0]+f2*n[1]+f3*n[2]
    a=-4
    b=f1*(4*e1-2*(e2+e3))+f2*(4*e2-2*(e1+e3))+f3*(4*e3-2*(e1+e2))
    c=f1*(2*e1*(e2+e3)-e2*e3)+f2*(2*e2*(e1+e3)-e1*e3)+f3*(2*e3*(e1+e2)-e1*e2)
    d=(f1+f2+f3)*(e1*e2*e3)
    e=root3(a,b,c,d)
    pn=[e[0]**(0.5),e[1]**(0.5),e[2]**(0.5)]
    distance=[abs(pn[0]-nguess), abs(pn[1]-nguess), abs(pn[2]-nguess)]
    return pn[distance.index(min(distance))]

def root3(a,b,c,d):
    #finds roots of eq ax^3+bx^2+cx+d=0
    a,b,c=b/a,c/a,d/a
    p=(-a**2)/3+b
    q=(2*a**3)/27-a*b/3+c
    u1=(-q/2+((q**2)/4+(p**3)/27)**(0.5))**(1.0/3)
    u,x=[],[]
    u.append(u1)#there are 3 roots for u
    u.append(u1*exp(2*pi/3*1j))
    u.append(u1*exp(4*pi/3*1j))
    for i in range(3):
        x.append(u[i]-p/(3*u[i])-a/3)
    return x

def CheckWaveRange(multilayer):
    #return wavelength range that can be computed
    lmin,lmax=[],[]
    for layer in multilayer.layer:
        n=layer.ln
        k=layer.lk
        ln,lk=[],[]
        for h in range(len(n)):
            ln.append(n[h][0])
        for h in range(len(k)):
            lk.append(k[h][0])
        lmin.append(max(min(ln),min(lk)))
        lmax.append(min(max(ln),max(lk)))
    multilayer.lmin=max(lmin)
    multilayer.lmax=min(lmax)

def Chi(R,T,Re,Te,st_dev):
    #compute chi square test
    Chi_R,Chi_T=0.0,0.0
    for i in range(len(R)):
        Chi_R=Chi_R+((Re[i]-R[i])/st_dev)**2
        Chi_T=Chi_T+((Te[i]-T[i])/st_dev)**2
    Chi_RT = (Chi_R+Chi_T)/(len(R)+len(T))
    Chi_R = Chi_R/len(R)
    Chi_T = Chi_T/len(T)

    return Chi_R,Chi_T,Chi_RT

def FileExists(filename):
    #check if a file exists
    try:
        in_file=open(filename,"r")
        in_file.close()
        answer=True
    except:
        answer=False
    return answer

def loadref(a):
    data=[]
    in_file=open(a,"r")
    line=in_file.readline()
    while True:
        line=in_file.readline()
        if len(line)<2:
            break
        line=line[:-1]
        v=line.split()
        data.append([float(v[0]),float(v[1])])
    in_file.close()
    return data
