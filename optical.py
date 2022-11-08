# This Python file uses the following encoding: utf-8

#Optical
#version="1.14-rc"

import os
from functions import *
from ScatteringMatrix import *

real_path = os.path.realpath(__file__)
real_dir = os.path.dirname(real_path)
IndexDir = os.path.split(real_dir)[0]+"\\optical\\n\\"

def create_layer(name,file1,file2="Not_selected",file3="Not_selected",fr1="100",fr2="0",fr3="0",thickness="0",incoherent="0",roughness="0",ln=[],lk=[]):
    # Create and return a new Layer object
    return Layer(name,file1,file2,file3,fr1,fr2,fr3,thickness,incoherent,roughness,ln=[],lk=[])

def create_multilayer(layers,top_air=True,bottom_air=True,start_l=200.,end_l=2000.):
    # Create new Multilayer object
    multilayer = Multilayer(filename="Untitled",name="Untitled",saved=True,lmin=start_l,lmax=end_l)
    if top_air:
        # Add "top" layer (air)
        multilayer.add(Layer(name="TopLayer",file1=IndexDir+"air.in3",fr1="100"))
    # Add all layers passed as list argument
    for layer in layers:
        multilayer.add(layer)
    if bottom_air:
        # Add "bottom" layer (air)
        multilayer.add(Layer(name="TopLayer",file1=IndexDir+"air.in3",fr1="100"))
    _ema(multilayer)
    return multilayer

def compute(multilayer,Fi,start_l=200.,end_l=2000.,PlotPoints=500):
    # Compute R,T for multilayer over lam wavelength range with Fi angle of incidence
    # wavelengths in nanometres
    # Fi in degrees

    # Calculate n and k for a mixed material using effective medium approximation (EMA)
    # or just use single material values from input file
    _ema(multilayer)
    CheckWaveRange(multilayer)
    start=max(start_l,multilayer.lmin)
    end=min(end_l,multilayer.lmax)

    if end > start:
        Fi=pi*Fi/180.0
        lam=[]
        for i in range(PlotPoints):
            lam.append(start+i*(end-start)/(PlotPoints-1))
        lam[-1]=end#to be sure last point is exactly end

        structure=PrepareList(multilayer,lam)
        R,T=ComputeRT(structure,lam,Fi)
        return lam,R,T
    else:
        print('No points to plot ! Check wavelength range.')


def _ema(multilayer):
    try:
        EMA(multilayer)
        CheckWaveRange(multilayer)
    except:
        for layer in multilayer.layer:
            test=Multilayer()
            test.add(layer)
            try:
                EMA(test)
            except:
                print("Layer "+layer.name+" has index files missing or corrupted, please check !")
