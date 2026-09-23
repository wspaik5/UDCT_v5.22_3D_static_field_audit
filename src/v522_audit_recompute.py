#!/usr/bin/env python3
"""UDCT v5.22 numerical audit: Newtonian and explicit QUMOND comparisons.

This script reports sampled galaxy/gas contrasts and true global maxima after
Gaussian smoothing.  It is deliberately separate from the nonlinear AQUAL
solver so the inexpensive FFT convergence matrix can be rerun independently.
"""
import argparse, csv, json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter

G = 4.30091e-6  # kpc (km/s)^2 Msun^-1
A0 = 1.082401e-10 * 3.085677581e19 / 1e6  # (km/s)^2/kpc

GEOMS = {
    "baseline": dict(mt=1e14, gas_scale=300., gal_scale=50., sep=500.),
    "favorable": dict(mt=3e14, gas_scale=500., gal_scale=70., sep=300.),
}

def fields(N, box, geom):
    p = GEOMS[geom]; dx = box/N
    c = (np.arange(N)-N//2)*dx
    x=c[:,None,None]; y=c[None,:,None]; z=c[None,None,:]
    def plum(m,a,x0):
        r2=(x-x0)**2+y*y+z*z
        return 3*m/(4*np.pi*a**3)*(1+r2/a**2)**-2.5
    rho=(plum(.85*p['mt'],p['gas_scale'],0)+
         plum(.075*p['mt'],p['gal_scale'],-p['sep'])+
         plum(.075*p['mt'],p['gal_scale'], p['sep']))
    k=2*np.pi*np.fft.fftfreq(N,d=dx)
    kx=k[:,None,None]; ky=k[None,:,None]; kz=k[None,None,:]
    k2=kx*kx+ky*ky+kz*kz
    rhok=np.fft.fftn(rho)
    phiNk=np.zeros_like(rhok)
    nz=k2>0
    phiNk[nz]=-4*np.pi*G*rhok[nz]/k2[nz]
    gx=np.fft.ifftn(-1j*kx*phiNk).real
    gy=np.fft.ifftn(-1j*ky*phiNk).real
    gz=np.fft.ifftn(-1j*kz*phiNk).real
    gN=np.sqrt(gx*gx+gy*gy+gz*gz)
    xarg=gN/A0
    # Explicit QUMOND choice: nu(x)=sqrt(1+1/x).  Zero-field cells are
    # regularised only for floating-point evaluation and are counted below.
    xfloor=1e-12
    affected=xarg<xfloor
    nu=np.sqrt(1+1/np.maximum(xarg,xfloor))
    qx=nu*gx; qy=nu*gy; qz=nu*gz
    divq=np.fft.ifftn(1j*kx*np.fft.fftn(qx)+1j*ky*np.fft.fftn(qy)+1j*kz*np.fft.fftn(qz)).real
    rhoq=-divq/(4*np.pi*G)
    return c,dx,rho,rhoq,float(affected.mean()),p['sep']

def diagnostics(c,dx,rho3,sep,smoothing):
    sigma=rho3.sum(axis=2)*dx
    sm=gaussian_filter(sigma,smoothing/dx,mode='wrap')
    j=len(c)//2
    ix0=np.argmin(abs(c)); ixm=np.argmin(abs(c+sep)); ixp=np.argmin(abs(c-sep))
    gas=sm[ix0,j]; gal=.5*(sm[ixm,j]+sm[ixp,j])
    imax=np.unravel_index(np.argmax(sm),sm.shape)
    return dict(contrast=float(gal/gas), argmax_x=float(c[imax[0]]),
                argmax_y=float(c[imax[1]]), gas=float(gas), galaxy_mean=float(gal))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',default='v522_qumond_newton_audit.csv')
    ap.add_argument('--json',default='v522_qumond_newton_audit.json')
    args=ap.parse_args()
    cases=[]
    for geom in GEOMS:
      for box,N in [(3000.,96),(3000.,128),(4500.,128),(6000.,128),(6000.,160)]:
        c,dx,rho,rhoq,ff,sep=fields(N,box,geom)
        for sm in (50,100,150,250):
          for theory,arr in [('Newtonian',rho),('QUMOND',rhoq)]:
            d=diagnostics(c,dx,arr,sep,sm)
            cases.append(dict(geometry=geom,theory=theory,N=N,box_kpc=box,
                              dx_kpc=dx,smoothing_kpc=sm,q_floor_fraction=ff,**d))
        print(geom,N,box,'dx',dx,'qfloor',ff,flush=True)
    keys=list(cases[0])
    with open(args.out,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys); w.writeheader(); w.writerows(cases)
    Path(args.json).write_text(json.dumps(cases,indent=2))
    print(args.out,args.json)

if __name__=='__main__': main()
