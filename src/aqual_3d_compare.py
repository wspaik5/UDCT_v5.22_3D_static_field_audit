import argparse, json
import numpy as np
from scipy.sparse.linalg import LinearOperator, cg
from scipy.ndimage import gaussian_filter

G=4.30091e-6  # kpc (km/s)^2 Msun^-1
A0=1.082401e-10*3.085677581e19/1e6  # (km/s)^2/kpc

def mu_of_g(g, floor):
    y=np.maximum(g/A0, floor)
    return 2*y/(np.sqrt(1+4*y*y)+1)

def solve_aqual(N=64, box=6000., floor=1e-5, outer_max=35,
                mt=1e14, gas_scale=300., gal_scale=50., sep=500.):
    dx=box/N; c=(np.arange(N)-N//2)*dx
    x=c[:,None,None]; y=c[None,:,None]; z=c[None,None,:]
    def plum(m,a,x0):
        r2=(x-x0)**2+y*y+z*z
        return 3*m/(4*np.pi*a**3)*(1+r2/a**2)**-2.5
    rho=plum(.85*mt,gas_scale,0)+plum(.075*mt,gal_scale,-sep)+plum(.075*mt,gal_scale,sep)
    r=np.sqrt(x*x+y*y+z*z+(.5*dx)**2)
    vinf2=np.sqrt(G*mt*A0)
    rb=np.sqrt(3)*(box/2)
    # Isolated monopole boundary; the additive constant is immaterial.
    boundary=vinf2*np.log(r/rb)-G*mt/r
    bmask=np.zeros((N,N,N),bool)
    bmask[[0,-1],:,:]=True; bmask[:,[0,-1],:]=True; bmask[:,:,[0,-1]]=True

    # Newtonian-like initial field with the same isolated boundary shape.
    phi=boundary.copy()
    phi[~bmask]=vinf2*np.log(r[~bmask]/rb)-G*mt/r[~bmask]
    shape=(N-2,N-2,N-2); n=np.prod(shape)
    history=[]
    mu=None
    for it in range(outer_max):
        gx,gy,gz=np.gradient(phi,dx,edge_order=2)
        target=mu_of_g(np.sqrt(gx*gx+gy*gy+gz*gz),floor)
        mu=target if mu is None else .5*mu+.5*target
        mx=.5*(mu[1:,:,:]+mu[:-1,:,:])
        my=.5*(mu[:,1:,:]+mu[:,:-1,:])
        mz=.5*(mu[:,:,1:]+mu[:,:,:-1])
        def apply_full(p):
            q=np.zeros_like(p)
            q[1:-1,1:-1,1:-1]=(
              mx[1:,1:-1,1:-1]*(p[1:-1,1:-1,1:-1]-p[2:,1:-1,1:-1])+
              mx[:-1,1:-1,1:-1]*(p[1:-1,1:-1,1:-1]-p[:-2,1:-1,1:-1])+
              my[1:-1,1:,1:-1]*(p[1:-1,1:-1,1:-1]-p[1:-1,2:,1:-1])+
              my[1:-1,:-1,1:-1]*(p[1:-1,1:-1,1:-1]-p[1:-1,:-2,1:-1])+
              mz[1:-1,1:-1,1:]*(p[1:-1,1:-1,1:-1]-p[1:-1,1:-1,2:])+
              mz[1:-1,1:-1,:-1]*(p[1:-1,1:-1,1:-1]-p[1:-1,1:-1,:-2]))/dx**2
            return q
        bfield=np.zeros_like(phi); bfield[bmask]=boundary[bmask]
        rhs=(-4*np.pi*G*rho-apply_full(bfield))[1:-1,1:-1,1:-1].ravel()
        diag=(mx[1:,1:-1,1:-1]+mx[:-1,1:-1,1:-1]+my[1:-1,1:,1:-1]+my[1:-1,:-1,1:-1]+mz[1:-1,1:-1,1:]+mz[1:-1,1:-1,:-1])/dx**2
        def mv(v):
            p=np.zeros_like(phi); p[1:-1,1:-1,1:-1]=v.reshape(shape)
            return apply_full(p)[1:-1,1:-1,1:-1].ravel()
        A=LinearOperator((n,n),matvec=mv,dtype=float)
        M=LinearOperator((n,n),matvec=lambda v:v/diag.ravel(),dtype=float)
        sol,info=cg(A,rhs,x0=phi[1:-1,1:-1,1:-1].ravel(),M=M,rtol=2e-7,atol=0,maxiter=500)
        new=phi.copy(); new[1:-1,1:-1,1:-1]=sol.reshape(shape); new[bmask]=boundary[bmask]
        rel=np.linalg.norm((new-phi)[1:-1,1:-1,1:-1])/np.linalg.norm(new[1:-1,1:-1,1:-1])
        # Nonlinear PDE residual using the newly updated field.
        phi=.35*phi+.65*new; phi[bmask]=boundary[bmask]
        history.append((it,rel,info,float(mu.min()),float(mu.max())))
        print('iter',history[-1],flush=True)
        if rel<2e-5 and it>5: break

    # Effective density inferred from the converged AQUAL potential.
    gx,gy,gz=np.gradient(phi,dx,edge_order=2)
    lap=sum(np.gradient(q,dx,axis=ax,edge_order=2) for ax,q in enumerate((gx,gy,gz)))
    rho_eff=lap/(4*np.pi*G)
    sigma=rho_eff.sum(axis=2)*dx
    j=N//2
    ratios={}; argmax={}
    for sm in (50,100,150,250):
        smap=gaussian_filter(sigma,sm/dx,mode='nearest'); p=smap[:,j]
        get=lambda pos:p[np.argmin(abs(c-pos))]
        ratios[sm]=.5*(get(-sep)+get(sep))/get(0)
        ij=np.unravel_index(np.argmax(smap),smap.shape)
        argmax[sm]=(float(c[ij[0]]),float(c[ij[1]]))
    gx,gy,gz=np.gradient(phi,dx,edge_order=2)
    raw_y=np.sqrt(gx*gx+gy*gy+gz*gz)/A0
    audit={'N':N,'box_kpc':box,'dx_kpc':dx,'floor':floor,
           'floor_fraction':float((raw_y<floor).mean()),'argmax':argmax}
    return history,ratios,audit

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--N',type=int,default=80)
    ap.add_argument('--box',type=float,default=6000); ap.add_argument('--floor',type=float,default=1e-5)
    ap.add_argument('--geometry',choices=['baseline','favorable'],default='baseline')
    ap.add_argument('--out',default='')
    a=ap.parse_args(); kw={}
    if a.geometry=='favorable': kw=dict(mt=3e14,gas_scale=500.,gal_scale=70.,sep=300.)
    h,r,audit=solve_aqual(N=a.N,box=a.box,floor=a.floor,**kw)
    result={'geometry':a.geometry,'last':h[-1],'ratios':r,'audit':audit}
    print('RESULT',json.dumps(result),flush=True)
    if a.out:
        with open(a.out,'w') as f: json.dump(result,f,indent=2)
