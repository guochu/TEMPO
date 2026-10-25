import numpy as np, itertools

SX = np.array([[0,1],[1,0]],complex)
SY = np.array([[0,-1j],[1j,0]],complex)
SZ = np.diag([1.,-1.]).astype(complex)
HS = .37*SX+.19*SY+.23*SZ
S = np.diag([-.7,.9]).astype(complex)
SV = np.diag(S).real
RHO = (np.eye(2)+.3*SX+.4*SY+.5*SZ)/2
OMEGA, G, BETA, T = 1.3, .35, 1.7, 1.6
NBAR = 1/np.expm1(BETA*OMEGA)

def prop(X, dt):
    """exp(-1j*dt*X) for Hermitian X."""
    w, v = np.linalg.eigh(X)
    return (v*np.exp(-1j*dt*w))@v.conj().T

def eta_analytic(dt,lag):
    a=OMEGA
    if lag==0:
        q=(1-1j*a*dt-np.exp(-1j*a*dt))/a**2
        return G**2*((NBAR+1)*q+NBAR*q.conjugate())
    q=2*(1-np.cos(a*dt))/a**2
    return G**2*q*((NBAR+1)*np.exp(-1j*a*lag*dt)+NBAR*np.exp(1j*a*lag*dt))

def if_eta(n,dt):
    p=list(itertools.product(range(2),repeat=n))
    out=np.zeros((len(p),len(p)),complex)
    for ip,pp in enumerate(p):
        for im,pm in enumerate(p):
            sp,sm=[SV[i] for i in pp],[SV[i] for i in pm]; z=0j
            for j in range(n):
                for k in range(j+1):
                    e=eta_analytic(dt,j-k)
                    z-=(sp[j]-sm[j])*(e*sp[k]-e.conjugate()*sm[k])
            out[ip,im]=np.exp(z)
    return out

def quapi_sum(f,n,dt,order='first',rho=RHO):
    p=list(itertools.product(range(2),repeat=n))
    u=prop(HS,dt); half=prop(HS,dt/2)
    r0=half@rho@half.conj().T if order=='second' else rho
    rt=np.zeros((2,2),complex)
    for ip,pp in enumerate(p):
        ap=np.prod([u[pp[j],pp[j-1]] for j in range(1,n)])
        for im,pm in enumerate(p):
            am=np.prod([u[pm[j],pm[j-1]] for j in range(1,n)])
            rt[pp[-1],pm[-1]]+=ap*am.conjugate()*r0[pp[0],pm[0]]*f[ip,im]
    last=half if order=='second' else u
    return last@rt@last.conj().T

def bath_matrices(nb):
    b=np.diag(np.sqrt(np.arange(1,nb)),1)
    hb=OMEGA*np.diag(np.arange(nb))
    bb=G*(b+b.T)
    p=np.exp(-BETA*OMEGA*np.arange(nb)); p/=p.sum()
    return hb,bb,np.diag(p)

def ed(nb,dt=None,order=None):
    hb,bb,rb=bath_matrices(nb)
    a=np.kron(HS,np.eye(nb))
    q=np.kron(np.eye(2),hb)+np.kron(S,bb)
    r=np.kron(RHO,rb)
    if order is None:
        u=prop(a+q,T)
    else:
        n=round(T/dt)
        if order=='first': step=prop(a,dt)@prop(q,dt)
        elif order=='second':
            h=prop(a,dt/2); step=h@prop(q,dt)@h
        u=np.linalg.matrix_power(step,n)
    x=u@r@u.conj().T
    return np.einsum('aibi->ab',x.reshape(2,nb,2,nb))

if __name__=='__main__':
    nb=32
    exact=ed(nb)
    print('trace(exact)=',np.trace(exact))
    for n in [4,8]:
        dt=T/n; f=if_eta(n,dt)
        qf=quapi_sum(f,n,dt,'first')
        qs=quapi_sum(f,n,dt,'second')
        rf=ed(nb,dt,'first'); rs=ed(nb,dt,'second')
        print(f'N={n}: |pathsum_first - ed_first|   = {np.linalg.norm(qf-rf):.3e}')
        print(f'       |pathsum_second - ed_second| = {np.linalg.norm(qs-rs):.3e}')
        print(f'       |pathsum_first - exact|  = {np.linalg.norm(qf-exact):.3e}')
        print(f'       |pathsum_second - exact| = {np.linalg.norm(qs-exact):.3e}')
        print(f'       err ratio first/second  = {np.linalg.norm(qf-exact)/np.linalg.norm(qs-exact):.1f}')
