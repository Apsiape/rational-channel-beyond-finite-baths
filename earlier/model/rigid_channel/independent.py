# Adapted from the verification package of the earlier report (Mghirbi and Douglas, Zenodo,
# doi:10.5281/zenodo.23093839; MIT licence, see LICENSE-MIT.txt), with paths set to this release.
"""Second exact arithmetic path for the finite data of the 43-generator presentation.

Does not import core.py.  A factor is (side, integer action, Laurent matrix).
Reduction scans adjacent factors, rather than using the compiler's append
algorithm.  Independence here is implementation independence, not a second
proof of the external infinite-dimensional theorems.
"""
from functools import lru_cache

ZERO=()
UNIT=((0,0,0,1),)
I=(1,0,0,0,1,0,0,0,1)
E=tuple(UNIT if k%4==0 else ZERO for k in range(9))


def need(ok,msg):
    if not ok: raise ValueError(msg)


def plus(*ps):
    d={}
    for p in ps:
        for i,j,k,n in p:
            e=i,j,k
            d[e]=(d.get(e,0)+n)%5
    return tuple((*e,n) for e,n in sorted(d.items()) if n)


def minus(p): return tuple((i,j,k,5-n) for i,j,k,n in p)


def times(p,q):
    if not p or not q: return ZERO
    d={}
    for i,j,k,n in p:
        for x,y,z,m in q:
            e=i+x,j+y,k+z
            d[e]=(d.get(e,0)+n*m)%5
    return tuple((*e,n) for e,n in sorted(d.items()) if n)


def determinant(m):
    a,b,c,d,e,f,g,h,i=m
    return plus(times(a,plus(times(e,i),minus(times(f,h)))),
                minus(times(b,plus(times(d,i),minus(times(f,g))))),
                times(c,plus(times(d,h),minus(times(e,g)))))


@lru_cache(maxsize=60000)
def inverse(m):
    need(determinant(m)==UNIT,'Laurent determinant is not one')
    a,b,c,d,e,f,g,h,i=m
    sub=lambda x,y,z,w:plus(times(x,y),minus(times(z,w)))
    return (sub(e,i,f,h),sub(c,h,b,i),sub(b,f,c,e),
            sub(f,g,d,i),sub(a,i,c,g),sub(c,d,a,f),
            sub(d,h,e,g),sub(b,g,a,h),sub(a,e,b,d))


@lru_cache(maxsize=60000)
def matrix_product(m,n):
    if m==E:return n
    if n==E:return m
    rows=[m[k:k+3] for k in (0,3,6)]
    cols=[n[k::3] for k in range(3)]
    return tuple(plus(*(times(a,b) for a,b in zip(r,c))) for r in rows for c in cols)


def det_int(a):
    return a[0]*a[4]*a[8]+a[1]*a[5]*a[6]+a[2]*a[3]*a[7]-a[2]*a[4]*a[6]-a[1]*a[3]*a[8]-a[0]*a[5]*a[7]


@lru_cache(maxsize=10000)
def inverse_int(m):
    need(det_int(m)==1,'integer determinant is not one')
    a,b,c,d,e,f,g,h,i=m
    return (e*i-f*h,c*h-b*i,b*f-c*e,f*g-d*i,a*i-c*g,c*d-a*f,d*h-e*g,b*g-a*h,a*e-b*d)


@lru_cache(maxsize=20000)
def product_int(a,b):
    return tuple(sum(x*y for x,y in zip(a[r:r+3],b[c::3])) for r in (0,3,6) for c in range(3))


@lru_cache(maxsize=60000)
def substitute(m,a):
    if a==I or m==E:return m
    result=[]
    for p in m:
        ts=[]
        for i,j,k,n in p:
            es=tuple(sum(x*y for x,y in zip(a[r:r+3],(i,j,k))) for r in (0,3,6))
            ts.append((*es,n))
        result.append(plus(tuple(ts)))
    return tuple(result)


def hfactor(f):
    _,a,m=f
    return a==I and all(min(i,j,k)>=0 for p in m for i,j,k,n in p)


def one(f):return f[1]==I and f[2]==E


def normal(seq):
    fs=[f for f in seq if not one(f)]
    while True:
        index=None
        for k in range(len(fs)-1):
            if fs[k][0]==fs[k+1][0] or hfactor(fs[k]) or hfactor(fs[k+1]):
                index=k;break
        if index is None:break
        k=index;left,right=fs[k:k+2]
        side=right[0] if hfactor(left) else left[0]
        act=product_int(left[1],right[1])
        mat=matrix_product(left[2],substitute(right[2],left[1]))
        f=side,act,mat
        fs[k:k+2]=[] if one(f) else [f]
    if len(fs)==1 and hfactor(fs[0]):fs[0]=(0,fs[0][1],fs[0][2])
    return tuple(fs)


def multiply(x,y):return normal(x+y)


def invert(x):
    out=[]
    for s,a,m in reversed(x):
        ai=inverse_int(a)
        out.append((s,ai,substitute(inverse(m),ai)))
    return normal(out)


def equal(x,y):return not multiply(invert(x),y)


def eval_word(w,gens):
    value=()
    for k in w:value=multiply(value,gens[k])
    return value


def generator_images(p):
    out={}
    for record in p['generators']:
        idx=record['id'];key=record['key']
        if key[0]=='witness':continue
        if key[0]=='p':
            _,f,r=key;side=0
        elif key[0]=='n':
            _,side,f,r=key
        elif key[0]=='t':
            _,side,i,j=key
            a=list(I);a[3*i+j]=1
            out[idx]=((side,tuple(a),E),)
            out[-idx]=invert(out[idx]);continue
        else:raise ValueError('unknown generator kind')
        root=((0,1),(1,2),(2,0))[f]
        exp=[0,0,0]
        if r:exp[abs(r)-1]=1 if r>0 else -1
        m=list(E);m[3*root[0]+root[1]]=((*exp,1),)
        out[idx]=normal(((side,I,tuple(m)),))
        out[-idx]=invert(out[idx])
    j=p['witness_generator']
    out[j]=eval_word(p['witness_word'],out);out[-j]=invert(out[j])
    return out


def decode(nf):
    out=[]
    for s in nf:
        a=tuple(s['action']);m=tuple(tuple(tuple(t) for t in p) for p in s['matrix'])
        need(len(a)==len(m)==9,'matrix length')
        need(det_int(a)==1 and determinant(m)==UNIT,'factor determinants')
        for p in m:
            need(all(len(t)==4 and all(type(v)==int for v in t) and 1<=t[3]<=4 for t in p),'polynomial coefficient encoding')
            need(plus(p)==p,'polynomial is not canonical')
        need(s['side'] in (0,1),'factor side')
        out.append((s['side'],a,m))
    return normal(out)
