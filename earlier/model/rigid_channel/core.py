# Adapted from the verification package of the earlier report (Mghirbi and Douglas, Zenodo,
# doi:10.5281/zenodo.23093839; MIT licence, see LICENSE-MIT.txt), with paths set to this release.
"""Exact Laurent-matrix amalgam and finite relator compiler for the 43-generator presentation.

Only Python integer arithmetic is used.  Polynomial tuples are sorted
(e1,e2,e3,coefficient mod 5); matrices are row-major 3 by 3 tuples.
Group products compose from right to left and [x,y]=xyx^-1y^-1.
The quotient used to decide the finite label equalities is
 (SL_3(F5[x1^+/-1,x2^+/-1,x3^+/-1]) rtimes SL_3(Z))
 amalgamated over SL_3(F5[x1,x2,x3]).
This module is not a test of the external centralizer theorem.
"""
from functools import lru_cache
from itertools import combinations, product

P = 5
Z = ()
ONE = ((0, 0, 0, 1),)
IMAT = tuple(ONE if i == j else Z for i in range(3) for j in range(3))
IACT = (1,0,0,0,1,0,0,0,1)
GE = (IMAT, IACT)
ROOTS = ((0,1),(1,2),(2,0))
PAIRS = tuple((i,j) for i in range(3) for j in range(3) if i != j)
LABELS = (0,1,2,3,-1,-2,-3)


def require(test, message):
    if not test:
        raise ValueError(message)


def poly(terms):
    out = {}
    for e1,e2,e3,c in terms:
        e = (e1,e2,e3)
        out[e] = (out.get(e,0)+c) % P
    return tuple((*e,c) for e,c in sorted(out.items()) if c)


def padd(a,b):
    return poly(a+b)


def pneg(a):
    return tuple((x,y,z,(-c)%P) for x,y,z,c in a)


def pmul(a,b):
    if not a or not b:
        return Z
    if a == ONE:
        return b
    if b == ONE:
        return a
    return poly((x+u,y+v,z+w,c*d) for x,y,z,c in a for u,v,w,d in b)


def monomial(label):
    if label == 0:
        return ONE
    e = [0,0,0]
    e[abs(label)-1] = 1 if label > 0 else -1
    return ((*e,1),)


def elementary(root, coefficient):
    i,j = root
    a = list(IMAT)
    a[3*i+j] = coefficient
    return tuple(a)


def act_elementary(i,j,sign=1):
    a = list(IACT)
    a[3*i+j] = sign
    return tuple(a)


@lru_cache(maxsize=100000)
def mmul(a,b):
    if a == IMAT: return b
    if b == IMAT: return a
    return tuple(poly(term for k in range(3)
                      for term in pmul(a[3*i+k],b[3*k+j]))
                 for i in range(3) for j in range(3))


def mdet(a):
    out = Z
    for s,inds in ((1,(0,4,8)),(1,(1,5,6)),(1,(2,3,7)),
                   (-1,(2,4,6)),(-1,(1,3,8)),(-1,(0,5,7))):
        t=pmul(pmul(a[inds[0]],a[inds[1]]),a[inds[2]])
        out=padd(out,t if s==1 else pneg(t))
    return out


@lru_cache(maxsize=50000)
def minv(a):
    require(mdet(a)==ONE,"matrix inverse requested outside SL3")
    out=[]
    for i in range(3):
        for j in range(3):
            rows=[r for r in range(3) if r!=j]
            cols=[c for c in range(3) if c!=i]
            d=padd(pmul(a[3*rows[0]+cols[0]],a[3*rows[1]+cols[1]]),
                   pneg(pmul(a[3*rows[0]+cols[1]],a[3*rows[1]+cols[0]])))
            out.append(d if (i+j)%2==0 else pneg(d))
    return tuple(out)


@lru_cache(maxsize=10000)
def amul(a,b):
    return tuple(sum(a[3*i+k]*b[3*k+j] for k in range(3))
                 for i in range(3) for j in range(3))


def adet(a):
    return a[0]*(a[4]*a[8]-a[5]*a[7])-a[1]*(a[3]*a[8]-a[5]*a[6])+a[2]*(a[3]*a[7]-a[4]*a[6])


@lru_cache(maxsize=10000)
def ainv(a):
    require(adet(a)==1,"action inverse requested outside SL3Z")
    out=[]
    for i in range(3):
        for j in range(3):
            rr=[r for r in range(3) if r!=j]
            cc=[c for c in range(3) if c!=i]
            x=a[3*rr[0]+cc[0]]*a[3*rr[1]+cc[1]]-a[3*rr[0]+cc[1]]*a[3*rr[1]+cc[0]]
            out.append(x if (i+j)%2==0 else -x)
    return tuple(out)


@lru_cache(maxsize=100000)
def transform(a,m):
    if a==IACT or m==IMAT: return m
    result=[]
    for entry in m:
        result.append(poly((a[0]*x+a[1]*y+a[2]*z,
                            a[3]*x+a[4]*y+a[5]*z,
                            a[6]*x+a[7]*y+a[8]*z,c)
                           for x,y,z,c in entry))
    return tuple(result)


@lru_cache(maxsize=100000)
def gmul(g,h):
    if g==GE: return h
    if h==GE: return g
    return (mmul(g[0],transform(g[1],h[0])),amul(g[1],h[1]))


@lru_cache(maxsize=50000)
def ginv(g):
    ai=ainv(g[1])
    return (transform(ai,minv(g[0])),ai)


@lru_cache(maxsize=100000)
def in_h(g):
    return g[1]==IACT and all(x>=0 and y>=0 and z>=0 for p in g[0] for x,y,z,c in p)


def reduce_amalgam(syllables):
    """Reduced alternating syllables; uniqueness is NOT asserted.

    Identity is exactly the empty tuple.  H-membership is checked in the
    actual SL3 polynomial subgroup, not in its elementary subgroup.
    """
    out=[]
    for side,g in syllables:
        require(side in (0,1),"invalid side")
        if g==GE: continue
        if in_h(g): side=0
        while out:
            old_side,old=out[-1]
            if old_side!=side and not in_h(old) and not in_h(g):
                break
            out.pop()
            side=old_side if not in_h(old) else side
            g=gmul(old,g)
            if g==GE: break
            if in_h(g): side=0
        if g!=GE: out.append((side,g))
    return tuple(out)


def dmul(a,b):
    return reduce_amalgam(a+b)


def dinv(a):
    return reduce_amalgam((side,ginv(g)) for side,g in reversed(a))


def dequal(a,b):
    if a==b: return True
    return not dmul(dinv(a),b)


def winv(w): return tuple(-x for x in reversed(w))


def free(w):
    out=[]
    for x in w:
        require(isinstance(x,int) and x!=0,"invalid free-group letter")
        if out and out[-1]==-x: out.pop()
        else: out.append(x)
    return tuple(out)


def comm(u,v): return free(tuple(u)+tuple(v)+winv(u)+winv(v))


def make_presentation():
    """Return literal numbered generators, relators and quotient images."""
    gens=[]; lookup={}; images={}
    def add(key,name,side,g):
        k=len(gens)+1
        lookup[key]=k
        gens.append({'id':k,'name':name,'key':list(key)})
        images[k]=reduce_amalgam(((side,g),))
    # Twelve shared positive-root coefficient symbols.
    for f in range(3):
        for r in (0,1,2,3):
            add(('p',f,r),f'{"abc"[f]}{r}',0,(elementary(ROOTS[f],monomial(r)),IACT))
    # Each copy has nine negative coefficients and six action roots.
    for side in (0,1):
        for f in range(3):
            for r in (-1,-2,-3):
                add(('n',side,f,r),f'{"abc"[f]}M{-r}_{side}',side,(elementary(ROOTS[f],monomial(r)),IACT))
        for i,j in PAIRS:
            add(('t',side,i,j),f'T{i+1}{j+1}_{side}',side,(IMAT,act_elementary(i,j)))
    def fword(side,f,r):
        return (lookup[('p',f,r)] if r>=0 else lookup[('n',side,f,r)],)
    def tword(side,i,j): return (lookup[('t',side,i,j)],)
    def prodword(side,f,r,s):
        return comm(comm(fword(side,f,r),fword(side,(f+1)%3,0)),
                    comm(fword(side,(f+2)%3,s),fword(side,f,0)))
    records=[]; seen=set(); raw_counts={}
    def rel(family,w):
        w=free(w)
        if not w: return
        raw_counts[family]=raw_counts.get(family,0)+1
        if w not in seen:
            seen.add(w); records.append({'family':family,'word':list(w)})
    for side in (0,1):
        for f in range(3):
            for r in LABELS:
                rel('triangle-order',fword(side,f,r)*5)
            for r,s in combinations(LABELS,2):
                rel('triangle-abelian',comm(fword(side,f,r),fword(side,f,s)))
        for f in range(3):
            other=(f+1)%3
            for r,s,u in product(LABELS,repeat=3):
                c=comm(fword(side,f,r),fword(side,other,s))
                rel('triangle-class-two',comm(c,fword(side,f,u)))
                rel('triangle-class-two',comm(c,fword(side,other,u)))
        for (i,j),(k,l) in combinations(PAIRS,2):
            if i!=l and j!=k:
                rel('action-Steinberg-commute',comm(tword(side,i,j),tword(side,k,l)))
        for i,j,k in product(range(3),repeat=3):
            if len({i,j,k})==3:
                rel('action-Steinberg-product',comm(tword(side,i,j),tword(side,j,k))+winv(tword(side,i,k)))
        for i,j in PAIRS:
            for sign in (1,-1):
                t=tuple(sign*x for x in tword(side,i,j))
                for f in range(3):
                    for r in LABELS:
                        old=fword(side,f,r)
                        if abs(r)==j+1:
                            other_label=(i+1)*(sign if r>0 else -sign)
                            rhs=prodword(side,f,other_label,r)
                        else: rhs=old
                        rel('action-conjugation',t+old+winv(t)+winv(rhs))
        for i in range(3):
            k=(i+1)%3
            z=(tword(side,i,k)+winv(tword(side,k,i))+tword(side,i,k))*2
            for f in range(3):
                rel('negative-name',winv(fword(side,f,-i-1))+z+fword(side,f,i+1)+winv(z))
    h=fword(0,0,2); tl=tword(0,0,1); tr=tword(1,0,1)
    witness=comm(h,tr+winv(tl))
    # j is the forty-third generator.
    k=len(gens)+1
    gens.append({'id':k,'name':'j','key':['witness']})
    images[k]=eval_word(witness,images)
    rel('witness-name',(-k,)+witness)
    for k in tuple(images): images[-k]=dinv(images[k])
    return {'format':'qi218-presentation-v1','characteristic':P,'generators':gens,
            'relators':records,'raw_relation_counts':raw_counts,
            'witness_generator':k,'witness_word':list(witness),
            'shared_generators':list(range(1,13))},images


def eval_word(w,images):
    out=()
    for k in w:
        g=images[k] if k in images else dinv(images[-k])
        out=dmul(out,g)
    return out


def serial_syllables(x):
    return [{'side':side,'matrix':[[list(t) for t in entry] for entry in g[0]],
             'action':list(g[1])} for side,g in x]


def parse_syllables(x):
    return tuple((s['side'],(tuple(tuple(tuple(t) for t in p) for p in s['matrix']),tuple(s['action']))) for s in x)


class Labels:
    def __init__(self):
        self.values=[()]; self.words=[()]; self.simple={():0}; self.multi=[]
    def add(self,val,word):
        if len(val)<=1:
            if val in self.simple: return self.simple[val]
        else:
            for k in self.multi:
                if dequal(val,self.values[k]): return k
        k=len(self.values); self.values.append(val); self.words.append(tuple(word))
        if len(val)<=1: self.simple[val]=k
        else: self.multi.append(k)
        return k


def compile_gadget(presentation,images):
    labels=Labels(); generator_labels={}
    for g in presentation['generators']:
        for k in (g['id'],-g['id']):
            generator_labels[k]=labels.add(images[k],(k,))
    triangles=[]; tset=set(); occurrences=0; trivial=0
    def triangle(x,y,z):
        nonlocal occurrences,trivial
        occurrences+=1
        if x==0 or y==0:
            require(z==(y if x==0 else x),"trivial triangle disagrees")
            trivial+=1;return
        t=(x,y,z)
        if t not in tset: tset.add(t);triangles.append(t)
    for rec in presentation['relators']:
        w=tuple(rec['word']); old=0; value=(); prefix=[]
        for pos,k in enumerate(w):
            prefix.append(k); value=dmul(value,images[k])
            new=labels.add(value,prefix)
            if pos: triangle(old,generator_labels[k],new)
            old=new
        require(old==0,"nontrivial quotient relator: "+rec['family'])
    for g in presentation['generators']:
        k=g['id']; triangle(generator_labels[k],generator_labels[-k],0)
    L=len(labels.values); T=len(triangles); d=L+2*T
    kraus=[[[g,g,5]] for g in range(L)]
    for t,(x,y,z) in enumerate(triangles):
        p=L+2*t
        kraus[0].append([p,p,3])
        kraus[y].append([p,p+1,4])
        kraus[x].append([p+1,p,4])
        kraus[z].append([p+1,p+1,-3])
    payload={'format':'qi219-rational-gadget-v1','dimension':d,'kraus_rank':L,
             'entry_denominator':5,'kraus':kraus,'triangles':[list(t) for t in triangles],
             'labels':[{'id':k,'word':list(w),'normal_form':serial_syllables(v)}
                       for k,(w,v) in enumerate(zip(labels.words,labels.values))],
             'generator_labels':[[k,generator_labels[k]] for k in sorted(generator_labels)],
             'raw_triangle_occurrences':occurrences,'trivial_triangle_occurrences':trivial}
    return payload
