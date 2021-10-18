class Vector():
  def __init__(self, vec):
    self._vec = list(vec)
    self._NF=[]
    self._const=[]
  def vec(self):
    return tuple(self._vec)
  def NF(self):
    return '*'.join(['sqrt(n__%s)' % (VECTOR_INDEX[i] + num2str(j)) for i,j in self._NF])
  def const(self):
    if len(self._const)>0:
      return '*'.join(['(%s)' % (i) for i in self._const if len(i)>0])
    else:return ''

def Creation(vec,i):
  while len(vec._vec)<i+1:
      vec._vec.append(0)
  vec._vec[i] +=1
  vec._NF.append((i,vec._vec[i]))

def Annihilation(vec,i):
  while len(vec._vec)<i+1:
      vec._vec.append(0)
  vec._NF.append((i,vec._vec[i]))
  vec._vec[i] -=1

def ksi_polinom(n:int)->list:
  ksi = [Annihilation,Creation]
  return [x[::-1] for x in product(ksi,repeat=n)]  

def CALC(bra,fksi,ket,ijk):
    const=''
    NF=''
    if isinstance(ket, Vector):
        const+=ket.const()
        NF+=ket.NF()
        ket=ket._vec
    if isinstance(bra, Vector):
        const+='*'+bra.const()
        NF+='*'+bra.NF()
        bra=bra._vec
    if isinstance(bra, tuple):bra=list(bra)
    if isinstance(ket, tuple):ket=list(ket)
    result = [Vector(ket) for i in range(len(fksi))]
    for vec,operations in zip(result,fksi):
        for op,i in zip(operations,ijk):
            op(vec,i)
    RESULT=[]
    for i in result:
        if bra==i._vec: RESULT.append('*'.join([c for c in ['1',i.NF(),const,i.const(),NF,] if c!='']))
    return '(%s)'% ('+'.join([i for i in RESULT if i!=''])) if len(RESULT)>0 else 0

from collections import namedtuple
Index = namedtuple('Index','p q betta gamma nu')
Index2 = namedtuple('Index2','p betta gamma')

def AV(ket=[0,0,0,0],ind=0):
  if ind==0:return [[tuple(ket),'']]
  N=[]
  for i in [Index(*i) for i in product(range(0,ind+1), repeat=5) if sum(i)==ind and i[1]%2==0 and i[0]>0]:
    KET=AV(ket,i.gamma)
    fksi=ksi_polinom(i.p+2)
    NF1='(%s/%s)'%(i.p,ind)
    for bra in VECTORS_m(3*(i.betta+i.gamma)+i.p+2,ket):
      if list(bra)==list(ket):continue
      if SelectionRule(bra,ket,i.p+2,i.betta,i.gamma)==0:continue
      DELTA=delta(bra,ket,i.q)
      if DELTA==0:continue
      BRA=AV(bra,i.betta)
      for ijk in INDEX(i.p,ket):
        A='(A__%s)'%(ijk)
        IJK =[ VECTOR_INDEX_MAP[x] for x in ijk]
        NF2=calc(BRA,KET,fksi,IJK,'s')
        if NF2==0:continue
        const='*'.join([k for k in [NF1,A,DELTA,NF2]])
        PP=[[V[0],'*'.join([k for k in [const,V[1]]if k!=''])] for V in AV(bra,i.nu)]
        N+=PP
  return N

def calc(BRA,KET,fksi,IJK,flag='s'):
  result=[]
  for i in BRA:
    for j in KET:
      A=CALC(i[0],fksi,j[0],IJK)
      if A==0:continue
      I='*'.join([k for k in [A,i[1],j[1]] if k!=''])
      if len(I)>0:
        I='(%s)' %(sy.simplify(I))
        result.append(I)
  if flag=='s': return '(%s)'%(sy.simplify('+'.join([k for k in result if k!='(0)']))) if len(result)>0 else 0
  if flag=='l': return [k for k in result if k!='()']

def AE(vec1=[0,0,0,0],vec2=[0,0,0,0],ind=0,flag='s'):
  if ind==0:
    E0=[]
    fksi=ksi_polinom(2)
    for ijk in INDEX(0,vec1):
      IJK =[ VECTOR_INDEX_MAP[x] for x in ijk]
      E=sy.simplify('0.5*(%s)*omega__%s'%(CALC(vec1,fksi,vec2,IJK),ijk[0]))
      if E!=0:
        E0.append(str(E))
    if len(E0)>0:
      return '+'.join(E0)
    else: return '0'
  if ind%2==1 and vec1==vec2: return '0'
  N=[]
  X=[]
  for i in [Index2(*i) for i in product(range(0,ind+1), repeat=3) if sum(i)==ind and i[0]>0]:#(p,бетта,гамма)
    KET=AV(vec1,i.gamma)
    BRA=AV(vec2,i.betta)
    fksi=ksi_polinom(i.p+2)
    NF1='(%s/%s)'%(i.p,ind)
    E_ijk=[]
    for ijk in INDEX(i.p,vec1):
      A='(A__%s)'%(ijk)
      IJK =[ VECTOR_INDEX_MAP[x] for x in ijk]
      E=calc(BRA,KET,fksi,IJK,'s')
      if E==0:continue
      E_ijk.append(str(sy.simplify(E+'*'+A)))
     # E_ijk.append('(%s)*%s'%(sy.simplify(E),A))
    if len(E_ijk)>0:
      E_ijk='+'.join(E_ijk)
      N.append('%s*(%s)'%(NF1,E_ijk))
  if len(N)==0: return '0'
  elif flag=='s': return '+'.join(N)
  elif flag=='l': return N

def SelectionRule(bra,ket,nksi,alfa,betta):
    k = sum(bra)-sum(ket)
    return 0 if (alfa+betta+nksi) % 2 != k % 2 else 1

def delta(bra,ket,q):
    if q==0:
        omega = ['(%s*omega__%s)' % (r-l,VECTOR_INDEX[i])  for i,(l,r) in enumerate(zip(bra,ket)) if l!=r]
        return '(1/(%s))' % '+'.join(omega) if omega else 0
    if q==2:
        A=delta(bra,ket,0)
        B='-'.join([AE(ket,ket,2,'s'),AE(bra,bra,2,'s')])
        return '*'.join([A,A,B]) if A else 0
    else: return 0
from itertools import product, combinations_with_replacement
import numpy as np

def num2str(n):
  return '%+d' % n if n else ''

VECTOR_INDEX = 'ijklmnopq'
VECTOR_INDEX_MAP = {k:i for i,k in enumerate(VECTOR_INDEX)}

def INDEX(n:int,vec)->list:
  return [''.join(sorted(i)) for i in combinations_with_replacement(VECTOR_INDEX[:len(vec)], n+2) ]    

def VECTORS_m(n:int,vec):
  rng = range(-n,n+1)
  el = [i for i in product(rng, repeat=len(vec)) if sum(map(abs,i)) <=n]
  el = np.array(el) + vec
  return list(map(tuple,el))

import sympy as sy
A__iiii, A__iiij, A__iiik, A__iiil, A__iijj, A__iijk, A__iijl, A__iikk, A__iikl, A__iill, A__ijjj, A__ijjk, A__ijjl, A__ijkk, A__ijkl, A__ijll, A__ikkk, A__ikkl, A__ikll, A__illl, A__jjjj, A__jjjk, A__jjjl, A__jjkk, A__jjkl, A__jjll, A__jkkk, A__jkkl, A__jkll, A__jlll, A__kkkk, A__kkkl, A__kkll, A__klll, A__llll = sy.symbols('A__iiii, A__iiij, A__iiik, A__iiil, A__iijj, A__iijk, A__iijl, A__iikk, A__iikl, A__iill, A__ijjj, A__ijjk, A__ijjl, A__ijkk, A__ijkl, A__ijll, A__ikkk, A__ikkl, A__ikll, A__illl, A__jjjj, A__jjjk, A__jjjl, A__jjkk, A__jjkl, A__jjll, A__jkkk, A__jkkl, A__jkll, A__jlll, A__kkkk, A__kkkl, A__kkll, A__klll, A__llll')
A__iii, A__iij, A__iik, A__iil, A__ijj, A__ijk, A__ijl, A__ikk, A__ikl, A__ill, A__jjj, A__jjk, A__jjl, A__jkk, A__jkl, A__jll, A__kkk, A__kkl, A__kll, A__lll = sy.symbols('A__iii, A__iij, A__iik, A__iil, A__ijj, A__ijk, A__ijl, A__ikk, A__ikl, A__ill, A__jjj, A__jjk, A__jjl, A__jkk, A__jkl, A__jll, A__kkk, A__kkl, A__kll, A__lll')
n__i, n__j, n__k, n__l=sy.symbols('n__i, n__j, n__k, n__l')
omega__i, omega__j, omega__k, omega__l=sy.symbols('omega__i, omega__j, omega__k, omega__l')

D__iiii, D__iiij, D__iiik, D__iiil, D__iijj, D__iijk, D__iijl, D__iikk, D__iikl, D__iill, D__ijjj, D__ijjk, D__ijjl, D__ijkk, D__ijkl, D__ijll, D__ikkk, D__ikkl, D__ikll, D__illl, D__jjjj, D__jjjk, D__jjjl, D__jjkk, D__jjkl, D__jjll, D__jkkk, D__jkkl, D__jkll, D__jlll, D__kkkk, D__kkkl, D__kkll, D__klll, D__llll = sy.symbols('D__iiii, D__iiij, D__iiik, D__iiil, D__iijj, D__iijk, D__iijl, D__iikk, D__iikl, D__iill, D__ijjj, D__ijjk, D__ijjl, D__ijkk, D__ijkl, D__ijll, D__ikkk, D__ikkl, D__ikll, D__illl, D__jjjj, D__jjjk, D__jjjl, D__jjkk, D__jjkl, D__jjll, D__jkkk, D__jkkl, D__jkll, D__jlll, D__kkkk, D__kkkl, D__kkll, D__klll, D__llll')
D__iii, D__iij, D__iik, D__iil, D__ijj, D__ijk, D__ijl, D__ikk, D__ikl, D__ill, D__jjj, D__jjk, D__jjl, D__jkk, D__jkl, D__jll, D__kkk, D__kkl, D__kll, D__lll = sy.symbols('D__iii, D__iij, D__iik, D__iil, D__ijj, D__ijk, D__ijl, D__ikk, D__ikl, D__ill, D__jjj, D__jjk, D__jjl, D__jkk, D__jkl, D__jll, D__kkk, D__kkl, D__kll, D__lll')

def Energy_Resonance(vec1,vec2,zamena):
  vec=[0 for i in range(len(vec1))]
  E0=[AE(vec,vec,0,'s'),AE(vec,vec,2,'s')]
    
  AA=[AE(vec1,vec1,0,'s'),AE(vec1,vec1,1,'s'),AE(vec1,vec1,2,'s')]
  AB=[AE(vec1,vec2,0,'s'),AE(vec1,vec2,1,'s'),AE(vec1,vec2,2,'s')]
  BA=[AE(vec2,vec1,0,'s'),AE(vec2,vec1,1,'s'),AE(vec2,vec1,2,'s')]
  BB=[AE(vec2,vec2,0,'s'),AE(vec2,vec2,1,'s'),AE(vec2,vec2,2,'s')]

  AA1=[sy.sympify(i).subs(zamena) for i in AA]
  AB1=[sy.sympify(i).subs(zamena) for i in AB]
  BA1=[sy.sympify(i).subs(zamena) for i in BA]
  BB1=[sy.sympify(i).subs(zamena) for i in BB]
  E0=[sy.sympify(i).subs(zamena) for i in E0]

  A=sum(AA1)-sum(E0)
  B=sum(AB1)
  C=sum(BA1)
  D=sum(BB1)-sum(E0)

  E1=(A+D+((A+D)**2-4*(A*D-C*B))**0.5)/2
  E2=(A+D-((A+D)**2-4*(A*D-C*B))**0.5)/2
  
  return E1,E2

