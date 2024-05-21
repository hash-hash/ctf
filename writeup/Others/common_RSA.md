# common_RSA

-  [attachment](https://github.com/hash-hash/crypto/blob/main/blog附件/祥云杯-common_RSA.py)
-  [paper](https://www.sciencedirect.com/science/article/pii/S0020025520305004#s0030)

CP-RSA系统的攻击，论文里面提到的方式还不少，先看下参数

$e≈N^{1-γ},d=N^δ,g=N^γ$

论文最新给出的攻击方案给的界

$δ<2-γ-\frac{1}{4}\sqrt{4γ^2-28γ+37}-ε$

可惜看不大懂他后面一些变量的描述，比较迷糊

这里我们用`Jochemsz-May’s asymptotic bounds`

$\delta<\frac{1}{4}(4+4γ-\sqrt{13+20γ+4γ^2})$

```python
import math
gama=0.314
delta=0.135
check= 1/4*(4+4*gama-math.sqrt(13+20*gama+4*gama**2))-delta
print(check)
#0.07010460367084215
```

由条件得到等式

$e^2d^2+e(kc_1+kc_2-2)d-kc_1kc_2(N-1)-kc_1-kc_2+1=0$

构造多项式

$f(x,y,z)=e^2x^2+ex(y+z-2)-(y+z-1)-(N-1)yz$

满足 $f(d,kc_1,kc_2)=0$

根据参数有 $X=N^δ,Y=N^{δ+\frac{1}{2}-γ},Z=N^{δ+\frac{1}{2}-γ}$

选择参数m，这个有关格的维度

模数 $R=WX^{max}Y^{max}Z^{max}$

M,S分别为 $f^ m,f^ {m-1}$ 的单项式集合

$g_{i,j,k}=x^iy^jz^kX^{2(m-1)+t-i}Y^{m-1-j}Z^{m-1-k}f(x,y,z)$   $x^iy^jz^k∈S$

$h_{i,j,k}=Rx^iy^jz^k$   $x^iy^jz^k∈M/S$

该题用Jochemsz-May’s asymptotic bounds是够的，因为是道陈题的缘故，没对参数做修改，所以没有用到这篇论文的改进结果

最后我本来是直接取的 $f,f_ 1,f_ 2$ 的理想的Gröbner基，结果跑了一小时还没整完，就换成gcd的板子了，感觉这个地方都差不多就没仔细研究了

-  exp:

```python
from Crypto.Util.number import *

def getpoly(coe,mon):
    h=0
    for idx,monomial in enumerate(mon):
        h+=coe[idx]//monomial(X,Y,Z)*monomial
    return h

N = 253784908428481171520644795825628119823506176672683456544539675613895749357067944465796492899363087465652749951069021248729871498716450122759675266109104893465718371075137027806815473672093804600537277140261127375373193053173163711234309619016940818893190549811778822641165586070952778825226669497115448984409
e = 31406775715899560162787869974700016947595840438708247549520794775013609818293759112173738791912355029131497095419469938722402909767606953171285102663874040755958087885460234337741136082351825063419747360169129165
c = 97724073843199563126299138557100062208119309614175354104566795999878855851589393774478499956448658027850289531621583268783154684298592331328032682316868391120285515076911892737051842116394165423670275422243894220422196193336551382986699759756232962573336291032572968060586136317901595414796229127047082707519

gama=0.314
delta=0.135
m=2
t=0
X=int(pow(N,delta))
Y=int(pow(N,(delta+1/2-gama)))
Z=int(pow(N,(delta+1/2-gama)))

R=(N-1)*Y*Z*X^2*Y*Z+1
PR.<x,y,z>=PolynomialRing(ZZ)
f=e^2*x^2+e*x*(y+z-2)-(y+z-1)-(N-1)*y*z

S,M=(f^(m-1)).monomials(),(f^m).monomials()

poly=[]
monomials=set()
for item in S:
    g=item*f*X^(2-item.degree(x))*Y^(1-item.degree(y))*Z^(1-item.degree(z))
    poly.append(g)
    monomials.add(item)
    
for item in M:
    if item not in S:
        poly.append(R*item)
        monomials.add(item)
        
A=Matrix(ZZ,len(poly),len(monomials))

for raw,shift in enumerate(poly):
    for col,monomial in enumerate(monomials):
        A[raw,col]=shift.monomial_coefficient(monomial)*monomial(X,Y,Z)
        

basis=A.LLL()
coe_f1=basis[0]
coe_f2=basis[1]
f1=getpoly(coe_f1,monomials).change_ring(ZZ)
f2=getpoly(coe_f2,monomials).change_ring(ZZ)
f=f.change_ring(ZZ)

load('poly_gcd.py')
pr = f.parent()
for roots in find_roots(pr, [f,f1,f2], method="resultants"):
    d=tuple(roots[xi] for xi in pr.gens())[0]
    break
    
print(long_to_bytes(int(pow(c,d,N))))
```

