# hardrsa

一种针对unbalanced CRT-RSA的攻击方式

翻了下N1之前的文档，2021 0CTF出过一道ezRSA+就是CRT-RSA with small dp,dq的攻击，这里的unbalanced的条件一定程度扩大了dp的范围，最新的改进里unbalanced的程度已经很接近正常取值了，只需要 $p < N^{0.468}$ 即可满足

相关参数 $e=N^\alpha,p=N^\beta,d_q=N^\delta$

```python
import math
alpha = 1
beta = 0.233
delta = 0.226
check=((1-beta)*(3+2*beta)-2*math.sqrt(beta*(1-beta)*(alpha*beta+3*alpha+beta)))/(3+beta)-delta
print(check)
#0.10940579921625301
```

由条件有 $ed_q=1+k(q-1)$

构造多项式 $f_q(x_q,y_q)=1+x_q(y_q-1)=0(mod~e)$

两边同乘p我们有 $ed_pp=p+k(N-p)=N+(k-1)(N-p)$

构造多项式 $f_p(x_p,y_p)=N+x_p(N-y_p)=0(mod~e)$

各变量的Bound

$X_p=N^{\alpha+\beta+\delta-1},X_q=N^{\alpha+\beta+\delta-1},Y_p=N^{\beta},Y_q=N^{1-\beta}$

计算 $\tau_p=\frac{1-2\beta-\delta}{2\beta},\tau_q=\frac{1-\beta-\delta}{1-\beta}$

```python
beta = 0.233
delta = 0.226
tp=(1-2*beta-delta)/(2*beta)
tq=(1-beta-delta)/(1-beta)
print(tp,tq)
```

构造划分

$S_1=\{i=0,...,m;j=0,...,m-i\}$

$S_2=\{i=0,...,m;j=1,...,\lceil\tau_pm\rceil\}$

$S_3=\{i=1,...,m;j=1,...,\lceil\tau_qi\rceil\}$

shift多项式定义

$g_{1[i,j]}(x_p,y_p)=x_p^jf_p^i(x_p,y_p)e^{m-i}$

$g_{2[i,j]}(x_p,y_p)=y_p^jf_p^i(x_p,y_p)e^{m-i}$

$g_{3[i,j]}(x_p,x_q,y_p,y_q)=f_p^{i-j}(x_p,y_p)f_q^j(x_q,y_q)e^{m-i}$

最后实际上我们可以四元简化成两元来降低求根的难度，做次代换即可

官方给的附件的参数给错了，看他wp写的0.233附件里给的0.223，然后估计是卡了下界，我用0.223的话求不出来

- exp

```python
from Crypto.Util.number import *

def getpoly(coe,mon):
    h=0
    for idx,monomial in enumerate(mon):
        h+=coe[idx]//monomial(Xp,Yp)*monomial
    return h

N=17898692915537057253027340409848777379525990043216176404521845629792286203459681133425615460580210961931628383718238208402935069434512008997422795968676635886265184398587211149645171148089263697198308448184434844310802022336492929706736607458307830462086477073132852687216229392067680807130235274969547247389
e=7545551830675702855400776411651853827548700298411139797799936263327967930532764763078562198672340967918924251144028131553633521880515798926665667805615808959981427173796925381781834763784420392535231864547193756385722359555841096299828227134582178397639173696868619386281360614726834658925680670513451826507
c=2031772143331409088299956894510278261053644235222590973258899052809440238898925631603059180509792948749674390704473123551866909789808259315538758248037806795516031585011977710042943997673076463232373915245996143847637737207984866535157697240588441997103830717158181959653034344529914097609427019775229834115

alpha = 1
beta = 0.233
delta = 0.226
Xp=int(RR(N)^(alpha+beta+delta-1))
Yp=int(RR(N)^beta)


tp=(1-2*beta-delta)/(2*beta)
tq=(1-beta-delta)/(1-beta)
m=5
PR.<xp,yp>=PolynomialRing(ZZ)
fq=yp+(xp+1)*(N-yp)
fp=N+xp*(N-yp)

poly=[]
monomials=set()
S1=[(i,j) for i in range(m+1) for j in range(m-i+1)]
S2=[(i,j) for i in range(m+1) for j in range(1,int(tp*m)+1)]
S3=[(i,j) for i in range(1,m+1) for j in range(1,int(tq*i)+1)]

for item in S1:
    i,j=item
    g=xp^j*fp^i*e^(m-i)
    poly.append(g)
    for mono in g.monomials():
        monomials.add(mono)
    
for item in S2:
    i,j=item
    g=yp^j*fp^i*e^(m-i)
    poly.append(g)
    for mono in g.monomials():
        monomials.add(mono)
    
for item in S3:
    i,j=item
    g=fp^(i-j)*fq^j*e^(m-i)
    poly.append(g)
    for mono in g.monomials():
        monomials.add(mono)
        
        
A=Matrix(ZZ,len(poly),len(monomials))

for raw,shift in enumerate(poly):
    for col,monomial in enumerate(monomials):
        A[raw,col]=shift.monomial_coefficient(monomial)*monomial(Xp,Yp)
        
basis=A.LLL()

coe_f1=basis[8]
coe_f2=basis[9]

f1=getpoly(coe_f1,monomials).change_ring(ZZ)
f2=getpoly(coe_f2,monomials).change_ring(ZZ)

load('poly_gcd.py') 
pr = f1.parent() 
for roots in find_roots(pr, [f1,f2], method="resultants"): 
    p=tuple(roots[xi] for xi in pr.gens())[1]
    break
    
q=N//p
phi=(p-1)*(q-1)
d=inverse_mod(e,phi)
m=int(pow(c,d,N))
print(long_to_bytes(m))
```

**[Reference]**

[Small CRT-Exponent RSA Revisited](https://www.iacr.org/archive/eurocrypt2017/10210359/10210359.pdf)
