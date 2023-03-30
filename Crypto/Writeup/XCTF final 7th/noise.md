```python
def gen(n, m, r, N):
    t1 = [ZZ.random_element(-2^15, 2^15) for _ in range(n*m)]
    t2 = [ZZ.random_element(N) for _ in range(r*n)]
    t3 = [ZZ.random_element(-2^20, 2^20) for _ in range(r*m)]
    B = matrix(ZZ, n, m, t1)
    L = IntegerLattice(B)
    A = matrix(ZZ, r, n, t2)
    X = matrix(ZZ, r, m, t3)
    C = (A * B + X) % N
    return L, C
```

HSSP基础上加了error，实际上没有任何影响，并且题目给了56组数据只需要8组即可利用`orthogonal lattice`还原B，全取的话规约速度极慢，八组7分钟就能出

$\vec{c}=\sum a_i·\vec{b_i}+\vec{e}(mod~N)$

$=>\vec{c}·\vec{v}=\sum a_i·\vec{b_i}·\vec{v}+\vec{e}·\vec{v}(mod~N)$

只需 $(a_1, ...,a_m,1)$ 正交格的最短向量量级在 $2^{20}$ 以上即可

取 $v$ 满足 $\vec{c}·\vec{v}=0(mod~N)$

有 $\vec{b_i}·\vec{v}=0$ 且 $\vec{e}·\vec{v}=0$

满足条件 $\vec{v}$ 的组数是m-n-r(m变量，n+r个方程)

故 $b_i$ 在 $\vec{v}$ 所构成的正交格上，取V的右核，做LLL约化

做的时候发现key在约化基的第一行或第二行，都试试就行了

```python
import time
import re
from Crypto.Cipher import AES
from hashlib import sha256

n = 75
m = 150
r = 8
N = 126633165554229521438977290762059361297987250739820462036000284719563379254544315991201997343356439034674007770120263341747898897565056619503383631412169301973302667340133958109


f=open('output.txt','r')

def s2row(s):
    return list(map(int, re.split(r'\s+', s[1:-1].strip())))

matA=[]
for _ in range(r):
    s = f.readline().strip()
    row = s2row(s)
    matA.append(row)


def orthogonal_lattice_N(matA):
    g=2^2048
    L=Matrix(ZZ,m+r,m+r)
    for i in range(m):
        L[i,i]=1
    for i in range(r):
        for j in range(m):
            L[j,m+i]=g*matA[i][j]
        L[m+i,m+i]=g*N
    return L.LLL()[:m-n-r]


C = Matrix(matA)

start = time.time()
LL=list(orthogonal_lattice_N(C))
print("LLL done")

for i in range(len(LL)):
    LL[i]=LL[i][:-r]


LL=Matrix(ZZ,LL)
sk=LL.right_kernel(basis="LLL").basis_matrix()[0]

ct = b'\xf1\xd5\xddA\xd8\xd7dK\xc2U\xa82\xa6\xbf\xea10\xecWA\xe6w\x14\x8d\x1b\x8c?i\x02\xc2\xb5\xda'
key = sha256(str(sk).encode()).digest()
aes = AES.new(key, AES.MODE_ECB)
pt = aes.decrypt(ct)

end = time.time()
print("cost {}".format(end-start))
print(pt)
```

