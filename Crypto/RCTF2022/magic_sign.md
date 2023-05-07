# magic_sign

- [attachment](https://github.com/hash-hash/crypto/tree/main/RCTF-2022/magic_sign)

虽然该题不算太难但我觉得确实是一道很不错的题目，后来看出题人的wp为了避免纯论文题还是特地改的局部，让它存在思维上更加简单的攻击方式

该题的方案存在一个十分脆弱的点，就是它的随机数发生器参数有问题，这个点出题人有暗示，所以特地试了下，发现参数问题大大降低了枚举的可能情况

```python
magic = Magic(137)          # most magical number
C = magic(C)
U = C.U
print(set(U))
#{4, 133, 135, 11, 78, 81, 18, 24, 26}
```

主要分析一下mix函数

```python
def mix(self, T, K):
 R = T+K
 e = self.U[0]
 for i in range(self.N ** 2):
  d = self.U[i+1]
  R.lst[d] = R.lst[d] + R.lst[e]
  e = d
 R = R+K
 return R
```

由于参数的问题U中的元素非常少，也就是 $T+K$ 中的大部分分量是不经过循环部分的，也就是这些分量只与T，K的对应位置有关系，所以在已知T和最终R的情况下我们枚举K的对应分量就行

而对于U中元素的部分我们采取穷搜的话需要 $2^{3*9}$ 次循环

但是即使这样不做任何改进直接穷搜的话还是不行，连上远程后等太长时间会断开，所以在搜索过程中加一步剪枝能稍微缩短一点时间，可能还是会连接断开，多试几次

官方wp中给了另外的一个方向，这个是通用的攻击方式不依赖于此处参数的漏洞，也是打破Xifrat0的攻击，利用了mix的一个特殊的结构，我们只需要枚举 $137*2^6$ 次即可，具体的可以参考官方wp给的那篇论文

- exp

```python
import itertools
from magic import *

C = '74032535155656707612174202154246122575411201740420517665613466343607106174320752346310145452517224771057452224514504103626243313340500141'
P1 = '26034574654564722553425220461341403675566726063433617616610164212572404636322263301253445543417275731360160736546535670526342042043423431'
P2 = '33623770257737645124556674275336430635672552177546635271132203077714551753677451066527734062615064741367336234266361430617305106626521751'
magic = Magic(137)
C = magic(C)
P1 = magic(P1)
P2 = magic(P2)
H = magic.shake(b'Never gonna give you flag~')
res = C*H*P2

N = 137
U = C.U
print(set(U))

'''
print(set(U))
{4, 133, 135, 11, 78, 81, 18, 24, 26}
'''

key = [11, 135, 81, 26, 4, 24, 18, 78, 133]

u = {}
for i in range(len(key)):
    u[key[i]] = 0

for i in range(len(U)):
    if U[i] in u.keys():
        u[U[i]] = i

v = [0]+list(u.values())
fake_S = [0]*137


for i in range(137):
    if i not in key:
        t = int(str(P1.lst[i]))
        s = int(str(res.lst[i]))
        for k in range(8):
            if Magic_Latin_Square[Magic_Latin_Square[t][k]][k] == s:
                fake_S[i] = k
                break


def check(item):
    item2 = [0]*9
    for i in range(9):
        item2[i] = Magic_Latin_Square[h[i]][item[i]]
    e = UU[0]
    for i in range(9):
        for j in range(v[i], v[i+1]):
            d = UU[j + 1]
            item2[d] = Magic_Latin_Square[item2[d]][item2[e]]
            e = d
        if Magic_Latin_Square[item2[d]][item[d]] != sign[d]:
            return False
    for i in range(9):
        item2[i] = Magic_Latin_Square[item2[i]][item[i]]
    print(item2)
    if item2 == sign:
        return True


sol = range(8)
h = [0]*9
sign = [0]*9
for i in range(9):
    h[i] = int(str(P1.lst[key[i]]))
    sign[i] = int(str(res.lst[key[i]]))

UU = []
for i in U:
    UU.append(key.index(i))

S = magic(fake_S)

print("start search")
for ele in itertools.product(sol, sol, sol, sol, sol, sol, sol, sol, sol):
    if check(ele):
        for i in range(9):
            fake_S[key[i]] = ele[i]
        S = magic(fake_S)
        print(ele)
        print(S)
        print(P1*S)
        print(res)
        break
```

**[Reference]**

[Xifrat Cryptanalysis - Compute the Mixing Function Without the Key](https://eprint.iacr.org/2021/487)