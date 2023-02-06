# S2DH

协议的整个流程和DH本质上没什么区别，就是换了个陷门的结构

其依赖于计算给定椭圆曲线之间同源的困难性，简单来说就是给你E和E的同态像 $E_{img}$ 理论上很难找到它们间的同源关系，由于我也是第一次见这种类型的密码协议所以对它困难程度的认知比较模糊。但是从最近的攻击来看这似乎与曲线的结构有关系，相关攻击利用了部分曲线已知的自同态结构，但是考虑到算法的效率问题是不大能将其替换为一般性的椭圆曲线。所以如何选取一个未知自同态环的 start curve 是当时初始版本的攻击出来后相关人员开始考虑修改的事情，但是几天时间又出了一篇对于一般性曲线的攻击，(xswl，求一波心理阴影面积)，至此为止 SIDH就被彻底打破了。看了一天论文，攻击倒是没怎么整明白但是一些工作者的吐槽倒是看乐了hhh。确实，一个自己研究的方向一夜之间被宣布打破内心也是挺麻木的了。

从相关的论文上可以看的出在这之前解决同源问题最优的算法复杂度是亚指数，经典算法复杂度是 $O(p^{\frac{1}{4}})$ ，量子算法是 $O(p^{\frac{1}{6}})$ ，而最新的Castryck Decru Attack则是给到了启发式多项式时间算法，可以说是基于数学手段上的巧妙且十分强大的攻击

稍微介绍一下SIDH的流程

1. 约定 $E_{start}:~E(F_{p^2})$ ，随机生成 $P_A,Q_A∈2^e$ -torsion group， $P_B,Q_B∈3^f$ -torsion group
2. A生成私钥 $r_A,s_A$ ，B生成私钥 $r_B,s_B$
3. A计算secrete kernel $G_A=r_AP_A+s_AQ_A$；B同理
4. A的公钥为 $\phi:E→E/G_A,\phi(P_B),\phi(Q_B)$；B同理
5. 最终的share key为 $(E/G_A)/\phi(G_B)$ 以及 $(E/G_B)/\psi(G_A)$ 的 j-invariant

算法的正确性很容易证明，由群同态定理我们可知上述曲线同构于 $\frac{E}{ker~f}$ ，由于他们映射的kernel是一致的所以这两曲线同构，而 j-invariant就是用来刻画同构关系的，所以最终结果一致

对于Castryck Decru Attack的那篇论文赛后又花一整天看了下，u1s1确实看不大明白，所以最后就匆匆浏览了一下攻击的主要思想，以我现在这数学水平死磕这东西也没啥意义，但这倒是足够解题了，因为这题就只是在现成的脚本上稍作修改

大致流程是基于97年提出的一个定理`Kani’s theorem`做猜测路径的判定条件，而在判断前要计算度为 $2^e-3^{f-k}$ 的辅助同源，而做这一步我们需要找初始曲线的一个满足要求的自同态 $2i$，后面通过KLPT算法可以很容易算得，额这个算法我也看，可惜wtcl完全看不懂

而这题就是要解决这个自同态的问题，当时改脚本的时候发现是这里自同态出了问题，然后还把论文2i这部分给漏掉了，搁这一直自己在造自同态人麻了，后面看和官方给的paper是同一篇，the auxiliary isogeny那块看的时候跳过了，悲):

- exp

```python
import public_values_aux
from public_values_aux import *

load('castryck_decru_shortcut.sage')

a = 43
b = 26

# Set the prime, finite fields and starting curve
# with known endomorphism
p = 2^a*3^b - 1
public_values_aux.p = p

Fp2.<i> = GF(p^2, modulus=x^2+1)
R.<x> = PolynomialRing(Fp2)

E = EllipticCurve(Fp2, [0,0,0,1,0])
E.set_order((p+1)^2) # Speeds things up in Sage

# Generation of the endomorphism 2i
#two_i = generate_distortion_map(E_start)

def two_i(A):
 x0,y0=A[0],A[1]
 x1=-x0
 y1=i*y0
 return 2*E(x1,y1)

# Generate public torsion points, for SIKE implementations
# these are fixed but to save loading in constants we can
# just generate them on the fly

P2 = E(20816113353953844596827139*i + 16418101434179547435831830, 9782287231195084940947894*i + 8305288838066432045414923)
Q2 = E(13022786448801065009926908*i + 21396754486749480260181021, 5027869541156315740937282*i + 8428382255806278677381816)
P3 = E(7582970089792232978539532*i + 6411668474015872447958400, 15459880436272725660545115*i + 7977012527121440514383975)
Q3 = E(10341548384598782389107676*i + 12525908271709247355078632, 6555843755802979256565190*i + 11595932163398809254591141)
Ea = EllipticCurve(Fp2, [4926878008530427712778566*i+8053083788709808436490360, 18771446501040649196825847*i+16306438728950797793375410])
Eb = EllipticCurve(Fp2, [18866222948911535725014127*i+21372353382532165741892023, 14780329017962693588095579*i+4731720677310255642021851])
Pb = Ea(2535790352220803985875373*i + 17699033710915047849396921, 2413558249712558899689063*i + 5157954648088691506046995)
Qb = Ea(16568070039544280994803013*i + 21423138055383385576701886, 5040448698696125071219900*i + 6672798507142407841550817)
Pa = Eb(3413055427164626562463192*i + 5176875496413372729075617, 17919859745180152815219510*i + 18120119720358642060676362)
Qa = Eb(18433160961475396600407402*i + 22312166252239187097449810, 10433258275941991434154560*i + 9029292514862239326241711)
enc = 243706092945144760206191226817331300960683091878992

'''
# Generate Bob's key pair
bob_private_key, EB, PB, QB = gen_bob_keypair(E_start, b, P2, Q2, P3, Q3)
solution = Integer(bob_private_key).digits(base=3)
'''

print(f"Running the attack against Baby SIDHp64 parameters, which has a prime: 2^{a}*3^{b} - 1")
#print(f"If all goes well then the following digits should be found: {solution}")

# ===================================
# =====  ATTACK  ====================
# ===================================


if __name__ == '__main__' and '__file__' in globals():
    recovered_key = CastryckDecruAttack(E, P2, Q2, Eb, Pa, Qa, two_i)
    K = Pb + recovered_key*Qb
    J = Ea.isogeny(K, algorithm='factored').codomain().j_invariant()
    m=enc^^((int(J[1]) << 84) + int(J[0]))
    print(bytes.fromhex(hex(m)[2:]))
    
# b'SIDH_isBr0ken_in_2O22'
```

**[Reference]**

[Supersingular isogeny key exchange](https://en.wikipedia.org/wiki/Supersingular_isogeny_key_exchange)

[AN EFFICIENT KEY RECOVERY ATTACK ON SIDH](https://eprint.iacr.org/2022/975.pdf)

[Breaking SIDH in polynomial time](https://eprint.iacr.org/2022/1038)

[A SageMath implementation of the Castryck-Decru Key Recovery attack on SIDH](https://github.com/jack4818/Castryck-Decru-SageMath)