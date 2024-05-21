# IS_THIS_LCG?

- [attachment](https://github.com/hash-hash/crypto/tree/main/RCTF-2022/IS_THIS_LCG)

需要解决三个部分

1. 截断的LCG，仅已知高位
2. 椭圆曲线上的随机数生成器
3. $Z_p$下矩阵LCG

***part1***

$h_{n+1}+x_{n+1}=a(h_n+x_n)+b(mod~m)$

$⇒ax_n+(ah_n+b-h_{n+1})=x_{n+1}(mod~m)$

造hnp的格取规约出x0

```python
m = 2 ** 1024
a = int.from_bytes(b'Welcome to RCTF 2022', 'big')
b = int.from_bytes(b'IS_THIS_LCG?', 'big')

x = [0xc65f1c882be27b574c70f10e155ed3d3792d037d3c7
    , 0x142e1a26667e31a70eb58fa1e2b296d31a09675fa687
    , 0x17f366e283147917cc044778bbce2816884577126a9c
    , 0x2a316775dda35ad9a0e8a038757c85f216e91516f1ce
    , 0x3ef873ee8fa84fd071777521c78cb10a929f92f10dc7
    , 0x14e228828cb5090361501acac3108f05096fa8976e9c
    , 0x2e664838384824369607284ad9950f839f23a85c1974
    , 0x11affcbdf3da150c318bcc7096d21e8eb4bdaf904b9e]

A = []
B1 = []
B1.append(a * (x[0] << 850) + b)
B = []
for i in range(7):
    A.append(a ^ (i + 1) % m)
    B.append((B1[i] - (x[i + 1] << 850)) % m)
    B1.append((a * B1[i] + b) % m)

L = Matrix(ZZ, 9, 9)
for i in range(7):
    L[i, i] = m
    L[7, i] = A[i]
    L[8, i] = B[i]

L[7, 7] = 1
L[8, 8] = 2 ^ 850
basis = list(L.LLL())

for i in range(len(basis)):
    if int(basis[i][-1] // (2 ^ 850)) == 1:
        x0 = (x[0] << 850) + basis[i][-2]
        for i in range(8):
            x0 = (a * x0 + b) % m
        print(next_prime(x0))
```

***part2***

比赛直接找的论文

[Predicting the elliptic curve congruential generator](https://link.springer.com/article/10.1007/s00200-016-0303-x)

整体实现还挺麻烦的，只是题目只要实现第一步就行了，当时Hermes没看论文直接线性变换就出了，tql

$x_{i-1}=x(W_i⨁(-G))$

$x_{i+1}=x(W_i⨁G)$

$⇒x_{i-1}+x_{i+1}=2\frac{y_i^2+y_G^2}{(x_i-x_G)^2}-2(x_i+x_G)$

$⇒(2x_i^2+2x_i(x_{i-1}+x_{i+1}))x_G+(2x_i-(x_{i-1}+x_{i+1}))x_G^2+2x_iA+2(B+y_G^2-x_G^3)=(x_{i-1}+x_{i+1})x_i^2(mod~p)$

换成矩阵的形式即

$A*($ $x_G$   $x_G^2$   $A$   $B+y_G^2-x_G^3$ $)=B(mod~p)$

改写成

$(A|B)*($ $x_G$   $x_G^2$   $A$   $B+y_G^2-x_G^3$   $-1$ $)=0(mod~p)$

注意到A|B是方阵，由克莱姆法则知 $det(A|B)=0(mod~p)$ ，和N做gcd即可

实际上也不大需要论文这题，主要是比较懒不想推式子，直接找的最后结论抄了

```python
N = 0x614d9a106993a792c144715b0269a2726eb18a2e7b1ea7061bce1f6acb31af6289309d67ce6b28b3e88110c42785c0ca23833cc0e2aa4a30aadb16d25db7a74ef03b0898b7af47d56d4538b0f556b2779ed86e0600f821354d51f8551ccd23bbf8bf91eb9a9283a3d4d5248e3f404b4c6646a7dc805f29940a7e29d2f50343e1acc0d0067606606b331a64881bbafeafeb8ca44e736b41eab4608097216f587a1a4f74518614b46e91505e07c3a280b701ee88ca189e9903d601bc934584409d560027e5b34adb1f4949333ab5db34e95e49374e354d4ddc088855f1aae7a95e32ef195521b33f118169ae613e3fd5bf8d2942c2bde9ef506346698b0b5192c86b1efe24cffb907652afd5f0cb3966c7470195122ced63f5c40a4d9a3b6704e0b186ab7b9e3296b1299b6fa133d2455a8f8d8a9007a22bc61546b357ea314b0d369d72d22063c5ed6c14aa2a7edf31bdf93e63149818ef3724ca1cac367ac22b51260c793212ea221e062fcca68f28a4cd0b3bbeee03b9c73fd064c8298e775ab8a63c94db480a1eba918d09cba975304eed4fa5e874fc964e328547c23790e97102c6ad0bca9810dabb6285906f13d41798d3237333288b4498610d1a8fa79be85a522232a7cb904cd7c9b7fab995f39cd22a9758a5c2b6dcf44299df1e3e2ac360339b341ca6beb31eccba39ebd6f98dee127c6b5298db152fa6920b9703ab
x = [0x524456278d175edd6bcc3f2bbb8160a87dfe07092db7eedd1e4e3521e9cef7925e9c965a47ce9b7349456938fbf6d1d92095cfe7cdc06c8dbeac5284982d027179d8d363b1d1a9b95c2bb1334e589ac3c013d8cff1c904d0c2aed1f281e997be89abe3d0d2d668dc53adc4ae9870474a23ff993598bf2b51679179c8a1568619
,0x2f340fb1c6761e084b1465c5078f36e9caf7f9d6deecb969cc84fb5b85b1e4070157094c835333349f3d317e6a78a31a27d1ff0f8dfb103ead7444f26ac7b8be6b8ec346a8c8b4fe6f983db2729b6490ce0e1ea115b62f5e2888911d278153e3377a7456705c4f1a56588d8f727a91a8a401a852dd26573b2dc2ccf6a4af1de2
,0x92b2bbdf0c336be756ac47cc0b98fbc76b9ac679db96a5afd8fe500d16f4997503ee33d0508a59fff172042d6dcc4994a2d8220adcd8f5e591458b9409468c51b92dcacf73e793af3f793b9becb9cbb0704834861a43e1d1fd5cd5a9be14fafa8ec02df059fb3e1a3b0e7a8fb9969d42ffc13e2e3404fd539cd0d95b15f69f33
,0x795e49d6bc45ecf4d349d16058166f6422311344e3e6d8913a8b0a28225c92e203dfba92dd809a58a3630bfbb4cdf6d3118f172f6d6cb7c35cd4f9cf70947d091659e2ec4e248eaf2c456d58a149dd1fff7667630504cbb55cd82e3a2fe681f9b23de329d70a85f4badce87168dfb37b96b9edbeeb39a3d4ab28c130e9150140
,0x5a2bf69e31eef5ec1b990a2d2e3f8ccb08ba9996db2022775770b3b486909653b5347c15ceab62b167ad1dfa6a997efb56315fd6afde2e6c1b5af5a6e9b818556669992f148525c990bdb61e712339856dcf6e0f27ed8279bb32aba553bbab2ad3ac4accf3084638528a34434ec80df33705e381b39e9786593cff3e04a5b23d
,0x5baeb38339d662e8c16b1f16cd6129af38adebb264ffb197d6245f56df813c64b7ef28e60137b54d15a4227ca6ecd08f6ccbbbcb598bc94b1f326d8d488e13179d2999fb2c922165c9f27c2d7d0267e6924ce6395c33ec52a35776e88874877d8ccbbf8ca9ee214a7b73a8f7da23db02978f3b8bd145c2cca66b17638169f5e9
,0x4c2fc188b5cd7f4d19a4b120402946d7f8ca11c711e7771c39814e01c692160b7545edcff82a22d4634c416185eb58ff44adfdce5dc36a6d7c663f57eb19fc34f1a6c7e493518b094ad46fb8f9b6eb741c4666878ac91898116eb353a0a5aab9289322aaca6bed2ee104db17be339af54538635208f756da15bf46d18b0549a]

A = Matrix(ZZ,5,5)
for i in range(5):
    A[i,0]=2*x[i+1]^2+2*x[i+1]*(x[i]+x[i+2])
    A[i,1] = 2*x[i+1]-x[i]-x[i+2]
    A[i,2] = 2*x[i+1]
    A[i,3] = 2
    A[i,4] = (x[i]+x[i+2])*x[i+1]^2
    
p = GCD(A.det(),N)
print(p)
```

***part3***

矩阵实现加速，把通项求出来就行

```python
X=[0xc54aad8bd2b3233576847209ad1ade5f535622aab2a6279464832dea3dc88e7898a58130e36273143a90fcd4497079010e50658c2981e66e09ae86de089bf1f7123abb7d71fe68cf8d9eab3a2fc4792f1cb6444eff47c0f666995096c43ef8149fa78c061ca62809a2eadb00ac0dff81fb4163335c0a8014082e95b5007a2e2c
,0x3a3944bb3fd77217ab57358b174dbef9f704b844fb09f0d05bd4cfafc5a758f3b4d60c5cb584b1bb37f0c83bce8cba67bd04d11826433afba1717106da48a6cc22d571a0fa57fe63c29896783d6a9676f241cf4c9b1081aee364334ed3f80d680ad4c52d8a9e026fdfc97c1cda397a1f37c368420176e3270299efa21fa4c614
,0x309838246999c3a8920a9e8911f0c643eb614a9c522fb2cd5776bf582d7ad79796558b839e8ffc393e479aa0761d961df6860f9c44dea9b073a5006c2705128a7e7b139c407d15f430bd1a60d679d9f40deab664c84553fa8b9c1e8aeeb42e75c5c305d8b86e09debc9e193617f9fd619a0053017f71810cc3a48bb1fe89878
,0x38acf9569013ea3a32b18aef48ed6d0ad6557afe3e929c757d541039faefab0eeb53c5341a4ae5b9df610efcc66d09ae4238c569929d46409dd4f21d75a7bc97f3d8eed2dd397124d5a94946ccb8e8da8d030b4db4ac8821c313bdc87c8c25576050503891ca629b232e4f1b5c9bac4809979fc4dfb8f07260b3cdd62b2f45a6
,0xda4663505a3ce430f75fe908c34f96dfc8e3a997dcd378205274b1855804d069044558eeb09f0e36feeb34edd82ffec268095eef4acc795cacf4921bb33dab678f0ee930e7718839962511c49f91dddb4389cd9db61ba49baadd3a876952291d31b85b04ab2561a85542879d0e3287ad6f1c60b28daf05a56cab18955dd8d48a
,0x8e7250916637c65a685f5db8a3e5e84e223ebe59346f807048f16f5ee98ccf10679b3b1952e50ba32730906794d40c1aebdeccd059b775bce13186907ea883230160254254ebc4006a452826eb75361f92e5ae9b30f87e8c8abed2a90117eccf1b4e6aac455b1fc6a0983141dfe1df81b912612649e3bb48560eca66af9c9b76
,0x7ba1fb51f424a6257d85599cb596aa3bb0e83c94fa14ca716e5d933a507ba8cd1b6addc171d260ade722e01c7d69eaba0f5f3dfccccb2711b8407d0891e2179525577619f96735d55c98414f61042457059f93bb8613c81dd656885b4dbd5554a792c1e8226e0207ab3bae04e63bf5ab68190dc4915709e2eb2c6e3ddbf0b89a
,0xcf3b0d393c7ef1753e602e0b088fc15d0c06f949631cc9083ef7ab16c65148b47aa63eabb6151e39d85a5a339c065d9f1b4a33ab587f6093eb097fc6bab25a6b27cc8ae7d77775869b0864f6bdb7c1d8dcfa1a28dce4df346d95eaf90047020f4f8ad7e9496ed86e7c1bd840724348d88a308ca21174c61cf759ba106c548458
,0xe8e2ba97f5bdb287984e2a61b5a489ea4dc45c2c8e3601f151a20b92d1d6b7c0800712d07e4de5d2ca6f9cbfff25a64989e0779b98e56df1f4d8c301d3d743b86690d567c7f3a6bd74aa08b7df1970eb4b53ef2d5f8a7c3be585462dc3a972f99cd99b4ee1738a719476ebe70ba5a89447e020566e2a98ebab5747be0758a312
,0xcb1cd415bb82a8036035396806a37e28f23a709a51301fea6b0195e3da1a5ff7f71c6bd89387b955ff9e0d743f00e09286cf32520428791bac19368936f2e9bda4ffc4487a2bc999bb22249cffedc16dd686ac91d9a4cfe459e114ce38858f2b4972b09fd3c463c5b40cd553e640afe5803a390766842d2b6a74152923f329db]

n, m = 8, next_prime(2^16)

def num2mat(num):
    A=Matrix(Zmod(m),n,n)
    for i in range(n):
        for j in range(n):
            A[i,j]=num%m
            num=num//m
    return A

def mt2dec(X, n, m):
    x = 0
    for i in range(n):
        for j in range(n):
            x = x + int(X[i, j]) * (m ** (i * n + j))
    return x

S=[]
for i in range(10):
    S.append(num2mat(X[i]))

T=[]
for i in range(9):
    T.append(S[i]-S[i+1])
    
A=T[1]*T[0]^(-1)
B=S[1]-A*S[0]
X=S[0]
T=(A-matrix.identity(n))^(-1)*B
time=1337^1337
P=A^(time)*(S[0]+T)-T
print(next_prime(mt2dec(P,n,m)))
```

后面就简单的RSA解密了
