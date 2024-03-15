# Numerical Algebra

## Ch1 线性方程组的直接解法(Ax=b)

### 1.1 三角形方程的解法

- 下三角形——前代法

  - 因为第一行只有$a_{11} \neq 0$，所以先解出$x_1$；
  - 再将$x_1$代入第二行，解出$x_2$。重复以上步骤直至解出$x$。

- 上三角形——回代法
  - 因为最后一行只有$a_{nn} \neq 0$，所以先解出$x_n$；
  - 再将$x_n$代入倒数第二行，解出$x_{n-1}$。重复以上步骤直至解出$x$。

### 1.2 **Gauss**变换与**Gauss**消元法

- **Gauss变换**的形式

$$
L_k = \begin{bmatrix}
1 & & &\\
& \ddots & &\\
& & 1 &\\
& & -l_{k+1, k}\;\;\;\; 1\\
& & \vdots & \ddots \\
& & -l_{n,k}& & 1
\end{bmatrix}
$$

> 其实Gauss变换可以看作一个算子，一个特殊的对应的Gauss变换作用于矩阵A，可以将A的第k列转化为$[x_1, x_2, \cdots, x_k, 0, \cdots, 0]^T$.

- 矩阵$A$的三角分解：

  - 令$L = (L_{n-1}L_{n-2}\cdots L_1)^{-1}, U = A^{(n-1)}$，则$A = LU$；
  - 其中$L$是下三角矩阵，$U$是上三角矩阵

- 三角分解之后，$Ax= b \rightarrow LUx=b \rightarrow Ly=b,Ux=y$

### 1.3 选主元三角分解

- 全主元消元法
  - $U = L_{n-1}P_{n-1}L_{n-2}P_{n-2} \cdots L_1P_1AQ_1 \cdots Q_{n-1}$，其中$P_i, Q_i$均为置换矩阵，$L_i$为Gauss变换。
  - 也就是说，从正式开始消元之前，就已经选过一轮全主元了。

- 列主元消元法
  - 只在当前主元这一列的下方选，如现在主元是$a_{kk}$，那么就是从$\{a_{ik}: k\leq i \leq n\}$中选出绝对值最大的一个，然后进行行交换——也就是说，没有全主元方法中的右端的$Q_i$了。
  - 也是从第一步消元前就先选了一轮了。

> 选主元消去是完全建立在Gauss消元法的基础上的，它只不过是在每一次消去之后，通过对比主元与其他元间的大小进而交换对应的行与列，来确保A的所有顺序主子式均非零。

### 1.4 平方根法(Cholesky分解)

- $Ax=b \rightarrow A=LL^T \rightarrow LL^Tx=b \rightarrow Ly=b,L^Tx=y$

- 改进的平方根法：$A = LDL^T$，不需要开方运算。

### 1.5 拓展

- 三对角矩阵$\rightarrow$带状矩阵
  - 追赶法求解三对角矩阵
- Doolittle分解，Courant分解
- 求逆的方法：Gauss-Jordan消元，构建增广矩阵（将两个矩阵拼起来）。

## Ch2 敏度分析与误差分析

### 2.1 范数

#### 2.1.1 向量范数

- 性质：正定性、齐次性（数乘）与三角不等式；
- 经常被使用的是p-范数——p-范数的等价性

#### 2.1.2 矩阵范数

- 比向量范数多一个相容性；

- 矩阵范数可以由向量范数诱导（大部分都是），但也有不是的（Frobenius范数（所有元素平方和开根））;
  - $||\cdot||_1$：列绝值和范数，$||\cdot||_{\infty}$：行绝对值和范数
  - $||\cdot||_2$：谱范数$= \sqrt{\lambda_{\max}(A^TA)}$

- 与矩阵条件数相关的性质与定理有很多，有意思的！

### 2.2 条件数

最常用的是矩阵条件数：$cond_A = ||A||\cdot||A^{-1}||$

### 2.3 敏度分析

敏度分析是说在求解$Ax= b$时，对$A$或$b$有一定的扰动，$x$的变换可以如何表示？或者$x$的变换是可以被控制的吗？

### 2.4 精度估计

- 相对误差与绝对误差的分析
- 舍入误差分析

## Ch3 最小二乘问题(Least Square)

### 3.1 定义与解的性质

- 定义：

$$
x^* = \mathop{\arg\min}_\limits{y\in \mathbb{R}^n} ||Ay - b||_2\\
\Longrightarrow x^* \in \mathcal{X}_{LS}
$$

- 解的存在唯一性（以下假设$A \in \mathbb{R}^{m\times n}$）

   - $A$的值域$\mathcal{R}(A) = \{y\in \mathbb{R}^m: y=Ax, x\in\mathbb{R}^n \}$；

     > $\mathcal{R}(A) = span(a_1, a_2, \cdots, a_n), a_i$为$A$的列向量

   - $A$的零空间$\mathcal{N}(A) = \{x\in \mathbb{R}^n: Ax=0\}$；

   - $S \in \mathbb{R}^n$，其正交补$S^{\perp} = \{y\in\mathbb{R}^n:y^Tx=0, \forall x \in S \}$。

   - LS问题的解总是存在的。解唯一 $\Longleftrightarrow$ $\mathcal{N}(A) = \{0\}$。

- 求解by正则化方程组/法方程组：将求解$x^*$转化为求解$A^TAx = A^Tb$.

### 3.2 Householder变换

1. 定义$H = I- 2ww^T$，其中$w \in \mathbb{R}^{n\times 1}, ||w||_2=1$
2. $H$的性质：
   1. 对称性：$H^T = H$;
   2. 正交性：$H^TH=I$;
   3. 反射性：$Hx$是$x$关于$w$的垂直超平面（$span\{w\}^{\perp}$）的镜像反射
3. $H$的求解：
   1. $v = x \pm ||x||_2e_1$;
   2. $w = \frac{v}{||v||_2}$;
   3. $H = I - 2ww^T = I - \frac{2vv^T}{v^Tv} = I - \beta vv^T, \beta=\frac{2}{v^Tv}$

### 3.3 Givens变换

> 亦称为平面旋转变换，可以选择性地将一些元素化为0.

1. G原来是一个单位阵，但第i行第i列与第k列第k行进行了一些操作。$Gy$可以使$y$的某一个分量变为0——利用三角函数性质

2. $G$的计算关键在于理解

$$
\begin{bmatrix}
cos \;\; sin\\
-sin \;\; cos
\end{bmatrix}
\begin{bmatrix}
a\\
b
\end{bmatrix}=
\begin{bmatrix}
r\\
0
\end{bmatrix}
$$

### 3.4 正交变换法求解LS问题

#### 3.4.1 QR分解定理

- 定理叙述：

   设$A \in \mathbb{R}^{m\times n}(m \geq n)$，则$A$有QR分解：

$$
A = Q\begin{bmatrix}R \\ 0\end{bmatrix},
$$

   其中$Q \in \mathbb{R}^{m\times m}$为正交阵，$R\in \mathbb{R}^{n\times n}$是具有非负对角元的上三角矩阵。且当$m = n$与$A$可逆时，上述分解唯一。

- 求解LS问题：

$$
\begin{aligned}
||Ax-b||_2^2 =& ||Q^TAx - Q^Tb||_2^2\\
=& ||\begin{bmatrix}R \\ 0\end{bmatrix}x - Q^Tb||_2^2\\
=& ||\begin{bmatrix}Rx \\ 0\end{bmatrix} - \begin{bmatrix}c_1 \\ c_2\end{bmatrix}||_2^2\\
&= ||Rx - c_1||_2^2 + ||c_2||_2^2
\end{aligned}
$$

   则$x\in \mathcal{X}_{LS} \Longleftrightarrow Rx=c_1$。

#### 3.4.2 利用Householder变换实现QR分解

- Householder变换可以将一个列向量（无论几维）变换为第一个元素非零而其他元素均为0的列向量。

- 有$H_r H_{r-1}\cdots H_2 H_1 A = \begin{bmatrix}R \\ 0\end{bmatrix}$，则$Q = H_1 H_2 \cdots H_{r-1}H_r$.

## Ch4 古典迭代法

### 4.1 单步线性定常迭代法

#### 4.1.1 一般定义

1. 新的近似解$x_k$是已知近似解$x_{k-1}$的线性函数，并且只与$x_{k-1}$有关；
2. 单步线性定常迭代法有如下形式：

$$
x_k = Tx_{k-1} + g \tag{4.1}
$$


   其中$T$是迭代矩阵，$g$是常数项

#### 4.1.2 Jacobi迭代与Gauss-Seidel迭代

考虑$Ax = b$，令$A = D -L -U$，则

- Jacobi: 

$$
T_{J} = D^{-1}(L+U), g_J = D^{-1}b \tag{4.2}
$$

- Weighted Jacobi: 

$$
x_{*} = T_{J}x^{(k)} + g_J \Longrightarrow x^{(k+1)} = (1-\omega)x^{(k)} + \omega x_{*}
$$

$$
T_{\omega} = (1-\omega)I + \omega D^{-1}(L+U) \tag{4.3}
$$



- Gauss-Seidel: 

$$
T_{GS} = (D-L)^{-1}U, g_{GS} = (D-L)^{-1}b \tag{4.4}
$$



- Backward G-S: 

$$
T_{BGS} = (D-U)^{-1}L, g_{BGS} = (D-U)^{-1}b \tag{4.5}
$$


注意在具体实现时，两种方法的分量计算顺序不尽相同！！！

Jacobi迭代法是每次完全只用第$k-1$次的，而G-S迭代法是使用混杂的分量：

- Jacobi：

$$
a_{ii}x_i^{(k+1)} = - \sum_{j\neq i;\;j=1}^{n} a_{ij}x_j^{(k)} + b_i \tag{4.6}
$$

- G-S：

$$
a_{ii}x_i^{(k+1)} = - \sum_{j=1}^{i-1} a_{ij}x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij}x_j^{(k)} + b_i \tag{4.7}
$$

- Backward G-S：

$$
a_{ii}x_i^{(k+1)} = - \sum_{j=1}^{i-1} a_{ij}x_j^{(k)} - \sum_{j=i+1}^{n} a_{ij}x_j^{(k+1)} + b_i \tag{4.8}
$$

#### 4.1.3 收敛性理论

- （一般性理论）单步线性定常迭代法收敛的充要条件是$M^k \rightarrow 0$，
  - 进一步，充要条件是$\rho(M) < 1$。

- Jacobi迭代法：
  - 若$A$对称，且其对角元均大于0，则Jacobi收敛的充要条件是$A$与$2D-A$均正定。
- G-S迭代法：
  - 若$A$对称正定，则G-S收敛。
- 若$A$严格对角占优或不可约对角占优，则Jacobi和G-S都收敛。

#### 4.1.4 收敛速度

$R_{\infty}(M) = -\ln\rho(M)$

### 4.2 超松弛迭代法

#### 4.2.1 迭代形式

分量形式：$x_i^{(k+1)} = (1-\omega)x^{(k)}_i + \omega x_i^{GS}, \; where\; x_i^{GS} \;equals \;x_i^{(k+1)}\; in\; (4.7).$

总形式：

$$
x_{k+1} = T_{SOR}x_k + \omega(D - \omega L)^{-1}b \tag{4.9}
$$

其中$T_{SOR} = (D - \omega L)^{-1}\left[(1-\omega)D + \omega U \right]$称为松弛迭代法的迭代矩阵，$\omega$称为松弛因子：

- $\omega =1$：G-S迭代法；$\omega>1$：超松弛（SOR）；$\omega < 1$：低松弛.

#### 4.2.2 收敛性分析

- SOR迭代法收敛的充要条件是$\rho(L_{\omega}) < 1$；

- SOR迭代法收敛的必要条件是$0 < \omega < 2$；

- 若$A$严格对角占优或不可约对角占优，且$\omega \in (0,1 )$，则SOR收敛；
- 若$A$是实对称矩阵，则当$\omega \in (0,2)$时，SOR收敛。

#### 4.2.3 最佳松弛因子

随着$\omega$从0开始增加，$\rho(T_{SOR})$逐渐减小，直至

$$
\omega = \omega_{opt} = \frac{2}{1 + \sqrt{1 - \rho^2(B)}}
$$

此时$\rho(T_{SOR})$达到极小：

$$
\rho(T_{opt}) = \frac{1 - \sqrt{1 - \rho^2(B)}}{1 + \sqrt{1 - \rho^2(B)}}
$$

此后$\omega$增加，$\rho$也随之增加。

## Ch5 共轭梯度法

### 5.1 最速下降法

亦称梯度下降法

- 定义二次泛函$\psi(x) = x^TAx - 2b^Tx$，若A对称正定，则求解$Ax = b$的解等价于求解二次泛函$\psi(x)$的极小值点。

- 求解思路：给定$x^{(0)}$，依次求$x^{(1)}, x^{(2)}, \cdots$，s.t. $\psi(x^{(k+1)} < \psi(x^{(k)})$

  $\rightarrow$给定$x^{(0)}$，在方向$x^{(0)} + \alpha p^{(0)}$上求$\psi(x)$的极小值点$x^{(1)}$；

  $\rightarrow \cdots$

  - 问题转化为两要点：
    - 如何确定方向序列$\{ p^{(k)}\}$；
    - 如何求得步长$\alpha$。

- 算法：

  先要取$x^{(0)} \in \mathbb{R}^n, r^{(0)} = b - Ax^{(0)}$，之后按如下方式做迭代：

$$
\left\{
\begin{aligned}
\alpha_k =& \frac{\langle r^{(k)}, r^{(k)}\rangle}{\langle r^{(k)}, Ar^{(k)}\rangle}\\
x^{(k+1)} =& x^{(k)} + \alpha_k r^{(k)}\\
r^{(k+1)} =& b - Ax^{(k+1)}
\end{aligned}
\right.
$$

### 5.2 共轭梯度法（C-G）

#### 5.2.1 与最速下降法作比较

1. 相同点：第一步仍然要取$x^{(0)} \in \mathbb{R}^n, r^{(0)} = b - Ax^{(0)}$；

2. 不同点：在第$k+1(k \geq 1)$步中，”下山“方向不再取负梯度方向。

   1. 而是在过点$x^{(k)}$由向量$r^{(k)}$和$p^{(k-1)}$所张成的平面内选取$\psi(x)$下降最快的方向：

$$
\pi = \{x = x^{(k)} + \xi r^{(k)} + \eta p^{(k-1)}: \xi, \eta \in \mathbb{R} \}
$$

#### 5.2.2 算法

1. $\alpha_k = \frac{\langle r^{(k)}, p^{(k)}\rangle}{\langle p^{(k)}, Ap^{(k)}\rangle}$
2. $x^{(k+1)} = x^{(k)} + \alpha p^{(k)}$
3. $r^{(k+1)} = b - Ax^{k+1}$
4. $\beta_k = -\frac{\langle r^{(k+1)}, Ap^{(k)}\rangle}{\langle p^{(k)}, Ap^{(k)}\rangle}$
5. $p^{(k+1)} = r^{(k+1)} + \beta_k p^{(k)}$

初始值：$\forall x^{(0)}, p^{(0)} = r^{(0)} = b - Ax^{(0)}, \alpha_0 = \frac{\langle r^{(0)}, r^{(0)}\rangle}{\langle p^{(0)}, Ap^{(0)}\rangle}, x^{(1)} = x^{(0)} + \alpha_0 p^{(0)}, r^{(1)} = b - Ax^{(1)}$。

#### 5.2.3 “共轭”的体现

1. 由C-G法得到的向量组$\{r^{(i)}\}$与$\{p^{(i)} \}$有以下性质：
   1. $\langle p^{(i)}, r^{(j)}\rangle = 0, 0 \leq i < j \leq k$；
   2. $\langle r^{(i)}, r^{(j)} \rangle = 0, 0 \leq i,j \leq k, i\neq j$；
   3. $\langle p^{(i)}, Ap^{(j)} \rangle = 0, 0 \leq i,j \leq k, i\neq j$；
   4. $span\{r^{(0)}, \cdots, r^{(k)} \} = span\{p^{(0)}, \cdots, p^{(k)} \} = \mathcal{K}(A, r^{(0)}, k+1)$，其中$\mathcal{K}(A, r^{(0)}, k+1) = span\{r^{(0)}, Ar^{(0)}, \cdots, A^k r^{(0)} \}$，通常称之为$Krylov$子空间。

## Ch6 非对称特征值问题的计算方法

### 6.1 代数基本知识与概念

#### 6.1.1 矩阵概念与性质

1. 特征值、特征多项式，谱集

   1. 特征值：$A \in \mathbb{C}^{n\times n},$ 存在$\lambda \in \mathbb{C}, s.t. Ax = \lambda x(x \neq 0)$；

   2. 特征多项式：$P_{A}(\lambda) = det(\lambda I - A) = (\lambda - \lambda_1)^{n_1} \cdots(\lambda - \lambda_r)^{n_r}$，其中$n_1 + \cdots + n_r = n$，$\lambda_1, \cdots , \lambda_r$互异，且称$\lambda_A = \{\lambda_1, \cdots, \lambda_r \}$为A的谱集。

      $n_i$称为$\lambda_i$的代数重数（简称重数），而$m_i = n - rank(\lambda_i I - A)$称为$\lambda_i$的几何重数（有$m_i \leq n_i$）。

   3. $n_i = 1$的$\lambda_i$称为单特征值，否则称为重特征值。如果$n_i = m_i$，则称$\lambda_i$为A的一个半单特征值。

      显然，单特征值必是半单特征值。

   4. 如果A所有特征值半单，则称A是非亏损的。

      A非亏损$\Longleftrightarrow$$A$有n个线性无关的特征向量（即A可对角化）。

2. 相似（变换）

   1. 相似的矩阵具有相同的特征值；
   2. A与B相似：$B = P^{-1}AP$。且若$x$为$A$的特征向量，则$Px$是$B$的特征向量


#### 6.1.2 相关定理

1. **Jordan分解定理**

   设$A \in \mathbb{C}^{n\times n}$，有$r$个互不相同的特征值$\{\lambda_1, \cdots,\lambda_r \}$，其重数分别为$n(\lambda_1), \cdots, n(\lambda_r)$，则必存在一个非奇异矩阵$P\in \mathbb{C}^{n\times n}$，使得
   
$$
P^{-1}AP = \begin{bmatrix}J(\lambda_1)& & & &\\
& J(\lambda_2) & & \\
& & \ddots &\\
& & & J(\lambda_r)
\end{bmatrix}
$$

其中$J(\lambda_i) = diag(J_1(\lambda_i), \cdots, J_{k_i}(\lambda_i)) \in \mathbb{C}^{n(\lambda_i) \times n(\lambda_i)}, i = 1, \cdots, r$，其中

$$
J_j(\lambda_i) = \begin{bmatrix}\lambda_i \quad1 & &\\
& \lambda_i \quad\ddots &\\
& & \ddots 1\\
& & & \lambda_i
\end{bmatrix} \in\mathbb{C}^{n_j(\lambda_i) \times n_j(\lambda_i)}, j = 1, \cdots, k_i
$$

$$
n_1(\lambda_i) + \cdots + n_{k_i}(\lambda_i) = n(\lambda_i), i = 1, \cdots, r;
$$

上述定理中的矩阵J称为A的Jordan标准形，其中每个子矩阵$J_i(\lambda_i)$称为Jordan块。

2. **Schur分解定理**

   设$A \in \mathbb{C}^{n \times n}$，则存在酉矩阵$U \in \mathbb{C}^{n\times n}$，使得$U^* AU = T$，其中$T$是上三角阵。

   且适当选取U，可以T的对角元按任意指定序列排序。

3. 实数域上的**Schur分解定理**

   设$A \in \mathbb{R}^{n \times n}$，则存在正交阵$Q \in \mathbb{R}^{n\times n}$，使得$Q^* AQ = \begin{bmatrix}R_{11}& R_{12}& \cdots& R_{1m}&\\ & \ddots & & \vdots \\ & & & R_{mm} \end{bmatrix}$，其中$R_{ii}$或是一个实数，或是一个具有一对共轭复数特征值的二阶方阵。

### 6.2 幂法

### 6.3 计算矩阵特征值的QR方法

#### 6.3.1 QR迭代：QR分解+反分解

1. 回顾QR分解：见**3.4.1**

2. 基本迭代与收敛性

   1. 令$A_0 =A$，作分解$A_0 = Q_1R_1$，令$A_1 = R_1Q_1$，则$A_0$与$A_1$正交相似（$A_1 = Q_1^T Q_1 R_1Q_1 = Q_1^TA_0Q_1$）；
   2. 作分解$A_1 = Q_2R_2$，令$A_2 = R_2 Q_2$；
   3. $\cdots$
   4. 作分解$A_m = Q_{m+1}R_{m+1}$，令$A_{m+1} = R_{m+1} Q_{m+1}$；

   得到矩阵序列$\{A_m\}：两两正交相似。

   令$\hat{Q}_m = Q_1 \cdots Q_m, \hat{R}_m = R_m \cdots R_1$，则有$\hat{Q}_m \hat{R}_m = A_m$。

3. 上述的$A_m = [\alpha_{ij}^{(m)}]$的对角线以下的元素趋于0，而对角元趋于A的各个特征值（特征值由大到小从左上方开始排列）。

#### 6.3.2 计算A的特征值的实现

##### 6.3.2.1 将A先转化为上Hessenberg矩阵

实际计算时，为了减少每次迭代所需的运算量，总是将原矩阵A经相似变换约化为一个**准上三角阵**，再对约化后的矩阵进行QR迭代。

1. 上Hessenberg矩阵：

$$
H = \begin{bmatrix}h_{11}& h_{12}& h_{13}& \cdots& h_{1,n-1}& h_{1n}&\\
h_{21}& h_{22}& h_{23}& \cdots& h_{2,n-1}& h_{2n}&\\
& h_{32}& h_{33}& \cdots& h_{3,n-1}& h_{3n}&\\
& & \ddots& \ddots& \vdots& \vdots&\\
& & & \ddots& \ddots& \vdots&\\
& & & & h_{n,n-1}& h_{nn}
\end{bmatrix}
$$

2. 上Hessenberg分解：

   现令$Q_0 = H_1H_2 \cdots H_{n-2}$，则有$Q_0^TAQ_0 = H$，称这个式子为A的上Hessenberg分解。

   其中$Q_0$是多个Householder变换的乘积，故$Q_0$是正交阵。

3. $Q_0$的具体计算：
   1. 首先$H_1A$会使得A的第一列只剩$a_{11} \neq 0$，但需要注意的是，我们需要进行右乘，即$H_1AH_1$，这一步操作会使A的第一列0元素变为非零元素——这不是我们想要的。
   2. 所以我们取$H_1 = \begin{bmatrix}1& 0\\ 0& \hat{H_1} \end{bmatrix}$ ，则$H_1AH_1 = \begin{bmatrix}a_{11}& a_2^T\hat{H_1}\\ \hat{H_1}a_1& \hat{H_1}A_{22}\hat{H_1} \end{bmatrix}$，其中$a_1^T = [a_{21}, a_{31}, \cdots, a_{n1}], a_2^T = [a_{12}, a_{13}, \cdots , a_{1n}]$。
   3. 最佳方案是取$\hat{H_1}, s.t. \hat{H_1}a_1 = pe_1, p\in \mathbb{R}$。（其实就是使$H_1AH_1$的第一列只有前两个非零元素）。
   4. 对$A_{22}$做类似的操作。

##### 6.3.2.2 上Hessenberg的QR迭代

因为H是一个准上三角阵（零元比较多），所以它的QR分解主要依靠Givens变换——也就是说，它的QR分解中的Q特指Givens变换。别忘了反分解（i.e.，要乘回去）。
##### 6.3.2.3 总结

1. 对A作上Hessenberg分解，$Q_0^T A Q_0 = H$；
2. 令$H_0 = H$，对$H_0$进行QR迭代——利用Givens变换。

### 6.4 实用的QR方法

#### 6.4.1 带原点位移方法

#### 6.4.2 双重步位移QR方法

## Ch7 对称特征值问题的计算方法

### 7.1 对称QR方法

### 7.2 奇异值分解定理
