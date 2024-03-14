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
\alpha_k =& \frac{(r^{(k)}, r^{(k)})}{(r^{(k)}, Ar^{(k)})}\\
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

1. $\alpha_k = \frac{(r^{(k)}, p^{(k)})}{(p^{(k)}, Ap^{(k)})}$
2. $x^{(k+1)} = x^{(k)} + \alpha p^{(k)}$
3. $r^{(k+1)} = b - Ax^{k+1}$
4. $\beta_k = -\frac{(r^{(k+1)}, Ap^{(k)})}{(p^{(k)}, Ap^{(k)})}$
5. $p^{(k+1)} = r^{(k+1)} + \beta_k p^{(k)}$

初始值：$\forall x^{(0)}, p^{(0)} = r^{(0)} = b - Ax^{(0)}, \alpha_0 = \frac{(r^{(0)}, r^{(0)})}{(p^{(0)}, Ap^{(0)})}, x^{(1)} = x^{(0)} + \alpha_0 p^{(0)}, r^{(1)} = b - Ax^{(1)}$。

#### 5.2.3 “共轭”的体现

## Ch6 非对称特征值问题的计算方法

### 6.1 代数基本知识与概念

#### 6.1.1 矩阵概念与性质

1. 特征多项式，谱集
2. 代数重数，几何重数
3. 相似（变换）

#### 6.1.2 相关定理

1. **Jordan分解定理**
2. **Schur分解定理**

### 6.2 幂法

### 6.3 QR方法

## Ch7 对称特征值问题的计算方法
