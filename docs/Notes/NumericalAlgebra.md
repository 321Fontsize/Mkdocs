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

1. 新的近似解$x_k$是已知近似解$x_{k-1}$的线性函数，并且只与$x_{k-1}$有关

2. 有如下形式：
   $$
   x_k = Mx_{k-1} + g
   $$
   其中$M$是迭代矩阵，$g$是常数项

#### 4.1.2 Jacobi迭代与Gauss-Seidel迭代

考虑$Ax = b$，令$A = D -L -U$，则

- Jacobi: $M = D^{-1}(L+U), g = D^{-1}b$；
- Gauss-Seidel: $M = (D-L)^{-1}U, g = (D-L)^{-1}b$

注意在具体实现时，两种方法的分量计算顺序不尽相同！！！

- Jacobi迭代法是每次完全只用第$k-1$次的，而G-S迭代法是使用混杂的分量！！！

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

$$
x_{k+1} = L_{\omega}x_k + \omega(D - \omega L)^{-1}b
$$

其中$L_{\omega} = (D - \omega L)^{-1}\left[(1-\omega D) + \omega U \right]$称为松弛迭代法的迭代矩阵，$\omega$称为松弛因子：

- $\omega =1$：G-S迭代法；$\omega>1$：超松弛（SOR）；$\omega < 1$：低松弛。

#### 4.2.2 收敛性分析

- SOR迭代法收敛的充要条件是$\rho(L_{\omega}) < 1$；

- SOR迭代法收敛的必要条件是$0 < \omega < 2$；

- 若$A$严格对角占优或不可约对角占优，且$\omega \in (0,1 )$，则SOR收敛；
- 若$A$是实对称矩阵，则当$\omega \in (0,2)$时，SOR收敛。

#### 4.2.3 最佳松弛因子^^^
