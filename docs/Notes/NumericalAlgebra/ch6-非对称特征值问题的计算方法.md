# Ch6 非对称特征值问题的计算方法

## 6.1 代数基本知识与概念

### 6.1.1 矩阵概念与性质

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
   2. A与B相似：$B = P^{-1}AP$。且若$x$为$A$的特征向量，则$Px$是$B$的特征向量。


### 6.1.2 相关定理

- **Jordan分解定理**

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

An Example:

$$
\begin{aligned}
J = \begin{bmatrix}
2 & & & \\
& 2& 1& \\
& & 2& \\
& & & 3&
\end{bmatrix}
\end{aligned},
$$

则 $A$ 共有2个互不相同的特征值$\lambda_1 = 2, \lambda_2 = 3$，其重数分别为 $3, 1$. 且 $J(\lambda_1) = \begin{bmatrix}
   2 & &  \\
& 2& 1 \\
& & 2 \\
\end{bmatrix}, J(\lambda_2) = [3] $.

- **Schur分解定理**

  设$A \in \mathbb{C}^{n \times n}$，则存在酉矩阵$U \in \mathbb{C}^{n\times n}$，使得$U^* AU = T$，其中$T$是上三角阵。

  且适当选取U，可以T的对角元按任意指定序列排序。

- <a name='Schur'>实数域上的**Schur分解定理**</a>

  设$A \in \mathbb{R}^{n \times n}$，则存在正交阵$Q \in \mathbb{R}^{n\times n}$，使得$Q^* AQ = \begin{bmatrix}R_{11}& R_{12}& \cdots& R_{1m}&\\ & \ddots & & \vdots \\ & & & R_{mm} \end{bmatrix}$，其中$R_{ii}$或是一个实数，或是一个具有一对共轭复数特征值的二阶方阵。称上述形式为A的实Schur标准形。

## 6.2 幂法

幂法是计算一个矩阵的模最大特征值和对应的特征向量的一种迭代方法。

### 6.2.1 幂法的思想

1. 设$A \in \mathbb{C}^{n\times n}$可对角化，有如下分解

$$
A = X\Lambda X^{-1}
$$

   其中$\Lambda = diag(\lambda_1, \cdots, \lambda_n), X = [x_1, \cdots, x_n] \in \mathbb{C}^{n\times n}$非奇异，再假定

$$
|\lambda_1| > |\lambda_2| \geq \cdots \geq |\lambda_n|.
$$

2. 现任取$u_0 \in \mathbb{C}^n$，由于X的列向量构成$\mathbb{C}^n$的一组基，故
   
$$
   u_0 = \alpha_1x_1 + \cdots + \alpha_n x_n, \alpha_i \in \mathbb{C}
$$

   这样，我们有

$$
\begin{aligned}
A^k u_0 =& \sum_{j=1}^{n} \alpha_j A^k x_j\\
=& \sum_{j=1}^{n}\alpha_j \lambda_j^{k} x_j\\
=& \lambda_1^k \left[\alpha_1x_1 + \sum_{j=2}^n \alpha_j(\frac{\lambda_j}{\lambda_1})^k x_j \right]
\end{aligned}
$$

   令$k \rightarrow +\infty$，有$\frac{A^k u_0}{\lambda_1^k} \rightarrow \alpha_1x_1$。这表明当$\alpha_1 \neq 0$且k充分大时，向量$u_k = \frac{A^k u_0}{\lambda_1^k}$是A的一个很好的近似特征向量。

### 6.2.2 实用的幂法

- **6.2.1**中幂法的弊端：

  - 我们事先不知道$A$的特征值$\lambda_1$；

  - 对充分大的$k$，计算$A^k$的工作量太大

- 解决思路：

  - $\lambda_1^k$并不改变向量$A^k u_0$的方向，而只改变它的大小，因此我们不必非用$\lambda_1^k$来约化$A^k u_0$，可以用其他合理的常数（注意，约化是必须的：防止溢出）；

  - 计算$A^k u_0$并非一定要计算出$A^k$，可以选择使用迭代。

- 实用的迭代格式：

$$
\begin{aligned}
& y_k = Au_{k-1}\\
& \mu_k = \zeta_j^{(k)}, \zeta_j^{(k)}是y_k的模最大分量 \\
& u_k = y_k / \mu_k
\end{aligned}
$$

其中$u_0 \in \mathbb{C}^{n}$是可以任意给定的初始向量。通常要求$||u_0||_{\infty} = 1$。

- 若$\lambda_1$是半单的（即几何重数与代数重数相等）且模最大，取$u_0$：$u_0$在$\lambda_1$的特征子空间投影不等于0，那么按上述实用迭代格式得到的数值序列$\{\mu_k\}$收敛到模最大特征值$\lambda_1$，向量序列$\{u_k\}$收敛到其对应的特征向量。

### 6.2.3 幂法拓展

幂法一般只能用来求模最大特征值及其对应的特征向量，如果想求第二大、第i大的特征值，则需要”降阶“：

在已知$\lambda_1$及其对应的特征向量的情况下，寻找酉矩阵P，将A化为如下形式（利用$Ax_1 = \lambda_1x_1$与$Px_1 = \alpha e_1$两个等式）：

$$
P^* A P =  \begin{bmatrix}\lambda_1& *&\\
& B
\end{bmatrix}
$$

继续对矩阵$B$做幂法即可。

## 6.3 计算矩阵特征值的QR方法

### 6.3.1 QR迭代：QR分解+反分解

1. 回顾QR分解：见[**3.4.1**](./ch3-LeastSquare.md)

2. 基本迭代与收敛性

   1. 令$A_0 =A$，作分解$A_0 = Q_1R_1$，令$A_1 = R_1Q_1$，则$A_0$与$A_1$正交相似（$A_1 = Q_1^T Q_1 R_1Q_1 = Q_1^TA_0Q_1$）；
   2. 作分解$A_1 = Q_2R_2$，令$A_2 = R_2 Q_2$；
   3. $\cdots$
   4. 作分解$A_m = Q_{m+1}R_{m+1}$，令$A_{m+1} = R_{m+1} Q_{m+1}$；

   得到矩阵序列$\{A_m\}$：两两正交相似——*也就是说，我们实施QR迭代其实是在对A一步步做正交相似变换*。

   令$\hat{Q}_m = Q_1 \cdots Q_m, \hat{R}_m = R_m \cdots R_1$，则有$\hat{Q}_m \hat{R}_m = A_m$。

3. 上述的$A_m = [\alpha_{ij}^{(m)}]$的对角线以下的元素趋于0，而对角元趋于A的各个特征值（特征值由大到小从左上方开始排列）。

4. 不难发现，上述的$A_m$最终会逼近于$A$的[实Schur标准形](#Schur)。

**故我们的目的即是计算$A$的实Schur标准形！！！**

### 6.3.2 计算A的特征值的实现

- **将A先转化为上Hessenberg矩阵**

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
   1. 首先$H_1A$会使得A的第一列只剩$a_{11} \neq 0$。但需要注意的是在左乘操作结束后，我们需要进行右乘，即$H_1AH_1$，这一步操作会使A的第一列0元素变为非零元素——这不是我们想要的——我们想保留尽可能多的零元素。
   2. 所以我们取$H_1 = \begin{bmatrix}1& 0\\ 0& \hat{H_1} \end{bmatrix}$ ，则$H_1AH_1 = \begin{bmatrix}a_{11}& a_2^T\hat{H_1}\\ \hat{H_1}a_1& \hat{H_1}A_{22}\hat{H_1} \end{bmatrix}$，其中$a_1^T = [a_{21}, a_{31}, \cdots, a_{n1}], a_2^T = [a_{12}, a_{13}, \cdots , a_{1n}]$。
   3. 最佳方案是取$\hat{H_1}, s.t. \hat{H_1}a_1 = pe_1, p\in \mathbb{R}$。（其实就是使$H_1AH_1$的第一列只有前两个非零元素）。
   4. 对$A_{22}$做类似的操作。

- **上Hessenberg的QR迭代**
  - 因为H是一个准上三角阵（零元比较多），所以它的QR分解主要依靠Givens变换——也就是说，它的QR分解中的Q特指Givens变换：
    - 对于一般的n阶上Hessenberg矩阵，我们可以确定n-1个平面旋转变换$G_{12}, G_{23}, \cdots , G_{n-1, n}$，使得$G_{n-1,n}\cdots G_{12}H = R$（R是上三角阵）。
    - 则QR分解中$Q = (G_{n-1,n}\cdots G_{12})^T$；
  - 计算$\hat{H} = RQ$。
  - 至此H的一次QR迭代结束。

- **总结**

1. 对A作上Hessenberg分解，$Q_0^T A Q_0 = H$，其中$Q_0 = H_1H_2 \cdots H_{n-2}$；
2. 令$H_0 = H$，对$H_0$进行QR迭代——利用Givens变换。

## 6.4 带原点位移方法

- 迭代格式如下：

$$
\begin{aligned}
H_m -\mu_mI = Q_mR_m\\
H_{m+1} = R_mQ_m + \mu_mI
\end{aligned}
$$

$m$从0开始，$H_0$参见**6.3.2**。

- 讨论位移$\mu_m$的选取：$\mu_m = h_{nn}^{(m)}$，也就是第m次迭代中$H_m$的右下角元素。
- 通过原点位移，特征值的渐近收敛速度从线性收敛加速而变成二次收敛。

## 6.5 双重步位移QR方法

**6.4**中的带原点位移的QR迭代存在严重的缺点：若A具有复共轭特征值，则实位移一般不能起到加速作用。

为了克服这一缺点，我们使用**双重步位移的QR迭代**，其基本思想是将两步带原点位移的QR迭代合并为一步，避免复数运算。

### 6.5.1 理论分析

- 考虑$H_k$的尾部$2 \times 2$子矩阵

$$
Z_k = \begin{bmatrix}h_{mm}^{(k)} &h_{mn}^{(k)}\\
h_{nm}^{(k)} &h_{nn}^{(k)}
\end{bmatrix}, \quad m = n-1
$$

有一对复共轭特征值$\mu_1, \mu_2$，此时我们就不能期望$h_{nn}^{(k)}$能收敛到A的某个特征值。我们想取$\mu_1$或$\mu_2$作为位移，但是这就引入了复数运算。

- 为了避免复数运算，我们用$\mu_1 \& \mu_2$连续作两次位移：

$$
\begin{aligned}
H-\mu_1I = U_1R_1, & \quad H_1 = R_1U_1 + \mu_1I,\\
H_1-\mu_2I = U_2R_2, & \quad H_2 = R_2U_2 + \mu_2I.
\end{aligned}
$$

这里的$H$即是$H_k$。对上面的迭代进行一些简单的运算——计算$H_2$：

1. 令$M = H^2 - sH + tI$，其中$s = \mu_1 + \mu_2 = h_{mm}^{(k)} + h_{nn}^{(k)} \in \mathbb{R}, \quad t = \mu_1\mu_2 = \det(Z_k) \in \mathbb{R}$.
2. 计算的QR分解：$M = (H-\mu_1I)(H-\mu_2I) \Longrightarrow M = QR, \quad Q = U_1U_2, R = R_1R_2$.
3. 计算$H_2 = Q^THQ$.

### 6.5.2 Francis双重步位移QR迭代算法

上述理论分析可行，但是在第一步——计算$M$——时就出现了$O(n^3)$的计算量，这是不能接受的。

我们的目标是$H \rightarrow H_2 , O(n^2)$：

- 首先，由$M=QR$与$R$是上三角阵可知：$M$与$Q$的第一列是共线的（其实Q的第一列即是M的第一列单位化而得）。又由$M = H^2 - sH + tI$，计算得出$Me_1 = (\xi_1, \xi_2, \xi_3, 0, \cdots, 0)^T$，其中

$$
\begin{aligned}
\xi_1 &= (h_{11}^{(k)})^2 + h_{12}^{(k)}h_{21}^{(k)} - sh_{11}^{(k)} + t,\\
\xi_2 &= h_{21}^{(k)}h_{11}^{(k)} + h_{22}^{(k)}h_{21}^{(k)} - sh_{21}^{(k)},\\
\xi_3 &= h_{32}^{(k)}h_{21}^{(k)}.
\end{aligned}
$$

- 其次，如果有Householder变换$P_0$将$Me_1$变换为$\alpha e_1$，则$P_0$的第一列与$Me_1$共线，从而$P_0$第一列就可作为$Q$的第一列，即$P_0e_1 = Qe_1$。由Householder变换的求解可知，$P_0 = diag(\hat{P_0}, I_{n-3})$，其中

$$
\begin{aligned}
\hat{P_0} = I_3 - \beta vv^T, \quad& v = \begin{bmatrix}\xi_1-\alpha \\ \xi_2\\ \xi_3 \end{bmatrix},\\
\alpha=(\xi_1^2+\xi_2^2+\xi_3^2)^{\frac{1}2}, \quad & \beta = \frac{2}{v^Tv}.
\end{aligned}
$$

令$B = P_0HP_0$$\rightarrow$

- 约化$B$上Hessenberg矩阵：
  - $P_k = diag(I_k, \hat{P_k}, I_{n-k-3}), \quad k = 1, \cdots, n-3$，其中$\hat{P_k}$是三阶Householder变换；
  - $P_{n-2} = diag(I_{n-2}, \hat{P}_{n-2})$，其中$\hat{P}_{n-2}$是二阶Householder变换；
  - 则$P \triangleq P_0P_1 \cdots P_{n-3}P_{n-2}$，最后得出$\hat{H_2} = P^THP$。
  - 我们想要的是$H_2 = Q^THQ$：由于H为不可约上Hessenberg矩阵且P、Q第一列相同，则$\hat{H_2} = H_2$。

## 6.6 隐式QR算法

[Next-->Ch7](./ch7-对称特征值问题的计算方法.md)

[Numerical Algebra's homepage](../NumericalAlgebra.md)