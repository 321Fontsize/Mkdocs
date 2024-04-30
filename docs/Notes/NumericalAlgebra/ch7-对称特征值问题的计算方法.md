# Ch7 对称特征值问题的计算方法

## 7.1 对称QR方法

对称QR方法就是求解对称特征值问题的QR方法，是将QR方法应用于对称矩阵，并且充分利用其对称性得到的。

### 7.1.1 三对角化

若$A$是实对称阵，并假定$A$的上Hessenberg分解为$Q^TAQ = T$，则容易验证$T$必定是对称三对角阵。

具体实现：

1. 上Hessenberg分解与[Ch6-6.3.2](./ch6-非对称特征值问题的计算方法.md)中的计算方法别无二样：

   - 定义$T = \begin{bmatrix}\alpha_1 &\beta_1 & & \\ \beta_1& \alpha_2 &\ddots & \\& \ddots &\ddots& \beta_{n-1} \\ & & \beta_{n-1}& \alpha_{n} \end{bmatrix}$，

   - 令$Q = H_1 H_2 \cdots H_{n-2}, \quad H_k =diag(I_k, \hat{H}_k)$，则$Q^T A Q = T$。

2. 在第k步约化时，主要工作量是计算$\hat{H}_{k}A_{k-1}\hat{H}_k$：

   - 设$\hat{H_k} = I -\beta v v^T, \quad v \in \mathbb{R}^{n-k}$，利用$A_{k-1}$的对称性，有$\hat{H}_{k}A_{k-1}\hat{H}_k = A_{k-1} - vw^T - wv^T$，其中$w = u - \frac{1}{2}\beta(v^Tu)v, \quad u = \beta A_{k-1}v$.

### 7.1.2 对T作QR迭代——隐式对称QR迭代

由于此时A的特征值均为实数，故没有必要使用双重步位移迭代，直接使用带原点位移迭代。

我们当然可以像[Ch6-6.4](./ch6-非对称特征值问题的计算方法.md)中那样选取位移为$\mu = T^{(k)}(n, n)$，但是由于A的对称性，我们有更好的选择方法：

- $Wilkinson$位移：
  - $\mu$选取为矩阵$T^{(k)}(n-1:n, n-1:n) = \begin{bmatrix}\alpha_{n-1}& \beta_{n-1}\\ \beta_{n-1}& \alpha_n \end{bmatrix}$的两个特征值中靠近$\alpha_n$的那一个；
  - 即$\mu = \alpha_n + \delta - sgn(\delta)\sqrt{\delta^2 + \beta_{n-1}^2},\quad \delta = \frac{\alpha_{n-1}-\alpha_n}{2}$.
  - 上述两种位移取法均是三次收敛速度，但后者比前者好。
- 选好了位移，我们再来考虑如何进行一次漂亮的QR迭代：$T_k - \mu I =Q_kR_k, \quad \hat{T}_k = R_kQ_K + \mu I$：
  - 利用Givens变换实现$T_k - \mu I$的QR分解。
- 带$Wilkinson$位移的隐式对称QR迭代：

### 7.1.3 隐式对称QR算法

## 7.2 奇异值分解定理

> SVD分解定理

设$A \in \mathbb{R}^{m\times n}$，则存在正交矩阵$U \in \mathbb{R}^{m\times m}$和$Q \in \mathbb{R}^{n\times n}$，使得
$$
U^TAV = \begin{bmatrix}\Sigma_r& 0\\ 0&0\end{bmatrix}
$$
其中$\Sigma_r = diag(\sigma_1, \cdots, \sigma_r), \quad \sigma_1 \geq \cdots \geq \sigma_r > 0$.

设$A$具有如上的奇异值分解，那么我们称数

$$
\sigma_1 \geq \cdots \geq \sigma_r > \sigma_{r+1} = \cdots = \sigma_n= 0
$$

为A的奇异值，V的列向量称为A的右奇异向量，U的列向量称为A的左奇异向量。

[Next-->Ch1](./ch1-线性方程组的直接解法(Ax=b).md)

[Numerical Algebra's homepage](../NumericalAlgebra.md)
