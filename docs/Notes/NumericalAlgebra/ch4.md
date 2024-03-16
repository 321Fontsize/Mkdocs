# Ch4 古典迭代法

## 4.1 单步线性定常迭代法

### 4.1.1 一般定义

1. 新的近似解$x_k$是已知近似解$x_{k-1}$的线性函数，并且只与$x_{k-1}$有关；
2. 单步线性定常迭代法有如下形式：

$$
x_k = Tx_{k-1} + g \tag{4.1}
$$


   其中$T$是迭代矩阵，$g$是常数项

### 4.1.2 Jacobi迭代与Gauss-Seidel迭代

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

### 4.1.3 收敛性理论

- （一般性理论）单步线性定常迭代法收敛的充要条件是$M^k \rightarrow 0$，
  - 进一步，充要条件是$\rho(M) < 1$。

- Jacobi迭代法：
  - 若$A$对称，且其对角元均大于0，则Jacobi收敛的充要条件是$A$与$2D-A$均正定。
- G-S迭代法：
  - 若$A$对称正定，则G-S收敛。
- 若$A$严格对角占优或不可约对角占优，则Jacobi和G-S都收敛。

### 4.1.4 收敛速度

$R_{\infty}(M) = -\ln\rho(M)$

## 4.2 超松弛迭代法

### 4.2.1 迭代形式

分量形式：$x_i^{(k+1)} = (1-\omega)x^{(k)}_i + \omega x_i^{GS}, \; where\; x_i^{GS} \;equals \;x_i^{(k+1)}\; in\; (4.7).$

总形式：

$$
x_{k+1} = T_{SOR}x_k + \omega(D - \omega L)^{-1}b \tag{4.9}
$$

其中$T_{SOR} = (D - \omega L)^{-1}\left[(1-\omega)D + \omega U \right]$称为松弛迭代法的迭代矩阵，$\omega$称为松弛因子：

- $\omega =1$：G-S迭代法；$\omega>1$：超松弛（SOR）；$\omega < 1$：低松弛.

### 4.2.2 收敛性分析

- SOR迭代法收敛的充要条件是$\rho(L_{\omega}) < 1$；

- SOR迭代法收敛的必要条件是$0 < \omega < 2$；

- 若$A$严格对角占优或不可约对角占优，且$\omega \in (0,1 )$，则SOR收敛；
- 若$A$是实对称矩阵，则当$\omega \in (0,2)$时，SOR收敛。

### 4.2.3 最佳松弛因子

随着$\omega$从0开始增加，$\rho(T_{SOR})$逐渐减小，直至

$$
\omega = \omega_{opt} = \frac{2}{1 + \sqrt{1 - \rho^2(B)}}
$$

此时$\rho(T_{SOR})$达到极小：

$$
\rho(T_{opt}) = \frac{1 - \sqrt{1 - \rho^2(B)}}{1 + \sqrt{1 - \rho^2(B)}}
$$

此后$\omega$增加，$\rho$也随之增加。

[Next-->Ch5](./ch5.md)

[Numerical Algebra's homepage](../NumericalAlgebra.md)