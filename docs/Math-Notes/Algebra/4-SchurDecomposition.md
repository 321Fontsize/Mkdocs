# Schur 分解

## 1. 定理

### 引理 1.1

设$T\in\mathbb{C}^{n\times n}$ 被划分如下，

$$
T = \begin{bmatrix}
    T_{11}&T_{12}\\0&T_{22}
\end{bmatrix},
$$

则$\lambda(T) = \lambda(T_{11})\cup \lambda(T_{22})$。

***证明***

假设存在$x \neq 0$ 使得$Tx = \begin{bmatrix}T_{11}&T_{12}\\0&T_{22} \end{bmatrix}\begin{bmatrix}x_1\\x_2 \end{bmatrix} = \lambda \begin{bmatrix}x_1\\x_2\end{bmatrix}$。

- 如果$x_2\neq 0$，则$T_{22}x_2 = \lambda x_2$，这意味着$\lambda\in\lambda(T_{22})$。
- 如果$x_2 = 0$，则$T_{11}x_1 = \lambda x_1$，这意味着$\lambda\in\lambda(T_{11})$。

因此$\lambda(T) \subset \lambda(T_{11})\cup \lambda(T_{22})$。但由于$|\lambda(T)| = n = |\lambda(T_{11}) \cup \lambda(T_{22})|$，则$\lambda(T) = \lambda(T_{11})\cup \lambda(T_{22})$。  $\blacksquare$

### 引理 1.2

设$A\in\mathbb{C}^{n\times n}$，$B\in\mathbb{C}^{p\times p}$ 和$X\in\mathbb{C}^{n\times p}$ 且$\text{rank}(X)=p$ 满足

$$
AX = XB,
$$

则存在一个酉矩阵$Q\in\mathbb{C}^{n\times n}$ 使得

$$
Q^H A Q = T = \begin{bmatrix}T_{11}&T_{12}\\ 0&T_{22} \end{bmatrix}
$$

其中$T_{11}\in\mathbb{C}^{p\times p}$ 和$T_{22}\in\mathbb{C}^{(n-p)\times (n-p)}$，且$\lambda(T_{11}) = \lambda(A) \cap \lambda(B)$。

***证明***

设$X = Q\begin{bmatrix}R_1 \\0 \end{bmatrix}, Q\in\mathbb{C}^{n\times n}, R_1\in\mathbb{C}^{p\times p}$ 是非奇异的。

则

$$
\begin{aligned}
&AQ\begin{bmatrix}
    R_1\\0
\end{bmatrix} = Q\begin{bmatrix}
    R_1\\0
\end{bmatrix}B,\\
\Longrightarrow & \begin{bmatrix}
    T_{11}&T_{12}\\T_{21}&T_{22}
\end{bmatrix}\begin{bmatrix}
    R_1\\0
\end{bmatrix} = \begin{bmatrix}
    R_1 B\\0
\end{bmatrix},
\end{aligned}
$$

其中$T:= \begin{bmatrix}T_{11}&T_{12}\\T_{21}&T_{22}\end{bmatrix} = Q^H A Q$。

由于$R_1$ 是非奇异的且$T_{21}R_1 = 0$，则$T_{21}=0$。

由于$T_{11}R_1 = R_1 B$，则$\lambda(T_{11}) = \lambda(B)$。又因为$A = QTQ^H$ 且$T_{21}=0$，则$\lambda(A) = \lambda(T) = \lambda(T_{11})\cup\lambda(T_{22})$，这意味着$\lambda(A)\cap \lambda(B) = \lambda(T_{11})$。 $\blacksquare$

### 定理 1.3

**(Schur 分解)。** 如果$A \in \mathbb{C}^{n\times n}$，则存在一个酉矩阵$Q\in\mathbb{C}^{n\times n}$ 使得$Q^HAQ = T := D+N$，其中$D = \text{diag}(\lambda_1,\cdots, \lambda_n)$ 且$N$ 是严格上三角的。

***证明***

当$n=1$ 时，结论成立。

假设对于所有$n-1$ 阶矩阵结论成立。如果$Ax = \lambda x$ 且$x\neq 0$，由引理 1.2（其中$B = (\lambda)\in\mathbb{C}^{n\times n}$ 在引理 1.2 中），存在一个酉矩阵$U$ 使得

$$
U^H A U = \begin{bmatrix}
    \lambda& w^H\\ 0 & C
\end{bmatrix}.
$$

通过归纳，存在一个酉矩阵$\tilde{U}$ 使得$\tilde{U}^H C \tilde{U}$ 是上三角的。

因此，如果$Q = U\cdot\text{diag}(1, \tilde{U})$，则$Q^HAQ$ 是上三角的。 $\blacksquare$

---

如果$Q = [q_1, \cdots, q_n]$，则$q_i$ 被称为*Schur 向量*。由$AQ = QT$，我们有

$$
Aq_k = \lambda_k q_k + \sum_{i=1}^{k-1}n_{ik}q_i,\quad k=1:n.
$$

由此，我们可以得出子空间

$$
S_k = \text{span}\{q_1, \cdots, q_k \}, \quad k=1:n
$$

对于$A$ 是不变的（即$\forall x\in S_k, Ax\in S_k$）。此外，如果$Q_k = [q_1, \cdots, q_k]$，则$\lambda(Q_k^H AQ_k) = \{\lambda_1, \cdots, \lambda_k\}$。

## 2. 应用

### 正规矩阵

定义：如果$A^HA = AA^H$，则$A$ 是正规的。

#### 定理 1.4

$A \in\mathbb{C}^{n\times n}$ 是正规的$\Longleftrightarrow$ 存在一个酉矩阵$Q\in\mathbb{C}^{n\times n}$ 使得$Q^H A Q = \text{diag}(\lambda_1, \cdots, \lambda_n)$。

***证明***

$\Longleftarrow$：

设$D = \text{diag}(\lambda_1 ,\cdots, \lambda_n)$，则$A^HA = Q D^H Q^H Q D Q^H = Q D^H D Q^H$ 且$AA^H = Q D D^H Q^H$。由于$D^H D = D D^H$，则$A^HA = AA^H$。

$\Longrightarrow$：

如果$A$ 是正规的，设$T = Q^H A Q$，其中$Q$ 是*Schur 分解*中的酉矩阵，则$T^H T = T T^H$，这意味着$T$ 是正规的。

由于$T$ 是上三角的，我们有

$$
\begin{bmatrix}
    \overline{t_{11}}& &\\
    \vdots& \ddots&\\
    \overline{t_{1n}}& \cdots& \overline{t_{nn}}
\end{bmatrix}\begin{bmatrix}
    {t_{11}}& \cdots &t_{1n} \\
    & \ddots& \vdots\\
    & & {t_{nn}}
\end{bmatrix} = \begin{bmatrix}
    {t_{11}}& \cdots &t_{1n} \\
    & \ddots& \vdots\\
    & & {t_{nn}}
\end{bmatrix}\begin{bmatrix}
    \overline{t_{11}}& &\\
    \vdots& \ddots&\\
    \overline{t_{1n}}
& \cdots& \overline{t_{nn}}
\end{bmatrix}
$$

对于(1,1)位置，我们有

$$
\text{LHS} = |t_{11}|^2 = \sum_{k=1}^{n} |t_{1k}|^2 = \text{RHS},
$$

这意味着$t_{1k} = 0,\; k=2:n$，然后对于(2,2)位置，我们有

$$
\text{LHS} = |t_{12}|^2 + |t_{22}|^2 = |t_{22}|^2 = \sum_{k=2}^{n}|t_{2k}|^2 = \text{RHS},
$$

这意味着$t_{2k}=0,\; k=3:n$。如上所述，最终我们可以得出$T = \text{diag}(t_{11}, \cdots, t_{nn})$。 $\blacksquare$
