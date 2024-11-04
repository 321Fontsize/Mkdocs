# Schur Decomposition

## 1. Theorems

### Lemma 1.1

$T\in\mathbb{C}^{n\times n}$ is partitioned as follow, 

$$
T = \begin{bmatrix}
    T_{11}&T_{12}\\0&T_{22}
\end{bmatrix},
$$

then $\lambda(T) = \lambda(T_{11})\cup \lambda(T_{22})$.

***Proof.***

Suppose $x \neq 0$ such that $Tx = \begin{bmatrix}T_{11}&T_{12}\\0&T_{22} \end{bmatrix}\begin{bmatrix}x_1\\x_2 \end{bmatrix} = \lambda \begin{bmatrix}x_1\\x_2\end{bmatrix}$.

- If $x_2\neq 0$, then $T_{22}x_2 = \lambda x_2$, which implies $\lambda\in\lambda(T_{22})$.
- If $x_2 = 0$, then $T_{11}x_1 = \lambda x_1$, which implies $\lambda\in\lambda(T_{11})$.

Hence $\lambda(T) \subset \lambda(T_{11})\cup \lambda(T_{22})$. But since both $|\lambda(T)| = n = |\lambda(T_{11}) \cup \lambda(T_{22})|$, then $\lambda(T) = \lambda(T_{11})\cup \lambda(T_{22})$.  $\blacksquare$    

### Lemma 1.2

$A\in\mathbb{C}^{n\times n}$, $B\in\mathbb{C}^{p\times p}$ and $X\in\mathbb{C}^{n\times p}$ with $\text{rank}(X)=p$ satisfy

$$
AX = XB,
$$

then there exists a unitary $Q\in\mathbb{C}^{n\times n}$ s.t.

$$
Q^H A Q = T = \begin{bmatrix}T_{11}&T_{12}\\ 0&T_{22} \end{bmatrix}
$$

with $T_{11}\in\mathbb{C}^{p\times p}$ and $T_{22}\in\mathbb{C}^{(n-p)\times (n-p)}$, together with $\lambda(T_{11}) = \lambda(A) \cap \lambda(B)$.

***Proof.***

Let $X = Q\begin{bmatrix}R_1 \\0 \end{bmatrix}, Q\in\mathbb{C}^{n\times n}, R_1\in\mathbb{C}^{p\times p}$ is nonsingular.

Then

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

where $T:= \begin{bmatrix}T_{11}&T_{12}\\T_{21}&T_{22}\end{bmatrix} = Q^H A Q$. 

Since $R_1$ is nonsingular and $T_{21}R_1 = 0$, then $T_{21}=0$.

Since $T_{11}R_1 = R_1 B$, then $\lambda(T_{11}) = \lambda(B)$. And since $A = QTQ^H$ and $T_{21}=0$ , then $\lambda(A) = \lambda(T) = \lambda(T_{11})\cup\lambda(T_{22})$, which implies that $\lambda(A)\cap \lambda(B) = \lambda(T_{11})$. $\blacksquare$

### Thm 1.3

**(Schur Decomposition).** If $A \in \mathbb{C}^{n\times n}$, then there exists a unitary $Q\in\mathbb{R}^{n\times n}$ such that $Q^HAQ = T := D+N$, where $D = \text{diag}(\lambda_1,\cdots, \lambda_n)$ and $N$ is strictly upper triangular.

***Proof.***

The conclusion holds for $n=1$.

Suppose it holds for all matrices of order $n-1$. If $Ax = \lambda x$ with $x\neq 0$, by Lemma 1.2(with $B = (\lambda)\in\mathbb{C}^{n\times n}$ in Lemma 1.2), there exists a unitary $U$ such that

$$
U^H A U = \begin{bmatrix}
    \lambda& w^H\\ 0 & C
\end{bmatrix}.
$$

By induction, there is a unitary $\tilde{U}$ such that $\tilde{U}^H C \tilde{U}$ is upper triangular.

Thus, if $Q = U\cdot\text{diag}(1, \tilde{U})$, then $Q^HAQ$ is upper triangular. $\blacksquare$

---

If $Q = [q_1, \cdots, q_n]$, then $q_i$ are referred to as *Schur vectors*. By $AQ = QT$, we have

$$
Aq_k = \lambda_k q_k + \sum_{i=1}^{k-1}n_{ik}q_i,\quad k=1:n.
$$

From this, we can conclude that the subspaces 

$$
S_k = \text{span}\{q_1, \cdots, q_k \}, \quad k=1:n
$$

are invariant for $A$ (i.e., $\forall x\in S_k, Ax\in S_k$). Moreover, if $Q_k = [q_1, \cdots, q_k]$, then $\lambda(Q_k^H AQ_k) = \{\lambda_1, \cdots, \lambda_k\}$.   

## 2. Applications

### Normal matrix

Defintion: $A$ is normal if $A^HA = AA^H$.

#### Thm 1.4

$A \in\mathbb{C}^{n\times n}$ is normal $\Longleftrightarrow$ there exists a unitary $Q\in\mathbb{C}^{n\times n}$ such that $Q^H A Q = \text{diag}(\lambda_1, \cdots, \lambda_n)$.

***Proof.***

$\Longleftarrow$:

Denote $D = \text{diag}(\lambda_1 ,\cdots, \lambda_n)$, then $A^HA = Q D^H Q^H Q D Q^H = Q D^H D Q^H$ and $AA^H = Q D D^H Q^H$. Since $D^H D = D D^H$, then $A^HA = AA^H$.

$\Longrightarrow$:

If $A$ is normal, set $T = Q^H A Q$, where $Q$ is the unitary in *Schur Decomposition*, then $T^H T = T T^H$, which implies $T$ is normal. 

Since $T$ is upper triangular, then we have

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
    \overline{t_{1n}}& \cdots& \overline{t_{nn}}
\end{bmatrix}
$$

For the (1,1) position, we have

$$
\text{LHS} = |t_{11}|^2 = \sum_{k=1}^{n} |t_{1k}|^2 = \text{RHS},
$$

which implies $t_{1k} = 0,\; k=2:n$, then for the $(2,2)$ position, we have

$$
\text{LHS} = |t_{12}|^2 + |t_{22}|^2 = |t_{22}|^2 = \sum_{k=2}^{n}|t_{2k}|^2 = \text{RHS},
$$

which implies $t_{2k}=0,\; k=3:n$. Do as the above, finally, we can get $T = \text{diag}(t_{11}, \cdots, t_{nn})$. $\blacksquare$ 