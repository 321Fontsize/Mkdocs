# QR Decomposition

## 1. Theorems

### Thm 1.1 QR Decomposition

If $A\in\mathbb{R}^{m\times n}$, then there exists an orthogonal $Q\in\mathbb{R}^{m\times m}$ and an upper triangular $R\in\mathbb{R}^{m\times n}$ such that $A = QR$.

***Proof.***

- For $n=1$, let $Q$ be a Householder matrix w.r.t. $A \in \mathbb{R}^{m\times 1}$: 
    
$$v = A \pm \|A\|_2e_1, Q = I - \frac{2vv^{\top}}{v^{\top}v}, 
$$

and note that $Q^{\top} = Q$, then $R = Q^{\top}A$ with $R(2:m)=0$.

- For general $n$, let $A = [A_1|u]$ with $u = A(:,n)$ and $A_1\in\mathbb{R}^{m\times(n-1)}$.

    By induction, there exists an orthogonal $Q_1\in\mathbb{R}^{m\times m}$ such that $R_1 = Q_1^{\top}A_1$ is upper triangular, which implies that if $R_1 = \begin{bmatrix}R_{11}\\R_{12}\end{bmatrix}$ with $R_{11}\in\mathbb{R}^{m\times m}$, then $R_{12}=0$.  

    Set $w = Q_1^{\top}u$, i.e., $u = Q_1 w$.

    Let $u(n:m) = Q_2R_2$ be the QR Decomposition of $u(n:m)$. 

    If $Q := Q_1\begin{bmatrix}I_{n-1}& 0\\ 0& Q_2\end{bmatrix}$, $R := \begin{bmatrix}R_{11}& w(1:n-1)\\ R_{12}&R_2 \end{bmatrix}$, then $Q$ is orthogonal and $R$ is upper triangular and  

$$
\begin{aligned}
QR &= Q_1 \begin{bmatrix}I_{n-1}& 0\\ 0& Q_2\end{bmatrix}\begin{bmatrix}R_{11}& w(1:n-1)\\ R_{12}&R_2 \end{bmatrix}\\
&= Q_1 \begin{bmatrix}R_{11}& w(1:n-1)\\ 0 & Q_2R_2=w(n:m) \end{bmatrix}\\
&= [A_1|u]= A,
\end{aligned}
$$

which completes the proof. $\blacksquare$

### Thm 1.2 Thin QR Decomposition

If $A\in\mathbb{R}^{m\times n}$ has full column rank, then there exists a unique QR Decomposition 

$$
A = Q_1 R_1
$$

with $Q\in\mathbb{R}^{m\times m}$ has orthonormal columns and $R\in\mathbb{R}^{m\times n}$ is upper triangular with positive diagonal entries.

## 2. Computing QR Decomposition

### 2.1 Householder QR

### 2.2 Gram-Schmidt QR