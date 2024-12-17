# Positive Definite

## 1. 定义

A is positive definite if $\forall x \neq 0\in \mathbb{R}^n, x^TAx > 0$.

## 2. 等价性

以下命题等价：

- $A\in\mathbb{C}^{n\times n}$ 是对称正定阵.

- $\lambda(A)\in \mathbb{R}^+$.

    ***Proof1.*** 假设 $(\lambda, x)$ 是 $A$ 的特征对，则 $x \neq 0, \Longrightarrow \lambda x^Tx = x^TAx > 0 \Longrightarrow \lambda > 0$.

    ***Proof2.*** 根据 Schur 分解定理，$A = U TU^H$，其中 $U$ 是酉矩阵而 $T$ 是上三角阵. 由于 $A = A^H$，故 $T$ 必定是对角的. 故 $A$ 正定 $\Longleftrightarrow T$ 的对角元均为正数.

- 由 $A$ 所定义的[半双线性形式](https://zh.wikipedia.org/wiki/%E5%8D%8A%E5%8F%8C%E7%BA%BF%E6%80%A7%E5%BD%A2%E5%BC%8F) $\langle x,y  \rangle _A = x^H A y$ 是 $\mathbb{C}^n$ 上的一个内积
   > 实际上，所有 $\mathbb{C}^n$ 上的内积都可认为是由此种方式定义得到.

- $A$ 的所有顺序主子式均为正，这也被称为 Sylvester's criterion.
  - ***Proof*** **ref.** [Sylvester's criterion](https://en.wikipedia.org/wiki/Sylvester%27s_criterion)
  - 类似的，我们有 $A$ 半正定 $\Longrightarrow A$ 的所有顺序主子式均为非负数，反之不然，如 $\begin{bmatrix}1& 1& 1\\ 1& 1& 1\\1&1&0 \end{bmatrix}$.
- 存在唯一的下三角矩阵 $L$，其主对角线上的元素全是正的，使得 $M = L^H L$. 这一分解被称为 [Cholesky decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition).
  - ***Proof*** **ref.** Decomposition of [this page](https://en.wikipedia.org/wiki/Definite_matrix)

对于实对称矩阵，只需将上述性质中的 $\mathbb{C}^n$ 改为 $\mathbb{R}^n$，并将“共轭转置”改为“转置”即可.

## 3. 性质

