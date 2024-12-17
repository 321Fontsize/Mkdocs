# Matrix Norm

## 1. 矩阵范数的定义与性质

### 定义 1.1-矩阵范数

设 $A \in \mathbb{C}^{m \times n}$，定义一个实值函数 $\| A \|$，它满足以下三个条件：

1. 非负性：当 $A \neq O$时，$\| A \| > 0$；当 $A = O$时，$\| A \| = 0$；

2. 齐次性：$\| \alpha A \| = | \alpha | \| A \|$（$\alpha \in \mathbb{C}$）；

3. 三角不等式：$\| A + B \| \leq \| A \| + \| B \|$（$B \in \mathbb{C}^{m \times n}$）。

则称 $\| A \|$为 $A$的广义矩阵范数。

若对 $\mathbb{C}^{m \times n}$，$\mathbb{C}^{n \times l}$及 $\mathbb{C}^{m \times l}$上的同类广义矩阵范数 $\| \cdot \|$，还应满足下面一个条件：

4. 相容性：$\|AB\| \leq \|A\|\cdot \|B\|, B\in\mathbb{R}^{n\times l}$。

则称 $\|A\|$为 $A$ 的矩阵范数。

### 定义 1.2-相容性

对于 $\mathbb{C}^{m \times n}$上的矩阵范数 $\| \cdot \|_M$和 $\mathbb{C}^m$与 $\mathbb{C}^n$上的同类向量范数 $\| \cdot \|_V$，如果

$$
\| Ax \|_V \leq \| A \|_M \| x \|_V \quad (\forall A \in \mathbb{C}^{m \times n}, \forall x \in \mathbb{C}^n)
$$

则称矩阵范数 $\| \cdot \|_M$与向量范数 $\| \cdot \|_V$是相容的。

### 定理 1.1 相容的矩阵范数与向量范数

设 $\|\cdot\|_{M}$ 是 $\mathbb{C}^{n\times n}$ 上的矩阵范数，任取 $\mathbb{C}^n$ 中的非零列向量 $y$，且固定 $y$，则函数

$$
\|x\|_{V}=\|x y^{H}\|_{M}\quad\left(\forall x\in \mathbb{C}^n\right)
$$

是 $\mathbb{C}^n$ 上的向量范数，且矩阵范数 $\|\cdot\|_M$ 与向量范数 $\|\cdot\|_V$ 相容。

**证明** 

非负性：当 $x\neq 0$ 时，$x y^{H}\neq 0$，从而 $\|x\|_{V}>0$；当 $x = 0$ 时，$x y^{H} = O$，从而 $\|x\|_{V} = 0$。

齐次性：对任意 $k \in \mathbb{C}$，有

$$
\|kx\|_{V} = \|kxy^{H}\|_{M} = |k| \|xy^{H}\|_{M} = |k| \|x\|_{V}
$$

三角不等式：对任意 $x_1, x_2 \in \mathbb{C}^n$，有

$$
\begin{aligned}
    \|x_1 + x_2\|_{V} = \|(x_1 + x_2)y^{H}\|_{M} = \|x_1 y^{H} + x_2 y^{H}\|_{M} \\
    \leq \|x_1 y^{H}\|_{M} + \|x_2 y^{H}\|_{M} = \|x_1\|_{V} + \|x_2\|_{V}
\end{aligned}
$$

因此，$\|x\|_{V}$ 是 $\mathbb{C}^n$ 上的向量范数。当 $A \in \mathbb{C}^{n \times n}, x \in \mathbb{C}^n$ 时

$$
\|Ax\|_{V} = \|(Ax)y^{H}\|_{M} = \|A(xy^{H})\|_{M} \leq \|A\|_{M} \|xy^{H}\|_{M} = \|A\|_{M} \|x\|_{V}
$$

所以矩阵范数 $\|\cdot\|_M$ 与向量范数 $\|\cdot\|_V$ 相容。证毕。$\blacksquare$

## 2. 具体范数

### 2.1 Frobenius 范数

**定义2.1** $\|A\|_F = \sqrt{\sum_{i,j}a_{ij}^2}$.

$\|A\|_F$ 有一特点，现以定理给出于下。

**定理 2.1** 设 $A \in \mathbb{C}^{m \times n}$，且 $P \in \mathbb{C}^{m \times m}$与 $Q \in \mathbb{C}^{n \times n}$都是酉矩阵，则

$$
\|P A\|_F = \|A\|_F = \|A Q\|_F
$$

即给 $A$左乘或右乘以酉矩阵后，其 $\|\cdot\|_F$ 值不变（在 $A \in \mathbb{R}^{m \times n}$时，$P$和 $Q$都是正交矩阵）。

**证明** 若记 $A$的第 $j$列为 $a_j (j=1,2,\cdots, n)$，则有

$$
\begin{align*}
\|P A\|_F^2 &= \left\|P\left(a_1, a_2,\cdots, a_n\right)\right\|_F^2 \\
&= \left\|\left(P a_1, P a_2,\cdots, P a_n\right)\right\|_F^2 \\
&= \sum_{j=1}^n \left\|P a_j\right\|_2^2 \\
&= \sum_{j=1}^n \left\|a_j\right\|_2^2 \\
&= \|A\|_F^2
\end{align*}
$$

即 $\|P A\|_F = \|A\|_F$。于是

$$
\begin{align*}
\|A Q\|_F &= \|(A Q)^H\|_F = \left\|Q^H A^H\right\|_F \\
&= \|A^H\|_F = \|A\|_F
\end{align*}
$$

证毕。$\blacksquare$

### 2.2 诱导范数

**定理 2.2** 已知 $\mathbb{C}^m$ 和 $\mathbb{C}^n$ 上的同类向量范数 $\|\cdot\|$，设 $A \in \mathbb{C}^{m \times n}$，则函数

\[
\|A\| = \max_{\|x\|=1} \|Ax\|
\]

是 $\mathbb{C}^{m \times n}$ 上的矩阵范数，且与已知的向量范数相容。

**证明** 由向量范数是其分量的连续函数的性质可知，对每一个矩阵 $A$ 而言，这个最大值都是可以达到的，也就是说，能够找到这样的向量 $x_0$，使得 $\|x_0\| = 1$，而 $\|Ax_0\| = \|A\|$。

非负性：当 $A \neq O$ 时，存在 $x_0 \in \mathbb{C}^n$ 满足 $\|x_0\| = 1$，使得 $Ax_0 \neq 0$，从而

\[
\|A\| \geqslant \|Ax_0\| > 0
\]

当 $A = O$ 时，$\|A\| = \max_{\|x\|=1} \|Ox\| = 0$。

齐次性：设 $\alpha \in \mathbb{C}$，则有

\[
\|\alpha A\| = \max_{\|x\|=1} \|\alpha A x\| = |\alpha| \max_{\|x\|=1} \|A x\| = |\alpha| \|A\|
\]

三角不等式：设 $B \in \mathbb{C}^{m \times n}$，对于矩阵 $A+B$，存在 $x_1 \in \mathbb{C}^n$ 满足 $\|x_1\|=1$，使得

\[
\|A+B\| = \|(A+B)x_1\|
\]

于是

\[
\begin{align*}
\|A+B\| &= \|Ax_1 + Bx_1\| \leqslant \|Ax_1\| + \|Bx_1\| \leqslant \|A\| + \|B\|
\end{align*}
\]

下面证明，对于任意的 $y \in \mathbb{C}^n$ 及 $A \in \mathbb{C}^{m \times n}$，有

\[
\|Ay\| \leqslant \|A\| \|y\|
\]

当 $y=0$ 时，结论显然成立；当 $y \neq 0$ 时，令 $y_0 = \frac{1}{\|y\|} y$，则

\[
\|y_0\| = 1, \text{且有} \|Ay_0\| \leqslant \|A\|, \text{于是}
\]

\[
\begin{align*}
\|Ay\| &= \|A(\|y\| y_0)\| = \|y\| \|Ay_0\| \leqslant \|y\| \|A\|
\end{align*}
\]

即结论亦成立。

最后证明，对于任意的 $A \in \mathbb{C}^{m \times n}$ 及 $B \in \mathbb{C}^{n \times l}$，有

\[
\|AB\| \leqslant \|A\| \|B\|
\]

对于矩阵 $AB$，存在 $x_2 \in \mathbb{C}^l$ 满足 $\|x_2\|=1$，使得

\[
\|AB\| = \|(AB)x_2\|
\]

利用 $\|Ay\| \leqslant \|A\| \|y\|$，可得

\[
\begin{align*}
\|AB\| &= \|A(Bx_2)\| \leqslant \|A\| \|Bx_2\| \leqslant \|A\| \|B\| \|x_2\| = \|A\| \|B\|
\end{align*}
\]

即 $\|A\|$ 是 $A$ 的矩阵范数。证毕。$\blacksquare$

**定理 2.3** 设 \( A = (a_{ij})_{m \times n} \in \mathbb{C}^{m \times n} \)，\( x = (\xi_1, \xi_2, \ldots, \xi_n)^T \in \mathbb{C}^n \)，则从属于向量 \( x \) 的三种范数 \( \|x\|_1 \)，\( \|x\|_2 \)，\( \|x\|_\infty \) 的矩阵范数依次是：

(1) \( \|A\|_1 = \max_j \sum_{i=1}^m |a_{ij}| \)；

(2) \( \|A\|_2 = \sqrt{\lambda_1} \)，\( \lambda_1 \) 为 \( A^H A \) 的最大特征值；

(3) \( \|A\|_\infty = \max_i \sum_{j=1}^n |a_{ij}| \)。

通常称 \( \|A\|_1 \)，\( \|A\|_2 \) 及 \( \|A\|_\infty \) 依次为列和范数，谱范数及行和范数。

**证明** 

(1) 设 \( \|x\|_1 = 1 \)，则

\[
\begin{aligned}
\|Ax\|_1 = \sum_{i=1}^m \left| \sum_{j=1}^n a_{ij} \xi_j \right| \leq \sum_{i=1}^m \sum_{j=1}^n |a_{ij}| |\xi_j| \\
= \left( \max_j \sum_{i=1}^m |a_{ij}| \right) \sum_{j=1}^n |\xi_j| = \max_j \sum_{i=1}^m |a_{ij}|
\end{aligned}
\]

因此有

\[
\|A\|_1 = \max_{\|x\|_1 = 1}\|Ax\|_1 \leqslant \max_j\sum_{i=1}^{m}|a_{ij}|.
\]

选取 $k$ 使得

\[
\sum_{i=1}^m |a_{ik}| = \max_j\sum_{i=1}^{m}|a_{ij}|.
\]

则 $Ae_k = (a_{1k}, a_{2k}, \ldots, a_{mk})^T$，且

\[
\|A\|_1 = \max_{\|x\|_1=1}\|Ax\|_1 \geqslant \|Ae_k\|_1 = \sum_{i=1}^m |a_{ik}| =  \max_j\sum_{i=1}^{m}|a_{ij}|
\]

因此(1)式成立。

(2) 因为 \( A^H A \) 是 Hermite 矩阵，且由

\[
x^H (A^H A) x = (Ax)^H (Ax) = \|Ax\|_2^2 \geq 0
\]

知 \( A^H A \) 是半正定的，从而它的特征值都是非负实数，设为 \( \lambda_1 \geq \lambda_2 \geq \ldots \geq \lambda_n \geq 0 \)。由于 \( A^H A \) 是 Hermite 矩阵，因此它具有 \( n \) 个互相正交的且 \( L_2 \) 范数为 1 的特征向量 \( x_1, x_2, \ldots, x_n \)，并设它们依次属于特征值 \( \lambda_1, \lambda_2, \ldots, \lambda_n \)。于是，任何一个范数 \( \|x\|_2 = 1 \) 的向量 \( x \)，可以用这些特征向量线性表示，即有

\[
x = \xi_1 x_1 + \xi_2 x_2 + \ldots + \xi_n x_n,
\]

由于

\[
A^HAx = \sum_{i=1}^n \xi_i A^HAx_i = \sum_{i=1}^n \lambda_i \xi_i x_i,
\]

因此有

\[
\|Ax\|_2^2 = (x, A^H A x) = \left( \sum_{i=1}^n \xi_i x_i, \sum_{i=1}^n \lambda_i \xi_i x_i \right) = \sum_{i=1}^n \lambda_i |\xi_i|^2 \leq \lambda_1 \sum_{i=1}^n |\xi_i|^2 = \lambda_1
\]

从而有

\[
\|A\|_2 = \max_{\|x\|_2 = 1} \|Ax\|_2 \leq \sqrt{\lambda_i}.
\]

另一方面，由于 $\|x\|_2=1$，而且 \(\|Ax\|_2^2 = \lambda_1 \)，所以 \( \|A\|_2 = \max_{\|x\|_2=1}\|Ax\|_2 \geq  \sqrt{\lambda_1} \)，因此(2)式成立。

(3) 设 $\|x\|_{\infty}=1$，则

\[
\|A x\|_{\infty}=\max_i|\sum_{j=1}^n a_{ij}\xi_j| \\
\leqslant \max_i\sum_{j=1}^n|a_{ij}||\xi_j| \leqslant \max_i\sum_{j=1}^n|a_{ij}|
\]

从而有 $\|A\|_{\infty}=\max_{\|x\|_{\infty}=1}\|A x\|_{\infty} \leqslant \max_{i}\sum_{j=1}^n\left|a_{ij}\right|$。

选取 $k$，使得

\[
\sum_{j=1}^n\left|a_{kj}\right|=\max_i\sum_{j=1}^n\left|a_{ij}\right|
\]

\[
y=\left[\begin{array}{l}\eta_1\\ \vdots\\ \eta_n\end{array}\right],\quad\eta_j=\left\{\begin{array}{l} 1\quad(a_{kj}=0)\\ \frac{\left|a_{kj}\right|}{a_{kj}}\quad(a_{kj}\neq 0)\end{array}\right.
\]

则有 $\|y\|_{\infty}=1$，且

\[
Ay=(*,\cdots,*,\sum_{j=1}^n\left|a_{kj}\right|,*,\cdots,*)^T
\]

从而

\[
\|A\|_{\infty}=\max_{\|x\|_{\infty}=1}\|A x\|_{\infty} \geqslant \|A y\|_{\infty} \geqslant \sum_{j=1}^n\left|a_{kj}\right| = \max_i\sum_{j=1}^n\left|a_{ij}\right|.
\]

从而(3)式成立。证毕。$\blacksquare$