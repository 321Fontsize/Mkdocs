# Ch5 共轭梯度法

## 5.1 最速下降法

亦称梯度下降法

- 定义二次泛函$\psi(x) = x^TAx - 2b^Tx$，若A对称正定，则求解$Ax = b$的解等价于求解二次泛函$\psi(x)$的极小值点。

- 求解思路：给定$x^{(0)}$，依次求$x^{(1)}, x^{(2)}, \cdots$，s.t. $\psi(x^{(k+1)} < \psi(x^{(k)})$

  $\rightarrow$给定$x^{(0)}$，在方向$x^{(0)} + \alpha p^{(0)}$上求$\psi(x)$的极小值点$x^{(1)}$；

  $\rightarrow \cdots$

  - 问题转化为两要点：
    - 如何确定方向序列$\{ p^{(k)}\}$；
    - 如何求得步长$\alpha$。

- 算法：

  先要取$x^{(0)} \in \mathbb{R}^n, p^{(0)} = r^{(0)} = b - Ax^{(0)}$，之后按如下方式做迭代：

$$
\left\{
\begin{aligned}
\alpha_k =& \frac{\langle r^{(k)}, r^{(k)}\rangle}{\langle r^{(k)}, Ar^{(k)}\rangle}\\
x^{(k+1)} =& x^{(k)} + \alpha_k r^{(k)}\\
r^{(k+1)} =& b - Ax^{(k+1)}
\end{aligned}
\right.
$$

## 5.2 共轭梯度法(C-G)

### 5.2.1 与最速下降法作比较

1. 相同点：第一步仍然要取$x^{(0)} \in \mathbb{R}^n, r^{(0)} = Ax^{(0)} - b$；

2. 不同点：在第$k+1(k \geq 1)$步中，”下山“方向不再取负梯度方向。

   1. 而是在过点$x^{(k)}$由向量$r^{(k)}$和$p^{(k-1)}$所张成的平面内选取$\psi(x)$下降最快的方向：

$$
\pi = \{x = x^{(k)} + \xi r^{(k)} + \eta p^{(k-1)}: \xi, \eta \in \mathbb{R} \}
$$

### 5.2.2 算法

我们的目标函数是$f(x) = \frac{1}{2}x^T A x  - b^T x$，故梯度为$(Ax - b)$，则下降方向要取负梯度$p = r = b - Ax$。

初始值：$\forall x^{(0)}, r^{(0)} = b - Ax^{(0)}, p^{(0)} = r^{(0)}, k =0$：

1. $\alpha_k = \frac{\langle r^{(k)}, p^{(k)}\rangle}{\langle p^{(k)}, Ap^{(k)}\rangle}$
2. $x^{(k+1)} = x^{(k)} + \alpha p^{(k)}$
3. $r^{(k+1)} = b - Ax^{(k+1)}$
4. $\beta_k = -\frac{\langle r^{(k+1)}, Ap^{(k)}\rangle}{\langle p^{(k)}, Ap^{(k)}\rangle}$
5. $p^{(k+1)} = r^{(k+1)} + \beta_k p^{(k)}$

第一步示例：
$$
\begin{aligned}
\alpha_0 &= \frac{\langle r^{(0)}, r^{(0)}\rangle}{\langle p^{(0)}, Ap^{(0)}\rangle}, \\
x^{(1)} &= x^{(0)} + \alpha_0 p^{(0)}, \\
r^{(1)} &= b - Ax^{(1)},\\
\beta_1 &= -\frac{\langle r^{1}, Ap^{(0)} \rangle}{\langle p^{(0)}, Ap^{(0)}\rangle},\\
p^{(1)} &= r^{(1)} + \beta_1 p^{(0)},\\
k &= k + 1.
\end{aligned}
$$

### 5.2.3 “共轭”的体现

1. 由C-G法得到的向量组$\{r^{(i)}\}$与$\{p^{(i)} \}$有以下性质：
   1. $\langle p^{(i)}, r^{(j)}\rangle = 0, 0 \leq i < j \leq k$；
   2. $\langle r^{(i)}, r^{(j)} \rangle = 0, 0 \leq i,j \leq k, i\neq j$；
   3. $\langle p^{(i)}, Ap^{(j)} \rangle = 0, 0 \leq i,j \leq k, i\neq j$；
   4. $span\{r^{(0)}, \cdots, r^{(k)} \} = span\{p^{(0)}, \cdots, p^{(k)} \} = \mathcal{K}(A, r^{(0)}, k+1)$，其中$\mathcal{K}(A, r^{(0)}, k+1) = span\{r^{(0)}, Ar^{(0)}, \cdots, A^k r^{(0)} \}$，通常称之为$Krylov$子空间。

[Next-->Ch6](./ch6.md)

[Numerical Algebra's homepage](../NumericalAlgebra.md)