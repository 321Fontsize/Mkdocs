# Ch5 梯度法

## 5.1 最速下降法

亦称梯度下降法

- 定义二次泛函$\psi(x) = x^TAx - 2b^Tx$，若A对称正定，则求解$Ax = b$的解等价于求解二次泛函$\psi(x)$的极小值点。

- 求解思路：给定$x_{0}$，依次求$x_{1}, x_{2}, \cdots$，s.t. $\psi (x_{k+1}) < \psi(x_{k})$

  $\rightarrow$给定$x_{0}$，在方向$x_{0} + \alpha p_{0}$上求$\psi(x)$的极小值点$x_{1}$；

  $\rightarrow \cdots$

  - 问题转化为两要点：
    - 如何确定方向序列$\{ p_{k}\}$；
    - 如何求得步长序列$\{\alpha_k\}$。这里步长的确定一般依靠线搜索方法。

算法：

**Set** $x_{0} \in \mathbb{R}^n, p_{0} = r_{0} = b - Ax_{0}, k = 0,$

**while $r_{k}\neq 0$, repeat**

1. $\alpha_k = \frac{ r_{k}^Tr_{k}}{ r_{k}^TAr_{k}}$

2. $x_{k+1} = x_{k} + \alpha_k r_{k}$

3. $r_{k+1} = b - Ax_{k+1}$

4. $k \leftarrow k+1$

## 5.2 共轭梯度法(Conjugate Gradient)

### 与最速下降法的异同

1. 相同点：第一步仍然要取$x_{0} \in \mathbb{R}^n, r_{0} = b- Ax_{0}$；

2. 不同点：在第$k+1(k \geq 1)$步中，“下山”方向不再取负梯度方向。共轭梯度法的选取策略如下：在过点$x_{k}$由向量$r_{k}$和$p_{k-1}$所张成的平面内选取$\psi(x)$下降最快的方向：$\pi = \{x = x_{k} + \xi r_{k} + \eta p_{k-1}: \xi, \eta \in \mathbb{R} \}$

### 标准共轭梯度法(Standard CG)

目标函数为$f(x) = \frac{1}{2}x^T A x  - b^T x$，故梯度为$\nabla f(x) = Ax - b$，则下降方向要取负梯度$p = r = b - Ax$。

---

算法：

初始值：$\forall x_{0}$, **set** $r_{0} = b - Ax_{0}, p_{0} = r_{0}, k =0$, 

**while $r_{k}\neq 0$, repeat**

1. $\alpha_k = \frac{ r_{k}^Tp_{k}}{ p_{k}^TAp_{k}}$
2. $x_{k+1} = x_{k} + \alpha_{k} p_{k}$
3. $r_{k+1} = b - Ax_{k+1}$
4. $\beta_k = -\frac{ r_{k+1}^TAp_{k}}{ p_{k}^TAp_{k}}$
5. $p_{k+1} = r_{k+1} + \beta_k p_{k}$
5. $k \leftarrow k+1$

---

其中第一步示例如下

$$
\begin{aligned}
&\alpha_0 = \frac{ r_{0}^Tp_{0}}{ p_{0}^TAp_{0}}, \\
&x_{1} = x_{0} + \alpha_0 p_{0}, \\
&r_{1} = b - Ax_{1},\\
&\beta_1 =-\frac{ r_{1}^TAp_{0} }{ p_{0}^TAp_{0}},\\
&p_{1} = r_{1} + \beta_1 p_{0},\\
&k = k + 1.
\end{aligned}
$$

Standard CG算法中部分步骤计算的简化

- $\alpha_k = \frac{r_k^T r_k}{p_k^T A p_k}$
- $r_{k+1} = r_k - \alpha_k A p_k$
- $\beta_{k+1} = \frac{r_{k+1}^T r_{k+1}}{r_k^T r_k}$

### “共轭”的体现

由C-G法得到的向量组$\{r^{(i)}\}$与$\{p^{(i)} \}$有以下性质：
- $ p^{(i)}, r^{(j)} = 0, 0 \leq i < j \leq k$；
- $ r^{(i)}, r^{(j)}  = 0, 0 \leq i,j \leq k, i\neq j$；
- $ p^{(i)}, Ap^{(j)}  = 0, 0 \leq i,j \leq k, i\neq j$；
- $span\{r_{0}, \cdots, r_{k} \} = span\{p_{0}, \cdots, p_{k} \} = \mathcal{K}(A, r_{0}, k+1)$，其中$\mathcal{K}(A, r_{0}, k+1) = span\{r_{0}, Ar_{0}, \cdots, A^k r_{0} \}$，通常称之为$Krylov$子空间。

## 5.3 Preconditioned CG

### 想法来源



[Next-->Ch6](./ch6-非对称特征值问题的计算方法.md)

[Numerical Algebra's homepage](../NumericalAlgebra.md)