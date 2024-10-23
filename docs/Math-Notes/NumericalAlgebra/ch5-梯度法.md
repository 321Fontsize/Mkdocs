# Ch5 梯度法

## 5.1 最速下降法

亦称梯度下降法

### 理论分析

- 定义二次泛函$\psi(x) = x^TAx - 2b^Tx$，若A对称正定，则求解$Ax = b$的解等价于求解二次泛函$\psi(x)$的极小值点。

- 求解思路：给定任意向量$x_{0}$，依次求$x_{1}, x_{2}, \cdots$，s.t. $\psi (x_{k+1}) < \psi(x_{k})$

  $\rightarrow$给定$x_{0}$，确定一个下降方向 $p_0$，在直线$x_{0} + \alpha p_{0}$上求$\psi(x)$的极小值点$x_{1} = x_0 + \alpha_0 p_0$，i.e.，确定 $\alpha_0 = \mathop{\arg \min }\limits_{\alpha \in \mathbb{R}}\psi(x_0 + \alpha p_0)$；

  $\rightarrow$给定 $x_1$，确定一个下降方向 $p_1$，$\cdots$

- 问题转化为两要点：
  - 如何确定方向序列$\{ p_{k}\}$；
  - 如何求得步长序列$\{\alpha_k\}$。这里步长的确定一般依靠线搜索方法。

### 算法

**Set** $x_{0} \in \mathbb{R}^n, p_{0} = r_{0} = b - Ax_{0}, k = 0,$

**while $r_{k}\neq 0$, repeat**

  &emsp;&emsp;$\alpha_k = \frac{ r_{k}^Tr_{k}}{ r_{k}^TAr_{k}}$

  &emsp;&emsp;$x_{k+1} = x_{k} + \alpha_k r_{k}$

  &emsp;&emsp;$r_{k+1} = b - Ax_{k+1}$

  &emsp;&emsp;$k \leftarrow k+1$
  
**end(while)**

## 5.2 共轭梯度法(Conjugate Gradient)

### 与最速下降法的异同

- 相同点：第一步仍然要取$x_{0} \in \mathbb{R}^n, r_{0} = b- Ax_{0}$；
- 不同点：在第$k+1(k \geq 1)$步中，“下降”方向不再取负梯度方向。CG的下降方向的选取策略如下：在过点$x_{k}$由向量$r_{k}$和$p_{k-1}$所张成的平面内选取$\psi(x)$下降最快的方向：$\pi = \{x = x_{k} + \xi r_{k} + \eta p_{k-1}: \xi, \eta \in \mathbb{R} \}$

### 标准共轭梯度法(Standard CG)

目标函数为$f(x) = \frac{1}{2}x^T A x  - b^T x$，故梯度为$\nabla f(x) = Ax - b$，则下降方向要取负梯度$p = r = b - Ax$。

### 算法

$\forall x_{0}$, **set** $r_{0} = b - Ax_{0}, p_{0} = r_{0}, k =0,$

**while $r_{k}\neq 0$, repeat**

&emsp;&emsp;$\alpha_k = \frac{ r_{k}^Tp_{k}}{ p_{k}^TAp_{k}}$

&emsp;&emsp;$x_{k+1} = x_{k} + \alpha_{k} p_{k}$

&emsp;&emsp;$r_{k+1} = b - Ax_{k+1}$

&emsp;&emsp;$\beta_k = -\frac{ r_{k+1}^TAp_{k}}{ p_{k}^TAp_{k}}$

&emsp;&emsp;$p_{k+1} = r_{k+1} + \beta_k p_{k}$

&emsp;&emsp;$k \leftarrow k+1$

**end(while)**


其中第一步示例如下.

$$
\begin{aligned}
&\alpha_0 = \frac{ r_{0}^Tp_{0}}{ p_{0}^TAp_{0}}, \\
&x_{1} = x_{0} + \alpha_0 p_{0}, \\
&r_{1} = b - Ax_{1},\\
&\beta_0 =-\frac{ r_{1}^TAp_{0} }{ p_{0}^TAp_{0}},\\
&p_{1} = r_{1} + \beta_0 p_{0},\\
&k = k + 1.
\end{aligned}
$$

Standard CG算法中部分步骤计算的简化

- $\alpha_k = \frac{r_k^T r_k}{p_k^T A p_k}$
- $r_{k+1} = r_k - \alpha_k A p_k$
- $\beta_{k} = \frac{r_{k+1}^T r_{k+1}}{r_k^T r_k}$

### “共轭”的体现-共轭基本性质

由C-G法得到的向量组$\{r^{(i)}\}$与$\{p^{(i)} \}$有以下性质：

(1) $p_{i}^T r_{j} = 0,\; 0 \leq i < j \leq k$；

(2) $r_{i}^T r_{j}  = 0,\; 0 \leq i,j \leq k,\; i\neq j$；

(3) $p_{i}^T A p_{j}  = 0,\; 0 \leq i,j \leq k,\;  i\neq j$；

(4)$\text{span}\{r_{0}, \cdots, r_{k} \} = \text{span}\{p_{0}, \cdots, p_{k} \} = \mathcal{K}(A, r_{0}, k+1)$，其中$\mathcal{K}(A, r_{0}, k+1) = \text{span}\{r_{0}, Ar_{0}, \cdots, A^k r_{0} \}$，通常称之为$Krylov$子空间。

**证明** 用数学归纳法. 当 $k=1$时，因为

$$
\begin{aligned}
&p_0 = r_0,\quad r_1 = r_0 - \alpha_0 A p_0,\quad p_1 = r_1 + \beta_1 p_0,
\end{aligned}
$$

故

$$
\begin{aligned}
&r_1^T r_0 = r_0^T r_0 - \alpha_0 r_0^T A r_0 = r_0^T r_0 - \frac{r_0^T r_0}{r_0^T A r_0}r_0^T A r_0 =  0,\\
&p_1^T A p_0 = (r_1 + \beta_0 r_0)^T A r_0 = r_1^T A r_0 - \frac{r_1^T A r_0}{r_0^T A r_0}r_0^T A r_0 = 0.
\end{aligned}
$$

所以定理的结论对 $k=1$ 成立. 现假设定理的结论对 $k(k\geq 1)$ 成立，下面证明其对 $k+1$ 也成立.

(1) 利用等式 $r_{k+1} = r_k - \alpha_k A p_k$ 及归纳法假设，有

$$
p_i^T r_{k+1} = p_i^T r_k  - \alpha_k p_i^T A p_k = 0 - 0 = 0,\quad 0\leq i \leq k-1.
$$

又由于

$$
p_{k}^T r_{k+1} = p_k^T r_k - \alpha_k p_k^T A p_k = p_k^T r_k - \frac{p_k^T r_k}{p_k^T A p_k} p_k^T A p_k = 0.
$$

故定理的结论(1)对 $k+1$ 亦成立.

(2) 利用归纳假设有

$$
\text{span}\{r_0, \cdots, r_k\} = \text{span}\{p_0,\cdots, p_k\},
$$

而由(1)所证可知，$r_{k+1}$与上述子空间正交，从而定理的结论(2)对 $k+1$ 也成立.

(3) 由等式

$$
r_{i+1} = r_i - \alpha_i A p_i, \quad i = 0, 1,\cdots, k-1.
$$

有

$$
p_i^T A = \frac{1}{\alpha_i}(r_{i} - r_{i+1})^T
$$

利用上述等式、等式 $p_{k+1} = r_{k+1} + \beta_k p_k$ 与归纳假设以及(2)所证的结论，有

$$
\begin{aligned}
p_i^T A p_{k+1} &= \frac{1}{\alpha_i}r_{k+1}^T(r_i - r_{i+1}) + \beta_k p_i^T A p_k \\
&= 0 - 0  + 0 = 0, \quad i = 0, 1, \cdots, k-1.
\end{aligned}
$$

而由 $\beta_k$ 的定义可得

$$
\begin{aligned}
p_{k+1}^T A p_k &= (r_{k+1} + \beta_{k+1}p_k)^T A p_k\\
&= r_{k+1}^T A p_k - \frac{r_{k+1}^T A p_k}{p_k^T A p_k} p_k^T A p_k  = 0 - 0 = 0. 
\end{aligned}
$$

故定理的结论(3)对 $k+1$ 也成立.

(4) 由归纳假设知

$$
r_k, \; p_k \in \mathcal{k}(A, r_0, k+1) = \text{span}\{r_0, Ar_0, \cdots, A^kr_0\},
$$

于是

$$
\begin{aligned}
& r_{k+1} = r_k - \alpha_k p_k \in \mathcal{K}(A, r_0, k+2) = \text{span}\{r_0, Ar_0, \cdots, A^{k+1}r_0\},\\
& p_{k+1} = r_{k+1} + \beta_{k}p_k \in \mathcal{K}(A, r_0, k+2) = \text{span}\{r_0, Ar_0, \cdots, A^{k+1}r_0\}.
\end{aligned}
$$

再注意到 (2) 和 (3) 所证的结论表明，向量组 $r_0, \cdots, r_{k+1}$ 是线性无关的，$p_0, \cdots, p_{k+1}$ 是线性无关的，因此定理的结论 (4) 对 $k+1$
同样成立.

综上所述，由归纳法原理知定理得证.

## 5.3 预优共轭梯度法(Preconditioned CG)

### 想法来源

当线性方程组 $Ax = b$ 的系数矩阵 $A$ 仅有几个少数几个互不相同的特征值或者非常良态的时候，共轭梯度法就会收敛得非常之快. 这就促使我们在应用共轭梯度法时，首先应将原方程转化为另一个方程 $\tilde{A}\tilde{x} = \tilde{b}$，其中 $\tilde{A}$ 满足上述讨论的性质.

这个想法的基础是下面这个定理：

用共轭梯度法求得的 $x_k$ 有如下的误差估计：

$$
||x_k - x^* ||_A \leq 2(\frac{\sqrt{\kappa_2(A)} - 1}{\sqrt{\kappa_2(A)} + 1})^k ||x_0 - x^* ||_A,
$$

其中 $||x||_A = \sqrt{x^TA x}, \;\kappa_2(A) = ||A||_2 ||A^{-1}||_2$.

### 理论分析

选择对称正定阵 $C$，使得 $\tilde{A} = C^{-T} A C^{-1}, \tilde{x} = Cx, \tilde{b} = C^{-T}b$，且我们希望 $\tilde{A}$ 具有良好的性质. 

应用Standard CG于上述方程组，常规的做法是首先计算 $\tilde{A} = C^{-T} A C^{-1}, \tilde{x} = Cx, \tilde{b} = C^{-T}b$，并且在获得近似解 $\tilde{x}_k$ 之后需通过变换 $x_k = C^{-1}\tilde{x}_k$ 获取原方程的解. 

实际上这些都是不必要的

令

$$
\forall x_0 \in \mathbb{R}^n,\; \tilde{x}_0 = C x_0, \;\;\text{Set}\;\; \tilde{r}_0 = \tilde{b} - \tilde{A}\tilde{x}_0 = C^{-T}r_0,\; \tilde{p}_0 = \tilde{r}_0,
$$

故

$$
\alpha_0 = \frac{r_0^T C^{-1} C^{-T}r_0}{r_0^T C^{-1}C^{-T} A C^{-1}C^{-T}r_0} = \frac{r^T_0 y_0}{p_0 A^T p_0},
$$

其中 $My_0 = r_0, p_0 = y_0$. 则

$$
\begin{aligned}
y_0 &= M^{-1}r_0\\
&= C^{-1}C^{-T}r_0,\\
&\Longrightarrow M = C^TC = C^2.
\end{aligned}
$$

在获得预优矩阵 $M$ 的表达式之后，我们可以设计如下的算法

### 算法

$\forall x_{0}$,  **set** $r_{0} = b - Ax_{0},\; p_{0} = y_{0} = M^{-1}r_0,\; k =0,$

**while $r_{k}\neq 0$, repeat**

&emsp;&emsp;$\alpha_k = \frac{ r_{k}^Ty_{k}}{ p_{k}^TAp_{k}}$

&emsp;&emsp;$x_{k+1} = x_{k} + \alpha_{k} p_{k}$

&emsp;&emsp;$r_{k+1} = b - Ax_{k+1} = r_k - \alpha_k A p_k$

&emsp;&emsp;**Solve** $My_{k+1}  = r_{k+1}$

&emsp;&emsp;$\beta_k = \frac{ r_{k+1}^T y_{k+1}}{ r_{k}^T y_{k}}$

&emsp;&emsp;$p_{k+1} = y_{k+1} + \beta_k p_{k}$

&emsp;&emsp;$k \leftarrow k+1$

**end(while)**

[Homepage of NLA](./index.md)