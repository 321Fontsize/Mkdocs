# Ch5 梯度法

## 5.1 最优化问题

定义二次泛函 $\phi(x) = x^TAx - 2b^Tx$，若 $A$ 对称正定，则求解 $Ax = b$ 的解等价于求解二次泛函 $\phi(x)$ 的极小值点。

***Proof.***

$\phi(x)$是凸函数，其梯度为 $\nabla \phi(x) = Ax - b$. 设 $x_*$ 为方程 $Ax = b$ 的解，$x_c$ 为函数 $\phi$ 的最小化点，定义 $A$-范数为

$$
\|v\|_A = \sqrt{v^T A v}.
$$

由于

$$
\phi(x_c) = \frac{1}{2} x_c^T A x_c - x_c^T b = \frac{1}{2} (x_c - x_*)^T A (x_c - x_*) - \frac{1}{2} b^T A^{-1} b
$$

且

$$
\phi(x_*) = -\frac{1}{2} b^T A^{-1} b,
$$

因此可得

$$
\phi(x_c) = \frac{1}{2} \|x_c - x_*\|_A^2 + \phi(x_*).
$$

故产生一系列逐次更优的 $\phi$ 的近似最小化点的迭代过程，亦即产生一系列逐次更优的 $Ax = b$ 的近似解的迭代过程，这些近似解以 $A$-范数衡量。$\blacksquare$


## 5.2 最速下降法

亦称梯度下降法

### 分析

改进的近似最小化点 $x_+$ 由下式给出

\[
x_+ = x_c - \mu_cg_c,
\]

其中 $g_c = Ax_c - b$ 是当前梯度，$\mu_c$ 解决 

\[
\min_{\mu\in\mathbb{R}}\phi(x_c - \mu g_c).
\]

$\mu_c$ 的表达式是 $\frac{g_c^Tg_c}{g_c^T A g_c}$

***Proof.***

由于

\[
\begin{aligned}
  \phi(x_c-\mu g_c) &=\frac{1}{2}(x_c-\mu g_c)^T A (x_c-\mu g_c) - (x_c - \mu g_c)^T b\\
  &= \frac{1}{2}\mu^2g_c^TAg_c - \mu g_c^T[Ax_c - b]\\
  &= \frac{1}{2}\mu^2g_c^TAg_c - \mu g_c^Tg_c,
\end{aligned}
\]

且 $\frac{\partial \phi}{\partial \mu} = \mu g_c^T A g_c - g_c^T g_c$，令 $\frac{\partial \phi}{\partial \mu} = 0$，则 $\mu_c = \frac{g_c^T g_c}{g_c^T A g_c}$。$\blacksquare$

现在我们有

\[
\begin{aligned}
\phi(x+) &= \phi(x_c - \mu_c g_c)\\
&= \frac{1}{2}x_c^T Ax_c - x_c^T b + \frac{1}{2}\mu^2 g_c^T A g_c - \mu g_c^T g_c\\
&= \phi(x_c) - \frac{1}{2}\frac{(g_c^T g_c)^2}{g_c^T A g_c}.
\end{aligned}
\]

因此，如果 $g_c \neq 0$，目标函数是减小的。为了建立该方法的全局收敛性，定义

\[
\kappa_c = \frac{g_c^T A g_c}{g_c^T g_c} \cdot \frac{g_c^T A^{-1} g_c}{g_c^T g_c}
\]

并观察到 $g_c^T A^{-1} g_c = 2\phi(x_c) + b^T A^{-1} b$ 和

\[
\phi(x_+) = \phi(x_c) - \frac{1}{2} \frac{1}{\kappa_c} g_c^T A^{-1} g_c = \phi(x_c) - \frac{1}{\kappa_c} \left( \phi(x_c) + \frac{1}{2} b^T A^{-1} b \right).
\]

如果 $\lambda_{\max}(A)$ 和 $\lambda_{\min}(A)$ 是 $A$ 的最大和最小特征值，那么我们有

\[
\kappa_c = \frac{g_c^T A g_c}{g_c^T g_c} \cdot \frac{g_c^T A^{-1} g_c}{g_c^T g_c} \leq \frac{\lambda_{\max}(A)}{\lambda_{\min}(A)} = \kappa_2(A).
\]

如果我们从 (11.3.8) 两边减去 $\phi(x_*) = -(b^T A^{-1} b)/2$ 并使用 (11.3.5)，那么我们得到

\[
\| x_+ - x_* \|_A^2 \leq \left( 1 - \frac{1}{\kappa_2(A)} \right) \| x_c - x_* \|_A^2.
\]

由此通过归纳可知，具有精确线搜索的最速下降法是全局收敛的。

### 算法

**Set** $x_{0} \in \mathbb{R}^n, r_{0} = Ax_0 - b, k = 0$

**while $r_{k}\neq 0$, repeat**

  &emsp;&emsp;$\mu_k = \frac{ r_{k}^Tr_{k}}{ r_{k}^TAr_{k}}$

  &emsp;&emsp;$x_{k+1} = x_{k} - \mu_k r_{k}$

  &emsp;&emsp;$r_{k+1} =  Ax_{k+1} - b$

  &emsp;&emsp;$k \leftarrow k+1$
  
**end(while)**

## 5.3 共轭梯度法(Conjugate Gradient)

### 理论分析 - A subspace strategy

定义 $v + S = \{x: x=v+s, s\in S\}$（我们可以发现最速下降法在第 $k$ 步对仿射空间 $x_k + \text{span}\{\nabla\phi(x_k)\}$ 进行优化）

给定 $x_0: Ax_0 \approx b$，我们的目标是找到 $\begin{cases}S_1\subset S_2\subset S_3\subset\cdots\\ \text{dim}(S_k)=k\end{cases}$ 来解决每一步的 $\min_{x\in x_0 + S_k} \phi(x)$。

如果 $x_k$ 是第 $k$ 步的最小化点，由于嵌套性质，则有

\[
\phi(x_1)\geq \phi(x_2) \geq \cdots \geq \phi(x_n) = \phi(x_*).
\]

接下来我们希望找到一个子空间，以促进 $\phi$ 值的快速减小。由于 $\phi$ 在负梯度方向上减小最快，因此在 $x_k$ 处，我们选择包含 $x_k$ 和 $\nabla\phi(x_k) = Ax_k-b$ 的 $S_k$，则有

\[
\phi(x_{k+1}) \triangleq \min_{x\in x_0 + S_{k+1}}\phi(x) \leq \min_{\mu\in\mathbb{R}}\phi(x_k-\mu g_k)
\]

如果 $x_0$ 是一个初始猜测，并且我们定义 $g_0 = Ax_0 - b$，那么由于 $\nabla\phi(x_k) \in \text{span}\{x_k, Ax_k\}$ 可以得出满足这一要求的唯一方法是设置 

\[
S_k = \mathcal{K}(A,g_0,k) = \text{span}\{g_0, Ag_0, \cdots, A^{k-1}g_0\}.
\]

### 标准共轭梯度法(Standard CG)算法

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

## 5.4 预优共轭梯度法(Preconditioned CG)

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