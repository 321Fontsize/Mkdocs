# Perturbations

## 1. 定理

### 引理 1.1

如果$F \in \mathbb{R}^{n\times n}$ 且$\|F\|_p < 1$，则$I-F$ 是非奇异的，且

$$
(I-F)^{-1} = \sum_{k=0}^{\infty}F^k
$$

同时

$$
\|(I-F)^{-1}\|_p \leq \frac{1}{1-\|F\|_p}.
$$

***证明***

假设$I-F$ 是奇异的。这意味着存在某个非零向量$x$ 使得$(I-F)x = 0$。但这样$\|x\|_p = \|Fx\|_p \leq \|F\|_p \|x\|_p$ 意味着$\|F\|_p \geq 1$，这与假设矛盾。因此，$I-F$ 是非奇异的。

为了得到其逆的表达式，考虑恒等式

$$
\left(\sum_{k=0}^N F^k\right)(I-F)=I - F^{N+1}.
$$

由于$\|F\|_p < 1$，可以得出$\lim\limits_{k\rightarrow \infty}F^k = 0$ 因为$\|F^k\|_p \leq \|F\|_p^k$。因此，

$$
\left(\lim_{N\rightarrow \infty}\sum_{k=0}^{N}F^k \right)(I-F) = I,
$$

这意味着$(I-F)^{-1} = \sum_{k=0}^{\infty}F^k$。由此可以很容易地证明

$$
\|(I-F)^{-1}\|_p \leq \sum_{k=0}^{\infty}\|F\|_p^k = \frac{1}{1 - \|F\|_p}
$$

完成了证明.$\blacksquare$


#### 引理 1.1 的推论

若 $I + F$ 或 $I - F$ 奇异，则 $\|F\|_p \geq 1$.  

### 引理 1.2

如果$A$ 是非奇异的且$r\triangleq \|A^{-1}E\|_p < 1$，则$A + E$ 是非奇异的且

$$
\|(A+E)^{-1} - A^{-1} \|_p \leq \frac{\|E\|_p \|A^{-1}\|_p^2 }{1-r}.
$$

***证明***

注意到$A + E = (I + EA^{-1})A$ 且$\|EA^{-1}\|_p < 1$，则由**引理 1.1**，我们有$I + EA^{-1}$ 是非奇异的且$\|(I + EA^{-1})^{-1}\|_p \leq \frac{1}{1-r}$。

因此$A+E$ 是非奇异的且

$$
\begin{aligned}
\|(A+E)^{-1} - A^{-1} \|_p &= \| A^{-1}\left(A - (A+E)  \right) (A+E)^{-1}\|_p \\
&= \|-A^{-1}EA^{-1}(I+EA^{-1})\|_p \\
&\leq \frac{\|E\|_p \|A^{-1}\|_p^2 }{1-r}.
\end{aligned}
$$

$\blacksquare$

## 2. 应用

### 2.1 Gershgorin 圆盘定理

如果$X^{-1}AX = D + F$ 其中$D = \text{diag}(d_1, \cdots, d_n)$ 且$F$ 有零对角线元素，则

$$
\lambda(A)\subset \bigcup_{i=1}^{n}D_i,
$$

其中$D_i = \{z\in\mathbb{C}: |z - d_i|\leq \sum_{j=1}^n|f_{ij}| \}$。

***证明***

假设$\lambda\in\lambda(A)$ 并假设不失一般性地$\lambda \neq d_i$ 对于$i=1:n$（因为$d_i \in \cup_{i=1}^n D_i$ 是显而易见的）。

由于$(D-\lambda I) + F = X^{-1}(A-\lambda I)X$ 是奇异的且$D-\lambda I$ 是非奇异的，则$I + (D-\lambda I)^{-1}F = (D-\lambda I)^{-1}(D-\lambda I + F)$ 是奇异的。

然后由**引理 1.1 的推论**，我们有

$$
1 \leq \|(D-\lambda I)^{-1}F\| = \sum_{j=1}^n \frac{f_{kj}}{|d_k -\lambda|},
$$

这完成了证明.$\blacksquare$

### 2.2 Bauer-Fike 定理

如果$\mu$ 是$A+E \in\mathbb{C}^{n\times n}$ 的特征值且$X^{-1}AX = D = \text{diag}(\lambda_1, \cdots, \lambda_n)$，则

$$
\mathop{\min}\limits_{\lambda \in \lambda(A)} |\lambda - \mu| \leq \kappa_p(X)\|E\|_p, \quad 1\leq p \leq +\infty.
$$

***证明***

如果$\mu\in\lambda(A)$，则定理显然是正确的。

不失一般性，假设$\mu\notin\lambda(A)$。由于$A + E -\mu I = X(D+X^{-1}EX - \mu I )X^{-1}$ 是奇异的且$D-\mu I$ 是非奇异的，则$I+(D-\mu I)^{-1}(X^{-1}EX) = (D-\mu I)^{-1}(D-\mu I + X^{-1}EX)$ 是奇异的。

然后由**引理 1.1 的归纳**，我们有

$$
\begin{aligned}
1 &\leq \|(D-\mu I)^{-1}(X^{-1}EX)\|_p\\
&\leq \|(D-\mu I)^{-1}\|_p\cdot  \kappa_p(X)\|E\|_p\\
&= \frac{\kappa_p(X)\|E\|_p}{\mathop{\min}\limits_{\lambda \in \lambda(A)} |\lambda - \mu|},
\end{aligned}
$$

这完成了证明.$\blacksquare$
