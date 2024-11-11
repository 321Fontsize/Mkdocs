# Perturbations

## 1. Theorems

### Lemma 1.1

If $F \in \mathbb{R}^{n\times n}$ and $\|F\|_p < 1$, then $I-F$ is nonsingular and 

$$
(I-F)^{-1} = \sum_{k=0}^{\infty}F^k
$$

with

$$
\|(I-F)^{-1}\|_p \leq \frac{1}{1-\|F\|_p}.
$$

***Proof.***

Suppose $I-F$ is singular. It follows that $(I-F)x = 0$ for some nonzero $x$. But then $\|x\|_p = \|Fx\|_p \leq \|F\|_p \|x\|_p$ implies $\|F\|_p \geq 1$, a contradiction. Thus, $I-F$ is nonsingular.

To obtain an expression for its inverse, consider the identity

$$
\left(\sum_{k=0}^N F^k\right)(I-F)=I - F^{N+1}.
$$

Since $\|F\|_p < 1$, it follows that $\lim\limits_{k\rightarrow \infty}F^k = 0$ because $\|F^k\|_p \leq \|F\|_p^k$. Thus, 

$$
\left(\lim_{N\rightarrow \infty}\sum_{k=0}^{N}F^k \right)(I-F) = I,
$$

which implies that $(I-F)^{-1} = \sum_{k=0}^{\infty}F^k$. From this it's easy to show that

$$
\|(I-F)^{-1}\|_p \leq \sum_{k=0}^{\infty}\|F\|_p^k = \frac{1}{1 - \|F\|_p}
$$

completing the proof. $\blacksquare$

**Induction**

If $I + F$ or $I - F$ is singular, then $\|F\|_p \geq 1$.  

### Lemma 1.2

If $A$ is nonsingular and $r\triangleq \|A^{-1}E\|_p < 1$, then $A + E$ is nonsingular and

$$
\|(A+E)^{-1} - A^{-1} \|_p \leq \frac{\|E\|_p \|A^{-1}\|_p^2 }{1-r}.
$$

***Proof.***

Note that $A + E = (I + EA^{-1})A$ with $\|EA^{-1}\|_p < 1$, then by **Lemma 1.1**, we have $I + EA^{-1}$ is nonsingular and $\|(I + EA^{-1})^{-1}\|_p \leq \frac{1}{1-r}$.

Thus $A+E$ is nonsingular and 

$$
\begin{aligned}
\|(A+E)^{-1} - A^{-1} \|_p &= \| A^{-1}\left(A - (A+E)  \right) (A+E)^{-1}\|_p \\
&= \|-A^{-1}EA^{-1}(I+EA^{-1})\|_p \\
&\leq \frac{\|E\|_p \|A^{-1}\|_p^2 }{1-r}.
\end{aligned}
$$

$\blacksquare$

## 2. Applications

### 2.1 Gershgorin Circle Theorem

If $X^{-1}AX = D + F$ where $D = \text{diag}(d_1, \cdots, d_n)$ and $F$ has zero diagonal entries, then 

$$
\lambda(A)\subset \bigcup_{i=1}^{n}D_i,
$$

where $D_i = \{z\in\mathbb{C}: |z - d_i|\leq \sum_{j=1}^n|f_{ij}| \}$.

***Proof.***

Suppose $\lambda\in\lambda(A)$ and assume W.L.O.G. that $\lambda \neq d_i$ for $i=1:n$ (since $d_i \in \cup_{i=1}^n D_i$ trivially).

Since $(D-\lambda I) + F = X^{-1}(A-\lambda I)X$ is singular and $D-\lambda I$ is nonsingular, then $I + (D-\lambda I)^{-1}F = (D-\lambda I)^{-1}(D-\lambda I + F)$ is singular. 

Then by the **Induction from Lemma 1.1**, we have

$$
1 \leq \|(D-\lambda I)^{-1}F\| = \sum_{j=1}^n \frac{f_{kj}}{|d_k -\lambda|},
$$

which completes the proof. $\blacksquare$

### 2.2 Bauer-Fike

If $\mu$ is an eigenvalue of $A+E \in\mathbb{C}^{n\times n}$ and $X^{-1}AX = D = \text{diag}(\lambda_1, \cdots, \lambda_n)$, then 

$$
\mathop{\min}\limits_{\lambda \in \lambda(A)} |\lambda - \mu| \leq \kappa_p(X)\|E\|_p, \quad 1\leq p \leq +\infty.
$$

***Proof.***

If $\mu\in\lambda(A)$, then the theorem is obviously true.

W.L.O.G., assume $\mu\notin\lambda(A)$. Since $A + E -\mu I = X(D+X^{-1}EX - \mu I )X^{-1}$ is singular and $D-\mu I$ is nonsingular, then $I+(D-\mu I)^{-1}(X^{-1}EX) = (D-\mu I)^{-1}(D-\mu I + X^{-1}EX)$ is singular.

Then by the **Induction from Lemma 1.1**, we have

$$
\begin{aligned}
1 &\leq \|(D-\mu I)^{-1}(X^{-1}EX)\|_p\\
&\leq \|(D-\mu I)^{-1}\|_p\cdot  \kappa_p(X)\|E\|_p\\
&= \frac{\kappa_p(X)\|E\|_p}{\mathop{\min}\limits_{\lambda \in \lambda(A)} |\lambda - \mu|},
\end{aligned}
$$

which completes the proof. $\blacksquare$ 