# Dense Set

## 1. Definition

- 在拓扑学及数学的其他相关领域，给定拓扑空间 $X$ 及其子集 $A$，如果对于 $X$ 中任一点 $x$，$x$ 的任一邻域与 $A$ 的交集不为空，则称 $A$ 在 $X$ 中**稠密**.
  - 直观上，如果 $X $ 中的任意一点 $x$ 可以被 $A$ 中的点很好地逼近，则称 $A $ 在 $X $ 中**稠密**.
- 等价地说，$A$ 在 $X$ 中**稠密**当且仅当 $A$ 在 $X$ 中的闭包是 $X$.

## 2. 度量空间中的稠密集

给定度量空间 $(E,d)$，当 $X$ 的拓扑由该度量空间给定时，有

$$
\overline{A} = A\cup \{\lim_n a_n\colon \forall n\geq 0,a_n \in A \},
$$

$\overline{A}$ 称为 $A$ 在 $X$ 中的闭包. 

当 $\overline{A} = X$ 时，$A$ 在 $X$ 中**稠密**，且注意到 $A \subset \{\lim_n a_n\colon \forall n\geq 0,a_n \in A \}$.

## 3. Examples

- 每一拓扑空间是其自身的稠密集.
- 有理数域和无理数域是实数域中的稠密集.