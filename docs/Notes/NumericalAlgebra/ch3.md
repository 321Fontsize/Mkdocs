# Ch3 最小二乘问题(Least Square)

## 3.1 定义与解的性质

- 定义：

$$
x^* = \mathop{\arg\min}_\limits{y\in \mathbb{R}^n} ||Ay - b||_2\\
\Longrightarrow x^* \in \mathcal{X}_{LS}
$$

- 解的存在唯一性（以下假设$A \in \mathbb{R}^{m\times n}$）

   - $A$的值域$\mathcal{R}(A) = \{y\in \mathbb{R}^m: y=Ax, x\in\mathbb{R}^n \}$；

     > $\mathcal{R}(A) = span(a_1, a_2, \cdots, a_n), a_i$为$A$的列向量

   - $A$的零空间$\mathcal{N}(A) = \{x\in \mathbb{R}^n: Ax=0\}$；

   - $S \in \mathbb{R}^n$，其正交补$S^{\perp} = \{y\in\mathbb{R}^n:y^Tx=0, \forall x \in S \}$。

   - LS问题的解总是存在的。解唯一 $\Longleftrightarrow$ $\mathcal{N}(A) = \{0\}$。

- 求解by正则化方程组/法方程组：将求解$x^*$转化为求解$A^TAx = A^Tb$.

## 3.2 Householder变换

1. 定义$H = I- 2ww^T$，其中$w \in \mathbb{R}^{n\times 1}, ||w||_2=1$
2. $H$的性质：
   1. 对称性：$H^T = H$;
   2. 正交性：$H^TH=I$;
   3. 反射性：$Hx$是$x$关于$w$的垂直超平面（$span\{w\}^{\perp}$）的镜像反射
3. $H$具体形式的求解：
   1. $v = x \pm ||x||_2e_1$;
   2. $w = \frac{v}{||v||_2}$;
   3. $H = I - 2ww^T = I - \frac{2vv^T}{v^Tv} = I - \beta vv^T, \beta=\frac{2}{v^Tv}$

## 3.3 Givens变换

> 亦称为平面旋转变换，可以选择性地将一些元素化为0.

1. G原来是一个单位阵，但第i行第i列与第k列第k行进行了一些操作。$Gy$可以使$y$的某一个分量变为0——利用三角函数性质

2. $G$具体形式的求解关键在于理解

$$
\begin{bmatrix}
cos \;\; sin\\
-sin \;\; cos
\end{bmatrix}
\begin{bmatrix}
a\\
b
\end{bmatrix}=
\begin{bmatrix}
r\\
0
\end{bmatrix}
$$

## 3.4 正交变换法求解LS问题

### <a name='3.4.1'>3.4.1 QR分解定理</a>

- 定理叙述：

   设$A \in \mathbb{R}^{m\times n}(m \geq n)$，则$A$有QR分解：

$$
A = Q\begin{bmatrix}R \\ 0\end{bmatrix},
$$

   其中$Q \in \mathbb{R}^{m\times m}$为正交阵，$R\in \mathbb{R}^{n\times n}$是具有非负对角元的上三角矩阵。且当$m = n$与$A$可逆时，上述分解唯一。

- 求解LS问题：

$$
\begin{aligned}
||Ax-b||_2^2 =& ||Q^TAx - Q^Tb||_2^2\\
=& ||\begin{bmatrix}R \\ 0\end{bmatrix}x - Q^Tb||_2^2\\
=& ||\begin{bmatrix}Rx \\ 0\end{bmatrix} - \begin{bmatrix}c_1 \\ c_2\end{bmatrix}||_2^2\\
&= ||Rx - c_1||_2^2 + ||c_2||_2^2
\end{aligned}
$$

   则$x\in \mathcal{X}_{LS} \Longleftrightarrow Rx=c_1$。

### 3.4.2 利用Householder变换实现QR分解

- Householder变换可以将一个列向量（无论几维）变换为第一个元素非零而其他元素均为0的列向量。

- 有$H_r H_{r-1}\cdots H_2 H_1 A = \begin{bmatrix}R \\ 0\end{bmatrix}$，则$Q = H_1 H_2 \cdots H_{r-1}H_r$.

[Next-->Ch4](./ch4.md)

[Numerical Algebra's homepage](../NumericalAlgebra.md)