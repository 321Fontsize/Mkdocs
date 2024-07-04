# 单纯形法

## 问题概述

待求解的问题具有如下的一般形式：

$$
\begin{aligned}
&\mathop{\min}\limits_{x \in \mathbb{R}^n} c^T x\\
&\text{s.t. } Ax = b, \\
&\quad \;\;\;x \geq 0.
\end{aligned}
$$

## 单纯形算法

给定任意基本指标集合 $\mathcal{B} \subset [n]$ 与非基本指标集合 $\mathcal{N} = [n] -  \mathcal{B}$.

**Repeat**

矩阵 $B$ 由系数矩阵 $A$ 中的部分列向量组成，这些列向量的下标即是指标集 $\mathcal{B}$，矩阵 $N$ 同理. 向量 $c_B$ 是由 $c$ 中的部分元素组成的列向量，$c_N$ 同理.

计算 $x_B = B^{-1}b,\; x_N = \vec{0}$. 

求解 $B^T \lambda = c_B$，计算 $s_N=c_N - N^T \lambda$.

**if** $s_N \geq \vec{0}$

&emsp;&emsp; **stop;**(找到了最优解为 $x_B + x_N$ )

**else**

&emsp;&emsp; **choose** $q\in \mathcal{N}$ randomly **with** $s_q < 0$

&emsp;&emsp; **Solve** $Bd = A_q$

&emsp;&emsp; **if** $d \leq \vec{0}$

&emsp;&emsp;&emsp;&emsp; **stop;**(问题无解)

&emsp;&emsp; **else**

&emsp;&emsp;&emsp;&emsp; **choose** $p = \mathop{\min}\limits_{i\;| \; d_i >0} \frac{(x_B)_i}{d_i}$ $\quad$ (Here, $p\in \mathcal{B}$!!!)

&emsp;&emsp;&emsp;&emsp; **add** $q \rightarrow\mathcal{B}$ and **add** $p\rightarrow\mathcal{N}$.

&emsp;&emsp; **end(if)**

**end(if)**

## Example

$$
\begin{aligned}
&\min -5x_1 - x_2 \\
\text{st.} \quad & x_1 + x_2 \leq 5 \\
& 2x_1 + x_2 = 8 \\
& x_1, x_2 \geq 0
\end{aligned}
$$

(1) 增加非限制变量，使成为标准型
(2) 利用单纯型法求解

解：(1) 原问题可转化为

$$
\begin{aligned}
\min &\quad c^T x \\
\text{st.} \quad & A x = b, x \geq 0
\end{aligned}
$$

其中

$$
\begin{aligned}
&c = [-5, -1, 0, 0]^T, \quad x = [x_1, x_2, x_3, x_4]^T, \\
&A = \begin{bmatrix} 1 & 1 & 1 & 0 \\ 2 & 1 & 0 & 1 \end{bmatrix}, \quad b = \begin{bmatrix} 5 \\ 8 \end{bmatrix}
\end{aligned}
$$

(2) 采用单纯形法求解上述问题如下

a. 初始时，选取$B = \{ 3, 4 \} \Rightarrow N = \{ 1, 2 \}$.

$$
\begin{aligned}
&\Rightarrow x_B = [5, 8]^T, \; x_N = [0, 0]^T \\
&\Rightarrow 初始可行解 x = [0, 0, 5, 8]^T, 且此时目标函数值为0
\end{aligned}
$$

又有$\lambda = B^{-T} c_B = [0, 0]^T, \; S_N = c_N - N^T \lambda = [-5, -1]^T < 0$

b. 选取 $q = 1\in\mathcal{N}$，此时 $A_1 = [1, 2]^T$，计算 $Bd = A_1$，得 $d = [1, 2]^T$。通过计算 $\frac{(x_B)_1}{d_1} = 5, \; \frac{(x_B)_2}{d_2} = 4$，则 $p = 4\in \mathcal{B}$。此时 $\mathcal{B} = \{1, 3\}, \; \mathcal{N} = \{2, 4\}.$ 对应矩阵 $B = \begin{bmatrix} 1 & 1 \\ 2 & 0 \end{bmatrix}, \; N = \begin{bmatrix} 1 & 0 \\ \frac{1}{2} & 1 \end{bmatrix}$，且$c_B = [-5, 0]^T, \; c_N = [-1, 0]^T$

计算得

$$
\begin{aligned}
& x_B = B^{-1} b = [4, 1]^T \; \Rightarrow x_N = [0, 0]^T \\
& \lambda = B^{-T} c_B = [0, -\frac{5}{2}]^T \\
& S_N = c_N - N^T \lambda = [\frac{1}{4}, \frac{5}{2}]^T \Rightarrow \text{stop}
\end{aligned}
$$

此时目标函数值为

$$
c^T [4, 0, 1, 0] = -20.
$$