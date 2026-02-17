# 共形场论与弦论 - 综合学习笔记

> **来源**: 
> - Polchinski, J. *String Theory Vol. 1*, Chapter 2
> - Becker, K., Becker, M., & Schwarz, J.H. *String Theory and M-Theory*, Chapter 3
> - Green, M.B., Schwarz, J.H., & Witten, E. *Superstring Theory Vol. 1*, Chapter 15
> **整理**: Susie 🩷  
> **日期**: 2026-02-17

---

## 目录

1. [引言与背景](#1-引言与背景)
2. [共形群的一般理论 (D维)](#2-共形群的一般理论-d维)
3. [二维共形场论基础](#3-二维共形场论基础)
4. [算符乘积展开 (OPE)](#4-算符乘积展开-ope)
5. [能量-动量张量与Virasoro代数](#5-能量-动量张量与virasoro代数)
6. [顶点算符与物理态](#6-顶点算符与物理态)
7. [Virasoro表示理论](#7-virasoro表示理论)
8. [高级主题：最小模型与模不变性](#8-高级主题最小模型与模不变性)
9. [总结与物理图像](#9-总结与物理图像)

---

## 1. 引言与背景

### 1.1 为什么需要共形场论

弦论的世界面理论具有**局域对称性**：
- **微分同胚不变性** (Diffeomorphism invariance)
- **Weyl对称性** (Weyl symmetry)

固定这些规范后，剩下的理论是**共形不变场论**。

**CFT在弦论中的核心作用**：

| 应用 | 说明 |
|------|------|
| **散射振幅** | 通过顶点算符关联函数计算 |
| **背景场** | 超越平坦26维Minkowski时空 |
| **弦场论** | 场算符创造和消灭整个弦 |
| **紧致化** | 额外维度的几何反映在CFT中 |

### 1.2 Wick旋转与欧几里得世界面

**物理设定**: 弦的世界面原本是Lorentzian签名。

**Wick旋转**: $\tau \rightarrow -i\tau$

结果：
- 世界面度规 $h_{\alpha\beta}$ 变为正定
- 可用复分析工具
- 闵可夫斯基和欧几里得振幅通过解析延拓联系

---

## 2. 共形群的一般理论 (D维)

### 2.1 共形平坦流形

**定义**: D维流形称为**共形平坦**，如果线元可写成：

$$ds^2 = e^{\omega(x)} dx \cdot dx$$

其中：
- $dx \cdot dx = \eta_{\mu\nu}dx^\mu dx^\nu$ (Lorentzian)
- $dx \cdot dx = \delta_{\mu\nu}dx^\mu dx^\nu$ (Euclidean)

### 2.2 D维共形变换

**共形群** = 保持共形平坦性的微分同胚子群

**几何性质**: 
- 保持**角度**
- 扭曲**长度**

### 2.3 共形群的生成元

**无穷小变换**:

$$\delta x^\mu = a^\mu + \omega^\mu_{\,\,\nu} x^\nu + \lambda x^\mu + b^\mu x^2 - 2x^\mu(b \cdot x)$$

| 变换 | 参数 | 含义 |
|------|------|------|
| **平移** | $a^\mu$ (D个) | 空间平移 |
| **旋转** | $\omega^{\mu\nu} = -\omega^{\nu\mu}$ ($\frac{D(D-1)}{2}$个) | Lorentz变换 |
| **尺度** | $\lambda$ (1个) | 整体缩放 |
| **特殊共形** | $b^\mu$ (D个) | 反演-平移-反演 |

**总计**: $\frac{(D+2)(D+1)}{2}$个生成元

### 2.4 特殊共形变换的推导

**技巧**: 反演-平移-反演

**反演变换**:
$$x^\mu \rightarrow \frac{x^\mu}{x^2}$$

它映射：
$$dx \cdot dx \rightarrow \frac{dx \cdot dx}{(x^2)^2}$$

保持度规共形平坦。

**合成变换**: 反演 → 平移 $x^\mu \rightarrow x^\mu + b^\mu$ → 反演

**结果**:
$$x^\mu \rightarrow \frac{x^\mu + b^\mu x^2}{1 + 2b \cdot x + b^2 x^2}$$

**无穷小形式**:
$$\delta x^\mu = b^\mu x^2 - 2x^\mu(b \cdot x)$$

### 2.5 共形群的代数结构

**Lie代数**: 
- Lorentzian: **SO(D, 2)**
- Euclidean: **SO(D+1, 1)**

**关键观察**: D维共形群的生成元数目 = (D+2)维旋转群的生成元数目

**注意**: 
- 当$D > 2$时，上述代数生成整个共形群
- 反演不是无穷小生成的
- 群有两个不连通分支

---

## 3. 二维共形场论基础

### 3.1 二维的特殊性

当$D = 2$时：
- SO(2,2)或SO(3,1)只是**更大代数的子代数**
- 二维共形群是**无穷维的**

**原因**: 所有解析坐标变换都是共形的

$$z \rightarrow f(z), \quad \bar{z} \rightarrow \bar{f}(\bar{z})$$

### 3.2 复坐标与度规

**坐标定义** (Polchinski约定):
$$z = \sigma^1 + i\sigma^2, \quad \bar{z} = \sigma^1 - i\sigma^2$$

**导数**:
$$\partial_z = \frac{1}{2}(\partial_1 - i\partial_2), \quad \partial_{\bar{z}} = \frac{1}{2}(\partial_1 + i\partial_2)$$

**度规**:
$$g_{z\bar{z}} = g_{\bar{z}z} = \frac{1}{2}, \quad g^{z\bar{z}} = g^{\bar{z}z} = 2$$

**积分测度**:
$$d^2z = 2d\sigma^1 d\sigma^2$$

### 3.3 自由标量场的作用量

**Polyakov作用量** (平坦欧几里得度规):

$$S = \frac{1}{4\pi\alpha'} \int d^2\sigma (\partial_1 X^\mu \partial_1 X^\mu + \partial_2 X^\mu \partial_2 X^\mu)$$

**复坐标形式**:

$$\boxed{S = \frac{1}{2\pi\alpha'} \int d^2z \partial X^\mu \bar{\partial} X_\mu}$$

**经典运动方程**:
$$\partial\bar{\partial} X^\mu = 0$$

**解的性质**:
- $\partial X^\mu$ 是全纯的 (holomorphic = 左移)
- $\bar{\partial} X^\mu$ 是反全纯的 (antiholomorphic = 右移)

---

## 4. 算符乘积展开 (OPE)

### 4.1 OPE的基本思想

**核心问题**: 两个局域算符靠近时的行为

**OPE表述**: 

$$\boxed{A_i(z_1, \bar{z}_1) A_j(z_2, \bar{z}_2) = \sum_k c^k_{ij}(z_{12}, \bar{z}_{12}) A_k(z_2, \bar{z}_2)}$$

**关键性质**:
- 在期望值内成立
- 收敛半径 = 到最近其他算符的距离
- 系数函数与期望值中的其他算符无关

**类比**: 类似泰勒级数，但系数可以是奇异的

### 4.2 自由标量场的OPE

**基本传播子**:

$$\langle X^\mu(z_1, \bar{z}_1) X^\nu(z_2, \bar{z}_2) \rangle = -\frac{\alpha'}{2}\eta^{\mu\nu}\ln|z_{12}|^2$$

**完整OPE**:

$$X^\mu(z_1, \bar{z}_1) X^\nu(z_2, \bar{z}_2) = -\frac{\alpha'}{2}\eta^{\mu\nu}\ln|z_{12}|^2 + :X^\nu X^\mu(z_2): + \cdots$$

**导数OPE** (最常用的形式):

$$\boxed{\partial X^\mu(z) \partial X^\nu(w) \sim -\frac{\alpha'}{(z-w)^2}\eta^{\mu\nu}}$$

$$\bar{\partial} X^\mu(\bar{z}) \bar{\partial} X^\nu(\bar{w}) \sim -\frac{\alpha'}{(\bar{z}-\bar{w})^2}\eta^{\mu\nu}$$

**混合导数**:
$$\partial X^\mu(z) \bar{\partial} X^\nu(\bar{w}) \sim 0$$

### 4.3 正规序 (Normal Ordering)

**递归定义**:

$$:X^{\mu_1}(z_1) \cdots X^{\mu_n}(z_n): = X^{\mu_1}(z_1) \cdots X^{\mu_n}(z_n) + \text{(减缩项)}$$

**减缩规则**: 对所有配对方式，每对用$\frac{\alpha'}{2}\eta^{\mu_i\mu_j}\ln|z_{ij}|^2$替代

**紧凑表达式**:
$$:F: = \exp\left(\frac{\alpha'}{4}\int d^2z_1 d^2z_2 \ln|z_{12}|^2 \frac{\delta}{\delta X^\mu(z_1)}\frac{\delta}{\delta X_\mu(z_2)}\right)F$$

### 4.4 格林函数恒等式

**重要公式**:

$$\partial\bar{\partial}\ln|z|^2 = 2\pi\delta^2(z, \bar{z})$$

这是OPE计算的基础。

---

## 5. 能量-动量张量与Virasoro代数

### 5.1 能量-动量张量的定义

**Noether流**: 对于共形变换，守恒流为

$$j(z) = iv(z)T(z), \quad \bar{j}(\bar{z}) = i\bar{v}(\bar{z})\bar{T}(\bar{z})$$

**能量-动量张量**:

$$\boxed{T(z) = -\frac{1}{\alpha'} :\partial X^\mu \partial X_\mu:}$$

$$\boxed{\bar{T}(\bar{z}) = -\frac{1}{\alpha'} :\bar{\partial} X^\mu \bar{\partial} X_\mu:}$$

**性质**:
- $T(z)$ 是全纯的 (权重$(2,0)$)
- $\bar{T}(\bar{z})$ 是反全纯的 (权重$(0,2)$)

### 5.2 TT OPE与中心荷

**基本结果**:

$$\boxed{T(z)T(w) \sim \frac{c/2}{(z-w)^4} + \frac{2T(w)}{(z-w)^2} + \frac{\partial T(w)}{z-w}}$$

对于D个自由标量场：

$$\boxed{c = D}$$

**物理诠释**: 
- 中心荷$c$是**量子反常**
- 标志经典共形对称性的量子破缺
- 弦论中$c = 26$保证一致性

### 5.3 Schwarzian导数

**有限变换**: 在$z \rightarrow w(z)$下，

$$T'(w) = (\partial_w z)^2 \left[T(z) - \frac{c}{12}S(w,z)\right]$$

**Schwarzian导数**:

$$\boxed{S(w,z) = \frac{2(\partial w)(\partial^3 w) - 3(\partial^2 w)^2}{2(\partial w)^2}}$$

**性质**: 当$w$是Möbius变换时，$S(w,z) = 0$

### 5.4 Virasoro生成元

**坐标选择**:

1. **$w$坐标**: $w = \sigma^1 + i\sigma^2$ (圆柱)
   - 时间 = $\text{Im}\, w$的平移

2. **$z$坐标**: $z = \exp(-iw)$ (平面)
   - 时间 = 径向
   - 原点 = 无限过去

**Laurent展开**:

$$T(z) = \sum_{m=-\infty}^{\infty} \frac{L_m}{z^{m+2}}$$

**Virasoro生成元**:

$$\boxed{L_m = \oint_C \frac{dz}{2\pi i} z^{m+1} T(z)}$$

### 5.5 Virasoro代数

**基本对易子**:

$$\boxed{[L_m, L_n] = (m-n)L_{m+n} + \frac{c}{12}(m^3-m)\delta_{m+n,0}}$$

**关键特性**:
- $c = 0$时退化为经典代数
- $m = 0, \pm 1$时中心荷项消失
- $L_0$生成时间平移（径向缩放）
- $L_{\pm 1}$生成特殊共形变换

**Hamilton量**:
$$H = L_0 + \tilde{L}_0 - \frac{c + \tilde{c}}{24}$$

### 5.6 经典Virasoro代数

**生成元**:
$$\ell_n = -z^{n+1}\partial, \quad \bar{\ell}_n = -\bar{z}^{n+1}\bar{\partial}$$

**代数**:
$$[\ell_m, \ell_n] = (m-n)\ell_{m+n}$$

**有限子群**: 
- $\ell_{0,\pm 1}$和$\bar{\ell}_{0,\pm 1}$生成**限制共形群**
- 同构于$SL(2,\mathbb{C})/\mathbb{Z}_2 = SO(3,1)$

---

## 6. 顶点算符与物理态

### 6.1 态-算符对应

**基本思想**: 世界面上每个局域算符对应Hilbert空间中的一个态

**对应关系**:

$$|\Phi\rangle = \lim_{z \rightarrow 0} \Phi(z)|0\rangle$$

**物理图像**: $z = 0$对应欧几里得时间的无限过去

### 6.2 最高权态 (Primary States)

**定义**: 权重为$h$的态满足

$$L_0|h\rangle = h|h\rangle, \quad L_n|h\rangle = 0 \text{ for } n > 0$$

**顶点算符**: 权重$h$的初级场$\Phi(z)$对应最高权态$|h\rangle$

**模式展开**:
$$\Phi(z) = \sum_{n=-\infty}^{\infty} \frac{\Phi_n}{z^{n+h}}$$

**真空条件**:
$$\Phi_n|0\rangle = 0 \text{ for } n > -h, \quad \Phi_{-h}|0\rangle = |\Phi\rangle$$

### 6.3 快子顶点算符

**态**: $|0; k\rangle$，动量为$k^\mu$的基态

**顶点算符**:
$$\boxed{V_k = :e^{ik \cdot X}:}$$

**权重**:
$$h = \tilde{h} = \frac{\alpha'k^2}{4}$$

**物理条件** ($h = \tilde{h} = 1$):
$$\alpha'k^2 = 4 \Rightarrow m^2 = -\frac{4}{\alpha'}$$

这是**快子** (tachyon)，质量平方为负。

### 6.4 矢量顶点算符 (光子)

**极化矢量**: $\zeta_\mu(k)$，满足$k^\mu \zeta_\mu = 0$（横向）

**顶点算符**:
$$\boxed{V_\zeta = \zeta_\mu :\partial X^\mu e^{ik \cdot X}:}$$

**权重**:
$$h = 1 + \frac{\alpha'k^2}{4}$$

**无质量条件** ($h = 1$):
$$k^2 = 0$$

对应**光子** (无质量矢量粒子)。

### 6.5 物理态条件

**开弦条件**:
$$(L_0 - 1)|\phi\rangle = 0, \quad L_n|\phi\rangle = 0 \text{ for } n > 0$$

**闭弦条件**:
$$(L_0 - 1)|\phi\rangle = (\tilde{L}_0 - 1)|\phi\rangle = 0, \quad L_n|\phi\rangle = \tilde{L}_n|\phi\rangle = 0 \text{ for } n > 0$$

**说明**:
- $L_0 - 1 = 0$是**质量壳条件**
- $L_n = 0$ ($n > 0$)是**Virasoro约束**
- 对应权重$h = \tilde{h} = 1$的初级场

---

## 7. Virasoro表示理论

### 7.1 最高权表示

**构造**: 从最高权态$|h\rangle$出发，用 raising operators $L_{-n}$ ($n > 0$)作用

**descendant态**:
$$L_{-k_1} L_{-k_2} \cdots L_{-k_l} |h\rangle, \quad k_1 \geq k_2 \geq \cdots \geq k_l \geq 1$$

**表示**: 最高权态 + 所有descendants = **共形族** (conformal family)

### 7.2 Kac行列式与零态

**内积矩阵**: 在每一能级$N$，定义矩阵

$$M^N_{\{k\},\{k'\}}(c,h) = \langle h, \{k\}|h, \{k'\}\rangle$$

**Kac行列式公式**:

$$\det M^N(c,h) = K_N \prod_{1 \leq rs \leq N} (h - h_{r,s})^{P(N-rs)}$$

**零点位置**:

$$h_{r,s} = \frac{c-1}{24} + \frac{(r\alpha_+ + s\alpha_-)^2}{4}$$

其中：
$$\alpha_\pm = \frac{\sqrt{1-c} \pm \sqrt{25-c}}{\sqrt{24}}$$

### 7.3 酉性条件

**要求**: 所有态的范数必须非负

**允许的$(c,h)$区域**：
1. **连续区域**: $c \geq 1, h \geq 0$
2. **离散点** ($0 \leq c < 1$):
   $$c = 1 - \frac{6}{m(m+1)}, \quad m = 2, 3, \ldots$$
   $$h = h_{r,s} = \frac{[r(m+1) - sm]^2 - 1}{4m(m+1)}$$
   其中$1 \leq r \leq m-1, 1 \leq s \leq m$

**弦论**: $c = 26$，在连续区域，需要$\alpha'k^2 \leq 0$保证$h \geq 0$

### 7.4 零态与退化表示

**零态**: 与所有态正交的态（在Kac行列式零点处出现）

**性质**:
- 零态也是最高权态
- 形成Verma模的不变子空间
- 可"mod out"得到不可约表示

**例子**: 在$h = 0$，总有$L_{-1}|0\rangle = 0$

---

## 8. 高级主题：最小模型与模不变性

### 8.1 最小模型

**定义**: 所有共形族都是退化的CFT

**条件**: 当$\alpha_-/\alpha_+ = -p/q$为有理数时，算符代数在有限集合上闭合

**参数**:
$$c = 1 - \frac{6(p-q)^2}{pq}$$
$$h_{r,s} = \frac{(rq-sp)^2 - (p-q)^2}{4pq}$$
其中$1 \leq r \leq p-1, 1 \leq s \leq q-1$

**酉性**: 只有当$q = p+1$时才是酉的

**例子**: $m = 3$时$c = 1/2$，对应**自由费米子**

### 8.2 融合规则 (Fusion Rules)

**基本规则**:
$$O_{r_1,s_1} O_{r_2,s_2} = \sum_{r,s} [O_{r,s}]$$

其中：
- $r = |r_1-r_2|+1, |r_1-r_2|+3, \ldots, \min(r_1+r_2-1, 2p-1-r_1-r_2)$
- $s = |s_1-s_2|+1, |s_1-s_2|+3, \ldots, \min(s_1+s_2-1, 2q-1-s_1-s_2)$

**物理意义**: OPE中出现的共形族

### 8.3 简单电流与离散对称性

**简单电流**: $J$满足$J \cdot [O_i] = [O_{i'}]$（单一项）

**性质**:
- 有确定的monodromy
- 产生$\mathbb{Z}_N$对称性
- OPE系数满足选择定则

**例子**: 在最小模型中，$O_{p-1,1}$是简单电流，产生$\mathbb{Z}_2$对称性

---

## 9. 总结与物理图像

### 9.1 关键公式汇总

| 概念 | 公式 |
|------|------|
| **复坐标** | $z = \sigma^1 + i\sigma^2$ |
| **作用量** | $S = \frac{1}{2\pi\alpha'}\int d^2z \partial X^\mu \bar{\partial} X_\mu$ |
| **基本OPE** | $\partial X^\mu(z)\partial X^\nu(w) \sim -\frac{\alpha'}{(z-w)^2}\eta^{\mu\nu}$ |
| **应力张量** | $T(z) = -\frac{1}{\alpha'}:\partial X^\mu \partial X_\mu:$ |
| **TT OPE** | $T(z)T(w) \sim \frac{c/2}{(z-w)^4} + \frac{2T(w)}{(z-w)^2} + \frac{\partial T(w)}{z-w}$ |
| **Virasoro代数** | $[L_m,L_n] = (m-n)L_{m+n} + \frac{c}{12}(m^3-m)\delta_{m+n,0}$ |
| **中心荷** | $c = D$ (自由标量) |
| **快子顶点** | $V_k = :e^{ik\cdot X}:$, $h = \frac{\alpha'k^2}{4}$ |
| **物理态条件** | $(L_0-1)|\phi\rangle = L_n|\phi\rangle = 0$ ($n > 0$) |
| **临界维数** | $c = D = 26$ |

### 9.2 物理图像

**世界面 = Riemann面**:
- 弦的世界面是二维曲面
- Wick旋转后可用复分析
- 时间是径向的，原点是无限过去

**CFT = 固定规范后的弦论**:
- 微分同胚 + Weyl对称性被固定
- 剩下共形对称性
- 无穷维对称性使理论可解

**态-算符对应**:
- 弦的每个量子态 ↔ 世界面上的局域算符
- 散射振幅 = 顶点算符的关联函数
- 通过OPE和共形对称性可精确计算

**临界维数 D = 26**:
- 保证异常抵消
- 确保Lorentz不变性
- 是协变量子化的必要条件

**谱与约束**:
- 快子: $m^2 = -4/\alpha'$ (不稳定)
- 光子: $m^2 = 0$ (无质量矢量)
- 更高激发态: 有质量粒子

### 9.3 后续发展

基于此框架，弦论继续发展：
1. **相互作用**: 通过世界面拓扑变化（管状图）
2. **超弦**: 加入费米场，$c = 10$，消除快子
3. **紧致化**: 额外维度的几何反映在CFT中
4. **对偶性**: 不同几何描述对应相同CFT
5. **AdS/CFT**: 引力/规范对偶的精确表述

---

## 参考文献

1. **Polchinski, J.** (1998). *String Theory, Vol. 1: An Introduction to the Bosonic String*. Cambridge University Press. (Chapter 2 - 基础CFT)

2. **Becker, K., Becker, M., & Schwarz, J. H.** (2007). *String Theory and M-Theory: A Modern Introduction*. Cambridge University Press. (Chapter 3 - 共形群与物理应用)

3. **Green, M. B., Schwarz, J. H., & Witten, E.** (1987). *Superstring Theory, Vol. 1: Introduction*. Cambridge University Press. (Chapter 15 - 高级CFT与表示理论)

4. **Di Francesco, P., Mathieu, P., & Sénéchal, D.** (1997). *Conformal Field Theory*. Springer. (全面的CFT教科书)

5. **Ginsparg, P.** (1988). "Applied Conformal Field Theory." *arXiv:hep-th/9108028*. (经典的CFT讲义)

---

*综合三份经典弦论教材的CFT章节*  
*完整覆盖：基础数学、物理应用、高级主题 🩷*