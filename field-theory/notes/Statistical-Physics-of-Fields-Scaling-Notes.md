# Statistical Physics of Fields - Chapter 4: The Scaling Hypothesis

> **来源**: Statistical Physics of Fields (Chapter 4)  
> **主题**: 标度假设与重整化群理论  
> **阅读日期**: 2026-02-07  
> **阅读者**: rimo  
> **笔记整理**: Susie 🩷

---

## 1. 引言：为什么需要标度假设？

### 1.1 临界现象的复杂性

在连续相变（continuous phase transition）附近，热力学量表现出**奇异性**（singular behavior）。我们用**临界指数**（critical exponents）$\alpha, \beta, \gamma, \delta, \nu, \eta, \ldots$ 来描述这些奇异行为。

**问题**: 鞍点近似（saddle-point approximation）给出的临界指数值由于**涨落**（fluctuations）的重要性而**不可靠**。

**关键观察**: 各种热力学量是相互关联的，因此临界指数**不能相互独立**。

**标度假设的目标**: 
1. 发现临界指数之间的**关系**
2. 找到描述临界点所需的**最小独立指数数目**

### 1.2 临界点附近的热力学

在 $(t, h)$ 平面中（$t = (T-T_c)/T_c$ 是约化温度，$h$ 是外场），临界点位于 $t = h = 0$。

**共存线**（coexistence line）: $t < 0$ 且 $h = 0$，终止于临界点。

**自由能**（在鞍点近似下）：

$$
f(t, h) = \min_m \left[ \frac{t}{2}m^2 + \frac{u}{4}m^4 - hm \right] = 
\begin{cases}
-\frac{t^2}{16u} & \text{for } h=0, t\lt0 \\
-\frac{3}{4^{4/3}} \frac{h^{4/3}}{u^{1/3}} & \text{for } h=0, t=0
\end{cases}
$$

---

## 2. 齐次性假设（The Homogeneity Assumption）

### 2.1 自由能的齐次函数形式

**基本假设**: 自由能的奇异部分可以用**齐次函数**（homogeneous function）描述。

**定义**: 函数 $f(x_1, x_2, \ldots)$ 是齐次的，如果对于任意标度因子 $b$:

$$
f(b^{p_1}x_1, b^{p_2}x_2, \ldots) = b^p f(x_1, x_2, \ldots)
$$

**自由能的标度形式**:

$$
\boxed{f(t, h) = t^{2-\alpha} g_f\left(\frac{h}{t^\Delta}\right)}
$$

其中：
- $\Delta$ 称为**间隙指数**（gap exponent）
- $g_f$ 是仅依赖于组合变量 $x \equiv h/t^\Delta$ 的标度函数

### 2.2 从鞍点结果确定间隙指数

比较鞍点近似结果和标度形式：

**$h = 0$ 极限**: 需要 $\lim_{x \to 0} g_f(x) \sim \text{const}$

**$t = 0$ 极限**: 需要 $\lim_{x \to \infty} g_f(x) \sim x^{4/3}/u^{1/3}$

后者给出：

$$
f \sim t^{2-\alpha} \frac{(h/t^\Delta)^{4/3}}{u^{1/3}} = \frac{h^{4/3}}{u^{1/3}} t^{2-\alpha-4\Delta/3}
$$

由于沿 $t = 0$ 自由能不能依赖于 $t$，必须有：

$$
2 - \alpha - \frac{4\Delta}{3} = 0 \quad \Rightarrow \quad \boxed{\Delta = \frac{3(2-\alpha)}{4}}
$$

对于鞍点近似（$\alpha = 0$），得到 $\Delta = 3/2$。

### 2.3 齐次性假设的广义形式

**标度假设**: 超越鞍点近似，奇异部分仍保持齐次形式：

$$
f_{\text{sing}}(t, h) = t^{2-\alpha} g_f\left(\frac{h}{t^\Delta}\right)
$$

其中实际的指数 $\alpha$ 和 $\Delta$ 取决于所考虑的临界点。

**选择 $t$ 依赖的原因**: 重现 $h = 0$ 处的热容奇异性 $C_{\text{sing}} \sim t^{-\alpha}$。

---

## 3. 热力学量的导数与标度形式

### 3.1 能量密度

从自由能计算奇异部分能量：

$$
E_{\text{sing}} \sim \frac{\partial f}{\partial t} \sim t^{1-\alpha} \left[ (2-\alpha)g_f(x) - \Delta x g_f'(x) \right] \equiv t^{1-\alpha} g_E(x)
$$

其中 $x = h/t^\Delta$。

**结论**: 齐次函数的导数仍是齐次函数。

### 3.2 热容

二阶导数给出热容：

$$
C_{\text{sing}} \sim -\frac{\partial^2 f}{\partial t^2} \sim t^{-\alpha} g_C(x)
$$

重现 $C_{\text{sing}} \sim t^{-\alpha}$（当 $h \to 0$ 时）。

### 3.3 $t > 0$ 和 $t < 0$ 的统一描述

**问题**: 能否对 $t > 0$ 和 $t < 0$ 使用不同的函数和指数？

**回答**: **不能**。原因如下：

假设：

$$
C_\pm(t, h) = |t|^{-\alpha_\pm} g_\pm\left(\frac{h}{|t|^{\Delta_\pm}}\right)
$$

**解析性条件**: 在 $t = 0$ 且 $h \neq 0$ 处，$C$ 必须是解析的，可展开为泰勒级数：

$$
C(t, h) = c_0(h) + t c_1(h) + t^2 c_2(h) + \cdots
$$

从标度形式展开：

$$
C_\pm = |t|^{-\alpha_\pm} \left[ A_\pm \left(\frac{h}{|t|^{\Delta_\pm}}\right)^{p_\pm} + B_\pm \left(\frac{h}{|t|^{\Delta_\pm}}\right)^{q_\pm} + \cdots \right]
$$

匹配泰勒级数要求：
- $p_\pm \Delta_\pm = -\alpha_\pm$
- $q_\pm \Delta_\pm = -1 + \alpha_\pm$

这导致：

$$
C_\pm(t, h) = A_\pm h^{-\alpha_\pm/\Delta_\pm} + B_\pm h^{-1+\alpha_\pm/\Delta_\pm} t + \cdots
$$

**连续性条件**: 在 $t = 0$ 处连续要求：

$$
\frac{\alpha_+}{\Delta_+} = \frac{\alpha_-}{\Delta_-}, \quad 1 + \frac{\alpha_+}{\Delta_+} = 1 + \frac{\alpha_-}{\Delta_-}
$$

因此：

$$
\boxed{\alpha_+ = \alpha_- \equiv \alpha, \quad \Delta_+ = \Delta_- \equiv \Delta}
$$

**重要结论**: 奇异部分的指数在临界点上下必须相同。

---

## 4. 其他热力学量的标度形式

### 4.1 磁化强度

$$
m(t, h) \sim \frac{\partial f}{\partial h} \sim t^{2-\alpha-\Delta} g_m\left(\frac{h}{t^\Delta}\right)
$$

**极限行为**：

1. **$x \to 0$ 极限**（$h = 0$）:

$$
m(t, h=0) \sim t^{2-\alpha-\Delta} \quad \Rightarrow \quad \boxed{\beta = 2 - \alpha - \Delta}
$$

2. **$x \to \infty$ 极限**（$t = 0$）:

若 $g_m(x) \sim x^p$，则：

$$
m(t=0, h) \sim t^{2-\alpha-\Delta} \left(\frac{h}{t^\Delta}\right)^p
$$

由于此极限不能依赖于 $t$，必须有 $p\Delta = 2 - \alpha - \Delta = \beta$。

因此：

$$
m(t=0, h) \sim h^{(2-\alpha-\Delta)/\Delta} = h^{\beta/\Delta}
$$

定义 $\delta$ 通过 $m \sim h^{1/\delta}$，得到：

$$
\boxed{\frac{1}{\delta} = \frac{\beta}{\Delta} = \frac{2-\alpha-\Delta}{\Delta}}
$$

### 4.2 磁化率

$$
\chi(t, h) \sim \frac{\partial m}{\partial h} \sim t^{2-\alpha-2\Delta} g_\chi\left(\frac{h}{t^\Delta}\right)
$$

在 $h = 0$：

$$
\chi(t, h=0) \sim t^{2-\alpha-2\Delta} \quad \Rightarrow \quad \boxed{\gamma = 2\Delta - 2 + \alpha}
$$

---

## 5. 临界指数恒等式（Exponent Identities）

### 5.1 Rushbrooke 恒等式

从上述关系：

$$
\alpha + 2\beta + \gamma = \alpha + 2(2-\alpha-\Delta) + (2\Delta-2+\alpha) = 2
$$

$$
\boxed{\alpha + 2\beta + \gamma = 2} \quad \text{[Rushbrooke's identity]}
$$

### 5.2 Widom 恒等式

从 $\delta$ 的定义：

$$
\delta = \frac{\Delta}{\beta} = \frac{2\Delta - 2 + \alpha}{2 - \alpha - \Delta} = \frac{\gamma}{\beta}
$$

$$
\boxed{\beta(\delta + 1) = 2 - \alpha + \gamma} \quad \text{[Widom's scaling relation]}
$$

或等价地：

$$
\beta\delta = \beta + \gamma
$$

### 5.3 指数独立性

**关键结论**: 所有（体）临界指数可以用**两个独立指数**表示，例如 $\alpha$ 和 $\Delta$（或 $\alpha$ 和 $\nu$ 等）。

**实验验证**: 这些恒等式与可靠的实验数据以及 $d=3$ 的理论估计（Ising模型、XY模型、Heisenberg模型）以及 $d=2$ 的精确解完全一致。

**临界指数表**（$d=3$）:

| 模型 | $\alpha$ | $\beta$ | $\gamma$ | $\delta$ | $\nu$ | $\eta$ |
|------|-----------|---------|-----------|-----------|--------|---------|
| $n=1$ (Ising) | 0.11 | 0.32 | 1.24 | 4.9 | 0.63 | 0.04 |
| $n=2$ (XY) | -0.01 | 0.35 | 1.32 | 4.7 | 0.67 | 0.04 |
| $n=3$ (Heisenberg) | -0.11 | 0.36 | 1.39 | 4.9 | 0.70 | 0.04 |
| $d=2$, $n=1$ (精确) | 0 | 1/8 | 7/4 | 15 | 1 | 1/4 |

---

## 6. 相关长度的发散（Divergence of Correlation Length）

### 6.1 相关长度的标度

齐次性假设涉及自由能及其导出量，但**不涉及关联函数**。为了得到包含相关长度指数 $\nu$ 的恒等式，我们做两个假设：

**假设 1**: 相关长度 $\xi$ 是齐次函数：

$$
\xi(t, h) \sim |t|^{-\nu} g_\xi\left(\frac{h}{|t|^\Delta}\right)
$$

**假设 2**: 在临界点附近，相关长度 $\xi$ 是系统中**唯一重要的长度尺度**，独自负责热力学量的奇异贡献。

### 6.2 自由能的体积依赖

由于 $\ln Z(t,h)$ 是广延量且无量纲，它必须取形式：

$$
\ln Z = \left(\frac{L}{\xi}\right)^d \times g_s + \cdots + \left(\frac{L}{a}\right)^d \times g_a
$$

其中：
- $g_s$: 奇异函数（非解析）
- $g_a$: 解析函数（来自微观尺度 $a$）
- $L$: 系统尺寸

**物理解释**: 将系统分成大小为 $\xi$ 的独立单元，每个单元贡献常数因子到临界自由能。单元数目 $\sim (L/\xi)^d$。

**奇异自由能**: 

$$
f_{\text{sing}}(t, h) \sim \frac{\ln Z}{L^d} \sim \xi^{-d} \sim |t|^{\nu d} g_f\left(\frac{h}{|t|^\Delta}\right)
$$

### 6.3 Josephson 恒等式（超标度关系）

比较自由能的两种表达式：

$$
|t|^{2-\alpha} \sim |t|^{\nu d}
$$

因此：

$$
\boxed{2 - \alpha = \nu d} \quad \text{[Josephson's identity / hyperscaling]}
$$

**重要注释**: 
- 此关系涉及空间维度 $d$
- 称为**超标度关系**（hyperscaling relations）
- 与表中指数一致
- **失效**: 对于 $d > 4$，鞍点值 $\alpha = 0$, $\nu = 1/2$ 不满足此关系

**物理解释**: 对于 $d > 4$，涨落不重要，平均场理论成立，超标度失效。

---

## 7. 临界关联函数与自相似性

### 7.1 临界点的幂律衰减

恰好处于临界点时，相关长度无限大，没有长度尺度（除样品尺寸外）来截断关联函数的衰减。因此所有关联都按**幂律**衰减。

**磁化强度关联**:

$$
G_{c,mm}(\mathbf{x}) \equiv \langle m(\mathbf{x})m(0) \rangle - \langle m \rangle^2 \sim \frac{1}{|\mathbf{x}|^{d-2+\eta}}
$$

**能量-能量关联**:

$$
G_{c,EE}(\mathbf{x}) = \langle \epsilon(\mathbf{x})\epsilon(0) \rangle - \langle \epsilon \rangle^2 \sim \frac{1}{|\mathbf{x}|^{d-2+\eta'}}
$$

### 7.2 Fisher 恒等式

磁化率可从关联函数的积分得到：

$$
\chi \sim \int d^d x \, G_{c,mm}(\mathbf{x}) \sim \int d^d x \frac{1}{|\mathbf{x}|^{d-2+\eta}} \sim \xi^{2-\eta} \sim |t|^{-\nu(2-\eta)}
$$

与 $\chi \sim |t|^{-\gamma}$ 比较：

$$
\boxed{\gamma = \nu(2-\eta)} \quad \text{[Fisher's identity]}
$$

类似地，对于热容：

$$
C \sim \int d^d x \, G_{c,EE}(\mathbf{x}) \sim \xi^{2-\eta'} \sim |t|^{-\nu(2-\eta')}
$$

因此：

$$
\alpha = \nu(2-\eta')
$$

### 7.3 自相似性与分形

在临界点，关联函数具有**额外的膨胀对称性**（dilation symmetry）：

$$
G_{\text{critical}}(\mathbf{x}) = \lambda^p G_{\text{critical}}(\lambda \mathbf{x})
$$

这意味着**尺度不变性**或**自相似性**: 如果临界系统的快照被放大因子 $\lambda$，除了对比度改变（乘以 $\lambda^p$），结果快照在统计上与原始快照相似。

**分形几何**: 这种统计自相似性是**分形几何**的标志。如Mandelbrot讨论的，许多自然形式（云、海岸线、流域等）都表现出这种行为。

**Landau-Ginzburg 概率**: 基于旋转不变性等局域对称性构造。如果能加上膨胀对称性的要求，所得概率就能描述临界点。不幸的是，直接看这如何约束有效哈密顿量并不容易。例外是在 $d=2$，其中膨胀对称性意味着**共形不变性**，可以通过构造共形不变理论获得大量信息。

**替代路径**: 遵循重整化群程序中膨胀操作对有效能量的影响。

---

## 8. 重整化群：概念基础

### 8.1 Kadanoff 的洞察

标度理论成功预测各种指数恒等式，强烈支持这样的假设：在临界点附近，**相关长度 $\xi$ 是唯一重要的长度尺度**，微观长度尺度是不相关的。临界行为由高达尺度 $\xi$ 的**自相似涨落**主导。

**Kadanoff 建议**: 利用涨落的自相似性，**逐渐消除**长度尺度 $x \ll \xi$ 上的关联自由度，直到剩下尺度 $\xi$ 上相对简单、无关联的自由度。

### 8.2 重整化群的三个步骤

**步骤 1: 粗粒化（Coarse grain）**

系统有一个隐含的短距离截断 $a$（晶格间距或Landau-Ginzburg哈密顿量的粗粒化尺度）。RG的第一步是将此最小尺度改为 $ba$（$b > 1$）。粗粒化磁化强度：

$$
m'(\mathbf{x}) = \frac{1}{b^d} \int_{\text{以 }\mathbf{x}\text{ 为中心的尺寸为 }b\text{ 的单元}} d^d x' \, m(\mathbf{x}')
$$

**步骤 2: 重标度（Rescale）**

由于分辨率改变，粗粒化"图像"比原始图像更粗糙。通过将所有长度尺度减少因子 $b$ 来恢复原始分辨率：

$$
\mathbf{x}_{\text{new}} = \frac{\mathbf{x}_{\text{old}}}{b}
$$

**步骤 3: 重整化（Renormalize）**

涨落变化一般与原始图像不同，即两幅图像之间存在**对比度**差异。这可以通过引入因子 $\zeta$ 的对比度改变来补救，通过定义重整化磁化强度：

$$
m_{\text{new}}(\mathbf{x}_{\text{new}}) = \frac{1}{\zeta b^d} \int_{\text{以 }b\mathbf{x}_{\text{new}}\text{ 为中心的单元}} d^d x \, m(\mathbf{x})
$$

### 8.3 参数空间的流动

通过这些步骤，对每个构型 $m_{\text{old}}(\mathbf{x})$，我们生成重整化构型 $m_{\text{new}}(\mathbf{x})$。方程 (4.29) 可视为从一组随机变量到另一组随机变量的映射，可用于构造概率分布或权重 $\mathcal{W}_b[m_{\text{new}}(\mathbf{x})] \propto \exp(-\mathcal{H}_b[m_{\text{new}}(\mathbf{x})])$。

**关键假设**: 由于 $x \ll \xi$ 上重整化构型在统计上与原始构型相似，它们可能由与原始构型"接近"的哈密顿量 $\mathcal{H}_b$ 分布。特别是，原始哈密顿量通过调节两个参数 $t$ 和 $h$ 到零而变得临界：此时原始构型在统计上与重标度系统的构型相似。临界哈密顿量因此在重标度和重整化下应该**不变**。

**简化假设**: 原始和重整化哈密顿量接近临界性的程度可以用两个参数 $t$ 和 $h$ 描述。

**RG 变换对参数的影响**:

由于变换只涉及最短长度尺度的改变，它不会导致任何奇异性。重整化参数必须是原始参数的解析函数，因此可展开为：

$$
\begin{cases}
t' \equiv t_b = A_b t + B_b h + \cdots \\
h' \equiv h_b = C_b t + D_b h + \cdots
\end{cases}
$$

由于 $t = h = 0$ 时系统处于临界点，必须有 $t' = h' = 0$（无常数项）。

**旋转对称性**: 由于旋转对称性被 RG 保持，在联合变换 $m(\mathbf{x}) \to -m(\mathbf{x})$，$h \to -h$，$t \to t$ 下构型权重不变。因此系数 $B$ 和 $C$ 必须为零：

$$
\begin{cases}
t' = A_b t + \cdots \\
h' = D_b h + \cdots
\end{cases}
$$

**群性质**: 剩余的系数 $A_b$ 和 $D_b$ 依赖于（任意的）重标度因子 $b$，显然 $A_1 = D_1 = 1$（$b=1$ 时）。由于变换可以**顺序执行**，且 $b_1$ 和 $b_2$ 重标度的净效应是尺度改变 $b_1 b_2$，RG 程序有时被称为**半群**（semi-group）。

这导致：

$$
A_b = b^{y_t}, \quad D_b = b^{y_h}
$$

因此：

$$
\boxed{\begin{cases} t' = b^{y_t} t + \cdots \\ h' = b^{y_h} h + \cdots \end{cases}}
$$

其中 $y_t$ 和 $y_h$ 是**反常维度**（anomalous dimensions）。

**关键性质**: 如果 $\mathcal{H}_{\text{old}}$ 稍微偏离临界性，它由大但有限的相关长度 $\xi_{\text{old}}$ 描述。RG 变换后，由于步骤 2 中的重标度，新的相关长度小了因子 $b$。因此重整化哈密顿量**不那么临界**，RG 程序将参数进一步推离原点。因此指数 $y_t$ 和 $y_h$ 必须为正。

---

## 9. 从 RG 推导标度形式

### 9.1 自由能

RG 变换是从原始构型到新构型的**多对一映射**。由于新构型的权重 $\mathcal{W}'[m']$ 是旧构型权重 $\mathcal{W}[m]$ 的和，配分函数保持不变：

$$
Z = \int \mathcal{D}m \, \mathcal{W}[m] = \int \mathcal{D}m' \, \mathcal{W}'[m'] = Z'
$$

因此 $\ln Z = \ln Z'$，相应的自由能满足：

$$
V f(t, h) = V' f(t', h')
$$

在 $d$ 维中，重标度体积小了因子 $b^d$，因此：

$$
f(t, h) = b^{-d} f(b^{y_t} t, b^{y_h} h)
$$

**选择 $b$**: 使 $b^{y_t}$ 为常数（如 1），即 $b = t^{-1/y_t}$，给出：

$$
f(t, h) = t^{d/y_t} f(1, h/t^{y_h/y_t}) \equiv t^{d/y_t} g_f(h/t^{y_h/y_t})
$$

**识别指数**:

与标度形式 $f(t, h) = t^{2-\alpha} g_f(h/t^\Delta)$ 比较：

$$
\boxed{2 - \alpha = \frac{d}{y_t}, \quad \Delta = \frac{y_h}{y_t}}
$$

### 9.2 相关长度

所有长度尺度在 RG 变换中都减少因子 $b$。这也适用于相关长度：$\xi' = \xi/b$，意味着：

$$
\xi(t, h) = b \cdot \xi(b^{y_t} t, b^{y_h} h) = t^{-1/y_t} \xi(1, h/t^{y_h/y_t}) \sim |t|^{-\nu}
$$

这识别出：

$$
\boxed{\nu = \frac{1}{y_t}}
$$

使用 $\alpha$ 的表达式，**Josephson 恒等式** $2 - \alpha = \nu d$ 被**恢复**。

### 9.3 磁化强度

从 $Z$、$V$ 和 $h$ 的 RG 结果，我们可以直接得出：

$$
m(t, h) = -\frac{1}{V} \frac{\partial \ln Z}{\partial h} = -\frac{1}{b^d V'} \frac{\partial \ln Z'}{b^{-y_h} \partial h'} = b^{y_h - d} m(b^{y_t} t, b^{y_h} h)
$$

选择 $b = t^{-1/y_t}$：

$$
m(t, h) = t^{(y_h - d)/y_t} m(1, h/t^{y_h/y_t})
$$

识别出 $\beta = (y_h - d)/y_t$ 和 $\Delta = y_h/y_t$（与之前一致）。

### 9.4 一般量的齐次形式

显然，任何量 $X$ 的奇异部分具有齐次形式：

$$
X(t, h) = b^{y_X} X(b^{y_t} t, b^{y_h} h)
$$

对于任何**共轭变量对**（在哈密顿量中以 $\int d^d x F \cdot X$ 形式出现），标度维度满足：

$$
y_X = y_F - d
$$

其中 $F' = b^{y_F} F$ 在 RG 下。

---

## 10. 重整化群：形式理论

### 10.1 参数空间中的流动

**步骤 1**: 从最一般的哈密顿量开始（受对称性约束）。例如，在旋转对称性下：

$$
\mathcal{H} = \int d^d x \left[ \frac{t}{2} m^2 + u m^4 + v m^6 + \cdots + \frac{K}{2} (\nabla m)^2 + \frac{L}{2} (\nabla^2 m)^2 + \cdots - h m \right]
$$

具有这种对称性的特定系统因此完全由（无限维）参数空间 $S \equiv (t, u, v, \ldots, K, L, \ldots)$ 中的点指定。

**步骤 2**: 应用构型空间中的重整化三步：(i) $b$ 粗粒化；(ii) 重标度，$\mathbf{x}' = \mathbf{x}/b$；(iii) 重整化，$m' = m/\zeta$。

这完成变量变换：

$$
m'(\mathbf{x}') = \frac{1}{\zeta b^d} \int_{\text{尺寸为 }b\text{ 且以 }b\mathbf{x}'\text{ 为中心的单元}} d^d x \, m(\mathbf{x})
$$

给定原始构型的概率 $\mathcal{W}[m(\mathbf{x})] \propto \exp(-\mathcal{H}[m(\mathbf{x})])$，我们可以使用上述变量变换构造新构型的相应概率 $\mathcal{H}'[m'(\mathbf{x}')]$。自然这是程序中最困难的步骤。

**步骤 3**: 由于旋转对称性被 RG 程序保持，重标度哈密顿量也必须由方程 (4.41) 的参数空间中的点描述：

$$
\mathcal{H}'[m'(\mathbf{x}')] \equiv -\ln \mathcal{W}'[m'(\mathbf{x}')] = \int d^d x' \left[ \frac{t'}{2} m'^2 + u' m'^4 + v' m'^6 + \cdots + \frac{K'}{2} (\nabla' m')^2 + \frac{L'}{2} (\nabla'^2 m')^2 + \cdots \right]
$$

重整化参数是原始参数的函数：$t' = t_b(t, u, \ldots)$，$u' = u_b(t, u, \ldots)$ 等，定义参数空间中的**映射** $\mathbf{S}' = \mathbf{R}_b(\mathbf{S})$。

**步骤 4**: 操作 $\mathbf{R}_b$ 描述膨胀对系统哈密顿量的影响。描述统计自相似构型的哈密顿量必须对应**不动点** $\mathbf{S}^*$，使得 $\mathbf{R}_b(\mathbf{S}^*) = \mathbf{S}^*$。

由于相关长度 $\xi$（哈密顿量参数的函数）在 RG 操作下减少 $b$（即 $\xi' = \xi/b = \xi(\mathbf{R}_b(\mathbf{S}))$），不动点的相关长度必须为零或无限：

- $\xi^* = 0$: 描述每点独立涨落，对应完全无序（无限温度）或完全有序（零温度）
- $\xi^* = \infty$: 描述临界点（$T = T_c$）

### 10.2 线性化与稳定性分析

在不动点 $\mathbf{S}^*$ 附近线性化递推关系：

在 RG 下，点 $\mathbf{S}^* + \delta\mathbf{S}$ 变换为：

$$
\mathbf{S}^* + \delta\mathbf{S}' = \mathbf{S}^* + \mathbf{L}_b \cdot \delta\mathbf{S} + \cdots
$$

其中 $\mathbf{L}_b \equiv \frac{\partial \mathbf{R}_b}{\partial \mathbf{S}}|_{\mathbf{S}^*}$ 是线性化矩阵。

**对角化**: 我们本征化矩阵 $\mathbf{L}_b$ 得到**本征矢量** $\mathbf{e}_i$ 和相应**本征值** $\Lambda_i$：

$$
\mathbf{L}_b \cdot \mathbf{e}_i = \Lambda_i \mathbf{e}_i
$$

**群性质**: 由于 $\mathbf{L}_{b_1} \mathbf{L}_{b_2} = \mathbf{L}_{b_1 b_2}$，本征值满足：

$$
\Lambda_i(b) = b^{y_i}
$$

其中 $y_i$ 是**反常维度**。

**分类**:

- **$y_i > 0$**: $g_i$ 在标度下增加，$\mathbf{e}_i$ 是**相关算符**（relevant operator）
- **$y_i < 0$**: $g_i$ 在标度下减小，$\mathbf{e}_i$ 是**无关算符**（irrelevant operator）
- **$y_i = 0$**: $g_i$ 称为**边缘算符**（marginal operator），需要高阶项追踪其行为

### 10.3 普适性的解释

由无关算符张成的子空间是**不动点的吸引域**（basin of attraction）。由于 $\xi$ 总是在 RG 下减少，且 $\xi(\mathbf{S}^*) = \infty$，那么对于吸引域上的任何点，$\xi$ 也是无限的。

对于 $\mathbf{S}^*$ 附近的一般点，相关长度满足：

$$
\xi(g_1, g_2, \ldots) = b \cdot \xi(b^{y_1} g_1, b^{y_2} g_2, \ldots)
$$

对于足够大的 $b$，所有无关算符标度到零。$\xi$ 的**主导奇异性**然后由剩余的**相关算符**集合决定。

**普适性的 RG 解释**: 

想象不动点 $\mathbf{S}^*$ 描述零外场下磁体的临界点。当温度（或某个其他控制参数）改变时，有效哈密顿量的系数改变，点 $\mathbf{S}$ 沿参数空间中的轨迹移动。

除了临界点处的一个点外，磁体具有有限的相关长度。如果 $\mathbf{S}$ 取的轨迹仅在一点与 $\mathbf{S}^*$ 的吸引域相交，就能实现这一点。为此，吸引域必须有**余维数一**，即不动点 $\mathbf{S}^*$ 必须**有且仅有一个相关算符**。

这解释了**普适性**: 系统的非常多个微观细节组成由无关算符构成的巨大空间（吸引域）。

**磁场的情况**: 存在磁场时，必须调整两个系统参数才能达到临界点（$T = T_c$ 和 $h = 0$）。因此磁场对应于 $\mathbf{S}^*$ 处的**额外相关算符**。其他"奇"相互作用，如 $m^3, m^5, \ldots$，不应该导致任何其他相关算符。

---

## 11. 总结与展望

### 11.1 标度假设的核心内容

1. **齐次性**: 临界点附近自由能和相关量的齐次函数形式
2. **指数关系**: 从齐次性导出临界指数之间的关系（Rushbrooke、Widom、Fisher、Josephson 等）
3. **超标度**: 涉及维度 $d$ 的关系，与平均场理论在 $d > 4$ 时失效有关

### 11.2 重整化群的核心思想

1. **粗粒化**: 逐步消除短程自由度
2. **重标度**: 恢复原始分辨率
3. **重整化**: 调整对比度

**关键洞察**: 临界点是 RG 的不动点。相关算符对应物理参数（$t, h$），无关算符对应微观细节（普适性的来源）。

### 11.3 未解决的问题

1. **如何实现 RG 变换**: Kadanoff 的原始表述是概念性的，缺乏具体实施方案
2. **无限维参数空间**: 对称性允许无限多个相互作用，参数空间不便处理
3. **不动点的存在性**: 如何先验地知道 RG 变换存在不动点？
4. **线性化假设**: 为什么可以在不动点附近线性化？
5. **相关算符稀少**: 为什么相关算符很少？

### 11.4 Wilson 的贡献

Wilson 展示了如何在 **Landau-Ginzburg 模型**中（至少微扰地）实施这些步骤。这是下一章的内容。

---

## 关键公式总结

| 概念 | 公式 | 说明 |
|-----|------|------|
| **齐次性假设** | $f(t, h) = t^{2-\alpha} g_f(h/t^\Delta)$ | 自由能标度形式 |
| **Rushbrooke 恒等式** | $\alpha + 2\beta + \gamma = 2$ | 热力学关系 |
| **Widom 恒等式** | $\beta(\delta + 1) = 2 - \alpha + \gamma$ | 或 $\beta\delta = \beta + \gamma$ |
| **Fisher 恒等式** | $\gamma = \nu(2-\eta)$ | 关联函数与响应 |
| **Josephson 恒等式** | $2 - \alpha = \nu d$ | 超标度关系 |
| **RG 递推** | $t' = b^{y_t} t$, $h' = b^{y_h} h$ | 相关长度指数 $y_t, y_h > 0$ |
| **指数-反常维度关系** | $2-\alpha = d/y_t$, $\Delta = y_h/y_t$, $\nu = 1/y_t$ | RG 与标度联系 |
| **临界关联** | $\langle m(x)m(0) \rangle \sim |x|^{-(d-2+\eta)}$ | 幂律衰减 |

---

## 思考问题

1. **标度假设与平均场理论的关系**: 标度假设如何包含平均场理论作为特例？在什么条件下标度假设失效？

2. **超标度关系的物理**: 为什么 $2 - \alpha = \nu d$ 称为"超"标度？它在高维度为什么失效？

3. **RG 的群结构**: 为什么 RG 被称为"半群"而不是完整的群？信息在哪里丢失？

4. **相关 vs 无关算符**: 从物理上解释为什么 $t$ 和 $h$ 是相关算符，而 $u, v, K, L$ 等是无关或边缘算符？

5. **自相似与分形**: 临界系统的自相似性与分形几何的联系是什么？可以用什么量来刻画临界涨落的"分形维数"？

6. **从标度到 RG**: 标度假设是"现象学"的，而 RG 提供了"微观基础"。Wilson 的 RG 如何具体实现 Kadanoff 的三步骤？

---

## 参考文献

- Kadanoff, L. P. (1966). "Scaling laws for Ising models near $T_c$." *Physics Physics Physique Fizika*, 2(6), 263.
- Widom, B. (1965). "Equation of state in the neighborhood of the critical point." *The Journal of Chemical Physics*, 43(11), 3898-3905.
- Wilson, K. G. (1971). "Renormalization group and critical phenomena. I. Renormalization group and the Kadanoff scaling picture." *Physical Review B*, 4(9), 3174.
- Wilson, K. G., & Fisher, M. E. (1972). "Critical exponents in 3.99 dimensions." *Physical Review Letters*, 28(4), 240.
- Cardy, J. (1996). *Scaling and Renormalization in Statistical Physics*. Cambridge University Press.

---

*笔记完成于 2026-02-07*  
*统计场论的核心：标度、重整化群、普适性 🩷*
