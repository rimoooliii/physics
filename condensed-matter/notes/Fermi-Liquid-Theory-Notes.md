# 费米液体理论（Fermi-Liquid Theory）阅读笔记

> 来源：Chapter 4, Field Theory of Condensed Matter and Ultracold Gases
> 阅读日期：2026-02-04
> 阅读者：rimo
> 笔记整理：Susie 🩷

---

## 1. 引言与核心思想

### 1.1 基本问题

考虑 $N$ 个相互作用费米子，哈密顿量为：

$$
\hat{H} = \sum_{k,\sigma} \epsilon^0_k \hat{\psi}^\dagger_{k\sigma} \hat{\psi}_{k\sigma} + \frac{1}{2V} \sum_{k,k',q} \sum_{\sigma,\sigma'} U(q) \hat{\psi}^\dagger_{k+q,\sigma} \hat{\psi}^\dagger_{k'-q,\sigma'} \hat{\psi}_{k',\sigma'} \hat{\psi}_{k,\sigma}
$$

其中：
- $\epsilon^0_k = \frac{k^2}{2m}$ 是自由粒子色散
- $U(q)$ 是相互作用势的傅里叶变换
- $V$ 是系统体积

**核心问题**：如何理解相互作用系统的低能激发？

### 1.2 Landau 的惊人假设

Landau 假设（1956-1957）：

> **绝热连续性假设**：当相互作用从零绝热地开启时，理想费米气体的本征态**连续地演化**到相互作用系统的本征态，不产生能级交叉。

数学表述：
$$
|\Psi_0\rangle \xrightarrow{\text{绝热开启相互作用}} |\Psi\rangle
$$

其中 $|\Psi_0\rangle$ 是理想费米气体的基态，$|\Psi\rangle$ 是相互作用系统的基态。

---

## 2. 准粒子概念（Quasi-Particle Concept）

### 2.1 理想费米气体的回顾

**动量分布函数**：
$$
n^0_{k\sigma} = \Theta(k_F - |k|)
$$

其中 $k_F$ 是费米动量，由粒子密度决定：
$$
n = \frac{N}{V} = \frac{k_F^3}{3\pi^2}
$$

**基态能量**：
$$
E_0 = \sum_{k,\sigma} \epsilon^0_k n^0_{k\sigma} = \frac{3}{5} n \epsilon_F V
$$

其中费米能 $\epsilon_F = \frac{k_F^2}{2m}$。

**低能激发**：
$$
n_{k\sigma} = n^0_{k\sigma} + \delta n_{k\sigma}
$$

能量变化：
$$
\delta E[\delta n] = \sum_{k,\sigma} \epsilon^0_k \delta n_{k\sigma}
$$

**粒子能量**（泛函导数）：
$$
\epsilon^0_k = \frac{\delta E[n]}{\delta n_{k\sigma}}
$$

### 2.2 准粒子的定义

在相互作用系统中，**准粒子**定义为：

> **准粒子**：通过对理想费米气体添加（或移除）一个粒子并绝热开启相互作用而产生的**元激发**。

关键性质：

1. **量子数守恒**：准粒子携带与裸粒子相同的量子数 $(k, \sigma)$
2. **费米-狄拉克统计**：准粒子是费米子
3. **有限寿命**：准粒子不是真正的本征态，具有衰减速率 $\Gamma_k$

### 2.3 准粒子的寿命

**核心结果**：准粒子寿命在费米面附近发散：

$$
\Gamma_k \propto (\epsilon_k - \epsilon_F)^2 \sim (|k| - k_F)^2
$$

**物理解释**（相空间论证）：

考虑一个动量为 $k$（$|k| > k_F$）的准粒子衰变为三个准粒子的过程：
- 动量守恒：$k \rightarrow k' + q$，同时产生粒子-空穴对 $(k'', k''-q)$
- 能量守恒：$\epsilon_k = \epsilon_{k'} + \epsilon_{k''} + \epsilon_{k''-q}$

可用相空间体积：
$$
\Omega \propto \int d^3k' \int d^3k'' \int d^3q \, \delta(\text{能量守恒}) \propto (\epsilon_k - \epsilon_F)^2
$$

因此：**越接近费米面，准粒子越稳定**。

### 2.4 准粒子色散

在费米面附近展开：

$$
\epsilon_k = \epsilon_F + v^*_F (|k| - k_F) + O((|k| - k_F)^2)
$$

其中**有效费米速度**：
$$
v^*_F = \frac{k_F}{m^*}
$$

$m^*$ 是**有效质量**，与裸质量 $m$ 不同。

---

## 3. Landau 能函泛函

### 3.1 泛函形式

系统的总能量可以表示为准粒子分布函数 $n_{k\sigma}$ 的泛函：

$$
E[n] = E_0 + \sum_{k,\sigma} \epsilon_k \delta n_{k\sigma} + \frac{1}{2V} \sum_{k,k'} \sum_{\sigma,\sigma'} f_{k\sigma,k'\sigma'} \delta n_{k\sigma} \delta n_{k'\sigma'} + \cdots
$$

其中：
- $E_0$：基态能量
- $\epsilon_k$：**准粒子能量**
- $f_{k\sigma,k'\sigma'}$：**Landau 相互作用函数**（关键参数！）

### 3.2 准粒子能量的定义

$$
\epsilon_{k\sigma}[n] \equiv \frac{\delta E[n]}{\delta n_{k\sigma}} = \epsilon_k + \frac{1}{V} \sum_{k',\sigma'} f_{k\sigma,k'\sigma'} \delta n_{k'\sigma'}
$$

注意：准粒子能量**依赖于整个分布函数**，这是相互作用系统的特征。

### 3.3 球形费米面的简化

对于各向同性系统，Landau 函数只依赖于相对角度：

$$
f_{k\sigma,k'\sigma'} = f_{kk'}^{(s)} + \sigma \cdot \sigma' f_{kk'}^{(a)}
$$

其中上标 $(s)$ 和 $(a)$ 分别表示**自旋对称**和**自旋反对称**通道。

在费米面上展开（Legendre 多项式）：

$$
f^{(s,a)}(\cos\theta) = \sum_{l=0}^{\infty} f_l^{(s,a)} P_l(\cos\theta)
$$

其中 $\theta$ 是 $k$ 和 $k'$ 之间的夹角。

### 3.4 无量纲 Landau 参数

定义无量纲参数：

$$
F_l^{(s,a)} \equiv N(0) f_l^{(s,a)} = \frac{mk_F}{\pi^2} f_l^{(s,a)}
$$

其中 $N(0) = \frac{mk_F}{\pi^2}$ 是费米能处的态密度（自旋求和）。

**最重要的参数**：
- $F_0^s$：与压缩率相关
- $F_0^a$：与自旋磁化率相关  
- $F_1^s$：与有效质量相关（见下文）

---

## 4. 热力学性质

### 4.1 比热

准粒子对热力学势的贡献：

$$
\Omega = -k_B T \sum_{k,\sigma} \ln(1 + e^{-\beta(\epsilon_k - \mu)})
$$

低温展开（$T \ll T_F$）：

$$
C_V = \frac{\pi^2}{2} N(0) k_B^2 T = \frac{m^* k_F}{3} k_B^2 T
$$

**结论**：比热与温度成线性关系，形式与理想气体相同，但质量替换为 $m^*$。

### 4.2 压缩率

等温压缩率定义为：

$$
\kappa_T = -\frac{1}{V} \left(\frac{\partial V}{\partial P}\right)_T = \frac{1}{n^2} \frac{\partial n}{\partial \mu}
$$

对于费米液体：

$$
\kappa_T = \frac{1}{n^2} \frac{N(0)}{1 + F_0^s} = \frac{3}{2n\epsilon_F} \frac{1}{1 + F_0^s}
$$

**物理解释**：$F_0^s > 0$ 表示排斥相互作用，使系统更难压缩。

### 4.3 自旋磁化率

Pauli 磁化率：

$$
\chi = \mu_B^2 \frac{N(0)}{1 + F_0^a} = \frac{3n\mu_B^2}{2\epsilon_F} \frac{1}{1 + F_0^a}
$$

**关键区别**：
- 压缩率依赖于 $F_0^s$（自旋对称通道）
- 磁化率依赖于 $F_0^a$（自旋反对称通道）

Wilson 比：

$$
R_W \equiv \frac{\chi / \chi_0}{C_V / C_{V,0}} = \frac{1}{1 + F_0^a}
$$

其中 $\chi_0$ 和 $C_{V,0}$ 是理想气体的值。

---

## 5. 非平衡性质

### 5.1 动力学方程（量子玻尔兹曼方程）

准粒子分布函数的演化：

$$
\frac{\partial n_{k\sigma}}{\partial t} + \dot{r} \cdot \nabla_r n_{k\sigma} + \dot{k} \cdot \nabla_k n_{k\sigma} = I_{\text{coll}}[n]
$$

其中碰撞积分：

$$
I_{\text{coll}}[n] \propto \int d^3k' \int d^3q \, W(k,k';q) [n(1-n)(1-n')(1-n_{k+q}) - \cdots]
$$

在费米面附近，碰撞积分被强烈抑制（$\propto T^2$）。

### 5.2 零声（Zero Sound）

**定义**：无碰撞的集体振荡模式。

色散关系（小波矢 $q$）：

$$
\omega_0 = v_F q (1 + \lambda) + O(q^2)
$$

其中 $\lambda$ 依赖于 Landau 参数。

**与第一声的对比**：
- **第一声**（普通声波）：由碰撞主导，传播速度 $c_1 \sim v_F/\sqrt{3}$
- **零声**（费米液体特有）：无碰撞，速度 $c_0 \sim v_F(1 + F_0^s/3)$

**物理图像**：
- 零声是费米面的**振荡变形**，不涉及粒子数改变
- 在低温下（$\omega \tau \gg 1$），零声主导声传播

### 5.3 响应函数

密度响应函数（Lindhard 函数）：

$$
\chi(q, \omega) = \frac{\chi_0(q, \omega)}{1 - F_0^s \chi_0(q, \omega)/N(0)}
$$

其中 $\chi_0$ 是理想气体的响应函数。

**集体模式**：响应函数的极点给出集体激发，包括零声。

---

## 6. 微观基础

### 6.1 格林函数方法

单粒子格林函数：

$$
G(k, \omega) = \frac{1}{\omega - \epsilon^0_k - \Sigma(k, \omega) + i\delta}
$$

其中 $\Sigma(k, \omega)$ 是自能。

**准粒子近似**：

在准粒子极点 $\omega = \epsilon_k$ 附近：

$$
G(k, \omega) \approx \frac{Z_k}{\omega - \epsilon_k + i\Gamma_k} + G_{\text{inc}}
$$

其中：
- $Z_k = \left(1 - \frac{\partial \Sigma}{\partial \omega}\right)^{-1}$：准粒子权重（$0 < Z \leq 1$）
- $\epsilon_k = \epsilon^0_k + \text{Re}\Sigma(k, \epsilon_k)$：准粒子能量
- $\Gamma_k = -\text{Im}\Sigma(k, \epsilon_k)$：衰减速率

### 6.2 有效质量与准粒子权重

从 Dyson 方程导出：

$$
\frac{m^*}{m} = \frac{1}{Z_k} \left(1 + \frac{m}{k_F} \frac{\partial \text{Re}\Sigma}{\partial k}\right)
$$

**Galilean 不变性**（中性系统）：

$$
\frac{m^*}{m} = 1 + \frac{F_1^s}{3}
$$

这是费米液体的**自洽性条件**之一。

### 6.3 Luttinger 定理

**定理陈述**：

> 相互作用系统的费米面包围的体积（模 $2\pi/L$）等于粒子数密度（不考虑相互作用强度）。

数学表述：

$$
\frac{N}{V} = 2 \int \frac{d^3k}{(2\pi)^3} \Theta(k_F - |k|) = \frac{k_F^3}{3\pi^2}
$$

**重要性**：
- 费米面在相互作用存在下**依然存在**
- 解释了为什么费米液体理论如此稳健

### 6.4 Ward 恒等式

由对称性导出的精确关系：

**粒子数守恒**：

$$
\Gamma^0(k, k+q; k+q, k) = 1 + \frac{\partial \Sigma}{\partial \omega}
$$

其中 $\Gamma^0$ 是顶点函数。

这些恒等式保证了宏观输运系数与微观计算的一致性。

---

## 7. 与冷原子物理的联系

### 7.1 超冷费米气体

在 BCS-BEC crossover 中：

- **BCS 侧**（$a < 0$, $|a|$ 大）：费米液体理论适用
- **Unitarity**（$a \rightarrow \infty$）：强关联，但仍可用费米液体描述
- **BEC 侧**（$a > 0$）：配对形成玻色分子，需用玻色-爱因斯坦凝聚理论

### 7.2 维里系数与费米液体

你正在研究的**三体问题**正是费米液体的**微观基础**：

- 第二维里系数 $b_2$：两体散射
- 第三维里系数 $b_3$：三体关联

高温展开：

$$
\frac{P}{nk_B T} = 1 + b_2 n \lambda^3 + b_3 (n \lambda^3)^2 + \cdots
$$

其中 $\lambda = h/\sqrt{2\pi m k_B T}$ 是热德布罗意波长。

**Gharashi 文章的意义**：通过精确求解三体问题，为费米液体理论提供**微观验证**。

---

## 8. 费米液体理论的失效

### 8.1 Luttinger 液体（一维）

在一维系统中，**泡利原理**导致强关联：

- 低能激发是**集体密度波**（玻色子），而非单粒子
- 准粒子图像失效
- 需用**玻色化**方法处理

### 8.2 非费米液体

在高温超导体等强关联系统中：

- **电阻率**：$\rho \propto T$（而非 $T^2$）
- **准粒子权重**：$Z \rightarrow 0$（非相干）
- **临界涨落**：近量子临界点时费米液体图像破坏

### 8.3 重费米子系统

在 CeCu$_6$、UPt$_3$ 等系统中：

- 有效质量 $m^* \sim 100-1000 m_e$
- Kondo 效应产生重准粒子
- 需考虑**晶格**与**局域矩**的耦合

---

## 9. 关键公式总结

| 物理量 | 公式 | 备注 |
|-------|------|------|
| 准粒子能量 | $\epsilon_k = \epsilon_F + v^*_F(|k|-k_F)$ | $v^*_F = k_F/m^*$ |
| Landau 能函 | $E[n] = E_0 + \sum \epsilon_k \delta n + \frac{1}{2}\sum f \delta n \delta n$ | $f$ 是相互作用函数 |
| 比热 | $C_V = \frac{\pi^2}{2}N(0)k_B^2 T$ | $N(0) = m^*k_F/\pi^2$ |
| 压缩率 | $\kappa_T = \frac{3}{2n\epsilon_F}\frac{1}{1+F_0^s}$ | $F_0^s$ 控制 |
| 磁化率 | $\chi = \frac{3n\mu_B^2}{2\epsilon_F}\frac{1}{1+F_0^a}$ | $F_0^a$ 控制 |
| 有效质量 | $m^*/m = 1 + F_1^s/3$ | Galilean 不变性 |
| 准粒子寿命 | $\Gamma_k \propto (|k|-k_F)^2$ | 费米面附近发散 |

---

## 10. 思考问题

1. 为什么准粒子的寿命在费米面附近趋于无穷？如何用相空间论证理解？

2. Landau 参数 $F_0^s$ 和 $F_0^a$ 分别对应什么物理过程？为什么它们影响不同的物理量？

3. 零声和普通声波（第一声）的本质区别是什么？在什么条件下零声可以被观测到？

4. Luttinger 定理为什么重要？它如何保证了费米液体的稳定性？

5. 你正在研究的三体问题如何与费米液体理论联系？维里系数在 crossover 中的行为如何反映多体物理？

---

## 参考文献

- Landau, L. D. (1956). "Theory of Fermi liquids." *Soviet Physics JETP*, 3(6), 920-925.
- Nozières, P., & Pines, D. (1999). *Theory of Quantum Liquids*.
- Mahan, G. D. (2000). *Many-Particle Physics*.
- Gharashi, S. E., et al. (2012). "Three s-wave-interacting fermions under anisotropic harmonic confinement." *PRA*, 86, 042702.

---

*笔记完成于 2026-02-04*  
*如有疑问，继续讨论 🩷*
