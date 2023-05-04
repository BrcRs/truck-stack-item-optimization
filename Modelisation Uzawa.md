# ROADEF 2022 Modelisation

Minimize nb of trucks used and inventory of the plants.


## The full model

**Minimize:**

$$\alpha_T c^\top_{{T}_1} \zeta^T + \alpha_E c^\top_{{T}_2} \zeta^E + \alpha_Ic^\top_{I} (IDL - TI^\top \times TDA)$$

**Subject to:**


$$-\zeta^T \geq -\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]\quad \bold{(\zeta^T_1)}$$  

$$-\zeta^T \geq -TI^T\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]\quad \bold{(\zeta^T_2)}$$


$$-\zeta^E \geq -\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]\quad \bold{(\zeta^E_1)}$$  

$$-\zeta^E \geq -TI^E\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]\quad \bold{(\zeta^E_2)}$$


$$TI \leq TR\quad \bold{(TI_1)}$$

$$TI^\top\times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] \leq \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right]\quad \bold{(TI_2)}$$

$$TI^\top \times TP = IP\quad \bold{(TI_3)}\quad\text{(dominated by TR)}$$


$$R \leq \Theta M^{I4}\quad \bold{(R_1)}\quad\text{(dominated by TR)}$$  
$$ -M^{I4}(1-\Theta) \leq R - (TI^\top\times TU) \leq M^{I4}(1 - \Theta) \quad \bold{(R_2)}\quad\text{(dominated by TR)}$$  


$$R\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]=IU\quad \bold{(R_3)}\quad\text{(dominated by TR)}$$  

$$\Theta\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \geq \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(\Theta_1)}\quad\text{(dominated by TR)}$$

$$IDE \leq TI^\top \times TDA \leq IDL\quad \bold{(TI_4)}\quad\text{(dominated by TR)}$$

$$Z\leq S\times M^{Z}\quad \bold{(Z_1)}$$  


$$\Psi = \left [ \Omega[\dots, 1, \dots]\cdot \left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ] \Omega[\dots, 2, \dots]\cdot \left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ]\dots  \right ]\quad \bold{(\Psi_1)}$$

$$\Omega[\dots, j, \dots]\leq S\cdot M^\Omega\quad \forall j\quad \bold{(\Omega_1)}$$

$$-(1-S)M^\Omega \leq \Omega[i,\dots,\dots] - TI^\top \leq (1-S)M^\Omega\quad \forall i\quad \bold{(\Omega_2)}$$

$$\Psi \leq ST \cdot M^\Psi\quad \bold{(\Psi_2)}$$

$$-(1-ST)M^\Psi \leq \Psi - \left [S\cdot\left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ]\dots\right ] \leq (1-ST)M^\Psi\quad \bold{(\Psi_3)}$$



$$S\times IS = Z\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(S_1)}$$

$$-M^{Z}(1-S) \leq Z - \left [  SS  \cdots  SS  \right ] \leq M^{Z}(1-S)\quad \bold{(Z_2)}$$

$$Q\leq S\times M^Q\quad \bold{(Q_1)}$$  
$$S\times IU = Q\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(Q_2)}$$

$$-M^Q(1-S) \leq Q - \left [  SU  \cdots  SU  \right ] \leq M^Q(1-S)\quad \bold{(Q_3)}$$

$$H\leq S\times M^H\quad \bold{(H_1)}$$  
$$S\times IP = H\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(H_2)}$$

$$-M^H(1-S) \leq H - \left [  SP  \cdots  SP  \right ] \leq M^H(1-S)\quad \bold{(H_3)}$$

$$V\leq S\times M^V\quad \bold{(V_1)}$$  
$$S\times IK = V\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(V_2)}$$

$$-M^V(1-S) \leq V - \left [  SK  \cdots  SK  \right ] \leq M^V(1-S)\quad \bold{(V_3)}$$

$$W\leq S\times M^{W}\quad \bold{(W_1)}$$  
$$\displaystyle S\times IPD = W\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(W_2)}$$

$$-M^{W}(1-S) \leq W - \left [  SPD  \cdots  SPD  \right ] \leq M^{W}(1-S)\quad \bold{(W_3)}$$

$$G^l\leq S\times M^G\quad \bold{(G^l_1)}$$  
$$G^r \leq S\times M^G\quad \bold{(G^r_1)}$$  
$$G^l\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = G^r\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(G^lG^r)}$$

$$-M^G(1-S) \leq G^r - \left [  SO  \cdots  SO  \right ] \leq M^G(1-S)\quad \bold{(G^r_2)}$$  
$$-M^G(1-S) \leq G^l - \left [  IOV  \cdots  IOV  \right ] \leq M^G(1-S)\quad \bold{(G^l_2)}$$

$$SX^e - SX^o = SL + SO\cdot M^{TL}\quad \bold{(SX_1)}$$

$$SY^e - SY^o = SW + SO\cdot M^{TW}\quad \bold{(SY_1)}$$

$$SX^e - SX^o = SW + (1-SO)\cdot M^{TW}\quad \bold{(SX_2)}$$

$$SY^e - SY^o = SL + (1-SO)\cdot M^{TL}\quad \bold{(SY_2)}$$

$$SZ^o = 0\quad \bold{(\zeta^T_1)}$$

$$SZ^e = S\cdot IH\quad \bold{(SZ_1)}$$

$$SX^e \leq  ST \times TL\quad \bold{(SX_3)}$$  
$$SY^e \leq  ST \times TW\quad \bold{(SY_3)}$$  
$$SZ^e \leq  ST \times TH\quad \bold{(SZ_2)}$$

$$\left [ \begin{matrix} 1 & 0 & \dots & 0\\ & I & &  \end{matrix} \right ]\times SX^o \leq  SX^o\quad \bold{(SX_4)}$$

$$\Xi^2 SX^o - \Xi^1 SX^{e} - \beta^- + \beta^+ =  - 0.0001\quad \bold{(a)}\quad \bold{(SX_5)}$$

$$\beta^- \leq \lambda M^\lambda\quad \bold{(\beta^-_1)}$$  
$$\beta^+ \leq (1-\lambda)M^\lambda\quad \bold{(\beta^+_1)}$$  

$$(1 - \mu) \leq \beta^- M^\mu\quad \bold{(\mu_1)}$$

$$\left ( \Xi^2 ST - \Xi^1 ST \right )\times \left [ \begin{matrix} 1 & 0 & 0 & \dots \\ 0 & 2 & 0& \dots \\ 0 & 0 & 3& \dots\\ etc \end{matrix} \right ] = \nu \quad \quad \bold{(b)}\quad \bold{(\nu_1)}$$


$$\tau - \phi \leq (\nu-|trucks|)M^\tau\quad \bold{(\tau_1)}$$  
$$\displaystyle \tau \geq \frac{(\nu-|trucks|)M^\tau}{10}\quad \bold{(\tau_2)}$$  
<!-- $$\tau - \phi \leq \eta^1 M$$  
$$\phi - \tau \leq (1-\eta^1)M$$   -->
$$\tau \leq \eta M^\eta\quad \bold{(\tau_3)}$$  
$$\phi \leq (1-\eta)M^\eta\quad \bold{(\phi_1)}$$  

$$\Xi^1SY^e \leq \Xi^2SY^o +  \xi M^{TW} + (\tau + \phi)M^{TW} + (1-\mu)M^{TW}\quad \bold{(SY_4)}$$

$$\Xi^2SY^e \leq \Xi^1SY^o + (1 - \xi)M^{TW} +  (\tau + \phi)M^{TW} + (1-\mu)M^{TW}\quad \bold{(SY_5)}$$

$$\Xi^1SU\cdot TE \leq \Xi^2 SU \cdot TE + (\tau + \phi)M^{TE}\quad \bold{(SU_1)}$$

$$\Xi^1SU - \Xi^2SU \geq \chi\epsilon - rM^{TE} - (\tau + \phi)M^{TE} - (1 - \sigma^1)M^{TE}\quad \bold{(SU_2)}$$  
$$\Xi^2SU - \Xi^1SU \geq (1 - \chi)\epsilon - rM^{TE} - (\tau + \phi)M^{TE} - (1 - \sigma^1)M^{TE}\quad \bold{(SU_2)}$$

$$\Xi^2SK \cdot TKE \geq \Xi^1SK \cdot TKE - (1-r)M^{TKE} - (\tau + \phi)M^{TKE}\quad \bold{(SK_1)}$$

$$\Xi^1SK \cdot TKE - \Xi^2SK \cdot TKE \geq \chi\epsilon - (\tau + \phi)M^{TKE} - (1 - \sigma^2)M^{TKE}\quad \bold{(SK_2)}$$  
$$\Xi^2SK \cdot TKE - \Xi^1SK \cdot TKE \geq (1 - \chi)\epsilon - (\tau + \phi)M^{TKE} - (1 - \sigma^2)M^{TKE}\quad \bold{(SK_3)}$$

$$\Xi^2SG\cdot TGE \geq \Xi^1SG\cdot TGE - (\tau + \phi)M^{TGE} - (1 - \sigma^3)M^{TGE}\quad \bold{(SG_1)}$$

$$\sigma^1 + \sigma^2 + \sigma^3 \geq 1\quad \bold{(\sigma_1)}$$  


**Variables:**

$\zeta^T \geq 0$

$\zeta^T \in \bold{R^+}^{|PlannedTrucks|}$

$\zeta^E \geq 0$

$\zeta^E \in \bold{R^+}^{|ExtraTrucks|}$

$SS \in \bold{R^+}^{|stacks|}$

$SP \in \bold{R^+}^{|stacks|}$

$SK \in \bold{R^+}^{|stacks|}$

$SPD \in \bold{R^+}^{|stacks|}$

$SU \in \bold{R^+}^{|stacks|}$

$SO\in \bold{R^+}^{|stacks|}$

$SX^e \in \bold{R^+}^{|stacks|}$

$SX^o \in \bold{R^+}^{|stacks|}$

$SY^e \in \bold{R^+}^{|stacks|}$

$SY^o \in \bold{R^+}^{|stacks|}$

$SZ^e \in \bold{R^+}^{|stacks|}$

$\beta^- \in \bold{R^+}^{|stacks|}$

$\beta^+ \in \bold{R^+}^{|stacks|}$

$\nu \in \bold{R^+}^{|stacks|}$

$\tau\in \bold{R^+}^{|stacks|}$

$\phi\in \bold{R^+}^{|stacks|}$

$SG \in \bold{R^+}^{|stacks| \times plants}$

$TI :$ truck $-$ item matrix

$TI \in \{0,1\}^{|trucks| \times N}$

$R \in \{0,1\}^{N \times |suppliers|}$

$\Theta \in \{0,1\}^{N \times |suppliers|}$

$S \in \{0,1\}^{|stacks| \times N}$

$Z \in \{0,1\}^{|stacks| \times N}$

$\Omega \in \{0,1\}^{|stacks| \times |trucks| \times |items|}$

$ST\in \{0,1\}^{|stacks| \times |trucks|}$

$IOV\in \{0, 1\}^{|items|}$

$\mu \in \{0, 1\}^{|stacks|}$

$\eta \in \{0, 1\}^{|stacks|}$

$\xi \in \{0, 1\}^{|stacks|}$

$\chi\in \{0, 1\}^{|stacks|}$

$r\in \{0, 1\}^{|stacks|}$

$\sigma^1 \in \{0, 1\}^{|stacks|}$

$\sigma^2 \in \{0, 1\}^{|stacks|}$

$\sigma^3 \in \{0, 1\}^{|stacks|}$

$\lambda \in  \{0, 1\}^{|stacks|\times (|stacks|+1)/2}$

$\Psi \in \bold{N}^{|stacks| \times |trucks|}$

$Q \in \bold{N}^{|stacks| \times |items|}$

$H \in \bold{N}^{|stacks| \times |items|}$

$V \in \bold{N}^{|stacks| \times |items|}$

$W \in \bold{N}^{|stacks| \times |items|}$

$G^l \in \bold{N}^{|stacks| \times |items|}$

$G^r \in \bold{N}^{|stacks| \times |items|}$



**Parameters:**

$\displaystyle M^{I4}_{i,j} = \max_i TU_{i,j}\quad \forall j \in  \left \{1,\; Col(TI^\top)\right \}$

$M^Z = \max IS + 1$

$IS \in \bold{N}^{|items|}$

$M^Q = \max IU + 1$

$M^H = \max IP + 1$

$M^V = \max IK + 1$

$M^W = \max{IPD} + 1$

$G^r = \max SO + 1$

$G^l = \max IO + 1$


$\displaystyle M^\eta = |trucks|\frac{M^\tau}{10}$

$M^{TE} = \max TE$

$M^{TKE} = \max TKE$

$M^{TGE} = \max TGE$

$M^{TW} = \max TW + 1$

$M^\Psi = |items|$

$IU \in \bold{N}^{|items|}$

$IP \in \bold{N}^{|items|}$


$IK \in \bold{N}^{|items|}$

$IPD  \in \bold{N}^{|items|}$


$TL  \in \bold{R}^{|trucks|}$

$TW\in \bold{R}^{|trucks|}$

$TH\in \bold{R}^{|trucks|}$

$TE \in \bold{N}^{|trucks| \times |suppliers|}$

$TKE \in \bold{N}^{|trucks| \times |supplierDocks|}$

$TGE \in \bold{N}^{|trucks| \times |plantDocks|}$

$M^\lambda = 2TL + 1$

**Constants:**

$SZ^o = 0$

$M^{ST} = 2$

$M^\Omega = 2$

$M^\mu = 2$

$M^\tau = 10$  

$\Xi^1 : \left [ \begin{matrix}1 & 0 & 0 & 0&\dots\\1 & 0& 0& 0&\dots\\1 & 0 & 0 & 0&\dots\\ \vdots\\ 0&1&0&0&\dots\\0&1&0&0&\dots\\etc \end{matrix} \right ]$

$\Xi^2 : \left [ \begin{matrix}0 & 1 & 0 & 0&\dots\\0 & 0& 1& 0&\dots\\0 & 0 & 0 & 1&\dots\\ \vdots\\ 0&0&1&0&\dots\\0&0&0&1&\dots\\etc \end{matrix} \right ]$

$\epsilon = 0.001$


## Resolution

Solving the problem using a Benders decomposition has proven to be difficult, because of the large memory footprint it implied.

We choose now to see the problem as a stochastic 2-stage problem. The different placement subproblems associated with each trucks acn be viewed as scenarios depending on the first stage decision of placing which items in which truck. This kind of problem can be solved thanks to an algorithm "à la Uzawa". Each subproblem is given a separate $TI$ variable. Iteratively, the subproblems are penalized proportionaly to the difference of their particular $TI$ in regards to the mean $TI$. The algorithm eventually converges and all the $TI$ first stage variables end up being equal and optimal.

The problem must be able to be rewritten in the following form:

$$\max_{{\{\kappa^s\}}_{s\in\bold{S}}} \sum_{s\in\bold{S}} \pi^s \min_{(x^s, y^s)}\left ( j^s(x^s,y^s)+(\kappa^s-\sum_{s'\in\bold{S}}\pi^{s'}\kappa^{s'})x^s\right)$$

$$\begin{align*} x^s &\in \bold{X}\\
y^s &\in \cal{Y}^s(x^s)
\end{align*}$$

where $x^s$ is first stage decision for scenario $s$, $y^s$ is second stage decision, $\kappa^s$ is penalization of scenario $s$. $\pi$ is just a cost coefficient.

Scheme of the scenario decomposition algorithm:

> P. Carpentier, J.-P. Chancelier, G. Cohen, M. De Lara, Stochastic Multi-Stage Optimization, Springer Cham, 2015.

Input: step $\delta > 0$, initial multipliers $\{\kappa_{(0)}\}_{s\in\bold{S}}$ and first decision $\overline{x}_{(0)}$.
Output: optimal first decision $x$

```Julia
while x[s, k+1] - sum([pi[s] * x[s, k+1] for s in bold_S]) != 0
    for s in bold_S
        # Solve the deterministic minimization problem for scenario s with
        # a penalization + kappa[s, k] * (x[s, k+1] - xbar[k])
        # and obtain optimal first decision x[s, k+1]
        ...
    end
    # Update the mean first decisions:
    xbar[k+1] = sum([pi[s] * x[s, k+1] for s in bold_S])

    # Update the multipliers by
    for s in bold_S
        kappa[s, k+1] = kappa[s, k] + delta * (x[s, k+1] - xbar[k+1])
    end
```

In our case $\bold{S}$ corresponds to our set of planned trucks + extra trucks.

## Solving the placement subproblem

Given a set of items and a truck, place each item in a stack, so that:
1. The stacks don't overlap
2. The stacks stay in the truck
3. The items of a stack have the same stackability code, supplier, supplier dock, plant dock, orientation
4. The loading orders are respected

The following program places stacks in the truck optimally and taking into account the Uzawa penalization:

**Minimize:**

$$\alpha_T c^\top_{{T}_1} \zeta^T + \alpha_E c^\top_{{T}_2} \zeta^E + \alpha_Ic^\top_{I} (IDL - TI^{t\top} \times TDA) + \kappa^t(TI^t - \overline{TI})$$

**Subject to:**


$$-\zeta^{tT} \geq -1\quad \bold{(\zeta^T_1)}$$  

$$-\zeta^{tT} \geq -TI^{tT}\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]\quad \bold{(\zeta^T_2)}$$


$$-\zeta^{tE} \geq -1\quad \bold{(\zeta^E_1)}$$  

$$-\zeta^{tE} \geq -TI^{tE}\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]\quad \bold{(\zeta^E_2)}$$

$$TI^t \leq TR\quad \bold{(TI_1)}$$

$$TI^{t\top}\times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] = \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right]\quad \bold{(TI_2)}$$

The following one ensures that if an item is in a stack, the item is in the truck, because all stacks of the subproblem are in the truck:
$$S^\top \cdot \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] = TI^t[t]$$

$$Z\leq S\times M^{Z}\quad \bold{(Z_1)}$$  

<!-- And this one ensures items of the truck are in a stack:

$$-\left (\left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] - TI[t]\right )M^S\leq S^\top \times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] - \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] \leq \left (\left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] - TI[t]\right )M^S$$  -->

$$S\times IS = Z\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(S_1)}$$


<!-- Omega variables become useless when we know for sure that the stacks we are dealing with are the ones of the truck -->
<!-- $$\Omega[\dots, t, \dots]\leq S\cdot M^\Omega\quad \bold{(\Omega_1)}$$ -->

<!-- $$-(1-S)M^\Omega \leq \Omega[i,t,\dots] - TI^t[t]^\top \leq (1-S)M^\Omega\quad \forall i\quad \bold{(\Omega_2)}$$ -->

<!-- $$\Omega[\dots, t, \dots] \leq ST \cdot M^{\Omega^t}\quad \bold{(\Psi_2)}$$ -->

<!-- $$-(1-ST)M^{\Omega^t} \leq \Omega[\dots, t, \dots] - \left [S\cdot\left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ]\dots\right ] \leq (1-ST)M^{\Omega^t}\quad \bold{(\Psi_3)}$$ -->

$$-M^{Z}(1-S) \leq Z - \left [  SS  \cdots  SS  \right ] \leq M^{Z}(1-S)\quad \bold{(Z_2)}$$

<!-- $$Q\leq S\times M^Q\quad \bold{(Q_1)}$$  
$$S\times IU = Q\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(Q_2)}$$ -->

$$Q\leq SU\times M^Q\quad \bold{(Q_1)}$$  
<!-- $$S\times IU = Q\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(Q_2)}$$ -->
<!-- $$S\times IU = SU^\top \cdot\left (S\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\right ) \quad \bold{(Q_2)}$$ -->
$$S\times IU = SU \circ \left [S\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \dots S\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\right ] \quad \bold{(Q_2)}$$
$$S\times IU = Q \quad \bold{(Q_2)}$$


$$-M^Q(1-SU^\top) \leq Q - \left [S\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \dots S\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\right ] \leq M^Q(1-SU^\top)\quad \bold{(Q_3)}$$ 

<!-- The following are unnecessary because the stacks naturally have the same plant as the truck, and since items of the stack must be in the truck, the item also have the same plant -->
<!-- $$H\leq S\times M^H\quad \bold{(H_1)}$$  
$$S\times IP = H\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(H_2)}$$

$$-M^H(1-S) \leq H - \left [  SP  \cdots  SP  \right ] \leq M^H(1-S)\quad \bold{(H_3)}$$ -->

$$V\leq S\times M^V\quad \bold{(V_1)}$$  
$$S\times IK = V\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(V_2)}$$

$$-M^V(1-S) \leq V - \left [  SK  \cdots  SK  \right ] \leq M^V(1-S)\quad \bold{(V_3)}$$

$$W\leq S\times M^{W}\quad \bold{(W_1)}$$  
$$\displaystyle S\times IPD = W\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(W_2)}$$

$$-M^{W}(1-S) \leq W - \left [  SPD  \cdots  SPD  \right ] \leq M^{W}(1-S)\quad \bold{(W_3)}$$


$$G^l\leq S\times M^G$$  
$$G^r \leq S\times M^G$$  
$$G^l\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = G^r\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$-M^G(1-S) \leq G^r - \left [  SO  \cdots  SO  \right ] \leq M^G(1-S)$$  
$$-M^G(1-S) \leq G^l - \left [  IOV  \cdots  IOV  \right ] \leq M^G(1-S)$$

Define $SL$ and $SW$

<!-- $$SL\times S\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = S\times IL$$
$$SW\times S\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = S\times IW$$ -->

$$D^L \leq S \times M^{D^L}$$

$$-M^{D^L}(1 - S) \leq D^L - SL \leq M^{D^L}(1 - S)$$

$$D^L\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = S\times IL$$

$$D^W \leq S \times M^{D^W}$$

$$-M^{D^W}(1 - S) \leq D^W - SW \leq M^{D^W}(1 - S)$$

$$D^W\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = S\times IW$$


$$SX^e - SX^o = SL + SO\cdot M^{TL}$$

$$SY^e - SY^o = SW + SO\cdot M^{TW}$$

$$SX^e - SX^o = SW + (1-SO)\cdot M^{TW}$$

$$SY^e - SY^o = SL + (1-SO)\cdot M^{TL}$$

$$SZ^o = 0$$

$$SZ^e = S\cdot IH$$

$$SX^e \leq  TL[t]$$  

<!-- $$\rho \geq SX^e - TL$$ -->

$$SY^e \leq  TW[t]$$  
$$SZ^e \leq  TH[t]$$

$$\left [ \begin{matrix} 1 & 0 & \dots & 0\\ & I & &  \end{matrix} \right ]\times SX^o \leq  SX^o$$

$$\Xi^2 SX^o - \Xi^1 SX^{e} - \beta^- + \beta^+ =  - 0.0001\quad \bold{(\Xi_a)}$$

$$\beta^- \leq \lambda M^\lambda$$  
$$\beta^+ \leq (1-\lambda)M^\lambda$$  

$$(1 - \mu) \leq \beta^- M^\mu$$

$$\Xi^1SY^e \leq \Xi^2SY^o +  \xi M^{TW} + (1-\mu)M^{TW}\quad \bold{(\Xi_b)}$$

$$\Xi^2SY^e \leq \Xi^1SY^o + (1 - \xi)M^{TW} + (1-\mu)M^{TW}\quad \bold{(\Xi_c)}$$

$$\Xi^1SU\cdot TE \leq \Xi^2 SU \cdot TE\quad \bold{(\Xi_d)}$$

$$\Xi^1SU - \Xi^2SU \geq \chi\epsilon - rM^{TE} - (1 - \sigma^1)M^{TE}\quad \bold{(\Xi_e)}$$  
$$\Xi^2SU - \Xi^1SU \geq (1 - \chi)\epsilon - rM^{TE} - (1 - \sigma^1)M^{TE}\quad \bold{(\Xi_f)}$$

$$\Xi^2SK \cdot TKE \geq \Xi^1SK \cdot TKE - (1-r)M^{TKE}\quad \bold{(\Xi_g)}$$

$$\Xi^1SK \cdot TKE - \Xi^2SK \cdot TKE \geq \chi\epsilon - (1 - \sigma^2)M^{TKE}\quad \bold{(\Xi_h)}$$  
$$\Xi^2SK \cdot TKE - \Xi^1SK \cdot TKE \geq (1 - \chi)\epsilon - (1 - \sigma^2)M^{TKE}\quad \bold{(\Xi_i)}$$

$$\Xi^2SG\cdot TGE \geq \Xi^1SG\cdot TGE - (1 - \sigma^3)M^{TGE}\quad \bold{(\Xi_j)}$$

$$\sigma^1 + \sigma^2 + \sigma^3 \geq 1$$  


**Variables:**

$\rho \in \bold{R^+}^{|stacks|}$

$SU \in \bold{R^+}^{|stacks|} \quad \text{(Entirely infered)}$

$SO\in \bold{R^+}^{|stacks|}$

$SX^e \in \bold{R^+}^{|stacks|}$

$SX^o \in \bold{R^+}^{|stacks|}$

$SY^e \in \bold{R^+}^{|stacks|}$

$SY^o \in \bold{R^+}^{|stacks|}$

$SZ^e \in \bold{R^+}^{|stacks|} \quad \text{(Entirely infered from \textit{S})}$

$\beta^- \in \bold{R^+}^{|stacks|} \quad \text{(Entirely infered)}$

$\beta^+ \in \bold{R^+}^{|stacks|} \quad \text{(Entirely infered)}$

$SG \in \bold{R^+}^{|stacks| \times plants}$

$G^l \in \bold{R^+}^{|stacks| \times |items|} \quad \text{(Entirely infered from \textit{S})}$

$G^r \in \bold{R^+}^{|stacks| \times |items|} \quad \text{(Entirely infered from \textit{S})}$

$IOV\in \{0, 1\}^{|items|} \quad \text{(Partially determined)}$

$\mu \in \{0, 1\}^{|stacks|} \quad \text{(Partially infered from } \beta^-)$

$\xi \in \{0, 1\}^{|stacks|}$

$\chi\in \{0, 1\}^{|stacks|}$

$r\in \{0, 1\}^{|stacks|}$

$\sigma^1 \in \{0, 1\}^{|stacks|}$

$\sigma^2 \in \{0, 1\}^{|stacks|}$

$\sigma^3 \in \{0, 1\}^{|stacks|}$

$\lambda \in  \{0, 1\}^{|stacks|\times (|stacks|+1)/2}$



**Parameters:**

$S \in \{0,1\}^{|stacks| \times N}$

$M^{G^r} = \max SO + 1$

$M^{G^l} = \max IOV + 1$


$M^{TE} = \max TE$

$M^{TKE} = \max TKE$

$M^{TGE} = \max TGE$

$M^{TW} = \max TW + 1$

$TL  \in \bold{R}^{|trucks|}$

$TW\in \bold{R}^{|trucks|}$

$TH\in \bold{R}^{|trucks|}$

$TE \in \bold{N}^{|trucks| \times |suppliers|}$

$TKE \in \bold{N}^{|trucks| \times |supplierDocks|}$

$TGE \in \bold{N}^{|trucks| \times |plantDocks|}$

$M^\lambda = 2TL + 1$


**Constants:**

$M^S = 2$

$SZ^o = 0$

$M^\mu = 2$

$\Xi^1 : \left [ \begin{matrix}1 & 0 & 0 & 0&\dots\\1 & 0& 0& 0&\dots\\1 & 0 & 0 & 0&\dots\\ \vdots\\ 0&1&0&0&\dots\\0&1&0&0&\dots\\etc \end{matrix} \right ]$

$\Xi^2 : \left [ \begin{matrix}0 & 1 & 0 & 0&\dots\\0 & 0& 1& 0&\dots\\0 & 0 & 0 & 1&\dots\\ \vdots\\ 0&0&1&0&\dots\\0&0&0&1&\dots\\etc \end{matrix} \right ]$

$\epsilon = 0.001$

This is a Mixed Integer Linear Program. We could try to solve it with B&B, with a smart branching on the binary variables which are not determined by other variables. The use of heuristic algorithms which ignore orientation for instance could prove useful. We could even use Benders decomposition.

## In Summary

```
1. Make first solution by distributing items in planned trucks
2. Instanciate as many subproblems than there are trucks. In each subproblem, only create constraints related to the corresponding truck.
while TI[t, k+1] - sum([pi[t] * TI[t, k+1] for t in bold_T]) != 0
    for t in bold_T
        # Solve the deterministic minimization problem for truck t with
        # a penalization + kappa[t, k] * (TI[t, k+1] - TIbar[k])
        # and obtain optimal first decision TI[t, k+1]
        ...
    end
    # Update the mean first decisions:
    TIbar[k+1] = sum([pi[t] * TI[t, k+1] for t in bold_T])

    # Update the multipliers by
    for t in bold_T
        kappa[t, k+1] = kappa[t, k] + delta * (TI[t, k+1] - TIbar[k+1])
    end
```