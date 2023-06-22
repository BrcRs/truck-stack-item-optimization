# Comparative Analysis of Different Solving Algorithms for the Truck-Stack-Item Optimization problem

This document aims to provide insights on the time and space complexity of different solving algorithms for the Truck-Stack-Item problem described in the ROADEF Challenge 2022.

## Truck Based Column Generation

The truck based column generation algorithm is an iterative algorithm initially solving the problem with a minimal number of trucks, and iteratively adding trucks until adding new trucks does not improve the solution anymore.

### Description

```
function columngeneration(solvefun!, problem, args...; eps=0.1)

    Let be m the maximal number of trucks

    Let be n the number of items

    Consider a list of trucks useless for the problem:
    deadtrucks is an empty list

    Consider a number of trucks chosen for the problem:
    chosentrucks is a minimal set of trucks

    Solve the problem with the initial set of trucks

    while adding trucks has proven useful

        Update chosen trucks based on a heuristic:
        E.g.: Add one new truck for each truck that have at least 1 item

        Solve the problem with chosen trucks

        Optional: remove trucks that do not participate in solution
        E.g.: remove empty unused useless trucks (update deadtrucks)

    end

    return solution
end
```

### Time Complexity

In the worst case, the algorithm adds one truck to the problem at each iteration. We would have thus $m$ iterations, $m$ being the maximum number of trucks.
The maximum number of trucks can be expressed in function of the number of items. Indeed, in the worst case we have one truck per item.

The number of iterations is $n$. But the size of the instance solved will increase with each iteration:

$$\sum_{k = 1}^{n}\phi(n, k)$$

This algorithm has thus a time complexity of $O(\sum_{k = 1}^{n}\phi(n, k))$, $\phi$ being the time complexity of the solving algorithm used internaly.

However in practice, we can presume that the actual number of trucks used $\hat{m}$ is significantly smaller than $n$. We could say that $p\times\hat{m} \leq n$.

We then have a time complexity of $O(\sum_{k = 1}^{n/p}\phi(n, k))$.

### Space Complexity

The space complexity of this algorithm will greatly depend on the solver used internaly. On its own its space complexity can be approximated to be $O(n)$. The utility of this method is to reduce the memory footprint of the resolution by a factor of $n$.

### Pros and Cons

This algorithm allows to solve the problem with a minimal number of trucks, at the cost of solving it multiple times.

A better approximation of the number of trucks necessary to the resolution could avoid solving many smaller problems in which the number of trucks is not sufficient.

Calculating the cumulative volume of candidate items for a given planned truck could give an approximation of the number of corresponding necessary extra trucks. This operation would be done in $O(n)$.

## Uzawa Inspired Algorithm

The problem can be subdivided in multiple sub-problems consisting in placing the items optimally within each truck. The Uzawa inspired algorithm makes use of this subdivision by solving separately the sub-problems. However solving separately the sub-problems means we get different values for coupling variables *i.e.* which items go in which truck. To avoid this effect, the algorithm adds a penalization to the objective functions of each sub-problems so that they converge towards a unique global value for the coupling variables.

### Description

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

### Time Complexity

This algorithm solves each subproblem an unkown finite number of times $r$.

$$O(r\cdot\sum_{t=1}^m\psi_t(n))$$

Where $\psi_t$ is the complexity of solving subproblem $t$.

Combining this algorithm with the truck based column generation, we get a complexity of:

$$O(\sum_{k = 1}^{n/p}\left [r\cdot\sum_{t=1}^k\psi_t(n)\right ])$$

$$= O(r\cdot \frac{n}{p}\cdot\sum_{k = 1}^{n/p}\left [\sum_{t=1}^k\psi_t(n)\right ])$$
$$\approx O(r\cdot \frac{n}{p}\cdot\sum_{k = 1}^{n/p}k \cdot\psi(n))$$

$$\approx O(r\cdot \frac{n}{p}\cdot\psi(n)\cdot\sum_{k = 1}^{n/p}k )$$

$$\approx O(r\cdot\psi(n)\cdot(n^3/p^3 + n^2/p^2)\cdot\frac{1}{2} )$$
$$\approx O(r\cdot\psi(n)\cdot\frac{n^3}{p^3})$$

### Space Complexity

$O(m\cdot n^2)$ variables.

$O(n^2)$ constraints.

### Pros and Cons

The advantage of this method is the smaller memory footprint, and the time gained by solving smaller and "independant" problems.

However the convergence time is unknown and we get a solution only at the end of the iterations.

## Benders Decomposition

A Benders Decomposition can leverage the subdivision of the problem in smaller independant placing problems.

### Description

```
1. Preprocess the instance and formulate a linear relaxation of the original MILP.
2. Setup the master problem with the following objective function:
    min     z
    subject to
    cuts = {z >= 0}
3. Set TI, the variables assigning each item to a truck, to random values, or the optimal value if there weren't placement constraints.
4. Make sure (How?) that the value chosen for TI is feasible for all truck-wise sub-problems.
5. Solve separately (How?) each subproblem consisting in placing the items of a truck into stacks and placing those stacks into the truck.
6. Get the solution for X\{TI} (all variables except TI).
7. Compute the objective function of the dual of the linear relaxation of the original MILP. Add penalities for items with no stacks.
8. Generate an optimality cut and add it to `cuts`.
9. Solve master problem and obtain new TI.
10. Unless the optimal solution is found, return to 4.
```

### Time Complexity

Solving all subproblems would be done in $O(\sum_{t=1}^m\pi_t(n_t))$ with $\sum_1^m n_t = n$ and $\pi$ the function solving a subproblem. Formulating the dual objective function is easily done beforehand, and retrieving the dual values is done in $O(n^2)$?**TODO** Solving the master problem will be done with a standard simplex implementation which in practice is linear.

We thus have an iterative algorithm of time complexity $O(s\cdot\sum_{t=1}^m\pi_t(n_t))$ with $s$ the number of iterations necessary to converge to the optimal solution.

Combining the Benders Decomposition with the truck based column generation, we get 
$$O(\sum_{k = 1}^{n/p}s\cdot\sum_{t=1}^k\pi_t(n_t)) $$
$$O(s\cdot\sum_{k = 1}^{n/p}\sum_{t=1}^k\pi_t(n_t)) $$




### Space Complexity

Variables and constraints are in $O(n^2)$. However we face memory issues when it comes to storing the objective function of the dual which is in space complexity $O(m\cdot n^2)$ which can be impractical without a column generation method.

### Pros and Cons

In this method the subproblems are easier to solve and of reduced size compared to the Uzawa algorithm. However memory issues have to be dealt with.

## Placement subproblem analysis

The placement subproblem consists in optimizing the placement of items into stacks and stacks into a truck, satisfying specific constraints on volume, load order, etc...


### Description

Formulation:

**Minimize:**

$$\alpha_T c^\top_{{T}_1} \zeta^T + \alpha_E c^\top_{{T}_2} \zeta^E + \alpha_Ic^\top_{I} (IDL - TI^{t\top} \times TDA) + \kappa^t(TI^t - \overline{TI})$$

**Subject to:**

$$-(1 - \zeta^T) \times M^\zeta + 1 \leq TI^T \times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]\leq \zeta^T \times M^\zeta$$


$$-(1 - \zeta^E) \times M^\zeta + 1 \leq TI^E \times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]\leq \zeta^E \times M^\zeta$$

$$TI^t \leq TR\quad \bold{(TI_1)}$$

$$TI^{t\top}\times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] = \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right]\quad \bold{(TI_2)}$$

$$S^\top \cdot \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] = TI^t[t]$$

$$Z\leq S\times M^{Z}\quad \bold{(Z_1)}$$  

$$S\times IS = Z\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\quad \bold{(S_1)}$$


$$-M^{Z}(1-S) \leq Z - \left [  SS  \cdots  SS  \right ] \leq M^{Z}(1-S)\quad \bold{(Z_2)}$$


$$Q\leq SU\times M^Q\quad \bold{(Q_1)}$$  
$$S\times IU = Q \quad \bold{(Q_2)}$$


$$-M^Q(1-SU) \leq Q - \left [S\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \dots S\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]\right ] \leq M^Q(1-SU)\quad \bold{(Q_3)}$$ 


$$V\leq SK\times M^V\quad \bold{(V_1)}$$  
$$S\times IK = V\quad \bold{(V_2)}$$

$$-M^V(1-SK) \leq V - \left [  S \times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \cdots  S \times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \right ] \leq M^V(1-SK)\quad \bold{(V_3)}$$

$$W\leq SG\times M^{W}\quad \bold{(W_1)}$$  
$$\displaystyle S\times IPD = W\quad \bold{(W_2)}$$

$$-M^{W}(1-SG) \leq W - \left [  S \times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \cdots  S \times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \right ] \leq M^{W}(1-SG)\quad \bold{(W_3)}$$


$$G^l\leq S\times M^G$$  
$$G^r \leq S\times M^G$$  
$$G^l\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = G^r\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$-M^G(1-S) \leq G^r - \left [  SO  \cdots  SO  \right ] \leq M^G(1-S)$$  
$$-M^G(1-S) \leq G^l - \left [  \begin{matrix}IOV^\top \\ \vdots \\ IOV^\top \end{matrix} \right ] \leq M^G(1-S)$$

Define $SL$ and $SW$


$$D^L \leq S \times M^{D^L}$$

$$-M^{D^L}(1 - S) \leq D^L - \left [ SL \cdots SL\right ] \leq M^{D^L}(1 - S)$$

$$D^L\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = S\times IL$$

$$D^W \leq S \times M^{D^W}$$

$$-M^{D^W}(1 - S) \leq D^W - \left [ SW \cdots SW\right ] \leq M^{D^W}(1 - S)$$

$$D^W\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = S\times IW$$


$$SO\cdot M^{TL} \leq SX^e - SX^o - SL \leq SO\cdot M^{TL}$$

$$SO\cdot M^{TW} \leq SY^e - SY^o - SL \leq SO\cdot M^{TW}$$

$$-(1-SO)\cdot M^{TW} \leq SX^e - SX^o - SW \leq (1-SO)\cdot M^{TW}$$

$$-(1-SO)\cdot M^{TL} \leq SY^e - SY^o - SW \leq (1-SO)\cdot M^{TL}$$

$$SZ^o = 0$$

$$SZ^e = S\cdot IH$$

$$SX^e \leq  TL[t]$$  


$$SY^e \leq  TW[t]$$  
$$SZ^e \leq  TH[t]$$

$$\left [ \begin{matrix} 1 & 0 & \dots & 0\\ & I & &  \end{matrix} \right ]\times SX^o \leq  SX^o$$

$$\Xi^2 SX^o - \Xi^1 SX^{e} - \beta^- + \beta^+ =  - \epsilon\quad \bold{(\Xi_a)}$$

$$\beta^- \leq \lambda M^\lambda$$  
$$\beta^+ \leq (1-\lambda)M^\lambda$$  

$$(1 - \mu) \leq \beta^- M^\mu$$

$$\Xi^1SY^e \leq \Xi^2SY^o +  \xi M^{TW} + (1-\mu)M^{TW}\quad \bold{(\Xi_b)}$$

$$\Xi^2SY^e \leq \Xi^1SY^o + (1 - \xi)M^{TW} + (1-\mu)M^{TW}\quad \bold{(\Xi_c)}$$

$$\Xi^1SU\cdot TE[t] \leq \Xi^2 SU \cdot TE[t]\quad \bold{(\Xi_d)}$$

$$\Xi^1SU - \Xi^2SU \geq \chi\epsilon - rM^{TE} - (1 - \sigma^1)M^{TE}\quad \bold{(\Xi_e)}$$  
$$\Xi^2SU - \Xi^1SU \geq (1 - \chi)\epsilon - rM^{TE} - (1 - \sigma^1)M^{TE}\quad \bold{(\Xi_f)}$$

$$\Xi^2SK \cdot TKE[t] \geq \Xi^1SK \cdot TKE[t] - (1-r)M^{TKE}\quad \bold{(\Xi_g)}$$

$$\Xi^1SK \cdot TKE[t] - \Xi^2SK \cdot TKE[t] \geq \chi\epsilon - (1 - \sigma^2)M^{TKE}\quad \bold{(\Xi_h)}$$  
$$\Xi^2SK \cdot TKE[t] - \Xi^1SK \cdot TKE[t] \geq (1 - \chi)\epsilon - (1 - \sigma^2)M^{TKE}\quad \bold{(\Xi_i)}$$

$$\Xi^2SG\cdot TGE[t] \geq \Xi^1SG\cdot TGE[t] - (1 - \sigma^3)M^{TGE}\quad \bold{(\Xi_j)}$$

$$\sigma^1 + \sigma^2 + \sigma^3 \geq 1$$  


### Analysis

### Pros and Cons

## Comparison
