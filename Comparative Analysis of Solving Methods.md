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
### Description
### Analysis
### Pros and Cons
## Comparison