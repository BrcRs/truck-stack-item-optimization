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

### Analysis

In the worst case, the algorithm adds one truck to the problem at each iteration. We would have thus $m$ iterations, $m$ being the maximum number of trucks.
The maximum number of trucks can be expressed in function of the number of items. Indeed, in the worst case we have one truck per item.

The number of iterations is $n$. But the size of the instance solved will increase with each iteration:

$$\sum_{k = 1}^{n}\phi(n, k)$$

This algorithm has thus a time complexity of $O(\sum_{k = 1}^{n}\phi(n, k))$, $\phi$ being the time complexity of the solving algorithm used internaly.

However in practice, we can presume that the actual number of trucks used $\hat{m}$ is significantly smaller than $n$. We could say that $p\times\hat{m} \leq n$

We then have a time complexity of $O(\sum_{k = 1}^{n/p}\phi(n, k)) \approx O(n\cdot \phi(n, n))$.

### Pros and Cons
## Uzawa Inspired Algorithm
### Description
### Analysis
### Pros and Cons
## Benders Decomposition
### Description
### Analysis
### Pros and Cons
## Branch & Bound for the placement subproblem
### Description
### Analysis
### Pros and Cons
## Comparison