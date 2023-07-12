# Placement Algorithms

The 2023 ROADEF challenge presents an affectation problem for the logistic chain of Renault.
One key component of the problem is to be able to place optimaly stacks made up of items
into trucks which have predetermined dimensions and constraints.

This document reviews some ideas of algorithms for the resolution of such problem, which
could be formalized as follow:
Given a truck's width, and a list of stacks or items, place the stacks/items as 
to satisfy dimension, weight, loading order and other constraints and so, minimizing 
the total length of the objects, in order to make them fit into the dimensions of the truck.

## Corner based algorithm

This is an online algorithm which places boxes going from the cabin to the opening of the truck.

```Text

Let be O a list of available positions. The list will be updated and is made up 
of corners created by adjacent stacks.
Initially, there is only one corner considered (or two), in position (0, 0) 
(or (0, W)), which is adjactent to the cabin and to the left side of the truck.

As long as there are stacks to place:
    Put the stack in the corner the closest to the cabin and orient the stack in order to minimize the total length of items.
    Placing the stack creates new available corners and removes covered corner.

```

With no other constraint than no overlapping of two boxes, the algorithm can be proved to have a 2-OPT optimality guarantee.
Indeed, consider an instance with a truck of width $W$ and length $L$. Let have $L = 2W$. Consider 3 stacks:
The first and second stack are almost squares of dimensions $W$. Indeed, they have one side shorter by an $\epsilon$. The third stack has a width of $\epsilon$ and a length of $L$.
The optimal configuration of the placement is to have the two almost squares aligned, touching one side of the truck, and oriented length wise as to leave a space of
width $\epsilon$ and length $L$ adjacent to the other side of the truck. This space can be taken by the third stack. The truck is completely filled with total length of $L$.

The corner based algorithm takes the first almost square and place it width wise as to minimize total length. It then takes the second almost square and does the same.
At this point, the partial solution has a length of $2W - 2\epsilon = L - 2\epsilon$, but there is no room for the third stack which is forced to be added length wise to the solution. We consequently have a solution of value $L - 2\epsilon + L$ which is roughly equal to $2L$ if we ignore $\epsilon$. The algorithm is thus at least 2-OPT in the worst case.

It would be interesting to be able to prove that the algorithm can't do worse than 2-OPT. We would need to prove that there is no instance for which the algorithm provides a solution worse than $2L$. We would need to prove that:

$!(\exist i \in I, algo(i) > 2L)$

or $\forall i \in I, algo(i) \leq 2L$

or prove that and instance i with value > 2 is impossible.

Even if we can prove this, this doesn't tell us if the algorithm adapts well with new constraints such as those of loading order, etc.

Proof:
Let's suppose that we indeed have an instance for which the algorithm finds a solution of value greater than 2-OPT. It means that we are in one of two cases:

- Case 1: the algorithm returns a solution in which there is one item of greatest length that is greater than 2-OPT. It is not possible because the optimum is of value at least the length of this item.
- Case 2: the algorithm returns a solution in which there are a number m of items adjacent aligned in the length axis which have a sumed size in the length axis greater than 2-OPT.
  - Case 2.1: the m items are all oriented length wise. In this case, the algorithm wasn't able to orient them width wise. It means the remaining space when these items were placed did not allow for orienting these width wise.
    - Case 2.1.1: the m items have a length greater than W.
      - Case 2.1.1.1: the m items have a width that allows some of them to be placed side to side width wise.
        - Case 2.1.1.1.1: the width of the truck W is less than the sum of the widths of two of these items.
          - Case 2.1.1.1.1.1: the width of the truck is exactly equal to the width of one item of the queue. In this case the optimal solution is equal to the solution found by the algrithm because there is no other way of placing the stacks.
          - Case 2.1.1.1.1.2: the width of the truck allows some space along side these items. But these items will always be stacked lengthwise anyways because none can be placed side to side. So the optimal solution has the same value than the one returned by the algorithm.
        - Case 2.1.1.1.2: at least two of these items can be placed side to side. Since the algorithm did not do so, it means some other items were already taking that space. These items have a width which is less than the width of the trucks minus the width of one long item. The algorithm always places stacks as to minimize the total length. It means the difference in length between the long queue and the blocking queue can't be greater than the length of the last item added to the long queue. Since the total length of this solution is more than two times the length of the optimal solution... etc

Or we could prove that the instance we found earlier is the worse.


If we adapt this algorithm to the loading order constraint, it requires to order the input list of stacks first by supplier, then by supplier dock for each supplier, and finally by plant dock for each supplier dock. The algorithm then builds a valid solution because it always tries to place stacks in positions the closest to the cabin. To be proved.

What about weight constraints?

The algorithm is not deterministic. Actually, the order in which stacks are provided can change the solution. Similarly, the way new corner positions are added to the list can affect the result of the algorithm. We can thus consider a tree that represents the different paths the algorithm can take to build its solution. The different branching possibilities are means to adapt a solution in regards to the weight constraints.
We can extend the algorithm to, when given a current stack, first explore the different available corner options and if none satisfies the constraint, try to place the following stack in the input list, hoping the constraint will be satisfied. If no similar stack is better to consider, rollback to the last stack placed, remove it and the previous one and try to place them in inverted order. If it still doesn't work, repeat the rollback process.

What it does concretely is exhausting all permutations of the list within what is allowed to satisfy the loading order. The permutations thus only affect stacks of same supplier, supplier dock and plant dock. Let's say there are k groups of similar stacks in terms of supplier, supplier dock and plant dock, for n stacks. In the worst case there are k-1 groups of exactly one stack and one group containing n-k-1 stacks. We thus have (n-k-1)! combinations to examine in the worst case.

Does the algorithm allows all possible solutions? Can the algorithm get stuck on incorrect solutions while not exploring valid solutions due to its design? This question did not occur until we tackled the weight constraints, although it can also be detrimental to optimality. Is a stack always in a corner? Well normaly, placing stacks in corners does not prevent the algorithm from testing feasible conditions. However, all corners should be considered, not only the (0, 0) one but also the (0,W) one and all similar corners which correspond not to the origin of the stack but also to its upper left corner. There are also cases in which a stack could be placed adjacent to a side, but with no side touching above, and a stack touching below but which have a starting coordinate greater than the current stack. In this particular case, the current stack is not placed in a corner, but this configuration could be part of a solution. So we need to adapt the algorithm to take into account other types of corners, like this "floating" one.

> Bottom-Left Placement Theorem for Rectangle Packing, W. Huang, T. Ye, D. Chen, 22 juil. 2011, arXiv:1107.4463v1

The corner based algorithm resembles a lot the Bottom-Left algorithm presented in this article.

The Bottom-Left algorithm has the benefits of providing bottom-left stable packings, ensuring stacks are always adjacent to stacks or side horizontally.

## Local search algorithm

I can't think of any local search based algorithm for this problem considering our constraints.

## Most constraints first

Put the stacks with most constraints and then play around as you add new stacks.

## The ARC Project

[The ARC project](https://intranet.csc.liv.ac.uk/~epa/surveyhtml.html)
contains a survey on two dimensional packing algorithms.

- First-Fit Decreasing Height (FFDH) algorithm
- Next-Fit Decreasing Height (NFDH) algorithm
- Best-Fit Decreasing Height (BFDH) algorithm
- Bottom-Left (BL) Algorithm which resembles a lot our own algorithm
- Baker's Up-Down (UD) algorithm 
- Reverse-fit (RF) algorithm 
- Steinberg's algorithm 
- Split-Fit algorithm (SF)  
- Sleator's algorithm  

However these algorithms might be difficult to adapt to the context of the loading order constraints and the weight constraints.

## Best fitting space

The idea is, at each iteration, to segment the remaining available space into largest rectangles, and then choosing the item to place into it which fits the best.

## Those algorithms

The algorithms seen so far don't consider modifying stacks composition to accomodate weight constraints.