# Challenge ROADEF 2022

The 2022 ROADEF challenge presents an assignment problem for the logistic chain of Renault. As a young graduate of a computer science master degree with a specialization in Artificial Intelligence and Operational Research as of 2023, I present in this project my work to analyse, model, and solve the problem formulated by Renault. However, I am not a competing participant, I only work on this project for training as well as exploratory purposes.

$\rightarrow$ [Website of the challenge](https://www.roadef.org/challenge/2022/en/index.php).

## Brief summary of the subject

The objective is to pack a set of items from suppliers into stacks and to pack the stacks into trucks which deliver the plants, in order to minimize (a) the number of the trucks used and (b) the inventory in the plants due to early deliveries. Items have various properties (a time window of delivery, weight, dimensions...). 

Renault hired trucks from transporters to deliver the items from the suppliers to the plants. These hired trucks are defined in an annual planning and are called "planned" trucks. Trucks have properties which translate into constraints (maximum weight on axles, dimensions). All items must be delivered in time to their plants. If necessary, "extra" trucks can be called from the transporter to deliver items which could not fit into planned trucks, at an extra cost.

## Application

The ultimate goal of this project is to implement a julia application allowing to solve efficiently the truck loading optimization problem as described by Renault.

## Status of the project

The most up-to-date branch is [the feature branch adding weight constraints to the model for solving the single-truck placement problem](https://github.com/BrcRs/truck-stack-item-optimization/tree/wip-30-satisfying-weight), still in progress. The difficulty lies in adapting the placement algorithm to continue outputting interesting solutions while accomodating for weight constraints.

The most up-to-date working branch, the [develop branch](https://github.com/BrcRs/truck-stack-item-optimization/tree/develop), hosts a working algorithm (along with necessary framework) for solving the problem of packing efficiently items into stacks and stacks into a single truck as to satisfy most of stack constraints and supplier/plant order constraints.

## Current content of the main branch

The main branch's `src/` code is outdated. However it contains some interesting documentation work:

- _Comparative Analysis of Solving Methods_ aims to work on and compile analyses of some resolution methods for the logistic problem.

- _Placement Algorithms_ presents different algorithms to solve the placement of stacks into a singular truck, with time and space analysis, aswell as optimality bounds analysis.

- _Modelisation_ presents a modelisation of the problem as a mixed-integer problem.



## Roadmap

### Part 1: A Working Single-Truck Packing Algorithm

Work in progress, the first iteration of the project have the following objectives:

- Implement an instance generator for the packing problem into a single truck. Generated instances have optimal solutions which fill up a truck entirely.
- Have a working packing algorithm which packs as many items as possible into stacks and then into a single truck. The efficiency of the algorithm is measured by allowing items to overflow out of the back of the truck and comparing the length of the overflow to the length of the truck. The algorithm can be proved to have an experimental optimality bound of 3-OPT on 1000 instances. This means solutions output by the algorithm have a length of at most 3 times that of the truck.
- Feature a simple packing visualizer.

### Part 2: A Working Multi-Truck Assignment Algorithm (no deadline expected)

- Solve all instances provided by Renault in reasonable time.
- Provided constraint checkers pass successfully.

## Contributors

I (@BrcRs) am the sole contributor of this project.
