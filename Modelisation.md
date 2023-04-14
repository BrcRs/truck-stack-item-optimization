# ROADEF 2022 Modelisation

Minimize nb of trucks used and inventory of the plants.

We will simplify the subproblem consisting in optimizing the stacks within the trucks.

We will only consider a number of sets of items, ordered in each truck according to the loading order between the picked-up suppliers.

We will consider the volume used by each items combined within each truck (and mass) to represent the capacity constraints of each truck.

Each truck $t$ has a maximum volume $TV_t$ which is equal to $TL_t \times TWt \times TH_t$.

N Items  
M Suppliers  
P Plants  

To each item corresponds exactly one plant.

## First Formulation

$x_{{T}_1} \in \{0, 1\}^{|T|}$  
$x_{{T}_2} \in \bold{N}^{|T|}$  
$y \in \bold{N}^{N}$

$x_{{T}_1}$ : planned trucks used or not  
$x_{{T}_2}$ : number of extra trucks used  
$y$ : arrival time of each item  
$d$ : latest arrival time authorized for each item  

$$ min\quad \alpha_T c^\top_{{T}_1} x_{{T}_1} + \alpha_E c^\top_{{T}_2} x_{{T}_2} + \alpha_Ic^\top_{I} (d - y) $$

How to indicate if truck $t$ contains item $i$? I see 2 solutions:

The first solution consists in considering a $T \times N$ matrix made of 0's and 1's, a 0 at $(t,i)$ meaning truck $t$ does not contain item $i$ and inversely if there is a 1.

The second solution consists in considering a vector of size $N$, indicating the index of the truck that contains item $i$.

The problem is: the amount of trucks used is variable... especially when it comes to extra ones.

Or you could prepare the matrix with extra trucks from the start even though a majority of those won't be used. You'd have then to manage the great number of potentially useless variables added (delayed column generation?).

Let's note $TI$ the truck $-$ item matrix.  

Another thing is not every item can be loaded in every truck. So the TI matrix would be multiplied by a "*filter*" matrix. Or maybe it is another constraint instead.

$TI \leq F$

with $F$ the candidate matrix.
Also, an item can't be loaded into 2 trucks:


$TI^\top\times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] \leq \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right]$

Explanations:  
$TI^\top : \begin{matrix}\text{item 1}\\ \text{item 2}
\end{matrix}\left [\begin{matrix} 1 & 0\\ 0 & 1\\ \end{matrix}\right ] \left [ \begin{matrix} 1 \\ 1 \end{matrix} \right ] = \left [ \begin{matrix} 1 \\ 1 \end{matrix} \right ]$

$TI^\top : \begin{matrix}\text{item 1}\\ \text{item 2}
\end{matrix}\left [ \begin{matrix} 1 & 0 \\ 1 & 1 \end{matrix} \right ]\left [ \begin{matrix} 1 \\ 1 \end{matrix} \right ] = \left [ \begin{matrix} 1 \\ 2 \end{matrix} \right ]$

### Item constraints

---

#### **(I1)** Stack constraints

---

Stack constraints are ignored and replaced by volume (and weight?) constraints.

The sum of the volumes of contained items by truck $t$ should not exceed volume of truck $t$.

$TI\times IV \leq TV$  
$TI\times IW \leq TW$

with $TW = EM^{mr} + EM^{mm} + EM^{mh}$


#### **(I2)** Arrival constraints

---

If Truck $t$ is loaded with $i$, it arrives at the item's plant $IP_i$. $TP$ is a vector.

$TI^\top \times TP = IP$

$TI^\top : \begin{matrix}\text{item 1}\\ \text{item 2} \\\text{item 3}
\end{matrix}\left [\begin{matrix} 1 & 0\\ 0 & 1\\ 1 & 0 \end{matrix}\right ] \left [ \begin{matrix} 5 \\ 2 \end{matrix} \right ] = IP$

$\left [\begin{matrix} 1 & 0\\ 0 & 1\\ 1 & 0 \end{matrix}\right ] \left [ \begin{matrix} 5 \\ 2 \end{matrix} \right ] = \left [ \begin{matrix} 5 \\ 2 \\ 5\end{matrix} \right ]$

#### **(I3)** Pick-up constraints

---

$TI \leq F$

$F$ is infered from candidate products.

#### **(I4)** Stop by suppliers constraints

---

$TI^\top \times TU = IU$ ?

$TU$ : $truck \times suppliers$

$TI^\top : \begin{matrix}\text{item 1}\\ \text{item 2} \\\text{item 3}
\end{matrix}\left [\begin{matrix} 1 & 0\\ 0 & 1\\ 1 & 0 \end{matrix}\right ] \left [ \begin{matrix} 5 & 4 \\ 2 & 8 \end{matrix} \right ] = IU$

$\left [\begin{matrix} 1 & 0\\ 0 & 1\\ 1 & 0 \end{matrix}\right ] \left [ \begin{matrix} 5 & 4 \\ 2 & 8 \end{matrix} \right ] = \left [\begin{matrix} 5 & 4\\ 2 & 8\\ 5 & 4 \end{matrix}\right ]$

$\left [\begin{matrix} 5\\ 8\\ 4 \end{matrix}\right ]$

$\left [\begin{matrix} 5y_{1,1} & 4y_{1,2}\\ 2 y_{2,1} & 8y_{2,2}\\ 5 y_{3,1} & 4y_{3,2} \end{matrix}\right ]\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = \left [\begin{matrix} 5\\ 8\\ 4 \end{matrix}\right ]$

$\left ( (TI^\top\times TU)  \circ \Theta\right ) \times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = IU$
and  
$\Theta\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \geq \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$

> (Impossible formulation because we multiply variables? TI and \Theta?).

With $y \in \{0, 1\}$. This formulation implies that for each row, at least one of the constraints corresponding to one column (supplier) is satisfied.

The $TI\circ B$ notation corresponds to the Hadamard matrix product (element-wise product).

$a^\top_{j,1}tu_{1, u}y_{j,u} + a^\top_{j,2}tu_{2, u}y_{j,u} + \dots + a^\top_{j,T}tu_{T, u}y_{j,u}$

$\displaystyle \sum_u\left (y_{j,u}\sum_t a^\top_{j, t}tu_{t,u}\right ) = iu_j,\quad \forall j$

Linearization:

$r_{j,u} \leq y_{j,u}M^{I4}$  
$(0 - M^{I4})(1 - y_{j,u}) \leq r_{j,u} - \sum_t a^\top_{j, t}tu_{t,u} \leq M^{I4}(1 - y_{j,u})$  
$y_{j,u}=0\implies-M^{I4} \leq r_{j,u} - \sum_t a^\top_{j, t}tu_{t,u}\leq M^{I4} \quad\land\quad r_{j,u} \leq 0$  
$y_{j,u}=1\implies 0 \leq r_{j,u} - \sum_t a^\top_{j, t}tu_{t,u}\leq 0$

The constraint becomes:

$\displaystyle \sum_ur_{j,u} = iu_j,\quad \forall j$

Matrix formulation:

$R \leq \Theta M^{I4}$  
$ -M^{I4}(1-\Theta) \leq R - (TI^\top\times TU) \leq M^{I4}(1 - \Theta) $  
<!-- $ M^{I4} = \max(R) = \max(TI^\top\times TU) = \left [ \begin{matrix}1 & 1 & \dots\\1 & 1 & \dots\\ \vdots&  \vdots&\ddots\end{matrix}\right ]\cdot TU$  

or -->

$\displaystyle M^{I4}_{i,j} = \max_i TU_{i,j}\quad \forall j \in  \left \{1,\; Col(TI^\top)\right \}$

$R\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]=IU$  

$\Theta\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \geq \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$

<!-- Is there another (simpler) way?

$\left [\begin{matrix} 5 & 4\\ 2  & 8\\ 5  & 4 \end{matrix}\right ] - \left [\begin{matrix} 5& & 5\\ 8 & \cdots & 8\\ 4 & & 4 \end{matrix}\right ] = \left [\begin{matrix} 0& -1\\ -6 & 0\\ 1 & 0 \end{matrix}\right ]$ -->

#### **(I5)** Time Window

---

$IDE \leq TI^\top \times TDA \leq IDL$

#### Stack constraints

---
All stack **(S\*)** constraints are ignored.

#### Placement constraints

---
Ignored. Loading orders aren't significant when the placement constraints are reduced to overall weight and volume.


Maybe the placement constraints are the actual heart of the problem. If I want a simpler problem to solve, maybe instead of simplifying the placement constraints I should to the contrary: consider only the subproblem consisting in placing optimally items in a truck and ignore everything else.

$\begin{matrix}\\ \min \\ \small items \end{matrix}\; \begin{matrix}\\ \min \\ \small items\;coord \end{matrix}\quad placement(items, items\; coord, truck )$

## The heart of the problem: placement

In this section, let's consider the subproblem consisting in, given a set of items and a truck (and all its properties), how to optimally place the items as to satisfy the placement constraints. Benders could then be used with the master problem having to determine which items go in each truck and the subproblem optimizing the placement. The placement problem has no objective, which makes it only a constraint satisfaction problem, much simpler? to solve. Or have an objective function to minimize the number of items which could not be included.
However I feel like the placement depends too much on the items chosen in each truck... Benders can't be used... or maybe I am thinking too much.

### Stacks

Which item goes into which stack? A matrix (similar to the one used for trucks/items) could be used:

$S : |stacks| \times N$ binary matrix

Similarly, we won't know in advance the number of stacks necessary.

We also need to define a similar matrix that indicates which plant docks are concerned by stack S:

$\widetilde{SG} : |stacks| \times | docks|$ binary matrix.

But that stack thing might change the way we write the problem.

Docks are unique by plant, but not in the whole problem.

<!-- S: number of stacks


Parameters:
A Truck (dimensions, suppliers and supplier docks served, plant and plant docks served), a number of items. -->

Each item is packed in exactly one stack:

$S^\top\times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] = \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]$

Each stack is loaded into exactly one truck:  
$ST:stacks\times trucks$, a binary matrix

$ST^\top\times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] = \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]$

Do we need to link the item's truck with the stack's? Else the program could put the item in a truck but place its stack in another...

<!-- $ST : \left [ \begin{matrix}1 & 0\\0 & 1 \\ 1 & 0 \end{matrix} \right ] = TI^\top : \left [ \begin{matrix}1 & 0\\0 & 1 \\ 1 & 0 \end{matrix} \right ]$ -->

We have to infer the truck of a stack from the items the stack contains.
The truck of a stack is no longer a variable.

It is easily done by the following multiplication:

$ST\circ \left ( S\times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]\right ) = S\times TI^\top$

















$\left [ \begin{matrix} 1 & 0 & 1 & 0\\ 0 & 1 & 0 & 0\\ 0& 0 & 0 & 1 \end{matrix}\right ]\times \left [ \begin{matrix} 1 & 0 & 0 & 0\\ 0 & 0 & 1 & 0\\ 1& 0 & 0 & 0 \\0&1&0&0\end{matrix}\right ] = \left [ \begin{matrix} 2 & 0 & 0 & 0\\ 0 & 0 & 1 & 0\\ 0& 1 & 0 & 0 \end{matrix}\right ] = \left [ \begin{matrix} 1 & 0 & 0 & 0\\ 0 & 0 & 1 & 0\\ 0& 1 & 0 & 0 \end{matrix}\right ] \circ  \left [ \begin{matrix} 2\\ 1\\ 1\end{matrix}\right ]$

But we can't multiply variables! Which means... we can't use this...


Let's linearize:

$st_{i,j}\times\sum_{j'} s_{i, j'} = s_{i,1}ti_{j, 1} + s_{i,2}ti_{j, 2} + \dots + s_{i,N}ti_{j, N} \quad \forall i \in \{1, |stacks|\}$

$\psi_{i,j} = \omega_{i,j, 1} + \omega_{i,j, 2} + \dots + \omega_{i,j, N} \quad \forall i \in \{1, |stacks|\}$

$\omega_{i,j,\iota} \leq s_{i,\iota}M^\omega\quad\implies\quad s_{i,\iota}=0\implies\omega_{i,j,\iota}=0$

$-(1-s_{i,\iota})M^\omega \leq \omega_{i,j,\iota} - ti_{j,\iota} \leq (1-s_{i,\iota})M^\omega$

$\psi_{i,j} \leq st_{i,j}M^\psi$

$-(1-st_{i,j})M^\psi \leq \psi_{i,j} - \sum_{j'} s_{i, j'} \leq (1-st_{i,j})M^\psi$


$\Omega[i,j,\iota] : stacks-truck-item$ matrix 

Matrix formulation:

<!-- $ST = \Omega\overset{items}{\cdot} \left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ]$ -->

$\Psi = \left [ \Omega[\dots, 1, \dots]\cdot \left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ] \Omega[\dots, 2, \dots]\cdot \left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ]\dots  \right ]$

$\Omega[\dots, j, \dots]\leq S\cdot M^\Omega\quad \forall j$

$-(1-S)M^\Omega \leq \Omega[i,\dots,\dots] - TI^\top \leq (1-S)M^\Omega\quad \forall i$

$\Psi \leq ST \cdot M^\Psi$

$-(1-ST)M^\Psi \leq \Psi - \left [S\cdot\left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ]\dots\right ] \leq (1-ST)M^\Psi$

$M^\Omega = 2$

$M^\Psi = |items|$

<!-- Linearization (wrong):

$ST \leq S\cdot M^{ST}$

$-(1-S)M^{ST} \leq ST - TI^\top \leq (1-S)M^{ST}$

$M^{ST} = 2$ -->

All the items packed in a stack must share the same supplier, plant, stackability code and supplier dock.

We do have access to the item $-$ stack matrix, so this shouldn't be too difficult.

$\displaystyle S\times IS = \left ( S \times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] \right )^\top \times IS$

Doesn't work because of the left hand side.

We need a new variable for the stackability code of each stack.
A simple vector:

$SS$

Then:


$\displaystyle S\times IS = S \times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]  \circ SS$

Valid

$\displaystyle \left [ \begin{matrix}1 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1 \end{matrix} \right ]\times \left [ \begin{matrix}1\\1\\3\\4 \end{matrix} \right ] =\left  [ \begin{matrix}1 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1 \end{matrix} \right ] \times \left [ \begin{matrix}1\\1\\1\\1 \end{matrix} \right ] \circ \left [ \begin{matrix}1\\3\\4 \end{matrix} \right ]$

$\iff \displaystyle \left [ \begin{matrix}2\\3\\4 \end{matrix} \right ] = \left [ \begin{matrix}2\\1\\1 \end{matrix} \right ] \circ \left [ \begin{matrix}1\\3\\4 \end{matrix} \right ]$

Invalid

$\displaystyle \left [ \begin{matrix}1 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1 \end{matrix} \right ]\times \left [ \begin{matrix}1\\2\\3\\4 \end{matrix} \right ] = \left [ \begin{matrix}1 & 1 & 0 & 0\\0 & 0 & 1 & 0\\0 & 0 & 0 & 1 \end{matrix} \right ] \times \left [ \begin{matrix}1\\1\\1\\1 \end{matrix} \right ] \circ \left [ \begin{matrix}1\\3\\4 \end{matrix} \right ]$

$\iff \displaystyle \left [ \begin{matrix}3\\3\\4 \end{matrix} \right ] \neq \left [ \begin{matrix}2\\1\\1 \end{matrix} \right ] \circ \left [ \begin{matrix}1\\3\\4 \end{matrix} \right ]$

Buuut this formulation is impossible because we multiply $S$ with $SS$ which are both variable matrices...

$\displaystyle s_{i,1}is_1 + s_{i,2}is_2 + \dots + s_{i,N}is_N = \sum_{1\leq j \leq N} s_{i,j}ss_i \quad\quad \forall i \in \{1, \dots, |stacks|\}$

$\displaystyle s_{i,1}is_1 + s_{i,2}is_2 + \dots + s_{i,N}is_N = s_{i,1}ss_i + s_{i,2}ss_i + \dots + s_{i,N}ss_i \quad\quad \forall i \in \{1, \dots, |stacks|\}$

Linearization:

$\displaystyle s_{i,1}is_1 + s_{i,2}is_2 + \dots + s_{i,N}is_N = z_{i,1} + z_{i,2} + \dots + z_{i,N} \quad\quad \forall i \in \{1, \dots, |stacks|\}$
$z_{i,j} \leq s_{i,j}M^{Z}$  
$(0 - M^{Z})(1 - s_{i,j}) \leq z_{i,j} - ss_i \leq M^{Z}(1 - s_{i,j})$  
$s_{i,j}=0\implies-M^{Z} \leq z_{i,j} - ss_i\leq M^{Z}$  
$s_{i,j}=1\implies 0 \leq z_{i,j} - ss_i\leq 0$

Matrix formulation:

$Z\leq S\times M^{Z}$  
$S\times IS = Z\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$

$-M^{Z}(1-S) \leq Z - \left [  SS  \cdots  SS  \right ] \leq M^{Z}(1-S)$

$M^Z = \max IS + 1$

Is it possible, with those constraints, that a stack doesn't have the same stackability code than the items? Answer: no. Proof:

$\displaystyle s_{i,1}is_1 + s_{i,2}is_2 + s_{i,3}is_3 = s_{i,1}ss_i + s_{i,2}ss_i + s_{i,3}ss_i \quad\quad \forall i \in \{1, \dots, |stacks|\}$  
$\implies is_j = ss_i \quad\forall j\quad \square$

All the items packed in a stack must share the same supplier, plant and supplier dock.

We reason similarly. Let SU be a variable-vector attributing a supplier to each stack.

$\displaystyle S\times IU = S \times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]  \circ SU$

Matrix formulation of the linearization:

$Q\leq S\times M^Q$  
$S\times IU = Q\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$

$-M^Q(1-S) \leq Q - \left [  SU  \cdots  SU  \right ] \leq M^Q(1-S)$

$M^Q = \max IU + 1$

All the items packed in a stack must share the same plant and supplier dock.

$\displaystyle S\times IP = S \times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]  \circ SP$

Matrix formulation of the linearization:

$H\leq S\times M^H$  
$S\times IP = H\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$

$-M^H(1-S) \leq H - \left [  SP  \cdots  SP  \right ] \leq M^H(1-S)$

$M^H = \max IP + 1$

All the items packed in a stack must share the same supplier dock.

$\displaystyle S\times IK = S \times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]  \circ SK$

Matrix formulation of the linearization:

$V\leq S\times M^V$  
$S\times IK = V\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$

$-M^V(1-S) \leq V - \left [  SK  \cdots  SK  \right ] \leq M^V(1-S)$

$M^V = \max IK + 1$

For any stack $s$ packed into truck $t$, if $(TF_t = \text{no})$, then all the items of stack $s$ must share the same plant dock **(S2)**.

We have a similar thing going on.

$\displaystyle S\times IPD = S \times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]  \circ SPD$ if $TF=0$

becomes

$-M^{IPD}\times TF \leq \displaystyle S\times IPD - S \times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]  \circ SPD \leq M^{IPD}\times TF$

Matrix formulation of the linearization:

$W\leq S\times M^{W}$  
$-M^{IPD}\times TF \leq \displaystyle S\times IPD - W\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \leq M^{IPD}\times TF$


$-M^{W}(1-S) \leq W - \left [  SPD  \cdots  SPD  \right ] \leq M^{W}(1-S)$


For any stack $s$ packed into truck $t$, if ($TF_t = \text{yes}$), $s$ may contain items with 2 plant docks with consecutive loading orders **(S3)**.
There is a restriction: for every stackability code SC present among the stacks of a truck $t$, only one stack with the stackability code SC can contain items with 2 plant docks.
> Lets ignore this for now.

Since we ignored $TF$, the previous constraints become:

$\displaystyle S\times IPD = S \times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]  \circ SPD$


Matrix formulation of the linearization:

$W\leq S\times M^{W}$  
$\displaystyle S\times IPD = W\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$

$-M^{W}(1-S) \leq W - \left [  SPD  \cdots  SPD  \right ] \leq M^{W}(1-S)$

$M^W = \max{IPD} + 1$

If one item $i$ of a stack has a forced orientation $IO_i$, then all the items of stack $s$ must share the same orientation ($so_s = IO_i$) **(S4)**. Consequently, there cannot be 2 different forced orientations of items in the same stack.

$SO$: orientation of stacks (variable).  
$IO$: forced orientation of items (parameter).  
$IOV$: orientation of items (variable).

Let's make sure that every items in a stack have the same orientation:

$\displaystyle S\times IOV = S \times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ]  \circ SO$

Matrix formulation of the (both sides) linearization:

$G^l\leq S\times M^G$  
$G^r \leq S\times M^G$  
$G^l\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = G^r\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$

$-M^G(1-S) \leq G^r - \left [  SO  \cdots  SO  \right ] \leq M^G(1-S)$  
$-M^G(1-S) \leq G^l - \left [  IOV  \cdots  IOV  \right ] \leq M^G(1-S)$

$G^r = \max SO + 1$

$G^l = \max IOV + 1$

Then add constraints for items with forced orientations to take the value of their respective forced orientations.

For any stack $s$ packed into truck $t$, the total weight of the items packed above the bottom item associated with product $\gamma$ must not exceed the maximal weight $TMM_{t\gamma}$ :  
$$\displaystyle\sum_{i\in \widetilde{SI_s}\text{ and }i\neq bottom\_item} IM_i \leq TMM_{t\gamma}$$
**(S5)**.
> Ignored for now

The number of items packed into a stack $s$ must not exceed the smallest "max stackability" $ISM_i$ of the items $i$ present in stack $s$. **(S6)**
In other words, for each item $i$ of a stack $s$, the number of items of $s$ must not exceed $ISM_i$.
> Ignored for now

The density of a stack $s$ must not exceed the maximal stack density defined for the truck $t$ into which the stack $s$ is loaded : $$\frac{sm_s}{sl_s\times sw_s} \leq TEM_t$$ **(S7)**.
> Ignored

### Placement constraints in details

The placement of a stack $s$ into a truck $t$ must not exceed the truck’s dimensions :
$\forall t \in \widetilde{T}, \forall s \in \widetilde{TS_t},\quad sx^e_s \leq TL_t $ and $sy^e_s \leq TW_t$ and $sz^e_s \leq TH_t$ **(P1)**.

$SX^e$ : vector of variable x position of stacks  
$SY^e$ : vector of variable y position of stacks  
$SZ^e$ : vector of variable z position of stacks

Reminder: $ST$, is the stack $-$ truck matrix.

The constraints are:

$SX^e \leq  ST \times TL$  
$SY^e \leq  ST \times TW$  
$SZ^e \leq  ST \times TH$

Easy stuff

$SX^e - SX^o = SL + SO\cdot M$

$SY^e - SY^o = SW + SO\cdot M$

$SX^e - SX^o = SW + (1-SO)\cdot M$

$SY^e - SY^o = SL + (1-SO)\cdot M$

$SZ^o = 0$

$SZ^e = S\cdot IH$

> (Ignoring nesting heights).

The stacks packed into a truck $t$ cannot overlap :
$\forall t \in \widetilde{T}, \forall s_1, s_2 \in \widetilde{TS_t}$ with $$sx^o_{s_1} \leq sx^o_{s_2}$$
, if $$(sx^o_{s_2} < sx^e_{s_1})$$ then $$sy^o_{s_2} \geq sy^e_{s_1}$$
or $$sy^e_{s_2} \leq sy^o_{s_1}$$
**(P2)**
If $s_1$ and $s_2$ overlap in the X axis, then they cannot overlap in the Y axis.

We are looking for a truck $-$ stack coordinate matrix

$\alpha M$

$SX^o_{s_1} \leq SX^o_{s_1}$

$ST^\top\times SX^o$

3 trucks, 4 stacks  
$\displaystyle \left [ \begin{matrix}0 & 0 & 1 & 1\\1 & 0& 0& 0\\0 & 1 & 0 & 0 \end{matrix} \right ] \times \left [ \begin{matrix}0.5 & 0 & 0 & 0 \\0 & 9.7 & 0& 0\\0&0&4.9&0\\0&0&0&1.4 \end{matrix} \right ] = \left [ \begin{matrix}0 & 0 & 4.9 & 1.4 \\0.5 & 0 & 0 & 0\\0 & 9.7 & 0 & 0 \end{matrix} \right ]$

Set an order between stacks to simplify (without loss of generality).

$SX^o_{s_1} \leq SX^o_{s_2}$  
$SX^o_{s_2} \leq SX^o_{s_3}$  
$\dots$

$\left [ \begin{matrix} 1 & 0 & \dots & 0\\ & I & &  \end{matrix} \right ]\times SX^o \leq  SX^o$

Does it introduce problems between stacks of different trucks?
But it introduces a coupling constraint between all stacks! Which is pretty bad!

> Not between all stacks, but some stacks from different trucks are coupled. If we could identify those... and abolish the constraint when the two stacks are in different trucks! No...

$s_1 \leq s_2 \leq s_3 \leq s_4 \leq s_5 \leq s_6 \leq s_7$  


$\begin{matrix} s_1 & s_3 \\ s_2 & s_4 \\ s_5 & s_6 \\ s_7 \end{matrix}$


$s_1 \leq s_2 \quad s_3 \leq s_4 \quad s_5 \quad s_6 \quad s_7$  

Things would go better if only contiguous stacks were in each truck. Maybe we could ignore this for the moment and look for solution when resolving.

Now that there is an order between stacks, we know for which stacks we have $sx^o_{s_1} \leq sx^o_{s_2}$.

There would be $$\frac{(|stacks| - 1)\times |stacks|}{2}$$

constraints of the "no overlap" type, filtered by truck.

$\left [ \begin{matrix}1 & 1 & 0 & 0&\dots\\1 & 0& 1& 0&\dots\\1 & 0 & 0 & 1&\dots\\ \vdots\\ 0&1&1&0&\dots\\0&1&0&1&\dots\\etc \end{matrix} \right ]$

Let's note $\Xi^1$ : $\left [ \begin{matrix}1 & 0 & 0 & 0&\dots\\1 & 0& 0& 0&\dots\\1 & 0 & 0 & 0&\dots\\ \vdots\\ 0&1&0&0&\dots\\0&1&0&0&\dots\\etc \end{matrix} \right ]$

Let's note $\Xi^2$ : $\left [ \begin{matrix}0 & 1 & 0 & 0&\dots\\0 & 0& 1& 0&\dots\\0 & 0 & 0 & 1&\dots\\ \vdots\\ 0&0&1&0&\dots\\0&0&0&1&\dots\\etc \end{matrix} \right ]$


$\Xi^2 SX^o - \Xi^1 SX^{e} - \beta^- + \beta^+ = (\nu - | trucks|) M^{TL} - 0.0001\quad \bold{(a)}$

$M^{TL} = \max TL + 1$

$\Xi^2 SX^o - \Xi^1 SX^{e} - \beta^- + \beta^+ =  - 0.0001\quad \bold{(a)}$

With $\beta^- \geq 0$ and $\beta^+ \geq 0$.

If $sx^o_{s_2} - sx^e_{s_1} > 0$, then it would imply $\beta^- = sx^o_{s_2} - sx^e_{s_1}$ and $\beta^+ = 0$ for $\bold{(a)}$ to be satisfied.

If $sx^o_{s_2} - sx^e_{s_1} = 0$, then it would imply $\beta^- = \beta^+ = 0$ for $\bold{(a)}$ to be satisfied.

If $sx^o_{s_2} - sx^e_{s_1} < 0$, then it would imply $\beta^+ = sx^o_{s_2} - sx^e_{s_1}$ and $\beta^- = 0$ for $\bold{(a)}$ to be satisfied.

<!-- since $sx^o_{s_1} \leq sx^o_{s_2}$ and that $sx^e_{s_2} = sx^o_{s_2} + sw_s$ or $sx^e_{s_2} = sx^o_{s_2} + sl_s$ then $sx^e_{s_2} \geq sx^o_{s_2}$ and thus $sx^o_{s_1} \leq sx^e_{s_2}$.

This means that $sx^o_{s_1} - sx^e_{s_2} \leq 0$ and thus $sx^e_{s_2} - sx^o_{s_1} \geq 0$. -->

$\nu$ is a binary variable which would be equal to $|trucks|$ when the two considered stacks are from the same truck. $\nu \geq 0$ by definition (see further).

$(1 - \mu) \leq \beta^- M^\mu\quad\quad(\beta^- = 0)\implies (\mu=1)$

$M^\mu = 2$

$\mu$ is a binary variable which is equal to 1 when the $\bold{(a)}$ constraint is satisfied.

How to know if the two stacks are from the same truck?

We have $ST:stacks\times trucks$, a binary matrix.

Example with 4 stacks and 3 trucks

$\left [ \begin{matrix}0 & 1 & 0 & 0\\0 & 0& 1& 0\\0 & 0 & 0 & 1\\ 0&0&1&0\\0&0&0&1 \end{matrix} \right ]\times \left [ \begin{matrix}0 & 1 & 0\\0 & 0& 1\\1 & 0 & 0 \\ 0&0&1\end{matrix} \right ] = \left [ \begin{matrix}0 & 0 & 1\\1 & 0& 0\\0 & 0 & 1 \\ 1&0&0\\0&0&1\end{matrix} \right ]$

$\left [ \begin{matrix}1 & 0 & 0 & 0\\1 & 0& 0& 0\\1 & 0 & 0 & 0\\ 0&1&0&0\\0&1&0&0 \end{matrix} \right ]\times \left [ \begin{matrix}0 & 1 & 0\\0 & 0& 1\\1 & 0 & 0 \\ 0&0&1\end{matrix} \right ] = \left [ \begin{matrix}0 & 1 & 0\\0 & 1& 0\\0 & 1 & 0 \\ 0&0&1\\0&0&1\end{matrix} \right ]$



$\left [ \begin{matrix}0 & 0 & 1\\1 & 0& 0\\0 & 0 & 1 \\ 1&0&0\\0&0&1\end{matrix} \right ] - \left [ \begin{matrix}0 & 1 & 0\\0 & 1& 0\\0 & 1 & 0 \\ 0&0&1\\0&0&1\end{matrix} \right ] = \left [ \begin{matrix}0 & -1 & 1\\1 & -1& 0\\0 & -1 & 1 \\ 1&0&-1\\0&0&0\end{matrix} \right ]$


$\left [ \begin{matrix}0 & -1 & 1\\1 & -1& 0\\0 & -1 & 1 \\ 1&0&-1\\0&0&0\end{matrix} \right ]\times \left [ \begin{matrix} 1 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 3 \end{matrix} \right ] = \left [ \begin{matrix} 0 & -2 & 3 \\ 1 & -2 & 0 \\ 0 & -2  & 3 \\ 1 & 0 & -3 \\ 0 & 0 & 0 \end{matrix} \right ]$

$\left [ \begin{matrix} 0 & -2 & 3 \\ 1 & -2 & 0 \\ 0 & -2  & 3 \\ 1 & 0 & -3 \\ 0 & 0 & 0 \end{matrix} \right ]\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] + 3 = \left [ \begin{matrix} 4\\ 2\\4\\1\\3 \end{matrix} \right ]$

$\left ( \Xi^2 ST - \Xi^1 ST \right )\times \left [ \begin{matrix} 1 & 0 & 0 & \dots \\ 0 & 2 & 0& \dots \\ 0 & 0 & 3& \dots\\ etc \end{matrix} \right ] = \nu \quad \quad \bold{(b)}$

$sy^o_{s_2} \geq sy^e_{s_1}$ or $sy^e_{s_2} \leq sy^o_{s_1}$

$\Xi^1SY^e \leq \Xi^2SY^o +  \xi M^{TW} + f(\nu)M^{TW}$

$\Xi^2SY^e \leq \Xi^1SY^o + (1 - \xi)M^{TW} +  f(\nu)M^{TW}$

$\xi \in \{0, 1\}$

$M^{TW} = \max TW + 1$

$s_1$ and $s_2$ $\in$ same truck $\implies (\nu-|trucks| = 0)$ else $(\nu-|trucks|) \in \bold{R}$

<!-- $\tau - \phi \leq (\nu-|trucks|)M$  
$\tau \geq 0$  
$\phi \geq 0$

$(\nu-|trucks|) > 0 \implies \tau \text{ and } \phi \text{ free}$

$(\nu-|trucks|) = 0 \implies \tau - \phi = 0$

$(\nu-|trucks|) < 0 \implies \tau = 0\quad \land \quad \phi \geq  - (\nu-|trucks|)M$


$\tau \leq (\nu-|trucks|)M + \phi$


$\phi \leq (\nu-|trucks|)M + \tau$

$\displaystyle \tau \geq \frac{(\nu-|trucks|)M}{10}$

$(\nu-|trucks|) > 0 \implies \phi \geq 0,\; \tau \geq (\nu-|trucks|)M/10 > 0$

$(\nu-|trucks|) = 0 \implies \tau - \phi = \phi - \tau = 0$

$(\nu-|trucks|) < 0 \implies \phi = - (\nu-|trucks|)M\quad \land \quad \tau = - (\nu-|trucks|)M$ Doesn't work -->

$\tau - \phi \leq (\nu-|trucks|)M^\tau$  
$\displaystyle \tau \geq \frac{(\nu-|trucks|)M^\tau}{10}$  
$\tau \geq 0$  
$\phi \geq 0$  
<!-- $\tau - \phi \leq \eta^1 M$  
$\phi - \tau \leq (1-\eta^1)M$   -->
$\tau \leq \eta M^\eta$  
$\phi \leq (1-\eta)M^\eta$  
$\eta \in \{0, 1\}$

$M^\tau = 10$  
$\displaystyle M^\eta = |trucks|\frac{M^\tau}{10}$

$(\nu-|trucks|) > 0 \implies \phi = 0,\; \tau \geq (\nu-|trucks|)M^\tau/10 > 0$

$(\nu-|trucks|) = 0 \implies ((\tau = \phi)\; \land\; (\tau = 0\; \lor\; \phi = 0)) \implies \tau = \phi = 0$

$(\nu-|trucks|) < 0 \implies \tau = 0\quad \land \quad \phi \geq  - (\nu-|trucks|)M^\tau > 0$

$\Xi^1SY^e \leq \Xi^2SY^o +  \xi M^{TW} + (\tau + \phi)M^{TW} + (1-\mu)M^{TW}$

$\Xi^2SY^e \leq \Xi^1SY^o + (1 - \xi)M^{TW} +  (\tau + \phi)M^{TW} + (1-\mu)M^{TW}$

> That was a hard one!

Any stack must be adjacent to another stack on its left on the X axis, or if there is a single stack in the truck, the unique stack must be placed at the front of the truck (adjacent to the truck driver) :  
$\forall t \in \widetilde{T}, \forall s_1 \in \widetilde{TS_t}$ with $sx^o_{s_1} > 0, \exists s_2 \in \widetilde{TS_t}$ with $(sx^e_{s_2} = sx^o_{s_1})$ and $(sy^o_{s_2} \in [sy^o_{s_1}, sy^e_{s_1}]$ or $sy^e_{s_2} \in [sy^o_{s_1}, sy^e_{s_1}])$ **(P3)**.

> Ignoring it for the moment

Stacks must be placed in an increasing fashion from the front to the rear of the truck, (1) according to the suppliers' pickup order, and among the stacks of the same supplier, stacks must be placed in an increasing fashion (2) according to the supplier dock loading order, and among the stacks with the same supplier and supplier dock, stacks must be placed in an increasing fashion (3) according to the plant dock loading order. We do not try to optimize the unloading of the stacks, other than the satisfaction of the 3 precited orders. **(P4)**

$\forall t \in \widetilde{T}, \forall s_1 \in \widetilde{TS_t}, \forall s_2 \in \widetilde{TS_t}$,

(1) if $TE_{us_1} < TE_{us_2}$ then $sx^o_{s_1} \leq sx^o_{s_2}$

(2) if ($u_{s_1} = u_{s_2}$ and $KE_{k_{s_1}} < KE_{k_{s_2}}$) then $sx^o_{s_1} \leq sx^o_{s_2}$

(3) if ($u_{s_1} = u_{s_2}$ and $KE_{k_{s_1}} = KE_{k_{s_2}}$ and $GE_{g_{s_1}} < GE_{g_{s_2}}, \forall g_{s_1} \in \widetilde{SG_{s_1}}, \forall g_{s_2} \in \widetilde{SG_{s_2}}$) then $sx^o_{s_1} \leq sx^o_{s_2}$

(1)

We can set an order similarly as for **(P2)**:

$\left [ \begin{matrix} 1 & 0 & \dots & 0\\ & I & &  \end{matrix} \right ]\times SU \leq  SU\quad \text{(Wrong)}$

This implies that the stacks the most to the front of the cabin have their suppliers first in the supplier order.

$\begin{matrix} a1 & a2 \\ a3 & b4 \\ a5 & c6 \\ b7 & c8 \\ b9 & c10 \\ c11 & c12\end{matrix}$

It would work well for 1 truck but it will cause disruption in all trucks $\rightarrow$ both orders can't be satisfied simultaneously!

(1) can be re-written as

if $sx^o_{s_1} > sx^o_{s_2}$ then $TE_{us_1} \geq TE_{us_2}$

Well, we know every $s_1$ and $s_2$ which satisfy $sx^o_{s_1} > sx^o_{s_2}$.

Since we know well that

$\Xi^1SX^o \leq \Xi^2SX^o$

then we can add the following constraint (with needed modification that they both need to be in the same truck):

$\Xi^1SU\cdot TE \leq \Xi^2 SU \cdot TE + (\tau + \phi)M^{TE}$

$M^{TE} = \max TE$

> Is the strict inequality significant?  
> Can we re-use the same $\tau$ and $\phi$ here? Since we are comparing the same stacks, then I am tempted to say yes.


(2) Is similar

(2) can be re-written as

$sx^o_{s_1} > sx^o_{s_2} \implies (u_{s_1} \neq u_{s_2}\; \lor \; KE_{k_{s_1}} \geq KE_{k_{s_2}})$

(And $s_1$ and $s_2$ can be exchanged freely).

$\Xi^1SU \neq \Xi^2SU$ can be written as:

$\Xi^1SU - \Xi^2SU \geq \chi\epsilon - rM^{TE} - (\tau + \phi)M^{TE} - (1 - \sigma^1)M^{TE}$  
$\Xi^2SU - \Xi^1SU \geq (1 - \chi)\epsilon - rM^{TE} - (\tau + \phi)M^{TE} - (1 - \sigma^1)M^{TE}$


$\chi \in \{0, 1\}$ and where $\epsilon$ is a small positive constant to ensure strict inequality.

$r \in \{0, 1\}$

$\Xi^2SK \cdot TKE \geq \Xi^1SK \cdot TKE - (1-r)M^{TKE} - (\tau + \phi)M^{TKE}$

$M^{TKE} = \max TKE$

(3) can be re-written as:

$sx^o_{s_1} > sx^o_{s_2} \implies u_{s_1} \neq u_{s_2}\;\lor\;KE_{k_{s_1}} \neq KE_{k_{s_2}}\;\lor\; GE_{g_{s_1}} \geq GE_{g_{s_2}}, \forall g_{s_1} \in \widetilde{SG_{s_1}}, \forall g_{s_2} \in \widetilde{SG_{s_2}}$

The first term can be integrated in the previous constraint of (2) with adding $\sigma^1$.

Second term:

$\Xi^1SK \cdot TKE - \Xi^2SK \cdot TKE \geq \chi\epsilon - (\tau + \phi)M^{TKE} - (1 - \sigma^2)M^{TKE}$  
$\Xi^2SK \cdot TKE - \Xi^1SK \cdot TKE \geq (1 - \chi)\epsilon - (\tau + \phi)M^{TKE} - (1 - \sigma^2)M^{TKE}$

Since we simplified the problem by imposing that a stack could only have one plant dock, then the last term is simpler to formulate.

$\Xi^2SG\cdot TGE \geq \Xi^1SG\cdot TGE - (\tau + \phi)M^{TGE} - (1 - \sigma^3)M^{TGE}$

$M^{TGE} = \max TGE$

Plant dock loading order may be equal to 0 for several plant docks in some trucks, in which case, there is no constraint on the loading order for these plant docks in these trucks.
I understand this as these plant docks then have the greatest loading order automatically... Same for TKE : Supplier docks may be unknown for some trucks, in which case there is no constraint on the loading order of supplier docks in these trucks.

$\begin{matrix} 9 \\ 4 \\ 0 \\ 7 \\ 4\\ 3 \\ 2 \\ 0 \end{matrix}$

$\rightarrow\begin{matrix} 2 \\ 3 \\ 4 \\ 4 \\ 7\\ 9 \\ 0 \\ 0 \end{matrix}$



$\sigma^1 + \sigma^2 + \sigma^3 \geq 1$  
$\sigma^1, \sigma^2, \sigma^3 \in \{0, 1\}$

## Weight constraints

> Ignored

## How to take into account extra trucks?

For each planned truck $t$, add as many extra trucks that there are candidate items for truck $t$. Problem solved.

## Rewriting the objective function

$x_{{T}_1} \in \{0, 1\}^{|T|}$  
$x_{{T}_2} \in \bold{N}^{|T|}$  
$y \in \bold{N}^{N}$

$TI^T \in \{0, 1\}^{|T|\times N}$

$\displaystyle TI^E \in \{0, 1\}^{(\sum TRI)\times N}$

$TRI_t:$ Candidate items that can be picked up by truck $t$. 


$\displaystyle TI = \left [ \frac{TI^T}{TI^E} \right ]$

$$ \min\quad \alpha_T c^\top_{{T}_1} x_{{T}_1} + \alpha_E c^\top_{{T}_2} x_{{T}_2} + \alpha_Ic^\top_{I} (d - y) $$


$$ \min\quad \alpha_T c^\top_{{T}_1} \min(1,\,TI^T\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]) + \alpha_E c^\top_{{T}_2} \min(1,\,TI^E\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]) + \alpha_Ic^\top_{I} (IDL - TI^\top \times TDA) $$


$$ \min\quad \alpha_T c^\top_{{T}_1} \zeta^T + \alpha_E c^\top_{{T}_2} \zeta^E + \alpha_Ic^\top_{I} (IDL - TI^\top \times TDA) $$

$-\zeta^T \geq -\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$  

$-\zeta^T \geq -TI^T\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$


$-\zeta^E \geq -\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$  

$-\zeta^E \geq -TI^E\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$

## The full model

**Minimize:**

$$\alpha_T c^\top_{{T}_1} \zeta^T + \alpha_E c^\top_{{T}_2} \zeta^E + \alpha_Ic^\top_{I} (IDL - TI^\top \times TDA)$$

**Subject to:**


$$-\zeta^T \geq -\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$$  

$$-\zeta^T \geq -TI^T\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$$


$$-\zeta^E \geq -\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$$  

$$-\zeta^E \geq -TI^E\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$$


$$TI \leq F$$

$$TI^\top\times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] \leq \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right]$$

$$TI^\top \times TP = IP$$


$$R \leq \Theta M^{I4}$$  
$$ -M^{I4}(1-\Theta) \leq R - (TI^\top\times TU) \leq M^{I4}(1 - \Theta) $$  


$$R\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]=IU$$  

$$\Theta\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \geq \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$IDE \leq TI^\top \times TDA \leq IDL$$

$$Z\leq S\times M^{Z}$$  


$$\Psi = \left [ \Omega[\dots, 1, \dots]\cdot \left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ] \Omega[\dots, 2, \dots]\cdot \left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ]\dots  \right ]$$

$$\Omega[\dots, j, \dots]\leq S\cdot M^\Omega\quad \forall j$$

$$-(1-S)M^\Omega \leq \Omega[i,\dots,\dots] - TI^\top \leq (1-S)M^\Omega\quad \forall i$$

$$\Psi \leq ST \cdot M^\Psi$$

$$-(1-ST)M^\Psi \leq \Psi - \left [S\cdot\left [ \begin{matrix} 1\\\vdots\\1\end{matrix} \right ]\dots\right ] \leq (1-ST)M^\Psi$$



$$S\times IS = Z\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$-M^{Z}(1-S) \leq Z - \left [  SS  \cdots  SS  \right ] \leq M^{Z}(1-S)$$

$$Q\leq S\times M^Q$$  
$$S\times IU = Q\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$-M^Q(1-S) \leq Q - \left [  SU  \cdots  SU  \right ] \leq M^Q(1-S)$$

$$H\leq S\times M^H$$  
$$S\times IP = H\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$-M^H(1-S) \leq H - \left [  SP  \cdots  SP  \right ] \leq M^H(1-S)$$

$$V\leq S\times M^V$$  
$$S\times IK = V\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$-M^V(1-S) \leq V - \left [  SK  \cdots  SK  \right ] \leq M^V(1-S)$$

$$W\leq S\times M^{W}$$  
$$\displaystyle S\times IPD = W\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$-M^{W}(1-S) \leq W - \left [  SPD  \cdots  SPD  \right ] \leq M^{W}(1-S)$$

$$G^l\leq S\times M^G$$  
$$G^r \leq S\times M^G$$  
$$G^l\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = G^r\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$-M^G(1-S) \leq G^r - \left [  SO  \cdots  SO  \right ] \leq M^G(1-S)$$  
$$-M^G(1-S) \leq G^l - \left [  IOV  \cdots  IOV  \right ] \leq M^G(1-S)$$

$$SX^e - SX^o = SL + SO\cdot M^{TL}$$

$$SY^e - SY^o = SW + SO\cdot M^{TW}$$

$$SX^e - SX^o = SW + (1-SO)\cdot M^{TW}$$

$$SY^e - SY^o = SL + (1-SO)\cdot M^{TL}$$

$$SZ^o = 0$$

$$SZ^e = S\cdot IH$$

$$SX^e \leq  ST \times TL$$  
$$SY^e \leq  ST \times TW$$  
$$SZ^e \leq  ST \times TH$$

$$\left [ \begin{matrix} 1 & 0 & \dots & 0\\ & I & &  \end{matrix} \right ]\times SX^o \leq  SX^o$$

$$\Xi^2 SX^o - \Xi^1 SX^{e} - \beta^- + \beta^+ =  - 0.0001\quad \bold{(a)}$$

$$\beta^- \leq \lambda M^\lambda$$  
$$\beta^+ \leq (1-\lambda)M^\lambda$$  

$$(1 - \mu) \leq \beta^- M^\mu$$

$$\left ( \Xi^2 ST - \Xi^1 ST \right )\times \left [ \begin{matrix} 1 & 0 & 0 & \dots \\ 0 & 2 & 0& \dots \\ 0 & 0 & 3& \dots\\ etc \end{matrix} \right ] = \nu \quad \quad \bold{(b)}$$


$$\tau - \phi \leq (\nu-|trucks|)M^\tau$$  
$$\displaystyle \tau \geq \frac{(\nu-|trucks|)M^\tau}{10}$$  
$$\tau \geq 0$$  
$$\phi \geq 0$$  
<!-- $$\tau - \phi \leq \eta^1 M$$  
$$\phi - \tau \leq (1-\eta^1)M$$   -->
$$\tau \leq \eta M^\eta$$  
$$\phi \leq (1-\eta)M^\eta$$  

$$\Xi^1SY^e \leq \Xi^2SY^o +  \xi M^{TW} + (\tau + \phi)M^{TW} + (1-\mu)M^{TW}$$

$$\Xi^2SY^e \leq \Xi^1SY^o + (1 - \xi)M^{TW} +  (\tau + \phi)M^{TW} + (1-\mu)M^{TW}$$

$$\Xi^1SU\cdot TE \leq \Xi^2 SU \cdot TE + (\tau + \phi)M^{TE}$$

$$\Xi^1SU - \Xi^2SU \geq \chi\epsilon - rM^{TE} - (\tau + \phi)M^{TE} - (1 - \sigma^1)M^{TE}$$  
$$\Xi^2SU - \Xi^1SU \geq (1 - \chi)\epsilon - rM^{TE} - (\tau + \phi)M^{TE} - (1 - \sigma^1)M^{TE}$$

$$\Xi^2SK \cdot TKE \geq \Xi^1SK \cdot TKE - (1-r)M^{TKE} - (\tau + \phi)M^{TKE}$$

$$\Xi^1SK \cdot TKE - \Xi^2SK \cdot TKE \geq \chi\epsilon - (\tau + \phi)M^{TKE} - (1 - \sigma^2)M^{TKE}$$  
$$\Xi^2SK \cdot TKE - \Xi^1SK \cdot TKE \geq (1 - \chi)\epsilon - (\tau + \phi)M^{TKE} - (1 - \sigma^2)M^{TKE}$$

$$\Xi^2SG\cdot TGE \geq \Xi^1SG\cdot TGE - (\tau + \phi)M^{TGE} - (1 - \sigma^3)M^{TGE}$$

$$\sigma^1 + \sigma^2 + \sigma^3 \geq 1$$  


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

## Analysis

The decision variables which constrain the problem the most are $TI$, $S$ and $ST$. In particular, $ST$ is the only passage to the placement constraints. With $ST$ and $ZE^e$ fixed, we obtain two independant problems, multiplied by the number of trucks.

I think we could probably fix $ST$ because we will define as many stacks as there are candidate items for which truck, which means each stack is actually already tied to a truck.

However we can't really fix $ST$ because the truck of each stack we depend on their position because of the order we imposed.

With $TI$ fixed, the problem becomes much simpler to solve since it becomes only a constraint satisfaction problem.


## Resolution

We are going to use the Benders decomposition as main method, with $TI$ being the variables of the master problem.

We represent all the decision variables excluding $TI$ by $X$, the constraints of the problem as $K$. We represent the constraints coefficients concerning the variable $x$ as $K(x)$, and the cost associated with $x$ by $C(x)$.

The problem can be re-written as:

**Objective function:**

$$\min\quad C(X)\cdot X + C(TI)\cdot TI$$

**Subject to:**

$$K(X)\cdot X + K(TI)\cdot TI \geq K$$

$$X \in \cal{X}$$

$$TI \in \cal{TI}$$

Residual problem:

**Objective function:**

$$\min\quad C(X)\cdot X + C(TI)\cdot \overline{TI}$$

**Subject to:**

$$K(X)X \geq K - K(TI)\overline{TI}$$

We will not solve the dual problem and will avoid having to generate feasibility cuts. However, without the dual problem we can't generate optimality cuts either...

The $TI$ fixed, what are the constraints that tie the subproblems consisting in placing the stacks optimally in each truck? We can already notice that the objective function is entirely determined by $TI$, which means that the remaining subproblems are not optimization problems. The stacks are globally ordered. But it is not necessarily an issue. When solving the subproblems of each truck, we can relabel each stack as to satisfy the global order.

We can probably make use of the fact that the objective function depends on so few variables. The thing is, the objective value will be the same wether the placement constraints are satisfied or not, as long as the items are placed in each truck... but not for the dual objective function...

Other approach: we can use the dual of the linear relaxation of the problem. We will not solve the dual anyway, we will find an integer solution, so the optimality cuts we will be able to produce will be correct.

However, when solving the master problem, we will likely find a fractional solution for $TI$. Will rounding allow for convergence? There is also the question of finding a feasible solution. Typically, if the $TI$ found is unfeasible for some trucks, then it means that some items could not have been put into a stack.

We could soften the constraint stating that each item should be in a stack, but heavily penalize it in the objective function. Heavily means that a single missing item should cost more that what the worst solution to the original problem costs. I guess one of the worst solutions would be to use a single truck per item. We could easily compute the value of this 'solution', and use it as the cost per missing item. This way, there always is a feasible solution for each truck: leaving it empty of stacks. That being said, there are other constraints linked with $TI$ which might be violated.

```
General approach:

1. Preprocess the instance and formulate a linear relaxation of the original MILP.
2. Setup the master problem with the following objective function:
    min     z
    subject to
    cuts = {z >= 0}
3. Set TI to a random value, or the optimal value if there wasn't placement constraints.
4. Make sure (How?) that the value chosen for TI is feasible for all truck-wise sub-problems.
5. Solve separately (How?) each subproblem consisting in placing the items of a truck into stacks and placing those stacks into the truck.
6. Get the solution for X\{TI}
7. Compute the objective function of the dual of the linear relaxation of the original MILP. Add penalities for items with no stacks.
8. Generate an optimality cut and add it to `cuts`
9. Solve master problem and obtain new TI
10. Unless the optimal solution is found, return to 4.
```

## Solving the placement subproblem

Given a set of items and a truck, place each item in a stack, so that:
1. The stacks don't overlap
2. The stacks stay in the truck
3. The items of a stack have the same stackability code, supplier, supplier dock, plant dock, orientation
4. The loading orders are respected

```
Naive algorithm:

1. Group items per supplier and plant loading orders.
2. Compare each item with each other as to make stacks of similar items within each group.
3. Place stacks into truck based on supplier and plant loading order. This might need solving a linear program.
```

Grouping is done in $O(n)$ where $n$ is the number of items assigned to the truck. Comparison is done in $O(n!)$. Placing stacks is done in $O(k)$ where $k$ is the number of stacks and logically $k\leq n$ so $O(n)$. The naive algorithm has thus a complexity of $O(n!)$, which is pretty bad, because of the second step.

However we can assign items to stacks on the fly in $O(n)$ instead of comparing every item to each other. Which makes it an algorithm of $O(n)$.

Placing could still require to be done via optimization. We assume placing items into stacks can be done beforehand. The following program places stacks in the truck as to minimize overflow in the X axis:

**Minimize:**

$$\rho$$

**Subject to:**


$$G^l\leq S\times M^G$$  
$$G^r \leq S\times M^G$$  
$$G^l\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] = G^r\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$-M^G(1-S) \leq G^r - \left [  SO  \cdots  SO  \right ] \leq M^G(1-S)$$  
$$-M^G(1-S) \leq G^l - \left [  IOV  \cdots  IOV  \right ] \leq M^G(1-S)$$

$$SX^e - SX^o = SL + SO\cdot M^{TL}$$

$$SY^e - SY^o = SW + SO\cdot M^{TW}$$

$$SX^e - SX^o = SW + (1-SO)\cdot M^{TW}$$

$$SY^e - SY^o = SL + (1-SO)\cdot M^{TL}$$

$$SZ^o = 0$$

$$SZ^e = S\cdot IH$$

$$(SX^e \leq  TL)\quad (relaxed)$$  

$$\rho \geq SX^e - TL$$

$$SY^e \leq  TW$$  
$$SZ^e \leq  TH$$

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

$SZ^o = 0$

$M^\mu = 2$

$\Xi^1 : \left [ \begin{matrix}1 & 0 & 0 & 0&\dots\\1 & 0& 0& 0&\dots\\1 & 0 & 0 & 0&\dots\\ \vdots\\ 0&1&0&0&\dots\\0&1&0&0&\dots\\etc \end{matrix} \right ]$

$\Xi^2 : \left [ \begin{matrix}0 & 1 & 0 & 0&\dots\\0 & 0& 1& 0&\dots\\0 & 0 & 0 & 1&\dots\\ \vdots\\ 0&0&1&0&\dots\\0&0&0&1&\dots\\etc \end{matrix} \right ]$

$\epsilon = 0.001$

This is a Mixed Integer Linear Program. We could try to solve it with B&B, with a smart branching on the binary variables which are not determined by other variables. The use of heuristic algorithms which ignore orientation for instance could prove useful. We could even use Benders decomposition.

Once  the solution is found, the stacks overflowing from the truck are deleted and we deduce the items which could not be added to stacks that way.

Problem: $\Xi$ constraints are in an exponential number. ERRATUM: $\Xi$ constraints are actually as many as $\frac{|stacks|(|stacks| + 1)}{2}$, which means $O(n^2)$. Following next are thoughts on how to tackle the 'high' number of $\Xi$ constraints before I corrected myself.

### Using branch and cut to solve subproblems? 

How to detect a violated $\Xi$ constraint? Would it even be useful? In the subproblems, maybe some constraints can be logically satisfied. For instance, let's say we have $s_1, s_2, s_3, s_4$, four stacks which are ordered by $SX^o$. Let's consider the following $\bold{(\Xi_a)}$ constraint:

$$\Xi^2 SX^o - \Xi^1 SX^{e} - \beta^- + \beta^+ =  - 0.0001$$

There are two cases:  
Case 1: $SX^o_{s_j} \geq SX^e_{s_i}\quad \forall i<j$

In this case, $\beta^- = 0$, and thus $\bold{(\Xi_b)}$ and $\bold{(\Xi_c)}$ are unecessary because relaxed.

Case 2: $SX^o_{s_j} < SX^e_{s_i}\quad \forall i<j$ in which case the contrary happens.

Now, what can the relationship between $s_1$ and $s_2$ can teach us about the one between $s_1$ and $s_3$?

If $SX^o_{s_2} \geq SX^e_{s_1}$, it just means that $s_1$ and $s_2$ aren't overlapping on the X axis. Since we know that $SX^o_{s_3} \geq SX^o_{s_2}$, then it is logical that $SX^o_{s_3} \geq SX^e_{s_1}$. It means that the satisfaction of some constraints is implied by the satisfaction of some other constraints. A Branch & Cut algorithm would thus prove useful.

The same concept applies to $\bold{(\Xi_b)}$.

$\bold{(\Xi_d)}$ can be replaced by a transitive constraint.

How to find a violated constraint? It will depend on the type of each constraint.

$\bold{(\Xi_a)}$: if both $\beta^-$ and $\beta^+$ are positive, we know for sure that $\bold{(\Xi_a)}$ is violated.

$\left [ \begin{matrix} 1 & 0 & \dots & 0\\ & I & &  \end{matrix} \right ]\times SX^o -  SX^e - \beta^- + \beta^+ =  - 0.0001\quad \bold{\Xi_a^{trans}}$

Let's prove that if all $\bold{\Xi_a^{trans}}$ are satisfied then all $\bold{\Xi_a}$ are aswell.
Suppose we have:

$\forall j,i\quad j=i+1,\; SX^o_j - SX^e_i-\beta^-_{i,j}+\beta^+_{i,j}=-0.0001$

Let's prove that $\forall h > j$

$SX^o_h - SX^e_i-\beta^-_{i,h}+\beta^+_{i,h}=-0.0001$

$\forall j,i\quad j=i+1,\; SX^o_j - SX^e_i-\beta^-_{i,j}+\beta^+_{i,j}=-0.0001$

$\iff \forall i\quad SX^o_{i+1} - (SX^o_i+\alpha_i)-\beta^-_{i,i+1}+\beta^+_{i,i+1}=-0.0001$

Thus

$\sum_{k=i}^{h-1}(SX^o_{k+1} - SX^o_k-\alpha_k-\beta^-_{k,k+1}+\beta^+_{k,k+1}) = -0.0001(h-i)$

$\overset{rec.}{\iff}SX^o_{h} - SX^o_{i}-\alpha_i-\sum_{k=i+1}^{h-1}(\alpha_k)-\sum_{k=i}^{h-1}(\beta^-_{k,k+1})+\sum_{k=i}^{h-1}\beta^+_{k,k+1} = -0.0001(h-i)$

$\iff SX^o_{h} - SX^e_{i}-\sum_{k=i+1}^{h-1}(\alpha_k)-\sum_{k=i}^{h-1}(\beta^-_{k,k+1})+\sum_{k=i}^{h-1}\beta^+_{k,k+1} = -0.0001(h-i)$

$\iff SX^o_{h} - SX^e_{i}-\sum_{k=i+1}^{h-1}(\alpha_k)-\sum_{k=i}^{h-1}(\beta^-_{k,k+1})+\sum_{k=i}^{h-1}\beta^+_{k,k+1}+0.0001(h-i-1) = -0.0001$

If we take 

$\beta^-_{i,h} = \sum_{k=i}^{h-1}(\beta^-_{k,k+1})-0.0001(h-i-1)$ 

and 

$\beta^+_{i,h} = \sum_{k=i}^{h-1}\beta^+_{k,k+1} - \sum_{k=i+1}^{h-1}(\alpha_k)$

Then we have

$SX^o_h - SX^e_i-\beta^-_{i,h}+\beta^+_{i,h}=-0.0001\quad \square$

However, there is no guarantee that at least one between $\beta^+_{i,h}$ and $\beta^-_{i,h}$ is equal to zero. So now the question is: does-it matter? Let's suppose that these are both positive. We thus have $\mu$ which is free in $\{0, 1\}$. Which means that $\bold{\Xi_b}$ and $\bold{\Xi_c}$ can be either relaxed or strict. We saw earlier that if $SX^o_h \geq SX^e_i$ then $\bold{(\Xi_b)}$ and $\bold{(\Xi_c)}$ should be relaxed. If we replace  those two constraints with their transitive counterparts, then this means these are automatically relaxed since they won't appear for the couple $(s_i,s_h)$. Now, if $SX^o_h < SX^e_i$...

I have found something interesting: if a truck can only have $n$ stacks side to side (overlapping on the $X$ axis) then it means that every stack $s_i$ only have to satisfy placement constraints in regards to their $n$ predecessors on the $SX^o$ order. This finding will be of no help when the stacks will be extremely tiny in $Y$. But it can help determine separating algorithms, or even placing algorithms. It is a rectangle packing problem, which is NP-hard in the general case. But here we can't rotate the rectangles freely. These are either rotated length-wise of width-wise.

```
separate_i
A separating algorithm to determine if stack i violates constraints with its adjacent predecessors:

1. Consider s[i-1]. Determine the range of s[i-1] in the X axis.
2. Make a list L of s[i-1] and the predecessors of s[i-1] which overlap with s[s-1]. The list should be ordered by decreasing order of Xo position.
3. For each predecessor p in L:
    4. For every \Xi constraint between p and s[i]:
        5. Check if it is satisfied:
            6. If not satisfied, terminate and return it.
7. Leaving the loop, no unsatisfied constraint has been found, return None.
```

Complexity:
- Step 1. is done in $O(1)$.
- Step 2. is done in $n$, with $n$ being the total number of stacks.
- Step 3. : There are at most $\displaystyle \left \lfloor \frac{TW}{\underset{\{s_j\in S-\{s_i\}|j < i\}}{\min}(SL_s, SW_s)} \right \rfloor \leq n-1$ predecessors in L, which means that in worst case we do go through every $\Xi$ constraints which gives $O(n)$.

Although in the general case the complexity of this separating algorithm seems mediocre, we can show that in practice, it can be very efficient:

```
separate_Xi
A separating algorithm to determine if stack i violates constraints with its adjacent predecessors, for all i:

1. Get an list of the stacks S order by increasing Xo.
2. For each stack s:
    3. If separate_i(s) returns a constraint c:
        4. return c
    5. else:
        6. continue
7. No violated constraint found: return None.
```

Complexity of ``separate_Xi``:
- Step 1. is done in worst case in $O(n\cdot log(n))$
- The separation of $s_1$ is done in $O(1)$ because no predecessors.
- For $s_i$:
    * Making a list of overlapping predecessors only consists in checking if $s_{i-1}$ overlaps with $s_{i-2}$. If it does, the list of predecessors becomes $L_i = L_{i-1} + \{s_{i-1}\}$. It is done in $O(1)$.
    * There are at most $(i-1)$ predecessors.


### Discarding those constraints in the dual relaxation hoping it doesn't affect the approximation?

We can keep those constraints because these are as many as $O(n^2)$ in contrary to what I previously thought.

### Expressing the penality per no-stack item

Any empty column in $S$, the stack-item matrix, indicates an item with no stack. The goal is to count those. In order to do so, we need to sum $S$ over the lines, substract 1 (at this point, we only have -1's and 0's) and multiply by -1. Finally, sum the vector. We achieve this with the following formula:

$$\left [ \begin{matrix}1 & \dots & 1\end{matrix} \right ]\cdot\left (- S^\top \cdot \left [ \begin{matrix}1 \\ \vdots \\ 1\end{matrix} \right ] + 1\right )$$

The cost of using one single truck per item can easily be calculated in advance. The objective function is:

$$\alpha_T c^\top_{{T}_1} \zeta^T + \alpha_E c^\top_{{T}_2} \zeta^E + \alpha_Ic^\top_{I} (IDL - TI^\top \times TDA)$$

with

$$-\zeta^T \geq -\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$$  

$$-\zeta^T \geq -TI^T\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$$


$$-\zeta^E \geq -\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$$  

$$-\zeta^E \geq -TI^E\times\left [ \begin{matrix} 1\\\vdots\\1 \end{matrix}\right ]$$

$\zeta$ is only a mean of knowing which truck is used or not. In this case, we know which trucks are used. We have $n$ items, $m$ planned trucks. We assume $n \geq m$. The $m$ planned trucks are necessarily used, which leaves $n - m$ items to dispatch in extra trucks. In the worst case scenario, only the costliest extra trucks are used. In the worst case all items arrive late. The value of an upper bound of the objective function of the original problem is thus:

$$\alpha_T c^\top_T \left [\begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] + (n-m)\alpha^E\max c_E + \alpha_I c_I^\top(IDL - TDE)$$

## Valid $TI$

We could formulate a problem in which we want to find a feasible solution for $TI$ but we want to minimize the distance from the fractional $\widehat{TI}$.

**Minimize:**

$$\delta$$


**Subject to:**

$$\delta \geq d_{i, j} \quad \forall i,j$$

$$d \geq TI - \widehat{TI}$$

$$d \geq \widehat{TI} - TI$$

$$TI \leq F$$

$$TI^\top\times \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right ] \leq \left [ \begin{matrix}1\\\vdots\\1 \end{matrix} \right]$$

$$TI^\top \times TP = IP$$


$$R \leq \Theta M^{I4}$$  
$$ -M^{I4}(1-\Theta) \leq R - (TI^\top\times TU) \leq M^{I4}(1 - \Theta) $$  


$$R\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]=IU$$  

$$\Theta\times \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ] \geq \left [ \begin{matrix} 1\\\vdots\\1 \end{matrix} \right ]$$

$$IDE \leq TI^\top \times TDA \leq IDL$$

**Variables:**

$TI :$ truck $-$ item matrix

$TI \in \{0,1\}^{|trucks| \times N}$

$R \in \{0,1\}^{N \times |suppliers|}$

$\Theta \in \{0,1\}^{N \times |suppliers|}$

$d \in \bold{R^+}^{|trucks| \times N}$

$\delta \in \bold{R^+}$

**Parameters:**

$\displaystyle M^{I4}_{i,j} = \max_i TU_{i,j}\quad \forall j \in  \left \{1,\; Col(TI^\top)\right \}$

$IU \in \bold{N}^{|items|}$

$IP \in \bold{N}^{|items|}$

Good news is we have very 'few' variables (arguable), but these are all integers. $\widehat{TI}$ is the fractional solution found by solving the master problem of the Benders decomposition.

Actually, this might be easier than I thought. We can compute beforehand which items are compatible with which trucks, an $O(nm)$ operation (n: nb of items, m: nb of trucks). There are at most $n$ trucks, which makes it a $O(n^2)$ operation. We then filter $\widehat{TI}$ based on the result to get rid of impossible affectations in $O(n^2)$. Then, for each item, we take the truck with the maximum affectation value and replace with 1, the rest of the line becomes 0, in $O(n^2)$. In the worst case, there are some trucks with too many items. When that happens, the value of the objective function is penalized in the sub-problems and the optimality cut will guide towards a better solution.

## In Summary

```
Detailed general approach:

1. Setup a linear relaxation of the original MILP for the given instance
2. Compute the value of the worst case of the not relaxed original MILP.
3. Setup the master problem with the following objective function:
    min     z
    subject to
    cuts = {z >= 0}
4. Set TI to a random value.
5. Use an algorithm to determine an integer solution close to the one found.
6. Solve separately each subproblem consisting in placing the items of a truck into stacks and placing those stacks into the truck, first by forming stacks logically and then solving via B&B the placing problem.
7. Get the solution for X\{TI}
8. Compute the number of missing items.
9. Re-label the stacks so as to satisfy the global order in the X axis for the relaxation of the original MILP.
10. Compute the objective function of the dual of the linear relaxation of the original MILP. Add penalities for items with no stacks.
11. Generate an optimality cut and add it to `cuts`
12. Solve master problem and obtain new potentially fractional TI
13. Unless the optimal solution is found, return to 5.
```