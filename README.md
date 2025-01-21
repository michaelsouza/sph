# SPH

SPH is an app to solve Molecular Distance Geometry Problems (MDGP's). SPH is written in C++ and depends on libgsl-dev library.

## Creating instances
1. Download 1GPV protein geometry from PDB
2. Using some bioinfo library, extract the atoms of the chain A
3. Save the geometry to the "1GPV_A.txt" file with the following columns: ATOM_SERIAL_NUMBER, ATOM_NAME, RESIDUE_NAME, RESIDUE_SEQUENCE_IDENTIFIER, CHAIN_ID, X, Y, Z
4. Define $V$ as the set of the first $N$ atoms (sorted by (RESIDUE_SEQUENCE_IDENTIFIER, ATOM_SERIAL_NUMBER)) in the chain A and extract this set of atoms.
5. Define $K$ as the set of residues with at least one atom in $V$.
6. Define $E \subset V \times V$ as the set of pairs of atoms in the same or consecutive residues.
7. Define $l_{ij} = (1 - \epsilon)d_{ij}$ and $u_{ij} = (1 + \epsilon)d_{ij}, \forall (i,j) \in E$, where $d_{ij}$ is the distance between atom $i$ and atom $j$.
7. For each $\epsilon \in \{0.04, 0.08, 0.12, 0.16\}$ and $N \in \{100, 200\}$, extract the first $N$ atoms from the geometry file and create an instance file with columns "ATOM_SERIAL_NUMBER_I, ATOM_SERIAL_NUMBER_J, DISTANCE_IJ, LOWER_BOUND_IJ, UPPER_BOUND_IJ". Name the file as "instance_N_<N>_eps_<eps>.txt".

## Functions definitions
Let 
- $G=(V,E)$ be an undirected graph with $E \subset V \times V$ and $V=\{1, \dots, N\}$, 
- $l_{i,j}, u_{i,j} \in \mathbb{R}$ for each $(i,j) \in E$, $\tau > 0$ and $\lambda > 0$,
- $\phi_{\lambda, \tau}(y) = \lambda y + \sqrt{\lambda^2 y^2 + \tau^2},$ where $y \in \mathbb{R}.$
- $\theta_{\tau}(x) = \sqrt{\tau^2 + x_{1}^2 + x_{2}^2 + x_{3}^2},$ where $x \in \mathbb{R}^3.$
- $\theta^{i,j}_{\tau}(X) = \theta_{\tau}(x_i - x_j)$, where $X \in \mathbb{R}^{N \times 3}$ and $x_i, x_j \in \mathbb{R}^3$ are the $i$-th and $j$-th rows of $X$.
- $f_{\lambda, \tau}(X) = \sum_{(i, j)\in E} \phi_{\lambda, \tau}(l_{i,j} - \theta^{i,j}_{\tau}(X)) + \phi_{\lambda, \tau}(\theta^{i,j}_{\tau}(X) - u_{i,j})$

## Derivatives
Let's calculate the gradient of $f_{\lambda,\tau}(X)$ w.r.t. the $k$-th row $x_k$. 

First, note that $f_{\lambda, \tau}(X)$ is given by

$$
f_{\lambda, \tau}(X) 
\;=\; \sum_{(i,j)\in E} \phi_{\lambda, \tau}\bigl(l_{i,j} - \theta_{\tau}^{i,j}(X)\bigr) 
\;+\; \phi_{\lambda, \tau}\bigl(\theta_{\tau}^{i,j}(X) - u_{i,j}\bigr),
$$

where

$$
\theta_{\tau}^{i,j}(X) \;=\; \sqrt{\tau^2 + \|x_i - x_j\|^2},
\quad
\phi_{\lambda,\tau}(y) 
\;=\; \lambda\,y \;+\; \sqrt{\lambda^2\,y^2 + \tau^2}.
$$


### 1. Partial derivatives of the building blocks

1. **Derivative of $\theta_{\tau}^{i,j}(X)$**

   For $\theta_{\tau}^{i,j}(X) = \sqrt{\tau^2 + \|x_i - x_j\|^2}$, we have
   $$
   \nabla_{x_i}\,\theta_{\tau}^{i,j}(X)
   \;=\;
   \frac{x_i - x_j}{\sqrt{\tau^2 + \|x_i - x_j\|^2}}
   \;=\;
   \frac{x_i - x_j}{\theta_{\tau}^{i,j}(X)}.
   $$
   By symmetry,
   $$
   \nabla_{x_j}\,\theta_{\tau}^{i,j}(X)
   \;=\;
   \frac{x_j - x_i}{\theta_{\tau}^{i,j}(X)}.
   $$
   If $k\neq i$ and $k\neq j$, then $\frac{\partial}{\partial x_k}\,\theta_{\tau}^{i,j}(X) = 0.$

2. **Derivative of $\phi_{\lambda,\tau}(y)$ wrt $y$**

   $$
   \phi_{\lambda,\tau}(y) 
   \;=\; \lambda\,y + \sqrt{\lambda^2\,y^2 + \tau^2}.
   $$
   Its derivative w.r.t. $y$ is
   $$
   \phi'_{\lambda,\tau}(y)
   \;=\; 
   \frac{d}{dy} \bigl[\lambda\,y + \sqrt{\lambda^2\,y^2 + \tau^2}\bigr]
   \;=\; 
   \lambda 
   \;+\;
   \frac{\lambda^2\,y}{\sqrt{\lambda^2\,y^2 + \tau^2}}.
   $$

### 2. Putting it all together for $\nabla_{x_k} f_{\lambda, \tau}(X)$

Because $f_{\lambda,\tau}(X)$ is a sum over edges $(i,j)\in E$, only edges that include node $k$ contribute to $\nabla_{x_k} f_{\lambda,\tau}(X)$. For an edge $(k,j)\in E$, the contribution to $f$ has two terms:

1. $\phi_{\lambda,\tau}\bigl(l_{k,j} - \theta_{\tau}^{k,j}(X)\bigr)$
2. $\phi_{\lambda,\tau}\bigl(\theta_{\tau}^{k,j}(X) - u_{k,j}\bigr)$

### 2.1 Derivative of the first term

$$
\frac{\partial}{\partial x_k}
\,\phi_{\lambda,\tau}\!\bigl(l_{k,j} - \theta_{\tau}^{k,j}(X)\bigr)
\;=\;
\phi'_{\lambda,\tau}\bigl(l_{k,j} - \theta_{\tau}^{k,j}(X)\bigr)
\;\cdot\;
\frac{\partial}{\partial x_k}\bigl(l_{k,j} - \theta_{\tau}^{k,j}(X)\bigr).
$$

Because $l_{k,j}$ is constant w.r.t. $x_k$, we get

$$
\frac{\partial}{\partial x_k}\bigl(l_{k,j} - \theta_{\tau}^{k,j}(X)\bigr)
\;=\;
-\,\frac{\partial}{\partial x_k}\,\theta_{\tau}^{k,j}(X)
\;=\;
-\,
\frac{x_k - x_j}{\theta_{\tau}^{k,j}(X)}.
$$

Hence,

$$
\frac{\partial}{\partial x_k}
\,\phi_{\lambda,\tau}\!\bigl(l_{k,j} - \theta_{\tau}^{k,j}(X)\bigr)
\;=\;
-\;\phi'_{\lambda,\tau}\bigl(l_{k,j} - \theta_{\tau}^{k,j}(X)\bigr)
\;\frac{x_k - x_j}{\theta_{\tau}^{k,j}(X)}.
$$

### 2.2 Derivative of the second term

$$
\frac{\partial}{\partial x_k}
\,\phi_{\lambda,\tau}\!\bigl(\theta_{\tau}^{k,j}(X) - u_{k,j}\bigr)
\;=\;
\phi'_{\lambda,\tau}\bigl(\theta_{\tau}^{k,j}(X) - u_{k,j}\bigr)
\;\cdot\;
\frac{\partial}{\partial x_k}\bigl(\theta_{\tau}^{k,j}(X) - u_{k,j}\bigr).
$$

Since $u_{k,j}$ is constant w.r.t. $x_k$, 

$$
\frac{\partial}{\partial x_k}\bigl(\theta_{\tau}^{k,j}(X) - u_{k,j}\bigr)
\;=\;
\frac{\partial}{\partial x_k}\,\theta_{\tau}^{k,j}(X)
\;=\;
\frac{x_k - x_j}{\theta_{\tau}^{k,j}(X)}.
$$

Thus,

$$
\frac{\partial}{\partial x_k}
\,\phi_{\lambda,\tau}\!\bigl(\theta_{\tau}^{k,j}(X) - u_{k,j}\bigr)
\;=\;
\phi'_{\lambda,\tau}\bigl(\theta_{\tau}^{k,j}(X) - u_{k,j}\bigr)
\;\frac{x_k - x_j}{\theta_{\tau}^{k,j}(X)}.
$$

### 2.3 Summing over edges

Summing the contributions from these two terms for each edge $(k,j)$ that touches $k$:

$$
\frac{\partial}{\partial x_k} f_{\lambda,\tau}(X)
\;=\;
\sum_{j : (k,j)\in E}
\Bigl[
-\;\phi'_{\lambda,\tau}\bigl(l_{k,j} - \theta_{\tau}^{k,j}(X)\bigr)
\;+\;
\phi'_{\lambda,\tau}\bigl(\theta_{\tau}^{k,j}(X) - u_{k,j}\bigr)
\Bigr]
\;\frac{x_k - x_j}{\theta_{\tau}^{k,j}(X)}.
$$

Since $(k,j)$ and $(j,k)$ represent the same undirected edge, you just sum over neighbors $j$ of $k$.

### 3. Final formula

Putting it plainly, for each $k=1,\dots,N$,

$$
\nabla_{x_k} f_{\lambda,\tau}(X)
\;=\;
\sum_{j : (k,j)\in E}
\Bigl[
-\;\phi'_{\lambda,\tau}\bigl(l_{k,j} - \theta_{\tau}^{k,j}(X)\bigr)
\;+\;
\phi'_{\lambda,\tau}\bigl(\theta_{\tau}^{k,j}(X) - u_{k,j}\bigr)
\Bigr]
\;\frac{x_k - x_j}{\theta_{\tau}^{k,j}(X)}.
$$

where 
$$
\phi'_{\lambda,\tau}(y)
\;=\;
\lambda \;+\; \frac{\lambda^2\,y}{\sqrt{\lambda^2\,y^2 + \tau^2}},
\quad
\theta_{\tau}^{k,j}(X)
\;=\;
\sqrt{\tau^2 + \|x_k - x_j\|^2}.
$$

## References

Souza, Michael, et al. "Hyperbolic smoothing and penalty techniques applied to molecular structure determination." Operations Research Letters 39.6 (2011): 461-465.
