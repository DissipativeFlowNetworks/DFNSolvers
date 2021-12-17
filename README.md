# DFNSolvers
Companion codes to the paper [Delabays et al (2022)]. 
Solves the Dissipative Network Flow problem for general graph structures. 

## Summary of the files
- **acyclic\_algorithm.jl**: Loads the scripts for computing the DFN solution for acyclic networks. 
- **cyclic\_iterations.jl**: Loads the scripts for computing the DFN solutions for general, cyclic graphs. 
- **Documentation.pdf**: Detailed documentation of the code. 
- **examples\_rts96.jl**: Runs the code on two examples, an arbitrary spanning tree of the IEEE RTS-96 test case, and a modified version of the full IEEE RTS-96 test case network.
- **ntw\_data**: Contains the data files necessary for the examples.
- **toolbox.jl**: Loads a list of tools useful for the code and for more general purposes.

## List of functions
- [acyclic\_algorithm](#acyclic\_algorithm)
- [cohesiveness\_adj](#cohesiveness\_adj)
- [cohesiveness\_inc](#cohesiveness\_inc)
- [cycle\_proj](#cycle\_proj)
- [dcc](#dcc)
- [deindex](#deindex)
- [dichot](#dichot)
- [hs](#hs)
- [iterations](#iterations)
- [L2B](#L2B)
- [load\_ksakaguchi](#load\_ksakaguchi)
- [reindex](#reindex)
- [retro\_function](#retro\_function)
- [run\_acyclic\_algorithm](#run\_acyclic\_algorithm)
- [Sδ](#Sδ)
- [targets](#targets)
- [winding](#winding)

## Detailed documentation for the functions

### acyclic\_aglorithm
*./acyclic\_algorithm.jl*

- `acyclic\_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, H::Vector{Function}, ω::Vector{Float64}, hγ::Vector{Tuple{Float64,Float64}}, ϵ::Float64=1e-10)`
	
Recursive algorithm determining the existence of a unique solution for the Dissipative Flow Problem [Delabays, Jafarpour, and Bullo (2021)] on an acyclic graph. The algorithm requires the nodes and edges to be indexed as described in the original paper.

**INPUT**:\
`B`: Incidence matrix of the graph considered. The nodes and edges are assumed to ordered according to the requirements of the original reference [Delabays et al. (2022)]. \
`ω`: Vector of natural frequencies of the oscillators.\
`H`: Vector (of dimension 2m) of the transfer functions, relating the flows over the two orientations of the edges. For `e = 1:m`, `H[e](f) = h[e](h[e+m]^{-1}(f))` and for `e = m+1:2m`, `H[e](f) = h[e](h[e-m]^{-1}(f))`. \
`hγ`: Vector of tuples of the lower and upper bounds of the domain of the transfer functions.\
`ϵ`: Correction parameter to avoid evaluating the transfer functions exactly at the boundary of their domain.

**OUTPUT**:\
`exists`: Returns true if a solution exists. If it does, then it is unique [Delabays et al. (2022)].\
`ff`: Vector (of dimension 2m) of flows on the edges, solving the Dissipative Flow Problem. For `e = 1:m` (resp. `e = m+1:2m`), `ff[e]` is the flow over the edge defined by `B0[:,e]` (resp. `-B0[:,e]`).\
`φ`: Synchronous frequency corresponding the solution.

---

### cohesiveness\_adj
*./toolbox.jl*

- `cohesiveness_adj(θ::Vector{Float64}, A::Matrix{Float64})`

Computes the maximal (in absolute value) angular difference over the edges, usgin the adjacency (or Laplacian) matrix of the graph. 

**INPUT**:\
`θ`: Vector of angles. \
`A`: Adjacency or Laplacian matrix. 

**OUTPUT**:\
`dθ`: Maximal angle difference over the edges. 

--- 

### cohesiveness\_inc
*./toolbox.jl*

`cohesiveness_inc(θ::Vector{Float64}, B::Matrix{Float64})`

Computes the maximal (in absolute value) angular difference over the edges, using the incidence matrix of the graph `B`.

**INPUT**:\
`θ`: Vector of angles. \
`B`: Incidence matrix. 

**OUTPUT**:
`dθ`: Maximal angle difference over the edges. 

---

### cycle\_proj
*./toolbox.jl*

- `cycle_proj(B::Matrix{Float64}, w::Vector{Float64}=Float64[])`

Returns the cycle projection matrix `P`, possibly weighted by the weight vector `w`. See [Jafarpour et al., SIAM Review (2021)] for details.

**INPUT**:\
`B`: Incidence matrix of the (undirected) graph.\
`w`: Weight vector. If none is given, unity weights are assumed.

**OUTPUT**:\
`P`: Cycle projection matrix, defined in [Jafarpour et al., SIAM Review (2021)].

---

### dcc
*./toolbox.jl*

- `dcc(x::Union{Float64,Vector{Float64},Matrix{Float64}})`

Takes the modulo 2π of `x`, in the interval [-π,π). Is applied elementwise. 

**INPUT**:\
`x`: Value(s) whose modulo is to be taken. 

**OUTPUT**:\
`d`: Result of the modulo. 

---

### deindex
*./toolbox.jl*

- `deindex(id::Vector{Int64}, ed::Vector{Int64})`

Reverses the re-numbering operated by `reindex`.

**INPUT**:\
`id`: Vector of the original node indexing, with their new position in the output matrix of `reindex`. This is the second output of `reindex`.\
`ed`: Vector of the original edge indexing, with their new position in the output matrix of `reindex`. This is the third output of `reindex`.

**OUTPUT**:\
`di`: Inverse indexing of nodes to recover the original ordering of nodes, before applying `reindex`.\
`de`: Inverse indexing of edges to recover the original ordering of edges, berore applying `reindex`.

---

### dichot
*./toolbox.jl*

- `dichot(F::Function, l::Float64, u::Float64, tol::Float64=1e-6)`

Computes the zero of the function `F` in the interval [`l`,`u`] with tolerance `tol`, if it exits. The function `F` is assumed strictly monotone and well-defined on [`l`,`u`]. The algorithm proceeds by dichotomy.

**INPUT**:\
`F`: Function whose zero needs to be found. It is assumed strictly monotone.\
`l`: Lower bound of the interval where the zero needs to be found.\
`u`: Upper bound of the interval where the zero needs to be found. \
`tol`: Tolerance on the error to the actual solution.

**OUTPUT**:\
`z`: Zero of the function `F` over [`l`,`u`], if it exists. 

---

### hs
*./cyclic\_iterations.jl*

- `hs(x::Float64, h::Function, γ::Tuple{Float64,Float64}, s::Float64=1.)`

Extended coupling function, extending `h` to the whole real axis. Matches `h` on [`γ[1]`,`γ[2]`], is continuous, and has slope `s` outside of [`γ[1]`,`γ[2]`]. 

**INPUT**:\
`x`: Argument of the extended coupling function. \
`h`: Coupling function to be extended. \
`γ`: Tuple of the two bounds of the domain of `h`. \
`s`: Slope of the extended coupling function outside of the domain of `h`. \

**OUTPUT**:\
`hx`: h_s(x).

**Variations**:

- `hs(x::Vector{Float64}, h::Vector{Function}, γ::Vector{Tuple{Float64,Float64}}, s::Union{Float64,Vector{Float64}})`

Computes the value of the extended coupling function elementwise on `x`, `h`, `γ`, and `s`. 

- `hs(x::Vector{Float64}, h::Function, γ::Tuple{Float64,Float64}, s::Float64)`

Computes the value of the extended coupling function elementwise on `x`.

---

### iterations
*./cyclic\_iterations.jl*

- `iterations(Δ0::Vector{Float64}, B::Matrix{Float64}, C::Matrix{Float64}, u::Vector{Int64}, ω::Vector{Float64}, h::Union{Function,Vector{Function}}, γ::Union{Tuple{Float,Float64},Vector{Tuple{Float64,Float64}}}, δ::Float64, s::Union{Float64,Vector{Float64}}=1., max_iter::Int64=100, tol::Float64=1e-6, verb::Bool=false)`

Runs the iteration scheme described in [Delabays et al. (2022)] for cyclic networks of diffusively coupled oscillators. Starts at initial conditions `Δ0`, which is projected on the afine subspace ker(B)^perp + 2π*pinv(`C`)*`u`. The script runs for at most `max_iter` iterations or until the correction is smaller than `tol`.

**INPUT**:\
`Δ0`: Initial conditions of the iterations. Each component should ideally be bounded by the corresponding components of `γ`. \
`B`: Incidence matrix of the (undirected) graph. \
`C`: Cycle-edge incidence matrix associated to the cycle basis of the graph (see [Delabays et al. (2022)]. \
`u`: Winding vector of the cell where the solution is searched. \
`ω`: Vector of natural frequencies. \
`h`: Vector of coupling functions over the (bidirected) edges of the graph. If a single function is given, the couplings are assumed homogenous. \
`γ`: Vector of tuples, composed of the lower (1st comp.) and upper (2nd comp.) bounds on the domain of `h`, such that it is strictly increasing. \
`δ`: Scaling parameter. \
`s`: Slope of the extended coupling functions. Should have the same dimension as `h`. \
`max_iter`: Maximal number of iterations allowed. \
`tol`: Minimal correction allowed between two iterations. \
`verb`: If true, enumerates the iterations.

**OUTPUT**:\
`Δ`: Final state at the end of the iterations. \
`Δs`: Sequence of states along the iterations. 

---

### L2B
*./toolbox.jl*

- `L2B(L::Union{Matrix{Float64},Matrix{ComplexF64},SparseMatrixCSC{Float64,Int64},SparseMatrixCSC{ComplexF64,Int64}})`

Compute an incidence matrix associated to a Laplacian matrix.

**INPUT**:\
`L`: Laplacian matrix (dense or sparse).

**OUTPUT**:\
`B`: Incidence matrix of the graph (same sparsity as the input).\
`w`: Vector of edge weights.\
`Bt`: Transpose of B.

--- 

### load\_ksakaguchi
*./toolbox.jl*

- `load_ksagaguchi(as::Vector{Float64},ϕs::Vector{Float64})`

Loads the Kuramoto-Sakaguchi coupling functions, ready to use in the acyclic_algorithm script.

**INPUT 1**:\
`as`: Vector of coupling weights, ordered according to the edge indexing in the incidence matrix. The vector `as` is of dimension 2m, as the two orientations of each edge are distiguished. If the edge e=(i,j), then the index of edge (j,i) is e+m. \
`ϕs`: Vector of phase frustrations. Ordering follows the same rules as `as`.

**OUTPUT 1**:\
`h`: Vector of the coupling functions, following the ordering of `as` and `ϕs`.\
`hi`: Vector of the inverse of the coupling functions, following the ordering of `as` and `ϕs`.\
`γ`: Vector of tuples composed of the lower (1st component) and upper (2nd component) bounds on the coupling functions, so that they are strictly increasing on (γ[1],γ[2]).

- `load_ksakaguchi(Y::Matrix{ComplexF64})`

Loads the coupling functions, their inverse, their domain bounds, and the incidence matrix of the network of Kuramoto oscillators, based on the associated admittance matrix. 

**INPUT 2**:\
`Y`: Admittance matrix of the system under consideration.

**OUTPUT 2**:\
`h`: Vector of the coupling functions, following the ordering of edges in the incidence matrix. Amplitudes are sqrt{B^2 + G^2} and phase frustrations are arctan(G/B).\
`hi`: Vector of the inverse of the coupling functions.\
`γ`: Vector of tuples composed of the lower (1st component) and upper (2nd component) bounds on the coupling functions, so that they are strictly increasing on (γ[1],γ[2]).\
`B`: Incidence matrix of the underlying network (undirected). \
`as`: List of coupling weights. \
`ϕs`: List of coupling frustrations.

---

### reindex
*./toolbox.jl*

- `reindex(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, id::Vector{Int64}=Int64[], ed::Vector{Int64}=Int64[])`

Re-orders the nodes and edges of the graph to satisfy the conditions for Theorem XXX [Delabays, Jafarpour, and Bullo (2021)] to apply. 

**INPUT**:\
`B`: Incidence matrix of the _bidirected_ graph. \
`id`: Indices of the nodes. If none is given, the indexing is assumed to be sequential, i.e., id = 1:n.\
`ed`: Indices of the edges. If none is given, the indexing is assumed to be sequential, i.e., ed = 1:n-1.

**OUTPUT**:\
`B2`: Re-ordered incidence matrix of the _bidirected_ graph. The graph structure is preserved, only the numbering is modified. \
`id2`: Vector of the original node indexing, with their new position in matrix B2.\
`ed2`: Vector of the original edge indexing, with their new position in matrix B2.

---

### retro\_function
*./toolbox.jl*

- `retro_function(i::Int64, m::Int64, ω::Vector{Float64}, H::Vector{Function}, et::Dict{Int64,Vector{Int64}})`

Recursively defines the flow over an edge as a function of the synchronous frequency `φ`.

**INPUT**:\
`i`: Index of the edge whose flow function has to be defined. Note that, according to indexing, `i` is also the index of the node at the source of edge `i`.\
`m`: Number of edges. \
`ω`: Vector of natural frequencies of the oscillators. \
`H`: Vector of the transfer functions over all (directed) edges. \
`et`: Dictionary associating to each node i, the list of nodes with lower index, connected to it. This is the second output of `targets`.

**OUTPUT**:\
`F`: Flow function over edge `i`, with respect to the synchronous frequency `φ`.

--- 

### run\_acyclic\_algorithm
*./acyclic\_algorithm.jl*

- `run_acyclic_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, ω::Vector{Float64}, h::Function, hi::Function, γ::Tuple{Float64,Float64})`
- `run_acyclic_algorithm(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, ω::Vector{Float64}, h::Vector{Function}, hi::Vector{Function}, γ::Vector{Tuple{Float64,Float64}})`

Runs the algorithm to decide if a solution exist for the Dissipative Flow Problem [Delabays et al. (2022)] on an acyclic graph. The nodes and edges of the graph are first reordered in order to match requirements described in the paper. Then the original numbering is retrieved.

**INPUT**:\
`B`: Incidence matrix of the graph considered, with n vertices and m edges. The graph is considered undirected and unweighted (weights can be incorporated in the coupling functions).\
`ω`: Vector of natural frequencies of the oscillators. \
`h`: Coupling function(s) over the edges of the interaction graph. If a single function is given, the couplings are assumed to identical. Otherwise, a (2m)-vector of functions needs to be given. Edge indexing is such that for e = 1:m, the orientation of e corresponds to the one given by the e-th column of B, and for e = m+1:2m, the orientation of e is the opposite of edge e-m.\
`hi`: Inverse(s) of the coupling function(s). \
`γ`: Tuple(s) of the lower and upper bounds on the argument of the coupling function(s) `h`, such that they are strictly increasing. Dimension is the same as `h`.

**OUTPUT**:\
`exists`: Returns true only if a solution exists. If the solution exists, then it is unique [Delabays, Jafarpour, and Bullo (2021)].\
`θ`: Solution (if it exists) of the Dissipative Flow Network problem, defined up to a constant angle shift. \
`ff`: Vector (of dimension 2m) of flows on the edges. For e = 1:m (resp. e = m+1:2m), ff[e] is the flow over the edge defined by B0[:,e] (resp. -B0[:,e]).\
`φ`: Synchronous frequency corresponding to the solution. 

---

### Sδ
*./cyclic\_iterations.jl*

- `Sδ(Δ::Vector{Float64}, ω::Vector{Float64}, B::Matrix{Float64}, Bout::Matrix{Float64}, P::Matrix{Float64}, W::Matrix{Float64}, δ::Float64, h::Union{Function,Vector{Function}}, γ::Union{Tuple{Float64,Float64},Vector{Tuple{Float64,Float64}}}, s::Union{Float64,Vector{Float64}}=1.)`

Iteration functions whose fixed points are solutions to the Dissipative Network Flow problem [Delabays et al. (2022)]. 

**INPUT**:\
`Δ`: Argument of the interation function. \
`ω`: Vector of natural frequencies. \
`B`: Incidence matrix of the (undirected) graph. \
`Bout`: Out-incidence matrix of the bidirected graph. \
`P`: Cycle projection matrix (see [Delabays et al. (2022)]. \
`W`: Weight matrix to be tuned. In [Delabays et al. (2022)], we take the pseudoinverse of the Laplacian. \
`δ`: Scaling parameter which, if small enough, guarantees `Sδ` to be contracting. \
`h`: Vector of the (directed) coupling functions. If a single `Function` is given, the couplings are assumed homogeneous. \
`γ`: Vector of tuples composed of the lower (1st comp.) and upper (2nd comp.) bounds on the domain of `h`, such that it is stricly increasing. Should have the same dimension as `h`. \
`s`: Slope of the extended coupling functions, outside of their domain. 

**OUTPUT**:\
`Δ2`: Updated value of the state `Δ`.

---

### targets
*./toolbox.jl*

- `targets(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}})`

Given an incidence matrix for an acyclic graph, with indexing satisfying the assumptions of Theorem XXX in [Delabays, Jafarpour, and Bullo (2021)], returns, for each node i: (1) The unique node j > i to which it is connected; (2) The list of nodes j < i that are connected to it.

**INPUT**:\
`B`: Incidence matrix of the graph. The graph is assumed connected and acyclic, and its numbering is assumed to satisfy the assumptions of Theorem XXX in [Delabays, Jafarpour, and Bullo (2021)].

**OUTPUT**:\
`te`: Dictionary associating to each node i = 1:n-1, the unique node with larger index to which it is connected.\
`et`: Dictionary associating to each node i = 1:n, the list of nodes with smaller index to which it is connected.

--- 

### winding
*./toolbox.jl*

- `winding(θ::Vector{Float64}, Σ::Vector{Int64})`

Computes the winding number associated to the state `θ` around the cycle `σ`. 

**INPUT**:\
`θ`: Vector of angles.  \
`σ`: Sequence of node indices forming a cycle. 

**OUTPUT**:\
`q`: Winding number.

- `winding(θ::Vector{Float64}, Σ::Vector{Vector{Int64}})`

Applies `winding` on each component of `Σ`.

---


