using SparseArrays, DelimitedFiles

include("acyclic_algorithm.jl")
include("cyclic_iterations.jl")
include("toolbox.jl")


@info "----------------------------------------------------------------"
@info "================================================================"
@info "Running the acyclic algorithm on a spanning tree of RTS-96."
@info "----------------------------------------------------------------"

Gsp = readdlm("ntw_data/myRTS96_spantree_G.csv",',')
G = sparse(Int64.(Gsp[:,1]),Int64.(Gsp[:,2]),Gsp[:,3])
Bsp = readdlm("ntw_data/myRTS96_spantree_B.csv",',')
B = sparse(Int64.(Bsp[:,1]),Int64.(Bsp[:,2]),Bsp[:,3])
Yst = Matrix(G + im*B)

id_st = Int64.(vec(readdlm("ntw_data/myRTS96_ids_spantree.csv",',')))

ω = vec(readdlm("ntw_data/myRTS96_om.csv",','))

h,hi,γ,Bst,as,ϕs = load_ksakaguchi(Yst)
Bstd = pinv(Bst)

exists,θ,ff,φ = run_acyclic_algorithm(Bst,.8*ω,h,hi,γ)

if exists
	dθ = cohesiveness_inc(θ,Bst)
	@info "The unique solution has been found. With maximal angular difference Δ = $(dθ)."
else
	@info "No solution as been found."
end



@info "================================================================"
@info "Running the cyclic iterations on RTS-96, in three winding cells (u1, u2, and u3)."
@info "----------------------------------------------------------------"

Gsp = readdlm("ntw_data/myRTS96_G.csv",',')
G = sparse(Int64.(Gsp[:,1]),Int64.(Gsp[:,2]),Gsp[:,3])
Bsp = readdlm("ntw_data/myRTS96_B.csv",',')
B = sparse(Int64.(Bsp[:,1]),Int64.(Bsp[:,2]),Bsp[:,3])
Y = Matrix(G + im*B)

include("ntw_data/myRTS96_cycles.jl")

u1 = Int64.(vec(readdlm("ntw_data/myRTS96_u1.csv",',')))
u2 = Int64.(vec(readdlm("ntw_data/myRTS96_u2.csv",',')))
u3 = Int64.(vec(readdlm("ntw_data/myRTS96_u3.csv",',')))

ω = vec(readdlm("ntw_data/myRTS96_om.csv",','))

h,hi,γ,B,as,ϕs = load_ksakaguchi(Y)

Bb = [B -B]
Bout = Bb.*(Bb .> 1e-2)
n,m = size(B)

δ = .01
s = 1.
max_iter = 1000
tol = 1e-5

Δ0 = [((γ[i][2]-γ[i][1])*rand() + γ[i][1]) for i in 1:m]


@info "First case:"
Δ1,Δ1s = iterations(Δ0,B,C,u1,ω,h,γ,δ,s,max_iter,tol)
θ1 = Bstd'*Δ1[id_st]
q1 = winding(θ1,Σ)
f1 = [[h[i](Δ1[i]) for i in 1:m]; [h[i+m](-Δ1[i]) for i in 1:m]]
res1 = ω - Bout*f1

if q1 == u1
	@info "Winding vector matches."
else
	@info "Winding vector does not match."
end
if maximum(res1) - minimum(res1) > 1e-2
	@info "Did not find a synchronous solution."
else
	@info "Found a synchronous solution, stored in θ1."
end
@info "----------------------------------------------------------------"


@info "Second case:"
Δ2,Δ2s = iterations(Δ0,B,C,u2,ω,h,γ,δ,s,max_iter,tol)
θ2 = Bstd'*Δ2[id_st]
q2 = winding(θ2,Σ)
f2 = [[h[i](Δ2[i]) for i in 1:m]; [h[i+m](-Δ2[i]) for i in 1:m]]
res2 = ω - Bout*f2

if q2 == u2
	@info "Winding vector matches."
else
	@info "Winding vector does not match."
end
if maximum(res2) - minimum(res2) > 1e-2
	@info "Did not find a synchronous solution."
else
	@info "Found a synchronous solution, stored in θ2."
end
@info "----------------------------------------------------------------"

@info "Third case:"
Δ3,Δ3s = iterations(Δ0,B,C,u3,ω,h,γ,δ,s,max_iter,tol)
θ3 = Bstd'*Δ3[id_st]
q3 = winding(θ3,Σ)
f3 = [[h[i](Δ3[i]) for i in 1:m]; [h[i+m](-Δ3[i]) for i in 1:m]]
res3 = ω - Bout*f3

if q3 == u3
	@info "Winding vector matches."
else
	@info "Winding vector does not match."
end
if maximum(res3) - minimum(res3) > 1e-2
	@info "Did not find a synchronous solution."
else
	@info "Found a synchronous solution, stored in θ3."
end
@info "================================================================"




