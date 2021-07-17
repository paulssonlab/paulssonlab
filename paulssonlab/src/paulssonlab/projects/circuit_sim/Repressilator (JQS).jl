### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 54918c32-e4cf-11eb-04d2-35c34dcb1016
begin
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add([
		Pkg.PackageSpec(name="Plots", version="1"),
        Pkg.PackageSpec(name="ModelingToolkit", version="5"),
		Pkg.PackageSpec(name="Catalyst", version="6"),
        Pkg.PackageSpec(name="OrdinaryDiffEq", version="5"),
        Pkg.PackageSpec(name="Latexify", version="0.15"),
        Pkg.PackageSpec(name="DiffEqJump", version="6"),
		Pkg.PackageSpec(url="https://github.com/augustinas1/MomentClosure.jl"),
		Pkg.PackageSpec(name="Turing", version="0.16"),
        Pkg.PackageSpec(name="Distributions", version="0.25"),
        Pkg.PackageSpec(name="DifferentialEquations", version="6"),
        Pkg.PackageSpec(name="MCMCChains", version="4"),
        Pkg.PackageSpec(name="StatsPlots", version="0.14"),
    ])
    using Plots, ModelingToolkit, Catalyst, LinearAlgebra, OrdinaryDiffEq, Latexify, DiffEqJump, MomentClosure, Turing, Distributions, DifferentialEquations, MCMCChains, StatsPlots
end

# ╔═╡ 26c8214c-13e6-4520-a2bc-13da053961cc
begin
	using Random
	Random.seed!(14);
end

# ╔═╡ 41af0c95-a13c-466f-8c50-2ca7c123b96e
md"""# Repressilator."""

# ╔═╡ 3f200780-bcdf-4be5-86e5-5162498ea382
repressilator = @reaction_network begin
    hillr(P₃,α,K,n), ∅ --> m₁
    hillr(P₁,α,K,n), ∅ --> m₂
    hillr(P₂,α,K,n), ∅ --> m₃
    (δ,γ), m₁ ↔ ∅
    (δ,γ), m₂ ↔ ∅
    (δ,γ), m₃ ↔ ∅
    β, m₁ --> m₁ + P₁
    β, m₂ --> m₂ + P₂
    β, m₃ --> m₃ + P₃
    μ, P₁ --> ∅
    μ, P₂ --> ∅
    μ, P₃ --> ∅
end α K n δ γ β μ;

# ╔═╡ 0a6ada38-9477-4c39-81f4-ab219ba87bc3
latexify(repressilator)

# ╔═╡ 571ba43d-c07b-40e9-8c88-df5368b7da41
odesys = convert(ODESystem, repressilator)

# ╔═╡ 90d627a9-ac70-46b7-a3f2-87ee81bea86d
begin
	# parameters [α,K,n,δ,γ,β,μ]
	p = (.5, 40, 2, log(2)/120, 5e-3, 20*log(2)/120, log(2)/60)

	# initial condition [m₁,m₂,m₃,P₁,P₂,P₃]
	u₀ = [0.,0.,0.,20.,0.,0.]

	# time interval to solve on
	tspan = (0., 10000.)

	# create the ODEProblem we want to solve
	oprob = ODEProblem(repressilator, u₀, tspan, p)

	sol = solve(oprob, Tsit5())
end

# ╔═╡ 440ffeb0-3d72-4790-a867-fb479710bcaa
plot(sol.t, sol')

# ╔═╡ a389e105-3bb7-4527-ad44-c0c7465c42c0
md"""# Compositional repressilator."""

# ╔═╡ e565eefb-115b-41dc-ac5b-5b72e8ab993d
begin
	@parameters t α₀ α K n δ β μ
	@variables m(t) P(t) R(t)
	rxs = [
	       Reaction(α₀, nothing, [m]),
	       Reaction(α / (1 + (R/K)^n), nothing, [m]),
	       Reaction(δ, [m], nothing),
	       Reaction(β, [m], [m,P]),
	       Reaction(μ, [P], nothing)
	      ]

	specs = [m,P,R]
	pars  = [α₀,α,K,n,δ,β,μ]
	@named rs = ReactionSystem(rxs, t, specs, pars)

	# using ODESystem components
	@named os₁ = convert(ODESystem, rs; include_zero_odes=false)
	@named os₂ = convert(ODESystem, rs; include_zero_odes=false)
	@named os₃ = convert(ODESystem, rs; include_zero_odes=false)
	connections = [os₁.R ~ os₃.P,
	               os₂.R ~ os₁.P,
	               os₃.R ~ os₂.P]
	@named connected = ODESystem(connections, t, [], [], systems=[os₁,os₂,os₃])
	repressilator_c = structural_simplify(connected)
end

# ╔═╡ f15ee055-440a-447c-8cbd-113c1ee75374
begin
	pvals_c = [os₁.α₀ => 5e-4,
	         os₁.α => .5,
	         os₁.K => 40.0,
	         os₁.n => 2,
	         os₁.δ => (log(2)/120),
	         os₁.β => (20*log(2)/120),
	         os₁.μ => (log(2)/600),
	         os₂.α₀ => 5e-4,
	         os₂.α => .5,
	         os₂.K => 40.0,
	         os₂.n => 2,
	         os₂.δ => (log(2)/120),
	         os₂.β => (20*log(2)/120),
	         os₂.μ => (log(2)/600),
	         os₃.α₀ => 5e-4,
	         os₃.α => .5,
	         os₃.K => 40.0,
	         os₃.n => 2,
	         os₃.δ => (log(2)/120),
	         os₃.β => (20*log(2)/120),
	         os₃.μ => (log(2)/600)]
	u₀_c    = [os₁.m => 0.0, os₁.P => 20.0, os₂.m => 0.0, os₂.P => 0.0, os₃.m => 0.0, os₃.P => 0.0]
	tspan_c = (0.0, 100000.0)
	oprob_c = ODEProblem(repressilator_c, u₀_c, tspan_c, pvals_c)
	sol_c = solve(oprob_c, Tsit5())
end

# ╔═╡ f3cef13d-366b-4218-aa1d-3f754fa1dd63
plot(sol_c.t, sol_c')

# ╔═╡ Cell order:
# ╠═54918c32-e4cf-11eb-04d2-35c34dcb1016
# ╠═26c8214c-13e6-4520-a2bc-13da053961cc
# ╟─41af0c95-a13c-466f-8c50-2ca7c123b96e
# ╠═3f200780-bcdf-4be5-86e5-5162498ea382
# ╠═0a6ada38-9477-4c39-81f4-ab219ba87bc3
# ╠═571ba43d-c07b-40e9-8c88-df5368b7da41
# ╠═90d627a9-ac70-46b7-a3f2-87ee81bea86d
# ╠═440ffeb0-3d72-4790-a867-fb479710bcaa
# ╟─a389e105-3bb7-4527-ad44-c0c7465c42c0
# ╠═e565eefb-115b-41dc-ac5b-5b72e8ab993d
# ╠═f15ee055-440a-447c-8cbd-113c1ee75374
# ╠═f3cef13d-366b-4218-aa1d-3f754fa1dd63
