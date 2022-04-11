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

# ╔═╡ e565eefb-115b-41dc-ac5b-5b72e8ab993d
rn = @reaction_network begin
    α, S + I --> 2I
    β, I --> R
end α β

# ╔═╡ acc1c643-fe73-4072-9896-641316150ad4
begin
	p     = [.1/1000, .01]           # [α,β]
	tspan = (0.0,250.0)
	u0    = [999.0,1.0,0.0]          # [S,I,R] at t=0
	Δt 	  = 1
	prob1 = ODEProblem(rn, u0, tspan, p)
	sol   = solve(prob1, Tsit5())       # use Tsit5 ODE solver
end

# ╔═╡ 75be9104-8023-41b5-8b53-8f37929765ec
plot(sol, lw=2)

# ╔═╡ f95497db-7bcb-479b-9ef7-ed20ad87ddd9
begin
	σ = 50
	skip = 20
	sol1 = solve(prob1, Tsit5(), saveat=Δt)
	data_t = sol1.t[1:skip:end]
	data_nonoise = Array(sol1)[:,1:skip:end]
	data = data_nonoise + σ * randn(size(data_nonoise))
	plot(sol1, alpha = 0.3, legend = false); scatter!(data_t, data')
end

# ╔═╡ 9376ab45-046b-4bb7-b5b0-1a64afe68181
begin
	Turing.setadbackend(:forwarddiff)

	@model function fitlv(data, prob1)
	    σ ~ InverseGamma(3, 200)
	    α ~ truncated(Normal(p[1], p[1]), 0, p[1]*3)
	    β ~ truncated(Normal(p[2], p[2]), 0, p[2]*3)

	    params = [α,β]
	    prob = remake(prob1, p=params)
	    predicted = solve(prob, Tsit5(), saveat=data_t)

	    for i = 1:length(predicted)
	        data[:,i] ~ MvNormal(predicted[i], σ)
	    end
	end

	model = fitlv(data, prob1)

	# This next command runs 3 independent chains without using multithreading.
	chain = mapreduce(c -> sample(model, NUTS(.65), 1000), chainscat, 1:3)
end

# ╔═╡ f91121b2-4e54-47a5-a038-b6d5c5781e94
plot(chain)

# ╔═╡ 22fe2c96-6e74-487d-902a-eb0c9a6f1dd8
chain_array = Array(chain);

# ╔═╡ 8cddb384-9604-4a90-90dd-fdb4e24c1b24
begin
	pl2 = plot()
	local resol
	for k in 1:100
		params_sample = chain_array[rand(1:size(chain_array, 1)), 2:end]
	    resol = solve(remake(prob1, u0=u0, p=params_sample), Tsit5(), saveat=Δt/10)
	    plot!(pl2, resol, alpha=0.1, color = "#BBBBBB", legend = false)
	end
	plot!(pl2, sol1, w=1, legend = false)
	scatter!(pl2, data_t, data');
	pl2
end

# ╔═╡ Cell order:
# ╠═54918c32-e4cf-11eb-04d2-35c34dcb1016
# ╠═26c8214c-13e6-4520-a2bc-13da053961cc
# ╠═e565eefb-115b-41dc-ac5b-5b72e8ab993d
# ╠═acc1c643-fe73-4072-9896-641316150ad4
# ╠═75be9104-8023-41b5-8b53-8f37929765ec
# ╠═f95497db-7bcb-479b-9ef7-ed20ad87ddd9
# ╠═9376ab45-046b-4bb7-b5b0-1a64afe68181
# ╠═f91121b2-4e54-47a5-a038-b6d5c5781e94
# ╠═22fe2c96-6e74-487d-902a-eb0c9a6f1dd8
# ╠═8cddb384-9604-4a90-90dd-fdb4e24c1b24
