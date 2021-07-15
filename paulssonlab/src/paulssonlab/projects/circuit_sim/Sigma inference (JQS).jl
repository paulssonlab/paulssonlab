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
# rn = @reaction_network begin
# 	μ, ∅ --> S
# 	ν, ∅ --> R
# 	k, S --> S + G
# 	γ, S --> ∅
# 	γ, R --> ∅
# 	γ, G --> ∅
# 	η, S + R --> ∅
# 	hillr(S, α, K, n), S --> 2S
# end μ ν k γ η α K n
rn = @reaction_network begin
	μ, ∅ --> S
	#ν, ∅ --> R
	k, S --> S + G
	γ, S --> ∅
	#γ, R --> ∅
	γ, G --> ∅
	#η, S + R --> ∅
	hillr(S, α, K, n), S --> 2S
end μ k γ α K n

# ╔═╡ 3aba3645-72fe-4f62-be6b-fcde666f352e
begin
	odesys = convert(ODESystem, rn)
	latexify(odesys)
end

# ╔═╡ be7d4fd8-004b-420d-87c3-b9135b0a3949
[:S => 10, :R => 10, :G => 0]

# ╔═╡ 7e931d11-256c-43a6-9de0-2f6f2b1c68e7
paramsmap(odesys)

# ╔═╡ acc1c643-fe73-4072-9896-641316150ad4
begin
	#pmap     = [:μ => 5, :ν => 10, :k => 5, :γ => 1, :η => 100, :α => 20, :K => 3, :n => 4]
	#u0map = [:S => 10, :R => 10, :G => 0]
	#u0map = map((x,y) -> Pair(x,y), species(rs), u0)
	#pmap  = map((x,y) -> Pair(x,y), params(rs), p)
	##p 	  = [5, 10, 5, 1, 100, 20, 3, 4]
	##u0    = [10.0, 10.0, 0.0]
	p 	  = [5, 5, 1, 20, 3, 4]
	u0    = [1.0, 0.0]
	tspan = (0.0, 5)
	Δt 	  = 0.001
	prob1 = ODEProblem(odesys, u0, tspan, p)
	sol   = solve(prob1, Tsit5())
end

# ╔═╡ 75be9104-8023-41b5-8b53-8f37929765ec
plot(sol, lw=2)

# ╔═╡ f95497db-7bcb-479b-9ef7-ed20ad87ddd9
begin
	σ = 1.
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

	@model function fit(data, prob0)
		# σ ~ InverseGamma(3, 200)
		# μ ~ truncated(Normal(10, 10), 0, 20)
		# ν ~ truncated(Normal(10, 10), 0, 20)
		# k ~ truncated(Normal(10, 10), 0, 20)
		# γ ~ truncated(Normal(1, 1), 0, 20)
		# η ~ truncated(Normal(100, 100), 0, 200)
		# α ~ truncated(Normal(10, 10), 0, 40)
		# K ~ truncated(Normal(3, 3), 0, 5)
		# n ~ truncated(Normal(4, 4), 0, 5)
		#[5, 10, 5, 1, 100, 20, 3, 4]
		#σ ~ InverseGamma(2, 3)
		σ ~ truncated(Normal(1, 2), 0.1, Inf)
	    μ ~ truncated(Normal(5, 1), 0.1, Inf)
		#ν ~ truncated(Normal(10, 0.1), 0, Inf)
		k ~ truncated(Normal(5, 1), 0.1, Inf)
		γ ~ truncated(Normal(1, 0.5), 0.1, Inf)
		#η ~ truncated(Normal(100, 0.1), 0, Inf)
		α ~ truncated(Normal(20, 10), 1, Inf)
		K ~ truncated(Normal(3, 1), 1, Inf)
		n ~ truncated(Normal(4, 1), 1, Inf)

		##params = [μ, ν, k, γ, η, α, K, n]
		params = [μ, k, γ, α, K, n]
	    prob = remake(prob0, p=params)
	    predicted = solve(prob, Tsit5(), saveat=data_t)

		# is this good?
		if predicted.retcode != :Success
			Turing.@addlogprob! -Inf
			return
		end

	    for i = 1:length(predicted)
	        data[:,i] ~ MvNormal(predicted[i], σ)
	    end
	end

	model = fit(data, prob1)
end

# ╔═╡ 2e9633c6-0fdc-4985-aa39-0d71533ee899
begin
	init_params = [1, 1, 1, 1, 10, 2, 2]

	# This next command runs 3 independent chains without using multithreading.
	chain = mapreduce(c -> sample(model, NUTS(100, .65), 100, init_params=init_params, progress=true), chainscat, 1:3)
end

# ╔═╡ f91121b2-4e54-47a5-a038-b6d5c5781e94
plot(chain)

# ╔═╡ 22fe2c96-6e74-487d-902a-eb0c9a6f1dd8
chain_array = Array(chain);

# ╔═╡ 4d4a1916-b06b-4808-bcad-d9572aa07ba6
chain_array

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
# ╠═3aba3645-72fe-4f62-be6b-fcde666f352e
# ╠═be7d4fd8-004b-420d-87c3-b9135b0a3949
# ╠═7e931d11-256c-43a6-9de0-2f6f2b1c68e7
# ╠═acc1c643-fe73-4072-9896-641316150ad4
# ╠═75be9104-8023-41b5-8b53-8f37929765ec
# ╠═f95497db-7bcb-479b-9ef7-ed20ad87ddd9
# ╠═9376ab45-046b-4bb7-b5b0-1a64afe68181
# ╠═2e9633c6-0fdc-4985-aa39-0d71533ee899
# ╠═f91121b2-4e54-47a5-a038-b6d5c5781e94
# ╠═22fe2c96-6e74-487d-902a-eb0c9a6f1dd8
# ╠═4d4a1916-b06b-4808-bcad-d9572aa07ba6
# ╠═8cddb384-9604-4a90-90dd-fdb4e24c1b24
