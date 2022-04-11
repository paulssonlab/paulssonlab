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
		Pkg.PackageSpec(name="ANSIColoredPrinters", version="0.0.1"),
		Pkg.PackageSpec(name="BenchmarkTools", version="1.1.1"),
    ])
    using Plots, ModelingToolkit, Catalyst, LinearAlgebra, OrdinaryDiffEq, Latexify, DiffEqJump, MomentClosure, Turing, Distributions, DifferentialEquations, MCMCChains, StatsPlots, BenchmarkTools
end

# ╔═╡ 4c79939c-ee08-4c31-917d-cdeedd296094
begin
	using ANSIColoredPrinters

	function color_print(f)
		io = IOBuffer()
		f(IOContext(io, :color=>true))
		html_str = sprint(io2->show(io2, MIME"text/html"(),
						  HTMLPrinter(io, root_class="documenter-example-output")))
		HTML("$html_str")
	end

	HTML("""
<style>
html .content pre {
    font-family: "JuliaMono", "Roboto Mono", "SFMono-Regular", "Menlo", "Consolas",
        "Liberation Mono", "DejaVu Sans Mono", monospace;
}
html pre.documenter-example-output {
    line-height: 125%;
	font-size: 60%
}
html span.sgr1 {
    font-weight: bolder;
}
html span.sgr2 {
    font-weight: lighter;
}
html span.sgr3 {
    font-style: italic;
}
html span.sgr4 {
    text-decoration: underline;
}
html span.sgr7 {
    color: #fff;
    background-color: #222;
}
html.theme--documenter-dark span.sgr7 {
    color: #1f2424;
    background-color: #fff;
}
html span.sgr8,
html span.sgr8 span,
html span span.sgr8 {
    color: transparent;
}
html span.sgr9 {
    text-decoration: line-through;
}
html span.sgr30 {
    color: #111;
}
html span.sgr31 {
    color: #944;
}
html span.sgr32 {
    color: #073;
}
html span.sgr33 {
    color: #870;
}
html span.sgr34 {
    color: #15a;
}
html span.sgr35 {
    color: #94a;
}
html span.sgr36 {
    color: #08a;
}
html span.sgr37 {
    color: #ddd;
}
html span.sgr40 {
    background-color: #111;
}
html span.sgr41 {
    background-color: #944;
}
html span.sgr42 {
    background-color: #073;
}
html span.sgr43 {
    background-color: #870;
}
html span.sgr44 {
    background-color: #15a;
}
html span.sgr45 {
    background-color: #94a;
}
html span.sgr46 {
    background-color: #08a;
}
html span.sgr47 {
    background-color: #ddd;
}
html span.sgr90 {
    color: #888;
}
html span.sgr91 {
    color: #d57;
}
html span.sgr92 {
    color: #2a5;
}
html span.sgr93 {
    color: #d94;
}
html span.sgr94 {
    color: #08d;
}
html span.sgr95 {
    color: #b8d;
}
html span.sgr96 {
    color: #0bc;
}
html span.sgr97 {
    color: #eee;
}
html span.sgr100 {
    background-color: #888;
}
html span.sgr101 {
    background-color: #d57;
}
html span.sgr102 {
    background-color: #2a5;
}
html span.sgr103 {
    background-color: #d94;
}
html span.sgr104 {
    background-color: #08d;
}
html span.sgr105 {
    background-color: #b8d;
}
html span.sgr106 {
    background-color: #0bc;
}
html span.sgr107 {
    background-color: #eee;
}
</style>""")
end

# ╔═╡ 26c8214c-13e6-4520-a2bc-13da053961cc
begin
	using Random
	Random.seed!(14);
end

# ╔═╡ 4c37a1a0-486c-489a-844d-052599d37575
begin
	macro code_warntype_(args...)
		code = macroexpand(@__MODULE__, :(@code_warntype $(args...)))
		@assert code.head == :call
		insert!(code.args, 2, :io)
		esc(quote # non-hygenic :(
			color_print() do io
			    $code
			end
		end)
	end

	macro code_llvm_(args...)
		code = macroexpand(@__MODULE__, :(@code_llvm $(args...)))
		@assert code.head == :call
		insert!(code.args, 2, :io)
		esc(quote # non-hygenic :(
			color_print() do io
			    $code
			end
		end)
	end

	macro code_native_(args...)
		code = macroexpand(@__MODULE__, :(@code_native $(args...)))
		@assert code.head == :call
		insert!(code.args, 2, :io)
		esc(quote # non-hygenic :(
			color_print() do io
			    $code
			end
		end)
	end
end

# ╔═╡ e565eefb-115b-41dc-ac5b-5b72e8ab993d
rn = @reaction_network begin
	μ, ∅ --> S
	ν, ∅ --> R
	k, S --> S + G
	γ, S --> ∅
	γ, R --> ∅
	γ, G --> ∅
	η, S + R --> ∅
	hillr(S, α, K, n), S --> 2S
end μ ν k γ η α K n

# ╔═╡ 3aba3645-72fe-4f62-be6b-fcde666f352e
begin
	odesys = convert(ODESystem, rn)
	latexify(odesys)
end

# ╔═╡ acc1c643-fe73-4072-9896-641316150ad4
begin
	#pmap     = [:μ => 5, :ν => 10, :k => 5, :γ => 1, :η => 100, :α => 20, :K => 3, :n => 4]
	#u0map = [:S => 10, :R => 10, :G => 0]
	#u0map = map((x,y) -> Pair(x,y), species(rs), u0)
	#pmap  = map((x,y) -> Pair(x,y), params(rs), p)
	p 	  = [5., 10., 5., 1., 2., 20., 3., 4.]
	u0    = [10.0, 10.0, 0.0]
	tspan = (0.0, 5.)
	Δt 	  = 0.01
	prob1 = ODEProblem(odesys, u0, tspan, p)
	sol   = solve(prob1, Tsit5())
end

# ╔═╡ 75be9104-8023-41b5-8b53-8f37929765ec
plot(sol, lw=2)

# ╔═╡ f95497db-7bcb-479b-9ef7-ed20ad87ddd9
begin
	σ = 0.1
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

	#@model function fit(ts::Array{T,1}, data::Array{T,2}, prob0) where {T,S}
	@model function fit(ts, data, prob0)
		σ ~ truncated(Normal(1., 10), 0.1, Inf)
	    μ ~ truncated(Normal(5., 10), 0.1, Inf)
		ν ~ truncated(Normal(5., 10), 0.1, Inf)
		#γ ~ truncated(Normal(1., 5), 0.1, Inf)
		η ~ truncated(Normal(1., 20), 0.1, Inf)
		#params = [μ, ν, γ, η]
		k = 5.
		γ = 1.
		η = 100.
		α = 20.
		K = 3.
		n = 4.
		params = [μ, ν, k, γ, η, α, K, n]
	    prob = remake(prob0, p=params)
		solution = solve(prob, Tsit5(), saveat=ts)
	    #predicted::Vector{Vector{T}} = solution.u
		predicted = solution.u

		# is this good?
		if solution.retcode != :Success
			Turing.@addlogprob! -Inf
			return
		end

	    for i = 1:length(predicted)
	        data[:,i] ~ MvNormal(predicted[i], σ)
	    end
		#data ~ MvNormal(Fill(predicted, length(predicted)), 0.2)
	end

	model = fit(data_t, data, prob1)
end

# ╔═╡ edd3b9bf-3b7f-492b-a977-ad069a9a050b
@code_warntype_ model.f(
    model,
    Turing.VarInfo(model),
    #Turing.SampleFromPrior(),
    Turing.DefaultContext(),
    model.args...,
)

# ╔═╡ 5f597d40-013a-4125-ad4a-ce306446d0e9
@benchmark model.f(
    model,
    Turing.VarInfo(model),
    #Turing.SampleFromPrior(),
    Turing.DefaultContext(),
    model.args...,
)

# ╔═╡ 2e9633c6-0fdc-4985-aa39-0d71533ee899
begin
	init_params = [0.1, 5., 5., 1.]

	# This next command runs 3 independent chains without using multithreading.
	chain = mapreduce(c -> sample(model, NUTS(100, .65), 200, init_params=init_params), chainscat, 1:3)
end

# ╔═╡ b690cb35-d99f-499e-bb0b-6340d94562e5
chain.value

# ╔═╡ 9a5d5be3-eecc-4231-bb02-84658bac794c
plot(chain.value[var=:lp,chain=2])

# ╔═╡ 308266b2-3f22-4df8-bba6-d33aa764a707
histogram(chain.value[var=:acceptance_rate,chain=3], bins=100)

# ╔═╡ f91121b2-4e54-47a5-a038-b6d5c5781e94
plot(chain)

# ╔═╡ fa89f0f0-040c-4fdd-b239-0ee06265bf57
chain.name_map

# ╔═╡ 160a25b4-1d06-480b-9ca7-87aa2033c931
chain

# ╔═╡ 22fe2c96-6e74-487d-902a-eb0c9a6f1dd8
chain_array = Array(chain);

# ╔═╡ 4d4a1916-b06b-4808-bcad-d9572aa07ba6
chain_array

# ╔═╡ 8cddb384-9604-4a90-90dd-fdb4e24c1b24
begin
	pl2 = plot()
	local resol
	for k in 1:30
		params_sample = chain_array[rand(1:size(chain_array, 1)), 2:end]
	    resol = solve(remake(prob1, u0=u0, p=params_sample), Tsit5(), saveat=Δt/10)
	    plot!(pl2, resol, alpha=0.1, color="#444", legend=false)
	end
	plot!(pl2, sol1, w=1, legend=false)
	scatter!(pl2, data_t, data', markerstrokewidth=0);
	pl2
end

# ╔═╡ Cell order:
# ╠═54918c32-e4cf-11eb-04d2-35c34dcb1016
# ╟─4c79939c-ee08-4c31-917d-cdeedd296094
# ╟─4c37a1a0-486c-489a-844d-052599d37575
# ╠═26c8214c-13e6-4520-a2bc-13da053961cc
# ╠═e565eefb-115b-41dc-ac5b-5b72e8ab993d
# ╠═3aba3645-72fe-4f62-be6b-fcde666f352e
# ╠═acc1c643-fe73-4072-9896-641316150ad4
# ╠═75be9104-8023-41b5-8b53-8f37929765ec
# ╠═f95497db-7bcb-479b-9ef7-ed20ad87ddd9
# ╠═9376ab45-046b-4bb7-b5b0-1a64afe68181
# ╠═edd3b9bf-3b7f-492b-a977-ad069a9a050b
# ╠═5f597d40-013a-4125-ad4a-ce306446d0e9
# ╠═b690cb35-d99f-499e-bb0b-6340d94562e5
# ╠═9a5d5be3-eecc-4231-bb02-84658bac794c
# ╠═308266b2-3f22-4df8-bba6-d33aa764a707
# ╠═2e9633c6-0fdc-4985-aa39-0d71533ee899
# ╠═f91121b2-4e54-47a5-a038-b6d5c5781e94
# ╠═fa89f0f0-040c-4fdd-b239-0ee06265bf57
# ╠═160a25b4-1d06-480b-9ca7-87aa2033c931
# ╠═22fe2c96-6e74-487d-902a-eb0c9a6f1dd8
# ╠═4d4a1916-b06b-4808-bcad-d9572aa07ba6
# ╠═8cddb384-9604-4a90-90dd-fdb4e24c1b24
