### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 1c20eab6-145e-4aac-9f19-be345740240e
import Pkg; Pkg.add("StaticArrays")

# ╔═╡ d463db70-31f4-11ef-15f5-ad6e77960b73
using DifferentialEquations, Plots, DataStructures, StaticArrays

# ╔═╡ 1a79db93-a1df-4a08-8ccc-5501ee03a7df
function nreed!(du,u,p,t)
    (μ,k,σ) = p
    du[1] = u[2] 
	du[2] = u[2]*(μ-u[2]^2)-k*u[1]*(1.0+σ*sqrt(k)*u[1]^2)
end

# ╔═╡ 614e6185-8426-414d-8610-8b71e42d5064
prob = ODEProblem(nreed!,[-1.;0],1280,[1.0 ,1.0, 0.18])

# ╔═╡ 4ab57402-989e-48f1-b0af-bc64296eedcc
@time sol = solve(prob,Tsit5(),saveat=0.02);

# ╔═╡ 92b9c457-5888-4fa5-a601-6b70a080ace5
@time sol2 = solve(prob,Tsit5());

# ╔═╡ e7d86a79-3965-4dbe-8a66-75df98cf59df
begin
	plot(sol,idxs=(0,1))
	plot!(sol,idxs=(0,2))
	scatter!(sol2.t,getindex.(sol2.u,1))
	scatter!(sol2.t,getindex.(sol2.u,2),xlims=(1200,1250))
end	

# ╔═╡ 49a82865-9073-4d03-b830-a1a14ba8b3d6
md"""
# With Remake
"""

# ╔═╡ fb6b3ab6-d27b-4f13-8f11-b90827e7223c
function solve_by_remake!(prob,buffer,dt)
	nbuf = length(buffer)
	u0 = prob.u0
	tmin = prob.tspan[1]
	tmax = tmin + dt*nbuf
	while tmax<prob.tspan[2]
		r_prob = remake(prob; u0=u0, tspan=(tmin,tmax), p=prob.p, saveat=dt)
		sol = solve(r_prob, Tsit5())
		for i in 1:nbuf
			for var in 1:2
				#push!(buffer,sol.u[i][var])
			end
		end	
		u0 = sol.u[end]
		tmin = sol.t[end]
		tmax = tmin + dt*nbuf
	end	
	return
end	

# ╔═╡ ae91bcc5-51f1-405f-aa86-e0a7f7b346e8
buffer = CircularBuffer{Float64}(100)

# ╔═╡ 83db9993-7371-4139-ad8e-8bf9e2b22029
@time solve_by_remake!(prob,buffer,0.02)

# ╔═╡ 6c724507-102a-41d9-aea0-0e0620a96e5a
md"""
# Integrator Interface
"""

# ╔═╡ 515894b0-4606-4508-aff4-c0d45b649386
function solve_by_steps(prob,dt)
	integrator = init(prob, Tsit5())
	t1 = prob.tspan[1]
	x = [integrator.u[1]]
	y = [integrator.u[2]]
	while t1<prob.tspan[2]
		step!(integrator)
		if integrator.t > t1+dt
			it = (t1+dt):dt:integrator.t
			append!(x,getindex.(integrator.sol(it).u,1))
			append!(y,getindex.(integrator.sol(it).u,2))
			t1 = it[end]
		end
	end	
	return x,y
end	

# ╔═╡ 2d856a6c-f1a2-4d68-9a05-9291e7023087
@time x,y = solve_by_steps(prob,0.02)

# ╔═╡ 98df0a47-9b23-4167-b516-5c1701d4de48
begin
	t2 = (0:length(x)-1)*0.05
	plot(sol,idxs=(0,1))
	plot!(sol,idxs=(0,2))
	plot!(t2,x)
	plot!(t2,y)
end	

# ╔═╡ 60d8722f-c582-40f0-89f1-ae5f1d400e65
function solve_by_steps!(buffer,prob,dt)
	integrator = init(prob, Tsit5())
	t = prob.tspan[1]
	while t<prob.tspan[2]
		step!(integrator)
		if integrator.t > t+dt
			it = (t+dt):dt:integrator.t
			append!(buffer,getindex.(integrator.sol(it).u,1))
			t = it[end]
		end
	end	
end	

# ╔═╡ c372b5fc-7601-43fc-8ab0-368185f44676
@time solve_by_steps!(buffer,prob,0.02)

# ╔═╡ 07a10700-ac8d-4fe9-beb1-ef925f49b33a
plot(buffer)

# ╔═╡ Cell order:
# ╠═1c20eab6-145e-4aac-9f19-be345740240e
# ╠═d463db70-31f4-11ef-15f5-ad6e77960b73
# ╠═1a79db93-a1df-4a08-8ccc-5501ee03a7df
# ╠═614e6185-8426-414d-8610-8b71e42d5064
# ╠═4ab57402-989e-48f1-b0af-bc64296eedcc
# ╠═92b9c457-5888-4fa5-a601-6b70a080ace5
# ╠═e7d86a79-3965-4dbe-8a66-75df98cf59df
# ╟─49a82865-9073-4d03-b830-a1a14ba8b3d6
# ╠═fb6b3ab6-d27b-4f13-8f11-b90827e7223c
# ╠═ae91bcc5-51f1-405f-aa86-e0a7f7b346e8
# ╠═83db9993-7371-4139-ad8e-8bf9e2b22029
# ╟─6c724507-102a-41d9-aea0-0e0620a96e5a
# ╠═515894b0-4606-4508-aff4-c0d45b649386
# ╠═2d856a6c-f1a2-4d68-9a05-9291e7023087
# ╠═98df0a47-9b23-4167-b516-5c1701d4de48
# ╠═60d8722f-c582-40f0-89f1-ae5f1d400e65
# ╠═c372b5fc-7601-43fc-8ab0-368185f44676
# ╠═07a10700-ac8d-4fe9-beb1-ef925f49b33a
