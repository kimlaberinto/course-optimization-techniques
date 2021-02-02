### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 1abc28f0-6516-11eb-2671-31a6fbf782e7
begin
	using LinearAlgebra
	using Plots
end

# ╔═╡ 468a2484-6516-11eb-14a6-c55cd458a893
function find_quadratic_coefs(xs ,ys)
	# A*coefs = y
	A = [xs[i]^power for i in 1:3, power in 2:-1:0]
	p, q, r = A \ ys
	return (p, q, r)
end

# ╔═╡ b8a947b6-6516-11eb-34b0-859b13081caa
function interpolated_quadratic(interval_begin, interval_end, xs, ys)
	p, q, r = find_quadratic_coefs(xs, ys)
	x = range(interval_begin, interval_end, length=30)
	y = @. p*x^2 + q*x + r
	x, y
end

# ╔═╡ 0ade54c8-6518-11eb-20cb-0b674f957dd3
begin
	xs = [0, 1, 2]
	ys = [0, 1, 4]
	x, y = interpolated_quadratic(-5, 5, xs, ys)
	
	plot(x, y, label="Quadratic Fitting")
	scatter!(xs, ys, label="Data")
end

# ╔═╡ Cell order:
# ╠═1abc28f0-6516-11eb-2671-31a6fbf782e7
# ╠═468a2484-6516-11eb-14a6-c55cd458a893
# ╠═b8a947b6-6516-11eb-34b0-859b13081caa
# ╠═0ade54c8-6518-11eb-20cb-0b674f957dd3
