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
	A = [x_val^power for x_val in xs, power in 2:-1:0]
	p, q, r = A \ ys
	return (p, q, r)
end

# ╔═╡ b8a947b6-6516-11eb-34b0-859b13081caa
function interpolated_quadratic(interval_begin, interval_end, xs, ys)
	p, q, r = find_quadratic_coefs(xs, ys)
	x = range(interval_begin, interval_end, length=30)
	y = @. p*x^2 + q*x + r
	x, y, p, q, r
end

# ╔═╡ a4c5b758-656a-11eb-151f-553e5ec5586d
function find_vertex_of_parabola(p, q, r)
	x_vertex = -q / (2*p)
	y_vertex = p*x_vertex^2 + q*x_vertex + r
	return (x_vertex, y_vertex)
end

# ╔═╡ 55cb058a-656b-11eb-1bc4-7d44883d115b


# ╔═╡ 0ade54c8-6518-11eb-20cb-0b674f957dd3
begin
	xs = [-2, 1, 2]
	ys = [5, 1, 4]
	x, y, p, q, r = interpolated_quadratic(-5, 5, xs, ys)

	x_vertex, y_vertex = find_vertex_of_parabola(p, q, r)
	
	plot(x, y, label="Quadratic Fitting")
	scatter!(xs, ys, label="Data")
	#scatter!(x_vertex, y_vertex, label="Vertex")
	scatter!([x_vertex], [y_vertex], label="Vertex")
end

# ╔═╡ Cell order:
# ╠═1abc28f0-6516-11eb-2671-31a6fbf782e7
# ╠═468a2484-6516-11eb-14a6-c55cd458a893
# ╠═b8a947b6-6516-11eb-34b0-859b13081caa
# ╠═a4c5b758-656a-11eb-151f-553e5ec5586d
# ╠═55cb058a-656b-11eb-1bc4-7d44883d115b
# ╠═0ade54c8-6518-11eb-20cb-0b674f957dd3
