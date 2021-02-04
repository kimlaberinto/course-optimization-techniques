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

# ╔═╡ 0ade54c8-6518-11eb-20cb-0b674f957dd3
begin
	xs = [-2, 1, 2]
	ys = [-5, -10, 4]
	x, y, p, q, r = interpolated_quadratic(-5, 5, xs, ys)

	x_vertex, y_vertex = find_vertex_of_parabola(p, q, r)
	
	plot(x, y, label="Quadratic Fitting")
	scatter!(xs, ys, label="Samples")
	#scatter!(x_vertex, y_vertex, label="Vertex")
	scatter!([x_vertex], [y_vertex], label="Minimum of Parabola")
end

# ╔═╡ 3ba6feca-6593-11eb-16c3-8f7bdfc7733f
function Himmelblau(x, y)
	# https://en.wikipedia.org/wiki/Himmelblau%27s_function
	
	(x^2 + y - 11)^2 + (x + y^2 - 7)^2
end

# ╔═╡ ffe0cb68-656b-11eb-20d5-75ef4f130bd7
#Visually looks unimodal
begin
	xs_func = -7:0.01:7
	ys_func = @. Himmelblau(xs_func,5.5)
	
	plot(xs_func, ys_func, legend=false)
	title!("Himmelblau(x, y=5.5)")
	xlabel!("x")
	ylabel!("Himmelblau function")
end

# ╔═╡ 533b2a96-6595-11eb-2c04-0fb5b472a717
begin
	samples_x = [-6, 0, 6]
	samples_y = @. Himmelblau(samples_x, 5.5)
	
	#Fit to a parabola
	x_fitted, y_fitted, p_fit, q_fit, r_fit = interpolated_quadratic(-6, 6, samples_x, samples_y)
	
	#Find the optimum of parabola
	vertex_fitted_x, vertex_fitted_y = find_vertex_of_parabola(p_fit, q_fit, r_fit)
	
	#Plot original
	plot(xs_func, ys_func, label="function")
	
	#Plot the samples
	scatter!(samples_x, samples_y, label="samples")
	
	plot!(x_fitted, y_fitted, label="parabola")
	scatter!([vertex_fitted_x], [vertex_fitted_y], label="guessed optimum")
end

# ╔═╡ Cell order:
# ╠═1abc28f0-6516-11eb-2671-31a6fbf782e7
# ╠═468a2484-6516-11eb-14a6-c55cd458a893
# ╠═b8a947b6-6516-11eb-34b0-859b13081caa
# ╠═a4c5b758-656a-11eb-151f-553e5ec5586d
# ╠═0ade54c8-6518-11eb-20cb-0b674f957dd3
# ╠═3ba6feca-6593-11eb-16c3-8f7bdfc7733f
# ╠═ffe0cb68-656b-11eb-20d5-75ef4f130bd7
# ╠═533b2a96-6595-11eb-2c04-0fb5b472a717
