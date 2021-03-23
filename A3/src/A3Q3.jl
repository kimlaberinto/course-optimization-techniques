using Plots
Plots.pyplot()

function objective(x::Array{T}) where T <: Real
    if  0 < x[1] < 4 &&
        0 < x[2] < 6 &&
        (x[1] + x[2] <= 8)
        return 2*x[1] + 5*x[2] 
    else
        return NaN
    end
end

x1s = -0.5:0.01:4.5
x2s = -0.5:0.01:6.5

plt = plot(; xlabel = "x1", ylabel = "x2", aspect_ratio = :equal)

# Plot contours of objective function
contour!(plt, x1s, x2s, (x1, x2) -> objective([x1, x2]), 
    levels = 30,
    c = :reds,
    lw = 3,
    label = "Contours of Objective")

# Plot variable bounds
vline!(plt,[0], label="A: x1 = 0")
vline!(plt,[4], label="B: x1 = 4")
hline!(plt,[0], label="C: x2 = 0")
hline!(plt,[6], label="D: x2 = 6")

# Plot more inequalities
# x1 + x2 <= 8
# x2 <= 8 - x1
plot!(plt, x1s, 8 .- x1s, label = "E: x1 + x2 = 8")

# Display
# Show maximized point
scatter!(plt, [2], [6], label = "P: LP Solution (2, 6)", markersize= 6)
xlims!(-0.5, 4.5)
ylims!(-0.5, 6.5)
plot!(size = (500, 600), title = "A3Q3 Graphical Solution")
savefig("assets/A3Q3_Plot.pdf")