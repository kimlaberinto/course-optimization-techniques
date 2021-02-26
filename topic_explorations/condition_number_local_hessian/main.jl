using ForwardDiff
using LinearAlgebra
using Plots

function rosenbrock_banana(x_input)
    return (1 - x_input[1])^2 + 0.1*100*(x_input[2] - x_input[1]^2)^2
end

function hessian_rosenbrock(x_input)
    return ForwardDiff.hessian(rosenbrock_banana, x_input)
end

function calculate_local_hessian_cond_number(x_input)
    hessian_matrix = hessian_rosenbrock(x_input)
    eigenvals, eigenvectors = eigen(hessian_matrix)
    eigenvals = @. abs(eigenvals)

    min_eigenval = minimum(eigenvals)
    max_eigenval = maximum(eigenvals)

    cond_number = abs(max_eigenval) / abs(min_eigenval)

    return cond_number
end

x_plot = -2.1:0.01:2.1
y_plot = -2.1:0.01:2.1

p1 = contour(x_plot, y_plot, (x, y) -> rosenbrock_banana([x, y]),  
    fill=false, 
    label = "Rosenbrock",
    c=:black,
    legend = :bottomright,
    levels = 0:5:500)
title!("Rosenbrock Function Contours")

p2 = contour(x_plot, y_plot, 
(x, y) -> log(calculate_local_hessian_cond_number([x, y])), 
fill=true,
c = cgrad(:hawaii, scale=:exp),
levels = -1:.5:60)

contour!(x_plot, y_plot, 
        (x, y) -> log(calculate_local_hessian_cond_number([x, y])), 
        levels = -1:.5:60,
        fill=false,
        c = :black)

title!("Log10 Condition Number of\nLocal Hessian Matrix")

layout_plot = @layout [a b];
plot(p1, p2, layout = layout_plot, size = (800, 500))

savefig("LocalLogConditionNumberPlot_Rosenbrock.svg")
savefig("LocalLogConditionNumberPlot_Rosenbrock.png")

plot!() #Show plot in VSCode