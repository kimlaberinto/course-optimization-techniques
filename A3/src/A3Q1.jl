using JuMP
using GLPK

model = Model(GLPK.Optimizer);

# Define the 11 variables representing the food servings
# At the same time, set the lower bound to be 0. (negative servings are non-physical)
num_food = 11
lower_bound = zeros(num_food)
@variable(model, x[i = 1:num_food], lower_bound = lower_bound[i])


# Define the Nutritional Value Constraints
proteins = [12; 9; 9; 1; 5; 4; 0; 2; 1; 3; 3]
carbs = [0; 1; 12; 0; 0; 14; 2; 20; 15; 12; 17]
fat = [12; 10; 5; 2; 3; 0; 0; 8; 0; 0; 4]
vitA = [119; 113; 147; 0; 0; 119; 112; 0; 28; 0; 18]
vitB1 = [0; 0; 0.10; 0.04; 0.23; 0.18; 0.03; 0.09; 0.11; 0.06; 0.05]
vitB2 = [0.14; 0.11; 0.43; 0.02; 0.07; 0.13; 0.05; 0.01; 0.05; 0.03; 0.08]
vitC = [0; 0; 2; 0; 0; 102; 11; 5; 70; 0; 0]
fibre = [0; 0; 0; 0; 0; 5; 0.9; 0; 2.6; 1.4; 1.8]

@assert length(proteins) == num_food
@assert length(carbs) == num_food
@assert length(fat) == num_food
@assert length(vitA) == num_food
@assert length(vitB1) == num_food
@assert length(vitB2) == num_food
@assert length(vitC) == num_food
@assert length(fibre) == num_food

protein_req = 60
carb_req = 300
fat_req = 40
vitA_req = 800
vitB1_req = 1.0
vitB2_req = 1.2
vitC_req = 60
fibre_req = 10


@constraint(model, sum(proteins[i]*x[i] for i in 1:num_food) >= protein_req)
@constraint(model, sum(carbs[i]*x[i] for i in 1:num_food) >= carb_req)
@constraint(model, sum(fat[i]*x[i] for i in 1:num_food) >= fat_req)
@constraint(model, sum(vitA[i]*x[i] for i in 1:num_food) >= vitA_req)
@constraint(model, sum(vitB1[i]*x[i] for i in 1:num_food) >= vitB1_req)
@constraint(model, sum(vitB2[i]*x[i] for i in 1:num_food) >= vitB2_req)
@constraint(model, sum(vitC[i]*x[i] for i in 1:num_food) >= vitC_req)
@constraint(model, sum(fibre[i]*x[i] for i in 1:num_food) >= fibre_req)

# Define calorie minimization objective
calories = [158; 132; 128; 25; 49; 64; 11; 158; 62; 61; 104]
@assert length(calories) == num_food

@objective(model, Min, sum(calories[i]*x[i] for i in 1:num_food))

# Print Model
println(model)

# Optimize Model
optimize!(model)

println("Termination status : ", termination_status(model))
println("Primal status      : ", primal_status(model))

obj_value = objective_value(model)
println("Lowest Calorie Found (objective value) : ", obj_value)

x_value = value.(x)
println("x vector found: ", x_value)

@assert sum(proteins[i]*x_value[i] for i in 1:num_food) >= protein_req
@assert sum(carbs[i]*x_value[i] for i in 1:num_food) >= carb_req
@assert sum(fat[i]*x_value[i] for i in 1:num_food) >= fat_req
@assert sum(vitA[i]*x_value[i] for i in 1:num_food) >= vitA_req
@assert sum(vitB1[i]*x_value[i] for i in 1:num_food) >= vitB1_req
@assert sum(vitB2[i]*x_value[i] for i in 1:num_food) >= vitB2_req
@assert sum(vitC[i]*x_value[i] for i in 1:num_food) >= vitC_req
@assert sum(fibre[i]*x_value[i] for i in 1:num_food) >= fibre_req

# Show everything nicely
println("\n\n")
using PrettyTables

let
    header = ["Variable" "Energy" "Protein" "Carbs" "Fat" "Vitamin A" "Vitamin B1" "Vitamin B2" "Vitamin C" "Fibre"
              ""         "[kcal]" "[g]"     "[g]"   "[g]" "[RE]"      "[mg]"       "[mg]"       "[mg]"      "[g]"];

    table_data = hcat(x, calories, proteins, carbs, fat, vitA, vitB1, vitB2, vitC, fibre)

    pretty_table(table_data, header, title = "Nutritional Values of Food") #Show Table of Nutritional Constants
end

let
    header = ["Protein" "Carbs" "Fat" "Vitamin A" "Vitamin B1" "Vitamin B2" "Vitamin C" "Fibre"
              "[g]"     "[g]"   "[g]" "[RE]"      "[mg]"       "[mg]"       "[mg]"      "[g]"];

    table_data = hcat(protein_req, carb_req, fat_req, vitA_req, vitB1_req, vitB2_req, vitC_req, fibre_req)

    pretty_table(table_data, header, title = "Daily Nutritional Requirements") #Show Table of Nutritional Constants
end

let
    header = ["Protein" "Carbs" "Fat" "Vitamin A" "Vitamin B1" "Vitamin B2" "Vitamin C" "Fibre"
              "[g]"     "[g]"   "[g]" "[RE]"      "[mg]"       "[mg]"       "[mg]"      "[g]"];

    table_data = hcat(sum(proteins[i]*x_value[i] for i in 1:num_food), 
        sum(carbs[i]*x_value[i] for i in 1:num_food), 
        sum(fat[i]*x_value[i] for i in 1:num_food),
        sum(vitA[i]*x_value[i] for i in 1:num_food), 
        sum(vitB1[i]*x_value[i] for i in 1:num_food), 
        sum(vitB2[i]*x_value[i] for i in 1:num_food), 
        sum(vitC[i]*x_value[i] for i in 1:num_food), 
        sum(fibre[i]*x_value[i] for i in 1:num_food))

    pretty_table(table_data, header, title = "Optimized Diet Nutrition") #Show Table of Nutritional Constants
end

let
    header = string.(x)

    table_data = x_value'

    pretty_table(table_data, header, title = "Optimized Diet Servings") #Show Table of Nutritional Constants
end

let

    header = ["" "Servings" "Energy" "Protein" "Carbs" "Fat" "Vitamin A" "Vitamin B1" "Vitamin B2" "Vitamin C" "Fibre";
              "" ""         "[kcal]" "[g]"     "[g]"   "[g]" "[RE]"      "[mg]"       "[mg]"       "[mg]"      "[g]"];

    x_value_rounded = round.(x_value; digits=2)
    food_rows = []
    for i in [6, 8, 9, 11]
        push!(food_rows, [x[i], 
            x_value_rounded[i], 
            calories[i]*x_value_rounded[i], 
            proteins[i]*x_value_rounded[i], 
            carbs[i]*x_value_rounded[i], 
            fat[i]*x_value_rounded[i], 
            vitA[i]*x_value_rounded[i], 
            vitB1[i]*x_value_rounded[i], 
            vitB2[i]*x_value_rounded[i], 
            vitC[i]*x_value_rounded[i], 
            fibre[i]*x_value_rounded[i]])
    end

    diet_row = ["Diet" "" "" sum(proteins[i]*x_value_rounded[i] for i in 1:num_food) sum(carbs[i]*x_value_rounded[i] for i in 1:num_food) sum(fat[i]*x_value_rounded[i] for i in 1:num_food) sum(vitA[i]*x_value_rounded[i] for i in 1:num_food)  sum(vitB1[i]*x_value_rounded[i] for i in 1:num_food) sum(vitB2[i]*x_value_rounded[i] for i in 1:num_food) sum(vitC[i]*x_value_rounded[i] for i in 1:num_food) sum(fibre[i]*x_value_rounded[i] for i in 1:num_food)]

    required_row = ["Daily Requirement" "" "" protein_req carb_req fat_req vitA_req vitB1_req vitB2_req vitC_req fibre_req]

    table_data = vcat(food_rows'..., diet_row, required_row)

    # @show size(table_data)
    # @show table_data
    pretty_table(table_data, header, title = "Optimized Diet - Nutrition Summary"; formatters = ft_printf("%5.1f"))
end