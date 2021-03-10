using Random

const MIN_MAX_VAL = 2;
function generate_random_inits()
    rng = Random.MersenneTwister(1234);

    random_array = (rand(rng, 4, 5) .* 2 .- 1) * MIN_MAX_VAL;
    round.(random_array; digits = 2)
end

generate_random_inits()