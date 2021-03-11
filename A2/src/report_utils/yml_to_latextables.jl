cd("assets")
#https://docs.julialang.org/en/v1/base/file/
#https://discourse.julialang.org/t/delete-all-files-with-a-specific-extension-in-a-particular-directory/38764/3

import YAML

begin 
    output_file = open("yml_to_tables_InitialVector1_Summary.txt", "w")
    for (root, dirs, files) in walkdir(".")
        println("YML Files in $root")
        png_files = filter(x->endswith(x, ".yml"), files)
        for file in png_files
            println(joinpath(root, file))
            write(output_file, join([file, "\n"]))

            dataInFile = YAML.load_file(joinpath(root, file))
            
            trialData = dataInFile["Initial Vector 1"]
            final_vector = trialData["final_vector"]
            final_loss = trialData["final_loss"]
            N_f_evals = trialData["N_f_evals"]
            N_grad_evals = trialData["N_grad_evals"]
            N_hessian_evals = trialData["N_hessian_evals"]
            N_linsys_solves = trialData["N_linsys_solves"]

            write(output_file, "final_vector = $final_vector\n")
            write(output_file, "final_loss = $final_loss\n")
            write(output_file, "N_f_evals = $N_f_evals\n")
            write(output_file, "N_grad_evals = $N_grad_evals\n")
            write(output_file, "N_hessian_evals = $N_hessian_evals\n")
            write(output_file, "N_linsys_solves = $N_linsys_solves\n")
            write(output_file, "\n\n")
        end
    end
    close(output_file)
end