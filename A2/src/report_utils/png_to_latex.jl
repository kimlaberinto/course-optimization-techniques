cd("assets")
#https://docs.julialang.org/en/v1/base/file/
#https://discourse.julialang.org/t/delete-all-files-with-a-specific-extension-in-a-particular-directory/38764/3

begin 
    output_file = open("png_to_latex.txt", "w")
    for (root, dirs, files) in walkdir(".")
        println("PNG Files in $root")
        png_files = filter(x->endswith(x, ".png"), files)
        for file in png_files
            file_clean = replace(file, "_"=>"\\_")
            file_path = joinpath(root, file)
            println(file)
            code_snippet = 
            """

            \\begin{figure}[H]
                \\centering
                \\includegraphics[width=0.5\\linewidth]{$file_path}
                \\caption{$file_clean}
                \\label{fig:$file}
            \\end{figure}
            """
            write(output_file, code_snippet)
        end
    end
    close(output_file)
end

