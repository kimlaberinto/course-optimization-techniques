cd("assets")
#https://docs.julialang.org/en/v1/base/file/
#https://discourse.julialang.org/t/delete-all-files-with-a-specific-extension-in-a-particular-directory/38764/3

begin 
    output_file = open("svg_to_latex.txt", "w")
    for (root, dirs, files) in walkdir(".")
        println("SVG Files in $root")
        svg_files = filter(x->endswith(x, ".svg"), files)
        for file in svg_files
            println(file)
            code_snippet = 
            """

            \\begin{figure}[H]
                \\centering
                \\includesvg[width=0.5\\linewidth]{$file}
                \\caption{$file}
                \\label{fig:$file}
            \\end{figure}
            """
            write(output_file, code_snippet)
        end
    end
    close(output_file)
end

