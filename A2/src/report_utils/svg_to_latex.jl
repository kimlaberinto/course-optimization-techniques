cd("assets")
#https://docs.julialang.org/en/v1/base/file/
#https://discourse.julialang.org/t/delete-all-files-with-a-specific-extension-in-a-particular-directory/38764/3
for (root, dirs, files) in walkdir(".")
    println("SVG Files in $root")
    svg_files = filter(x->endswith(x, ".svg"), files)
    for file in svg_files
        println(file) # path to files
    end
end