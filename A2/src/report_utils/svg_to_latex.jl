cd("assets")
for (root, dirs, files) in walkdir(".")
    println("SVG Files in $root")
    svg_files = filter(x->endswith(x, ".svg"), files)
    for file in svg_files
        println(file) # path to files
    end
end