push!(LOAD_PATH,"../src/")

using Documenter, DMSuite

makedocs(sitename="DMSuite.jl",
         pages = ["Intro" => "index.md",
                  "Library" => "api.md",
                  "References" => "references.md"]
        )