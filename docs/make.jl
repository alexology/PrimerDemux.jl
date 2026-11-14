using Documenter
using PrimerDemux

DocMeta.setdocmeta!(PrimerDemux, :DocTestSetup, :(using PrimerDemux); recursive=true)

makedocs(;
    modules=[PrimerDemux],
    sitename="PrimerDemux.jl",
    authors="Alex Laini",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://alexology.github.io/PrimerDemux.jl",
        edit_link="main",
    ),
    pages=[
        "Home" => "index.md",
        "Validation" => "validation.md",
        "API" => "api.md",
    ],
    checkdocs=:exports,
)

deploydocs(; repo="github.com/alexology/PrimerDemux.jl", devbranch="main")
