using GhilSellersEBM
using Documenter

DocMeta.setdocmeta!(GhilSellersEBM, :DocTestSetup, :(using GhilSellersEBM); recursive=true)

makedocs(;
    modules=[GhilSellersEBM],
    authors="Maximilian Gelbrecht <maximilian.gelbrecht@posteo.de and contributors",
    repo="https://github.com/maximilian-gelbrecht/GhilSellersEBM.jl/blob/{commit}{path}#{line}",
    sitename="GhilSellersEBM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://maximilian-gelbrecht.github.io/GhilSellersEBM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Model Description" => "model_description.md"
    ],
)

deploydocs(;
    repo="github.com/maximilian-gelbrecht/GhilSellersEBM.jl",
    devbranch="main",
)
