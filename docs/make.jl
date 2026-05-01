using ArnoldiVandermonde
using Documenter

DocMeta.setdocmeta!(ArnoldiVandermonde, :DocTestSetup, :(using ArnoldiVandermonde); recursive=true)

makedocs(;
    modules=[ArnoldiVandermonde],
    authors="Toby Driscoll <driscoll@udel.edu> and contributors",
    sitename="ArnoldiVandermonde.jl",
    format=Documenter.HTML(;
        canonical="https://tobydriscoll.github.io/ArnoldiVandermonde.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tobydriscoll/ArnoldiVandermonde.jl",
    devbranch="main",
)
