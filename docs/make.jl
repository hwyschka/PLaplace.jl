using PLaplace, MinFEM
using Documenter, DocStringExtensions

Home = "Home" => "index.md"

Algorithm = "Algorithm" => "algorithm.md"

Examples = "Examples" => [
    "A Dirichlet Problem" => "examples/dirichlet-square.md",
    "A Neumann Problem" => "examples/neumann-square.md",
    "A Vector-Valued Problem" => "examples/vector-neumann-cube.md"
]

Library = "Library" => [
    "Public"    => "lib/public.md",
    "Internal"  => "lib/internal.md"
]

License = "License" => "license.md"

PAGES = [
    Home,
    Algorithm,
    Examples,
    Library,
    License
]

FORMAT = Documenter.HTML(
    prettyurls = true,
    assets = ["assets/favicon.ico"]
)

makedocs(
    sitename = "PLaplace.jl",
    authors = "Henrik Wyschka",
    format = FORMAT,
    pages = PAGES
)

deploydocs(
    repo = "github.com/hwyschka/PLaplace.jl.git"
)
