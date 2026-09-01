using Pkg

# docs 环境默认只含 Documenter；主包 TEMPO（含本地路径依赖
# ExpExp / ImpurityModelBase / QuAPI）通过 develop 挂载进来。
Pkg.activate(@__DIR__)
if Base.find_package("TEMPO") === nothing
    Pkg.develop(PackageSpec(path = dirname(@__DIR__)))
end
Pkg.instantiate()

using Documenter
using TEMPO

DocMeta.setdocmeta!(TEMPO, :DocTestSetup, :(using TEMPO, ImpurityModelBase, LinearAlgebra); recursive = true)

makedocs(;
    modules = [TEMPO],
    sitename = "TEMPO.jl",
    authors = "Guo Chu <guochu604b@gmail.com>",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://guochu.github.io/TEMPO.jl",
        edit_link = "master",
        assets = String[],
        size_threshold = nothing,   # api.md 单页收录全部 docstring，体积较大
    ),
    pages = [
        "首页" => "index.md",
        "快速上手" => "quickstart.md",
        "使用手册" => "manual.md",
        "实践指南" => "practice.md",
        "实现细节" => "internals.md",
        "API 参考" => "api.md",
    ],
    checkdocs = :nothing,
    linkcheck = false,
)

deploydocs(;
    repo = "github.com/guochu/TEMPO.git",
    devbranch = "master",
)
