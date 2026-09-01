# run a notebook .ipynb in-process: evaluate each code cell in Main,
# capture stdout, embed the last value (plots -> png) as cell outputs
using JSON
using Base64

function run_notebook(path::String)
    nb = JSON.parsefile(path)
    exec = 0
    for cell in nb["cells"]
        cell["cell_type"] == "code" || continue
        exec += 1
        cell["execution_count"] = exec
        src = cell["source"] isa Vector ? join(cell["source"]) : cell["source"]
        outfile = tempname()
        result = nothing
        err = nothing
        open(outfile, "w") do f
            try
                redirect_stdout(f) do
                    result = Base.include_string(Main, src, "In[$exec]")
                end
            catch e
                err = e
            end
        end
        out = isfile(outfile) ? read(outfile, String) : ""
        rm(outfile, force=true)
        outputs = Any[]
        isempty(out) || push!(outputs, Dict("output_type"=>"stream", "name"=>"stdout", "text"=>out))
        if err !== nothing
            push!(outputs, Dict("output_type"=>"error", "ename"=>string(typeof(err)),
                                "evalue"=>sprint(showerror, err),
                                "traceback"=>[sprint(showerror, err)]))
            println("ERROR in cell $exec of $path: ", sprint(showerror, err))
            cell["outputs"] = outputs
            break
        end
        # embed returned plot object as png
        if result !== nothing && occursin("Plots", string(typeof(result)))
            pngbuf = IOBuffer()
            try
                Base.invokelatest(show, pngbuf, MIME("image/png"), result)
                b64 = base64encode(take!(pngbuf))
                push!(outputs, Dict("output_type"=>"display_data",
                                    "data"=>Dict("image/png"=>b64),
                                    "metadata"=>Dict()))
            catch e
                println("plot embedding failed: ", e)
            end
        end
        cell["outputs"] = outputs
    end
    open(path, "w") do f
        JSON.print(f, nb, 1)
    end
    println("executed: ", path)
end

for p in ARGS
    run_notebook(p)
end
