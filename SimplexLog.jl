function open_log(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64})
    fname = "Simplex.log"
    mode = isfile(fname) ? "a" : "w"
    stream = open(fname, mode)
    pwrite(stream, "=======================")
    pwrite(stream, "Comeco da Solucao do PL")
    pwrite(stream, "=======================")
    pwrite(stream, "Problema:")
    pwrite(stream, "A = $A")
    pwrite(stream, "b = $b")
    pwrite(stream, "c = $c")
    pwrite(stream, "")
    close(stream)
    nothing
end

function get_log(state::Int)
    fname = "SimplexFase2.log"
    stream = open(fname, "a")
    if state == 1
        pwrite(stream, "Simplex Fase 1")
        pwrite(stream, "--------------")
    else
        pwrite(stream, "Simplex Fase 2")
        pwrite(stream, "--------------")
    end
    return stream
end

function simplex_log(it::Int, x::Vector{Float64}, bidx::Vector{Int}, nidx::Vector{Int}, z::Float64, status::Int, stream::IOStream, debug=true)
    pwrite(stream, "iter $it:", debug)
    pwrite(stream, "x = $x", debug)
    pwrite(stream, "Base = $bidx", debug)
    pwrite(stream, "Nbase = $nidx", debug)
    pwrite(stream, "", debug)

    if status == 1
        pwrite(stream, "| Solucao otima obtida:", debug)
        pwrite(stream, "| ---------------------", debug)
        pwrite(stream, "| x = $x", debug)
        pwrite(stream, "| z = $z", debug)
        pwrite(stream, "| status = $status", debug)
        pwrite(stream, "", debug)
    elseif status == -1
        pwrite(stream, "| Solucao ilimitada obtida:", debug)
        pwrite(stream, "| -------------------------", debug)
        pwrite(stream, "| de = $x", debug)
        pwrite(stream, "| z = $z", debug)
        pwrite(stream, "| status = $status", debug)
        pwrite(stream, "", debug)
    elseif status == -2
        pwrite(stream, "| Solucao inviavel obtida:", debug)
        pwrite(stream, "| -------------------------", debug)
        pwrite(stream, "| z = $z", debug)
        pwrite(stream, "| status = $status", debug)
        pwrite(stream, "", debug)
    end
end

function pwrite(stream::IOStream, string::AbstractString, debug=true)
    if debug
        println(string)
    end
    write(stream, string * "\n")
end
