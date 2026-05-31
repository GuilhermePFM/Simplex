#
# Simplex Revisado
#

# Guilherme Pereira Freire Machado

using LinearAlgebra

include("SimplexLog.jl")

function SimplexFase2(A::Array, b::Array, c::Array, debug=true; norder::Array = [0]) 
    stream = get_log(2)

    # separa em basicas/ nao basicas
    m, n = size(A)
    nv = n - m # numero de variaveis basicas
    
    nidx = [i for i in 1:nv]
    bidx = [i for i in (nv+1):(n)]
    
    # resolve valores iniciais de x
    x = zeros(n)
    x[bidx] = A[:, bidx] \ b

    # inicia loop de solucao
    status = 3
    z = 0.0
    it = 0
    maxit = 100
    while status > 1
        it += 1
        
        # seleciona base
        B = A[:, bidx]
        N = A[:, nidx]
        
        # matriz direcao
        BAj = B \ N
        db = BAj
        
        # separa os custos em basicos/nao basicos
        cn = c[nidx]
        cb = c[bidx]
        
        xb = x[bidx]
        xn = x[nidx]
        
        # calcula custos reduzidos como vetor coluna
        cr = cn - db' * cb

        z = dot(cb, xb) + dot(cr, xn)

        # escreve o log
        simplex_log(it, x, bidx, nidx, z, status, stream, debug)
                
        # testa otimo global
        if maximum(cr) <= 0
            println("otimo")
            status = 1
            if length(norder) > 1
                println("norder")
                new_x = zero(x)
                for i in 1:length(x)
                    new_x[norder[i]] = x[i]
                end
               x = new_x
            end
            # escreve o log
            pwrite(stream, "Reordenando: ")
            simplex_log(it, x, bidx, nidx, z, status, stream, debug)
            break
        end
        
        if it > maxit
            status = -2
            pwrite(stream, "Maximum number of iterations exceded! Solution not found.")
            break
        end
        
        # seleciona proxima variavel basica
        j = argmax(cr)
        
        # seleciona a variavel basica que sai
        r = xb ./ db[:, j]
        i = 1
        if it == 1
            r[r .<= 0] .= NaN
            i = argmin(ifelse.(isnan.(r), Inf, r))
        else
            i = lexicographic_smallest(N, db[:, j])
        end

        theta = r[i]
        
        # testa ilimitado
        if isnan(theta)
            println("ilimitado")
            status = -1
            
            d = - db[:, j]
            d[nidx[j]] = 1
           
            # # escreve o log
            simplex_log(it+1, d, bidx, nidx, Inf, status, stream, debug)
            break
        end

        # atualiza x 
        x[nidx[j]] = theta # novo j
        x[bidx] = x[bidx] - theta * db[:, j]
       
        # atualiza indices das variaveis
        old_bidx = bidx[:]
        bidx[i] = nidx[j] 
        nidx[j] = old_bidx[i]
    end

    #x: o ponto ótimo (ou direção extrema se ilimitado),
    #z: função objetivo no ponto ótimo, e
    #status: 1 se ótimo e -1 se ilimitado.
    close(stream)
    return x, z, status, it
end

function SimplexFase1(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64}, debug=true) 
    stream = get_log(1)

    # pega n de variaveis basicas/ nao basicas originais
    m, n = size(A)
    nv = n - m # numero de variaveis basicas
    
    # elimina  linhas com b negativo
    A[b.<0,:] .= -A[b.<0,:]
    b[b.<0]   .= -b[b.<0]

    # adiciona a variavel de folga w
    Aw = [A Matrix{Float64}(I, m, m)]
    bw = b
    cw = [zeros(n); -ones(m)]
    
    # escolhendo como basicas as slacks das desigualdades    
    x = zeros(n+m)
    art_idx = collect((n+1):(n+m))
    nidx = [i for i in (1):(size(Aw)[2] - size(Aw)[1])]
    bidx = [i for i in 1:(n+m) if !(i in nidx) ]
    xn = zeros(length(nidx))
        
    status = 3 
    it = 0
    while status > 1
        it += 1

        # seleciona base
        B = Aw[:, bidx]
        N = Aw[:, nidx]
        
        # resolve as demais variaveis
        d = B \ b
        db = B \ N

        xb = d - db * xn
        x[bidx] = xb
        x[nidx] = xn

        # separa os custos em basicos/nao basicos
        cn = cw[nidx]
        cb = cw[bidx]

        # calcula custos reduzidos como vetor coluna
        cr = cn - db' * cb
        z = dot(cb, xb) + dot(cr, xn)
        simplex_log(it, x, bidx, nidx, z, status, stream, debug)

        # seleciona proxima variavel basica
        j = argmax(cr)

        # seleciona a variavel basica que sai
        r = xb ./ db[:, j]
        i = 1
        if it == 1
            r[r .<= 0] .= NaN
            r[r .== Inf] .= NaN
            i = argmin(ifelse.(isnan.(r), Inf, r))
        else
            i = lexicographic_smallest(N, db[:, j])
        end

        # testa otimo global
        if maximum(cr) <= 0 && it > 1
            status = 1

            x = zeros(n+m)
            x[bidx] = xb

            if dot(cb, xb) > 0
                print()
                status = -2
            end

            # escreve o log
            simplex_log(it, x, bidx, nidx, z, status, stream, debug)
            break
        end
       
        # descobre os verdadeiros indices das variaveis
        old_bidx = deepcopy(bidx)
        bidx[i] = nidx[j] 
        nidx[j] = old_bidx[i]
    end

    # fix artificial in basis
    bidx, nidx = fix_artificial_var_in_B(Aw, bidx, nidx, art_idx)

    # reorder original matrix
    bidx_out = [i for i in bidx if i<=(n)]
    nidx_out = [i for i in nidx if i<=(n)]

    A1 = zeros(size(A))
    A1[:, 1:length(nidx_out)] = A[:,nidx_out] 
    A1[:, (length(nidx_out)+1):end] = A[:, bidx_out] 
    c1 = zeros(n)
    c1[1:length(nidx_out)] = c[nidx_out]
    c1[(length(nidx_out)+1):end] = c[bidx_out]
    
    new_idx = ones(Int, n)
    new_idx[1:length(nidx_out)] = nidx_out
    new_idx[ (length(nidx_out)+1):n] = bidx_out

    if debug
        pwrite(stream, "Novas matrizes:")
        pwrite(stream, "A = $A1")
        pwrite(stream, "c = $c1")
    end

    close(stream)
    return A1, b, c1, status, new_idx
end

function lexicographic_smallest(N, Aj)
    # N = [1  0  5  3; 2  4  6  -1; 3  0  7  9]
    # Aj = [3 -1 9]
    
    r = zeros(size(N))
    for i in 1:size(N)[1]
        if Aj[i] != 0
            r[i,:] = N[i,:] ./ Aj[i]
        else
            r[i,:] .= NaN
        end
    end
    
    # inicia vetor — ignora a primeira linha se tiver NaN
    index = findfirst(i -> all(r[i,:] .>= 0) && all(isfinite.(r[i,:])), 1:size(r,1))
    index = isnothing(index) ? 1 : index
    smallest = r[index,:]
    for i in (index+1):size(r)[1] # elimina NaN caso exista na inicializacao
        if all(r[i,:] .>= 0) && all(isfinite.(r[i,:])) && Tuple(r[i,:]) < Tuple(smallest)
            smallest = r[i,:]
            index = i
        end
    end

    return index
end

function fix_artificial_var_in_B(Aw, bidx, nidx, art_idx)
    # eliminates redundant
    art_in_b = [i for i in art_idx if i in bidx] 
    for i in art_in_b
        B = Aw[:, bidx]
        test = B \ Aw

        l = findall(==(i), bidx)
        
        if all(test[l,:] .== 0) # removes row
            Aw[l,:] .= 0

        else # change basis

            candidates = findall(!iszero, test[l[1],:])
            candidates = [i for i in candidates if i in nidx && !(i in art_idx)]
            j = candidates[1]
            
            old_bidx = deepcopy(bidx)
            bidx[l[1]] = nidx[j] 
            nidx[j] = old_bidx[l[1]]
            end
    end
    return bidx, nidx
end

function Simplex(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64}, debug=true)
    # Init logger
    open_log(A, b, c)

    A1, b, c1, status, norder = SimplexFase1(A, b, c, debug)
    A=A1
    c=c1
    if status == 1
       return SimplexFase2(A1, b, c1, debug, norder=norder)
    else
        println("Unfeasible problem.")
    end
end

function problemas()
    cd(pwd())
    # 2)

    # a) Problema da Producao
    println("a) Problema da Producao")
    println("")
    A = Float64[2 1 1 0; 1 2 0 1]
    b = Float64[4; 4]
    c = Float64[4; 3; 0; 0]
    x,z,status = SimplexFase2(A, b, c)
    
    # b) Prob 2
    println("b) Problema ilimitado")
    println("")
    A = Float64[0.5 -1 1 0; -4 1 0 1]
    b = Float64[0.5; 1]
    c = Float64[1; 1; 0; 0]
    x,z,status = SimplexFase2(A, b, c)

    # c) Prob 3 - fase 1
    println("c) Problema fase 1")
    println("")
    A = Float64[2 1 1 0 0; 1 2 0 1 0; -1 -1 0 0 1]
    b = Float64[4; 4; -1]
    c = Float64[4; 3; 0; 0; 0]
    x,z,status = Simplex(A, b, c)
end

#problemas()
