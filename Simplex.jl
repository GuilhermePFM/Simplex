#
# Simplex Revisado
#

# Guilherme Pereira Freire Machado

function SimplexFase2(A::Array, b::Array, c::Array) 
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
    z = 0
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
        
        # calcula custos
        cr = cn' - cb' * db

        z = cb' * xb + (cr) * xn

        # escreve o log
        simplex_log(it, x, bidx, nidx, z, status, stream)
                
        # testa otimo global
        if maximum(cr) <= 0
            status = 1

            x = zeros(n)
            x[bidx] = xb

            # escreve o log
            simplex_log(it, x, bidx, nidx, z, status, stream)
            break
        end
        
        if it > maxit
            status = -2
            pwrite(stream, "Maximum number of iterations exceded! Solution not found.")
            break
        end
        
        # seleciona proxima variavel basica
        j = indmax(cr)
        
        # seleciona a variavel basica que sai
        r = xb ./ db[:, j]
        r[r .< 0] = NaN
        i = indmin(r)
        theta = r[i]
        
        # testa ilimitado
        if isnan(theta)
            status = -1
            
            d = - db[:, j]
            d[nidx[j]] = 1
           
            # # escreve o log
            simplex_log(it+1, d, bidx, nidx, Inf, status, stream)
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
    return x, z, status
end

function SimplexFase1(A::Array{Float64,2}, b::Array{Float64,1}, c::Array{Float64,1}) 
    stream = get_log(1)

    # pega n de variaveis basicas/ nao basicas originais
    m, n = size(A)
    nv = n - m # numero de variaveis basicas
    
    # adiciona a variavel de folga w
    Aw = -1 * ones(m, n + 1)
    Aw[:,1:n] = A
    bw = b
    cw = zeros(n+1)
    cw[n+1] = -1
    
    # escolhendo como basicas as slacks das desigualdades    
    x = zeros(n+1)
    bidx = [i for i in (nv+1):(n)]
    nidx = [i for i in 1:(n+1) if !(i in bidx) ]
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

        # calcula custos
        cr = cn' - cb' * db
        z = cb' * xb + (cr) * xn
        simplex_log(it, x, bidx, nidx, z, status, stream)

        # seleciona proxima variavel basica
        if it == 1
            j = length(nidx)
        else
            j = indmax(cr)
        end

        # seleciona a variavel basica que sai
        r = xb ./ db[:, j]
        r[r .< 0] = NaN
        i = indmin(r)
        
        # testa otimo global
        if maximum(cr) <= 0 && it > 1
            status = 1

            x = zeros(n+1)
            x[bidx] = xb

            # escreve o log
            simplex_log(it, x, bidx, nidx, z, status, stream)
            break
        end
       
        # descobre os verdadeiros indices das variaveis
        old_bidx = deepcopy(bidx)
        bidx[i] = nidx[j] 
        nidx[j] = old_bidx[i]
    end

    # reorder original matrix
    bidx_out = [i for i in bidx if i!=(n+1)]
    nidx_out = [i for i in nidx if i!=(n+1)]

    A1 = zeros(size(A))
    A1[:, 1:length(nidx_out)] = A[:,nidx_out] 
    A1[:, (length(nidx_out)+1):n] = A[:, bidx_out] 
    c1 = zeros(n)
    c1[1:length(nidx_out)] = c[nidx_out]
    c1[(length(nidx_out)+1):end] = c[bidx_out]

    close(stream)
    return A1, b, c1, status
end

function Simplex(A::Array{Float64,2}, b::Array{Float64,1}, c::Array{Float64,1})
    # Init logger
    open_log(A, b, c)

    A1, b, c1, status = SimplexFase1(A, b, c)
    
    if status == 1
       SimplexFase2(A1, b, c1)
    else
        println("Unfeasible problem.")
    end
end

function open_log(A::Array{Float64,2}, b::Array{Float64,1}, c::Array{Float64,1})
    fname = "SimplexFase2.log"
    if isfile(fname)
        stream = open(fname, "a")
        pwrite(stream, "=======================")
        pwrite(stream, "Comeco da Solucao do PL")
        pwrite(stream, "=======================")
        pwrite(stream, "Problema:")
        pwrite(stream, "A = $A")
        pwrite(stream, "b = $b")
        pwrite(stream, "c = $c")
        pwrite(stream, "")
        close(stream)
    else
        stream = open(fname, "w")
        pwrite(stream, "=======================")
        pwrite(stream, "Comeco da Solucao do PL")
        pwrite(stream, "=======================")
        pwrite(stream, "Problema:")
        pwrite(stream, "A = $A")
        pwrite(stream, "b = $b")
        pwrite(stream, "c = $c")
        pwrite(stream, "")
        close(stream)
    end
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

function simplex_log(it::Int, x::Array{Float64,1}, bidx::Array{Int,1}, nidx::Array{Int,1}, z::Float64, status::Int, stream::IOStream)
    pwrite(stream, "iter $it:")
    pwrite(stream, "x = $x")
    pwrite(stream, "Base = $bidx")
    pwrite(stream, "Nbase = $nidx")
    pwrite(stream, "")
    
    if status == 1
        pwrite(stream, "| Solucao otima obtida:")
        pwrite(stream, "| ---------------------")
        pwrite(stream, "| x = $x")
        pwrite(stream, "| z = $z")
        pwrite(stream, "| status = $status")
        pwrite(stream, "")
    elseif status == -1
        pwrite(stream, "| Solucao ilimitada obtida:")
        pwrite(stream, "| -------------------------")
        pwrite(stream, "| de = $x")
        pwrite(stream, "| z = $z")
        pwrite(stream, "| status = $status")
        pwrite(stream, "")
    end
end

function pwrite(stream::IOStream, string::AbstractString)
    println(string)
    write(stream, string * "\n")
end

# 2)

# a) Problema da Producao
A = float([2 1 1 0; 1 2 0 1])
b = float([4 ; 4])
c = float([4 ; 3; 0; 0])
x,z,status = SimplexFase2(A, b, c)

# b) Prob 2
A = float([0.5 -1 1 0; -4 1 0 1])
b = float([0.5 ; 1])
c = float([1 ; 1; 0; 0])
x,z,status = SimplexFase2(A, b, c)

# c) Prob 3 - fase 1
A = float([2 1 1 0 0; 1 2 0 1 0; -1 -1 0 0 1])
b = float([4 ; 4 ; -1])
c = float([4 ; 3; 0; 0; 0])
x,z,status = Simplex(A, b, c)

