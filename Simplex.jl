#
# Simplex Revisado
#

# Guilherme Pereira Freire Machado

function SimplexFase2(A, b, c) 
    stream = open_log()

    # separa em basicas/ nao basicas
    m, n = size(A)
    nv = n - m

    nidx = [i for i in 1:nv]
    bidx = [i for i in (nv+1):(n)]

    # inicia loop de solucao
    status = 3
    z = 0
    x = zeros(n)
    it = 0
    maxit = 100
    while status > 1
        it += 1

        # seleciona base
        B = A[:, bidx]
        N = A[:, nidx]
        
        # resolve as demais variaveis
        xn = zeros(length(nidx))
        d = B \ b
        db = B \ N

        xb = d - db * xn
        x[bidx] = xb
        x[nidx] = xn

        # separa os custos em basicos/nao basicos
        cn = c[nidx]
        cb = c[bidx]

        # calcula custos
        cr = cn' - cb' * db
        z = cb' * xb + (cr) * xn

        # escreve o log
       simplex_log(it, x, B, bidx, nidx, z, status, stream)

        # testa otimo global
        if maximum(cr) < 0
            status = 1

            x = zeros(n)
            x[bidx] = xb
            # escreve o log
            simplex_log(it, x, B, bidx, nidx, z, status, stream)
            break
        end
        
        # testa ilimitado
        for i in 1:nb
            if all(db[:, i] .<= zeros(m)) && cr[i] > 0
                status = -1

                x = zeros(n)
                x[nidx] = db[:, i]
                # escreve o log
                simplex_log(it, x, B, bidx, z, status, stream)
                break
            end
        end

        if it > maxit
            status = -2
            pwrite(stream, "Maximum number of iterations exceded! Solution not found.")
            break
        end

        # seleciona proxima variavel basica
       nvbidx = indmax(cr)

       # seleciona a variavel basica que sai
       nvnidx = indmin(xb ./ abs.(db[:,nvbidx]))

       # descobre os verdadeiros indices das variaveis
       old_bidx = bidx[:]
       bidx[nvnidx] = nidx[nvbidx] 
       nidx[nvbidx] = old_bidx[nvnidx]
    end

    #x: o ponto ótimo (ou direção extrema se ilimitado),
    #z: função objetivo no ponto ótimo, e
    #status: 1 se ótimo e -1 se ilimitado.
    close(stream)
    return x, z, status
end

function SimplexFase1(A, b, c)  
    status = 3 
    # pega n de variaveis basicas/ nao basicas originais
    m, n = size(A)
    nv = n - m # numero de variaveis nao basicas

    # adiciona a variavel de folga w
    Aw = -1 * ones(m, n + 1)
    Aw[:,1:n] = A
    bw = b
    cw = zeros(n+1)
    cw[n+1] = 1

    # escolhendo como basicas as slacks das desigualdades    
    x = zeros(n+1)
    bidx = [i for i in (nv+1):(n)]
    nidx = [i for i in 1:(n+1) if !(i in bidx)]

    it = 0
    while status > 1
        it += 1

        # seleciona base
        B = Aw[:, bidx]
        N = Aw[:, nidx]
        
        # resolve as demais variaveis
        xn = zeros(length(nidx))
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
        
        # seleciona proxima variavel basica
        if it == 1
            nvbidx = length(nidx)
        else
            nvbidx = indmax(cr)
        end
        # seleciona a variavel basica que sai
        nvnidx = indmin(xb ./ abs.(db[:,nvbidx]))
        
        # testa otimo global
        if minimum(cr) >= 0 && it > 1
            status = 1

            x = zeros(n)
            x[bidx] = xb
            # escreve o log
            # simplex_log(it, x, B, bidx, nidx, z, status, stream)
            break
        end

        # testa ilimitado
        for i in 1:m
            if all(db[:, i] .<= zeros(m)) && cr[i] > 0
                status = -1

                x = zeros(n)
                x[nidx] = db[:, i]
                # escreve o log
                # simplex_log(it, x, B, bidx, z, status, stream)
                break
                # return A, b, c, status
            end
        end

        # descobre os verdadeiros indices das variaveis
        old_bidx = deepcopy(bidx)
        bidx[nvnidx] = nidx[nvbidx] 
        nidx[nvbidx] = old_bidx[nvnidx]
    end

    # reorder original matrix
    orig_bidx = [i for i in bidx if i!=(n+1)]
    orig_nidx = [i for i in nidx if i!=(n+1)]

    A1 = zeros(size(A))
    A1[:, 1:length(orig_nidx)] = A[:,orig_nidx] 
    A1[:, (length(orig_nidx)+1):n] = A[:, orig_bidx] 
    c1 = zeros(n)
    c1[1:length(orig_nidx)] = c[orig_nidx]
    c1[(length(orig_nidx)+1):end] = c[orig_bidx]

    return A1, b, c1, status
end

function Simplex()
    A1, b, c1, status = SimplexFase1(A, b, c)
    A= A1
    c=c1
    if status == 1
       SimplexFase2(A1, b, c1)
    else
        println("Unbounded problem.")
    end

end
function open_log()
    fname = "SimplexFase2.log"
    if isfile(fname)
        stream = open(fname, "a")
        pwrite(stream, "Comeco da Solucao do PL")
        return stream
    else
        stream = open(fname, "w")
        pwrite(stream, "Comeco da Solucao do PL")
        return stream
    end
end

function simplex_log(it, x, B, bidx, nidx, z, status, stream)
    pwrite(stream, "iter $it:")
    pwrite(stream, "x = $x")
    pwrite(stream, "Base = $bidx")
    pwrite(stream, "Nbase = $nidx")
    pwrite(stream, "")
    
    if status == 1
        pwrite(stream, "Solucao otima obtida:")
        pwrite(stream, "x = $x")
        pwrite(stream, "z = $z")
        pwrite(stream, "status = $status")
        pwrite(stream, "")
    elseif status == -1
        pwrite(stream, "Solucao ilimitada obtida:")
        pwrite(stream, "x = $x")
        pwrite(stream, "z = $z")
        pwrite(stream,"status = $status")
        pwrite(stream, "")
    end
end

function pwrite(stream, string)
    println(string)
    write(stream, string * "\n")
end

# 2)

# a) Problema da Producao
A = [2 1 1 0; 1 2 0 1]
b = [4 ; 4]
c = [4 ; 3; 0; 0]
x,z,status = SimplexFase2(A, b, c)

# b) Prob 2
A = [0.5 -1 1 0; -4 1 0 1]
b = [0.5 ; 1]
c = [1 ; 1; 0; 0]
x,z,status = SimplexFase2(A, b, c)

# c) Prob 3 - fase 1
A = [2 1 1 0 0; 1 2 0 1 0; 1 1 0 0 -1]
b = [4 ; 4 ; 1]
c = [4 ; 3; 0; 0; 0]
x,z,status = Simplex(A, b, c)