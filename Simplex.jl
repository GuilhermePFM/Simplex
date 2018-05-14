#
# Simplex Revisado
#

# Guilherme Pereira Freire Machado

function SimplexFase2(A, b, c) 
    stream = open_log()

    # separa em basicas/ nao basicas
    m, n = size(A)
    nb = n - m

    nidx = [i for i in 1:m]
    bidx = [i for i in (m+1):(n)]

    # inicia loop de solucao
    status = 3
    z = 0
    x = zeros(n)
    it = 0
    while status > 1
        it += 1

        # seleciona base
        B = A[:, bidx]
        N = A[:, nidx]
        
        # resolve as demais variaveis
        xn = zeros(length(nidx))
        d = B \ b
        D = B \ N
        xb = d - D * xn

        # separa os custos em basicos/nao basicos
        cn = c[nidx]
        cb = c[bidx]

        # calcula custos
        cr = cn' - cb' * D
        z = cb' * xb + (cr) * xn
        
        # testa otimo global
        if maximum(cr) < 0
            status = 1

            x = zeros(n)
            x[bidx] = xb
        end
        
        # testa ilimitado
        for i in 1:m
            if all(D[:, i] .<= zeros(m))
                status = -1

                x = zeros(n)
                x[nidx] = D[:, i]
            end
        end

        # seleciona proxima variavel basica
       nvbidx = indmax(cr)

       # seleciona a variavel basica que sai
       nvnidx = indmax(N[:, nvbidx] ./ d)

       # descobre os verdadeiros indices das variaveis
       old_bidx = bidx[:]
       bidx[nvnidx] = nidx[nvbidx] 
       nidx[nvbidx] = old_bidx[nvnidx]

       # escreve o log
       simplex_log(it, x, B, bidx, z, status, stream)
    end

    #x: o ponto ótimo (ou direção extrema se ilimitado),
    #z: função objetivo no ponto ótimo, e
    #status: 1 se ótimo e -1 se ilimitado.
    close(stream)
    return x, z, status
end

function SimplexFase1()

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

function simplex_log(it, x, B, bidx, z, status, stream)
    pwrite(stream, "iter $it:")
    pwrite(stream, "x = $x")
    pwrite(stream, "Base = $bidx")
    pwrite(stream, "Nbase = $bidx")
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