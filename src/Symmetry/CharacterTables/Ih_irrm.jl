
function generate_irrm_I()
    symels = pg_to_symels("I")
    mtable = build_mult_table(symels)
    # Define C5 as z-axis, and C2 as the C2 axis on an edge of the z+-pentagon in the yz plane. The z+-pentagon points to the -y region
    sq3 = sqrt(3)
    sq5 = sqrt(5)
    α1p = (1+(1/sqrt(5)))/4
    α1m = (1-(1/sqrt(5)))/4
    bα1p = sqrt((5*α1p)/2)
    bα1m = sqrt((5*α1m)/2)
    T1_E = [Float64(i==j) for i=1:3, j=1:3]
    T2_E = [Float64(i==j) for i=1:3, j=1:3]
    G_E = [Float64(i==j) for i=1:4, j=1:4]
    H_E = [Float64(i==j) for i=1:5, j=1:5]
    T1_C5 = [sq5*α1m bα1p 0;-bα1p sq5*α1m 0;0 0 1]
    T1_C2 = [-1 0 0;0 -1/sq5 2/sq5; 0 2/sq5 1/sq5]
    T2_C5 = [1 0 0;0 -sq5*α1p -bα1m;0 bα1m -sq5*α1p]
    T2_C2 = [-1/sq5 0 -2/sq5;0 -1 0;-2/sq5 0 1/sq5]
    G_C5 = [sq5*α1m bα1p 0 0;-bα1p sq5*α1m 0 0;0 0 -sq5*α1p bα1m;0 0 -bα1m -sq5*α1p]
    G_C2 = [0 0 -1 0;0 2/sq5 0 1/sq5;-1 0 0 0;0 1/sq5 0 -2/sq5]
    H_C5 = [1 0 0 0 0;0 -sq5*α1p bα1m 0 0;0 -bα1m -sq5*α1p 0 0;0 0 0 sq5*α1m -bα1p;0 0 0 bα1p sq5*α1m]
    H_C2 = [-1/5 -2*sq3/5 0 2*sq3/5 0;-2*sq3/5 3/5 0 2/5 0;0 0 1/sq5 0 -2/sq5;2*sq3/5 2/5 0 3/5 0;0 0 -2/sq5 0 -1/sq5]
    # 9 = 14*25*25*25*25
    # 25 = 14*14*14*59
    strangs = Vector{Int64}[[1]]
    knowns = [9,14,25,59]
    for i = 2:60
        if i ∉ knowns
            r = el_to_prod(i,knowns,1,mtable,0,0,Int64[])
            push!(strangs, r)
            append!(knowns, i)
            m = mprod(mtable, r)
            if m != i
                throw(ErrorException("POOPOO"))
            end
        else
            push!(strangs, [i])
            m = i
            r = [i]
        end
        #println("i = $i, resulting prod.: $m, ", r)
    end
    A = [Float64[1] for i=1:60]
    T1 = Matrix{Float64}[Float64[i==j for i=1:3, j=1:3] for i = 1:60]
    T2 = Matrix{Float64}[Float64[i==j for i=1:3, j=1:3] for i = 1:60]
    G = Matrix{Float64}[Float64[i==j for i=1:4, j=1:4] for i = 1:60]
    H = Matrix{Float64}[Float64[i==j for i=1:5, j=1:5] for i = 1:60]
    
    T1[1] = T1_E
    T2[1] = T2_E
    G[1] = G_E
    H[1] = H_E

    T1[14] = T1_C5
    T2[14] = T2_C5
    G[14] = G_C5
    H[14] = H_C5
    
    T1[59] = T1_C2
    T2[59] = T2_C2
    G[59] = G_C2
    H[59] = H_C2

    T1[25] = (T1_C5^3)*T1_C2
    T2[25] = (T2_C5^3)*T2_C2
    G[25] = (G_C5^3)*G_C2
    H[25] = (H_C5^3)*H_C2
    
    T1[9] = T1_C5*(T1[25]^4)
    T2[9] = T2_C5*(T2[25]^4)
    G[9] = G_C5*(G[25]^4)
    H[9] = H_C5*(H[25]^4)
    
    skip_these = [1,9,14,25,59]
    for i = 1:60
        if i ∉ skip_these
            lil_strang = strangs[i]
            for j in lil_strang
                T1[i] *= T1[j]
                T2[i] *= T2[j]
                G[i] *= G[j]
                H[i] *= H[j]
            end
        end
    end
    irrep_mat = Dict("A"=>A,"T1"=>T1,"T2"=>T2,"G"=>G,"H"=>H)
    for irr in ["T1","T2","G","H"]
        for i = 1:60
            if abs(abs(LinearAlgebra.det(irrep_mat[irr][i]))-1) > 1e-12
                throw(ErrorException("Bad booty determinant in I"))
            end
        end
    end
    return irrep_mat
end

function generate_irrm_Ih()
    symels = pg_to_symels("Ih")
    mtable = build_mult_table(symels)
    irrm_I = generate_irrm_I()
    A = irrm_I["A"]
    T1 = irrm_I["T1"]
    T2 = irrm_I["T2"]
    G = irrm_I["G"]
    H = irrm_I["H"]
    Ag = [Float64[1] for i=1:120]
    T1g = Matrix{Float64}[Float64[i==j for i=1:3, j=1:3] for i = 1:120]
    T2g = Matrix{Float64}[Float64[i==j for i=1:3, j=1:3] for i = 1:120]
    Gg = Matrix{Float64}[Float64[i==j for i=1:4, j=1:4] for i = 1:120]
    Hg = Matrix{Float64}[Float64[i==j for i=1:5, j=1:5] for i = 1:120]
    Au = [Float64[1] for i=1:120]
    T1u = Matrix{Float64}[Float64[i==j for i=1:3, j=1:3] for i = 1:120]
    T2u = Matrix{Float64}[Float64[i==j for i=1:3, j=1:3] for i = 1:120]
    Gu = Matrix{Float64}[Float64[i==j for i=1:4, j=1:4] for i = 1:120]
    Hu = Matrix{Float64}[Float64[i==j for i=1:5, j=1:5] for i = 1:120]
    for i = 1:60
        Ag[i]  = A[i]
        T1g[i] = T1[i]
        T2g[i] = T2[i]
        Gg[i]  = G[i]
        Hg[i]  = H[i]
        Au[i]  = A[i]
        T1u[i] = T1[i]
        T2u[i] = T2[i]
        Gu[i]  = G[i]
        Hu[i]  = H[i]
    end
    for i = 1:60
        j = mtable[61,i]
        Ag[j]  = A[i]
        T1g[j] = T1[i]
        T2g[j] = T2[i]
        Gg[j]  = G[i]
        Hg[j]  = H[i]
        Au[j]  = -A[i]
        T1u[j] = -T1[i]
        T2u[j] = -T2[i]
        Gu[j]  = -G[i]
        Hu[j]  = -H[i]
    end
    return Dict("Ag"=>Ag,"T1g"=>T1g,"T2g"=>T2g,"Gg"=>Gg,"Hg"=>Hg,"Au"=>Au,"T1u"=>T1u,"T2u"=>T2u,"Gu"=>Gu,"Hu"=>Hu)
end

function mprod(mtable, V)
    r = V[1]
    for v = 2:length(V)
        r = mtable[r,V[v]]
    end
    return r
end

function el_to_prod(element,knowns,kidx,mtable,depth,counter,strang)
    if depth > 500
        throw(ErrorException("'They delved too greedily and too deep, and disturbed that from which they fled, Durin's Bane.'"))
    end
    if counter > 4
        new_kidx = ((kidx+1)%length(knowns))+1
        return el_to_prod(element,knowns,new_kidx,mtable,depth,0,strang)
    end
    cidx = 0
    ridx = knowns[kidx]
    row = mtable[ridx,:]
    for i = 1:length(row)
        if mtable[ridx,i] == element
            cidx = i
            break
        end
    end
    if cidx == 0
        throw(ErrorException("It broke at line 171"))
    end
    append!(strang,ridx)
    if cidx in knowns
        return append!(strang, cidx)
    else
        return el_to_prod(cidx,knowns,kidx,mtable,depth+1,counter+1,strang)
    end
end