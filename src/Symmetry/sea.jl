export buildD, findSEA, checkSEA

struct DistanceMatrix
    atoms::Vector{Any}
    distances::AbstractArray
end

mutable struct SEA{I,F}
    label::String
    set::Vector{I}
    paxis::Vector{F}
end

function distance(A::Atom, B::Atom)
    a = A.xyz
    b = B.xyz
    return (((a[1]-b[1])^2) + ((a[2]-b[2])^2) + ((a[3]-b[3])^2))^0.5
end

"""
Build interatomic distance matrix
"""
function buildD(mol::Vector{<:Atom})
    len = size(mol)[1]
    atoms = []
    arr = zeros(Float64, (len,len))
    for i = 1:len
        push!(atoms, mol[i])
    end
    for i = 1:len, j = i+1:len
        arr[i,j] = distance(mol[i], mol[j])
        arr[j,i] = arr[i,j]
    end
    return DistanceMatrix(atoms, arr)
end

"""
Check if two distance matrix vectors are equivalent under permutation 
within some tolerance δ
"""
function checkSEA(A::Vector{Float64}, B::Vector{Float64}, δ)::Bool
    a = sort(A)
    b = sort(B)
    z = a - b
    for i in z
        if abs(i) < 10.0^(-δ)
            continue
        else
            return false
        end
    end
    return true
end

"""
Find sets of Symmetry Equivalent Atoms
"""
function findSEA(D::DistanceMatrix, δ::Int)
    len = size(D.distances)[1]
    equiv_atom_pairs = Tuple[]
    # Find all pairs of equivalent Atoms
    for i = 1:len, j = i+1:len
        if cmp(D.atoms[i].mass, D.atoms[j].mass) == 0
            c = checkSEA(D.distances[:,i], D.distances[:,j], δ)
            if c
                push!(equiv_atom_pairs, (i,j))
            end
        end
    end
    
    # Find all sets of Symmetry Equivalent Atoms
    already_checked = Int[]
    seas = SEA[]
    for i = 1:len
        if i in already_checked
            continue
        else
            sea_idxs = Int[]
            push!(sea_idxs,i)
        end

        for k in equiv_atom_pairs
            if i in k
                if i == k[1]
                    push!(sea_idxs, k[2])
                    push!(already_checked, k[2])
                else
                    push!(sea_idxs, k[1])
                    push!(already_checked, k[1])
                end
            end
        end
    sea = SEA("None", sea_idxs, [0.0, 0.0, 0.0])
    push!(seas, sea)
    end

    return seas
end


