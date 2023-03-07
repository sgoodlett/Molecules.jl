struct Symel
    symbol::String
    rrep::Matrix{Float64}
end

struct CharacterTable
    name::String
    irreps::Vector{String}
    classes::Vector{String}
    class_orders::Vector{Int64}
    characters::Matrix{Float64}
    irrep_dims::Dict{String, Int64}
end

struct SymText
    pg::String
    symels::Vector{Molecules.Symmetry.CharacterTables.Symel}
    ctab::Molecules.Symmetry.CharacterTables.CharacterTable
    class_map::Vector{Int64}
    atom_map::Matrix{Int64}
    mult_table::Matrix{Int64}
    order::Int64
end

struct PG
    str
    family
    n
    subfamily
end

"""
Reduce n by i by dividing both by their greatest common divisor
"""
function reduce(n, i)
    g = gcd(n, i)
    return div(n, g), div(i, g)
end

"""
A quick implementation of the Euclid algorithm for finding the greatest common divisor
"""
function gcd(A, B)
    a = max(A,B)
    b = min(A,B)
    if a == 0
        return b
    elseif b == 0
        return a
    else
        r = a % b
        gcd(b, r)
    end
end

function Base.:(==)(A::Symel, B::Symel)
    if sum(A.rrep .- B.rrep) < tol
        return true
    else
        return false
    end
end

using Formatting
function string_repr(ctab::CharacterTable)
    l = length(ctab.classes) + 1
    longstrang = "-"^(l*10)*"\n"
    out = longstrang
    out *= "$(ctab.name) Character Table:\n"*(" "^15)
    for i = 1:length(ctab.classes)
        out *= "$(format("{:10s}",ctab.classes[i]))"
    end
    out *= "\n"
    for i = 1:length(ctab.irreps)
        out *= "$(format("{:10s}",ctab.irreps[i]))"
        for j = 1:length(ctab.characters[i,:])
            out *= "$(format("{:>10.3f}",ctab.characters[i,j]))"
        end
        out *= "\n"
    end
    out *= longstrang
    out *= "Irrep. dims.: $(ctab.irrep_dims)\n"
    return out
end

function string_repr(symels::Vector{Symel})
    out = ""
    for s in symels
        out *= "$(format("{:>8s}",s.symbol))"*": $(s.rrep)\n"
    end
    return out
end

function string_repr(symtext::SymText)
    out = "SymText for the $(symtext.pg) Point Group (order $(symtext.order))\n"
    out *= string_repr(symtext.ctab)
    out *= "\nSymmetry Elements (Symels)\n"
    out *= string_repr(symtext.symels)
    out *= "\nMultiplication table (row×column)\n"
    h = length(symtext.symels)
    out *= " "^(6)*"|"
    for i = 1:h
        out *= "$(format("{:>5d}",symtext.mult_table[1,i]))"
    end
    longstrang = "-"^6*"|"*"-"^(5*h)
    out *= "\n"*longstrang*"\n"
    for i = 1:h
        out *= "$(format("{:>5d}",symtext.mult_table[i,1])) |"
        for j = 1:h
            out *= "$(format("{:>5d}",symtext.mult_table[i,j]))"
        end
        out *= "\n"
    end
    out *= "\nClass Map\n$(symtext.class_map)\n\nAtom Map\n"
    c,d = size(symtext.atom_map)
    out *= " "^(9)*"|"
    for s = 1:d
        out *= "$(format("{:>8s}",symtext.symels[s].symbol))"
    end
    longstrang2 = "-"^(9)*"|"*"-"^(8*d)
    out *= "\n"*longstrang2*"\n"
    for a = 1:c
        out *= "$(format("{:>8d} |",a))"
        for b = 1:d
            out *= "$(format("{:>8d}",symtext.atom_map[a,b]))"
        end
        out *= "\n"
    end
    return out
end

function show(io::IO, ::MIME"text/plain", symels::Vector{Symel})
    print(io, string_repr(symels))
end

function show(io::IO, ::MIME"text/plain", ctab::CharacterTable)
    print(io, string_repr(ctab))
end

function show(io::IO, ::MIME"text/plain", symtext::SymText)
    print(io, string_repr(symtext))
end

function getindex(ctab::CharacterTable, irrep::String)
    return ctab.characters[findall(x->x==irrep, ctab.irreps)[1],:]
end