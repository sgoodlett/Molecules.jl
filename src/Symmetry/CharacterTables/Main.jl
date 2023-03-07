"""
Take in a point group string and return all symmetry elements for
that point group. For the groups Cn, Cnv, Cnh, S2n, Dn, Dnd, and Dnh
the nth order rotation axis is assumed to be the z-axis. Then, for Cnv,
the reflection planes are oriented such that the x-axis is contained in
one of the planes. For Dn, Dnd, and Dnh the x-axis is one of the n C2' axes.
The name of some symels will not match the typical conventions in Cotton's
character tables, for instance in the case of D2h, Cotton's C2'(y) is named
C2''. Similarly, σ(xz) is σv, σ(yz) is σd, and σ(xy) is σh.
The ordering of the symmetry elements is also not consistent with Cotton, 
though this issue is mostly corrected by the function `generate_symel_to_class_map`.
"""
function pg_to_symels(PG)
    arg_err_string = "An invalid point group has been given or unexpected parsing of the point group string has occured"
    pg = parse_pg_str(PG)
    symels = Vector{Symel}([Symel("E", [1 0 0; 0 1 0; 0 0 1])])
    σh = [1 0 0; 0 1 0; 0 0 -1]
    if pg.family == "C"
        if pg.subfamily == "h"
            push!(symels, Symel("σh", σh)) # σh
            if pg.n % 2 == 0
                push!(symels, Symel("i", i()))
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            symels = vcat(symels, cns, sns)
        elseif pg.subfamily == "v"
            cns = generate_Cn(pg.n)
            if pg.n % 2 == 0
                n = div(pg.n, 2)
                σds = generate_σd(n)
            else
                n = pg.n
                σds = Vector{Symel}([])
            end
            σvs = generate_σv(pg.n)
            symels = vcat(symels, cns, σvs, σds)
        elseif pg.subfamily == "s"
            push!(symels, Symel("σh", σh))
        elseif pg.subfamily == "i"
            push!(symels, Symel("i", i()))
        elseif isnothing(pg.subfamily)
            cns = generate_Cn(pg.n)
            symels = vcat(symels, cns)
        else
            throw(ArgumentError(arg_err_string))
        end
    elseif pg.family == "D"
        if pg.subfamily == "h"
            push!(symels, Symel("σh", σh))
            if pg.n % 2 == 0
                push!(symels, Symel("i", i()))
                n = pg.n >> 1
                σds = generate_σd(n)
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = vcat(c2ps, c2pps)
            else
                n = pg.n
                σds = Vector{Symel}([])
                c2s = generate_C2p(pg.n)
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n)
            σvs = generate_σv(pg.n)
            #c2s = generate_C2(pg.n)
            symels = vcat(symels, cns, c2s, sns, σvs, σds)
        elseif pg.subfamily == "d"
            if pg.n % 2 == 0
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = vcat(c2ps, c2pps)
            else
                c2s = generate_C2p(pg.n)
                push!(symels, Symel("i", i()))
            end
            cns = generate_Cn(pg.n)
            sns = generate_Sn(pg.n * 2, true)
            σds = generate_σd(pg.n)
            symels = vcat(symels, cns, sns, c2s, σds)
        elseif isnothing(pg.subfamily)
            cns = generate_Cn(pg.n)
            if pg.n % 2 == 0
                c2ps = generate_C2p(pg.n)
                c2pps = generate_C2pp(pg.n)
                c2s = vcat(c2ps, c2pps)
            else
                c2s = generate_C2p(pg.n)
            end
            symels = vcat(symels, cns, c2s)
        else
            throw(ArgumentError(arg_err_string))
        end
    elseif pg.family == "S"
        if isnothing(pg.subfamily) & (pg.n % 2 == 0)
            n = div(pg.n, 2)
            if n % 2 != 0
                push!(symels, Symel("i", i()))
            end
            cns = generate_Cn(n)
            sns = generate_Sn(pg.n, true)
            symels = vcat(symels, cns, sns)
        else
            throw(ArgumentError(arg_err_string))
        end
    else
        if pg.family == "T"
            if pg.subfamily == "h"
                Ths = generate_Th()
                symels = vcat(symels, Ths)
            elseif pg.subfamily == "d"
                Tds = generate_Td()
                symels = vcat(symels,Tds)
            else
                Ts = generate_T()
                symels = vcat(symels,Ts)
            end
        elseif pg.family == "O"
            if pg.subfamily == "h"
                Ohs = generate_Oh()
                symels = vcat(symels, Ohs)
            else
                Os = generate_O()
                symels = vcat(symels, Os)
            end
        elseif pg.family == "I"
            if pg.subfamily == "h"
                Ihs = generate_Ih()
                symels = vcat(symels, Ihs)
            else
                Is = generate_I()
                symels = vcat(symels, Is)
            end
        else
            throw(ArgumentError(arg_err_string))
        end
    end
    return symels
end

"""
Parses point group information from a string. Family indicates the large letter in
the Schönflies symbol (C,S,D,T,O, or I), the order is the natural number associated
with the group (3 in C3v), and the subfamily is the smaller letter following the order
(d in D3d, h in Th, etc.). If the Schönflies symbol has no associated order or subfamily,
`nothing` is returned.
"""
function parse_pg_str(s)
    re = r"([A-Z]+)(\d+)?([a-z]+)?"
    m = match(re, s)
    family, n, subfamily = m.captures
    if !isnothing(n)
        n = parse(Int, n)
    end
    if !isnothing(subfamily)
        subfamily = string(subfamily)
    end
    family = string(family)
    return PG(s, family, n, subfamily)
end

"""
Calculates the character table entries for a given Schönflies string.
As of now, complex valued characters are ignored and added with their associated irrep.
Returns a `CharacterTable` 
"""
function pg_to_chartab(PG)
    pg = parse_pg_str(PG)
    irreps = []
    if pg.family == "C"
        if pg.subfamily == "s"
            irreps = ["A'","A''"]
            classes = ["E", "σh"]
            chars = [1 1; 1 -1]
        elseif pg.subfamily == "i"
            irreps = ["Ag","Au"]
            classes = ["E", "i"]
            chars = [1 1; 1 -1]
        elseif pg.subfamily == "v"
            irreps, classes, chars = Cnv_irr(pg.n)
        elseif pg.subfamily == "h"
            irreps, classes, chars = Cnh_irr(pg.n)
        else
            irreps, classes, chars = Cn_irrmat(pg.n)
        end
    elseif pg.family == "D"
        if pg.subfamily == "d"
            irreps, classes, chars = Dnd_irr(pg.n)
        elseif pg.subfamily == "h"
            irreps, classes, chars = Dnh_irr(pg.n)
        else
            irreps, classes, chars = Dn_irr(pg.n)
        end
    elseif pg.family == "S"
        irreps, classes, chars = Sn_irr(pg.n)
    else
        cp3 = cos(pi/3)
        pr5 = 0.5*(1.0+sqrt(5.0))
        mr5 = 0.5*(1.0-sqrt(5.0))
        if pg.family == "T"
            if pg.subfamily == "h"
                irreps, classes, chars = (["Ag","Au","Eg","Eu","Tg","Tu"],
                 ["E", "4C_3", "4C_3^2", "3C_2", "i", "S_6", "S_6^5", "3σh"],
                 [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0;
                  1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0;
                  2.0  cp3  cp3  2.0  2.0  cp3  cp3  1.0;
                  2.0  cp3  cp3  2.0 -2.0 -cp3 -cp3 -1.0;
                  3.0  0.0  0.0 -1.0  1.0  0.0  0.0 -1.0;
                  3.0  0.0  0.0 -1.0 -1.0  0.0  0.0  1.0])
            elseif pg.subfamily == "d"
                irreps, classes, chars = (["A1","A2","E","T1","T2"],
                 ["E", "8C_3", "3C_2", "6S_4", "6σd"],
                 [1.0  1.0  1.0  1.0  1.0;
                  1.0  1.0  1.0 -1.0 -1.0;
                  2.0 -1.0  2.0  0.0  0.0;
                  3.0  1.0 -1.0  1.0 -1.0;
                  3.0 -1.0 -1.0 -1.0  1.0])
            else
                irreps, classes, chars = (["A","E","T"],
                 ["E", "4C_3", "4C_3^2", "3C_2"],
                 [1.0  1.0  1.0  1.0;
                  2.0  cp3  cp3  2.0;
                  3.0  0.0  0.0 -1.0])
            end
        elseif pg.family == "O"
            if pg.subfamily == "h"
                irreps, classes, chars = (["A1g","A2g","Eg","T1g","T2g","A1u","A2u","Eu","T1u","T2u"],
                 ["E", "8C_3", "6C_2", "6C_4", "3C_2", "i", "6S_4", "8S_6", "3σh", "6σd"],
                 [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0;
                  1.0  1.0 -1.0 -1.0  1.0  1.0 -1.0  1.0  1.0 -1.0;
                  2.0 -1.0  0.0  0.0  2.0  2.0  0.0 -1.0  2.0  0.0;
                  3.0  0.0 -1.0  1.0 -1.0  3.0  1.0  0.0 -1.0 -1.0;
                  3.0  0.0  1.0 -1.0 -1.0  3.0 -1.0  0.0 -1.0  1.0;
                  1.0  1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0 -1.0;
                  1.0  1.0 -1.0 -1.0  1.0 -1.0  1.0 -1.0 -1.0  1.0;
                  2.0 -1.0  0.0  0.0  2.0 -2.0  0.0  1.0 -2.0  0.0;
                  3.0  0.0 -1.0  1.0 -1.0 -3.0 -1.0  0.0  1.0  1.0;
                  3.0  0.0  1.0 -1.0 -1.0 -3.0  1.0  0.0  1.0 -1.0])
            else
                irreps, classes, chars = (["A1","A2","E","T1","T2"],
                 ["E", "6C_4", "3C_2", "8C_3", "6C_2"],
                 [1.0  1.0  1.0  1.0  1.0;
                  1.0 -1.0  1.0  1.0 -1.0;
                  2.0  0.0  2.0 -1.0  0.0;
                  3.0  1.0 -1.0  0.0 -1.0;
                  3.0 -1.0 -1.0  0.0  1.0])
            end
        elseif pg.family == "I"
            if pg.subfamily == "h"
                irreps, classes, chars = (["Ag","T1g","T2g","Gg","Hg","Au","T1u","T2u","Gu","Hu"],
                 ["E", "12C_5", "12C_5^2", "20C_3", "15C_2", "i", "12S_10", "12S_10^3", "20S_6", "15σ"],
                 [1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0  1.0;
                  3.0  pr5  mr5  0.0 -1.0  3.0  mr5  pr5  0.0 -1.0;
                  3.0  mr5  pr5  0.0 -1.0  3.0  pr5  mr5  0.0 -1.0;
                  4.0 -1.0 -1.0  1.0  0.0  4.0 -1.0 -1.0  1.0  0.0;
                  5.0  0.0  0.0 -1.0  1.0  5.0  0.0  0.0 -1.0  1.0;
                  1.0  1.0  1.0  1.0  1.0 -1.0 -1.0 -1.0 -1.0 -1.0;
                  3.0  pr5  mr5  0.0 -1.0 -3.0 -mr5 -pr5  0.0  1.0;
                  3.0  mr5  pr5  0.0 -1.0 -3.0 -pr5 -mr5  0.0  1.0;
                  4.0 -1.0 -1.0  1.0  0.0 -4.0  1.0  1.0 -1.0  0.0;
                  5.0  0.0  0.0 -1.0  1.0 -5.0  0.0  0.0  1.0 -1.0])
            else
                irreps, classes, chars = (["A","T1","T2","G","H"],
                 ["E", "12C_5", "12C_5^2", "20C_3", "15C_2"],
                 [1.0  1.0  1.0  1.0  1.0;
                  3.0  pr5  mr5  0.0 -1.0;
                  3.0  mr5  pr5  0.0 -1.0;
                  4.0 -1.0 -1.0  1.0  0.0;
                  5.0  0.0  0.0 -1.0  1.0])
            end
        else
            throw(ArgumentError("An invalid point group has been given or unexpected parsing of the point group string has occured"))
        end
    end
    class_orders = grab_class_orders(classes)
    irr_dims = Dict{String,Int}()
    for (irr_idx,irrep) in enumerate(irreps)
        irr_dims[irrep] = chars[irr_idx, 1]
    end
    return CharacterTable(PG, irreps, classes, class_orders, chars, irr_dims)
end

"""
Extract class orders into a vector
"""
function grab_class_orders(classes)
    ncls = length(classes)
    class_orders = zeros(Int64, ncls)
    for i = 1:ncls
        class_orders[i] = grab_order(classes[i])
    end
    return class_orders
end

"""
Determine class order from string
"""
function grab_order(class_str)
    re = r"^(\d+)"
    m = match(re, class_str)
    if !isnothing(m)
        return parse(Int64, m.captures[1])
    else
        return 1
    end
end

"""
A vector mapping symmetry elements to the classes they belong to.
"""
function generate_symel_to_class_map(symels::Vector{Symel}, ctab::CharacterTable)
    pg = parse_pg_str(ctab.name)
    if pg.n !== nothing
        ns = pg.n>>1 # pg.n floor divided by 2
    end
    ncls = length(ctab.classes)
    nsymel = length(symels)
    class_map = zeros(Int64, nsymel)
    class_map[1] = 1 # E is always first
    if pg.family == "C"
        if pg.subfamily == "s" || pg.subfamily == "i"
            class_map[2] = 2
        elseif pg.subfamily == "h"
            if pg.n % 2 == 0
                class_map[4:pg.n+2] .= 2:pg.n # C_n
                class_map[3] = pg.n+1 # i
                class_map[2] = pg.n+ns+1 # σh
                for i = pg.n+3:2*pg.n # S_n
                    if i > 3*ns+1
                        class_map[i] = i-ns
                    else
                        class_map[i] = i+ns-1
                    end
                end
            else
                for i = 2:pg.n+1 # C_n
                    class_map[i] = i-1
                end
                class_map[2] = pg.n+1 # σh
                for i = pg.n+2:2*pg.n # S_n
                    class_map[i] = i
                end
            end
        elseif pg.subfamily == "v"
            # The last class is σv (and then σd if n is even), and the last symels are also these!
            cn_class_map!(class_map, pg.n, 0, 0)
            if pg.n % 2 == 0
                class_map[end-pg.n+1:end-ns] .= ncls-1 # σv
                class_map[end-ns+1:end] .= ncls # σd
            else
                class_map[end-pg.n+1:end] .= ncls # σv
            end
        else
            class_map[2:end] .= 2:nsymel
        end
    elseif pg.family == "S"
        if pg.n % 4 == 0
            for i = 2:pg.n
                if i <= ns
                    class_map[i] = 2*i-1
                else
                    class_map[i] = 2*(i-ns)
                end
            end
        else
            class_map[2] = ns+1 # i
            class_map[3:ns+1] .= 2:ns # C_n
            for i = ns+2:pg.n # S_n
                if i > ns+(pg.n>>2)+1
                    class_map[i] = i-(ns-1)>>1
                else
                    class_map[i] = i + (ns-1)>>1
                end
            end
        end
    elseif pg.family == "D"
        if pg.subfamily == "h"
            if pg.n % 2 == 0
                class_map[2] = ncls-2 # σh
                class_map[3] = (ncls>>1)+1 # i
                cn_class_map!(class_map, pg.n, 2, 0) # Cn
                class_map[pg.n+3:3*ns+2] .= ns+2 # C2'
                class_map[3*ns+3:2*pg.n+2] .= ns+3 # C2''
                for i = 2*pg.n+3:3*pg.n+1 # Sn
                    if i > 3*pg.n-ns+1
                        class_map[i] = i-2*pg.n+3
                    else
                        class_map[i] = 3*pg.n+6-i
                    end
                end
                # The result of C2'×i changes depending if pg.n ≡ 0 (mod 4), 
                # but also D2h doesn't need to be flipped because I treated it special
                if pg.n % 4 == 0 || pg.n == 2
                    class_map[end-pg.n+1:end-ns] .= ncls-1 # σv
                    class_map[end-ns+1:end] .= ncls # σd
                else
                    class_map[end-pg.n+1:end-ns] .= ncls # σv
                    class_map[end-ns+1:end] .= ncls-1 # σd
                end
            else
                class_map[2] = (ncls>>1)+1
                cn_class_map!(class_map, pg.n, 1, 0)
                class_map[pg.n+2:2*pg.n+1] .= ns+2
                cn_class_map!(class_map, pg.n, 2*pg.n, ns+2)
                class_map[end-pg.n+1:end] .= ncls
            end
        elseif pg.subfamily == "d"
            if pg.n % 2 == 0
                cn_class_map!(class_map, pg.n, 0, 0) # Cn
                class_map[2:pg.n] .= 2*class_map[2:pg.n].-1 # Reposition Cn
                cn_class_map!(class_map, pg.n+1, pg.n-1, 0) # Sn
                class_map[pg.n+1:2*pg.n] .= 2*(class_map[pg.n+1:2*pg.n].-1) # Reposition Sn
                class_map[end-2*pg.n+1:end-pg.n] .= ncls-1 # C2'
                class_map[end-pg.n+1:end] .= ncls # σd
            else
                class_map[2] = (ncls>>1)+1 # i
                cn_class_map!(class_map, pg.n, 1, 0) # Cn
                for i = pg.n+2:2*pg.n # Sn
                    if i > pg.n+1+ns
                        class_map[i] = i+2-pg.n
                    else
                        class_map[i] = 2*pg.n+4-i
                    end
                end
                class_map[end-2*pg.n+1:end-pg.n] .= ns+2
                class_map[end-pg.n+1:end] .= ncls # σd
            end
        else
            cn_class_map!(class_map, pg.n, 0, 0) # Cn
            if pg.n % 2 == 0
                class_map[end-pg.n+1:end-ns] .= ncls-1 # Cn'
                class_map[end-ns+1:end] .= ncls # Cn''
            else
                class_map[end-pg.n+1:end] .= ncls # Cn
            end
        end
    else
        if pg.family == "T"
            if pg.subfamily == "h"
                class_map = [1,2,3,2,3,2,3,2,3,4,4,4,5,6,7,6,7,6,7,6,7,8,8,8]
            elseif pg.subfamily == "d"
                class_map = [1,2,2,2,2,2,2,2,2,3,3,3,5,5,5,5,5,5,4,4,4,4,4,4]
            else
                class_map = [1,2,3,2,3,2,3,2,3,4,4,4]
            end
        elseif pg.family == "O"
            if pg.subfamily == "h"
                class_map = [1,4,5,4,4,5,4,4,5,4,2,2,2,2,2,2,2,2,3,3,3,3,3,3,6,7,9,7,7,9,7,7,9,7,8,8,8,8,8,8,8,8,10,10,10,10,10,10]
            else
                class_map = [1,2,3,2,2,3,2,2,3,2,4,4,4,4,4,4,4,4,5,5,5,5,5,5]
            end
        elseif pg.family == "I"
            if pg.subfamily == "h"
                class_map = [1,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
                             6,7,8,8,7,7,8,8,7,7,8,8,7,7,8,8,7,7,8,8,7,7,8,8,7,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10]
            else
                class_map = [1,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,2,3,3,2,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5]
            end
        else
            throw(ArgumentError("Unrecognized point group: $(ctab.name)"))
        end
    end
    return class_map
end

"""
Auxilliary function for the `symel_to_class_map`, puts the Cn primary rotations into classes
of two except for C_2 (i.e. {C_6,C_6^5}, {C_3,C_3^2}, {C_2}).
"""
function cn_class_map!(class_map, n, idx_offset, cls_offset)
    for i = 2:n
        if i > (n>>1)+1
            class_map[i+idx_offset] = n-i+2+cls_offset
        else
            class_map[i+idx_offset] = i+cls_offset
        end
    end
end

"""
Given primary and secondary axes, rotate the molecule to the symmetry element frame
"""
function rotate_mol_to_symels(mol, paxis, saxis)
    φ, θ, χ = get_euler_angles(paxis, saxis)
    dc = direction_cosine_matrix(φ, θ, χ)
    new_mol = Molecules.transform(mol, dc)
    return new_mol
end

"""
Given primary and secondary axes, rotate the symels to the molecule frame. Rotating in this
way will render the irreducible representation matrices we provide incorrect.
"""
function rotate_symels_to_mol(symels, paxis, saxis)
    φ, θ, χ = get_euler_angles(paxis, saxis)
    dc = direction_cosine_matrix(φ, θ, χ)
    dct = deepcopy(dc)
    dct = transpose!(dct, dc)
    new_symels = Vector{Symel}([])
    for s in symels
        push!(new_symels, Symel(s.symbol, dct*s.rrep*dc))
    end
    return new_symels
end

"""
Determine the Euler angles for rotating (x,y,z) onto (saxis, paxis × saxis , paxis)
"""
function get_euler_angles(paxis, saxis)
    x = [1.0;0.0;0.0]
    y = [0.0;1.0;0.0]
    z = [0.0;0.0;1.0]
    ynew = paxis×saxis
    zxp = normalize!(z×paxis)
    if isnan(zxp[1])
        φ = 0.0
    else
        xproj = zxp⋅x
        if xproj <= 0
            φ = acos(y⋅zxp)
        else
            φ = 2*π-acos(y⋅zxp)
        end
    end
    rφ = Molecules.rotation_matrix(z,φ)
    yN = rφ*y
    xN = rφ*x
    θ = acos(z⋅paxis)
    rθ = Molecules.rotation_matrix(yN,θ)
    x2N = rθ*xN
    Z = rθ*z
    χ = acos(x2N⋅saxis)
    rχ = Molecules.rotation_matrix(Z,χ)
    X = rχ*x2N
    Y = rχ*yN
    return φ, θ, χ
end

"""
Form the direction cosine matrix for rotating by Euler angles φ, θ, χ
"""
function direction_cosine_matrix(φ, θ, χ)
    sp = sin(φ)
    cp = cos(φ)
    st = sin(θ)
    ct = cos(θ)
    sc = sin(χ)
    cc = cos(χ)
    direction_cosine = [cp*ct*cc-sp*sc sp*ct*cc+cp*sc -cc*st; -cp*ct*sc-sp*cc -sp*ct*sc+cp*cc sc*st; st*cp st*sp ct]
    return direction_cosine
end

"""
Determine where the symels send an atom. Returns an Natom by Nsymel array, where element (i,j)
is where the ith atom is sent under the jth symel.
"""
function get_atom_mapping(mol, symels)
    "symels after transformation"
    amap = zeros(Int, length(mol), length(symels))
    for (a, atom) in enumerate(mol)
        for (s, symel) in enumerate(symels)
            w = where_you_go(mol, atom, symel)
            if w !== nothing
                amap[a,s] = w
            else
                throw(ErrorException("Atom $(atom) not mapped to another atom under symel $(symel)"))
            end
        end
    end
    return amap
end

"""
A function that checks where an atom is sent under a symel.
"""
function where_you_go(mol, atom, symel)
    ratom = Atom(atom.Z, atom.mass, symel.rrep*atom.xyz)
    len = size(mol,1)
    for i = 1:len
        if isapprox(mol[i].xyz, ratom.xyz, atol=tol)
            return i
        end
    end
    return nothing
end

"""
Given a molecule file, determine the symtext for that molecule
"""
function symtext_from_file(fn)
    mol = Molecules.parse_file(fn)
    return symtext_from_mol(mol)
    #mol = Molecules.translate(mol, Molecules.center_of_mass(mol))
    #pg, paxis, saxis = Molecules.Symmetry.find_point_group(mol)
    #symels = pg_to_symels(pg)
    #symels = rotate_symels_to_mol(symels, paxis, saxis)
    #ctab = pg_to_chartab(pg)
    #class_map = generate_symel_to_class_map(symels, ctab)
    #atom_map = get_atom_mapping(mol, symels)
    #return SymText(pg, symels, ctab, class_map, atom_map)
end

"""
Given a molecule, determine the symtext for that molecule
"""
function symtext_from_mol(mol)
    mol = Molecules.translate(mol, Molecules.center_of_mass(mol))
    pg, paxis, saxis = Molecules.Symmetry.find_point_group(mol)
    symels = pg_to_symels(pg)
    #symels = rotate_symels_to_mol(symels, paxis, saxis)
    mol = rotate_mol_to_symels(mol, paxis, saxis)
    #mol = Molecules.transform(mol, Molecules.Cn([0 0 1],3)^2)
    ctab = pg_to_chartab(pg)
    class_map = generate_symel_to_class_map(symels, ctab)
    atom_map = get_atom_mapping(mol, symels)
    mtable = build_mult_table(symels)
    return mol, SymText(pg, symels, ctab, class_map, atom_map, mtable, length(symels))
end

"""
A function that returns an index so that irreps can be sorted into canonical order
"""
function irrep_sort_idx(irrep_str)
    bunyon = 0
    # g and ' always go first
    gchk = r"g"
    ppchk = r"''"
    gm = occursin(gchk, irrep_str)
    pm = occursin(ppchk, irrep_str)
    if gm
        bunyon += 0
    elseif pm
        bunyon += 10000
    else
        bunyon += 1000
    end
    bub = irrep_str[1] # the letter
    numbro = r"(\d+)"
    mn = match(numbro, irrep_str)
    if mn != nothing
        bunyon += parse(Int,mn.captures[1])
    end
    if bub == 'A'
        bunyon += 0
    elseif bub == 'B'
        bunyon += 100
    elseif bub == 'E'
        bunyon += 200
    elseif bub == 'T'
        bunyon += 300
    elseif bub == 'G'
        bunyon += 400
    elseif bub == 'H'
        bunyon += 500
    else
        throw(ArgumentError("Oh shit the trout population!"))
    end
    return bunyon
end