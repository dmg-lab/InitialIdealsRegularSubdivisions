
function subsets_lex(S::Vector{Int64}, k::Int64)
    return sort!(subsets(S,k))
end

function subsets_revlex(S::Vector{Int64}, k::Int64)
    Sk = sort!([reverse!(a) for a in subsets(S,k)])
    return [reverse!(a) for a in Sk]
end

function lineality_space_H_rep(I)
    Grob = groebner_basis(I, complete_reduction = true)
    W_perp = Vector{Vector{Int64}}()
    for f in Grob
        exponents_f = collect(exponents(f))
        for i in 1:length(exponents_f) - 1, j in i+1:length(exponents_f)
            push!(W_perp, exponents_f[i] - exponents_f[j])
        end
    end
    n, W_perp = rref(matrix(QQ, W_perp))
    return W_perp[1:n,:]
end

function lineality_space_V_rep(I)
    H_rep = lineality_space_H_rep(I)
    V_rep = transpose(kernel(H_rep, side = :right))
    n, V_rep = rref(V_rep)
    return V_rep
end

function point_configuration(I)
    M = lineality_space_V_rep(I)
    A = [M[:,i] for i in (1:ncols(M))]
    return A
end

function stratum(I, B)
    R = base_ring(I)
    x = gens(R)
    xNB = [x[i] for i in 1:length(x) if !(x[i] in B)]
    return eliminate(I + ideal(R, xNB), xNB)
end

function ideals_of_max_cells(I, w, Delta=point_configuration(I))
    x = gens(base_ring(I))
    subd = subdivision_of_points(Delta, w)
    return [stratum(I, x[B]) for B in maximal_cells(subd)]
end

function ideal_w(I, w, Delta=point_configuration(I))
    return sum(ideals_of_max_cells(I, w, Delta))
end

function ideal_up_w(I, w, Delta=point_configuration(I))
    R = base_ring(I)
    x = gens(R)
    subd = subdivision_of_points(Delta, w)
    Iw=ideal(R, [R(1)])
    for C in maximal_cells(subd)
        NV = [gens(R)[i] for i in setdiff(1:length(Delta), C)]
        SC, x = polynomial_ring(QQ, "x" => sort(C))
        f = hom(SC, R, [gens(R)[i] for i in sort(C)])
        CI = preimage(f, I)
        G = ([R(f(g)) for g in gens(CI)])
        append!(G,NV)
        Iw_new = ideal(R,G)
        Iw = intersect(Iw,ideal(R,G))
    end
    return Iw
end

nu_t = tropical_semiring_map(QQ, min);

function Omega_fan(I, Delta, Sec, outside = false)

    if typeof(Sec) == PolyhedralFan{QQFieldElem}
        rays_Delta = rays_modulo_lineality(Sec)[1];
        cones_Delta = cones(Sec);
    else
        rays_Delta, cones_Delta = Sec
    end

    n_cones_Delta, n_rays_Delta = size(cones_Delta)
    reps = [sum(cones_Delta[i,:]) == 0 ? Rational[0 for _ in 1 : dim(Sec)] : sum(rays_Delta[j] for j in 1:n_rays_Delta if cones_Delta[i,j])  
                               for i in 1:n_cones_Delta]

    tests = []
    for w in reps
        w = lcm(denominator.(w))*w
        init_w_I = initial(I, nu_t, w)
        I_w = ideal_w(I, w, Delta)
        push!(tests, init_w_I == I_w)
    end
    
    if outside
        outside_Omega = [i for i in 1:length(reps) if !tests[i]]
        return [rays_Delta, cones_Delta[outside_Omega, :]]
    else
        inside_Omega = [i for i in 1:length(reps) if tests[i]];
        if typeof(Sec) == PolyhedralFan{QQFieldElem}
            ls = lineality_space(Sec)
        else 
            ls = lineality_space_V_rep(I)
        end
        return polyhedral_fan(cones_Delta[inside_Omega, :], rays_Delta, ls)
    end   
end

function Omega_star_fan(I, Delta, Sec, outside = false)

    if typeof(Sec) == PolyhedralFan{QQFieldElem}
        rays_Delta = rays_modulo_lineality(Sec)[1];
        cones_Delta = cones(Sec);
    else
        rays_Delta, cones_Delta = Sec
    end

    n_cones_Delta, n_rays_Delta = size(cones_Delta)
    reps = [sum(cones_Delta[i,:]) == 0 ? Rational[0 for _ in 1 : dim(Sec)] : sum(rays_Delta[j] for j in 1:n_rays_Delta if cones_Delta[i,j])  
                               for i in 1:n_cones_Delta]

    tests = []
    for w in reps
        w = lcm(denominator.(w))*w
        init_w_I = initial(I, nu_t, -w)
        I_w = ideal_up_w(I, w, Delta)
        push!(tests, init_w_I == I_w)
    end
    
    if outside
        outside_Omega = [i for i in 1:length(reps) if !tests[i]]
        return [rays_Delta, cones_Delta[outside_Omega, :]]
    else
        inside_Omega = [i for i in 1:length(reps) if tests[i]];
        return polyhedral_fan(incidence_matrix(cones_Delta[inside_Omega, :]), rays_Delta, lineality_space(Sec))
    end
    
end