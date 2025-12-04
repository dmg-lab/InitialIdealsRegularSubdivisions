function parse_tropical_file(path::String, cone_orbits_given = false, negate_rs = false)
    lines = readlines(path)
    rays = Vector{Vector{Int}}()
    cones = Vector{Vector{Int}}()
    cone_orbits = Vector{Vector{Int}}()
    section = ""

    for line in lines
        
        if isempty(line)
            section = ""
            continue
        elseif startswith(line, "RAYS")
            section = "RAYS"
            continue
        elseif startswith(line, "CONES_ORBITS")
            section = "CONES_ORBITS"
            continue
        elseif startswith(line, "CONES")
            section = "CONES"
            continue
        end

        if section == "RAYS"
            line = split(line, "#")[1]
            line = parse.(Int, split(line))
            #rays = [rays; line]
            push!(rays, line)
        elseif cone_orbits_given && section == "CONES_ORBITS"
            line = split(line, "#")[1]
            line = replace(line, ['{', '}']=>"")
            push!(cone_orbits, parse.(Int, split(line)))

        elseif section == "CONES"
            line = split(line, "#")[1]
            line = replace(line, ['{', '}']=>"")
            push!(cones, parse.(Int, split(line)))
        end
    end
    
    if cone_orbits_given
        cone_orbits = pm.to_one_based_indexing(cone_orbits)
    else
        cones = pm.to_one_based_indexing(cones)
    end
    if negate_rs
        rays = -rays
    end
    if cone_orbits_given
        return rays, incidence_matrix(cone_orbits)
    else
        return rays, incidence_matrix(cones)
    end
end

function point_to_string(pt)
    return "(1, " * join([string(a) for a in pt], ", ") * ")"
end

function point_configuration_to_string(pts)
    rs, cs = size(pts)
    return "{" * join([point_to_string(pts[i,:]) for i in 1:rs], ", ") * "}"
end

function pad_variable_numbers(s::AbstractString)
    regex = r"([A-Za-z]+)(\d+)"
    matches = collect(eachmatch(regex, s))
    isempty(matches) && return s

    nums = parse.(Int, getfield.(matches, :captures) .|> x -> x[2])
    width = length(string(maximum(nums)))

    # replacement: match again on the substring to get captures
    result = replace(s, regex => (m -> begin
        mm = match(regex, String(m))         # mm is a RegexMatch with .captures
        prefix = mm.captures[1]
        numstr = mm.captures[2]
        prefix * lpad(numstr, width, '0')
    end))

    return result
end

function ring_to_gfan(n)
    y = join(["y$(i)" for i in 1:n], ",")
    z = pad_variable_numbers(y)
    return "Q[" * z * "]"
end

function polynomial_to_gfan(f_str, var_to_str)
    reduce((acc, kv) -> replace(acc, kv), var_to_str, init=f_str)
end

function ideal_to_gfan(I, make_file = false, file_name = "test.dat")
    if is_zero(I)
        return ""
    end
    I_gens = gens(I)
    
    I_gens_str = [string(f) for f in I_gens]
    R = base_ring(I)
    x = gens(R)
    var_to_str = Dict(string(x[i]) => "p$(lpad(i, 2, '0'))" for i in 1:length(x))

    gfan_ring = ring_to_gfan(length(x))
    gfan_ideal = "{\n" * join(I_gens_str, ",\n") * "\n}"

    gfan_ideal = pad_variable_numbers(gfan_ideal)
    #gfan_ideal = "{\n" * join([polynomial_to_gfan(f_str, var_to_str) 
    #                          for f_str in I_gens_str], ",\n") * "\n}"

    output = gfan_ring * "\n\n" * gfan_ideal
    if make_file
        write(file_name, output)
    end
    return output
end