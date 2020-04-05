export build_edge_element_map, build_normal_length_list, build_nbr_edge_list, get_mesh_information

# get element neighbors information

"""
Return <element_index><edge_index>[nbr_element_index, nbr_edge_index].
The order of "element_index" and "edge_index" are the same as returned by gmsh API.
"""
function build_nbr_edge_list(::Type{Tri3}, eletag2index, edge2ele)
    elenbrlist = [[Vector{Int}(undef, 2) for _ ∈ 1: 3] for _ ∈ 1: length(keys(eletag2index))]

    @inbounds for (edge, t) in edge2ele
        @assert length(t) ≤ 2 "Each edge has no more than 2 neighbors."
        if length(t) == 1 # boundary edge
            eleindex = eletag2index[t[1].eletag]
            elenbrlist[eleindex][t[1].ithedge][1] = Int(-1) # denote no neighbor element on this edge
            elenbrlist[eleindex][t[1].ithedge][2] = Int(0) # placeholder for boundary type
        else length(t) == 2 # interior edge
            eleindex1, eleindex2 = eletag2index[t[1].eletag], eletag2index[t[2].eletag]
            elenbrlist[eleindex1][t[1].ithedge][1] = eleindex2
            elenbrlist[eleindex1][t[1].ithedge][2] = t[2].ithedge
            elenbrlist[eleindex2][t[2].ithedge][1] = eleindex1
            elenbrlist[eleindex2][t[2].ithedge][2] = t[1].ithedge
        end
    end

    return elenbrlist
end

build_edge_element_map(::Type{Tri3}, eletags, edgenodes) = build_edge_element_map(eletags, edgenodes, 2, 3)

function build_edge_element_map(eletags, edgenodes, nnodesperedge, nedges)
    nele = length(eletags)
    nedgenodes = nnodesperedge * nedges
    @assert length(edgenodes) == nele * nedgenodes

    dummy = [(eletag=UInt(0), ithedge=0)]
    edge2ele = Dict{Tuple{UInt, UInt}, typeof(dummy)}()

    let ithedge = 1
        @inbounds for i in 1: nnodesperedge: length(edgenodes)
            edge = TupleTools.sort(ntuple(x -> edgenodes[i - 1 + x], Val(nnodesperedge)))
            t = eletags[i ÷ nedgenodes + Int(i % nedgenodes ≠ 0)]
            if edge ∉ keys(edge2ele)
                edge2ele[edge] = [(eletag=t, ithedge=ithedge)]
            else
                push!(edge2ele[edge], (eletag=t, ithedge=ithedge))
            end
            ithedge = ithedge % nedges + 1
        end
    end
    edge2ele
end

# get element outward normal vector and edge length

function build_normal_length_list(::Type{Tri3}, nodecoords, nodetag2index, elenodes)
    nele = length(elenodes) ÷ 3
    dummy = [zeros(Float64, 3) for _ in 1: 3]
    elenormallist = Vector{typeof(dummy)}(undef, nele)
    dummy = zeros(Float64, 3)
    elelengthlist = Vector{typeof(dummy)}(undef, nele)

    @inbounds for i in 1: nele
        nd = (nodetag2index[k] for k in elenodes[3 * (i-1) + 1: 3 * i])
        p1, p2, p3 = [nodecoords[3x-2: 3x] for x in nd]

        p12 = p2 - p1
        p23 = p3 - p2
        p31 = p1 - p3

        p12n = normalize([-p12[2], p12[1], 0.0]) # gmsh return Float64 coordinates
        p23n = normalize([-p23[2], p23[1], 0.0])
        p31n = normalize([-p31[2], p31[1], 0.0])

        p12n ⋅ p31 < 0 && (p12n *= -1)
        p23n ⋅ p12 < 0 && (p23n *= -1)
        p31n ⋅ p23 < 0 && (p31n *= -1)

        elenormallist[i] = [p12n, p23n, p31n]
        elelengthlist[i] = [norm(p12), norm(p23), norm(p31)]
    end

    return elenormallist, elelengthlist
end

function get_mesh_information(eletype::Int64, intorder::Integer=1)
    mi = Dict{String, Any}()
    get_node_info!(mi)
    get_element_info!(mi, eletype)
    get_edge_node_info!(mi, eletype)
    get_jacobian_info!(mi, eletype, intorder)
    return mi
end

function get_node_info!(mi::AbstractDict)
    nodetags, nodecoords, _ = gmsh.model.mesh.getNodes() # no parallel
    nodetag2index = Dict(nodetags[i] => i for i in eachindex(nodetags))
    mi["nodetags"] = nodetags
    mi["nodecoords"] = nodecoords
    mi["nodetag2index"] = nodetag2index
    return nodetags, nodecoords, nodetag2index
end

function get_element_info!(mi::AbstractDict, eletype::Int64)
    eletags, elenodes = gmsh.model.mesh.getElementsByType(eletype)
    eletag2index = Dict(eletags[i] => i for i in eachindex(eletags))
    mi["eletags"] = eletags
    mi["elenodes"] = elenodes
    mi["eletag2index"] = eletag2index
    return eletags, elenodes, eletag2index
end

function get_edge_node_info!(mi::AbstractDict, eletype::Int64)
    edgenodes = gmsh.model.mesh.getElementEdgeNodes(eletype, -1, 1) # only primary nodes
    mi["edgenodes"] = edgenodes
    edgenodes
end

function get_jacobian_info!(mi::AbstractDict, eletype::Int64, intorder::Integer=1)
    points, weights = gmsh.model.mesh.getIntegrationPoints(eletype, "Gauss" * string(intorder))
    jacobians, determinants, intpoints = gmsh.model.mesh.getJacobians(eletype, points)
    np = length(weights)
    elestride = 1: np: length(determinants)
    _detJ = determinants[elestride]
    _jac = [reshape(jacobians[9i-8: 9i], 3, 3) for i ∈ elestride]
    _intp = [intpoints[3np * (i-1) + 1: 3np * i] for i ∈ 1: length(elestride)]
    _intp2 = [[x[3i-2: 3i] for i ∈ 1: np] for x ∈ _intp]
    _points = [points[3i-2: 3i] for i ∈ 1: np]

    mi["jacobians"] = _jac
    mi["detJ"] = _detJ
    mi["intpoints"] = _intp2
    mi["weights"] = weights
    mi["points"] = _points
    return _jac, _detJ, _intp2, _points, weights
end
