using Test
using InteractiveUtils
using Distributed
using Base.Threads
using Logging
using RecursiveArrayTools
using DifferentialEquations
using WriteVTK
using Formatting

logger = SimpleLogger(stdout, Logging.Debug)
logger2 = global_logger(logger)
const msh = "/Users/spc/Cirilla/Arbitary/dr/fault.msh"
gmsh.clear()
gmsh.model.occ.addRectangle(-50, -50, 0, 100, 100)
meshsize = 5
@addOption begin
    "Mesh.CharacteristicLengthMax", meshsize
    "Mesh.CharacteristicLengthMin", meshsize
end
gmsh.model.occ.synchronize()
gmsh.model.mesh.generate(2)
gmsh.write("mesh.msh")
# gmsh.open(msh)

mi = get_mesh_information(2)
edge2ele = build_edge_element_map(Tri3, mi["eletags"], mi["edgenodes"])
elenbrlist = build_nbr_edge_list(Tri3, mi["eletag2index"], edge2ele)
elenorlist, elelenlist = build_normal_length_list(Tri3, mi["nodecoords"], mi["nodetag2index"], mi["elenodes"])
bdy = Dict(1 => Absorbing, 2 => Absorbing, 3 => Absorbing, 4 => Absorbing)

function apply_boundary_type!(elenbrlist, eletag2index, edge2ele, edgetype, bcdict)
    if edgetype == 1 # 2-point-line
        nnodesperedge = 2
    else
        return
    end
    count_edge = 0
    for (k, bt) in bcdict
        edgenodes = gmsh.model.mesh.getElementEdgeNodes(edgetype, k)
        for i ∈ 1: nnodesperedge: length(edgenodes)
            edgetagtuple = TupleTools.sort(ntuple(x -> edgenodes[i - 1 + x], Val(nnodesperedge)))
            eleinfo = edge2ele[edgetagtuple]
            @assert length(eleinfo) == 1 "Boundary edge should only associate with one element."
            eletag, ithedge = eleinfo[1]
            eleindex = eletag2index[eletag]
            @assert elenbrlist[eleindex][ithedge][1] == -1 "Not a boundary edge."
            if elenbrlist[eleindex][ithedge][2] ≠ 0
                @warn  "Not a placeholder any more."
            end
            elenbrlist[eleindex][ithedge][2] = Int(bt)
        end
        count_edge += length(edgenodes)
    end
end

function check_unset_boundary(elenbrlist)
    for nbrs ∈ elenbrlist
        for nbr ∈ nbrs
            if nbr[1] == -1 && nbr[2] == 0
                return false
            end
        end
    end
    return true
end

check_unset_boundary(elenbrlist)
apply_boundary_type!(elenbrlist, mi["eletag2index"], edge2ele, 1, bdy)
check_unset_boundary(elenbrlist)


function build_kernels(geotype, mi)
    eletags = mi["eletags"]
    jac = mi["jacobians"]
    detJ = mi["detJ"]

    edge2ele = build_edge_element_map(Tri3, mi["eletags"], mi["edgenodes"])
    elenbrlist = build_nbr_edge_list(Tri3, mi["eletag2index"], edge2ele)
    elenorlist, elelenlist = build_normal_length_list(Tri3, mi["nodecoords"], mi["nodetag2index"], mi["elenodes"])
    apply_boundary_type!(elenbrlist, mi["eletag2index"], edge2ele, 1, bdy)
    if !check_unset_boundary(elenbrlist)
        @error "There are unset boundaries."
    end

    mass = mass_matrix(geotype)
    factormass = qr(mass)
    stiff = stiff_matrix(geotype)
    fluxs = flux_matrix(geotype)
    fluxn = flux_matrix_nbr(geotype)
    reflect = free_surface_reflect_matrix(geotype)
    kernels = Vector{Elasticity2D}(undef, length(eletags))

    @inbounds for i ∈ eachindex(eletags)
        A = A_matrix(geotype, λ, μ, ρ)
        B = B_matrix(geotype, λ, μ, ρ)
        A⁺, A⁻ = upwinding_matrix(A)
        J⁻¹ = inv(jac[i]) # 9 entries for each
        Aˣ = Aˣ_matrix(geotype, A, B, J⁻¹)
        Bˣ = Bˣ_matrix(geotype, A, B, J⁻¹)
        TA⁺T⁻¹S = [similar(A) for _ ∈ 1: 3]
        TA⁻T⁻¹S = [similar(A) for _ ∈ 1: 3]
        Fⁿ = Vector{typeof(mass)}(undef, 3)
        for j ∈ 1: 3 # triangle
            T = T_matrix(geotype, elenorlist[i][j][1], elenorlist[i][j][2])
            TA⁺T⁻¹S[j] = T * A⁺ * inv(T) * elelenlist[i][j]
            if elenbrlist[i][j][1] == -1 # boundary edge
                Fⁿ[j] = fluxs[j] # self edge flux
                if elenbrlist[i][j][2] == Int(Absorbing)
                    TA⁻T⁻¹S[j] = zeros(5, 5) # no incoming flux from neighbor
                elseif elenbrlist[i][j][2] == Int(FreeSurface)
                    TA⁻T⁻¹S[j] = T * A⁻ * reflect * inv(T) * elelenlist[i][j] # reflect
                else
                    # TODO
                    @error "Unrecognizable boundary type $(elenbrlist[i][j][2])."
                    @debug i, j
                end
            else # interior edge
                Fⁿ[j] = fluxn[elenbrlist[i][j][2]][j]
                TA⁻T⁻¹S[j] = T * A⁻ * inv(T) * elelenlist[i][j]
            end
        end
        ker = Elasticity2D(
            factormass, # M
            fluxs, # F⁰
            Fⁿ, # Fⁿ
            stiff[1], # Kξ
            stiff[2], # Kη
            TA⁺T⁻¹S, # TA⁺T⁻¹S
            TA⁻T⁻¹S, # TA⁻T⁻¹S
            Aˣ, # Aˣ
            Bˣ, # Bˣ
            detJ[i], # detJ
        )
        kernels[i] = ker
    end
    kernels
end

function local_solve!(kernel::Elasticity2D, du, u, unbr, p, t, rhs)
    uᵀ = u'
    fill!(rhs, 0)
    for i ∈ 1: 3
        rhs .-= kernel.F⁰[i] * (kernel.TA⁺T⁻¹S[i] * uᵀ)'
        rhs .-= kernel.Fⁿ[i] * (kernel.TA⁻T⁻¹S[i] * unbr[i]')'
    end
    rhs .+= kernel.Kξ * (kernel.Aˣ * uᵀ * kernel.detJ)'
    rhs .+= kernel.Kη * (kernel.Bˣ * uᵀ * kernel.detJ)'
    # no source term for now
    @. rhs /= kernel.detJ
    ldiv!(du, kernel.M, rhs)
end

function global_solve!(kers, du, u, unbrs, p, t, rhss, elelenlist)
    @show t
    for i ∈ eachindex(kers)
        for (j, nbr) ∈ enumerate(elenbrlist[i])
            if nbr[1] == -1 # boundary
                unbrs[i][j] = u[i] # replicate itself
            else # interior
                unbrs[i][j] = u[nbr[1]] # get nbr
            end
        end
        local_solve!(kers[i], du[i], u[i], unbrs[i], p, t, rhss[i])
    end
end

function L2project(f::Function, kers, intpoints, ϕs, points, weights) where T
    ϕw = [f(x) for (x, f) ∈ Iterators.product(points, ϕs)] .* weights
    uinit = [Matrix{Float64}(undef, 3, 5) for _ ∈ 1: length(kks)]
    for (i, u) ∈ enumerate(uinit)
        u_exact = [f(p) for p ∈ intpoints[i]]
        u_exact = hcat(u_exact...)
        rhs = (u_exact * ϕw)'
        ldiv!(u, kers[i].M, rhs)
    end
    return uinit
end

λ, μ, ρ = 2.0, 1.0, 1.0
cp, cs = √((λ + 2μ) / ρ), √(μ / ρ)
geotype = Tri3
@time kks = build_kernels(geotype, mi)
n = [1, 1]
Rp2n = [-2μ*n[1]*n[2], 2μ*n[1]*n[2], μ*(n[1]^2-n[2]^2), -n[2]*cs, n[1]*cs]
Rp5n = [λ+2μ*n[1]^2, λ+2μ*n[2]^2, 2μ*n[1]*n[2], -n[1]*cp, -n[2]*cp]

get_jacobian_info!(mi, 2, 3)
points = mi["points"]
intpoints = mi["intpoints"]
weights = mi["weights"]
ϕ = OrthoBasisFunctions(Tri3, 1)

finit(x) = (Rp2n + Rp5n) * sin(2π/25 * (x[1] + x[2]))
uinit = L2project(finit, kks, intpoints, ϕ, points, weights)
uv = VectorOfArray(uinit)
rhss_cache = similar(uv).u
unbrs_cache = [[Matrix{Float64}(undef, 3, 5) for _ ∈ 1: 3] for _ ∈ eachindex(uinit)]

f(du, u, p, t) = global_solve!(kks, du, u, unbrs_cache, p, t, rhss_cache, elenbrlist)
prob = ODEProblem(f, uv, (0, 100*√2))
sol = solve(prob, SSPRK22(), dt=meshsize/cp/100)

centroid, _ = gmsh.model.mesh.getIntegrationPoints(2, "Gauss1")
ϕi = [f(centroid) for f ∈ ϕ]

nodes = mi["nodecoords"]
conn = mi["elenodes"]
eletags = mi["eletags"]
nn = gmsh.model.mesh.getElementProperties(2)[4]
conn = reshape(conn, Int(nn), :)
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, view(conn, :, i)) for i ∈ eachindex(eletags)]
vtk_points = reshape(nodes, 3, :)

file = "/Users/spc/Downloads/testmovie/test.pvd"
fmt = "%0$(ndigits(length(sol.t)))d"
t, u = sol.t, sol.u
paraview_collection(file) do pvd
    for i = 1: length(t)
        vtk_grid(file * sprintf1(fmt, i), vtk_points, cells) do vtk
            uu = [x' * ϕi for x ∈ u[i]]
            vtk
            vtk["sxx", VTKCellData()] = [x[1] for x ∈ uu]
            vtk["syy", VTKCellData()] = [x[2] for x ∈ uu]
            vtk["sxy", VTKCellData()] = [x[3] for x ∈ uu]
            vtk["v", VTKCellData()] = [x[4] for x ∈ uu]
            vtk["w", VTKCellData()] = [x[5] for x ∈ uu]
            collection_add_timestep(pvd, vtk, t[i])
        end
    end
end

sol.u[1][1]




sol.u[500][1]
