using Test

function calculate_matrix(geo, order)
    # maybe sufficient using Gauss 2order
    points, weights = gmsh.model.mesh.getIntegrationPoints(element_code(geo), "Gauss" * string(2order))

    u, v, w = points[1: 3: end], points[2: 3: end], points[3: 3: end]
    coords = zip(u, v, w)
    bf = OrthoBasisFunctions(geo, order)
    basispairs = Iterators.product(bf, bf)

    # mass matrix
    M = [[f₁(q) * f₂(q) for q in coords] ⋅ weights for (f₁, f₂) in basispairs]

    # stiffness matrix
    Kξ = [[ReverseDiff.gradient(f₁, collect(q))[1] * f₂(q) for q in coords] ⋅ weights for (f₁, f₂) in basispairs]
    Kη = [[ReverseDiff.gradient(f₁, collect(q))[2] * f₂(q) for q in coords] ⋅ weights for (f₁, f₂) in basispairs]
    Kζ = [[ReverseDiff.gradient(f₁, collect(q))[3] * f₂(q) for q in coords] ⋅ weights for (f₁, f₂) in basispairs]

    # flux matrix
    points, weights = gmsh.model.mesh.getIntegrationPoints(element_edge_code(geo), "Gauss" * string(np))
    u, v, w = points[1: 3: end], points[2: 3: end], points[3: 3: end]
    coords = zip(u, v, w)

    paramaps = paramatricalmap(Tri3)
    paramaps_nbr = paramatricalmap_nbr(Tri3)
    paramap_pair = Iterators.product(paramaps, paramaps_nbr)

    # within the element
    F0 = [[[f₁(pm(q)) * f₂(pm(q)) for q in coords] ⋅ weights for (f₁, f₂) in basispairs] for pm in paramaps]
    # across the interior boundary
    Fn = [[[f₁(pmn(q)) * f₂(pms(q)) for q in coords] ⋅ weights for (f₁, f₂) in basispairs] for (pms, pmn) in paramap_pair]

    return (M, Kξ, Kη, Kζ, F0, Fn)
end

@testset "Tri3 Order 1" begin
    M, Kξ, Kη, Kζ, F⁰, Fⁿ = calculate_matrix(Tri3, 1)
    @test M ≈ diagm([1/2, 1/12, 1/4]) atol=1e-8
    @test Kξ ≈ [0 0 0; 1 0 0; 0 0 0] atol=1e-8
    @test Kη ≈ [0 0 0; 1/2 0 0; 3/2 0 0] atol=1e-8
    @test F⁰[1] ./2 ≈ [1 0 -1; 0 1/3 0; -1 0 1] atol=1e-8
    @test F⁰[2] ./2 ≈ [1 1/2 1/2; 1/2 1/3 0; 1/2 0 1] atol=1e-8
    @test F⁰[3] ./2 ≈ [1 -1/2 1/2; -1/2 1/3 0; 1/2 0 1] atol=1e-8
    @test Fⁿ[1,1] ./2 ≈ [1 0 -1; 0 -1/3 0; -1 0 1] atol=1e-8
    @test Fⁿ[2,1] ./2 ≈ [1 1/2 1/2; 0 1/6 -1/2; -1 -1/2 -1/2] atol=1e-8
    @test Fⁿ[3,1] ./2 ≈ [1 -1/2 1/2; 0 1/6 1/2; -1 1/2 -1/2] atol=1e-8
    @test Fⁿ[1,2] ./2 ≈ [1 0 -1; 1/2 1/6 -1/2; 1/2 -1/2 -1/2] atol=1e-8
    @test Fⁿ[2,2] ./2 ≈ [1 1/2 1/2; 1/2 1/6 1/2; 1/2 1/2 -1/2] atol=1e-8
    @test Fⁿ[3,2] ./2 ≈ [1 -1/2 1/2; 1/2 -1/3 0; 1/2 0 1] atol=1e-8
    @test Fⁿ[1,3] ./2 ≈ [1 0 -1; -1/2 1/6 1/2; 1/2 1/2 -1/2] atol=1e-8
    @test Fⁿ[2,3] ./2 ≈ [1 1/2 1/2; -1/2 -1/3 0; 1/2 0 1] atol=1e-8
    @test Fⁿ[3,3] ./2 ≈ [1 -1/2 1/2; -1/2 1/6 -1/2; 1/2 -1/2 -1/2] atol=1e-8
end

@testset "Tri3 Order 2" begin
    M, Kξ, Kη, Kζ, F⁰, Fⁿ = calculate_matrix(Tri3, 2)
    @test M ≈ diagm([1/2, 1/12, 1/30, 1/4, 1/18, 1/6]) atol=1e-8
    @test Kξ ≈ [0 0 0 0 0 0;
                1 0 0 0 0 0;
                0 1/2 0 0 0 0;
                0 0 0 0 0 0;
                2/3 0 0 5/6 0 0;
                0 0 0 0 0 0] atol=1e-8
    @test Kη ≈ [0 0 0 0 0 0;
                1/2 0 0 0 0 0;
                1/3 1/4 0 -1/12 0 0;
                3/2 0 0 0 0 0;
                1/3 5/12 0 5/12 0 0;
                -2/3 0 0 5/3 0 0] atol=1e-8
    @test F⁰[1] ./2 ≈ [1 0 0 -1 0 1;
                       0 1/3 0 0 -1/3 0;
                       0 0 1/5 0 0 0;
                       -1 0 0 1 0 -1;
                       0 -1/3 0 0 1/3 0;
                       1 0 0 -1 0 1] atol=1e-8
    @test F⁰[2] ./2 ≈ [1 1/2 1/3 1/2 1/3 1/3;
                       1/2 1/3 1/4 0 1/12 0;
                       1/3 1/4 1/5 -1/12 0 0;
                       1/2 0 -1/12 1 5/12 2/3;
                       1/3 1/12 0 5/12 1/3 0;
                       1/3 0 0 2/3 0 1] atol=1e-8
    @test F⁰[3] ./2 ≈ [1 -1/2 1/3 1/2 -1/3 1/3;
                       -1/2 1/3 -1/4 0 1/12 0;
                       1/3 -1/4 1/5 -1/12 0 0;
                       1/2 0 -1/12 1 -5/12 2/3;
                       -1/3 1/12 0 -5/12 1/3 0;
                       1/3 0 0 2/3 0 1] atol=1e-8
    @test Fⁿ[1,1] ./2 ≈ [1 0 0 -1 0 1;
                         0 -1/3 0 0 1/3 0;
                         0 0 1/5 0 0 0;
                         -1 0 0 1 0 -1;
                         0 1/3 0 0 -1/3 0;
                         1 0 0 -1 0 1] atol=1e-8
    @test Fⁿ[2,1] ./2 ≈ [1 1/2 1/3 1/2 1/3 1/3;
                         0 1/6 1/6 -1/2 -1/6 -1/3;
                         0 0 1/30 0 -1/6 1/3;
                         -1 -1/2 -1/3 -1/2 -1/3 -1/3;
                         0 -1/6 -1/6 1/2 1/6 1/3;
                         1 1/2 1/3 1/2 1/3 1/3] atol=1e-8
    @test Fⁿ[3,1] ./2 ≈ [1 -1/2 1/3 1/2 -1/3 1/3;
                         0 1/6 -1/6 1/2 -1/6 1/3;
                         0 0 1/30 0 1/6 1/3;
                         -1 1/2 -1/3 -1/2 1/3 -1/3;
                         0 -1/6 1/6 -1/2 1/6 -1/3;
                         1 -1/2 1/3 1/2 -1/3 1/3] atol=1e-8
    @test Fⁿ[1,2] ./2 ≈ [1 0 0 -1 0 1;
                         1/2 1/6 0 -1/2 -1/6 1/2;
                         1/3 1/6 1/30 -1/3 -1/6 1/3;
                         1/2 -1/2 0 -1/2 1/2 1/2;
                         1/3 -1/6 -1/6 -1/3 1/6 1/3;
                         1/3 -1/3 1/3 -1/3 1/3 1/3] atol=1e-8
    @test Fⁿ[2,2] ./2 ≈ [1 1/2 1/3 1/2 1/3 1/3;
                         1/2 1/6 1/12 1/2 1/4 1/3;
                         1/3 1/12 1/30 5/12 1/6 1/3;
                         1/2 1/2 5/12 -1/2 -1/12 -1/3;
                         1/3 1/4 1/6 -1/12 1/6 -1/3;
                         1/3 1/3 1/3 -1/3 -1/3 1/3] atol=1e-8
    @test Fⁿ[3,2] ./2 ≈ [1 -1/2 1/3 1/2 -1/3 1/3;
                         1/2 -1/3 1/4 0 -1/12 0;
                         1/3 -1/4 1/5 -1/12 0 0;
                         1/2 0 -1/12 1 -5/12 2/3;
                         1/3 -1/12 0 5/12 -1/3 0;
                         1/3 0 0 2/3 0 1] atol=1e-8
    @test Fⁿ[1,3] ./2 ≈ [1 0 0 -1 0 1;
                         -1/2 1/6 0 1/2 -1/6 -1/2;
                         1/3 -1/6 1/30 -1/3 1/6 1/3;
                         1/2 1/2 0 -1/2 -1/2 1/2;
                         -1/3 -1/6 1/6 1/3 1/6 -1/3;
                         1/3 1/3 1/3 -1/3 -1/3 1/3] atol=1e-8
    @test Fⁿ[2,3] ./2 ≈ [1 1/2 1/3 1/2 1/3 1/3;
                         -1/2 -1/3 -1/4 0 -1/12 0;
                         1/3 1/4 1/5 -1/12 0 0;
                         1/2 0 -1/12 1 5/12 2/3;
                         -1/3 -1/12 0 -5/12 -1/3 0;
                         1/3 0 0 2/3 0 1] atol=1e-8
    @test Fⁿ[3,3] ./2 ≈ [1 -1/2 1/3 1/2 -1/3 1/3;
                         -1/2 1/6 -1/12 -1/2 1/4 -1/3;
                         1/3 -1/12 1/30 5/12 -1/6 1/3;
                         1/2 -1/2 5/12 -1/2 1/12 -1/3;
                         -1/3 1/4 -1/6 1/12 1/6 1/3;
                         1/3 -1/3 1/3 -1/3 1/3 1/3] atol=1e-8
end
