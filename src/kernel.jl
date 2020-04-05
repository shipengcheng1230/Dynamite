export Elasticity2D

@enum BoundaryType Absorbing=1 FreeSurface=2

abstract type AbstractKernel end

struct Elasticity2D <: AbstractKernel
    M
    F⁰ # self
    Fⁿ # nbr
    Kξ
    Kη
    TA⁺T⁻¹S
    TA⁻T⁻¹S
    Aˣ
    Bˣ
    detJ
    # S # source term
    # nbr # neighbor
    # sol # solution
end

solution_shape(::Elasticity2D) = (3, 5)

function A_matrix(::Type{Tri3}, λ, μ, ρ)
    [0.0 0 0 -(λ+2μ) 0;
     0 0 0 -λ 0;
     0 0 0 0 -μ;
     -1/ρ 0 0 0 0;
     0 0 -1/ρ 0 0]
end

function B_matrix(::Type{Tri3}, λ, μ, ρ)
    [0.0 0 0 0 -λ;
     0 0 0 0 -(λ+2μ);
     0 0 0 -μ 0;
     0 0 -1/ρ 0 0;
     0 -1/ρ 0 0 0]
end

function T_matrix!(::Type{Tri3}, dest, nx::T, ny::T) where T
    nx², ny² = nx^2, ny^2
    nxny = nx * ny

    # [nx² ny² -2nxny 0 0;
    #  ny² nx² 2nxny 0 0;
    #  nxny -nxny nx²-ny² 0 0;
    #  0 0 0 nx -ny;
    #  0 0 0 ny nx]

    fill!(dest, zero(T))
    dest[1,1] = nx²
    dest[1,2] = ny²
    dest[1,3] = -2nxny
    dest[2,1] = ny²
    dest[2,2] = nx²
    dest[2,3] = 2nxny
    dest[3,1] = nxny
    dest[3,2] = -nxny
    dest[3,3] = nx²-ny²
    dest[4,4] = nx
    dest[4,5] = -ny
    dest[5,4] = ny
    dest[5,5] = nx
end

function T_matrix(::Type{Tri3}, nx::T, ny::T) where T
    dest = Matrix{T}(undef, 5, 5)
    T_matrix!(Tri3, dest, nx, ny)
    return dest
end

function upwinding_matrix!(A⁺, A⁻, A)
    e = eigen(A)
    # remove allocation in the future
    AA = e.vectors * diagm(abs.(e.values)) * inv(e.vectors)
    A⁺ .= 1/2 * (A + AA)
    A⁻ .= 1/2 * (A - AA)
end

function upwinding_matrix(A)
    A⁺, A⁻ = similar(A), similar(A)
    upwinding_matrix!(A⁺, A⁻, A)
    return A⁺, A⁻
end

function Aˣ_matrix!(::Type{Tri3}, dest, A, B, J⁻¹)
    @. dest = A * J⁻¹[1,1] + B * J⁻¹[1,2]
end

function Aˣ_matrix(::Type{Tri3}, A, B, J⁻¹)
    dest = similar(A)
    Aˣ_matrix!(Tri3, dest, A, B, J⁻¹)
    return dest
end

function Bˣ_matrix!(::Type{Tri3}, dest, A, B, J⁻¹)
    @. dest = A * J⁻¹[2,1] + B * J⁻¹[2,2]
end

function Bˣ_matrix(::Type{Tri3}, A, B, J⁻¹)
    dest = similar(A)
    Bˣ_matrix!(Tri3, dest, A, B, J⁻¹)
    return dest
end

function free_surface_reflect_matrix(::Type{Tri3}, T=Float64)
    diagm([-one(T), one(T), -one(T), one(T), one(T)])
end
