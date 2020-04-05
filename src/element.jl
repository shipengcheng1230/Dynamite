export Seg2, Tri3

abstract type AbstractElementType end

struct Seg2 <: AbstractElementType end
struct Tri3 <: AbstractElementType end
struct Tri6 <: AbstractElementType end

# gmsh element code
element_code(::Type{Seg2}) = 1
element_code(::Type{Tri3}) = 2
element_code(::Type{Tri6}) = 9

element_edge_code(::Type{Tri3}) = 1

# (ξ, η) -> ((1 + t) / 2, (1 - t) / 2), must be couter-clockwise
paramatricalmap(::Type{Tri3}) = (
    t -> ((1 + t[1]) / 2, zero(t[1])),
    t -> ((1 - t[1]) / 2, (1 + t[1]) / 2),
    t -> (zero(t[1]), (1 - t[1]) / 2),
)

# t′ = -t
paramatricalmap_nbr(::Type{Tri3}) = (
    t -> ((1 - t[1]) / 2, zero(t[1])),
    t -> ((1 + t[1]) / 2, (1 - t[1]) / 2),
    t -> (zero(t[1]), (1 + t[1]) / 2),
)
