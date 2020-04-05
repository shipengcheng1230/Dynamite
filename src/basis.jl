export OrthoBasisFunctions

# This function subjects to change in the future.
function OrthoBasisFunctions(::Type{Tri3}, order::Integer=1, T=Float64)
    if order == 1
        (x -> oneunit(T),
         x -> 2x[1] + x[2] - oneunit(T),
         x -> 3x[2] - oneunit(T),)
    elseif order == 2
        (x -> oneunit(T),
         x -> 2x[1] + x[2] - oneunit(T),
         x -> 6x[1]^2 + x[2]^2 + 6x[1]*x[2] - 6x[1] - 2x[2] + oneunit(T),
         x -> 3x[2] - oneunit(T),
         x -> 5x[2]^2 + 10x[1]*x[2] - 2x[1] - 6x[2] + oneunit(T),
         x -> 10x[2]^2 - 8x[2] + oneunit(T),)
    else
        @error "Unsupported order: $order > 2."
    end
end
