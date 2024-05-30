module Queue
export FIFO

mutable struct FIFO{T}
    data::Vector{T}
    left::Int
    right::Int
    const N::Int
    FIFO{T}(N::Int) where {T} = new{T}(Vector{T}(undef, N), 1, 1, N)
end

function Base.push!(b::FIFO{T}, x::T) where {T}
    b.data[b.right] = x
    b.right += 1
    if b.right > b.N
        b.right = 1
    end
end

function Base.popfirst!(b::FIFO{T}) where {T}
    if isempty(b)
        throw(ArgumentError("FIFO is empty"))
    end
    x = b.data[b.left]
    b.left += 1
    if b.left > b.N
        b.left = 1
    end
    x
end

function Base.empty!(b::FIFO{T}) where {T}
    b.left = 1
    b.right = 1
    b
end

function Base.isempty(b::FIFO{T}) where {T}
    b.left == b.right
end
end