# module TSIMatrixOps

# export fillXi1!, identityMat, fillXi2!, ones

"""
    fillXi1!(Xi::BitArray)

Fill in place the `m * n` matrix Xi as to get the first part of a couple forming matrix.

`m` must be equal to `n * (n + 1) / 2`.

# Example
```julia-repl
julia> n = 5
5

julia> A = falses(convert(Integer, n*(n+1)/2), n)
15×5 BitMatrix:
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0

julia> fillXi1!(A)

julia> A
15×5 BitMatrix:
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  0  0
 1  0  0  0  0
 0  1  0  0  0
 0  1  0  0  0
 0  1  0  0  0
 0  1  0  0  0
 0  0  1  0  0
 0  0  1  0  0
 0  0  1  0  0
 0  0  0  1  0
 0  0  0  1  0
 0  0  0  0  1
```

See also [`fillXi2!`](@ref).
"""
function fillXi1!(Xi::BitArray)
    if size(Xi, 1) != size(Xi, 2) * (size(Xi, 2) + 1) / 2
        throw(ArgumentError("`m * n` matrix `Xi`: `m` must be equal to `n * (n + 1) / 2`."))
    end
    n = size(Xi)[2]
    m = 1
    i = 1
    while n > 0
        Xi[m:m+n-1, i] .= 1
        m = m + n
        n = n - 1
        i = i + 1
    end
end

"""Return the identity matrix of size `n`"""
function identityMat(n)
    mat = falses(n, n)
    for i in 1:n
        mat[i,i] = 1
    end

    return mat
end

"""
    fillXi2!(Xi::BitArray)

Fill in place the `m * n` matrix Xi as to get the second part of a couple forming matrix.

`m` must be equal to `n * (n + 1) / 2`.

# Example
```julia-repl
julia> n = 5
5

julia> A = falses(convert(Integer, n*(n+1)/2), n)
15×5 BitMatrix:
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0
 0  0  0  0  0

julia> fillXi2!(A)

julia> A
15×5 BitMatrix:
 1  0  0  0  0
 0  1  0  0  0
 0  0  1  0  0
 0  0  0  1  0
 0  0  0  0  1
 0  1  0  0  0
 0  0  1  0  0
 0  0  0  1  0
 0  0  0  0  1
 0  0  1  0  0
 0  0  0  1  0
 0  0  0  0  1
 0  0  0  1  0
 0  0  0  0  1
 0  0  0  0  1
```

See also [`fillXi1!`](@ref).
"""
function fillXi2!(Xi::BitArray)
    if size(Xi, 1) != size(Xi, 2) * (size(Xi, 2) + 1) / 2
        throw(ArgumentError("`m * n` matrix `Xi`: `m` must be equal to `n * (n + 1) / 2`."))
    end
    n = size(Xi)[2]
    m = 1
    i = 1
    while n > 0
        # Xi[m:m+n-1, i] .= 1
        # Xi[m:m+n-1, i:size(Xi)[2]] = diagm([1 for _ in 1:n])
        # Xi[m:m+n-1, i:size(Xi)[2]] = identityMat(n)
        Xi[m:m+n-1, i:size(Xi)[2]] = LinearAlgebra.I(n)
        m = m + n
        n = n - 1
        i = i + 1
    end
end

"""Return a vertical matrix of dimension 1, type `T` and length `n`"""
function vones(::Type{T}, n::Int) where {T<:Number}
    return Base.ones(T, n, 1)
end



# end
