using Iterators
#simple_partition(n::Int, k::Int) = push!(map(first, partition(1:n,floor(Int,n/k))),n+1)
function simple_partition(n::Int,k::Int)
  rval = ones(Int,k+1)
  m,r = divrem(n,k)
  fill!(rval, m)
  rval[1] = 1
  for i = 2:k+1
    rval[i] += r .> (i-2)
  end
  cumsum!(rval,rval)
end

import Base.LinAlg.At_mul_B, Base.LinAlg.At_mul_B!, Base.A_mul_B, Base.LinAlg.A_mul_B!

"""
A type to encapsulate the multithreaded matrix-vector product operation.

Usage
-----
~~~~
A = sprandn(10^5,10^5,100/10^5); # create a fairly dense sparse matrix.
x = randn(10^5);
M = MultithreadedMatVec(A, simple_partition(A.n, Threads.nthreads()));
y = M.'*x;
norm(y - A.'*x)
~~~~
"""
mutable struct MultithreadedMatVec{T,I} <: AbstractArray{T,2}
  A::SparseMatrixCSC{T,I}
  regions::Vector{Int}

  MultithreadedMatVec(A::SparseMatrixCSC) = MultithreadedMatVec(A,Threads.nthreads())
  MultithreadedMatVec(A::SparseMatrixCSC, k::Int) = MultithreadedMatVec(A, simple_partition(A.n,k))
  function MultithreadedMatVec(A::SparseMatrixCSC{T,I}, regions::Vector{Int}) where {T,I}
    new{T,I}(A, regions)
  end
end

import Base.size
size(M::MultithreadedMatVec,inputs...) = size(M.A,inputs...)

# Julia's built in multiplication operations are called with
# A_mul_B!
# Ac_mul_B!
# At_mul_B!
# which take in
# α::Number, A::SparseMatrixCSC, B::StridedVecOrMat, β::Number, C::StridedVecOrMat
# and compute
# βC += αA B
# βC += αA^* B
# βC += αA^T B
# respectively
# look in base/sparse/linalg.jl
# for their implementations

""" Run the internal loop """
function internal_loop_transmult(C,B,nzv,rv,colptr,i1,i2,α::Number)
  for k=1:size(C,2)
    @inbounds for col=i1:i2
      tmp = zero(eltype(C))
      for j=colptr[col]:(colptr[col+1]-1)
        tmp += nzv[j]*B[rv[j],k]
      end
      C[col,k] += α*tmp
    end
  end
  return
end

"""
we are going to make these work with MuilthreadedMatVec types
"""
function At_mul_B!(α::Number, M::MultithreadedMatVec, B::StridedVecOrMat, β::Number, C::StridedVecOrMat)
  M.A.m == size(B,1) || throw(DimensionMismatch())
  M.A.n == size(C,1) || throw(DimensionMismatch())
  size(B,2) == size(C,2) || throw(DimensionMismatch())
  nzv = M.A.nzval
  rv = M.A.rowval
  colptr = M.A.colptr
  if β != 1
    β != 0 ? scale!(C, β) : fill!(C, zero(eltype(C)))
  end
  # this is the parallel construction
  Threads.@threads for t=1:(length(M.regions)-1)
    internal_loop_transmult(C,B,nzv,rv,colptr,M.regions[t],M.regions[t+1]-1,α)
  end
  C
end

function At_mul_B(M::MultithreadedMatVec{TA,S}, x::StridedVector{Tx}) where {TA,S,Tx}
  T = promote_type(TA,Tx)
  At_mul_B!(one(T), M, x, zero(T), similar(x, T, M.A.n))
end

function internal_loop_mult(C,B,nzv,rv,colptr,i1,i2,α::Number)
  for k=1:size(C,2)
    @inbounds for col=i1:i2
      αxj = α*B[col,k]
      for j=colptr[col]:(colptr[col+1]-1)
        # need to be done atomically... (ideally...)
        Threads.@atomic C[rv[j],k] += nzv[j]*αxj
      end
    end
  end
  return
end

function A_mul_B!(α::Number, M::MultithreadedMatVec, B::StridedVecOrMat, β::Number, C::StridedVecOrMat)
  M.A.n == size(B,1) || throw(DimensionMismatch())
  M.A.m == size(C,1) || throw(DimensionMismatch())
  size(B,2) == size(C,2) || throw(DimensionMismatch())
  nzv = M.A.nzval
  rv = M.A.rowval
  colptr = M.A.colptr
  if β != 1
    β != 0 ? scale!(C, β) : fill!(C, zero(eltype(C)))
  end
  # this is the parallel construction
  Threads.@threads for t=1:(length(M.regions)-1)
    internal_loop_mult(C,B,nzv,rv,colptr,M.regions[t],M.regions[t+1]-1,α)
  end
  C
end

import Base.eltype
eltype(M::MultithreadedMatVec{T,I}) where {T,I} = T

function A_mul_B!(C::StridedVecOrMat, M::MultithreadedMatVec, B::StridedVecOrMat)
  T = promote_type(eltype(M), eltype(B))
  A_mul_B!(one(T), M, B, zero(T), C)
end
