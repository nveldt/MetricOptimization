# A Helper Library so that I can keep my other .jl files a little cleaner

# Builds the weights matrix and desired distance matrix D
# for a LambdaCC LP relaxation problem
function LamCC_DandW(A::SparseMatrixCSC{Float64,Int64})
  D = zeros(Float64,n,n)
  W = (1-lam)*ones(n,n)
  for i = 1:n-1
      for j = i+1:n
          if A[i,j] < .1
              D[j,i] = 1.0
              W[j,i] = lam
          end
      end
  end
  return D, W
end

# TriangleCheck
# Returns whether or not all triangle inequality constraints are satisfied
# to within the desired tolerance
function TriangleCheck(D::Matrix{Float64},tol::Float64)
  # Checks whether or not the triangle inequality constraints for matrix D
  # are satisfied to within the given tolerance
  n = size(D,1)
  @inbounds for i = 1:n-2
    for j = i+1:n-1
      for k = j+1:n
        a = D[j,i]
        b = D[k,i]
        c = D[k,j]

        if a - b - c > tol || b - c - a > tol || c - a - b > tol
          return false
        end

      end
    end
  end
  return true
end

# Check if the double loop constraints are satisfied to within tolerance
function DoubleCheck(E::Matrix{Float64},F::Matrix{Float64},tol::Float64)

n = size(E,1)
for i = 1:n-1
  for j = i+1:n
    eij = E[j,i]
    fij = F[j,i]
    if eij - fij > tol || -eij - fij > tol
      return false
    end
  end
end
return true

end

# FullTriangleCheck
# Returns the worst triangle violation in the whole matrix
function FullTriangleCheck(D::Matrix{Float64})
  n = size(D,1)
  maxi = 0.0
  @inbounds for i = 1:n-2
    for j = i+1:n-1
      a = D[j,i]
      for k = j+1:n
        b = D[k,i]
        c = D[k,j]

        #vio = maximum([a-b-c,b-a-c,c-a-b])
        if a-b > maxi && a-c > maxi && a-b-c > maxi
          maxi = a-b-c
        end
        if b - a > maxi && b-c > maxi && b-c-a > maxi
          maxi = b-c-a
        end
        if c-a > maxi && c-c > maxi && c-a-b > maxi
          maxi = c-a-b
        end

      end
    end
  end
  maxi
end

function FullConstraintCheck(X::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64})

  tri = FullTriangleCheck(X)
  n = size(X,1)
  maxi = 0.0

  for i = 1:n-1
    for j = i+1:n
      eij = E[j,i]
      fij = F[j,i]
      if eij - fij > maxi
        maxi = eij - fij
      end
      if -eij - fij > maxi
        maxi = -eij - fij
      end
    end
  end
  doublecheck = maxi

  return round(max(tri,doublecheck),4), round(tri,4)
end


# Evaluate the LambdaCC LP relaxed objective, given a distance matrix D
function LPcc_obj(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Float64},lam::Float64)
  n = size(A,1)
  # assert(issymmetric(D))
  # assert(issymmetric(A))
  numedges = countnz(A)/2
  lccBound = sum((A[j,i]-lam)*D[j,i] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
  return lccBound
end

function LPobj(F::Matrix{Float64},W::Matrix{Float64})
  n = size(F,1)
  obj = 0.0
  for i = 1:n-1
    for j = i+1:n
      obj += F[j,i]*W[j,i]
    end
  end
  return obj
end

function Wnorm(W::Matrix{Float64},E::Matrix{Float64},F::Matrix{Float64})

n = size(W,1)
out = 0.0
for i = 1:n-1
  for j = i+1:n
    out += (E[j,i]^2 + F[j,i]^2)*W[j,i]
  end
end
return out

end

function LR_obj(A::SparseMatrixCSC{Float64,Int64},X::Matrix{Float64})
  n = size(A,1)
  lr = sum((A[j,i]*X[j,i] for i=1:n-1 for j = i+1:n))
  return lr
end
