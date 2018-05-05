# TriangleCheck
# Returns whether or not all triangle inequality constraints are satisfied
# to within the desired tolerance
function TriangleCheck(D::Matrix{Float64},tol::Float64)
  # Checks whether or not the triangle inequality constraints for matrix D
  # are satisfied to within the given tolerance
  n = size(D,1)
  answer = true;
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

# Evaluate the LambdaCC LP relaxed objective, given a distance matrix D
function LPcc_obj(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Float64},lam::Float64)
  n = size(A,1)
  # assert(issymmetric(D))
  # assert(issymmetric(A))
  numedges = countnz(A)/2
  lccBound = sum((A[j,i]-lam)*D[j,i] for i=1:n-1 for j = i+1:n) + lam*(n*(n-1)/2 - numedges)
  return lccBound
end

function TripleLoop!(E::Matrix{Float64},Enew::Matrix{Float64},D::Matrix{Int8},W::Matrix{Float64},corrections::Array{Float64})

    n = size(D,1)
    numcon = n*(n-1)*(1+(n-2)/2)        # number of constraints
    lamt = 1/numcon
    epsi = 1e-12
    totalLam1 = 0
    index = 1
    violations = 0

    noChangeCount = 0
    Epart = zeros(n,n)

    for i = 1:n-2
        for j = i+1:n-1
            Dij = D[j,i]
            Wij = W[j,i]
            for k = j+1:n
                Dik = D[k,i]
                Djk = D[k,j]
                Wik = W[k,i]
                Wjk = W[k,j]
                denom = Wij*Wjk + Wik*Wij + Wjk*Wik

                # First constraint
                cor = corrections[index]

                # Corretion step
                Cij = E[j,i] + cor*(Wik*Wjk/denom)
                Cik = E[k,i] - cor*(Wij*Wjk/denom)
                Cjk = E[k,j] - cor*(Wik*Wij/denom)

                # Projection step
                b = Dik + Djk - Dij
                mu = (Cij - Cjk - Cik - b)
                if mu > epsi
                    Pij = Cij - mu*(Wik*Wjk/denom)
                    Pik = Cik + mu*(Wij*Wjk/denom)
                    Pjk = Cjk + mu*(Wik*Wij/denom)
                    violations+=1
                    corrections[index] = mu
                else
                    Pij = Cij
                    Pik = Cik
                    Pjk = Cjk
                    corrections[index] = 0
                end
                index = index+1

                # Update Enew step
                Enew[j,i] += lamt*(Pij)
                Enew[k,i] += lamt*(Pik)
                Enew[k,j] += lamt*(Pjk)
                totalLam1 += lamt

                # Second constraint
                cor = corrections[index]

                # Corretion step
                Cij = E[j,i] - cor*(Wik*Wjk/denom)
                Cik = E[k,i] + cor*(Wij*Wjk/denom)
                Cjk = E[k,j] - cor*(Wik*Wij/denom)

                # Projection step
                b = -Dik + Djk + Dij
                mu = (-Cij - Cjk + Cik - b)
                if mu > epsi
                    Pij = Cij + mu*(Wik*Wjk/denom)
                    Pik = Cik - mu*(Wij*Wjk/denom)
                    Pjk = Cjk + mu*(Wik*Wij/denom)
                    violations+=1
                    corrections[index] = mu
                else
                    Pij = Cij
                    Pik = Cik
                    Pjk = Cjk
                    corrections[index] = 0
                end

                index = index+1

                # Update Enew step
                Enew[j,i] += lamt*(Pij)
                Enew[k,i] += lamt*(Pik)
                Enew[k,j] += lamt*(Pjk)
                totalLam1 += lamt
                # third constraint
                cor = corrections[index]

                # Corretion step
                Cij = E[j,i] - cor*(Wik*Wjk/denom)
                Cik = E[k,i] - cor*(Wij*Wjk/denom)
                Cjk = E[k,j] + cor*(Wik*Wij/denom)

                # Projection step
                b = Dik - Djk + Dij
                mu = (-Cij + Cjk - Cik - b)
                if mu > epsi
                    Pij = Cij + mu*(Wik*Wjk/denom)
                    Pik = Cik + mu*(Wij*Wjk/denom)
                    Pjk = Cjk - mu*(Wik*Wij/denom)
                    violations+=1
                    corrections[index] = mu
                else
                    Pij = Cij
                    Pik = Cik
                    Pjk = Cjk
                    corrections[index] = 0
                end

                index = index+1

                # Update Enew step
                Enew[j,i] += lamt*(Pij)
                Enew[k,i] += lamt*(Pik)
                Enew[k,j] += lamt*(Pjk)
                totalLam1 += lamt
            end
        end
    end
#    @show totalLam1
#    @show violations
end

function DoubleLoop!(E::Matrix{Float64},Enew::Matrix{Float64},F::Matrix{Float64},Fnew::Matrix{Float64},W::Matrix{Float64},P::Matrix{Float64},Q::Matrix{Float64})
    epsi = 1e-12
    n = size(E,1)
    numcon = n*(n-1)*(1 + (n-2)/2)        # number of constraints
    lamt = 1/numcon
    totalLam2 = 0

    dvio = 0

    for i = 1:n-1
      for j = i+1:n

        cor = P[j,i]
        Eij = E[j,i] - cor
        Fij = F[j,i] - cor
        delta = -Eij - Fij
        if delta > epsi
            Enew[j,i] += lamt*(Eij + delta/2)
            Fnew[j,i] += lamt*(Fij + delta/2)
            P[j,i] = delta/2
            totalLam2 += lamt
            dvio += 1
        else
            Enew[j,i] += lamt*(Eij)
            Fnew[j,i] += lamt*(Fij)
            P[j,i] = 0.0
            totalLam2 += lamt
        end
    end
    end

    for i = 1:n-1
      for j = i+1:n
        cor = Q[j,i]
        Eij = E[j,i] - cor
        Fij = F[j,i] + cor
        delta = Eij - Fij
        if delta > epsi
          Enew[j,i] += lamt*(Eij - delta/2)
          Fnew[j,i] += lamt*(Fij + delta/2)
          Q[j,i] = -delta/2
          totalLam2 += lamt
          dvio += 1
        else
            Enew[j,i] += lamt*(Eij)
            Fnew[j,i] += lamt*(Fij)
            Q[j,i] = 0.0
            totalLam2 += lamt
        end
      end
    end
  #  @show totalLam2
  #  @show dvio
end

using MAT

function TriangleFixingStub(n)

    mat = matread("KarateA.mat")
    A = mat["A"]

    n = size(A,1)
    Dgraph = zeros(Int8,n,n)
    for i = 1:n-1
        for j = i+1:n
            if A[i,j] < .1
                Dgraph[j,i] = 1
            end
        end
    end
    D = Dgraph

    E = zeros(n,n);
    lam = .5
    W = lam*ones(n,n)
    for i = 1:n-1
        for j = i+1:n
            if A[i,j] > .1
                W[j,i] = 1-lam
            end
        end
    end
    gam = 5; Etol = .05; TriTol = .01
    F = -gam*ones(n,n); P = zeros(n,n); Q = zeros(n,n)
    Fnew = zeros(n,n)
    Enew = zeros(n,n)
    @show Threads.nthreads()

    numtricon = Int(n*(n-1)*(n-2)/2)
    corrections = zeros(Float64,numtricon)

    numcon = Int(n*(n-1)*(1 + (n-2)/2))        # number of constraints
    lamt = 1/numcon

    TripleLoop!(E,Enew,D,W,corrections)

    # We didn't update F for the projections, so we account for that now.
    # For the first large number of constraints, it stays the same
    Fnew = lamt*(n*(n-1)*(n-2)/2)*F

    #@show lamt*(n*(n-1)*(n-2)/2)

    ## Now we handle the constraints of the form Eij - Fij <= 0,
    # -Eij - Fij <= 0.

    P = zeros(n,n)   # correction variables for constraints -Eij-Fij <=0
    Q = zeros(n,n)   # correction variables for constraints Eij-Fij <=0
    DoubleLoop!(E,Enew,F,Fnew,W,P,Q)

   # Now we enter the main loop

   # We are ready to reset E and Enew, F and Fnew
   F = Fnew
   E = Enew

iter = 0
   while true && iter < 5
     iter +=1
       tempE = E[:,:]
       tempF = F[:,:]

       totalLam = 0
       Enew = zeros(n,n)
       tic()
       TripleLoop!(E,Enew,D,W,corrections)
       lasttime = toq()
       @show norm(E-Enew)
       Fnew = lamt*(n*(n-1)*(n-2)/2)*F

       ## Now we handle the constraints of the form Eij - Fij <= 0, -Eij - Fij <= 0.
       DoubleLoop!(E,Enew,F,Fnew,W,P,Q)

       # We are ready to reset E and Enew, F and Fnew
       F = Fnew
       E = Enew
       Fnew = zeros(n,n)
       Enew = zeros(n,n)

       change = norm(vec(E-tempE))
       tricheck = TriangleCheck(E+D,TriTol)
       objective = LPcc_obj(A,E+D,lam)
       ch = round(change,3)
       lt = round(lasttime,2)
       ob = round(objective,3)
       #if iter%1 == 1
         G = E+D
         tr = FullTriangleCheck(G)
        # println("Triangle inequality check: $tr")
       #end
       println("Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Tri $tr")
       if change < Etol && tricheck
         break
       end
   end  # end main loop

end # end TriangleFixingStub

#TriangleFixingStub(parse(Int64,ARGS[1]))
TriangleFixingStub(10)
