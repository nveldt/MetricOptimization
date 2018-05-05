include("Tricon_Helper.jl")

# Currently this code shouldn't work for weighted graphs, so beware.
# Either way it is impractical, as it takes forever to converge, even if the 3Loop
# isn't super slow on tiny graphs.
function cyclic_triangle_constraints!(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},
    W::Matrix{Float64},F::Matrix{Float64},gamma::Float64)

n = size(A,1)
@inbounds for i = 1:n-2
    for j = i+1:n-1
        Wij = W[j,i]
        for k = j+1:n
            Wik = W[k,i]
            Wjk = W[k,j]
            denom = Wik*Wjk + Wij*Wik + Wij*Wjk
            # i,j,k here satisfies i < j < k
            # There are three ways to order these, each with a different constraint
            ABCproj!(D,E,F,gamma,i,j,Wij,Wik,Wjk,n,denom,k,i,k,j)
            ABCproj!(D,E,F,gamma,i,k,Wik,Wij,Wjk,n,denom,j,i,k,j)
            ABCproj!(D,E,F,gamma,j,k,Wjk,Wij,Wik,n,denom,j,i,k,i)
        end
    end
end
end

function ABCproj!(D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},gamma::Float64,a::Int64,b::Int64,
  wab::Float64,wac::Float64,wbc::Float64,n::Int64,denom::Float64,maxAC::Int64,minAC::Int64,maxBC::Int64,minBC::Int64)

   # Get current iterates
   Eab = E[b,a]
   Eac = E[maxAC,minAC]
   Ebc = E[maxBC,minBC]
   Dab = D[b,a]
   Dac = D[maxAC,minAC]
   Dbc = D[maxBC,minBC]

   delta = Eab+Dab - Eac -Dac - Ebc - Dbc
   # We only change things if the constraint is violated
   if delta > 0

       EabNew = Eab - delta*wbc*wac/denom
       EacNew = Eac + delta*wab*wbc/denom
       EbcNew = Ebc + delta*wab*wac/denom

       v = (delta/denom)^2*(wbc^2*wac^2 + wab^2*wbc^2 + wab^2*wac^2)
       X = (Eab*(-wbc*wac)+Eac*(wab*wbc) + Ebc*(wab*wac))*delta/denom

       # mu is the hard thing to calculate
       uE = 0
       uF = 0
       for i = 1:n-1
         for j = i+1:n
           uE += E[j,i]^2
           uF += (F[j,i]+gamma)^2
         end
       end
       u = uE + uF
       p = u*v - X^2

       if p == 0 && X >= 0
         E[a,b] = Eab - delta*wbc*wac/denom
         E[maxAC,minAC] = Eac + delta*wab*wbc/denom
         E[maxBC,minBC] = Ebc + delta*wab*wac/denom
       end

       if p > 0
         #println("Something went horribly wrong with rho. Triple \n")
         #@show p, delta

       if X*v >= p
         println("Hard reset, didn't expect this \n")
         E = zeros(n,n)
         E[b,a] = -delta*wbc*wac/denom*(1 + X/v)
         E[maxAC,minAC] = delta*wab*wbc/denom*(1 + X/v)
         E[maxBC,minBC] = delta*wab*wac/denom*(1 + X/v)
         F = -gamma*tril(ones(n,n))
       else
         # This is what I expect to happen most often if there is a change
        #println("As expected, triple loop")
         del = X*v/p

         # Update E part and F part simultaneously
         for i = 1:n-1
           for j = i+1:n
            E[j,i] = (1-del)*E[j,i]
            F[j,i] -= del*(F[j,i] + gamma)
           end
         end
         E[b,a] -= u*v/p*(delta*wbc*wac/denom)
         E[maxAC,minAC] += u*v/p*(delta*wab*wbc/denom)
         E[maxBC,minBC] += u*v/p*(delta*wab*wac/denom)
       end
   end
 end
   # At the end of all this, E must be the correct new vector, even if it takes and extra O(n^2) time for this single projection

end

function doubleloop!(E::Matrix{Float64},F::Matrix{Float64},n::Int64,gamma::Float64)

  # Part 2: -E <= F and E <= F checking, the second set of constraints
  for i = 1:n-1
    for j = i+1:n

      Eij = E[j,i]
      Fij = F[j,i]
      delta = -Eij - Fij
      if delta > 0
        EijNew = Eij + delta/2
        FijNew = Fij + delta/2

        X = Eij*(delta/2) + (gamma+Fij)*(delta/2)
        v = delta^2/2

        uE = 0
        uF = 0
        for a = 1:n-1
          for b = a+1:n
            uE += E[b,a]^2
            uF += (F[b,a]+gamma)^2
          end
        end
        u = uE + uF
        p = u*v - X^2

        if p == 0 && X >= 0
          E[j,i] = EijNew
          F[j,i] = FijNew
        end

        if p > 0
          #println("Something went horribly wrong with rho. Double Loop 1 \n")
          #@show p, delta, X, v, u


        if X*v >= p
          println("Hard reset, didn't expect this. Double Loop 1 \n")
          E = zeros(n,n)
          E[j,i] = (1+X/v)*(EijNew - Eij)
          F = -gamma*tril(ones(n,n))
          F[j,i] += (1+X/v)*(FijNew - Fij)
        else
          # This is what I expect to happen most often if there is a change
        #  println("As expected, double loop")
          del = X*v/p

          # Update E part and F part simultaneously
          for a = 1:n-1
            for b = a+1:n
             E[b,a] = (1-del)*E[b,a]
             F[b,a] -= del*(F[b,a] + gamma)
            end
          end
          E[j,i] += u*v/p*(delta/2)
          F[j,i] += u*v/p*(delta/2)
        end
      end
    end
  end
  end

  for i = 1:n-1
    for j = i+1:n
      Eij = E[j,i]
      Fij = F[j,i]
      delta = Eij - Fij


      if delta > 0

        EijNew = Eij - delta/2
        FijNew = Fij + delta/2

        X = -Eij*(delta/2) + (gamma+Fij)*(delta/2)
        v = delta^2/2

        uE = 0
        uF = 0
        for a = 1:n-1
          for b = a+1:n
            uE += E[b,a]^2
            uF += (F[b,a]+gamma)^2
          end
        end
        u = uE + uF
        p = u*v - X^2

        if p == 0 && X >= 0
          E[j,i] = EijNew
          F[j,i] = FijNew
        end

        if p > 0
          #  println("Something went horribly wrong with rho. Double Loop 2\n")
          #  @show p, delta,

          if X*v >= p
            println("Hard reset, didn't expect this \n")
            E = zeros(n,n)
            E[j,i] = (1+X/v)*(EijNew - Eij)
            F = -gamma*tril(ones(n,n))
            F[j,i] += (1+X/v)*(FijNew - Fij)
          else
            # This is what I expect to happen most often if there is a change
  #println("As expected, double loop 2")
            del = X*v/p

            # Update E part and F part simultaneously
            for a = 1:n-1
              for b = a+1:n
               E[b,a] = (1-del)*E[b,a]
               F[b,a] -= del*(F[b,a] + gamma)
              end
            end
            E[j,i] -= u*v/p*(delta/2)
            F[j,i] += u*v/p*(delta/2)
          end
        end
      end

    end
  end
end

function Haugazeaus(A::SparseMatrixCSC{Float64,Int64},D::Matrix{Int8},E::Matrix{Float64},F::Matrix{Float64},W::Matrix{Float64},Etol::Float64,TriTol::Float64,lam::Float64,filename::String,gamma::Float64)
n = size(A,1)
@show filename
  open(filename, "w") do f
          write(f, "Writing output from calling the Triangle Fixing Algorithm\n")
          write(f, "Lambda = $lam, gamma = $gamma, Etol = $Etol, TriTol = $TriTol \n")
  end

# Tolerance for inequalty constraints
epsi = 1e-12
tr = 0


@time cyclic_triangle_constraints!(A,D,E,W,F,gamma)

@time doubleloop!(E,F,n,gamma)

# Done with the first time through the constraints

iter = 0
maxits = 1e10
while true && iter < maxits
    tempE = E[:,:]

    # Start again with the triangle inequality constraints
    tic()
    cyclic_triangle_constraints!(A,D,E,W,F,gamma)
    lasttime = toq()
    doubleloop!(E,F,n,gamma)


    # Print and update of results for this loop
    iter += 1
    change = norm(vec(E-tempE))
    tricheck = TriangleCheck(E+D,TriTol)
    objective = LPcc_obj(A,E+D,lam)
    ch = round(change,3)
    lt = round(lasttime,1)
    ob = round(objective,3)
    if iter%20 == 1
      G = E+D
      tr = FullTriangleCheck(G)
    end
    println("Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr")

    open(filename, "a") do f
      write(f,"Iteration $iter, E changed by $ch, 3Loop: $lt, Obj: $ob, Last TriCheck = $tr\n")
    end

    if change < Etol && tricheck
      break
    end

end #end while loop


return E+D

end
