# This Julia file constains functions for creating animations
# that illustrate the L2 metric nearness,specifically for dissimilarity
# matrices associated with a sparse graph.
# I.e., if A is the adjacency matrix for a sparse graph, the dissimilarity
# matrix is D = 1-A. This is the choice that leads to a relationship between
# correlation clustering and metric nearness. Here we're just solving L2
# metric nearness, so you can think of this as a quadratic programming relaxation
# of correlation clustering that isn't as tight as the LP relaxation. It's easier
# to visualize this L2 version than the L1 version.

using Plots
using MatrixNetworks
pyplot()

# Given a sparse adjacency matrix A for a graph, return a
# dissimilarity matrix D = 1 - A
function GetDfromA(A::SparseMatrixCSC{Float64,Int64})
  n = size(A,1)
  D = zeros(n,n)
  for i = 1:n-1
      for j = i+1:n
          if A[i,j]== 0
              D[j,i] = 1
              D[i,j] = 1
          end
      end
  end
  return D
end

# Given an n x n matrix, find the constraint violations at each triplet
# (x[i], y[i], z[i], V[i]) give the violation V[i] at the given triplet
function GetViolations(X::Matrix{Float64})
    n = size(X,1)
    # Get violation at each triplet
    x = Vector{Int64}()
    y = Vector{Int64}()
    z = Vector{Int64}()
    V = Vector{Float64}()
    maxvio = 0.0
    @inbounds for i = 1:n-2
        for j = i+1:n-1
          a = X[j,i]
          for k = j+1:n
            b = X[k,i]
            c = X[k,j]
            vio = maximum([a-b-c,b-a-c,c-a-b])
            if vio > maxvio
              maxvio = vio
            end
            if vio > 0
            push!(x,k)
            push!(y,j)
            push!(z,i)
            push!(V,vio)
            end
          end
        end
      end
      return x, y, z, V, Base.round(maxvio,4)
  end


# All functions start by plotting everything before we make any changes
# with Dykstra's algorithm
function FirstFrame(A::SparseMatrixCSC{Float64,Int64},xy::Matrix{Float64},anim::Animation)

  n = size(A,1)

  # Empty plots for a space-filler
  info = plot(leg=false, axis = false,grid = false)
  pass = 0
  i = 0; j = 0; k = 0;
  localvio = 0
  maxvio = 1.0
  annotate!([
            (.6,.8,text("Round = ",15,:right)),
            (.95,.8,text("$pass",15,:right)),
            # (.6,.6,text("(i,j,k) = ",15,:right)),
            # (.95,.6,text("(00,00,00)",15,:right)),
            (.6,.2,text("Max Violation = ",15,:right)),
            (.95,.2,text("$maxvio",15,:right))
            ])
  # Construct the initial dissimilarity matrix
  X = GetDfromA(A)

  # Plot the enties of X as a heatmap
  xp = plot(heatmap(X,color=:grays_r),yflip=true)

  # Get violation at each triplet
  x,y,z,V,maxvio = GetViolations(X)

  # Show the constraint violation at each triplet as a scatterplot
  C = scatter3d(x,y,z,alpha=sqrt.(V),zcolor=V,
    markerstrokewidth=0,xlabel = "",ylabel = "",zlabel = "",xticks=[],yticks=[],zticks=[], xlims = (0,n), ylims = (0,n), zlims = (0,n),
    framestyle=:semi,  label="", color=:blues,markersize = 5, colorbar = false)

  # Plot of the graph
  g = graphplot(A,xy)

  # Put all plots together
  #plot(xp,C,mt,g,mt2,layout = @layout [a{.4w} b; c d{.5w} e])
  plot(xp, C, g, info, layout = @layout [a{.4w} b; c d])
  frame(anim)

end

# Perform and plot all steps of the L2 Metric Nearness Problem.
#
# Itnum is the number of iterations
# X is the dissimilarity matrix that we are turning into a metric matrix
# triplet_corrections is the dictionary where we store nonzero dual variables for Dykstra's method
# A is the sparse adjacency matrix for the graph we're performing computations on
# xy stores the coordinates for drawing the graph A
function AllStepsMetricL2!(Itnum::Int64,A::SparseMatrixCSC{Float64,Int64},xy::Matrix{Float64})

    n = size(A,1)

    # Begin the animation
    anim = Animation()

    # Construct the initial dissimilarity matrix, set X_0 = D essentially
    X = GetDfromA(A)

    # Save the first plot where we haven't yet changed anything
    FirstFrame(A,xy,anim)

    triplet_corrections = Dict{Int64,Float64}()

    # Now perform Dykstra's method and plot intermediate steps
    for iters = 1:Itnum
    Xtemp = X[:,:]

    @inbounds for i = 1:n-2
            for j = i+1:n-1
                for k = j+1:n
                    # i,j,k here satisfies i < j < k
                    println("$i $j $k")

                    updated1, lvio1 = ABCproj!(X,triplet_corrections,i,j,n,(i-1)*n^2+(j-1)*n + k,k,i,k,j)
                    updated2, lvio2 = ABCproj!(X,triplet_corrections,i,k,n,(i-1)*n^2+(k-1)*n + j,j,i,k,j)
                    updated3, lvio3 = ABCproj!(X,triplet_corrections,j,k,n,(j-1)*n^2+(k-1)*n + i,j,i,k,i)

                    # If anything was updated, update the animation
                    if updated1 || updated2 || updated3
                    g = graphplot(A,xy,i,j,k)
                    xp = plot(heatmap(X,color=:grays_r),yflip=true)
                    cx, cy, cz, V, maxvio = GetViolations(X)
                    info = plot(leg=false, axis = false,grid = false)
                    localvio = maximum([lvio1,lvio2,lvio3])
                    annotate!([
                              (.6,.8,text("Round = ",15,:right)),
                              (.95,.8,text("$iters",15,:right)),
                              (.6,.6,text("(i,j,k) = ",15,:right)),
                              (.95,.6,text("($i,$j,$k)",15,:right)),
                              (.6,.4,text("Local Violation = ",15,:right)),
                              (.95,.4,text("$localvio",15,:right)),
                              (.6,.2,text("Max Violation = ",15,:right)),
                              (.95,.2,text("$maxvio",15,:right))
                              ])
                    C = scatter3d(cx,cy,cz,alpha=sqrt.(V),zcolor=V,
                      markerstrokewidth=0,xlabel = "",ylabel = "",zlabel = "",xticks=[],yticks=[],zticks=[],xlims = (0,n), ylims = (0,n), zlims = (0,n),
                      framestyle=:semi,  label="", color=:blues,markersize = 5,colorbar = false)
                      plot(xp,C,g,info,layout = @layout [a{.4w} b; c d])
                      frame(anim)
                    end
                end
            end
        end
      change = norm(X-Xtemp)
      vio = FullTriangleCheck(X)
      println("Iteration $iters, change = $change,  ConVio = $vio")
    end
    return anim
end


# This starts by plotting every single change at each triplet, and then switches
# to only showing many updates all at once--all updates that contain a single
# node i at a time.
#
# Variables are the same as in the previous function
function HybridMetricLoop!(Itnum::Int64,A::SparseMatrixCSC{Float64,Int64},xy::Matrix{Float64})

  n = size(A,1)

  # Begin the animation
  anim = Animation()

  # Construct the initial dissimilarity matrix, set X_0 = D essentially
  X = GetDfromA(A)

  # Save the first plot where we haven't yet changed anything
  FirstFrame(A,xy,anim)

  triplet_corrections = Dict{Int64,Float64}()

  for iters = 1:Itnum
  Xtemp = X[:,:]

  @inbounds for i = 1:n-2

      if i == 1 && iters == 1
          for j = i+1:n-1
              for k = j+1:n
                  println("$i $j $k")
                  updated1, lvio1 = ABCproj!(X,triplet_corrections,i,j,n,(i-1)*n^2+(j-1)*n + k,k,i,k,j)
                  updated2, lvio2 = ABCproj!(X,triplet_corrections,i,k,n,(i-1)*n^2+(k-1)*n + j,j,i,k,j)
                  updated3, lvio3 = ABCproj!(X,triplet_corrections,j,k,n,(j-1)*n^2+(k-1)*n + i,j,i,k,i)
                  if updated1 || updated2 || updated3
                  g = graphplot(A,xy,i,j,k)
                  xp = plot(heatmap(X,color=:grays_r),yflip=true)
                  cx, cy, cz, V, maxvio = GetViolations(X)
                  info = plot(leg=false, axis = false,grid = false)
                  localvio = maximum([lvio1,lvio2,lvio3])
                  annotate!([
                            (.6,.8,text("Round = ",15,:right)),
                            (.95,.8,text("$iters",15,:right)),
                            (.6,.6,text("(i,j,k) = ",15,:right)),
                            (.95,.6,text("($i,$j,$k)",15,:right)),
                            (.6,.4,text("Local Violation = ",15,:right)),
                            (.95,.4,text("$localvio",15,:right)),
                            (.6,.2,text("Max Violation = ",15,:right)),
                            (.95,.2,text("$maxvio",15,:right))
                            ])
                  C = scatter3d(cx,cy,cz,alpha=sqrt.(V),zcolor=V,
                    markerstrokewidth=0,xlabel = "",ylabel = "",zlabel = "",xticks=[],yticks=[],zticks=[],xlims = (0,n), ylims = (0,n), zlims = (0,n),
                    framestyle=:semi,  label="", color=:blues,markersize = 5,colorbar = false)
                    plot(xp,C,g,info,layout = @layout [a{.4w} b; c d])
                    frame(anim)
                  end
              end
          end

      else
            for j = i+1:n-1
                for k = j+1:n
                    # i,j,k here satisfies i < j < k
                    ABCproj!(X,triplet_corrections,i,j,n,(i-1)*n^2+(j-1)*n + k,k,i,k,j)
                    ABCproj!(X,triplet_corrections,i,k,n,(i-1)*n^2+(k-1)*n + j,j,i,k,j)
                    ABCproj!(X,triplet_corrections,j,k,n,(j-1)*n^2+(k-1)*n + i,j,i,k,i)
                end
            end
      end
        println("$i")
        g = graphplot(A,xy,i,i,i)
        xp = plot(heatmap(X,color=:grays_r),yflip=true)
        cx, cy, cz, V, maxvio = GetViolations(X)
        info = plot(leg=false, axis = false,grid = false)
        annotate!([
                  (.6,.8,text("Round = ",15,:right)),
                  (.95,.8,text("$iters",15,:right)),
                  (.6,.6,text("i = ",15,:right)),
                  (.95,.6,text("$i",15,:right)),
                  (.6,.2,text("Max Violation = ",15,:right)),
                  (.95,.2,text("$maxvio",15,:right))
                  ])
        C = scatter3d(cx,cy,cz,alpha=sqrt.(V),zcolor=V,
          markerstrokewidth=0,xlabel = "",ylabel = "",zlabel = "",xticks=[],yticks=[],zticks=[],xlims = (0,n), ylims = (0,n), zlims = (0,n),
          framestyle=:semi,  label="", color=:blues,markersize = 5,colorbar = false)
          plot(xp,C,g,info,layout = @layout [a{.4w} b; c d])
          frame(anim)
    end
    change = norm(X-Xtemp)
    vio = FullTriangleCheck(X)
    println("$change  $vio")

  end
  return anim
end

# This only plots a picture every pass through the full set of constraints
# It allows you to visualize how much change happens from iteration to iteration,
# but doesn't give you a close look at how changes happen at each triplet
function RoundMetricLoopAnimate!(Itnum::Int64,A::SparseMatrixCSC{Float64,Int64},xy::Matrix{Float64})
  n = size(A,1)

  # Plot stating information about the round and max violation
  info = plot(leg=false, axis = false,grid = false)
  # Begin the animation
  anim = Animation()

  # Construct the initial dissimilarity matrix, set X_0 = D essentially
  X = GetDfromA(A)

  # Save the first plot where we haven't yet changed anything
  FirstFrame(A,xy,anim)

  triplet_corrections = Dict{Int64,Float64}()

  for iters = 1:Itnum
    Xtemp = X[:,:]
    n = size(X,1)
    @inbounds for i = 1:n-2
              for j = i+1:n-1
                  for k = j+1:n
                      # Perform projections for each triplet
                      ABCproj!(X,triplet_corrections,i,j,n,(i-1)*n^2+(j-1)*n + k,k,i,k,j)
                      ABCproj!(X,triplet_corrections,i,k,n,(i-1)*n^2+(k-1)*n + j,j,i,k,j)
                      ABCproj!(X,triplet_corrections,j,k,n,(j-1)*n^2+(k-1)*n + i,j,i,k,i)
                  end
              end
        end
      change = norm(X-Xtemp)
      vio = FullTriangleCheck(X)
      g = graphplot(A,xy)
      xp = plot(heatmap(X,color=:grays_r),yflip=true)
      cx, cy, cz, V, maxvio = GetViolations(X)
      C = scatter3d(cx,cy,cz,alpha=sqrt.(V),zcolor=V,
        markerstrokewidth=0,xlabel = "",ylabel = "",zlabel = "",xticks=[],yticks=[],zticks=[], xlims = (0,n), ylims = (0,n), zlims = (0,n),
        framestyle=:semi,  label="", color=:blues,markersize = 5,colorbar = false)
        info = plot(leg=false, axis = false,grid = false)
        vio = Base.round(vio,4)
        annotate!([
                  (.6,.8,text("Round = ",15,:right)),
                  (.95,.8,text("$iters",15,:right)),
                  (.6,.2,text("Max Violation = ",15,:right)),
                  (.95,.2,text("$maxvio",15,:right))
                  ])
        plot(xp,C,g,info,layout = @layout [a{.4w} b; c d])
        frame(anim)
      println("$change  $vio")
  end
  return anim
end

# Text!
# This only plots a picture every pass through the full set of constraints
# It allows you to visualize how much change happens from iteration to iteration,
# but doesn't give you a close look at how changes happen at each triplet
function RoundMetricLoopAnimateOld!(Itnum::Int64,A::SparseMatrixCSC{Float64,Int64},xy::Matrix{Float64})
  n = size(A,1)

  # Empty plots for a space-filler
  mt = plot(leg=false, axis = false,grid = false)
  mt2 = plot(leg=false, axis = false,grid = false)

  # Begin the animation
  anim = Animation()

  # Construct the initial dissimilarity matrix, set X_0 = D essentially
  X = GetDfromA(A)

  # Save the first plot where we haven't yet changed anything
  FirstFrame(A,xy,anim)

  triplet_corrections = Dict{Int64,Float64}()

  for iters = 1:Itnum
    Xtemp = X[:,:]
    n = size(X,1)
    @inbounds for i = 1:n-2
              for j = i+1:n-1
                  for k = j+1:n
                      # Perform projections for each triplet
                      ABCproj!(X,triplet_corrections,i,j,n,(i-1)*n^2+(j-1)*n + k,k,i,k,j)
                      ABCproj!(X,triplet_corrections,i,k,n,(i-1)*n^2+(k-1)*n + j,j,i,k,j)
                      ABCproj!(X,triplet_corrections,j,k,n,(j-1)*n^2+(k-1)*n + i,j,i,k,i)
                  end
              end
        end
      change = norm(X-Xtemp)
      vio = FullTriangleCheck(X)
      g = graphplot(A,xy)
      xp = plot(heatmap(X,color=:grays_r),yflip=true)
      cx, cy, cz, V = GetViolations(X)
      C = scatter3d(cx,cy,cz,alpha=sqrt.(V),zcolor=V,
        markerstrokewidth=0,xlabel = "",ylabel = "",zlabel = "",xticks=[],yticks=[],zticks=[], xlims = (0,n), ylims = (0,n), zlims = (0,n),
        framestyle=:semi,  label="", color=:blues,markersize = 5,colorbar = false)
        plot(xp,C,mt,g,mt2,layout = @layout [a{.4w} b; c d{.5w} e])
        frame(anim)
      println("$change  $vio")
  end
  return anim
end

# This is code for performing correction and projection steps (and updating
# dual variables) at a single triplet (a,b,c) in the graph.
function ABCproj!(X::Matrix{Float64},triplet_corrections::Dict{Int64,Float64},
  a::Int64,b::Int64,n::Int64,ijkKey::Int64,ACmax::Int64,ACmin::Int64,BCmax::Int64,BCmin::Int64)

  updated = false
  localvio = 0
  corr = pop!(triplet_corrections,ijkKey,0.0)
  Xab = X[b,a]
  Xac = X[ACmax,ACmin]
  Xbc = X[BCmax,BCmin]
  # Correction step

  if corr > 0
      X[b,a] = Xab + corr/3
      X[ACmax,ACmin] = Xac - corr/3
      X[BCmax,BCmin] = Xbc - corr/3

      Xab = X[b,a]
      Xac = X[ACmax,ACmin]
      Xbc = X[BCmax,BCmin]
  end

  delta = Xab - Xac - Xbc
  if delta > 0
      X[b,a] = Xab - delta/3
      X[ACmax,ACmin] = Xac + delta/3
      X[BCmax,BCmin] = Xbc + delta/3
      triplet_corrections[ijkKey] = delta
  end

  if corr + delta > 0
      updated = true
      localvio = Base.round(delta,4)
  end
  return updated, localvio
end

# FullTriangleCheck
# Returns the worst triangle constraint violation in the whole matrix D
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

# Plot a graph. Code taken (and slightly edited) from Huda Nassar's Julia tutorial:
# https://github.com/nassarhuda/JuliaTutorials/blob/master/plotting.ipynb
function graphplot(A::SparseMatrixCSC{Float64,Int64},xy::Matrix{Float64},a::Int64=0,b::Int64=0,c::Int64=0)
  f = plot(leg=false, axis = false,grid = false)
  ei,ej,w = findnz(triu(A))
  lx = [xy[ei,1]';xy[ej,1]';NaN*ones(1,length(ei))]
  ly = [xy[ei,2]';xy[ej,2]';NaN*ones(1,length(ei))]
  for i = 1:length(w)
      plot!(f,lx[:,i],ly[:,i],color = :black, linewidth = 1)
  end
  scatter!(f,xy[:,1],xy[:,2],color = :black)
  if a > 0
      scatter!(f,[xy[a,1]],[xy[a,2]], color = :red,markersize = 5,linewidth = 0)
      scatter!(f,[xy[b,1]],[xy[b,2]], color = :red,markersize = 5,linewidth = 0)
      scatter!(f,[xy[c,1]],[xy[c,2]], color = :red,markersize = 5,linewidth = 0)
  end

  return f
end
