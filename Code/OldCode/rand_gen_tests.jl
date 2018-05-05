numtimes = 100000

tic()
RR = rand(1:n,numtimes,3)
for t = 1:numtimes
    i,j,k = sort(RR[t,:])
end
toc()

tic()
RR = rand(1:n,numtimes,3)
for t = 1:numtimes
    a,b,c = RR[t,:]
    i,j,k = sort([a,b,c])
end
toc()

tic()
for i = 1:numtimes
    a,b,c = rand(1:n,1,3)
    i,j,k = sort([a,b,c])
end
toc()
