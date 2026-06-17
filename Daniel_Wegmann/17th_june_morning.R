simulateWF <- function(twoN, f, G){
  p <- numeric(G + 1)
  p[1] <- f
  for(i in 1:G){
    p[i+1] <- rbinom(1, size = twoN, prob = p[i]) / twoN
  }
  return(p)
}
trajectories <- replicate(1000, simulateWF(twoN = 100, f = 0.1, G = 1000))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))

print(paste("Allele was lost in", sum(trajectories[1000,] == 0), "/", ncol(trajectories), "cases."))

trajectories <- replicate(1000, simulateWF(twoN = 100, f = 1/100, G = 1000))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))

simulateWFWithSelection <- function(twoN, f, G, v){
  p <- numeric(G + 1)
  p[1] <- f
  for(i in 1:G){
    # selection
    fA <- p[i]
    fa <- 1-fA
    fPrime <- (v[1]*fA*fA + v[2]*fA*fa)/(v[1]*fA*fA + v[2]*2*fA*fa + v[3]*fa*fa);
    # drift
    p[i+1] <- rbinom(1, size = twoN, prob = fPrime) / twoN
  }
  return(p)
}

s <- 0.01
trajectories <- replicate(100, simulateWFWithSelection(twoN = 10000, f = 0.1, G = 1000, v=c(1,1-s,(1-s)^2)))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))

print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))

s <- 0.05
trajectories <- replicate(100, simulateWFWithSelection(twoN = 10^3, f = 0.1, G = 1000, v=c(1,1,1-s)))
plot(0, type='n', ylim=c(0,1), xlim=c(0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))

print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))

##try with stronger selection
s <- 0.5
#With population size 10000, initital frequency 0.1, and 1000 generations
trajectories <- replicate(100, simulateWFWithSelection(twoN = 10^4, f = 0.1, G = 1000, v=c(1,1-s,(1-s)^2)))
plot(0, type='n', ylim=c(0,1), xlim=c(
  0, nrow(trajectories)))
invisible(apply(trajectories, 2, lines, type='l'))
print(paste("Allele was fixed in", sum(trajectories[1000,] == 1), "/", ncol(trajectories), "cases."))

