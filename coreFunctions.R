 ###############################################################
####### Functions used for simulating information flow ########
####### in multilayer, heterogeneous biological networks ######
###############################################################

# Regulatory factors
activation <- function(x, fenergy) {
  (1 + x * fenergy) / (1 + x)
}

repression <- function(x) {
  1 / (1 + abs(x))
}

# Transmission functions
regulation <- function(target, inputs, type, capacities) {
  #freg <- 4
  strength <- 0
  regs <- regulators[[type]][[target]]
  if (length(regulators[[type]][[target]]) == 0) {
    return(0)
  }
  #cat(paste0(regs,"_",target,"\n"))
  weights <- as.numeric(capacities[[type]][regs, target])
  return(sum(log2((1 + weights * (inputs[regs])) / (1 + inputs[regs]))))
  # pregs <- which(weights > 0)
  # nregs <- which(weights <= 0)
  # #pweights <- weights[pregs]
  # #nweights <- weights[nregs]
  # activity <- inputs[regs] * weights
  # pact <- activity[pregs]
  # nact <- activity[nregs]
  # if (length(pact) != 0) {
  #   peffects <- sapply(pact,activation,fenergy=freg) 
  # } else {
  #   peffects <- 1
  # }
  # if (length(nact) != 0) {
  #   neffects <- sapply(nact,repression)
  # } else {
  #   neffects <- 1
  # }
  # return(as.numeric((sum(log2(peffects))) + (sum(log2(neffects)))))
#   strength <- sum(capacities[[type]][regs, target] * inputs[regs])
#   for (regulator in regulators[[type]][[target]]) {
#     strength <- strength + as.numeric(inputs[regulator] * capacities[[type]][regulator, target])
#   }
#   return(exp(strength))
}

# Dummy function, to be implemented in future
reg.translation <- function() {
  return(null)
}

scaleScores <- function(vec) {
  l <- min(vec)
  u <- max(vec)
  if (abs(l - u) < 1e-50) {
    return(vec)
  }
  (vec - l) / (u - l)
}

scaleScores2 <- function(vec) {
  vec / max(abs(vec))
}

scaleScores3 <- function(vec,qt) {
  l <- min(vec)
  u <- max(vec)
  if (abs(l - u) < 1e-50) {
    return(vec)
  }
  temp <- vec
  temp[temp >= quantile(temp,probs=qt)] <- quantile(temp,probs=qt)
  u <- max(temp)
  (temp - l) / (u - l)
}

# Calculate updated node values for specified layer
# Node abundance values are in log2 space
computeNewNodeValues <- function(inputs, targets, input_types=NULL, target_type, deg_rate, capacities) {
  # Use hard-coded rules first
  if (target_type == "Protein") {
    mRNAs <- inputs$mRNAs
    effects <- as.numeric(mRNAs)
    names(effects) <- names(targets)
    # This part is a little weird...no protein degradation taken into consideration here?
    return(log2(2^targets + 2^effects))
  }
  if (target_type == "mRNA") {
    proteins <- inputs$proteins
    miRNAs <- scaleScores(inputs$miRNAs)
    # Protein values are on log2 scale
    pr.effects <- sapply(names(targets), regulation,
                         inputs=2^proteins, 
                         type="Pr2mRNA",
                         capacities=capacities)
    # mRNA / miRNAs are on log2 scale
    # mi.effects <- sapply(names(targets), regulation,
    #                      inputs=2^miRNAs, 
    #                      type="mi2mRNA",
    #                      capacities=capacities)
    mi.effects <- 0
    if (mi.effects == 0) {
      effects <- pr.effects
    } else {
      effects <- {0.9 * pr.effects + 0.1 * mi.effects} 
    }
    names(effects) <- names(targets)
    res <- targets + log2(1 - deg_rate) + effects
    #res[res < 0] <- 0
    #head(targets * effects)
    return(res)
#     return(targets * {pr.effects * mi.effects})
  }
  if (target_type == "miRNA") {
    proteins <- inputs$proteins
    #miRNAs <- scaleScores(inputs$miRNAs)
    pr.effects <- sapply(names(targets), regulation,
                         inputs=proteins, 
                         type="Pr2miRNA",
                         capacities=capacities)
    mi.effects <- 0
#     mi.effects <- sapply(targets, regulation,
#                          inputs=miRNAs, 
#                          type="mi2miRNA")
    effects <- pr.effects
    names(effects) <- names(targets)
    res <- targets + log2(1 - deg_rate) + effects
    #res[res < 0] <- 0
    #head(targets * effects)
    return(res)
#     return(targets * {pr.effects * mi.effects})
  }
}

# Update transition probabilities for random walk
# updateTransitionProb <- function(nvalues, weights) {
#   n <- length(nvalues)
#   Arow <- matrix(nvalues,ncol=n,nrow=1)
#   Acol <- matrix(nvalues,ncol=1,nrow=n)
#   A <- Acol %*% Arow
#   W <- A * weights
#   P <- t(apply(W,1,function(x){
#     temp <- sum(x)
#     if (sum(x) != 0) {
#       return(x/sum(x))
#     } else {
#       return(rep(0,ncol(W)))
#     }
#   }))
#   z <- which(rowSums(P)==0)
#   if (length(z)) {
#     diag(P[z,z]) <- 1 
#   }
# #   P <- t(apply(array(1:n),1,function(x){
# #     arr <- W[x,]
# #     s <- sum(arr)
# #     if (s != 0) {
# #       return(arr / s)
# #     }
# #     res <- rep(0,n)
# #     res[x] <- 1
# #     return(res)
# #   }))
#   return(P)
# }

updateTransitionProb <- function(nvalues, weights) {
    n <- length(nvalues)
    Arow <- matrix(nvalues,ncol=n,nrow=1)
    Acol <- matrix(nvalues,ncol=1,nrow=n)
    A <- Acol %*% Arow
    W <- A * weights
    rs <- rowSums(W)
    z <- which(rs==0)
    if (length(z)) {
        diag(W[z,z]) <- 1
        rs[z] <- 1
    }
    P <- as.matrix(W) / rs
    #   P <- t(apply(array(1:n),1,function(x){
    #     arr <- W[x,]
    #     s <- sum(arr)
    #     if (s != 0) {
    #       return(arr / s)
    #     }
    #     res <- rep(0,n)
    #     res[x] <- 1
    #     return(res)
    #   }))
    return(P)
}

updateTransitionProbExpRank <- function(nvalues, weights) {
  n <- length(nvalues)
  Arow <- matrix(nvalues,ncol=n,nrow=1)
  Acol <- matrix(nvalues,ncol=1,nrow=n)
  A <- Acol %*% Arow
  W <- A * weights
  # rs <- rowSums(W)
  # z <- which(rs==0)
  # if (length(z)) {
  #   diag(W[z,z]) <- 1
  #   rs[z] <- 1
  # }
  # P <- as.matrix(W) / rs
  #   P <- t(apply(array(1:n),1,function(x){
  #     arr <- W[x,]
  #     s <- sum(arr)
  #     if (s != 0) {
  #       return(arr / s)
  #     }
  #     res <- rep(0,n)
  #     res[x] <- 1
  #     return(res)
  #   }))
  return(W)
}


topologicalBiasMatrix <- function(P) {
  tbm <- P
  tbm[P!=0] <- 1
  tbm <- tbm / rowSums(tbm)
  return(tbm)
}

# The PageRank algorithm (i.e. random walk with restart)
PageRank <- function(p_init, P, u=0.7, maxiter=100, errtol=1e-7,tp=F) {
  n <- length(p_init)
  if (tp) {
    restart <- matrix(rep(1,n),ncol=1,nrow=n) %*% matrix(rep(1/n,n),ncol=n,nrow=1)
  } else {
    restart <- matrix(rep(1,n),ncol=1,nrow=n) %*% matrix(p_init,ncol=n,nrow=1)
  }
  #restart <- matrix(rep(1,n),ncol=1,nrow=n)
  Q <- t(restart * matrix(u,n,n) + P * matrix(1-u,n,n))
  p_now <- matrix(p_init,ncol=1,nrow=n)
  p_prev <- p_now
  error <- 1e6
  iter <- 1
  while (error >= errtol & iter <= maxiter) {
    p_now <- Q %*% p_prev
    # Error is quantified as Frobenius norm
    error <- norm(abs(p_now-p_prev),"F")
    p_prev <- p_now
    iter <- iter + 1
  }
  if (iter > maxiter) {
    cat("Warning: Algorithm did not converge before reaching max. number of iterations...")
    cat("\n")
  }
  return(list(PRScores=p_now,error=error))
}

# 2017-01-05
# Sharp approximation PageRank algorithm by Chung et al.
# beta is defined as 2*alpha/(1-alpha)
# Need a less 'sharp' approximation function first
# Potential problem is that the algorithm is built for undirected networks!!!
# set.seed(123456789)
# library(igraph)
# g <- erdos.renyi.game(10000,0.01,"gnp",directed=T)
# g <- set.edge.attribute(g,"weight",E(g),runif(length(E(g)),-10,10))
# V(g)$name <- as.character(c(1:length(V(g))))

# Another method that treats negative links is the PageTrust algorithm
# by Kerchove et al. (2008). This algorithm allows for directed networks.
# 'Memory after zapping' not considered in this version.

# pageTrust <- function(g, alpha, beta, restart, p_init, epsilon=1e-7) {
#   # First get positive and negative links
#   g.pos <- subgraph.edges(g, eids=which(E(g)$weight>0),
#                           delete.vertices = F)
#   adj.pos <- data.matrix(get.adjacency(g.pos))
#   d.pos <- degree(g.pos)
#   g.neg <- subgraph.edges(g, eids=which(E(g)$weight<0),
#                           delete.vertices = F)
#   adj.neg <- data.matrix(get.adjacency(g.neg))
#   d.neg <- degree(g.neg)
#   N <- length(V(g))
#   # Initialize
#   p_curr <- p_init
#   P_curr_tilde <- data.matrix(get.adjacency(g.neg))
#   P_curr <- P_curr_tilde
#   t <- 0
#   # Build Google matrix for the subgraph induced by positive edges
#   G <- t((alpha * 1 / d.pos) * (adj.pos)) + matrix(1 - alpha, N, N) * restart
#   p_next <- rep(0,N)
#   P_next <- matrix(0,N,N)
#   P_next_tilde <- matrix(0,N,N)
#   while (TRUE) {
#     for (i in c(1:N)) {
#       p_next[i] <- (1 - P_curr_tilde[i,i])^beta * (G[i,] %*% p_curr)
#       # Build transition matrix T(t)
#       # Use matrix operations wherever possible!
#       T_curr_num <- G %*% diag(p_curr)
#       neighborhood.scores <- sapply(c(1:N),function(v){
#         in.vertices <- as.numeric(neighborhood(g.pos,1,v,mode=c("in"))[[1]])
#         alpha * sum(p_curr[in.vertices] / d.pos[in.vertices]) + (1 - alpha) * restart[v]
#       })
#       T_curr <- 1 / neighborhood.scores * T_curr_num
#       for (j in c(1:N)) {
#         P_next_tilde[i,j] <- T_curr[i,] %*% P_curr[,j]
#         if (get.edge.ids(g.neg,c(i,j))>0) {
#           P_next[i,j] <- 1
#         } else if (i==j) {
#           P_next[i,j] <- 0
#         } else {
#           P_next[i,j] <- P_next_tilde[i,j]
#         }
#       }
#     }
#     p_next <- p_next / sum(p_next)
#     t <- t + 1
#     p_curr <- p_next
#     P_curr <- P_next
#     P_curr_tilde <- P_next_tilde
#     if (max(abs(p_next - p_curr)) <= epsilon) {
#       break
#     }
#     cat(".")
#   }
#   return(list(res=p_curr,distrust=P_curr))
# }
# system.time(res <- pageTrust(g,0.8,1,
#                              rep(1/length(V(g)),length(V(g))),
#                              rep(1/length(V(g)),length(V(g)))))


normalize <- function(x) {
  if (sum(x) != 0) {
    x / sum(x)
  } else {
    x 
  }
}

# Another alternative: 'exponential ranking' by Traag et al. (2010)
# There is a large number problem here!!
expRank <- function(p_init, A, mu, maxiter=100, errtol=1e-7,tp=F) {
  n <- length(p_init)
  p_now <- normalize(matrix(p_init,ncol=1,nrow=n))
  p_prev <- p_now
  error <- 1e6
  iter <- 1
  while (error >= errtol & iter <= maxiter) {
    p_now <- exp(1 / mu * (t(A) %*% p_prev))
    p_now <- normalize(p_now)
    # Error is quantified as Frobenius norm
    error <- norm(abs(p_now-p_prev),"F")
    p_prev <- p_now
    iter <- iter + 1
    #cat(".")
  }
  if (iter > maxiter) {
    cat("Warning: Algorithm did not converge before reaching max. number of iterations...")
    cat("\n")
  }
  return(list(PRScores=p_now,error=error))
}

expRank.test <- function(p_init, A, mu, maxiter=100, errtol=1e-7,tp=F) {
  n <- length(p_init)
  p_now <- normalize(matrix(p_init,ncol=1,nrow=n))
  p_prev <- p_now
  error <- 1e6
  iter <- 1
  while (error >= errtol & iter <= maxiter) {
    p_now <- exp(1 / mu * (t(A) %*% p_prev))
    p_now <- normalize(p_now)
    # Error is quantified as Frobenius norm
    error <- norm(abs(p_now-p_prev),"F")
    p_prev <- p_now
    iter <- iter + 1
    cat(as.character(error))
    cat("\n")
    #cat(".")
  }
  if (iter > maxiter) {
    cat("Warning: Algorithm did not converge before reaching max. number of iterations...")
    cat("\n")
  }
  return(list(PRScores=p_now,error=error))
}


model.pred <- function(tfvec, mivec, fvecmi, fvectf, const) {
  if (length(fvecmi)==0 | sum(fvecmi)==length(fvecmi)) {
    temp <- log2(const) + sum(log2((1 + fvectf * (tfvec^2)) / (1 + tfvec^2)))
    temp
  } else {
    temp <- log2(const) + sum(0.9 * log2((1 + fvectf * (tfvec)^2) / (1 + tfvec^2))) + sum(0.1 * log2((1 + fvecmi * (mivec)^2) / (1 + mivec^2)))
    temp
  }
}
