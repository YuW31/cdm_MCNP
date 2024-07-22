### all NonParametric functions
### Updated on September 21, 2023: function to find eta epc.generate()
### Updated on September 14, 2023

#===Calculate PAR====
PAR = function(x, y) {
  out = mean(1 * (rowSums(abs(x - y)) == 0))
  return(out)
}

PAR.class = function(true, est) {
  # calculate PARs for each latent class
  # input: matrix and the order matters!
  
  true = as.matrix(true)
  est = as.matrix(est)
  
  K = ncol(true)
  N = nrow(true)
  LatentClass = class.generate(K)
  
  out = c()
  class.n = c()
  for (c in 1:2^K) {
    group = which(rowSums(abs(true - t(replicate(N, LatentClass[c, ])))) == 0)
    par <- mean(1 * (rowSums(abs(true[group, , drop = FALSE] - est[group, , drop = FALSE])) == 0))
    out <- c(out, par)
    class.n = c(class.n,length(group))
  }
  
  out = cbind(class.n,out)
  
  rownames(out) <- apply(LatentClass, 1, function(x) {
    paste(x, collapse = "")
  })
  
  colnames(out) <- c("ClassSize","PAR")
  out = data.frame(out)
  return(out)
}

PAR.class.k <- function(true, est){
  # this function is used to calculate the PARs across groups that require the same number of attributes
  # input: matrix and the order matters!
  K = ncol(true)
  N = nrow(true)
  LatentClass = class.generate(K)
  
  class.par = PAR.class(true,est)
  out = data.frame()
  for (k in 0:K){
    class.k = class.par[which(rowSums(LatentClass)==k),,drop = FALSE]
    par = sum(class.k$ClassSize * class.k$PAR, na.rm = TRUE)
    n = sum(class.k$ClassSize)
    out = rbind(out,c(n,par/n))
  }
  
  colnames(out) =  c("ClassSize","PAR")
  rownames(out) =  paste0("K",0:K)
  return(out)
}



#===Calculate AAR====
AAR = function(x, y) {
  out = mean(1 * (x - y == 0))
  return(out)
}

### generate binary Q for small J due to the constraint of Condition II
### Nset_complete: the number of sets of complete Q-matrix for binary items. The default value is 1.
generate.Q <- function(K, J, q, Nset_complete=1) {
  # need to be changed bc q cannot equal to 0.7
  
  
  R = 1:K
  
  ### generate # of items requiring x attributes
  binom.p = dbinom(R, K, q)
  select.od = order(binom.p, decreasing = TRUE)
  select.od = select.od[-which(select.od == 1)]
  
  item.p = binom.p[-1]
  item.p = item.p / sum(item.p)
  item.p = c(0, item.p)
  item.p = item.p[select.od]

  num.item = c(K*Nset_complete, rep(0, K - 1))
  for (i in 1:(K - 2)) {
    det = sum(num.item)
    if (det >= J) {
      break
    }
    
    largest.J <- ceiling(max((J - det) * item.p) * 1.3)
    smallest.J <- floor(max((J - det) * item.p) * 0.8)
    
    x = J
    while (x > largest.J || x < smallest.J) {
      x = rbinom(1, J - det, item.p[i])
    }
    # print(x)
    
    num.item[select.od[i]] = x
    
    item.p = binom.p[-select.od[i]]
    item.p = item.p / sum(item.p)
    od = order(item.p, decreasing = TRUE)
    item.p  = item.p[od]
    
  }
  
  num.item[K] = max(J - sum(num.item), 0)
  # sum(num.item)
  
  # replicate: simplify = FALSE: list; or matrix combined by column
  con.1 = TRUE
  con.2 = TRUE
  
  while (con.1 || con.2) {
    q.matrix <- do.call("rbind",replicate(Nset_complete,diag(K),simplify = FALSE))
    wk.item <- num.item
    wk.item[1] = wk.item[1] - K
    
    for (i in 2:K) {
      if (wk.item[i] == 0) {
        next
      } else {
        test.m <- matrix(0, nrow = wk.item[i], ncol = K)
        sub.q <- apply(test.m, 1, function(x) {
          t = sample(1:K, i, replace = FALSE)
          x[t] = 1
          return(x)
        })
        sub.q = t(sub.q)
        q.matrix <- rbind(q.matrix, sub.q)
      }
    }
    con.1 <- 0 %in% colSums(q.matrix)
    con.2 <- max(colSums(q.matrix)) / min(colSums(q.matrix)) > 2
  }
  return(q.matrix)
  # return(num.item)
}


# O=3,4,5,6
# may need to change based on the Q-completeness
mcq.generate = function(O, Q) {
  mc.q <- c()
  key <- c()
  save.m <- c()
  
  J <- nrow(Q)
  K <- ncol(Q)
  q.matrix <- cbind(Q, rowSums(Q))
  
  for (i in 1:J) {
    num <- q.matrix[i, K + 1]
    
    if (O == 3) {
      p1 = c(0.8, 0.1, 0.1)
      p2 = c(0.2, 0.2, 0.6)
      p3 = c(0.2, 0.2, 0.2)
    } else if (O == 4) {
      p1 = c(0.7, 0.2, 0.1, 0)
      p2 = c(0.1, 0.2, 0.2, 0.5)
      p3 = c(0.1, 0.2, 0.2, 0.1)
    } else {
      # p1 = c(0.7, 0.1, 0.1, 0.1, 0)
      # p2 = c(0.1, 0.1, 0.1, 0.1, 0.6)
      # p3 = c(0.1, 0.1, 0.1, 0.1, 0.1)
      ## to acoomodate the large number of options
      p1 = rep(0.1,O)
      p1[1] = 0.7
      p2 = rep(0.1,O)
      p2[O] = 0.7
      p3 = rep(0.1,O)
    }
    
    if (num == 1) {
      p = p1
    } else if (num >= O) {
      p = p2
    } else {
      p = p3
      p[num] = 0.6
    }
    
    m <- sample(1:O, 1, prob = p)
    sub.q <- matrix(0, nrow = O, ncol = K)
    od.option <- sample(1:O, O)
    id.key <- od.option[1]
    sub.q[1, ] <- q.matrix[i, 1:K]
    
    r.q.key <- which(q.matrix[i, 1:K] == 1)
    not.q.key <- which(q.matrix[i, 1:K] == 0)
    cb <-
      expand.grid(lapply(r.q.key, function(x)
        (c(0, x))))[1:(2 ^ num - 1), , drop = FALSE]  ## keep dataframe
    cb2 <- expand.grid(lapply(not.q.key, function(x)
      (c(0, x))))
    
    
    # remaining m-1 coded options
    if (m == 2) {
      x1 = sample(nrow(cb), 1)
      x2 = ifelse(x1 == 1, sample(2:nrow(cb2), 1), sample(nrow(cb2), 1))
      # prob(nesting within the key) = 0.5: prob=c(0.5,rep(0.5/(nrow(cb2)-1),nrow(cb2)-1))
      
      v <- rep(0, K)
      loc <- unlist(c(cb[x1, ], cb2[x2, ]))
      v[loc] = 1
      sub.q[2, ] <- v
      
    } else if (m > 2 && m <= nrow(cb)) {
      ### two possibilities: 1) nested within the key; 2) nested within the first generated "large" q-vector
      ### No hybrid!!!
      
      nest <- sample(c(1, 0), 1)
      nest <- ifelse(num == K, 1, nest)
      
      m0 = m
      if (nest == 1) {
        ### nested within the key
        x1 = sample(2:nrow(cb), m - 1, replace = FALSE)
        
        v.matrix <- matrix(0, nrow = m - 1, ncol = K)
        for (ttt in 1:(m - 1)) {
          loc <- unlist(cb[x1[ttt], ])
          v.matrix[ttt, loc] = 1
        }
        
        largest.d <- v.matrix[which.max(rowSums(v.matrix)), ]
        r.large.d <- which(largest.d == 1)
        sub.q[2, ] <- largest.d
        
        full.q <- q.matrix[i, 1:K]
        
        abc = 1
        complement = 0
        sub.ava = length(r.large.d) > 1
        
        while (abc <= (m0 - 2) && complement == 0 && sub.ava) {
          complt.largest <- full.q - largest.d
          complement = sample(c(0, 1), 1) ## of select the complement, the selection of distractors ends
          
          if (complement == 1) {
            if (sum(complt.largest) == 0) {
              m = 2 + abc - 1
            } else {
              sub.q[2 + abc, ] <- complt.largest
              m = 2 + abc
            }
            
          } else {
            # r.large.d <- which(largest.d==1)
            # sub.ava = length(r.large.d)>1
            
            if (sub.ava) {
              full.q <- sub.q[abc + 1, ]
              cb3 <-
                expand.grid(lapply(r.large.d, function(x)
                  (c(0, x))))[2:(2 ^ length(r.large.d) - 1), , drop = FALSE]
              
              x1 = sample(2:nrow(cb3), m0 - 1 - abc, replace = TRUE)
              
              v.matrix <- matrix(0, nrow = m0 - 1 - abc, ncol = K)
              for (ttt in 1:nrow(v.matrix)) {
                loc <- unlist(cb3[x1[ttt], ])
                v.matrix[ttt, loc] = 1
              }
              
              largest.d <- v.matrix[which.max(rowSums(v.matrix)), ]
              r.large.d <- which(largest.d == 1)
              sub.q[2 + abc, ] <- largest.d
              
              m = 2 + abc
              abc = abc + 1
              sub.ava = length(r.large.d) > 1
              
            } else {
              m = 2 + abc - 1
            }
            
          }
          
        }
        
        
      } else {
        ### nested within the first generated "large" q-vector
        
        x2 = sample(2:nrow(cb2), m - 1, replace = TRUE)
        
        v.matrix <- matrix(0, nrow = m - 1, ncol = K)
        for (ttt in 1:(m - 1)) {
          loc <- unlist(cb2[x2[ttt], ])
          v.matrix[ttt, loc] = 1
        }
        
        largest.d <- v.matrix[which.max(rowSums(v.matrix)), ]
        r.large.d <- which(largest.d == 1)
        sub.q[2, ] <- largest.d
        
        full.q <- largest.d
        
        
        abc = 1
        complement = 0
        sub.ava = length(r.large.d) > 1
        
        while (abc <= (m0 - 2) && complement == 0 && sub.ava) {
          complt.largest <- full.q - largest.d
          # complement = sample(c(0,1),1) ## of select the complement, the selection of distractors ends
          
          if (complement == 1) {
            if (sum(complt.largest) == 0) {
              m = 2 + abc - 1
            } else {
              sub.q[2 + abc, ] <- complt.largest
              m = 2 + abc
            }
            
          } else {
            if (sub.ava) {
              full.q <- sub.q[abc + 1, ]
              cb3 <-
                expand.grid(lapply(r.large.d, function(x)
                  (c(0, x))))[2:(2 ^ length(r.large.d) - 1), , drop = FALSE]
              
              x1 = sample(2:nrow(cb3), m0 - 1 - abc, replace = TRUE)
              
              v.matrix <- matrix(0, nrow = m0 - 1- abc, ncol = K)
              for (ttt in 1:nrow(v.matrix)) {
                loc <- unlist(cb3[x1[ttt], ])
                v.matrix[ttt, loc] = 1
              }
              
              largest.d <- v.matrix[which.max(rowSums(v.matrix)), ]
              r.large.d <- which(largest.d == 1)
              sub.q[2 + abc, ] <- largest.d
              
              m = 2 + abc
              abc = abc + 1
              sub.ava = length(r.large.d) > 1
              complement = sample(c(0, 1), 1) ## of select the complement, the selection of distractors ends
              
            } else {
              m = 2 + abc - 1
            }
            
          }
          
        }
        
      }
      
    } else if (m > nrow(cb)) {
      ### one possibility: nested within the first generated "large" q-vector
      m0 = m
      
      x2 = sample(2:nrow(cb2), m - 1, replace = TRUE)
      
      v.matrix <- matrix(0, nrow = m - 1, ncol = K)
      for (ttt in 1:(m - 1)) {
        loc <- unlist(cb2[x2[ttt], ])
        v.matrix[ttt, loc] = 1
      }
      
      largest.d <- v.matrix[which.max(rowSums(v.matrix)), ]
      r.large.d <- which(largest.d == 1)
      sub.q[2, ] <- largest.d
      
      full.q <- largest.d
      # full.q <- rep(0,K)
      # full.q[unlist(cb2[nrow(cb2),])] = 1
      
      abc = 1
      complement = 0
      sub.ava = length(r.large.d) > 1
      
      while (abc <= (m0 - 2) && complement == 0 && sub.ava) {
        complt.largest <- full.q - largest.d
        
        if (complement == 1) {
          if (sum(complt.largest) == 0) {
            m = 2 + abc - 1
          } else {
            sub.q[2 + abc, ] <- complt.largest
            m = 2 + abc
          }
          
        } else {
          if (sub.ava) {
            full.q <- sub.q[abc + 1, ]
            cb3 <-
              expand.grid(lapply(r.large.d, function(x)
                (c(0, x))))[2:(2 ^ length(r.large.d) - 1), , drop = FALSE]
            
            x1 = sample(2:nrow(cb3), m0 - 1 - abc, replace = TRUE)
            
            v.matrix <- matrix(0, nrow = m0 - 1 - abc, ncol = K)
            for (ttt in 1:nrow(v.matrix)) {
              loc <- unlist(cb3[x1[ttt], ])
              v.matrix[ttt, loc] = 1
            }
            
            largest.d <- v.matrix[which.max(rowSums(v.matrix)), ]
            r.large.d <- which(largest.d == 1)
            sub.q[2 + abc, ] <- largest.d
            
            m = 2 + abc
            abc = abc + 1
            sub.ava = length(r.large.d) > 1
            complement = sample(c(0, 1), 1) ## of select the complement, the selection of distractors ends
            
            
          } else {
            m = 2 + abc - 1
          }
          
        }
        
      }
      
    }
    
    
    item.q <- cbind(rep(i, O), od.option, sub.q)
    
    dl <- which(rowSums(sub.q) == 0)
    if (length(dl) == 0) {
      item.q = item.q
    } else {
      item.q = item.q[-which(rowSums(sub.q) == 0), ]
    }
    
    mc.q <- rbind(mc.q, item.q)
    key <- c(key, id.key)
    save.m <- c(save.m, m)
    
  }
  
  save.m <- c()
  for (j in 1:J) {
    save.m <- c(save.m, length(which(mc.q[, 1] == j)))
  }
  
  rownames(mc.q) <- NULL
  colnames(mc.q) <- c("Item", "Option", paste0("Att", seq(1:K)))
  
  return(list(
    "mcq" = mc.q,
    "key" = key,
    "H_star" = save.m
  ))
}


# count # of q-vectors that measure each pattern
count.pattern <- function(Q, type) {
  if (type == "binary") {
    K = ncol(Q)
  } else if (type == "MC") {
    K = ncol(Q) - 2
    Q = Q[, -c(1:2)]
  } else {
    print("Not compatible for this Q-matrix for now.")
  }
  
  LatentClass = class.generate(K)
  J = nrow(Q)
  
  out = c()
  for (c in 1:2 ^ K) {
    group = which(rowSums(abs(Q - t(
      replicate(J, LatentClass[c, ])
    ))) == 0)
    out <- c(out, length(group))
  }
  
  names(out) <-
    apply(LatentClass, 1, function(x) {
      paste(x, collapse = "")
    })
  return(out)
}


class.generate = function(K) {
  M <- diag(K)
  for (i in 2:K) {
    M <- rbind(M, t(apply(combn(K, i), 2, function(x) {
      apply(M[x, ], 2, sum)
    })))
  }
  M <- rbind(0, M)
  return(M)
}

l.q <- function(Qj){
  ### requirement of input Qj
  # a matrix
  # the first row must be the key
  # no IDs
  
  # rules to order coded options
  # 1. The key is at the highest level
  # 2. coded distractors that require more skills are at the higher level
  # 3. coded distractors that require the same number of skills, 1 at the inconsistent entry......
  
  Qj_nokey <- Qj[-1,,drop = FALSE]
  Kj <- rowSums(Qj_nokey)
  Kj.order <- rank.Kj <- rank(Kj)
  if(any(duplicated(rank.Kj))){
    # currently, only one set of ties is considered
    # at most three coded distractors are considered
    loc.tie = c(which(duplicated(rank.Kj))[1]-1,which(duplicated(rank.Kj)))
    
    Qj_tie <- Qj_nokey[loc.tie,]
    tie.order <- order(apply(Qj_tie, 1, function(x) {
      as.integer(paste(x, collapse = ""))
    }))
    
    
    Kj.order[loc.tie] = mean(rank.Kj[loc.tie]-tie.order) + tie.order
    
  } else{
    # all coded distractors require different numbers of skills
    Kj.order <- rank(Kj)
  }
  
  
  Kj.order <- c(nrow(Qj),Kj.order)
  
  return(Kj.order)
  
}


epc.generate = function(mcq,O,key,LS = NULL){
  # for ideal responses
  
  Q <- mcq[, -c(1:2),drop = FALSE]
  Item.info <- mcq[,1:2]
  item.no <- mcq[,1]
  coded.op <- mcq[,2]
  num.coded <- tabulate(item.no)   # number of coded option  save.m
  
  J = length(unique(item.no))
  K = ncol(Q)
  if (is.null(LS)){
    Class <- class.generate(K)
    M = 2^K
  } else {
    Class <- LS
    M = nrow(LS)
  }
  
  no.options <- rep(O, J)
  
  # "scored" option
  eta.class <- matrix(0,J,M)
  # B.list = vector(length = J,mode = "list")
  for(j in 1:J){
    
    j.id = unique(item.no)[j]
    
    Qj <- Q[which(item.no==j.id),,drop=FALSE]  # won't change data type
    row.names(Qj) = NULL
    coded.op.j <- coded.op[which(item.no==j.id)]
    
    
    if (num.coded[j.id]>1){
      
      kj.order <- l.q(Qj)
      
      Qj <- Qj[order(kj.order),]
      A <- Class%*%t(Qj)
      B <- t(1*(A==(matrix(1,M)%*%t(rowSums(Qj)))))
      
      eta.j = apply(B,2,function(x){
        l = which(x==1)
        ifelse(length(l)==0,0,max(l))})
      
    }else{
      
      A <- Class%*%t(Qj)
      B <- rbind(0,t(1*(A==(matrix(1,M)%*%t(rowSums(Qj))))))
      
      eta.j <- apply(B,2,which.max)
      eta.j <- eta.j-1
    }
    
    eta.class[j,] <- eta.j
    
  }
  return(eta.class)
}

stu.generate = function(I, K, N) {
  Class = class.generate(K)
  
  if (I == 1) {
    #independent
    true.att <- c()
    for (i in 1:N) {
      att <- Class[sample(1:2 ^ K, 1), ]
      true.att <- rbind(true.att, att)
      rownames(true.att) <- c()
    }
  } else{
    # corr <- c(0,3,0.5)
    mu = rep(0, K)
    while (TRUE) {
      gene.sigma <- function() {
        # lower.sigma <- matrix(runif(K*K,corr[1],corr[2]),K)
        lower.sigma <- matrix(0.5, nrow = K, ncol = K)
        lower.sigma[lower.tri(lower.sigma)] <- 0
        trial.sigma <- t(lower.sigma) + lower.sigma
        diag(trial.sigma) <- 1
        return(trial.sigma)
      }
      sigma <- gene.sigma()
      if (all(eigen(sigma)$values > 0)) {
        break
      }
    }
    prob.att <- mvrnorm(N, mu, sigma)
    true.att <-
      apply(prob.att, c(1, 2), function(x) {
        if (x > 0) {
          x = 1
        } else {
          x = 0
        }
      })
  }
  return(true.att)
}

prob.generate = function(mcq, O, att, corr.q) {
  # att must be a matrix
  
  item.no <- unique(mcq[, 1])
  J = length(item.no)
  K = ncol(mcq) - 2
  
  N = nrow(att)
  
  LatentClass = class.generate(K)
  
  Ideal <- LatentClass %*% t(mcq[, 3:(2 + K)])
  Ideal.met <- 1 * (Ideal == (matrix(1, 2^K) %*% t(rowSums(mcq[, 3:(2 + K)]))))
  
  
  
  prob <- vector(length = 2^K, mode = "list")
  for (i in 1:2^K) {
    # 0.06; 0.82; 0.25: same as Jimmy
    #op <- mc.q[,1:2][which(Ideal.conj[i,]>0),,drop=FALSE]
    sub.r <- matrix(1 / O, nrow = J, ncol = O)
    
    i.eta <- Ideal.met[i, ]
    
    for (j in 1:J) {
      q = corr.q[j]
      
      work.q <- mcq[which(mcq[, 1] == item.no[j]), , drop = FALSE]
      work.eta <- i.eta[which(mcq[, 1] == item.no[j])]
      
      
      if (sum(work.eta) == 1 || work.eta[1] == 1) {
        sub.r[j, ] = (1 - q) / (O - 1)
        sub.r[j, work.q[which(work.eta == 1)[1], 2]] = q
      } else if (sum(work.eta) > 1) {
        possible.improper.Q = work.q[work.eta == 1, -(1:2)]
        
        ori.order <-
          order(rowSums(possible.improper.Q), decreasing = TRUE)
        all.include = all(apply(possible.improper.Q, 1, function(row)
          all(row <= possible.improper.Q[ori.order[1], ])))
        
        if (all.include) {
          work.eta[which(work.eta == 1)[-ori.order[1]]] = 0
        }
        
        sub.r[j, ] = (1 - q) / (O - 1)
        sub.r[j, work.q[which(work.eta == 1), 2]] = (q + (1 - q) * (length(which(
          work.eta == 1
        )) - 1) / (O - 1)) / length(which(work.eta == 1))
      }
      
    }
    
    prob[[i]] <- sub.r
  }
  
  
  prob.each <- vector(length = N, mode = "list")
  for (i in 1:N){
    prob.each[[i]] = prob[[which(apply(LatentClass,1,function(y){
      all.equal(att[i,],y)
    })==TRUE)]]
  }
  
  return(prob.each)
  
}

score.option = function(mcq, O) {
  J = length(unique(mcq[, 1]))
  Q <- mcq[,-c(1:2), drop = FALSE]
  K = ncol(Q)
  
  Item.info <- mcq[, 1:2, drop = FALSE]
  item.no <- mcq[, 1]
  coded.op <- mcq[, 2]
  num.coded <- tabulate(item.no)   # number of coded option  save.m
  
  op <- vector(mode = "list", length = J)
  
  mc.q = as.data.frame(mcq)
  colnames(mc.q) = c("Item", "Option", paste0("K", 1:K))
  
  save.m <- c()
  key <- c()
  for (i in 1:J) {
    j.id = unique(item.no)[i]
    
    key[i] <- mc.q$Option[min(which(mc.q$Item == j.id))]
    save.m <- c(save.m, length(which(mc.q$Item == j.id)))
  }
  
  
  for (i in 1:J) {
    j.id = unique(item.no)[i]
    
    Qj <- Q[which(item.no == j.id), , drop = FALSE]  # won't change data type
    
    ### wrong!!! the new index was updated
    # Kj <- rowSums(Qj)
    # kj.order <- order(Kj, decreasing = TRUE)  #get location;
    
    kj.order <- l.q(Qj)  #get location;
    coded.op.j <- coded.op[which(item.no == j.id)]
    
    sub.info <- Item.info[which(Item.info[, 1] == j.id), , drop = FALSE]
    
    if (num.coded[j.id] > 1) {
      key.loc <- which(coded.op.j == key[i])
      # if (key.loc != kj.order[1]) {
      #   kj.order <-
      #     c(key.loc, setdiff(kj.order, key.loc))  # make the answer goes first
      # }
      
      if (num.coded[j.id] == O) {
        g = c(O:1)
      } else {
        g = c(c(nrow(sub.info):1), rep(0, O - nrow(sub.info)))
      }
      
      
      ans = c(sub.info[, 2][order(kj.order,decreasing = TRUE)], setdiff(c(1:O), sub.info[, 2][order(kj.order,decreasing = TRUE)]))
      score <- cbind(ans, g)
      
      
    } else {
      ans = c(unname(sub.info[, 2]), setdiff(c(1:O), sub.info[, 2]))
      g = c(1, rep(0, (O - 1)))
      score <- cbind(ans, g)
    }
    
    op[[i]] <- score
    
  }
  
  return(op)
}

get.dis <- function(d, M) {
  J = ncol(M)
  v = c()
  for (i in 1:J) {
    v = c(v, M[d[i], i])
  }
  return(sum(v))
}

### input mc.q Observedresponse O
### Modified for CAT
### Cannot run for 1 item
algo_mc.npc <- function(mcq, dat, H, LS = NULL) {
  # input: matrix not a list
  # LS: a matrix
  
  K = ncol(mcq) - 2
  J = ncol(dat)
  N = nrow(dat)
  
  
  # mcq = as.data.frame(mcq)
  
  item.no <- mcq[, 1]
  
  save.m <- c()
  key <- c()
  for (i in 1:J) {
    j.id = unique(item.no)[i]
    
    key[i] <- mcq[, 2][min(which(mcq[, 1] == j.id))]
    save.m <- c(save.m, length(which(mcq[, 1] == j.id)))
  }
  
  if (is.null(LS)){
    LatentClass <- class.generate(K)
    M = 2^K
  } else {
    LatentClass <- LS
    M = nrow(LS)
  }
  
  
  eta.class <- epc.generate(mcq, H, key,LatentClass) # "scored" option
  
  
  score.op <- score.option(mcq, H)
  
  # number of coded options is larger than the noncoded options
  prob.j1 <- which(save.m < H & save.m > H / 2)
  
  prob.j2 <- which(save.m == H)
  
  
  w.class <- matrix(1, J, M)
  for (i in 1:M) {
    p = 1 / H
    w = seq(1, 0, by = -p)[-c(1, H, H + 1)]
    
    for (m in 2:(H - 1)) {
      ll1 <- intersect(which(eta.class[, i] == 0), which(save.m == m))
      w.class[ll1, i] = w[m - 1]
    }
    
    ll3 <- intersect(which(eta.class[, i] == 0), which(save.m == H))
    w.class[ll3, i] = 0
  }
  
  ### score examinees's responses
  score.response <- matrix(0,
                           nrow = N,
                           ncol = J,
                           byrow = TRUE)
  for (i in 1:J) {
    gl <- score.op[[i]]
    op <- c(1:H)
    for (j in 1:H) {
      score.response[, i][which(dat[, i] == j)] <-
        gl[, 2][which(gl[, 1] == j)]
    }
  }
  
  dis.ham <- vector(length = N, mode = "list")
  dis.order <- matrix(nrow = N, ncol = M)
  dis.all <- matrix(nrow = N, ncol = M)
  for (i in 1:N) {
    # cat("The estimation by the MC-NPC is in process: ",i,"/",N,"\n")
    d <- apply(t(eta.class), 1, function(x) {
      x - score.response[i, ]
    })
    d <- matrix(d, byrow = FALSE, ncol = M)
    d.hamming <- apply(d, c(1, 2), function(x) {
      ifelse(x != 0, 1, 0)
    })
    sum.d <- diag(t(d.hamming) %*% w.class)
    dis.all[i, ] <- sum.d
    
    # randomly order ties
    # sum.d_ties <- unique(sum.d)
    # sum.d_ties <- sum.d_ties[order(sum.d_ties, decreasing = FALSE)]
    # 
    # i.order <- c()
    # for (m in 1:length(sum.d_ties)){
    #   dis.class = which(sum.d==sum.d_ties[m])
    #   if(length(dis.class)==1){
    #     i.order = c(i.order,dis.class)
    #   } else {
    #     i.order = c(i.order,sample(dis.class))
    #   }
    # }
    # dis.order[i, ] <- i.order
    
    dis.order[i, ] <- order(sum.d, decreasing = FALSE)
    
    min.d <- which(sum.d == min(sum.d))
    x = ifelse(length(min.d) == 1, min.d, sample(min.d, 1))
    
    dis.ham[[i]] <- x
  }
  
  class.ham <- unlist(dis.ham)
  est.att.ham <- matrix(nrow = N, ncol = K)
  for (i in 1:length(class.ham)) {
    est.att.ham[i, ] <- LatentClass[class.ham[i], ]
  }
  
  est.att.ham = as.data.frame(est.att.ham)
  colnames(est.att.ham) = paste0("Att", 1:K)
  rownames(est.att.ham) = paste0("Examinee", 1:N)
  return(
    list(
      "est" = est.att.ham,
      "classID" = class.ham,
      "distance order" = dis.order,
      "distance" = dis.all,
      "latentclass" = LatentClass
    )
  )
}

# algo_mc.npc(mc.q,observed.response,4)

#######################################
######### NO Longer Used #############
#######################################

### before adding improper Q

### generate mc.q: no because not use lemma proposed in the paper
# O=3,4,5
free.mcq.generate = function(O, Q) {
  mc.q <- c()
  key <- c()
  save.m <- c()
  
  J <- nrow(Q)
  K <- ncol(Q)
  q.matrix <- cbind(Q, rowSums(Q))
  
  for (i in 1:J) {
    num <- q.matrix[i, K + 1]
    
    if (O == 3) {
      p1 = c(0.8, 0.1, 0.1)
      p2 = c(0.2, 0.2, 0.6)
      p3 = c(0.2, 0.2, 0.2)
    } else if (O == 4) {
      p1 = c(0.7, 0.2, 0.1, 0)
      p2 = c(0.1, 0.2, 0.2, 0.5)
      p3 = c(0.1, 0.2, 0.2, 0.1)
      
    } else {
      p1 = c(0.7, 0.1, 0.1, 0.1, 0)
      p2 = c(0.1, 0.1, 0.1, 0.1, 0.6)
      p3 = c(0.1, 0.1, 0.1, 0.1, 0.1)
    }
    
    if (num == 1) {
      p = p1
    } else if (num >= O) {
      p = p2
    } else {
      p = p3
      p[num] = 0.6
    }
    
    m <- sample(1:O, 1, prob = p)
    sub.q <- matrix(0, nrow = O, ncol = K)
    od.option <- sample(1:O, O)
    id.key <- od.option[1]
    sub.q[1, ] <- q.matrix[i, 1:K]
    
    
    # remaining m-1 coded options
    if (m > 1) {
      if (num == 1) {
        not.q.key <- which(q.matrix[i, 1:K] == 0)
        x = m - 1
        x = min(length(not.q.key), x)
        l <- sample(not.q.key, x)
        m = x + 1
        
        for (jj in 1:(m - 1)) {
          v <- rep(0, K)
          v[l[i = jj]] = 1
          sub.q[1 + jj, ] <- v
        }
        
      } else if (num == K) {
        r.q.key <- which(q.matrix[i, 1:K] == 1)
        cb <- expand.grid(lapply(r.q.key, function(x)
          (c(0, x))))
        m1 <- sample(2:(nrow(cb) - 1), m - 1)
        
        v1 <- cb[m1, ]
        for (ll in 1:nrow(v1)) {
          v <- rep(0, K)
          v[unlist(v1[ll, ][which(v1[ll, ] > 0)])] = 1
          sub.q[1 + ll, ] <- v
        }
        
      } else {
        r.q.key <- which(q.matrix[i, 1:K] == 1)
        not.q.key <- which(q.matrix[i, 1:K] == 0)
        cb <- expand.grid(lapply(r.q.key, function(x)
          (c(0, x))))
        #cb2 <- expand.grid(lapply(not.q.key,function(x)(c(0,x))))[-1,]
        
        
        # the probability that the option includes unrequired skill: 1 not include
        if (m - 1 > (nrow(cb) - 2)) {
          sel <- rbinom(nrow(cb) - 2, 1, 0.8)
          sel <- c(sel, rep(0, m - nrow(cb) + 1))
        } else {
          sel <- rbinom(m - 1, 1, 0.8)
        }
        
        
        m1 <- sample(2:(nrow(cb) - 1), length(which(sel == 1)))
        m2 <- sample(1:(nrow(cb) - 1), length(which(sel == 0)))
        
        n2 <-
          ifelse(length(not.q.key) > 1,
                 sample(not.q.key, min(
                   length(which(cb[m2, ] == 0)), length(not.q.key)
                 )),
                 not.q.key)
        
        v1 <- cb[m1, ]
        if (length(m1) > 0) {
          for (ll in 1:nrow(v1)) {
            v <- rep(0, K)
            v[unlist(v1[ll, ][which(v1[ll, ] > 0)])] = 1
            sub.q[1 + ll, ] <- v
          }
        }
        
        if (length(m2) > 0) {
          v2 <- cbind(cb[m2, ], n2)
          for (lll in 1:nrow(v2)) {
            v <- rep(0, K)
            v[unlist(v2[lll, ][which(v2[lll, ] > 0)])] = 1
            sub.q[1 + nrow(v1) + lll, ] <- v
          }
        }
      }
      
    }
    
    print(rbind(rep(i, K), sub.q))
    
    item.q <- cbind(rep(i, O), od.option, sub.q)
    
    dl <- which(rowSums(sub.q) == 0)
    if (length(dl) == 0) {
      item.q = item.q
    } else {
      item.q = item.q[-which(rowSums(sub.q) == 0), ]
    }
    
    mc.q <- rbind(mc.q, item.q)
    key <- c(key, id.key)
    save.m <- c(save.m, m)
    
  }
  
  rownames(mc.q) <- NULL
  colnames(mc.q) <- c("Item", "Option", paste0("Att", seq(1:K)))
  
  return(list(
    "mcq" = mc.q,
    "key" = key,
    "coded options" = save.m
  ))
}


epc.generate2 = function(mcq, O, key) {
  Q <- mcq[,-c(1:2)]
  Item.info <- mcq[, 1:2]
  item.no <- mcq[, 1]
  coded.op <- mcq[, 2]
  num.coded <- tabulate(item.no)   # number of coded option  save.m
  
  J = max(item.no)
  K = ncol(Q)
  no.options <- rep(O, J)
  Class <- class.generate(K)
  
  # "scored" option
  eta.class <- matrix(0, J, 2 ^ K)
  for (j in 1:J) {
    j.id = unique(item.no)[j]
    
    Qj <-
      Q[which(item.no == j.id), , drop = FALSE]  # won't change data type
    row.names(Qj) = NULL
    Kj <- rowSums(Qj)
    kj.order <- order(Kj, decreasing = TRUE)  #get location;
    coded.op.j <- coded.op[which(item.no == j.id)]
    
    if (num.coded[j.id] > 1) {
      key.loc <- which(coded.op.j == key[j])
      if (key.loc != kj.order[1]) {
        kj.order <-
          c(key.loc, setdiff(kj.order, key.loc))  # make the answer goes first
      }
      
      
      Qj <- Qj[kj.order, ]
      #et.label[[j]] <- c(apply(Qj,1,paste,collapse=""),paste(rep("*",ncol(Qj)),collapse = ""))
      A <- Class %*% t(Qj)
      B <- rbind(0, t(1 * (A == (
        matrix(1, 2 ^ K) %*% t(rowSums(Qj))
      ))))
      #matrix(unlist(lapply(apply(B,2,function(x){which(x==max(x))}),function(y){if (length(y)<O){c(y,rep(0,O-length(y)))}else{y=y}})),nrow=2^K,byrow=TRUE)
      
      #tmp <- eta(Qj)
      eta.j <- apply(B, 2, which.max)
      max.j <- max(eta.j)
      eta.j <- eta.j - 1
      eta.j[eta.j == 0] <- max.j
      eta.j <- num.coded[j.id] + 1 - eta.j
      
      
    } else{
      A <- Class %*% t(Qj)
      B <- rbind(0, t(1 * (A == (
        matrix(1, 2 ^ K) %*% t(rowSums(Qj))
      ))))
      
      eta.j <- apply(B, 2, which.max)
      eta.j <- eta.j - 1
    }
    
    eta.class[j, ] <- eta.j
    
  }
  return(eta.class)
}

eta.generate2 = function(mcq, att, J) {
  # except for noncoded option, others will be the original ID
  N = nrow(att)
  K = ncol(att)
  
  item.no <- mcq[, 1]
  num.coded <- tabulate(item.no)
  Item.info <- mcq[, 1:2]
  
  Ideal <- att %*% t(mcq[, 3:(2 + K)])
  Ideal.met <- 1 * (Ideal == (matrix(1, N) %*% t(rowSums(mcq[, 3:(2 + K)]))))
  Ideal.conj <- matrix(nrow = N, ncol = J)
  cum.save.m <- cumsum(num.coded)
  
  
  for (i in 1:J) {
    if (i == 1) {
      y <- Ideal.met[, 1:cum.save.m[i]]
    } else {
      y <- Ideal.met[, (cum.save.m[i - 1] + 1):cum.save.m[i]]
    }
    
    if (is.vector(y)) {
      # only one coded option
      p.y <- y
      p.y[which(y == 1)] <- Item.info[which(Item.info[, 1] == i), 2]
      Ideal.conj[, i] <- p.y
    } else {
      # more than one coded option
      work.q <- mcq[which(mcq[, 1] == i), ]
      sum.q <- rowSums(work.q[, 3:(2 + K)])
      #work.q <- rbind(work.q[1,],work.q[order(sum.q,decreasing=TRUE)[-1],])
      # decide which option is an ideal response
      p <-
        unlist(lapply(apply(y, 1, function(x) {
          which(x == 1)
        }), function(y2) {
          if (length(y2) == 0) {
            aaa = 0
          } else if (length(y2) == 1) {
            aaa = y2
          } else {
            aaa =  ifelse(y2[1] == 1, y2[1], y2[which.max(sum.q[y2])])
          }
          return (aaa)
        }))
      if (length(p) == 0) {
        p.y = rep(0, N)
      } else {
        p.y <- p
        p.y[which(is.na(p))] <- 0
        for (j in 1:ncol(y)) {
          p.y[which(p == j)] <- Item.info[which(Item.info[, 1] == i), 2][j]
        }
      }
      Ideal.conj[, i] <- p.y
    }
  }
  return(Ideal.conj)
}

prob.generate2 = function(O, eta, corr.q) {
  N = nrow(eta)
  J = ncol(eta)
  
  prob <- vector(length = N, mode = "list")
  for (i in 1:N) {
    # 0.06; 0.82; 0.25: same as Jimmy
    #op <- mc.q[,1:2][which(Ideal.conj[i,]>0),,drop=FALSE]
    sub.r <- matrix(1 / O, nrow = J, ncol = O)
    for (j in 1:J) {
      q = corr.q[j]
      if (eta[i, j] > 0) {
        sub.r[j, ] = (1 - q) / (O - 1)
        sub.r[j, eta[i, j]] = q
      }
    }
    prob[[i]] <- sub.r
  }
  
  return(prob)
}

### the output shouold be the ID of options (observed)
### designed for MC-NPS, fixed the first several observed responses as the ideal responses
etaID.generate = function(mcq,O,key){
  # for eta = 0, the option ID will be selected randomly
  
  Q <- mcq[, -c(1:2),drop = FALSE]
  Item.info <- mcq[,1:2]
  item.no <- mcq[,1]
  coded.op <- mcq[,2]
  num.coded <- tabulate(item.no)   # number of coded option  save.m
  
  J = length(unique(item.no))
  K = ncol(Q)
  Class = class.generate(K)
  no.options <- rep(O, J)
  
  # option ID
  etaID.class <- matrix(0,J,2^K)
  # B.list = vector(length = J,mode = "list")
  for(j in 1:J){
    
    j.id = unique(item.no)[j]
    
    Qj <- Q[which(item.no==j.id),,drop=FALSE]  # won't change data type
    row.names(Qj) = NULL
    coded.op.j <- coded.op[which(item.no==j.id)]
    
    
    if (num.coded[j.id]>1){
      
      kj.order <- l.q(Qj)
      h.id <- Item.info[which(item.no==j.id),2]
      
      Qj <- Qj[order(kj.order),]
      A <- Class%*%t(Qj)
      B <- t(1*(A==(matrix(1,2^K)%*%t(rowSums(Qj)))))
      
      etaID.j = apply(B,2,function(x){
        l = which(x==1)
        zero.ID = ifelse(length(c(1:O)[-h.id])==1,c(1:O)[-h.id],sample(c(1:O)[-h.id],1))
        ifelse(length(l)==0,zero.ID,h.id[order(kj.order)][max(l)])})
      
    }else{
      
      A <- Class%*%t(Qj)
      B <- rbind(0,t(1*(A==(matrix(1,2^K)%*%t(rowSums(Qj))))))
      
      eta.j <- apply(B,2,which.max)
      eta.j <- eta.j-1
      etaID.j <- ifelse(eta.j==1,key[j],sample(c(1:O)[-key[j]],1))
      
    }
    
    etaID.class[j,] <- etaID.j
    
  }
  return(etaID.class)
}
