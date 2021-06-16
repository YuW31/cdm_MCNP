### all NonParametric functions

#===Calculate PAR====
PAR=function(x,y){
  out=mean(1*(rowSums(abs(x-y))==0))
  return(out)
}
#===Calculate AAR====
AAR=function(x,y){
  out=mean(1*(x-y==0))
  return(out)
}

generate.Q <- function(K,J,q){
  R = 1:K
  
  ### generate # of items requiring x attributes
  binom.p = dbinom(R,K,q)
  select.od = order(binom.p,decreasing=TRUE)
  
  item.p = binom.p
  item.p = item.p/sum(item.p)
  item.p = item.p[select.od]
  
  num.item = rep(0,K)
  for (i in 1:(K-1)){
    det = sum(num.item)
    if (det == J){
      break
    }
    
    largest.J <- ceiling(max((J-det)*item.p)*1.3)
    smallest.J <- floor(max((J-det)*item.p)*0.8)
    
    x = J
    while (x > largest.J || x < smallest.J){
      x = rbinom(1,J-det,item.p[1])
    }
    
    num.item[select.od[i]] = x
    
    item.p = binom.p[-select.od[i]]
    item.p = item.p/sum(item.p)
    od = order(item.p,decreasing=TRUE)
    item.p  = item.p[od]
  }
  num.item[1] = ifelse(num.item[1]>K,num.item[1],K)
  num.item[K] = max(J-sum(num.item),0)
  num.item
  
  # replicate: simplify = FALSE: list; or matrix combined by column
  con.1 = TRUE
  con.2 = TRUE
  
  while (con.1 || con.2){
    q.matrix <- diag(K)
    wk.item <- num.item
    wk.item[1] = wk.item[1]-K
    
    for (i in 1:K){
      if (wk.item[i]==0){
        next
      } else {
        test.m <- matrix(0,nrow=wk.item[i],ncol=K)
        sub.q <- apply(test.m,1,function(x){
          t = sample(1:K,i,replace = FALSE)
          x[t] = 1
          return(x)
        }
        )
        sub.q = t(sub.q)
        q.matrix <- rbind(q.matrix,sub.q)
      }
    }
    con.1 <- 0 %in% colSums(q.matrix)
    con.2 <- max(colSums(q.matrix))/min(colSums(q.matrix))>2
  }
  return(q.matrix)
}

# O=3,4,5
mcq.generate = function(O,Q){
  mc.q <- c()
  key <- c()
  save.m <- c()
  
  J <- nrow(Q)
  K <- ncol(Q)
  q.matrix <-cbind(Q,rowSums(Q))
  
  for (i in 1:J){
    num <- q.matrix[i,K+1]
    
    if (O==3){
      p1 = c(0.8,0.1,0.1)
      p2 = c(0.2,0.2,0.6)
      p3 = c(0.2,0.2,0.2)
    } else if (O==4){
      p1 = c(0.7,0.2,0.1,0)
      p2 = c(0.1,0.2,0.2,0.5)
      p3 = c(0.1,0.2,0.2,0.1)
    } else {
      p1 = c(0.7,0.1,0.1,0.1,0)
      p2 = c(0.1,0.1,0.1,0.1,0.6)
      p3 = c(0.1,0.1,0.1,0.1,0.1)
    }
    
    if (num == 1) {
      p = p1 
    } else if (num >= O){
      p = p2
    } else {
      p = p3
      p[num] = 0.6
    }
    
    m <- sample(1:O,1,prob=p)
    sub.q <- matrix(0,nrow=O,ncol=K)
    od.option <- sample(1:O,O)
    id.key <- od.option[1]
    sub.q[1,] <- q.matrix[i,1:K]
    
    r.q.key <- which(q.matrix[i,1:K]==1)
    not.q.key <- which(q.matrix[i,1:K]==0)
    cb <- expand.grid(lapply(r.q.key,function(x)(c(0,x))))[1:(2^num-1),,drop=FALSE]  ## keep dataframe
    cb2 <- expand.grid(lapply(not.q.key,function(x)(c(0,x))))
    
    
    # remaining m-1 coded options
    if (m == 2){
      x1 = sample(nrow(cb),1)
      x2 = ifelse(x1==1,sample(2:nrow(cb2),1),sample(nrow(cb2),1))
      # prob(nesting within the key) = 0.5: prob=c(0.5,rep(0.5/(nrow(cb2)-1),nrow(cb2)-1))
      
      v <- rep(0,K)
      loc <- unlist(c(cb[x1,],cb2[x2,]))
      v[loc]=1
      sub.q[2,] <- v
      
    } else if (m > 2 && m <= nrow(cb)) {
      ### two possibilities: 1) nested within the key; 2) nested within the first generated "large" q-vector
      ### No hybrid!!!
      
      nest <- sample(c(1,0),1)
      nest <- ifelse(num==K,1,nest)
      
      if (nest==1){
        ### nested within the key
        x1 = sample(2:nrow(cb),m-1,replace = FALSE)
        
        v.matrix <- matrix(0,nrow=m-1,ncol=K)
        for (ttt in 1:(m-1)){
          loc <- unlist(cb[x1[ttt],])
          v.matrix[ttt,loc]=1
        }
        
        largest.d <- v.matrix[which.max(rowSums(v.matrix )),]
        r.large.d <- which(largest.d==1)
        sub.q[2,] <- largest.d
        
        full.q <- q.matrix[i,1:K]
        
        abc = 1
        complement = 0
        sub.ava = length(r.large.d)>1
        
        while (abc <= (m-2) && complement == 0){
          
          complt.largest <- full.q-largest.d
          complement = sample(c(0,1),1) ## of select the complement, the selection of distractors ends
          
          if (complement == 1){
            
            if (sum(complt.largest)==0){
              m=2+abc-1
            } else {
              sub.q[2+abc,] <- complt.largest
              m = 2+abc
            }
            
          } else {
            
            # r.large.d <- which(largest.d==1)
            # sub.ava = length(r.large.d)>1
            
            if (sub.ava){
              
              full.q <- sub.q[abc+1,]
              cb3 <- expand.grid(lapply(r.large.d,function(x)(c(0,x))))[2:(2^length(r.large.d)-1),,drop=FALSE]
              
              x1 = sample(2:nrow(cb3),m-1-abc,replace = TRUE)
              
              v.matrix <- matrix(0,nrow=m-1,ncol=K)
              for (ttt in 1:(m-1)){
                loc <- unlist(cb3[x1[ttt],])
                v.matrix[ttt,loc]=1
              }
              
              largest.d <- v.matrix[which.max(rowSums(v.matrix )),]
              r.large.d <- which(largest.d==1)
              sub.q[2+abc,] <- largest.d
              
              m = 2+abc
              abc = abc+1
              sub.ava = length(r.large.d)>1
              
            } else {
              m = 2+abc-1
            }
            
          }
          
        }
        
        
      } else {
        ### nested within the first generated "large" q-vector
        
        x2 = sample(2:nrow(cb2),m-1,replace = TRUE)
        
        v.matrix <- matrix(0,nrow=m-1,ncol=K)
        for (ttt in 1:(m-1)){
          loc <- unlist(cb2[x2[ttt],])
          v.matrix[ttt,loc]=1
        }
        
        largest.d <- v.matrix[which.max(rowSums(v.matrix )),]
        r.large.d <- which(largest.d==1)
        sub.q[2,] <- largest.d
        
        full.q <- largest.d
        
        
        abc = 1
        complement = 0
        sub.ava = length(r.large.d)>1
        
        while (abc <= (m-2) && complement == 0){
          
          complt.largest <- full.q-largest.d
          # complement = sample(c(0,1),1) ## of select the complement, the selection of distractors ends
          
          if (complement == 1){
            
            if (sum(complt.largest)==0){
              m=2+abc-1
            } else {
              sub.q[2+abc,] <- complt.largest
              m = 2+abc
            }
            
          } else {
            
            if (sub.ava){
              
              full.q <- sub.q[abc+1,]
              cb3 <- expand.grid(lapply(r.large.d,function(x)(c(0,x))))[2:(2^length(r.large.d)-1),,drop=FALSE]
              
              x1 = sample(2:nrow(cb3),m-1-abc,replace = TRUE)
              
              v.matrix <- matrix(0,nrow=m-1,ncol=K)
              for (ttt in 1:(m-1)){
                loc <- unlist(cb3[x1[ttt],])
                v.matrix[ttt,loc]=1
              }
              
              largest.d <- v.matrix[which.max(rowSums(v.matrix )),]
              r.large.d <- which(largest.d==1)
              sub.q[2+abc,] <- largest.d
              
              m = 2+abc
              abc = abc+1
              sub.ava = length(r.large.d)>1
              complement = sample(c(0,1),1) ## of select the complement, the selection of distractors ends
              
            } else {
              m = 2+abc-1
            }
            
          }
          
        }
        
      }
      
    } else if (m > nrow(cb)){
      ### one possibility: nested within the first generated "large" q-vector
      
      x2 = sample(2:nrow(cb2),m-1,replace = TRUE)
      
      v.matrix <- matrix(0,nrow=m-1,ncol=K)
      for (ttt in 1:(m-1)){
        loc <- unlist(cb2[x2[ttt],])
        v.matrix[ttt,loc]=1
      }
      
      largest.d <- v.matrix[which.max(rowSums(v.matrix )),]
      r.large.d <- which(largest.d==1)
      sub.q[2,] <- largest.d
      
      full.q <- largest.d
      # full.q <- rep(0,K)
      # full.q[unlist(cb2[nrow(cb2),])] = 1
      
      abc = 1
      complement = 0
      sub.ava = length(r.large.d)>1
      
      while (abc <= (m-2) && complement == 0){
        
        complt.largest <- full.q-largest.d
        
        if (complement == 1){
          
          if (sum(complt.largest)==0){
            m=2+abc-1
          } else {
            sub.q[2+abc,] <- complt.largest
            m = 2+abc
          }
          
        } else {
          
          if (sub.ava){
            
            full.q <- sub.q[abc+1,]
            cb3 <- expand.grid(lapply(r.large.d,function(x)(c(0,x))))[2:(2^length(r.large.d)-1),,drop=FALSE]
            
            x1 = sample(2:nrow(cb3),m-1-abc,replace = TRUE)
            
            v.matrix <- matrix(0,nrow=m-1,ncol=K)
            for (ttt in 1:(m-1)){
              loc <- unlist(cb3[x1[ttt],])
              v.matrix[ttt,loc]=1
            }
            
            largest.d <- v.matrix[which.max(rowSums(v.matrix )),]
            r.large.d <- which(largest.d==1)
            sub.q[2+abc,] <- largest.d
            
            m = 2+abc
            abc = abc+1
            sub.ava = length(r.large.d)>1
            complement = sample(c(0,1),1) ## of select the complement, the selection of distractors ends
            
            
          } else {
            m = 2+abc-1
          }
          
        }
        
      }
      
      
      
      
      # x2 = sample(2:nrow(cb2),m-1,replace = TRUE)
      # 
      # v.matrix <- matrix(0,nrow=m-1,ncol=K)
      # for (ttt in 1:(m-1)){
      #   loc <- unlist(cb2[x2[ttt],])
      #   v.matrix[ttt,loc]=1
      # }
      # 
      # largest.d <- v.matrix[which.max(rowSums(v.matrix )),]
      # sub.q[2,] <- largest.d
      # 
      # r.large.d <- which(largest.d==1)
      # 
      # if (length(r.large.d)>1){
      #   cb3 <- expand.grid(lapply(r.large.d,function(x)(c(0,x))))[2:(2^length(r.large.d)-1),,drop=FALSE]
      #   m = ifelse(nrow(cb3)>(m-2),m,nrow(cb3)+2)
      #   if(nrow(cb3)>(m-2)){
      #     x3 = sample(nrow(cb3),m-2)
      #   } else {
      #     x3 = sample(nrow(cb3),nrow(cb3))
      #   }
      #   
      #   for (xxx in 1:(m-2)){
      #     loc <- unlist(cb3[x3[xxx],])
      #     v <- rep(0,K)
      #     v[loc]=1
      #     sub.q[2+xxx,] <- v
      #   }
      # } else {
      #   m = 2
      # }
      
    }
    
    
    item.q <- cbind(rep(i,O),od.option,sub.q)
    
    dl <- which(rowSums(sub.q)==0)
    if (length(dl)==0){
      item.q = item.q
    } else {
      item.q = item.q[-which(rowSums(sub.q)==0),]
    }
    
    mc.q <- rbind(mc.q,item.q)
    key <- c(key,id.key)
    save.m <- c(save.m,m)
    
  }
  
  rownames(mc.q) <- NULL
  colnames(mc.q) <- c("Item","Option",paste0("Att",seq(1:K)))
  
  return(list("mcq"=mc.q,"key"=key,"coded options"=save.m))
}

class.generate = function(K){
  M <- diag(K)
  for (i in 2:K){
    M <- rbind(M,t(apply(combn(K,i),2,function(x){apply(M[x,],2,sum)})))
  }
  M <- rbind(0,M)
  return(M)
}

epc.generate = function(mcq,O,J,K,key,Class){
  Q <- mcq[, -c(1:2)]
  Item.info <- mcq[,1:2]
  item.no <- mcq[,1]
  coded.op <- mcq[,2]
  num.coded <- tabulate(item.no)   # number of coded option  save.m
  no.options <- rep(O, J)
  
  # "scored" option
  eta.class <- matrix(0,J,2^K)
  for(j in 1:J){
    
    Qj <- Q[which(item.no==j),,drop=FALSE]  # won't change data type
    Kj <- rowSums(Qj)
    kj.order <- order(Kj,decreasing = TRUE)  #get location;
    coded.op.j <- coded.op[which(item.no==j)]
    
    if (num.coded[j]>1){
      
      key.loc <- which(coded.op.j==key[j])
      if(key.loc!=kj.order[1]){
        kj.order <- c(key.loc,setdiff(kj.order,key.loc))  # make the answer goes first
      }
      
      
      Qj <- Qj[kj.order,]
      #et.label[[j]] <- c(apply(Qj,1,paste,collapse=""),paste(rep("*",ncol(Qj)),collapse = ""))
      A <- Class%*%t(Qj)
      B <- rbind(0,t(1*(A==(matrix(1,2^K)%*%t(rowSums(Qj))))))
      #matrix(unlist(lapply(apply(B,2,function(x){which(x==max(x))}),function(y){if (length(y)<O){c(y,rep(0,O-length(y)))}else{y=y}})),nrow=2^K,byrow=TRUE)
      
      #tmp <- eta(Qj)
      eta.j <- apply(B,2,which.max)
      max.j <- max(eta.j)
      eta.j <- eta.j - 1
      eta.j[eta.j==0] <- max.j
      eta.j <- num.coded[j]+1-eta.j
      
      
    }else{
      
      A <- Class%*%t(Qj)
      B <- rbind(0,t(1*(A==(matrix(1,2^K)%*%t(rowSums(Qj))))))
      
      eta.j <- apply(B,2,which.max)
      eta.j <- eta.j-1
    }
    
    eta.class[j,] <- eta.j
    
  }
  return(eta.class)
}

stu.generate = function(I,Class,N){
  K = ncol(Class)
  
  if (I==1){ #independent
    true.att <- c()
    for (i in 1:N){
      att <- Class[sample(1:2^K, 1),]
      true.att <- rbind(true.att,att)
      rownames(true.att) <- c()
    }
  }else{
    corr <- c(0,3,0.5)
    mu = rep(0,K)
    while (TRUE){
      gene.sigma <- function(){
        lower.sigma <- matrix(runif(K*K,corr[1],corr[2]),K)
        lower.sigma[lower.tri(lower.sigma)]<-0
        trial.sigma <- t(lower.sigma)+lower.sigma
        diag(trial.sigma) <- 1
        return(trial.sigma)
      }
      sigma <- gene.sigma()
      if(all(eigen(sigma)$values > 0)){
        break
      }
    }
    prob.att<-mvrnorm(N, mu, sigma)
    true.att <- apply(prob.att, c(1,2), function(x){if (x>0) {x=1} else {x=0} })
  }
  return(true.att)
}

eta.generate = function(mcq,att,J){
  N = nrow(att)
  K = ncol(att)
  
  item.no <- mcq[,1]
  num.coded <- tabulate(item.no)
  Item.info <- mcq[,1:2]
  
  Ideal <- att%*%t(mcq[,3:(2+K)])
  Ideal.met <- 1*(Ideal==(matrix(1,N)%*%t(rowSums(mcq[,3:(2+K)]))))
  Ideal.conj <- matrix(nrow=N,ncol = J)
  cum.save.m <- cumsum(num.coded)
  
  
  for (i in 1:J){
    
    if (i==1){
      y <- Ideal.met[,1:cum.save.m[i]]
    } else {
      y <- Ideal.met[,(cum.save.m[i-1]+1):cum.save.m[i]]
    }
    
    if (is.vector(y)){
      p.y <- y
      p.y[which(y==1)] <- Item.info[which(Item.info[,1]==i),2]
      Ideal.conj[,i] <- p.y
    } else {
      work.q <- mc.q[which(mc.q[,1]==i),]
      sum.q <- rowSums(work.q[,3:(2+K)])
      #work.q <- rbind(work.q[1,],work.q[order(sum.q,decreasing=TRUE)[-1],])
      p <- unlist(lapply(apply(y,1,function(x){which(x==1)}),function(y){
        if (length(y)==0){
          aaa = 0
        } else if (length(y)==1){
          aaa = y
        } else {
          aaa =  ifelse(y[1]==1,y[1],y[which.max(sum.q[y])])
        }
        return (aaa)}))
      if(length(p)==0){
        p.y = rep(0,N)
      } else {
      p.y <- p
      p.y[which(is.na(p))] <- 0
      for (j in 1:ncol(y)){
        p.y[which(p==j)] <- Item.info[which(Item.info[,1]==i),2][j]
      }
      }
      Ideal.conj[,i] <- p.y
    }
  }
  return(Ideal.conj)
}

prob.generate = function(O,eta,corr.q){
  N = nrow(eta)
  J = ncol(eta)
  
  prob <- vector(length = N, mode = "list")
  for (i in 1:N){
    # 0.06; 0.82; 0.25: same as Jimmy
    #op <- mc.q[,1:2][which(Ideal.conj[i,]>0),,drop=FALSE]
    sub.r <- matrix(1/O,nrow = J, ncol = O)
    for (j in 1:J){
      q = corr.q[j]
      if (eta[i,j]>0){
        sub.r[j,] = (1-q)/(O-1)
        sub.r[j,eta[i,j]] = q
      }
    }
    prob[[i]] <- sub.r
  }
  
  return(prob)
}

score.option = function(mcq,O){
  J = max(mcq[,1])
  Q <- mcq[, -c(1:2)]
  Item.info <- mcq[,1:2]
  item.no <- mcq[,1]
  coded.op <- mcq[,2]
  num.coded <- tabulate(item.no)   # number of coded option  save.m
  
  op <- vector(mode="list",length=J)
  
  
  for (i in 1:J){
    Qj <- Q[which(item.no==i),,drop=FALSE]  # won't change data type
    Kj <- rowSums(Qj)
    kj.order <- order(Kj,decreasing = TRUE)  #get location;
    coded.op.j <- coded.op[which(item.no==i)]
    
    sub.info <- Item.info[which(Item.info[,1]==i),,drop=FALSE]
    
    if (num.coded[i]>1){
      
      key.loc <- which(coded.op.j==key[i])
      if(key.loc!=kj.order[1]){
        kj.order <- c(key.loc,setdiff(kj.order,key.loc))  # make the answer goes first
      }
      
      if (num.coded[i]==O){
        g = c(O:1)
      } else {g = c(c(nrow(sub.info):1),rep(0,O-nrow(sub.info)))}
      
      
      ans = c(sub.info[,2][kj.order],setdiff(c(1:O),sub.info[,2][kj.order]))
      score <- cbind(ans,g)
      
    } else {
      ans =c(unname(sub.info[,2]),setdiff(c(1:O),sub.info[,2]))
      g = c(1,rep(0,(O-1)))
      score <- cbind(ans,g)
    }
    
    op[[i]] <- score
    
  }
  
  return(op)
}

get.dis <- function(d,M){
  J = ncol(M)
  v = c()
  for (i in 1:J){
    v = c(v,M[d[i],i])
  }
  return(sum(v))
}

mc.npc <- function(dat,eta.class,w.class){
  N = nrow(dat)
  dis.ham <- vector(length = N, mode = "list")
  for (i in 1:N){
    d <- apply(t(eta.class),1,function(x){x-dat[i,]})
    d.hamming <- apply(d,c(1,2),function(x){ifelse(x!=0,1,0)})
    
    sum.d <- diag(t(d.hamming) %*% w.class)
    min.d <- which(sum.d==min(sum.d))
    x = ifelse(length(min.d)==1,min.d,sample(min.d,1))
    dis.ham[[i]] <- x
  }
  
  class.ham <- unlist(dis.ham)
  est.att.ham <- matrix(nrow=N,ncol=K)
  for (i in 1:length(class.ham)){
    est.att.ham[i,] <- LatentClass[class.ham[i],]
  }
  
  return (est.att.ham)
}

nested.q <- function(parent.q){
  m = length(which(parent.q != 0))
  q.loc = which(parent.q == 1)
  cb = do.call("cbind",expand.grid(lapply(q.loc,function(x)(c(0,x))))[2:(2^m-1),,drop=FALSE])
  x2 = sample(1:nrow(cb),1,replace = FALSE)
  candidate.q = rep(0,length(parent.q))
  candidate.q[cb[x2,]] = 1
  return(candidate.q)
}
