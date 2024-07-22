# Jan 2, 2024: modify the JSD index: change posterior(alpha) is from the calibration sample

### functions in the MCNPS method
library(DescTools)

generate.Q_largeJ <- function(K, J, q, Nset_complete=1) {
  # need to be changed bc q cannot equal to 0.7
  
  R = 1:K
  
  ### generate # of items requiring x attributes
  binom.p = dbinom(R,K,q)
  select.od = order(binom.p,decreasing=TRUE)
  select.od = select.od[-which(select.od==1)]
  
  item.p = binom.p[-1]
  item.p = item.p/sum(item.p)
  item.p = c(0,item.p)
  item.p = item.p[select.od]
  
  # in this function, no need a complex function to generate Q for keys
  num.item = floor(c(K*Nset_complete,item.p*(J-K*Nset_complete)))
  num.item[K] = num.item[K] + J-sum(num.item)
  # By this generation, the item is relatively easy
  
  q.matrix <- do.call("rbind",replicate(Nset_complete,diag(K),simplify = FALSE))
  
  wk.item <- num.item
  for (i in 2:K){ #1:K, to run simulation I
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
  
  return(q.matrix)
  # return(num.item)
}

## need to be revised
mcq.complete.generate <- function(K,J,H,q,prop_B,Nset_complete){
  # This function may not be generalized to other situations
  # prop_B = 0, means no constraint on the proportion of binary items
  
  if(H >= 5){
    m0 = diag(K)
    m1 = matrix(nrow = floor(K/2)*3,ncol = (2+K))
    b1 = matrix(nrow = floor(K/2),ncol = K)
    for (j in 1: floor(K/2)){
      # m1[,1] = rep(1:floor(K/2), each=3)
      m1[,2] = rep(c(1,2,3),floor(K/2))
      
      e1 = m0[2*(j-1)+1,]
      e2 = m0[2*j,]
      
      m1[3*(j-1)+1,3:(2+K)] = e1+e2
      m1[3*(j-1)+2,3:(2+K)] = e1
      m1[3*j,3:(2+K)] = e2
      
      b1[j,] = e1+e2
    }
    
    complete.M = matrix(rep(t(m1),Nset_complete),ncol=ncol(m1),byrow=TRUE)
    complete.M[,1] = rep(1:(floor(K/2)*Nset_complete), each=3)
    
    binary.M = matrix(rep(t(b1),Nset_complete),ncol=ncol(b1),byrow=TRUE)
    
    J_complete = J-max(complete.M[,1])
    
    q.matrix = c()
    
    if (K%%2 == 1){
      q.matrix = matrix(rep(t(m0[K,]),Nset_complete),ncol=K,byrow=TRUE)
      J_complete = J_complete-Nset_complete
    }
    
    q.matrix = rbind(q.matrix,generate.Q_largeJ(K, J_complete, q))
    
    mcq2 = mcq.generate(H,q.matrix)
    
    mcq2$mcq[,1] = mcq2$mcq[,1] + max(complete.M[,1])
    all.mcq = rbind(complete.M,mcq2$mcq)
    all.key = c(rep(1,max(complete.M[,1])),mcq2$key)
    all.num.options = c()
    for (j in 1:J){
      all.num.options = c(all.num.options,nrow(all.mcq[all.mcq[,1]==j,,drop = FALSE]))
    }
    
    q.matrix = rbind(binary.M,q.matrix)
    
    binary.J = which(all.num.options==1)
    if(length(binary.J) < prop_B*J){ # In the current simulation, this will always hold.
      ToB = sample(which(all.num.options>1)[-c(1:max(complete.M[,1]))],prop_B*J-length(binary.J))
      
      for (j in 1:length(ToB)){
        # print(j)
        # print(which(all.mcq[,1]==ToB[j]))
        # all.mcq[which(all.mcq[,1]==ToB[j]),]
        
        all.mcq = all.mcq[-which(all.mcq[,1]==ToB[j])[-1],]
        
        # j = j+1
        # all.mcq
        # 
        # print(dim(all.mcq))
        all.num.options[ToB[j]] = 1
      }
    } else {
      print("The proportion of H*=1 is not met! Please revise the code.")
    }
    
  } else {
    q.matrix = generate.Q_largeJ(K, J, q, Nset_complete)
    
    mcq2 = mcq.generate(H,q.matrix)
    
    all.mcq = mcq2$mcq
    all.key = mcq2$key
    all.num.options = mcq2$H_star
    
    binary.J = which(all.num.options==1)
    if(length(binary.J) < prop_B*J){ # In the current simulation, this will always hold.
      ToB = sample(which(all.num.options>1),prop_B*J-length(binary.J))
      
      for (j in 1:length(ToB)){
        # print(j)
        # print(which(all.mcq[,1]==ToB[j]))
        # all.mcq[which(all.mcq[,1]==ToB[j]),]
        
        all.mcq = all.mcq[-which(all.mcq[,1]==ToB[j])[-1],]
        
        # j = j+1
        # all.mcq
        # 
        # print(dim(all.mcq))
        all.num.options[ToB[j]] = 1
      }
    } else {
      print("The proportion of H*=1 is not met! Please revise the code.")
    }
    
    
  }
  

  
  rownames(all.mcq) <- NULL
  colnames(all.mcq) <- c("Item","Option",paste0("Att",seq(1:K)))
  
  return(list("binaryQ"= q.matrix, "mcQ"=all.mcq,"key"=all.key,"H_star"=all.num.options))
  
}

## Select a subset of mcQ

# Q: mcQ for the whole test
# j.set: set of IDs of selected items
sub.mcQ = function(Q,j.set){
  ini.mcQ <- Q[which(Q[,1]==j.set[1]),,drop=FALSE]
  
  if (length(j.set)>1){
    for (i in 2:length(j.set)){
      ini.mcQ <- rbind(ini.mcQ,Q[which(Q[,1]==j.set[i]),])
    }
  }
  return(ini.mcQ)
}

Q.optimal.binary <- function(mcB,obs.x,key){
  ### obs.x should be a _length*2 matrix: the first column is the item ID 
  ### and the second column is the response
  
  
  K <- ncol(mcB)
  J <- nrow(mcB)
  LatentClass <- class.generate(K)
  
  item.id <- c(1:J)
  
  I.M <- diag(K)
  
  
  if (is.null(nrow(obs.x))){
    ### first item
    o.key = I.M[1,]
    
    s1 = which(apply(mcB,1,function(x){all.equal(x,o.key)})==TRUE)
    
    
  } else {
    
    blank = rep(0,K)
    for (j in 1:nrow(obs.x)){
      blank[j] = ifelse(obs.x[j,2] == key[obs.x[j,1]], 1, 0)
    }
    blank[nrow(obs.x)+1]=1
    
    num.absent = length(which(blank == 0))
    s1 = which(mcB[,nrow(obs.x)+1]==1)
    if (num.absent > 0){
      for (j in 1:num.absent){
        s1 = intersect(s1,which(mcB[,which(blank == 0)[j]]==0))
      }
    }

  } 
  
  
  next.J = ifelse(length(s1)>1,sample(s1,1),s1)
  
  return(next.J)
  
}

# 
# mcQ = mc.q
# key = key.all
# O = H
#### for the current simulations
Q.optimal.mc <- function(mcQ,obs.x,key,save.m,O){
  ### obs.x should be a _length*2 matrix: the first column is the item ID 
  ### and the second column is the response
  
  
  Q.mcQ <- mcQ[,-c(1:2)]
  item.id <- unique(mcQ[,1])
  
  J <- length(unique(item.id))
  K <- ncol(Q.mcQ)
  
  mcB = matrix(nrow = J, ncol = K)
  for (j in item.id){
    mcB[j,] = Q.mcQ[which(mcQ[,1]==j)[1],]
  }
  
  LatentClass <- class.generate(K)

  optimal.J.even <- floor(K/2)
  I.M <- diag(K)
  
  if(O==4){
    next.J = Q.optimal.binary(mcB,obs.x,key)
  } else {
    if (is.null(nrow(obs.x))){
      ### first item
      o.key = I.M[1,]+I.M[2,]
      o.Q = rbind(o.key,I.M[1,],I.M[2,])
      
      s1 = which(apply(mcB,1,function(x){all.equal(x,o.key)})==TRUE)
      s2 = which(save.m>2)
      s3 = intersect(s1,s2)
      
      good.J <- c()
      for (l in 1:length(s3)){
        q.J = Q.mcQ[which(mcQ[,1]==s3[l]),]
        good.J = c(good.J,sum(abs(q.J[1:3,] - o.Q))==0)
        
      }
      
      j1 = s3[good.J]
      
    } else if (floor(K/2)<=nrow(obs.x) && nrow(obs.x)<(K/2)){
      ### odd: the last item
      
      ### get the estimate of the previous response
      subQ <- sub.mcQ(mcQ,obs.x[,1])
      
      score.sub.op <- score.option(subQ,O)
      
      score.sub.response <- c()
      for (j in 1:nrow(obs.x)){
        gl <- score.sub.op[[j]]
        op <- c(1:O)
        score.sub.response <- c(score.sub.response,gl[,2][which(gl[,1]==obs.x[j,2])])
      }
      
      eta.sub.class <- epc.generate(subQ,O,key[obs.x[,1]])
      d.sub <- apply(t(eta.sub.class),1,function(x){x-score.sub.response})
      d.sub <- matrix(d.sub,byrow = FALSE,ncol=2^K)
      d.sub.hamming <- apply(d.sub,c(1,2),function(x){ifelse(x!=0,1,0)})
      
      class.sub <- LatentClass[which.min(colSums(d.sub.hamming)),]
      
      k_der = 1:(2*nrow(obs.x))
      k_un = (1:K)[-k_der]  # should be Kth entry
      
      must.zero = which(class.sub[k_der]==0)
      
      Q.mcQ_der = Q.mcQ[,k_der]
      Q.mcQ_un = Q.mcQ[,k_un,drop=FALSE]
      
      s1 = which(mcB[,k_un]==1) # item id

      good.J <- c()
      for (l in 1:length(s1)){

        if(!all(class.sub[k_der]==1)){
          # not master all 1:(K-1) attributes
          q.J = Q.mcQ_der[which(mcQ[,1]==s1[l]),,drop=FALSE]
          q.J2 = Q.mcQ_un[which(mcQ[,1]==s1[l]),,drop=FALSE]
          
          yes.J <- FALSE
          for (m in 1:nrow(q.J)){
            if (all(q.J[m,][must.zero]==0) && q.J2[m,]==1){
              yes.J = TRUE
            } else {
              yes.J = yes.J
            }
          }
          
          # yes.J = ifelse(mcB[s1[l],must.zero])
          
          
          if (yes.J){
            good.J = c(good.J,s1[l])
          } else {
            good.J = good.J
          }
        } else {
          # master all 1:(K-1) attributes
          good.J = c(good.J,s1[l])
        }
        
      }
      
      
      j1 = good.J
      
      
    } else {
      ### 2 to floor(K/2)
      subQ <- sub.mcQ(mcQ,obs.x[,1])
      
      score.sub.op <- score.option(subQ,O)
      
      score.sub.response <- c()
      for (j in 1:nrow(obs.x)){
        gl <- score.sub.op[[j]]
        op <- c(1:O)
        score.sub.response <- c(score.sub.response,gl[,2][which(gl[,1]==obs.x[j,2])])
      }
      
      eta.sub.class <- epc.generate(subQ,O,key[obs.x[,1]])
      d.sub <- apply(t(eta.sub.class),1,function(x){x-score.sub.response})
      d.sub <- matrix(d.sub,byrow = FALSE,ncol=2^K)
      d.sub.hamming <- apply(d.sub,c(1,2),function(x){ifelse(x!=0,1,0)})
      
      class.sub <- LatentClass[which.min(colSums(d.sub.hamming)),]
      
      k_der = 1:(2*nrow(obs.x))
      k_un = (1:K)[-k_der]
      
      Q.mcQ_der = Q.mcQ[,k_der]
      Q.mcQ_un = Q.mcQ[,k_un]
      
      v_test = I.M[k_un[1],]+I.M[k_un[2],]
      v_test_un = v_test[k_un]
      
      o.Q = rbind(v_test,I.M[k_un[1],],I.M[k_un[2],])
      o.Q_un = o.Q[,k_un]
      
      s1 = which(apply(mcB[,k_un],1,function(x){all.equal(x,v_test_un)})==TRUE)
      s2 = which(save.m>2)
      s3 = intersect(s1,s2)
      
      good.J <- c()
      for (l in 1:length(s3)){
        q.J = Q.mcQ_un[which(mcQ[,1]==s3[l]),]
        good.J = c(good.J,sum(abs(q.J[1:3,] - o.Q_un))==0)
      }
      j1 = s3[good.J]
      
      # Reduce(intersect, list(c(1,2),c(2,3),c(2,3,4)))
      ### optimal
      if(!all(class.sub[k_der]==1)){
        must.zero = which(class.sub[k_der]==0)
        
        good.J <- c()
        for (l in 1:length(j1)){
          q.J = Q.mcQ_der[which(mcQ[,1]==j1[l]),]
          good.J = c(good.J,!1 %in% q.J[,must.zero])
        }
        
        j1 = j1[good.J]
        
      } else {
        j1 = j1
      }
    }
    
    
    next.J = ifelse(length(j1)>1,sample(j1,1),j1)
  }
  
  
  
  
  return(next.J)
  
}



# etas: matrix of ideal responses J*M
# j1: administered J.set
# j2: remaining J.set
# v1, v2: two candidates
MCNPS <- function(v1,v2,etas,j2,O,Hj_star,ws = "",LS = NULL){
  # v1 = LatentClass[ est.mcnps$`distance order`[1],]; v2=LatentClass[ est.mcnps$`distance order`[2],];
  # etas=eta.class;
  # j2=IB[-j_mcnps[1:j.look]];
  # O = H; Hj_star = save.m;
  # ws = "Type II"
  
  # M0 = etas[j1,] administered J
  if (is.null(LS)){
    LatentClass <- class.generate(K)
    M_num = 2^K
  } else {
    LatentClass <- LS
    M_num = nrow(LS)
  }
  
  M = etas[j2,]
  
  
  J = nrow(M)
  K = length(v1)
  Hj_star = Hj_star[j2]
  
  
  c1 = which(apply(LatentClass,1,function(x){all.equal(x,v1)})==TRUE) #the group id
  c2 = which(apply(LatentClass,1,function(x){all.equal(x,v2)})==TRUE)
  
  d.set = numeric(J)
  for (j in 1:J){
    
    eta1 = M[j,c1]
    eta2 = M[j,c2]
    
    if (eta1==eta2){
      d.set[j] = 0
    } else if (eta1 && eta2 > 0){
      d.set[j] = 1
    } else {
      if (ws=="Type I"){
        w = 1
      } else if (ws=="Type II"){
        w = ifelse(Hj_star[j]==O,0,1-1/O) ### H_j* = H_J this item should not be selected
      } else {
        # w = 1-save.m[j2[j]]/O
        w = 1/Hj_star[j]-1/O
      }
      
      d.set[j] = w
    }
    
      
  }
  
  dist.max = max(d.set)
  candidate.J = which(d.set == dist.max)
  next.J = j2[ifelse(length(candidate.J)==1,candidate.J,sample(candidate.J))]
  
  return(next.J)
  
}


post.p <- function(P_j,post_t_1,dat,administer.J){
  J_t = length(administer.J)
  
  p = matrix(nrow = J_t,ncol = 2^K) #likelihood
  for (j in 1:J_t){
    p[j,] = P_j[[administer.J[j]]][dat[j],]
  }
  
  L = apply(p,2,prod)
  
  post_t = (post_t_1*L)/(sum(post_t_1*L))
  
  return(post_t)
}

# post_t_1: prior
JSD.func = function(P_j,post_t_1,dat,administer.J,remain.J){
  # P_j = item.par
  # post_t_1 = prior
  J_t = length(administer.J)
  M = length(post_t_1)
  
  p = matrix(nrow = J_t,ncol = M) #likelihood
  for (j in 1:J_t){
    p[j,] = P_j[[administer.J[j]]][dat[j],]
  }
  
  L = apply(p,2,prod)
  
  post_t = (post_t_1*L)/(sum(post_t_1*L))
  
  # ln.L = apply(log(p),2,sum)
  # ln.prior = log(prior)
  # which.max(ln.L + ln.prior)
  # exp(ln.L + ln.prior)
  
  
  est.att = which.max.randomtie(post_t)
  # if(post_t ==0){
  #   0.0001
  # }
  
  JSD = c()                                    
  for (j in remain.J){
    S1 = Entropy(P_j[[j]] %*% post_t,base = exp(1))
    
    S2 = 0
    for (c in 1:M){
      S2 = S2+post_t[c] * Entropy(P_j[[j]][,c], base = exp(1))
    }
    
    JSD = c(JSD,S1-S2)
    
    j.maxJSD = which.max.randomtie(JSD)
    
    next.J.JSD = remain.J[j.maxJSD]
    
  }
  post_t =  pmin(pmax(post_t, 0.00001), 1 - 0.00001)
  return(list("Next"=next.J.JSD,"Est" = est.att,"UpPrior" = post_t))
}



which.max.randomtie <- function(x,na.rm=TRUE){
  # from GDINA
  loc <- which(x==max(x,na.rm = na.rm))
  if(length(loc)>1){
    loc <- sample(loc,1)
  }
  return(loc)
}
