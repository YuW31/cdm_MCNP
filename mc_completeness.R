mc_completenes <- function(Q){ # input should be a data frame
  if (!is.data.frame(Q)){
    Q = data.frame(Q)
  }
  
  J.ID = unique(Q$Item)
  J = length(J.ID)
  K = ncol(Q)-2
  
  Q.body = Q[,-c(1:2)]
  
  Q_j = vector(length = J, mode = "list")
  
  for (j in 1:J){
    Q_j[[j]] = Q.body[Q$Item == J.ID[j],]
  }
  
  e_k_M = diag(K)
  
  # include all possible e_k
  condition1 = apply(e_k_M,1,function(x){
    any(apply(Q.body,1,function(y){
      all(y==x)
      }))
    })
  
  if(!all(condition1)){
    message(paste0("e_",which(condition1==FALSE)," does not exist."))
    stop("This Q-matrix is not complete because not all possible e_k are included!")
  }
  
  # e_k \preceq key or Lemma 3
  
  for (k in 1:K){
    e_k = e_k_M[k,]
    loc = apply(Q.body,1,function(y){
      all(y==e_k)
      })
    loc.J = Q$Item[loc]
    
    condition2.a = FALSE
    for (j1 in 1:length(loc.J)){
      condition2.a = condition2.a || all(Q_j[[loc.J[j1]]][1,] - e_k >= 0) # the key
    }
    
    # if e_k \preceq key is not satisfied, then proceed Lemma 3
    if (!condition2.a){
      
      condition2.b = TRUE
      mixed.g = rep(0,K)
      
      while(condition2.b){
        
        for (j1 in 1:length(loc.J)){
          mixed.g = mixed.g + Q_j[[loc.J[j1]]][1,]
        }
        
        mixed.g[mixed.g >1] = 1
        print(mixed.g)
        
        Q.e_k =  Q.body[Q.body[,1]==1,]
        
        loc = apply(Q.e_k[,-k],1,function(y){
          all(mixed.g[-k]-y >= 0)
        })
        loc.J = unique(Q$Item[Q.body[,1]==1][loc]) # find all Q_j
        
        if (is.na(loc.J)){
          stop("The Q-matrix is not complete because the examinees who have mastered Attribute ", k, "
               cannot be differentiated from those who have not mastered it.")
        }
        
        condition2.b2 = FALSE
        for (j1 in 1:length(loc.J)){
          condition2.b2= condition2.b2 || all(Q_j[[loc.J[j1]]][1,] - e_k >= 0) # the key requires e_k
        }
        
        condition2.b = !condition2.b2
        
      }
    }
     
  }
  
  print("The Q-matrix is complete!")
  
}
