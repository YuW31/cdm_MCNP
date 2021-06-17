lemma.func <- function(item.q){
  ## By default, the first row is always the key
  K = ncol(item.q)
  r = rowSums(item.q)
  distractor.q <- item.q[which(r>0),]
  distractor.q <- item.q[-1,]
  key.q <- item.q[1,]
  usefulness = any(apply(distractor.q,1,function(x) all(x==key.q)))
  
  any.two <- combn(nrow(distractor.q),2)
  all.proper <- FALSE

  for (i in 1:ncol(any.two)){
    loc = any.two[,i]
    q1 = distractor.q[loc[1],]
    q2 = distractor.q[loc[2],]
    candi.q = q1+q2
    candi.q[which(candi.q > 0)] = 1
    candi.q = data.frame(t(candi.q))
    
    all.par.q <- sapply(apply(item.q,1,function(x){x-candi.q}),function(y){all(y==0)})
    if(!any(all.par.q)) {
      all.proper <- TRUE
      }
  }
  
  if (usefulness) {cat("not useful\n")}
  if (all.proper) {cat("not proper\n")}
  if (!(usefulness || all.proper)) {cat("useful and proper\n")}
  
}


ex1 <- matrix(c(1,1,1,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,1,0),nrow = 4,ncol = 5, byrow = TRUE)
lemma.func(ex1)

ex2 <- matrix(c(1,1,1,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,1),nrow = 4,ncol = 5, byrow = TRUE)
lemma.func(ex2)

ex3 <- matrix(c(1,1,1,0,0,1,1,0,0,1,1,0,0,0,0,0,0,0,0,0),nrow = 4,ncol = 5, byrow = TRUE)
lemma.func(ex3)
                         

