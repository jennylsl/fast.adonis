#' @noRd
#' check the boot matrix
X.check.func <- function(X){
  X.boot.check <- apply(as.matrix(X[,-1]), 2, function(vec){
    length(unique(vec))>1
  })
  X.boot.check <- c(1,as.numeric(X.boot.check))==1
  return(X.boot.check)
}

#' @noRd
#' compute the conditional R2 sequentially
R2.select.fun <- function(ind.col, ind.covariate,
                          rhs, H.p2A, t.AK, weights, N){
  # if cleaned bootstrap sample do not contain specified varaibles
  if(all(!ind.covariate %in% ind.col)){
    R2.model <- NA
    # warnings: this time of bootstrap failed
  }
  else {
    model.X.n<- rhs[,c(TRUE,(ind.col %in% ind.covariate))]

    H.p1.model<- tcrossprod(model.X.n,solve(t(model.X.n*weights)%*%model.X.n))
    H.p2A.model <- H.p2A[c(TRUE, (ind.col %in% ind.covariate)),]
    t.HA.model <- sum(H.p1.model*weights*t(H.p2A.model))

    R2.model <- (t.HA.model-1/N*t.AK)/(-1/N*t.AK)
  }

  return(R2.model)
}

#' @noRd
#'  Sequential computing: create combination of variables put in the design matrix
num.list_gener <- function(nterms, num_orders = num_orders, order_list = order_list){
  num.list0 <- list(1)
  num.list1 <- lapply(1:  nterms, function(i) (1:nterms)[-i])[-nterms]
  num.list2 <- lapply(1:  nterms, function(i) 1:i)[-c(1,nterms)]
  num.list <- c(num.list0,num.list1,num.list2)

  # if there are additional orders
  if(!is.null(num_orders)){
    for(ind_order in 1:num_orders){
      # ind_order<-1
      num.list.temp <- lapply(1: nterms, function(i)order_list[[ind_order]][1:i])[-c(nterms)]
      # create new lists
      assign(paste0("num.list",(ind_order+2)),num.list.temp)
      num.list <- c(num.list,eval(parse(text=paste0("num.list",ind_order+2))))
    }
  }
  return(num.list)
}

#' @noRd
#' compute R2
R2.calc<-function( lhs, rhs, weights, nterms, ind.col,
                   t.AK, num.list, num_orders, SS=FALSE){

  # weights
  # rhs = X
  # lhs = A
  N <-  sum(weights)

  ## SSR
  # H=X(X^tX)^{-1}X^t
  # H.p=X(X^tX)^{-1}
  H.p1 <- tcrossprod(rhs,solve(t(rhs*(weights))%*%rhs))
  # warnins : X might be singular if large number of variables are categorical
  # H.P2A <-X^tA
  H.p2A<- tcrossprod(t(rhs*weights),lhs)

  # HA= X(X^tX)^{-1}X^tA
  t.HA <- sum(H.p1*weights*t(H.p2A))

  # HAK =X(X^tX)^{-1}X^tAK
  # tr(AK)=tr(HAK)=tr(HKA)=1/N*tr(HKAK)
  t.AK <- t.AK
  # SSTO: trace(G) = -1/N*tr(AK) = -1/N*sum(A)
  # R^2 = SSR/SSTO
  R2.o <- (t.HA-1/N*t.AK)/(-1/N*t.AK)

  R2_calc_result <- rep(R2.o,2)
  SS.RES <- -t.HA
  SS.TO <- sum(R2.o*(-1/N*t.AK)+SS.RES)

  if(SS==TRUE & nterms==1){
    SS.E <- R2.o*(-1/N*t.AK)
    result_r <- list(c(R2_calc_result,R2_calc_result),c(SS.E,SS.E,SS.RES, SS.TO))
  }

  # R2 computing
  if(nterms>1){
    # sequential
    # number of groups of variables in one time of boot (p single + p-2 sequential)
    var.total <- length(num.list)

    R2.single.boot <- lapply(num.list,function(ind.list) {
      R2.model <- R2.select.fun(ind.col, ind.list, rhs, H.p2A, t.AK, weights, N)
    } )
    SigSeq.R2 <- unlist( R2.single.boot)[1:(nterms*2-2)]

    # then conditional R2
    if(nterms>2){
      R2.condi.ori<- c(SigSeq.R2[c(1, nterms+1:(nterms-2))],R2.o)
    }
    if(nterms==2){
      R2.condi.ori<- c(SigSeq.R2[1],R2.o)
    }
    R2_condi_all_set <- R2.condi.ori[-1]-R2.condi.ori[-length(R2.condi.ori)]

    # return result if no additional orders
    ## total R2, marginal R2 for each variables from rhs, sequential R2
    R2_calc_result <- c(R2.o, SigSeq.R2[1], R2.o-SigSeq.R2[2:nterms], R2_condi_all_set)
    if(SS==TRUE){
      SS.E <- c(R2_calc_result*(-1/N*t.AK))
      result_r <- list(R2_calc_result,c(SS.E,SS.RES, SS.TO))
    }
    if(!is.null(num_orders)){
      R2_add_set <- matrix(NA, num_orders, nterms)
      for(ind_order in 1:num_orders){
        # ind_order<-1
        R2_condi_add <- numeric(nterms+1)
        R2_condi_add[2] <- unlist( R2.single.boot)[c(nterms*2-1+(ind_order-1)*(nterms-1))]
        R2_condi_add[-c(1:2,(nterms+1))]<-  unlist( R2.single.boot)[(2*nterms+(ind_order-1)*(nterms-1)):(2*nterms-1+nterms-2+(ind_order-1)*(nterms-1))]
        R2_condi_add[(nterms+1)]<- R2.o
        R2_add_set[ind_order,]<- R2_condi_add[-1]-R2_condi_add[-(nterms+1)]
      }
      R2_calc_result <- c(R2_calc_result, c(t(R2_add_set)))
      if(SS==TRUE){
        SS.E <- c(R2_calc_result*(-1/N*t.AK))
        result_r <- list(R2_calc_result,c(SS.E, SS.RES, SS.TO))
      }
    }
  }
  if(SS==FALSE){
    return(R2_calc_result)
  }
  if(SS==TRUE){
    return(result_r)
  }
  ## add on orders
}

#' @noRd
#' generate within cluster boot sample
WCB.sample.fun <- function(IND.matrix, weights, n){
  IND.WCB.matrix <- t(apply(IND.matrix, 1, function(x) {
    # x<- IND.WCB.matrix[1,]
    new_boot <- ddply(.data = data.frame(weights,temp_id=1:n)
                      ,.(weights),
                      function(subset){
                        # subset<-data.frame(weights,temp_id=1:n)[weights==1,]
                        ind_t <- dim(subset)[1]
                        boot_ind <- sample(1:ind_t,replace = T)
                        new_subset <- subset[boot_ind,]
                        new_subset
                      })
    # x<-(rep(new_boot$temp_id,ceiling(new_boot$w)))
    x<- new_boot$temp_id
  }))
  return(IND.WCB.matrix)
}
#' @noRd
#' generate sample by weights
AS.sample.fun <- function(IND.matrix, weights, n){
  IND.o.matrix <- t(apply(IND.matrix, 1, function(x){
    x <- sample(1:n,
                size=boot.sample.size,
                prob = weights/sum(weights),
                replace = T)
  }))
  return(IND.o.matrix)
}

#' @noRd
#' bootstrap function
boot.fun <- function(ind.boot, weights, rhs, lhs, ind.col, Ind.matrix,
                     nterms, num.list, num_orders){
  # ind.boot=1
  ind.temp.o <- Ind.matrix[ind.boot,]

  # summarize the boot sample frequencies
  ind.freq.table <- data.frame(table(ind.temp.o))
  ind.temp <- as.numeric(as.character(ind.freq.table[,1]))
  ind.weights <- ind.freq.table[,2]

  # the new weights= boot sample frequencies * sample weights
  weights.boot <- ind.weights*weights[ind.temp]

  # boot lhs
  lhs.boot <- lhs[ind.temp,ind.temp]

  # boot rhs
  rhs.boot <- as.matrix(rhs)[ind.temp,]
  delete.id<- X.check.func(rhs.boot)
  rhs.boot <- rhs.boot[,delete.id]
  ind.col.boot <- ind.col[delete.id[-1]]

  # check the singularity of boot design matrix
  invisible(tryCatch(solve(t(rhs.boot)%*%rhs.boot),
           error = function(e) {print("failed boot samples")}))
  t.AK.boot <- sum(colSums(lhs.boot*weights.boot)*weights.boot)
  boot.R2 <- R2.calc(lhs.boot, rhs.boot, weights.boot, nterms,
                     ind.col.boot, t.AK.boot, num.list, num_orders)

  return(boot.R2)
}
#' @noRd
#' @importFrom vegan getPermuteMatrix
pern<- function(permutations,n){
  p <- vegan::getPermuteMatrix(permutations, n, strata =NULL)
  return(p)
}

