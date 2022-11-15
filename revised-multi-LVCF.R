library("MASS")
library(tidyverse)
library("extraDistr")
library(nleqslv)
library(nloptr)
library(splines)
library(gaussquad)

NonHPP_gen <- function(lambdat, censor, lambda.bar) {
  nn <- 0
  while (nn <= 0) {
    nn <- rpois(1, censor * lambda.bar)
    if (nn > 0) {
      tt <- sort(runif(nn, min = 0, max = censor))
      
      ind <- (sapply(tt, lambdat) /  lambda.bar) > runif(nn)
      tt <- tt[ind]
    } else {
      tt <- numeric(0)
    }
    nn <- length(tt)
  }
  
  tt
}


simhomoPoipro <- function(rate, endTime) {
  num <- rtpois(1, rate * endTime, 0)
  return(sort(endTime * runif(num)))
}

# Simulate homogeneous poisson process
simhomoPoipro_zero <- function(rate, endTime) {
  num <- rpois(1, rate * endTime)
  return(sort(endTime * runif(num)))
}


trans_fun <- function(x, s) {
  if (s == 0) {
    res <- exp(x)
  } else {
    res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (1 / s))
  }
  return(res)
}


## old code ---------------
trans_fun_d <- function(x, s) {
  if (s == 0) {
    res <- exp(x)
  } else {
    res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (1 / s - 1))
  }
  return(res)
}

trans_fun_d1o1 <- function(x, s) {
  if (s == 0) {
    res <- rep(1, length(x))
  } else {
    res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (-1))
  }
  return(res)
}

trans_fun_d12o1 <- function(x, s) {
  if (s == 0) {
    res <- exp(x)
  } else {
    res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (1 / s - 2))
  }
  return(res)
}

trans_fun_d12o2 <- function(x, s) {
  if (s == 0) {
    res <- rep(1, length(x))
  } else {
    res = ifelse(x < -1 / s, .Machine$double.eps, (s * x + 1) ^ (-2))
  }
  return(res)
}

## new code ---------------

# trans_fun <- function(x, s) {
#   if (s == 0) {
#     res <- exp(x)
#   } else {
#     res = (s * x + 1) ^ (1 / s)
#   }
#   return(res)
# }



# trans_fun_d <- function(x, s) {
#   if (s == 0) {
#     res <- exp(x)
#   } else {
#     res =  (s * x + 1) ^ (1 / s - 1)
#   }
#   return(res)
# }

# trans_fun_d1o1 <- function(x, s) {
#   if (s == 0) {
#     res <- rep(1, length(x))
#   } else {
#     res = (s * x + 1) ^ (-1)
#   }
#   return(res)
# }

# trans_fun_d12o1 <- function(x, s) {
#   if (s == 0) {
#     res <- exp(x)
#   } else {
#     res = (s * x + 1) ^ (1 / s - 2)
#   }
#   return(res)
# }

# trans_fun_d12o2 <- function(x, s) {
#   if (s == 0) {
#     res <- rep(1, length(x))
#   } else {
#     res = (s * x + 1) ^ (-2)
#   }
#   return(res)
# }


kerfun <- function(xx){
  pmax((1-xx^2)*0.75,0)
}



simAsytransdata <-
  function(mu,
           mu_bar,
           alpha,
           beta,
           s,
           cen,
           nstep = 20) {
    p <- length(beta)
    #browser()
    # 1. Generate censoring times
    cen=runif(1,cen,1.5)
    cen <- min(cen,1)  
    
    
    
    # 2. Generate Z(t) and h(t) as a step function
    
    Sigmat_z <- exp(-abs(outer(1:nstep, 1:nstep, "-")) / nstep)
    
    z <-  2*(pnorm(c(mvrnorm(
      1, rep(0, nstep), Sigmat_z
    ))) - 0.5)
    
    left_time_points <- (0:(nstep - 1)) / nstep
    
    # 3. Generate observation time
    z_fun <- stepfun(left_time_points, c(0, z))
    # z_fun=function(x){
    #   z_star(x)
    # }
    z_cons<- beta[2] * rbernoulli(1)
    h_fun <- function(x) {
      beta[1] %*% z_fun(x) + 
        z_cons
    }
    
    # 4. Generate failure time
    if (s == 0) {
      lam_fun <- function(tt)
        exp(alpha(tt) + h_fun(tt))
    } else {
      lam_fun <- function(tt)
        (s * (alpha(tt) + h_fun(tt)) + 1) ^ (1 / s)
    }
    u <- runif(1)
    
    fail_time <- nleqslv(cen/2 , function(ttt)
      legendre.quadrature(lam_fun,
                          lower = 0,
                          upper = ttt,
                          lqrule64) + log(u))$x
    
    
    
    X <- min(fail_time, cen)
    
    # 1. Generate R_ik
    
    obs_times <-
      NonHPP_gen(mu,
                 cen, mu_bar)
    
    if (length(obs_times) == 0)
      obs_times <- cen
    
    covariates_obscov <- matrix(c(z_fun(obs_times),rep((z_cons/beta[2]),length(z_fun(obs_times)))),length(z_fun(obs_times)))
    
    # Results
    return(
      tibble(
        X = X,
        delta = fail_time < cen,
        covariates = covariates_obscov,
        obs_times = obs_times,
        censoring = cen
      ) # %>% filter(X >= obs_times)
    )
  }

# lqrule <- legendre.quadrature.rules(256)[[256]]
lqrule64 <- legendre.quadrature.rules(64)[[64]]

estproc_ori_d_mul <- function(data, n, nknots, norder, s, h, pl = 0) {
  #gammap <- nknots + norder+1
  gammap<- nknots+norder-1
  X <- data$X
  id <- data$id
  covariates <- data$covariates 
  obs_times <- data$obs_times
  delta <- data$delta
  kerval <- kerfun((X - obs_times) / h) / h * (X > obs_times)
  #kerval <- ifelse(kerval ==0,1,kerval)
  knots <- (1:nknots) / (nknots + 1)
  #knots <- data %>% filter(delta & (abs(X - obs_times) < h)) %>% pull(X) %>% quantile((1:nknots) / (nknots + 1))
  
  bsmat <-
    # bs(
    ns(
      X,
      knots = knots,
      # degree = 4,
      intercept = TRUE,
      Boundary.knots = c(0, 1)
    )
  
  loglik1_1 <- function(beta, gamma) {
    alphaX <- bsmat %*% gamma
    res <-
      sum(log(trans_fun(alphaX +  covariates %*% beta, s)) * delta * kerval)
    
    res
  }
  
  loglik1_1_d <- function(beta, gamma) {
    alphaBeta <- bsmat %*% gamma +covariates %*% beta
    temp1 <- trans_fun_d1o1(alphaBeta, s) * delta * kerval
    res <- as.vector(t(temp1) %*% cbind(covariates, bsmat))
    
    res
  }
  
  loglik2_inner_1 <- function(tt, beta, gamma) {
    #browser()
    dist <- outer(tt, obs_times, "-")
    kerval_tt <- kerfun(dist / h) / h * (dist > 0)
    alpha_tt <-
      # bs(
      ns(  tt,
           knots = knots,
           # degree = 4,
           intercept = TRUE,
           Boundary.knots = c(0, 1)  ) %*% gamma %>% as.vector()
    res <-
      trans_fun(outer(alpha_tt,c(covariates %*% beta), "+"), s) 
    res <- ifelse(kerval_tt==0,0,res)* kerval_tt
    res <- ifelse(outer(tt, X, "<"),res,0)
    rowSums(res) #/ rowSums(kerval_tt)
  }
  
  loglik2_inner_1_d <- function(tt, beta, gamma) {
    
    dist <- outer(tt, obs_times, "-")
    kerval_tt <- kerfun(dist / h) / h * (dist > 0)
    bsmat_tt <-
      #  bs(
      ns(
        tt,
        knots = knots,
        # degree = 4,
        intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
    
    alpha_tt <-
      bsmat_tt %*% gamma %>% as.vector()
    temp1 <- trans_fun_d(outer(alpha_tt,c(covariates %*% beta), "+"), s)
    temp1 <- ifelse(outer(tt, X, "<"),temp1,0)
    temp1 <- ifelse(kerval_tt==0,0,temp1) * kerval_tt
    res <- cbind(temp1 %*% covariates, rowSums(temp1) * bsmat_tt)
    
    res
  }
  
  loglik_1 <- function(beta, gamma) {
    res <-
      loglik1_1(beta, gamma) -
      legendre.quadrature(
        function(tt)
          loglik2_inner_1(tt, beta, gamma),
        lower = 0,
        upper = 1,
        lqrule64
      )
    
    res / n
    
  }
  
  loglik_1_d <- function(beta, gamma) {
    res1 <-
      loglik1_1_d(beta, gamma)
    res2 <- loglik2_inner_1_d(0.5 * lqrule64$x + 0.5, beta, gamma)
    res2 <- as.vector(0.5 * colSums(lqrule64$w * res2))
    (res1 - res2) / n
    
  }
  
  
  f <- function(xx) {
    -loglik_1(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))])
  }
  # estres <- optim(rep(1, 1 + gammap),f  )
  # browser()
  # estresl <- lbfgs(f, function(xx) - loglik_1_d(xx[1], xx[-1], s), rep(1, 1 + gammap))
  
  
  A <- function(beta, gamma) {
    
    # Get Z_bar
    dist_XX <- outer(X, obs_times, "-")
    kerval_XX <- kerfun(dist_XX / h) / h * (dist_XX > 0)
    bsmat_XX <-
      #bs(
      ns(
        X,
        knots = knots,
        # degree = 4,
        intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
    
    alpha_XX <-
      bsmat_XX %*% gamma %>% as.vector()
    inner <- outer(alpha_XX, covariates %*% beta, "+")
    S0 <- outer(X, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(X)) * kerval_XX
    S1 <- S0 %*% covariates/ rowSums(S0)
    outf=function(x){x %o% x}
    S2 <- S0 %*% ( t(apply(covariates,1,outf))) / rowSums(S0)
    outerprod <-  (S2 - t(apply(S1,1,outf))) * c(trans_fun_d1o1(alpha_XX +  covariates %*% beta, s)
    ) ^ 2 
    outerprod[is.na(outerprod)] <-0
    matrix(colSums( outerprod *  kerval * delta) / n,ncol(data$covariates))
    
  }
  
  B <- function(beta, gamma) {
    
    dist_XX <- outer(X, obs_times, "-")
    kerval_XX <- kerfun(dist_XX / h) / h * (dist_XX > 0)
    bsmat_XX <-
      #bs(
      ns(
        X,
        knots = knots,
        # degree = 4,
        intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
    
    alpha_XX <-
      bsmat_XX %*% gamma %>% as.vector()
    inner <- outer(alpha_XX,   covariates %*% beta, "+")
    S0 <- outer(X, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(X)) * kerval_XX
    S1 <- S0 %*% covariates / rowSums(S0)
    outf=function(x){x %o% x}
    outerproducinner1 <- (S1 - covariates) * c(trans_fun_d1o1(alpha_XX +  covariates %*% beta, s))
    outerproducinner1[is.na(outerproducinner1)] <-0
    bb1=outerproducinner1*  kerval * delta * ( (X-obs_times)>0 )
    bb1=tibble(a=bb1 ,id=as.numeric(id)) %>% group_by(id) %>%  summarise( colSums(a))
    nc=ncol(data$covariates)
    bbvec=NULL
    for(kk in 1:n)
    {
      bbvec=c(bbvec,c(bb1%>% filter(row_number() == kk))$`colSums(a)`)
    }
    bb1=matrix(bbvec,n)
    
    b_inner <- function(tt, beta, gamma) {
      dist <- outer(tt, obs_times, "-")
      kerval_tt <- kerfun(dist / h) / h * (dist > 0)
      bsmat_tt <-
        #bs(
        ns(
          tt,
          knots = knots,
          # degree = 4,
          intercept = TRUE,
          Boundary.knots = c(0, 1)
        )
      
      
      
      
      
      alpha_tt <-
        bsmat_tt %*% gamma %>% as.vector()
      inner <- outer(alpha_tt,   covariates %*% beta  , "+")
      S0_tt <- outer(tt, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(tt)) * kerval_tt
      S1_tt <- S0_tt %*% covariates / rowSums(S0_tt)
      res=NULL
      
      for( jj in 1:nc){
        temp1  <-
          outer(tt, X, "<=") * matrix(trans_fun_d(outer(alpha_tt, covariates %*% beta,"+"), s),length(tt))*
          kerval_tt     *outer(as.vector(S1_tt[,jj]),covariates[,jj],"-") 
        temp1[is.na(temp1)] <-0
        r1=tibble(a=t(temp1) ,id=as.numeric(id)) %>% group_by(id) %>%  summarise( colMeans(a))
        mr1=NULL
        for(i in 1:length(tt)){
          mr1 =c(mr1,c(r1%>% filter(row_number() == i))$`colMeans(a)`)
        }
        r1=matrix( mr1,n)
        res=cbind(res,r1)
      }
      res
    }
    
    binner <- b_inner (0.5 * lqrule64$x + 0.5, beta, gamma)
    bb2=NULL
    for( jj in 1:nc){
      b1 <- as.vector(0.5 * rowSums(lqrule64$w * binner[,(1+(nc-1)*64):(64*nc)]))
      bb2=cbind(bb2,b1)
    }
    bb=bb1-bb2
    matrix(colSums(t(apply(bb,1,outf))) / n,ncol(data$covariates))
  }
  
  # estres <- optim(rep(1, 1 + gammap),f  )
  #browser()
  # estresl <- lbfgs(f, function(xx) - loglik_1_d(xx[1], xx[-1], s), rep(1, 1 + gammap))
  # estres <- nloptr(
  #   x0 = rep(0, ncol(data$covariates) + gammap),
  #   eval_f = f,
  #   eval_grad_f = function(xx)
  #     - loglik_1_d(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))]),
  #   opts = list(
  #     "algorithm" = "NLOPT_LD_LBFGS",
  #     "maxeval" = 10000,
  #     "print_level" = 0
  #   )
  # )
  
  if(s==0) {
    estres <- nloptr(
      x0 = rep(0, ncol(data$covariates) + gammap),
      eval_f = f,
      eval_grad_f = function(xx)
        - loglik_1_d(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))]),
      opts = list(
        "algorithm" = "NLOPT_LD_SLSQP",
        "xtol_rel"=1.0e-6,
        "maxeval" = 10000,
        "print_level" = 0
      )
    )
  } else {
    #ineqmat <- as.matrix(expand_grid(covariates,unique(bsmat)))
    ineqmat <- cbind(covariates,bsmat)
    ineqmat <- ineqmat[(X>obs_times)& (abs(X-obs_times) <= h),]
    estres <- nloptr(
      x0 = rep(0, ncol(data$covariates) + gammap),
      eval_f = f,
      eval_grad_f = function(xx)
        - loglik_1_d(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))]),
      eval_g_ineq = function(xx) -ineqmat %*%xx-1/s+1e-6,
      eval_jac_g_ineq = function(xx) -ineqmat,
      opts = list(
        "algorithm" = "NLOPT_LD_SLSQP",
        "xtol_rel"=1.0e-6,
        "maxeval" = 10000,
        "print_level" = 0
      )
    )
  }
  
  
  # Sig_est <- Sigma1(estres$solution[1], estres$solution[-1])
  A_est <- A(estres$solution[1:ncol(data$covariates)], estres$solution[-(1:ncol(data$covariates))])
  B_est <- B(estres$solution[1:ncol(data$covariates)], estres$solution[-(1:ncol(data$covariates))])
  # browser()
  # list(est=estres$solution,var=solve(A_est) %*% B_est %*% solve(A_est)/n,A_est=A_est,B_est=B_est)
  
  se=sqrt(diag(solve(A_est) %*% B_est %*% solve(A_est)/n))
  list(est=estres$solution,se=se,A_est=A_est,B_est=B_est)
  
}
## Comment code ---------------------
# estproc_ori_d1 <- function(data, n, nknots, norder, s, h, pl = 0) {
#   gammap <- nknots + norder
#   X <- data$X
#   id <- data$idd
#   covariates <- data$covariates
#   obs_times <- data$obs_times
#   delta <- data$delta
#   kerval <- kerfun((X - obs_times) / h) / h * (X > obs_times)
#   knots <- (1:nknots) / (nknots + 1)
#   bsmat <-
#     bs(
#       X,
#       knots = knots,
#       degree = norder,
#       Boundary.knots = c(0, 1)
#     )
#   
#   loglik1_1 <- function(beta, gamma) {
#     alphaX <- bsmat %*% gamma
#     res <-
#       sum(log(trans_fun(alphaX + beta * covariates, s)) * delta * kerval)
#     # if (is.nan(res))
#     #   browser()
#     res
#   }
#   
#   loglik1_1_d <- function(beta, gamma) {
#     alphaBeta <- bsmat %*% gamma + beta * covariates
#     temp1 <- trans_fun_d1o1(alphaBeta, s) * delta * kerval
#     res <- as.vector(t(temp1) %*% cbind(covariates, bsmat))
#     res
#   }
#   
#   loglik2_inner_1 <- function(tt, beta, gamma) {
#     dist <- outer(tt, obs_times, "-")
#     kerval_tt <- kerfun(dist / h) / h * (dist > 0)
#     alpha_tt <-
#       bs(
#         tt,
#         knots = knots,
#         degree = norder,
#         Boundary.knots = c(0, 1)
#       ) %*% gamma %>% as.vector()
#     res <-
#       trans_fun(outer(alpha_tt, beta * covariates, "+"), s) * kerval_tt
#     rowSums(outer(tt, X, "<") * res) #/ rowSums(kerval_tt)
#   }
#   
#   loglik2_inner_1_d <- function(tt, beta, gamma) {
#     dist <- outer(tt, obs_times, "-")
#     kerval_tt <- kerfun(dist / h) / h * (dist > 0)
#     bsmat_tt <-
#       bs(
#         tt,
#         knots = knots,
#         degree = norder,
#         Boundary.knots = c(0, 1)
#       )
#     alpha_tt <-
#       bsmat_tt %*% gamma %>% as.vector()
#     temp1 <-
#       outer(tt, X, "<") * trans_fun_d(outer(alpha_tt, beta * covariates, "+"), s) * kerval_tt
#     
#     
#     
#     res <- cbind(temp1 %*% covariates, rowSums(temp1) * bsmat_tt)
#     res
#   }
#   
#   loglik_1 <- function(beta, gamma) {
#     res <-
#       loglik1_1(beta, gamma) -
#       legendre.quadrature(
#         function(tt)
#           loglik2_inner_1(tt, beta, gamma),
#         lower = 0,
#         upper = 1,
#         lqrule64
#       )
#     #if (is.nan(res))
#     # browser()
#     res / n
#     
#   }
#   
#   loglik_1_d <- function(beta, gamma) {
#     res1 <-
#       loglik1_1_d(beta, gamma)
#     res2 <- loglik2_inner_1_d(0.5 * lqrule64$x + 0.5, beta, gamma)
#     res2 <- as.vector(0.5 * colSums(lqrule64$w * res2))
#     (res1 - res2) / n
#     
#   }
#   
#   
#   f <- function(xx) {
#     -loglik_1(xx[1], xx[-1])
#   }
#   
#   estres <- nloptr(
#     x0 = rep(1, 1 + gammap),
#     eval_f = f,
#     eval_grad_f = function(xx)
#       - loglik_1_d(xx[1], xx[-1]),
#     opts = list(
#       "algorithm" = "NLOPT_LD_LBFGS",
#       "maxeval" = 10000,
#       "print_level" = 0
#     )
#   )
#   
#   estres$solution
#   
# }

# boot<-function( simdata,
#                 n,
#                 nknots,
#                 norder,
#                 s,
#                 h ){
#   
#   bootn=100
#   id=simdata$id
#   res=rep(0,bootn)
#   for  (i in 1:bootn){
#     data2 <-
#       replicate(
#         n,
#         simdata[which(id==sample(1:n,size=1)),],
#         simplify = F
#       ) %>% bind_rows(.id = "idd")
#     res[i]=estproc_ori_d1(
#       data2,
#       n,
#       nknots = 3,
#       norder = 3,
#       s,
#       h = n ^ (-0.7)
#     )[1]
#   }
#   return(sqrt(var(res)))
#   
# }

# CVgroup <- function(k,datasize,seed){
#   cvlist <- list()
#   set.seed(seed)
#   n <- rep(1:k,ceiling(datasize/k))[1:datasize]     
#   temp <- sample(n,datasize)    
#   x <- 1:k
#   dataseq <- 1:datasize
#   cvlist <- lapply(x,function(x) dataseq[temp==x])   
#   return(cvlist)
# }


# hcv<-function(simdata,n,K,s){
#   # browser()
#   test_idk <- lapply(split(sample(1:n,n),rep(1:K,n/K)),sort)
#   
#   nn=11
#   hdiff <- simdata %>% mutate(diff=X-obs_times) %>% filter(diff>0) %>% group_by(id) %>% mutate(diff1=min(diff)) %>% pull(diff1)
#   #hmin <- min(hdiff[hdiff!=min(hdiff)]) 
#   hmin <- simdataout %>% filter(delta) %>% filter(X==max(X)) %>% mutate(a=X-obs_times) %>% filter(a>0) %>% filter(a==min(a)) %>% pull(a)
#   hmin<-max(hmin,n^(-0.6))
#   hmax <- max(simdata %>% mutate(diff=X-obs_times) %>% filter(diff>0) %>% group_by(id) %>% mutate(diff1=max(diff)) %>% pull(diff1))
#   hmax= min(hmax,n^(-0.3)) 
#   hn <- exp(seq(log(hmin),log(hmax),length.out=nn+1)[-1])
#   
#   browser()
#   res <- foreach(hh=hn) %do% {
#     foreach(test_idx_one = test_idk) %do% {
#       foldkpar=estproc_ori_dh(simdata %>% 
#                                 filter(!(as.numeric(id) %in% test_idx_one)),n*(1-1/K), 3, 3, s, hh, pl = 0) 
#       testing <- simdata %>% filter(as.numeric(id) %in% test_idx_one)
#       DIres <- testing %>% mutate(kerval = kerfun((X - obs_times) / hh) / hh * (X >= obs_times),
#                                   H = trans_fun(ns(X, knots =  (1:3) / (3 + 1), intercept = TRUE, Boundary.knots = c(0, 1)) %*% foldkpar[-(1:2)] 
#                                                 + covariates %*% foldkpar[1:2], s)) %>% 
#         group_by(id) %>% 
#         summarise(dMI =  sum(kerval*(delta-H))/sum(kerval),
#                   delta=mean(delta),
#                   .groups = "drop") %>% 
#         mutate(DI = ifelse(is.na(dMI),0,sign(dMI)*sqrt(-2*(dMI+delta*log(delta-dMI)))))
#       DIres %>% mutate(bd=hh)
#     } %>% bind_rows(.id="fold")
#   }
#   #browser()
#   res %>% bind_rows() %>% group_by(bd) %>% summarise(cvloss=sum(DI^2)/n)
#   
# }


## unComment code ---------------------

estproc_LVCF_mul <- function(data, n, nknots, norder, s, h, initpara=NULL) {
  gammap<- nknots+norder-1
  X <- data$X
  id <- data$id
  covariates <- data$covariates 
  obs_times <- data$obs_times
  delta <- data$delta
  #kerval <- kerfun((X - obs_times) / h) / h * (X > obs_times)
  kervalinner=X - obs_times
  kervalinner[which(kervalinner<0)]=100
  kerval <- unsplit(sapply(split(kervalinner,as.factor(id)),
                           function(xx) as.numeric((xx==min(xx)) &(xx < Inf) )),
                    factor(simdata$id))
  
  #kerval <- ifelse(kerval ==0,1,kerval)
  knots <- (1:nknots) / (nknots + 1)
  #knots <- data %>% filter(delta & (abs(X - obs_times) < h)) %>% pull(X) %>% quantile((1:nknots) / (nknots + 1))
  
  bsmat <-
    # bs(
    ns(
      X,
      knots = knots,
      # degree = 4,
      intercept = TRUE,
      Boundary.knots = c(0, 1)
    )
  
  loglik1_1 <- function(beta, gamma) {
    alphaX <- bsmat %*% gamma
    res <-
      sum(log(trans_fun(alphaX +  covariates %*% beta, s)) * delta * kerval)
    
    res
  }
  
  loglik1_1_d <- function(beta, gamma) {
    alphaBeta <- bsmat %*% gamma +covariates %*% beta
    temp1 <- trans_fun_d1o1(alphaBeta, s) * delta * kerval
    res <- as.vector(t(temp1) %*% cbind(covariates, bsmat))
    
    res
  }
  ismin<-function(a){
    a[which(a<0)]=100
    t=ifelse((abs(a)==min(abs(a)))&(a<99),1,0)
    t[1]=ifelse(sum(t)==0,1,t[1])
    t
  }
  ismin_mat= function(a){
    if(dim(a)[1]==1){matrix(rep(1,dim(a)[2]),1)}else{apply(a, 2,ismin)}
  }
  loglik2_inner_1 <- function(tt, beta, gamma) {
    #browser()
    #dist <- outer(tt, obs_times, "-")
    #kerval_tt <- kerfun(dist / h) / h * (dist > 0)
    
    kerval_ttinner=outer(tt, obs_times, "-")
    kerval_ttinner[which(kerval_ttinner<0)]=100
    kerval_tt <- t(tibble(id=as.numeric(id),dist=t(kerval_ttinner))%>%group_by(id)%>%
                     summarise(indicator= ismin_mat(dist),.groups = "keep")%>%pull(indicator))  
    
  
    alpha_tt <-
      # bs(
      ns(  tt,
           knots = knots,
           # degree = 4,
           intercept = TRUE,
           Boundary.knots = c(0, 1)  ) %*% gamma %>% as.vector()
    res <-
      trans_fun(outer(alpha_tt,c(covariates %*% beta), "+"), s) 
    res <- ifelse(kerval_tt==0,0,res)* kerval_tt
    res <- ifelse(outer(tt, X, "<"),res,0)
    rowSums(res) #/ rowSums(kerval_tt)
  }
  
  loglik2_inner_1_d <- function(tt, beta, gamma) {
    
    #dist <- outer(tt, obs_times, "-")
    #kerval_tt <- kerfun(dist / h) / h * (dist > 0)
    
    kerval_ttinner=outer(tt, obs_times, "-")
    kerval_ttinner[which(kerval_ttinner<0)]=100
    kerval_tt <- t(tibble(id=as.numeric(id),dist=t(kerval_ttinner))%>%group_by(id)%>%
                     summarise(indicator= ismin_mat(dist),.groups = "keep")%>%pull(indicator))  
    
    
    bsmat_tt <-
      #  bs(
      ns(
        tt,
        knots = knots,
        # degree = 4,
        intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
    
    alpha_tt <-
      bsmat_tt %*% gamma %>% as.vector()
    temp1 <- trans_fun_d(outer(alpha_tt,c(covariates %*% beta), "+"), s)
    temp1 <- ifelse(outer(tt, X, "<"),temp1,0)
    temp1 <- ifelse(kerval_tt==0,0,temp1) * kerval_tt
    #browser()
    res <- cbind(temp1 %*% covariates, rowSums(temp1) * bsmat_tt)
    
    res
  }
  
  loglik_1 <- function(beta, gamma) {
    res <-
      loglik1_1(beta, gamma) -
      legendre.quadrature(
        function(tt)
          loglik2_inner_1(tt, beta, gamma),
        lower = 0,
        upper = 1,
        lqrule64
      )
    
    res / n
    
  }
  
  loglik_1_d <- function(beta, gamma) {
    res1 <-
      loglik1_1_d(beta, gamma)
    res2 <- loglik2_inner_1_d(0.5 * lqrule64$x + 0.5, beta, gamma)
    res2 <- as.vector(0.5 * colSums(lqrule64$w * res2))
    (res1 - res2) / n
    
  }
  
  
  f <- function(xx) {
    -loglik_1(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))])
  }
  # estres <- optim(rep(1, 1 + gammap),f  )
  # browser()
  # estresl <- lbfgs(f, function(xx) - loglik_1_d(xx[1], xx[-1], s), rep(1, 1 + gammap))
  
  
  A <- function(beta, gamma) {
    
    # Get Z_bar
    kerval_XXinner=outer(X, obs_times, "-")
    kerval_XXinner[which(kerval_XXinner<0)]=100
    kerval_XX <- t(tibble(id=as.numeric(id),dist=t(kerval_XXinner))%>%group_by(id)%>%
                     summarise(indicator= ismin_mat(dist),.groups = "keep")%>%pull(indicator))  

    bsmat_XX <-
      #bs(
      ns(
        X,
        knots = knots,
        # degree = 4,
        intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
    
    alpha_XX <-
      bsmat_XX %*% gamma %>% as.vector()
    inner <- outer(alpha_XX, covariates %*% beta, "+")
    S0 <- outer(X, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(X)) * kerval_XX
    S1 <- S0 %*% covariates/ rowSums(S0)
    outf=function(x){x %o% x}
    S2 <- S0 %*% ( t(apply(covariates,1,outf))) / rowSums(S0)
    outerprod <-  (S2 - t(apply(S1,1,outf))) * c(trans_fun_d1o1(alpha_XX +  covariates %*% beta, s)
    ) ^ 2 
    outerprod[is.na(outerprod)] <-0
    matrix(colSums( outerprod *  kerval * delta) / n,ncol(data$covariates))
    
  }
  
  B <- function(beta, gamma) {
    
    # dist_XX <- outer(X, obs_times, "-")
    # kerval_XX <- kerfun(dist_XX / h) / h * (dist_XX > 0)
    kerval_XXinner=outer(X, obs_times, "-")
    kerval_XXinner[which(kerval_XXinner<0)]=100
    kerval_XX <- t(tibble(id=as.numeric(id),dist=t(kerval_XXinner))%>%group_by(id)%>%
                     summarise(indicator= ismin_mat(dist),.groups = "keep")%>%pull(indicator))  
     
    
    bsmat_XX <-
      #bs(
      ns(
        X,
        knots = knots,
        # degree = 4,
        intercept = TRUE,
        Boundary.knots = c(0, 1)
      )
    
    alpha_XX <-
      bsmat_XX %*% gamma %>% as.vector()
    inner <- outer(alpha_XX,   covariates %*% beta, "+")
    S0 <- outer(X, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(X)) * kerval_XX
    S1 <- S0 %*% covariates / rowSums(S0)
    outf=function(x){x %o% x}
    outerproducinner1 <- (S1 - covariates) * c(trans_fun_d1o1(alpha_XX +  covariates %*% beta, s))
    outerproducinner1[is.na(outerproducinner1)] <-0
    bb1=outerproducinner1*  kerval * delta * ( (X-obs_times)>0 )
    bb1=tibble(a=bb1 ,id=as.numeric(id)) %>% group_by(id) %>%  summarise( colSums(a))
    nc=ncol(data$covariates)
    bbvec=NULL
    for(kk in 1:n)
    {
      bbvec=c(bbvec,c(bb1%>% filter(row_number() == kk))$`colSums(a)`)
    }
    bb1=matrix(bbvec,n)
    
    b_inner <- function(tt, beta, gamma) {
      #dist <- outer(tt, obs_times, "-")
      #kerval_tt <- kerfun(dist / h) / h * (dist > 0)
      kerval_ttinner=outer(tt, obs_times, "-")
      kerval_ttinner[which(kerval_ttinner<0)]=100
      kerval_tt <- t(tibble(id=as.numeric(id),dist=t(kerval_ttinner))%>%group_by(id)%>%
                       summarise(indicator= ismin_mat(dist),.groups = "keep")%>%pull(indicator))  
      
      
      
      bsmat_tt <-
        #bs(
        ns(
          tt,
          knots = knots,
          # degree = 4,
          intercept = TRUE,
          Boundary.knots = c(0, 1)
        )
      
      
      
      
      
      alpha_tt <-
        bsmat_tt %*% gamma %>% as.vector()
      inner <- outer(alpha_tt,   covariates %*% beta  , "+")
      S0_tt <- outer(tt, X, "<=") * matrix(trans_fun_d12o1(inner, s),length(tt)) * kerval_tt
      S1_tt <- S0_tt %*% covariates / rowSums(S0_tt)
      res=NULL
      
      for( jj in 1:nc){
        temp1  <-
          outer(tt, X, "<=") * matrix(trans_fun_d(outer(alpha_tt, covariates %*% beta,"+"), s),length(tt))*
          kerval_tt     *outer(as.vector(S1_tt[,jj]),covariates[,jj],"-") 
        temp1[is.na(temp1)] <-0
        r1=tibble(a=t(temp1) ,id=as.numeric(id)) %>% group_by(id) %>%  summarise( colMeans(a))
        mr1=NULL
        for(i in 1:length(tt)){
          mr1 =c(mr1,c(r1%>% filter(row_number() == i))$`colMeans(a)`)
        }
        r1=matrix( mr1,n)
        res=cbind(res,r1)
      }
      res
    }
    
    binner <- b_inner (0.5 * lqrule64$x + 0.5, beta, gamma)
    bb2=NULL
    for( jj in 1:nc){
      b1 <- as.vector(0.5 * rowSums(lqrule64$w * binner[,(1+(nc-1)*64):(64*nc)]))
      bb2=cbind(bb2,b1)
    }
    bb=bb1-bb2
    matrix(colSums(t(apply(bb,1,outf))) / n,ncol(data$covariates))
  }
  
 
  if(is.null(initpara)) initpara <- rep(0, ncol(data$covariates) + gammap)
  if(s==0) {
    estres <- nloptr(
      x0 = initpara,
      eval_f = f,
      eval_grad_f = function(xx)
        - loglik_1_d(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))]),
      opts = list(
        "algorithm" = "NLOPT_LD_SLSQP",
        "xtol_rel"=1.0e-6,
        "maxeval" = 10000,
        "print_level" = 0
      )
    )
  } else {
    #ineqmat <- as.matrix(expand_grid(covariates,unique(bsmat)))
    ineqmat <- cbind(covariates,bsmat)
    ineqmat <- ineqmat[(X>obs_times)& (abs(X-obs_times) <= h),]
    estres <- nloptr(
      x0 = initpara,
      eval_f = f,
      eval_grad_f = function(xx)
        - loglik_1_d(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))]),
      eval_g_ineq = function(xx) -ineqmat %*%xx-1/s+1e-8,
      eval_jac_g_ineq = function(xx) -ineqmat,
      opts = list(
        "algorithm" = "NLOPT_LD_SLSQP",
        "xtol_rel"=1.0e-8,
        "maxeval" = 10000,
        "print_level" = 0
      )
    )
  }
  
  A_est <- A(estres$solution[1:ncol(data$covariates)], estres$solution[-(1:ncol(data$covariates))])
  B_est <- B(estres$solution[1:ncol(data$covariates)], estres$solution[-(1:ncol(data$covariates))])
  
  se=sqrt(diag(solve(A_est) %*% B_est %*% solve(A_est)/n))
  list(est=estres$solution,se=se,A_est=A_est,B_est=B_est)
}





hcv<-function(simdata,n,K,s){
  # browser()
  test_idk <- lapply(split(sample(1:n,n),rep(1:K,n/K)),sort)
  
  nn=11
  #hdiff <- simdata %>% mutate(diff=X-obs_times) %>% filter(diff>0) %>% group_by(id) %>% mutate(diff1=min(diff)) %>% pull(diff1)
  hmin <- simdata %>% filter(delta) %>% filter(X>0.75) %>% mutate(a=X-obs_times) %>% filter(a>0) %>% filter(a==min(a)) %>% pull(a)
  hmin<-max(hmin,n^(-0.6))
  hmax <- max(simdata %>% mutate(diff=X-obs_times) %>% filter(diff>0) %>% group_by(id) %>% mutate(diff1=max(diff)) %>% pull(diff1))
  hmax= min(hmax,n^(-0.3)) 
  hn <- exp(seq(log(hmin),log(hmax),length.out=nn+1)[-1])
  
  # browser()
  res <- foreach(hh=hn) %do% {
    foreach(test_idx_one = test_idk) %do% {
      foldkpar=estproc_ori_dh(simdata %>% 
                                filter(!(as.numeric(id) %in% test_idx_one)),n*(1-1/K), 3, 3, s, hh) 
      testing <- simdata %>% filter(as.numeric(id) %in% test_idx_one)
      test_logll <- logll_val(foldkpar,testing,n/K, 3, 3, s, hh, pl = 0)
      tibble(test_logll,bd=hh)
    } %>% bind_rows(.id="fold")
  }
  #browser()
  res %>% bind_rows() %>% group_by(bd) %>% summarize(cvloss=mean(test_logll))
  
}

#### Calculate loglikelihood


logll_val <- function(para,data, n, nknots, norder, s, h, pl = 0) {
  #gammap <- nknots + norder+1
  gammap<- nknots+norder-1
  X <- data$X
  id <- data$id
  covariates <- data$covariates 
  obs_times <- data$obs_times
  delta <- data$delta
  kerval <- kerfun((X - obs_times) / h) / h * (X > obs_times)
  #kerval <- ifelse(kerval ==0,1,kerval)
  knots <- (1:nknots) / (nknots + 1)
  
  bsmat <-
    # bs(
    ns(
      X,
      knots = knots,
      # degree = 4,
      intercept = TRUE,
      Boundary.knots = c(0, 1)
    )
  
  loglik1_1 <- function(beta, gamma) {
    alphaX <- bsmat %*% gamma
    res <-
      sum(log(trans_fun(alphaX +  covariates %*% beta, s)) * delta * kerval,na.rm = T)
    
    res
  }
  
  
  loglik2_inner_1 <- function(tt, beta, gamma) {
    #browser()
    dist <- outer(tt, obs_times, "-")
    kerval_tt <- kerfun(dist / h) / h * (dist > 0)
    alpha_tt <-
      # bs(
      ns(  tt,
           knots = knots,
           # degree = 4,
           intercept = TRUE,
           Boundary.knots = c(0, 1)  ) %*% gamma %>% as.vector()
    res <-
      trans_fun(outer(alpha_tt,c(covariates %*% beta), "+"), s) * kerval_tt
    rowSums(outer(tt, X, "<") * res) #/ rowSums(kerval_tt)
  }
  
  
  loglik_1 <- function(beta, gamma) {
    res <-
      loglik1_1(beta, gamma) -
      legendre.quadrature(
        function(tt)
          loglik2_inner_1(tt, beta, gamma),
        lower = 0,
        upper = 1,
        lqrule64
      )
    
    res / n
    
  }
  
  f <- function(xx) {
    -loglik_1(xx[1:ncol(data$covariates)], xx[-(1:ncol(data$covariates))])
  }
  # estres <- optim(rep(1, 1 + gammap),f  )
  # browser()
  # estresl <- lbfgs(f, function(xx) - loglik_1_d(xx[1], xx[-1], s), rep(1, 1 + gammap))
  
  
  
  # browser()
  # list(est=estres$solution,var=solve(A_est) %*% B_est %*% solve(A_est)/n,A_est=A_est,B_est=B_est)
  f(para) 
  
}




library(tidyverse)
library(caret)
library(foreach)
library(doParallel)
cores=detectCores()
cl <- makeCluster(cores)
doParallel::registerDoParallel(cl)

s_vec <- c(0,0.25,0.5,0.75,1)
cenrate_vec <- c(20,30)
# cen_mat <- matrix(c(0.789,0.0526,0.789,0.158,0.895,0.105,0.842,0.0526,0.895,0.0526),nrow=2)

cen_mat <- matrix(c(0.39807293, 0.4714146, 0.5447563, 0.6180980, 0.6914397,
                    0.07700907, 0.1503508, 0.2236924, 0.2970341, 0.3703758),nrow=2,byrow=T)

true_coef <- c(1,-0.5)

res <-
  foreach(s_ind = c(1)) %do% {
    foreach(cenrate_ind = 1:2) %do% {
      foreach(nn = c(200, 400)) %do% {
        cen_val <- cen_mat[cenrate_ind, s_ind]
        cenrate <- cenrate_vec[cenrate_ind]
        s_val <- s_vec[s_ind]
        print(c(nn, s_val, cenrate))
        print(Sys.time())
        estres_all_h <-
          foreach(
            ii = 1:1000,
            .packages = c(
              "tidyverse",
              "extraDistr",
              "MASS",
              "nleqslv",
              "nloptr",
              "splines",
              "gaussquad",
              "foreach",
              "survival"
            )
          ) %do% {
            set.seed(ii)
            print(ii)
            
            simdata <-
              replicate(
                nn,
                simAsytransdata(
                  mu = function(tt)
                    8 * (0.75 + (0.5 - tt) ^ 2),
                  # 5,
                  mu_bar = 8,
                  #rep(10,length(tt)),
                  alpha = function(tt)
                    s_val + (tt * (1 - sin(
                      2 * pi * (tt - 0.25)
                    ))),
                  #s_val+(tt*(1-sin(2*pi*(tt-0.25))))/2,
                  #0.1+0.5*s_val+(tt*(1-sin(2*pi*(tt-0.25)))),
                  #function(tt) 0.5*s_val+tt*sin(2*pi*(tt-0.25))/2,
                  #(tt*(1-sin(2*pi*(tt-0.25)))+0.2+0.3*s_val),
                  beta = true_coef ,
                  s = s_val,
                  cen = cen_val  ,
                  nstep = 20
                ),
                simplify = F
              ) %>% bind_rows(.id = "id")
            # simdata %>% group_by(id) %>% filter(row_number() == 1)   %>% pull(delta) %>% mean()
            # unlist(estres_ori)
            # print(simdata)
            #  hres <- hcv(simdata,n,5,s)
            # browser()
            # print(t(res))
            
            estres_ori <- estproc_LVCF_mul(
              simdata,
              nn,
              nknots = 3,
              norder = 3,
              s_val,
              0.1)
            
            #browser()
            list(estres_ori,
                 n = nn,
                 s = s_val,
                 cenrate = cenrate)
            #browser()
            # list(hres,estres_ori)
            #list(estres_ori)
          }
        
        # filename <-
        #   paste(
        #     "/home/dsun30/asyntrans/simres/newLVCF6multi",
        #     "s",
        #     s_val,
        #     "n",
        #     nn,
        #     "cen",
        #     cenrate,
        #     "LVCF.Rdata",
        #     sep = "_"
        #   )
        # save(estres_all_h, file = filename)
        
        res_temp <- estres_all_h  %>%
          lapply(function(xxx)
            tibble(
              est = xxx[[1]]$est[1:2],
              se = xxx[[1]]$se[1:2],
              name = c("gamma1", "gamma2"),
              true_val = true_coef
            ))  %>%
          bind_rows(.id = "rep") %>% 
          mutate(
            lb = est - qnorm(0.975) * se,
            ub = est + qnorm(0.975) * se,
            CP = (lb < true_val) * (ub > true_val)
          ) %>% group_by(name) %>%
          summarise(
            bias = mean(est - true_val),
            sd = sd(est),
            se = mean(se),
            CP = mean(CP)
          ) %>% mutate(n = nn,
                       s = s_val,
                       cenrate = cenrate)
        print(res_temp)
        res_temp
      }
    }
  }


#save(res,file="/home/dsun30/asyntrans/multi6_fixed_all_beta_sum.RData")

###
# newmulti0200censor20bd4=estres_all_h
# save(newmulti0200censor20bd4,file="newmulti0200censor20bd4.RData")
# 
#  estres_all_h=newmulti1400censor40bd4
# estres_all_h  %>%
#   lapply(function(xx) tibble(est=xx[[1]]$est[1],se=xx[[1]]$se[1])) %>% 
#   bind_rows(.id="rep") %>% 
#   mutate(lb=est-qnorm(0.975)*se,
#          ub=est+qnorm(0.975)*se,
#          CP = (lb < 1)*(ub>1)) %>% 
#   summarise(bias=mean(est)-1,sd=sd(est),se=mean(se),CP=mean(CP))
# 
# estres_all_h %>%
#   lapply(function(xx) tibble(est=xx[[1]]$est[2],se=xx[[1]]$se[2])) %>% 
#   bind_rows(.id="rep") %>% 
#   mutate(lb=est-qnorm(0.975)*se,
#          ub=est+qnorm(0.975)*se,
#          CP = (lb < -0.5)*(ub>-0.5)) %>% 
#   summarise(bias=(mean(est)+0.5)/(-0.5),sd=sd(est),se=mean(se),CP=mean(CP))
