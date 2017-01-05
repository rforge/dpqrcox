
dpqr_cox <- function(x) {
    .Call('dpqrcox_dpqr_cox',x, PACKAGE = 'dpqrcox')
}

dcox <- function(x, probvec, lambdvec, log = FALSE){
  dpqrflag <- 1
  dpqrname <- "dcox"
  outy <- NA
  if(check_cox(dpqrname, x, probvec, lambdvec, log = log)){
    outy <- base_fun(dpqrflag, x, probvec, lambdvec, log = log)
  }
  if(log) outy <- log(outy);
  outy
}

pcox <- function(x, probvec, lambdvec, lower.tail = TRUE, log.p = FALSE){
  dpqrflag <- 2
  dpqrname <- "pcox"
  outy <- NA
  if(check_cox(dpqrname, x, probvec, lambdvec, lower.tail = lower.tail, log.p = log.p)){
    outy <- base_fun(dpqrflag, x, probvec, lambdvec, lower.tail = lower.tail, log.p = log.p)
  }
  if(!lower.tail) outy <- rep(1,length(outy))-outy;
  if(log.p) outy <- log(outy);
  outy
}

qcox <- function(p, probvec, lambdvec, lower.tail = TRUE, log.p = FALSE){
  dpqrflag <- 3
  dpqrname <- "qcox"
  outy <- NA
  if(check_cox(dpqrname, p, probvec, lambdvec, lower.tail = lower.tail, log.p = log.p)){
#    outy <- base_fun(dpqrflag, p, probvec, lambdvec, lower.tail = lower.tail, log.p = log.p)
    if(!lower.tail) p <- rep(1,length(p))-p;
    if(log.p) p <- exp(p);
    rsample <- rcox(10000, probvec, lambdvec) #生成抽样数据
    n <- length(p)
    outy <- rep(0,n)
    tmp <- quantile(rsample, probs = p) #计算抽样数据的分位数
    for(i in 1:n) outy[i] <- mean(tmp[i]) #消除R给出的分位数的格式
  }
  outy
}

rcox <- function(n, probvec, lambdvec){
  dpqrname <- "rcox"
  rcox_outy <- NA
  if(check_cox(dpqrname, n, probvec, lambdvec)){
    if(length(n)>1) n <- length(n)
    m <- length(lambdvec)
    rlambdvec <- rexp(n*m, rate = 1)/lambdvec
    dim(rlambdvec) <- c(m,n) 
    rprobvec <- runif(n, min = 0, max = 1)
    cdv <- rep(0,m)
    cdv[[1]] <- probvec[[1]]
    if(m>=3) for(i in 2:(m-1)) cdv[[i]] <- cdv[[i-1]] + probvec[[i]]
    cdv[[m]] <- Inf 
    for(i in 1:n){
      j <- 1
      while(rprobvec[i] > cdv[[j]]) j <- j+1
      rprobvec[i] <- j
    }
    rcox_outy <- 0*rep(0,n)
    for(j in 1:n)rcox_outy[[j]] <- sum((rlambdvec[,j])[1:rprobvec[[j]]])
    rcox_outy
  }
}

#hc 检查各个程序的参数是否合适
check_cox <- function(dpqrname, xdata, probvec, lambdvec, lower.tail = TRUE, log.p = FALSE, log = FALSE){
  check_result <- TRUE
  
 #hc 一般条件：xdata, probvec, lambdvec 全是数值量
  if(!(is.numeric(xdata)|is.numeric(probvec)|is.numeric(probvec))){
    check_result <- FALSE
    warning(paste(c(dpqrname,": ", "the first three arguments must be numeric"), collapse=""), call. = FALSE)
  }

 #hc 一般条件：lower.tail, log.p, log 全是长度为1的逻辑量
  if(length(lower.tail)!=1|length(log.p)!=1|length(log)!=1){
    check_result <- FALSE
    warning(paste(c(dpqrname,": ", "length of lower.tail, log.p, and log must be 1"), collapse=""), call. = FALSE)
  }
  if(!(is.logical(lower.tail)&is.logical(log.p)&is.logical(log))){
    check_result <- FALSE
    warning(paste(c(dpqrname,": ", "lower.tail, log.p, and log must be logical"), collapse=""), call. = FALSE)
  }

 #hc 一般条件：probvec与lambdvec长度相同，probvec介于0,1间，且probvec为正
 #hc probvec各元素和为1，注意检查除末元素之外所有的和是否小于1
  if(length(probvec)!=length(lambdvec)){
    check_result <- FALSE
    warning(paste(c(dpqrname,": ", "lengths of probvec and lambdvec must be equal"), collapse=""), call. = FALSE)
  }
  if(any(probvec<0|probvec>1)){
    check_result <- FALSE
    warning(paste(c(dpqrname,": ", "all elements in probvec must be in [0, 1]"), collapse=""), call. = FALSE)
  }
  if(sum(probvec)-probvec[[length(probvec)]] > 1.0){
    check_result <- FALSE
    warning(paste(c(dpqrname,": ", "summation of probvec is greater than 1"), collapse=""), call. = FALSE)
  }
  if(any(lambdvec<=0)){
    check_result <- FALSE
    warning(paste(c(dpqrname,": ", "all elements in lambdvec must be positive"), collapse=""), call. = FALSE)
  }

 #hc 条件：dcox,pcox中，xdata是非负的
  if(dpqrname=="dcox" |dpqrname== "pcox"){
    if(any(xdata<0)){
      check_result <- FALSE
      warning(paste(c(dpqrname,": ", "some element(s) in vector of position are negative"), collapse=""), call. = FALSE)
    }
  }

 #hc 条件：qcox中，xdata介于0,1间,或取对数后是非正的
  if(dpqrname=="qcox"){
    if(log.p){
      if(any(xdata>0)){
        check_result <- FALSE
        warning(paste(c(dpqrname,": ", "all element in p must be negative"), collapse=""), call. = FALSE)
      }
    } else {
      if(any(xdata<0|xdata>1)){
        check_result <- FALSE
        warning(paste(c(dpqrname,": ", "all elements in p must be in [0, 1]"), collapse=""), call. = FALSE)
      }
    }
  }

 #hc 条件：rcox中，xdata是长度大于1的向量或正整数
  if(dpqrname=="rcox"){
    if(!(length(xdata)>1|(is.double(xdata)&round(xdata)>0))){
      check_result <- FALSE
      warning(paste(c(dpqrname,": ", "n must be a vector or positive integer"), collapse=""), call. = FALSE)
    }
  }

  check_result
}

base_fun <- function(dpqrflag, xdata, probvec, lambdvec, lower.tail = TRUE, log.p = FALSE, log = FALSE){
  logv2double <- rep(0,3)
  if(lower.tail) logv2double[[1]] <- 1
  if(log.p) logv2double[[2]] <- 1
  if(log) logv2double[[3]] <- 1
  y <- Im(sort(1i*1:length(lambdvec)+lambdvec))
  inx <- c(dpqrflag, length(xdata), length(probvec), logv2double, xdata, probvec, lambdvec, y)
  outy <- result_cox(inx)
  outy 
}

result_cox <- function(inx){
  dpqr_cox(inx) 
}


