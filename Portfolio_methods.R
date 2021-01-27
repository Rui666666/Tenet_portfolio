markowitz = function(mu,cv,Er) {
  n = length(mu)
  wuns = matrix(1,n,1)
  A = t(wuns) %*% solve(cv) %*% mu
  B = t(mu) %*% solve(cv) %*% mu
  C = t(wuns) %*% solve(cv) %*% wuns
  D = B*C - A^2
  lam = (C*Er-A)/D
  gam = (B-A*Er)/D
  wts = lam[1]*(solve(cv) %*% mu) + gam[1]*(solve(cv) %*% wuns)
  g = (B[1]*(solve(cv) %*% wuns) - A[1]*(solve(cv) %*% mu))/D[1]
  h = (C[1]*(solve(cv) %*% mu) - A[1]*(solve(cv) %*% wuns))/D[1]
  wts = g + h*Er
}

markowitz_nss = function(n,mu,cv,Er,nss) {
  Bmat = matrix(0,n,n)    #No Short sales matrix
  diag(Bmat) = 1
  Amat <- matrix(c(mu,rep(1,n)),n,2) 
  if (nss==1) { Amat = matrix(c(Amat,Bmat),n,2+n) }
  bvec = matrix(c(Er,1),2,1)
  if (nss==1) { bvec = t(c(bvec,matrix(0,n,1))) }
  Dmat <-cv ## matrix in the quadratic function 2*cv 
  dvec <- matrix(0,n,1)  #vector in the quadratic function
  sol=solve.QP(Dmat,dvec,Amat,bvec,meq=2)
}

ltec = function(n,mu,l_mu,Er) {
  f.obj <- l_mu  # Set coefficients of the objective function
  # Set matrix corresponding to coefficients of constraints by rows
  f.con <- matrix(c(rep(1,n),colMeans(rdata[,c(2:(n+1))])), nrow = 2, byrow = TRUE)
  f.dir <- c("=",">=") # Set unequality signs
  f.rhs <- c(1,Er) # Set right hand side coefficients
  sol=lp("min", f.obj, f.con, f.dir, f.rhs)
}  


qtec = function(n,mu,l_mu,l_cv,Er,nss,gamma) {
  nss=0
  Bmat = matrix(0,n,n)    #No Short sales matrix
  diag(Bmat) = 1
  Amat <- matrix(c(mu,rep(1,n)),n,2) 
  if (nss==1) { Amat = matrix(c(Amat,Bmat),n,2+n) }
  bvec = matrix(c(Er,1),2,1)
  if (nss==1) { bvec = t(c(bvec,matrix(0,n,1))) }
  if (is.element("TRUE", c(l_cv < 0))){
    pd_D_mat <- nearPD(l_cv)
    l_cv=as.matrix(pd_D_mat$mat)
  }
  Dmat <- gamma*2*l_cv  
  dvec <- (1-gamma)*as.matrix(-l_mu) 
  sol=solve.QP(Dmat,dvec,Amat,bvec,meq=2)
  sol = abs(sol$solution)/sum(abs(sol$solution))
} 



