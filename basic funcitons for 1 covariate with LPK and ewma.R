# 核函数
kh_fun=function(tdiff,h){
  u=tdiff/h
  kh=exp(-abs(u)/2)/h
  return(kh)
}
# 归一化函数
tildefunc = function(datalist,coef.fixed,r.mat,sig.sq){
  t = datalist$tlist[[1]]
  z = datalist$ylist[[1]]
  x = datalist$xlist[[1]]
  
  t = round(as.numeric(t),0);time = round(time,0)
  id = match(t,time)
  
  # 去除含有固定效应的项，使z变为y
  y = z - x * coef.fixed[id]
  
  # y的协方差矩阵
  r.matsub = r.mat[id,id]
  y.sigma.sqrt = diag(x) %*% r.matsub %*% t(diag(x)) + diag(sig.sq, nrow(r.matsub))
  
  # 求y的协方矩阵y.sigma.sqrt的负二分之一：A
  library(MASS)
  lambda.eig <- ginv(eigen(y.sigma.sqrt)$vectors)%*%y.sigma.sqrt%*%(eigen(y.sigma.sqrt)$vectors)
  lambda.sqrt = abs(diag(sqrt(diag(lambda.eig))))
  B = (eigen(y.sigma.sqrt)$vectors)%*%lambda.sqrt%*%solve(eigen(y.sigma.sqrt)$vectors)
  A = solve(B)
  # y的独立正态变换
  y.tilde = A %*% y
  # mean(y.tilde);var(y.tilde);shapiro.test(y.tilde); acf(y.tilde)
  # mean(y);var(y);shapiro.test(y); acf(y)
  # 为了保证等式成立，还需要对x进行同等变换
  x.tilde = A %*% x
  return(list(x.tilde= x.tilde, y.tilde = y.tilde, t = t))
}
# 偏移函数
delta.func=function(t,x.tilde,y.tilde){
  ni = length(y.tilde)
  delta.m = c()
  for (mm in 1:ni) {
    # mm = 1
    tt = t[mm]
    tdiff = tt-t[1:mm]
    delta.m.term1 = sum(kh_fun(tdiff,h)*x.tilde[1:mm]^2)
    delta.m.term2 = sum(kh_fun(tdiff,h)*x.tilde[1:mm]*y.tilde[1:mm])
    delta.m = c(delta.m, delta.m.term2/delta.m.term1)
  }
  return(delta.m)
}

# 似然比统计量的常规计算
Hmfunc = function(y.tilde,y.tilde.fit){
  Hm = c()
  for (m in 1:ni) {
    # m = 1
    wjlist = c()
    for (j in 1:m) {
      wj = (1-lambda)^(t[m]-t[j])
      wjlist = c(wjlist,wj)
    }
    Hm.term1 = (y.tilde[1:m]^2)
    Hm.term2 = ((y.tilde[1:m]-y.tilde.fit[1:m])^2)
    Hm = c(Hm,sum((Hm.term1-Hm.term2)*wjlist))
  }
  return(Hm)
}
# 似然比统计量的递归计算
HmstatFunc = function(t,x.tilde,y.tilde){
  Umlist = Gmlist = deltalist =vector("list",ni)
  # 计算Um.1和Gm.1
  m = j = 1
  rj = exp(-(t[m]-t[j])/(2*h))/h
  Um = rj*x.tilde[j]*y.tilde[j]
  Gm = solve(rj*x.tilde[j]*x.tilde[j])
  Umlist[[1]] = Um; Gmlist[[1]] = Gm; deltalist[[1]] = Gm*Um
  # 递归
  for (m in 2:ni) {
    rj = exp(-(t[m]-t[m-1])/(2*h))/h
    r = exp(-(t[m]-t[m])/(2*h))/h
    Umlist[[m]] = h*rj*Umlist[[m-1]] + r*x.tilde[m]*y.tilde[m]
    Gmlist[[m]] = solve(h*rj*solve(Gmlist[[m-1]]) + r*x.tilde[m]*x.tilde[m])
    deltalist[[m]] = Gmlist[[m]]*Umlist[[m]]
  }
  # unlist(deltalist)
  Hmlist = vector("list",ni)
  m = 1
  Hmlist[[1]] = 2*y.tilde[m]*x.tilde[m]*deltalist[[m]]-(x.tilde[m]*deltalist[[m]])^2
  for (m in 2:ni) {
    wj = (1-lambda)^(t[m]-t[m-1])
    Hmlist[[m]] = (wj*Hmlist[[m-1]] + 2*y.tilde[m]*x.tilde[m]*deltalist[[m]]-(x.tilde[m]*deltalist[[m]])^2)
  }
  return(list(Umlist = Umlist, Gmlist = Gmlist, deltalist = deltalist, Hmlist = Hmlist))
}


# 统计量的期望
EHfunc = function(t,x.tilde,Gmlist){
  # EH_1
  m = 0; j = 1
  rm1 = exp(-(t[m+1]-t[j])/(2*h))/h
  qmj = x.tilde[m+1]*Gmlist[[m+1]]*x.tilde[j]
  r = exp(-(t[m+1]-t[m+1])/(2*h))/h
  EH = 2*r*x.tilde[m+1]*Gmlist[[m+1]]*x.tilde[m+1]-rm1^2*qmj^2
  
  EHlist = c()
  EHlist = c(EHlist,EH)  
  for (m in 1:(ni-1)) {
    w = (1-lambda)^(t[m+1]-t[m])
    EH.term1 = w*EH
    EH.term2 = 0
    for (j in 1:(m+1)) {
      rmj = exp(-(t[m+1]-t[j])/(2*h))/h
      qmj = x.tilde[m+1]*Gmlist[[m+1]]*x.tilde[j]
      EH.term2 = EH.term2 + rmj^2*qmj^2
    }
    EH.term3 = 2*r*x.tilde[m+1]*Gmlist[[m+1]]*x.tilde[m+1]
    EH = EH.term1 - EH.term2 + EH.term3
    EHlist = c(EHlist,EH)
  }
  return(EHlist)
}

EUUHfunc = function(t,x.tilde,Gmlist){
  m = 1; j = 1
  EUUH.term1 = 0
  
  r.m1m = exp(-(t[m+1]-t[m])/(2*h))/h
  rmj = exp(-(t[m]-t[j])/(2*h))/h
  qmj = x.tilde[m]*Gmlist[[m]]*x.tilde[j]
  EUUH.term2 = h^2*r.m1m^2*rmj^4*x.tilde[j]*x.tilde[j]*qmj^2*epsilon4
  
  r.mm = exp(-(t[m]-t[m])/(2*h))/h
  EUUH.term3 = 2*h^2*r.m1m^2*r.mm^3*x.tilde[m]*x.tilde[m]*Gmlist[[m]]*x.tilde[m]*x.tilde[m]*epsilon4
  
  EUUH.term4 = 0
  
  r.m1m1 = exp(-(t[m+1]-t[m+1])/(2*h))/h
  EUUH.term5 = r.m1m1^2*x.tilde[m+1]*x.tilde[m+1]*rmj^2*qmj^2
  
  EUUH.term6 = 2*r.m1m1^2*r.mm*x.tilde[m+1]*x.tilde[m]*Gmlist[[m]]*x.tilde[m]*x.tilde[m+1]
  
  EUUH = EUUH.term1 - EUUH.term2 + EUUH.term3 + EUUH.term4 - EUUH.term5 + EUUH.term6
  
  EUUHlist = c()
  EUUHlist = c(EUUHlist,EUUH)
  
  for (m in 2:(ni-1)) {
    r.m1m = exp(-(t[m+1]-t[m])/(2*h))/h
    w = (1-lambda)^(t[m]-t[m-1])
    EUUH.term1 = h^2*r.m1m^2*w*EUUHlist[m-1]
    
    EUUH.term2.1 = EUUH.term2.2 = EUUH.term2.3 = 0
    j = c(1:m)
    rmj = exp(-(t[m]-t[j])/(2*h))/h
    qmj = x.tilde[m]*as.numeric(Gmlist[[m]])*x.tilde[j]
    EUUH.term2.1 = sum(rmj^4*x.tilde[j]*x.tilde[j]*qmj^2*epsilon4)
    
    for (j in 1:m) {
      rmj = exp(-(t[m]-t[j])/(2*h))/h
      qmj = x.tilde[m]*Gmlist[[m]]*x.tilde[j]
      m_c = c(1:m)
      i = m_c[-j]
      rmi = exp(-(t[m]-t[i])/(2*h))/h
      qmi = x.tilde[m]*as.numeric(Gmlist[[m]])*x.tilde[i]
      EUUH.term2.2 = EUUH.term2.2 + rmj^2*x.tilde[j]*x.tilde[j]*sum(rmi^2*qmi^2)
      EUUH.term2.3 = EUUH.term2.3 + rmj^2*x.tilde[j]*qmj*sum(rmi^2*x.tilde[i]*qmi)
    }
    r.m1m = exp(-(t[m+1]-t[m])/(2*h))/h
    EUUH.term2 = h^2 * r.m1m^2 * (EUUH.term2.1 + EUUH.term2.2 + EUUH.term2.3)
    
    EUUH.term3.1 = EUUH.term3.2 = EUUH.term3.3 = 0
    j = c(1:(m-1))
    rmj = exp(-(t[m-1]-t[j])/(2*h))/h
    qmj = x.tilde[m]*as.numeric(Gmlist[[m]])*x.tilde[j]
    EUUH.term3.1 = sum(rmj^2*x.tilde[j]*qmj)*x.tilde[m]
    EUUH.term3.2 = sum(rmj^2*x.tilde[j]*x.tilde[j])
    EUUH.term3.3 = sum(rmj^2*x.tilde[j]*x.tilde[j])
    
    r.mm1 = exp(-(t[m]-t[m-1])/(2*h))/h
    r.mm = exp(-(t[m]-t[m])/(2*h))/h
    EUUH.term3.1 = h^2*r.mm1^2*r.mm*EUUH.term3.1
    EUUH.term3.2 = h^2*r.mm1^2*r.mm*x.tilde[m]*Gmlist[[m]]*x.tilde[m]*EUUH.term3.2
    EUUH.term3.3 = h^2*r.mm1^2*r.mm*x.tilde[m]*x.tilde[m]*Gmlist[[m]]*EUUH.term3.3
    EUUH.term3.4 = r.mm^3*x.tilde[m]*x.tilde[m]*Gmlist[[m]]*x.tilde[m]*x.tilde[m]*epsilon4
    r.m1m = exp(-(t[m+1]-t[m])/(2*h))/h
    EUUH.term3 = 2*h^2*r.m1m^2*(EUUH.term3.1 + EUUH.term3.2 + EUUH.term3.3 + EUUH.term3.4)
    
    r.m1m1 = exp(-(t[m+1]-t[m+1])/(2*h))/h
    EUUH.term4 = r.m1m1^2 * w * x.tilde[m+1]*x.tilde[m+1] * EHlist[m-1]
    
    EUUH.term5 = 0
    j = c(1:m)
    rmj = exp(-(t[m]-t[j])/(2*h))/h
    qmj = x.tilde[m]*as.numeric(Gmlist[[m]])*x.tilde[j]
    EUUH.term5 = sum(rmj^2*qmj^2)
    EUUH.term5 = r.m1m1^2*x.tilde[m+1]*x.tilde[m+1]*EUUH.term5
    
    EUUH.term6 = 2*r.m1m1^2*r.mm*x.tilde[m+1]*x.tilde[m]*Gmlist[[m]]*x.tilde[m]*x.tilde[m+1]
    
    EUUH = EUUH.term1 - EUUH.term2 + EUUH.term3 + EUUH.term4 - EUUH.term5 + EUUH.term6
    
    EUUHlist = c(EUUHlist,EUUH)
  }
  EUUHlist = c(0,EUUHlist)
  return(EUUHlist)
}


varALLfunc = function(t,x.tilde,Gmlist,EUUH.list){
  cov.AB.list = c(); var.B.list = c(); var.C.list = c(); cov.BC.list = c(); var.A.list = c(); Vm.list = c()
  m = 0
  covAB = 0
  cov.AB.list = c(cov.AB.list,covAB)
  
  j = 1
  rmj = exp(-(t[m+1]-t[j])/(2*h))/h
  qmj = x.tilde[m+1]*Gmlist[[m+1]]*x.tilde[j]
  var.B = (epsilon4-1)*rmj^4*qmj^4
  var.B.list = c(var.B.list,var.B)
  
  rmj = exp(-(t[m+1]-t[m+1])/(2*h))/h
  qmj = x.tilde[m+1]*Gmlist[[m+1]]*x.tilde[m+1]
  var.C = 4*(epsilon4-1)*rmj^2*qmj^2
  var.C.list = c(var.C.list,var.C)
  
  rmj = exp(-(t[m+1]-t[m+1])/(2*h))/h
  qmj = x.tilde[m+1]*Gmlist[[m+1]]*x.tilde[m+1]
  cov.BC = 2*(epsilon4-1)*rmj^3*qmj^3
  cov.BC.list = c(cov.BC.list,cov.BC)
  
  
  var.A = 0;var.A.list = c(var.A.list,var.A)
  Vm = var.A.list[m+1] + var.B.list[m+1] + var.C.list[m+1] - 2*cov.AB.list[m+1] - 2*cov.BC.list[m+1]
  Vm.list = c(Vm.list,Vm)
  
  for (m in 1:(ni-1)) {
    w = (1-lambda)^(t[m+1]-t[m])
    covAB.term1 = w*x.tilde[m+1]^2*Gmlist[[m+1]]^2*EUUH.list[m+1]
    j = c(1:(m+1))
    rmj = exp(-(t[m+1]-t[j])/(2*h))/h
    qmj = x.tilde[m+1]*as.numeric(Gmlist[[m+1]])*x.tilde[j]
    covAB.term2 = w*EHlist[m]*sum(rmj^2*qmj^2)
    covAB = covAB.term1 - covAB.term2
    cov.AB.list = c(cov.AB.list,covAB)
    
    var.B.term1 = sum(rmj^4*qmj^4*(epsilon4-1))
    var.B.term2 = 0
    m_c = c(1:(m+1))
    for (j in 1:(m+1)) {
      rmj = exp(-(t[m+1]-t[j])/(2*h))/h
      qmj = x.tilde[m+1]*as.numeric(Gmlist[[m+1]])*x.tilde[j]
      i = m_c[-j]
      rmi = exp(-(t[m+1]-t[i])/(2*h))/h
      qmi = x.tilde[m+1]*as.numeric(Gmlist[[m+1]])*x.tilde[i]
      var.B.term2 = var.B.term2 + rmj^2*qmj^2*sum(rmi^2*qmi^2)
    }
    var.B.term = var.B.term1+var.B.term2
    var.B.list = c(var.B.list, var.B.term)
    
    j = c(1:m)
    r.m1m = exp(-(t[m+1]-t[m])/(2*h))/h
    rmj = exp(-(t[m]-t[j])/(2*h))/h
    qmj = x.tilde[m+1]*as.numeric(Gmlist[[m+1]])*x.tilde[j]
    var.C.term1 = 4*h^2*r.m1m^2*sum(rmj^2*qmj^2)
    
    r.m1m1 = exp(-(t[m+1]-t[m+1])/(2*h))/h
    qmj = x.tilde[m+1]*as.numeric(Gmlist[[m+1]])*x.tilde[m+1]
    var.C.term2 = 4*r.m1m1^2*qmj^2*(epsilon4-1)
    var.C = var.C.term1+var.C.term2
    var.C.list = c(var.C.list, var.C)
    
    j = c(1:m)
    rmj = exp(-(t[m]-t[j])/(2*h))/h
    qmj = x.tilde[m+1]*as.numeric(Gmlist[[m+1]])*x.tilde[j]
    r.m1m1 = exp(-(t[m+1]-t[m+1])/(2*h))/h
    r.m1m = exp(-(t[m+1]-t[m])/(2*h))/h
    qm1m1 = x.tilde[m+1]*as.numeric(Gmlist[[m+1]])*x.tilde[m+1]
    cov.BC.term2 = 2*h^2*r.m1m1*r.m1m^2*qm1m1*sum(rmj^2*qmj^2)
    cov.BC.term1 = (epsilon4-1)*r.m1m1^3*qm1m1^3
    cov.BC = 2*(cov.BC.term1 + cov.BC.term2)
    cov.BC.list = c(cov.BC.list, cov.BC)
    
    
    w = (1-lambda)^(t[m+1]-t[m])
    var.A = w^2*Vm.list[m]
    var.A.list = c(var.A.list,var.A)
    Vm = var.A.list[m+1] + var.B.list[m+1] + var.C.list[m+1] - 2*cov.AB.list[m+1] - 2*cov.BC.list[m+1]
    Vm.list = c(Vm.list,Vm)
  }
  return(list(cov.AB.list = cov.AB.list, var.B.list = var.B.list, var.C.list = var.C.list, cov.BC.list = cov.BC.list,Vm.list = Vm.list)) 
}



