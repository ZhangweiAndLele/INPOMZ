
ICdata.model1 = function(rate,N,tc.end){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = xlist = vilist = ylist = vector("list",N)
  nilist = c()
  for (i in 1:N) {
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    vilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    xlist[[i]] =sin(0.01*pi*tlist[[i]])+rnorm(ni,0,0.1)
    beta = (1+tlist[[i]])^0.3
    ylist[[i]] = beta*xlist[[i]] + vilist[[i]]*xlist[[i]] + rnorm(ni,0,1)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist,xlist = xlist,ylist = ylist,nilist = nilist,dayslist = dayslist))
}


ICdata.model2 = function(rate,N,tc.end){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = x1list = v1ilist = x2list = v2ilist = ylist = vector("list",N)
  nilist = c()
  
  for (i in 1:N) {
    # i = 1
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    v1ilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    v2ilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    x1list[[i]] = sin(0.01*tlist[[i]])+rnorm(ni,0, 0.1)
    x2list[[i]] = cos(0.05*tlist[[i]])+rnorm(ni,0, 0.1)
    beta1 = (1+0.02*tlist[[i]])
    beta2 = (1+0.05*tlist[[i]])
    ylist[[i]] = (beta1+ v1ilist[[i]])*x1list[[i]] + (beta2+ v2ilist[[i]])*x2list[[i]] + rnorm(ni,0,0.1)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist,x1list = x1list,x2list = x2list,ylist = ylist,nilist = nilist,dayslist = dayslist))
}


OCdata.case1 = function(rate,N,tc.end,delta){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = xlist = vilist = ylist = vector("list",N)
  nilist = c()
  for (i in 1:N) {
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    vilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    xlist[[i]] =sin(0.01*pi*tlist[[i]])+rnorm(ni,0,0.1)
    beta = (1+tlist[[i]])^0.3*(1+delta)
    ylist[[i]] = beta*xlist[[i]] + vilist[[i]]*xlist[[i]] + rnorm(ni,0,1)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist,xlist = xlist,ylist = ylist,nilist = nilist,dayslist = dayslist))
}

OCdata.case2 = function(rate,N,tc.end,delta){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = xlist = vilist = ylist = vector("list",N)
  nilist = c()
  for (i in 1:N) {
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    vilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    xlist[[i]] =sin(0.01*pi*tlist[[i]])+rnorm(ni,0,0.1)
    drift = ((tc-tc[1])/50)^2
    beta = (1+tlist[[i]])^0.3*(1+drift*delta)
    ylist[[i]] = beta*xlist[[i]] + vilist[[i]]*xlist[[i]] + rnorm(ni,0,1)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist,xlist = xlist,ylist = ylist,nilist = nilist,dayslist = dayslist))
}

OCdata.case3 = function(rate,N,tc.end,delta){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = xlist = vilist = ylist = vector("list",N)
  nilist = c()
  for (i in 1:N) {
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    vilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    xlist[[i]] =sin(0.01*pi*tlist[[i]])+rnorm(ni,0,0.1)
    oscillation = sin((tc-tc[1])*0.1*pi)
    beta = (1+tlist[[i]])^0.3*(1+oscillation*delta)
    ylist[[i]] = beta*xlist[[i]] + vilist[[i]]*xlist[[i]] + rnorm(ni,0,1)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist,xlist = xlist,ylist = ylist,nilist = nilist,dayslist = dayslist))
}

OCdata.case4 = function(rate,N,tc.end,delta){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = x1list = v1ilist = x2list = v2ilist = ylist = vector("list",N)
  nilist = c()
  
  for (i in 1:N) {
    # i = 1
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    v1ilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    v2ilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    x1list[[i]] = sin(0.01*tlist[[i]])+rnorm(ni,0, 0.1)
    x2list[[i]] = cos(0.05*tlist[[i]])+rnorm(ni,0, 0.1)
    beta1 = (1+0.02*tlist[[i]])*(1+delta)
    beta2 = (1+0.05*tlist[[i]])*(1+delta)
    
    ylist[[i]] = (beta1+ v1ilist[[i]])*x1list[[i]] + (beta2+ v2ilist[[i]])*x2list[[i]] + rnorm(ni,0,0.5)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist,x1list = x1list,x2list = x2list,ylist = ylist,nilist = nilist,dayslist = dayslist))
}

OCdata.case5 = function(rate,N,tc.end,delta){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = x1list = v1ilist = x2list = v2ilist = ylist = vector("list",N)
  nilist = c()
  
  for (i in 1:N) {
    # i = 1
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    v1ilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    v2ilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    
    drift = ((tc-tc[1])/50)^2
    x1list[[i]] = sin(0.01*tlist[[i]])+rnorm(ni,0, 0.1)
    x2list[[i]] = cos(0.05*tlist[[i]])+rnorm(ni,0, 0.1)
    beta1 = (1+0.02*tlist[[i]])*(1+drift*delta)
    beta2 = (1+0.05*tlist[[i]])*(1+drift*delta)

    ylist[[i]] = (beta1+ v1ilist[[i]])*x1list[[i]] + (beta2+ v2ilist[[i]])*x2list[[i]] + rnorm(ni,0,0.5)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist,x1list = x1list,x2list = x2list,ylist = ylist,nilist = nilist,dayslist = dayslist))
}

OCdata.case6 = function(rate,N,tc.end,delta){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = x1list = v1ilist = x2list = v2ilist = ylist = vector("list",N)
  nilist = c()
  
  for (i in 1:N) {
    # i = 1
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    v1ilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    v2ilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    oscillation = sin((tc-tc[1])*0.01*pi)
    x1list[[i]] = sin(0.01*tlist[[i]])+rnorm(ni,0, 0.1)
    x2list[[i]] = cos(0.05*tlist[[i]])+rnorm(ni,0, 0.1)
    beta1 = (1+0.02*tlist[[i]])*(1+oscillation*delta)
    beta2 = (1+0.05*tlist[[i]])*(1+oscillation*delta)
    
    ylist[[i]] = (beta1+ v1ilist[[i]])*x1list[[i]] + (beta2+ v2ilist[[i]])*x2list[[i]] + rnorm(ni,0,0.5)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist,x1list = x1list,x2list = x2list,ylist = ylist,nilist = nilist,dayslist = dayslist))
}


ICdata.model3 = function(rate,N,tc.end){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = x1list = vilist = x2list = ylist = vector("list",N)
  nilist = c()
  
  for (i in 1:N) {
    # i = 1
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    vilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    x1list[[i]] = sin(0.001*tlist[[i]]*20*pi)+rnorm(ni,0, 0.01)
    x2list[[i]] = cos(0.001*tlist[[i]]*10*pi)+rnorm(ni,0, 0.01)
    beta1 = 1
    beta2 = 1.5
    ylist[[i]] = beta1*x1list[[i]] + beta2*x2list[[i]] + vilist[[i]] + rnorm(ni,0,0.01)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist, x1list = x1list, x2list = x2list, ylist = ylist, nilist = nilist, dayslist = dayslist))
}


OCdata.case7 = function(rate,N,tc.end,delta){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = x1list = vilist = x2list = ylist = vector("list",N)
  nilist = c()
  
  for (i in 1:N) {
    # i = 1
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    vilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    x1list[[i]] = sin(0.001*tlist[[i]]*20*pi)+rnorm(ni,0, 0.01)
    x2list[[i]] = cos(0.001*tlist[[i]]*10*pi)+rnorm(ni,0, 0.01)
    beta1 = 1*(1+delta)
    beta2 = 1.5*(1+delta)
    ylist[[i]] = beta1*x1list[[i]] + beta2*x2list[[i]] + vilist[[i]] + rnorm(ni,0,0.01)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist, x1list = x1list, x2list = x2list, ylist = ylist, nilist = nilist, dayslist = dayslist))
}

OCdata.case8 = function(rate,N,tc.end,delta){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = x1list = vilist = x2list = ylist = vector("list",N)
  nilist = c()
  
  for (i in 1:N) {
    # i = 1
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    vilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    x1list[[i]] = sin(0.001*tlist[[i]]*20*pi)+rnorm(ni,0, 0.01)
    x2list[[i]] = cos(0.001*tlist[[i]]*10*pi)+rnorm(ni,0, 0.01)
    drift = ((tc-tc[1])/50)^2
    beta1 = 1*(1+drift*delta)
    beta2 = 1.5*(1+drift*delta)
    ylist[[i]] = beta1*x1list[[i]] + beta2*x2list[[i]] + vilist[[i]] + rnorm(ni,0,0.01)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist, x1list = x1list, x2list = x2list, ylist = ylist, nilist = nilist, dayslist = dayslist))
}


OCdata.case9 = function(rate,N,tc.end,delta){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = x1list = vilist = x2list = ylist = vector("list",N)
  nilist = c()
  
  for (i in 1:N) {
    # i = 1
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    vilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    x1list[[i]] = sin(0.001*tlist[[i]]*20*pi)+rnorm(ni,0, 0.01)
    x2list[[i]] = cos(0.001*tlist[[i]]*10*pi)+rnorm(ni,0, 0.01)
    oscillation = sin((tc-tc[1])*0.01*pi)
    beta1 = 1*(1+oscillation*delta)
    beta2 = 1.5*(1+oscillation*delta)
    ylist[[i]] = beta1*x1list[[i]] + beta2*x2list[[i]] + vilist[[i]] + rnorm(ni,0,0.01)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist, x1list = x1list, x2list = x2list, ylist = ylist, nilist = nilist, dayslist = dayslist))
}

ICdata.model4 = function(rate,N,tc.end){
  gaussprocess <- function(from = 0,
                           to = tc.end,
                           K = function(s, t) {min(s, t)},
                           start = NULL,
                           m = ni) {
    t <- seq(from = from, to = to, length.out = m)
    Sigma <- sapply(t, function(s1) {
      sapply(t, function(s2) {
        K(s1, s2)
      })
    })
    path <- MASS::mvrnorm(mu = rep(0, times = m), Sigma = Sigma)
    if (!is.null(start)) {
      path <- path - path[1] + start  # Must always start at "start"
    }
    return(data.frame("t" = t, "xt" = path))
  }
  
  gaussian_data_generate = function(times = N, n = ni) {
    gaussian_data <- gaussprocess(K = function(s, t) {exp(-(s - t) ^ 2)})
    # for(i in 1:(times-1)){
    #     temp <- gaussprocess(K = function(s, t) {exp(-16 * (s - t) ^ 2)})[,2]
    #     gaussian_data <- cbind(gaussian_data, temp)
    # }
    return(gaussian_data[,-1])
  }
  
  dayslist = tlist = xlist = vilist = ylist = vector("list",N)
  nilist = c()
  for (i in 1:N) {
    # set.seed(i)
    tc = c(); tc.ini = 0; j = 1
    while(tc.ini<=tc.end){
      d = round(rexp(1,rate),7)
      if (d == 0) {
        d = round(rexp(1,rate),7)
      }else{
        tc.ini = tc.ini + d
        tc = c(tc,tc.ini)
        j = j + 1
      }
    }
    tc = head(tc,-1)
    ni = length(tc)
    vilist[[i]] = gaussian_data_generate(times = 1, n = ni)
    dayslist[[i]] = rep(i,ni)
    tlist[[i]] = tc
    xlist[[i]] =sin(0.01*pi*tlist[[i]])+rnorm(ni,0,1)
    beta = (1+tlist[[i]])^0.3
    ylist[[i]] = beta*xlist[[i]] + vilist[[i]]*xlist[[i]] +rt(ni,2.5)
    nilist = c(nilist,ni)
  }
  return(list(tlist = tlist,xlist = xlist,ylist = ylist,nilist = nilist,dayslist = dayslist))
}
