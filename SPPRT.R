#  optimal sequentially planned probability ratio test for Bernoulli observations


prRatio <- function(j, y, th0, th1){
  dbinom(y, j, th1) / dbinom(y, j, th0)
}


cost <- function(n){
  n
}


back_first <- function(z, j, l0, l1, th0, th1){
  s = 0
  for (y in (0:j))
    s = s + min(l0, l1 * z * prRatio(j, y, th0, th1)) * dbinom(y, j, th0)
  s
}


back_next <- function(z, j, vals, rho, l0, l1, th0, th1){

  incorp <- function(z){
    lz = log(z)
    if((lz > A) & (lz < B))
      approx(vals, rho, xout=lz)$y
    else
      min(l0, l1 * (z))
  }

  A = head(vals, 1)
  B = tail(vals, 1)
  s = 0
  for (y in (0:j))
  s = s + incorp(z * prRatio(j, y, th0, th1)) * dbinom(y, j, th0)
  s
}


design_test <- function(l0, l1, th0, th1, H, gsizes, gamma, h, step_fn=NULL){

  tmp <- function(x, z, valsold, rhoold)
    cost(x)*((1 - gamma) + gamma * z) + back_next(z, x, valsold, rhoold, l0, l1, th0, th1)

  tmp1 <- function(z, valsold, rhoold)
    min(Vectorize(tmp,vectorize.args="x")(gsizes, z, valsold, rhoold))

  tmp2 <- function(z, valsold, rhoold)
    which.min(Vectorize(tmp,vectorize.args="x")(gsizes, z, valsold, rhoold))

  step_eff <- function(z)
    min(l0, l1 * (z)) - min(Vectorize(tmp, vectorize.args="x")(gsizes, z, valsold, rhoold))

  b = uniroot(
    function(z)
      min(l0, l1 * (z)) - min(
        sapply(gsizes, function(x) cost(x) * ((1 - gamma) + gamma * z) + back_first(z, x, l0, l1, th0, th1))
      ),
    c(l0 / l1,l0 / l1 * 2),
    extendInt="downX")

  a = uniroot(
    function(z)
      min(l0, l1 * (z)) - min(
        sapply(gsizes, function(x) cost(x) * ((1 - gamma) + gamma * z) + back_first(z, x, l0, l1, th0, th1))
      ),
    c(0,l0/l1),
    extendInt="upX")

  A = log(a$root)
  B = log(b$root)
  nint = ceiling((B - A) / h)
  vals = seq(A, B, length.out=nint + 1)

  rho = sapply(
    exp(vals),
    function(z)
      min(sapply(gsizes, function(x) cost(x) * ((1 - gamma) + gamma * z) + back_first(z, x, l0, l1, th0, th1)))
  )

  test = list()
  test[[H]] = list(H=H, l0=l0, l1=l1, th0=th0, th1=th1, gsizes=gsizes, gamma=gamma, h=h)
  valsold = vals
  rhoold = rho
  test[[H - 1]] = list(H=H - 1, vals=vals, rho=rho)

  if(!is.null(step_fn))
    step_fn(H - 1)

  if(H > 2)
    for(i in seq(2, H - 1)){
      B = log(uniroot(step_eff, c(l0 / l1, l0 / l1 * 2), extendInt="downX")$root)
      A = log(uniroot(step_eff, c(0, l0 / l1), extendInt="upX")$root)
      nint = ceiling((B - A) / h)
      vals = seq(A, B, length.out=nint + 1)
      rho = Vectorize(tmp1, vectorize.args="z")(exp(vals), valsold, rhoold)

      test[[H - i]] = list(H=H - i, vals=vals, rho=rho)
      valsold = vals
      rhoold = rho

      if(!is.null(step_fn))
        step_fn(H - i)
    }

  test
}


operating_characteristic <- function(test, th, step_fn=NULL){

  back_first_operating <- function(z, j, th){
    s = 0
    for (y in (0:j))
    s = s + dbinom(y, j, th) * ifelse(l0 > l1 * z * prRatio(j, y, th0, th1), 1, 0)

    s
  }

  back_next_operating <- function(z, j, vals, rho, th){

    incorp <- function(z){
      lz = log(z)
      if((lz > A) & (lz < B))
        approx(vals, rho, xout=lz)$y
      else
        ifelse(l0 > l1 * (z), 1, 0)
    }

    A = head(vals, 1)
    B = tail(vals, 1)

    s = 0
    for (y in (0:j))
      s = s + incorp(z * prRatio(j, y, th0, th1)) * dbinom(y, j, th)

    s
  }

  control_first_fn <- function(z){
    tmp <- function(x)
      cost(x)*((1-gamma)+z*((gamma)))+back_first(z,x,l0,l1,th0,th1)

    gsizes[which.min(Vectorize(tmp)(gsizes))]
  }

  control_next_fn <- function(z, vals, rho){
    tmp <- function(x)
      cost(x)*((1-gamma)+z*((gamma)))+back_next(z,x,vals,rho,l0,l1,th0,th1)

    gsizes[which.min(Vectorize(tmp)(gsizes))]
  }


  H = length(test)
  gsizes = test[[H]]$gsizes
  gamma = test[[H]]$gamma
  l0 = test[[H]]$l0
  l1 = test[[H]]$l1
  th0 = test[[H]]$th0
  th1 = test[[H]]$th1


  rho = test[[H - 1]]$rho
  vals = test[[H - 1]]$vals
  A = head(vals, 1)
  B = tail(vals, 1)

  alpha = sapply(
    exp(vals),
    function(z) back_first_operating(z, control_first_fn(z), th)
  )

  valsold = vals
  rhoold = rho
  alphaold = alpha

  if(!is.null(step_fn))
    step_fn(H - 1)

  if(H > 2)
    for(n in seq(H - 2, 1)){
      rho = test[[n]]$rho
      vals = test[[n]]$vals
      A = head(vals, 1)
      B = tail(vals, 1)

      alpha = sapply(
        exp(vals),
        function(z)
          back_next_operating(z, control_next_fn(z, valsold, rhoold), valsold, alphaold, th)
      )

      valsold = vals
      rhoold = rho
      alphaold = alpha

      if(!is.null(step_fn))
        step_fn(n)
    }

  back_next_operating(1, control_next_fn(1, valsold, rhoold), valsold, alphaold, th)
}


sampling_report <- function(test, th, step_fn=NULL, accounting_fn=cost){

  back_first_loss <- function(z, j, th){
    0
  }

  back_next_loss <- function(z, j, vals, rho, th){
    incorp <- function(z){
      lz = log(z)
      if(lz > A & lz < B)approx(vals, rho, xout=lz)$y
      else 0
    }
    A = head(vals, 1)
    B = tail(vals, 1)
    s = 0
    for (y in (0:j))
      s = s + incorp(z * prRatio(j, y, th0, th1)) * dbinom(y, j, th)

    s
  }

  control_first_fn <- function(z){
    tmp <- function(x)
      cost(x)*((1-gamma)+z*((gamma)))+back_first(z, x, l0, l1, th0, th1)

    gsizes[which.min(Vectorize(tmp)(gsizes))]
  }

  control_next_fn <- function(z, vals, rho){
    tmp <- function(x)
      cost(x)*((1-gamma)+z*((gamma)))+back_next(z, x, vals, rho, l0, l1, th0, th1)

    gsizes[which.min(Vectorize(tmp)(gsizes))]
  }

  H = length(test)

  gsizes = test[[H]]$gsizes
  gamma = test[[H]]$gamma
  l0 = test[[H]]$l0
  l1 = test[[H]]$l1
  th0 = test[[H]]$th0
  th1 = test[[H]]$th1

  rho = test[[H - 1]]$rho
  vals = test[[H - 1]]$vals
  A = head(vals, 1)
  B = tail(vals, 1)

  loss = sapply(
    exp(vals),
    function(z){
        j = control_first_fn(z)
        accounting_fn(j) + back_first_loss(z, j, th)
    }
  )

  valsold = vals
  rhoold = rho
  lossold = loss

  if(!is.null(step_fn))
    step_fn(H - 1)

  if(H > 2)
    for(n in seq(H - 2, 1)){
      rho = test[[n]]$rh
      vals = test[[n]]$vals
      A = head(vals, 1)
      B = tail(vals, 1)

      tmp <- function(z){
        j = control_next_fn(z, valsold, rhoold)
        accounting_fn(j) + back_next_loss(z, j, valsold, lossold, th)
      }

      loss = Vectorize(tmp)(exp(vals))
      valsold = vals
      rhoold = rho
      lossold = loss

      if(!is.null(step_fn))
        step_fn(n)
    }

  j = control_next_fn(1, valsold, rhoold)
  s = accounting_fn(j) + back_next_loss(1, j, valsold, lossold, th)

  return(s)
}
