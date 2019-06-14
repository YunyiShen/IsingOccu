pPCA = function(cnt, dis, ...) {
  require(polycor)
  cnt = as.matrix(cnt)
  dis = as.matrix(dis)
  n = rep(0, 2)
  n[1] = ncol(cnt)
  n[2] = ncol(dis)
  
  cnt = apply(cnt, 2, function(kk) {
    (kk - mean(kk)) / sd(kk)
  }) #standardize
  
  cormat = matrix(0, nrow = sum(n)
                  , ncol = sum(n))
  cormat[1:n[1], 1:n[1]] = cor(cnt) # standard correlation for continuous
  cormat[1:n[2] + n[1], 1:n[2] + n[1]] = apply(as.matrix(1:n[2]), 1, function(i, dat, n2, ...) {
        apply(matrix(1:n2, 1, n2), 2, function(j, dat, i, ...) {
              if (i==j) return(1)
              polychor(dat[, i], dat[, j], ...) # use polychor for dis vars
            }, dat = dat, i = i, ...)
          }, dat = dis, n2 = n[2], ...)
  
  
  
  cormat[1:n[1], 1:n[2] + n[1]] = 
    apply(as.matrix(1:n[1]), 1, function(i, dat1, dat2, n2, ...) {
        apply(matrix(1:n2, 1, n2), 2, function(j, dat1, dat2, i, ...) {
                  polyserial(dat1[, i], dat2[, j], ...)# polyserial for cnt and dis
               }, dat1 = dat1, dat2 = dat2, i = i, ...)
           }, dat1 = cnt, dat2 = dis, n2 = n[2], ...)
  cormat[1:n[2] + n[1], 1:n[1]] = t(cormat[1:n[1], 1:n[2] + n[1]])
  
  diag(cormat) = 1
  
  eig = eigen(cormat, T) # corresponding eigen problem
  
  V = eig$vector
  X = cbind(cnt, dis)
  X_new = X %*% V
  w = eig$value / sum(eig$value)
  unexplained = lapply(1:ncol(X),function(i,ww){1-sum(ww[1:i])},ww = w)
  
  return(list(
    original_X = X,
    correlation_matrix = cormat,
    eigen_vectors = V,
    eigen_values = eig$value,
    weight = w,
    unexplained_sd = unlist( unexplained),
    newX = X_new
  ))
  
}