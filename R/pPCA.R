function(cnt,dis,...){
	
    require(polycor)
	n = rep(0,3)
	n[1] = ncol(cnt)
	n[2] = ncol(dis)
	
	cnt = apply(cnt,2,function(kk){(kk-mean(kk))/sd(kk)}) #standardize
	
	cormat = matrix(0,nrow = sum(n)
	                ,ncol = sum(n)
	cormat[1:n[1],1:n[1]] = cor(cnt) # standard correlation for continuous
	cormat[1:n[2]+n[1],1:n[2] + n[1]] = apply(as.matrix(1:n[2]),1,function(i,dat,n2,...){
	        apply(matrix(1:n2,1,n2),2,function(j,dat,i,...){
			    polychor(dat[,i],dat[,j],...) # use polychor for dis vars
			},dat = dat,i=i,...)
	},dat = dis,n2 = n[2],...)
    
	cormat[1:n[1],1:n[2] + n[1]] = apply(as.matrix(1:n[1]),1,function(i,dat1,dat2,n2,...){
	        apply(matrix(1:n2,1,n2),2,function(j,dat1,dat2,i,...){
			    polyserial(dat1[,i],dat2[,j],...)# polyserial for cnt and dis
			},dat1 = dat1,dat2=dat2,i=i,...)
	},dat1 = cts,dat2 = dis,n2 = n[2],...)
    cormat[1:n[2] + n[1],1:n[1]]=t(cormat[1:n[1],1:n[2] + n[1]])
    
    eig = eigen(cormat,T) # corresponding eigen problem
	
	V = eig$vector
	X = cbind(cnt,dis)
	X_new = X%*%V
	w = eig$value/sum(eig$value)
	
	return(list(eigen_vector = V,eigen_value = eig$value,weight = w,newX = X_new))

}