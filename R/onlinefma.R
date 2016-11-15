onlinefma <-
function(y,k,r,x.z=NULL,x.w=NULL,it=1,eps=0.0001,seed=4,scaling=FALSE,init=NULL,no_iter=100,batch_size=1000)
{

ptm <- proc.time()

set.seed(seed)

if (scaling) y<-scale(y)

numobs<-nrow(y)
p<-ncol(y)
ybar<- apply(y, 2, mean)
y<-scale(y, ybar, scale=FALSE) 
lik<--100000000000

if(!is.null(x.w)) x.w <- x.w[,1:99]
if(!is.null(x.z)) x.z <- x.z[,1:99]

x.z<-cbind(rep(1,numobs),x.z)
q.z<-ncol(x.z) 
x.w<-cbind(rep(1,numobs),x.w)
q.w<-ncol(x.w) 


## initialization of H e psi

if (is.null(init)) { stima=try(factanal(y,r,rotation="none"),silent=TRUE) 
                          if (is.character(stima)) {psi=0.1*diag(p)
                                           H=matrix(runif(p*r),p,r)}

                          if (!is.character(stima)) {psi=diag(stima$uniq)
                                                     H=stima$load}

                            } else {psi<-init$psi
                                    H=init$H}
                                               

if (!is.character(stima)) {
	z=matrix(factanal(y,r,scores="Bartlett")$scores,numobs)
	if (k>1) memb=kmeans(z,k)$cl else memb=rep(1,numobs)
} 
else {
	if (k>1) memb=kmeans(y,k)$cl else memb=rep(1,numobs)
}

for  (i in 1:k) 
	if ((table(memb)[i])<2) 
		memb[sample(1:numobs,2,replace=FALSE)]=i

if (is.null(init$phi)) {
	phi<-matrix(0,k,q.w) 
} 
else phi<-init$phi

if (is.null(init$w)) {  
	w<-table(memb)/sum(table(memb))
    w<-t(t(w))
} 
else 
	w<-init$w
                                

w<-matrix(w,nrow=k,ncol=numobs)
                      
sigma<-array(0,c(k,r,r))
Beta<-array(0,c(k,r,q.z))

if (q.z>1) if (is.null(init$Beta)) { for (i in 1:k) Beta[i,,]=runif(r*q.z,-1,1)} else Beta<-init$Beta
if (q.z==1) if (is.null(init$Beta)) { for (i in 1:k) Beta[i,,1]=colMeans(matrix(z[memb==i,],ncol=r))} else Beta<-init$Beta
if (is.null(init$sigma)) { for (i in 1:k) sigma[i,,]<-var(z[memb==i,])} else sigma<-init$sigma

rate <- -1
for (i in 1:no_iter) {

	gamma <- (i+1)^rate
	# print (batch_size)
	batch_indices <- c(sample(1:numobs, batch_size))
	y_batch <- y[batch_indices,]
	w_batch <- w[,batch_indices]
	# print (y_batch)
	if(!is.null(x.w)) x.w_batch <- x.w[batch_indices,]
	else x.w_batch <- NULL
	if(!is.null(x.z)) x.z_batch <- x.z[batch_indices,]
	else x.z_batch <- NULL
 	# print (dim(w))
 	# print (dim(w_batch))
	out<-fma.em.alg(y_batch,batch_size,r,k,p,x.z_batch,q.z,x.w_batch,q.w,phi,it,H,w_batch,Beta,sigma,eps,psi,lik)

	likelihood<-out$likelihood
	ph.y<-out$ph.y
	py.h<-out$py.h
	sigma<-(1-gamma)*sigma+gamma*(out$sigma)
	H<-(1-gamma)*H+gamma*(out$H)
	w_batch<-out$w
	Beta<-(1-gamma)*Beta+gamma*(out$Beta)
	psi<-(1-gamma)*psi+gamma*(out$psi)
	psi<-diag(diag(psi))
	phi<-(1-gamma)*phi+gamma*(out$phi)
	likelihood<-matrix(likelihood[!likelihood==0])

	h<-(k-1)*(q.w+(r*q.z)+r*(r+1)/2)+p*(r+1)   
	pen<-h*log(numobs)
	lik<-likelihood[length(likelihood)]
	bic<--2*lik+pen
	aic<--2*lik+2*h
	print(lik)

	l <- 1
	for (j in  batch_indices){
		w[,j] <- (1-gamma)*w[,j]+gamma*w_batch[,l]
		l <- l+1
	}
}

print ('finished')
if (k>1)  {  index<-(apply(ph.y, 2, order))[k,]} else index<-rep(k,numobs) 
z<-t(solve(t(H)%*%solve(psi)%*%H)%*%t(H)%*%solve(psi)%*%t(y))

for (i in 1:k) if (sum(index==i)==0) index[c(sample(1:numobs,q.z))]<-i

output<-list(H=H,lik=likelihood,w=w,Beta=Beta,psi=psi,sigma=sigma,ph.y=ph.y,index=index,z=z,phi=phi,bic=bic,elapsed=proc.time()-ptm,aic=aic) 
invisible(output)
}
