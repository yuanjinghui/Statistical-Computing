#Call data
data <- read.table("C://Onedrive/OneDrive - Knights - University of Central Florida/UCF/Courses/Statistical Computing/Project1/pb2.txt")

#Convert the data to matrix form
Data <- as.matrix(data, ncol=5)
Y <- Data[,1]
Y<-as.vector(Y)
X <- as.matrix(Data[,2:5],ncol=4)
N <- length(Y)
for (i in 1:N){
  if (Y[i] > 1){
    Y[i]<--1
  }
}

## Defining the Gaussian kernel
rbf_kernel <- function(x1,x2,gamma){
  K<-exp(-(1/gamma^2)*t(x1-x2)%*%(x1-x2))
  return(K)
}


## Generate the matrix A and B
gamma<-1.5
Dm<-matrix(0,N,N)
for(i in 1:N){
  for(j in 1:N){
    Dm[i,j]<-Y[i]*Y[j]*rbf_kernel(X[i,],X[j,],gamma)
  }
}
H<-Dm+diag(N)*(1/gamma)+diag(N)*1e-12 # adding a very small number to the diag, some trick
d2<-as.vector(rep(1,N))


## Conjugate Gradient Method to solving the problem Ax=B
conjugate_gradient_method<-function (A,B,N){
  i<-1
  x<-matrix(0,62,N)
  r<-matrix(0,62,N)
  p<-matrix(0,62,N)
  beta<-numeric(N)
  lamda<-numeric(N)
  r[,1]<-B
  while (t(r[,i])%*%r[,i]>0){
    i<-i+1
    if (i==2){
      p[,i]<-r[,i-1]
    }
    else {
      beta[i]<-t(r[,i-1])%*%r[,i-1]/t(r[,i-2])%*%r[,i-2]
      p[,i]<-r[,i-1]+beta[i]*p[,i-1]
    }
    lamda[i]<-t(r[,i-1])%*%r[,i-1]/t(p[,i])%*%A%*%p[,i]
    x[,i]<-x[,i-1]+lamda[i]*p[,i]
    r[,i]<-r[,i-1]-(lamda[i]*A)%*%p[,i]
  }
  return(x[,i])
}

## Training the LSSVM
lssvmtrain<-function(H, Y, d2, step){
  nta<-conjugate_gradient_method(H,Y,step)
  vta<-conjugate_gradient_method(H,d2,step)
  s<-t(Y)%*%nta
  b<-(t(nta)%*%d2)/s
  alpha<-vta-nta*b
  alpha<-as.vector(alpha)
  
  list(alpha=alpha, b=b)
}

## Calculate the alpha and b
model<-lssvmtrain(H, Y, d2, 5000)

model

### Predict the class of an object X
lssvmpredict <- function(Xv,Yv,x,model,N){
  alpha<-model$alpha
  b<-model$b
  gamma<-1.5
  ayK<-numeric(N)
  for(i in 1:N){
    ayK[i]<-alpha[i]*Yv[i]*rbf_kernel(Xv[i,],x,gamma)
  }
  result <- sign(sum(ayK)+b)
  return(result)
}

z <- c(16,13,16,14)
lssvmpredict(X,Y,z,model,62)

z <- c(18,17,33,26)
lssvmpredict(X,Y,z,model,62)



################################
###  Problem 1-2   #############
################################

QR_Com<-function(X,Y,gamma){
  A_up<-append(0,Y)
  N<-length(Y)
  omega<-matrix(0,N,N)
  for(i in 1:N){
    for(j in 1:N){
      omega[i,j]<-Y[i]*Y[j]*rbf_kernel(X[i,],X[j,],gamma)
    }
  }
  H<-omega+diag(N)*(1/gamma)+diag(N)*1e-12 # adding a very small number to the diag, some trick
  A_down<-cbind(Y,H)
  A<-rbind(A_up,A_down)
  A<-as.matrix(A)
  
  qr<-qr(A)
  Q<-qr.Q(qr)
  R<-qr.R(qr)
  list(Q=Q,R=R)
}

###Givens function
Givens<-function(a,b){
  if (b==0){
    c<-1
    s<-0
  }
  else {
    if (abs(b)>=abs(a)){
      t<-(-a/b)
      s<-1/sqrt(1+t^2)
      c<-s*t
    }
    else {
      t<-(-b/a)
      c<-1/sqrt(1+t^2)
      s<-c*t
    }
  }
  list(c=c,s=s)
}

###QR update
qrupdate<-function(X,Y,X_add,Y_add,gamma){
  n<-nrow(X)
  u_t<-matrix(0,1,n)
  for (i in 1:n) {
    u_t[i]<-Y_add*Y[i]*rbf_kernel(X_add,X[i,],gamma)
  }
  u<-append(Y_add,u_t)
  
  Q<-QR_Com(X,Y,gamma)$Q
  R<-QR_Com(X,Y,gamma)$R
  
  n_q<-n+1
  Q_1<-rbind(Q,matrix(0,1,n_q))
  Q_1<-cbind(Q_1,matrix(0,n_q+1,1))
  Q_1[n_q+1,n_q+1]<-1
  
  
  ##Update Q and R by adding a new row
  c<-numeric(n_q)
  s<-numeric(n_q)
  for (j in 1:n_q) {
    c[j]<-Givens(R[j,j],u[j])$c
    s[j]<-Givens(R[j,j],u[j])$s
    R[j,j]<-c[j]*R[j,j]-s[j]*u[j]
    ##Update jth row of R and u
    if (j<n_q){
      t1<-R[j,(j+1):n_q]
      t2<-u[(j+1):n_q]
      R[j,(j+1):n_q]=c[j]*t1-s[j]*t2
      u[(j+1):n_q]=s[j]*t1+c[j]*t2
    }
    ##Update Q
    t1_q<-Q_1[1:(n_q+1),j]
    t2_q<-Q_1[1:(n_q+1),(n_q+1)]
    Q_1[1:(n_q+1),j]<-c[j]*t1_q-s[j]*t2_q
    Q_1[1:(n_q+1),(n_q+1)]<-s[j]*t1_q+c[j]*t2_q
  }
  R_1<-rbind(R,matrix(0,1,n_q))
  ##Update Q and R by adding a new column
  u2<-append(u_t,Y_add*Y_add*rbf_kernel(X_add,X_add,gamma))
  u2<-append(Y_add,u2)
  u2<-t(Q_1)%*%u2
  
  m<-nrow(R_1)
  k<-ncol(R_1)+1
  
  c1<-numeric(k-m+1)
  s1<-numeric(k-m+1)
  for (i in m:k) {
    c1[i]<-Givens(u2[i-1],u2[i])$c
    s1[i]<-Givens(u2[i-1],u2[i])$s
    u2[i]<-c1[i]*u2[i-1]-s1[i]*R_1[(i-1),(i-1)]
    
    R_1<-cbind(R_1,matrix(u2,m,1))
    Q_1[1:m,(i-1):i]<-Q_1[1:m,(i-1):i]%*%matrix(c(c[j],-s[j],s[j],c[j]),2,2)
  }
  list(Q_1=Q_1,R_1=R_1)
}

##Testing the alghrithm to generate the update Q and R
X1<-head(X,5)
Y1<-head(Y,5)
X_add<-X[6,]
Y_add<-Y[6]
Q_1<-qrupdate(X1,Y1,X_add,Y_add,1.5)$Q_1
R_1<-qrupdate(X1,Y1,X_add,Y_add,1.5)$R_1


X2<-head(X,6)
Y2<-head(Y,6)
Q<-QR_Com(X2,Y2,1.5)$Q
R<-QR_Com(X2,Y2,1.5)$R

Q_diff<-Q_1-Q
R_diff<-R_1-R
Q_diff
R_diff

###Problem1-3
Solution<-function(X1,Y1,X_add,Y_add,gamma){
  Q_1<-qrupdate(X1,Y1,X_add,Y_add,gamma)$Q_1
  R_1<-qrupdate(X1,Y1,X_add,Y_add,gamma)$R_1
  Y_new<-append(Y1,Y_add)
  n<-length(Y_new)
  B<-t(Q_1)%*%append(0,rep(1,n))
  x<-solve(R_1,B)
  x<-as.vector(x)
  b<-x[1]
  alpha<-x[2:length(x)]
  list(b=b,alpha=alpha)
}

###Compare the results with the 1-1
X1<-head(X,61)
Y1<-head(Y,61)
X_add<-X[62,]
Y_add<-Y[62]
solution<-Solution(X1,Y1,X_add,Y_add,1.5)

solution$b
solution$alpha

alpha_diff<-solution$alpha-model$alpha
alpha_diff

b_diff<-solution$b-model$b
b_diff

################################
###    Problem 2   #############
################################

#Call data
data2 <- read.table("C://Onedrive/OneDrive - Knights - University of Central Florida/UCF/Courses/Statistical Computing/Project1/Boston.txt")
Data2<-data.frame(data2)
train.rows<-sample(nrow(Data2),300)
train<-Data2[train.rows,]
Y<-train[,14]
X <- as.matrix(train[,1:13],ncol=13)

require('quadprog')
## Defining the Gaussian kernel
rbf_kernel <- function(x1,x2,gamma){
  K<-exp(-(1/gamma^2)*t(x1-x2)%*%(x1-x2))
  return(K)
}

svrtrain <- function(X,Y,C=Inf, gamma=1.5,esp=1e-10){
  N<-length(Y)
  Dm<-matrix(0,N,N)
  X<-as.matrix(X)
  Y<-as.vector(Y)
  for(i in 1:N){
    for(j in 1:N){
      Dm[i,j]<-rbf_kernel(X[i,],X[j,],gamma)
    }
  }
  Dm<-Dm+diag(N)*1e-12 # adding a very small number to the diag, some trick
  
  dv<-t(Y)
  meq<-1
  Am<-cbind(matrix(rep(1,N),N),diag(N)) 
  if(C!=Inf){
    Am<-cbind(Am,-1*diag(N))
    bv<-append(0,rep(-C,2*N))
  }
  alpha_org<-solve.QP(Dm,dv,Am,meq=meq,bvec=bv)$solution
  indx<-which(alpha_org>esp,arr.ind=TRUE)
  alpha<-alpha_org[indx]
  nSV<-length(indx)
  if(nSV==0){
    throw("QP is not able to give a solution for these data points")
  }
  Xv<-X[indx,]
  Yv<-Y[indx]
  Yv<-as.vector(Yv)
  ## choose one of the support vector to compute b. Instead of using an arbitrary Support 
  ##Vector xs, it is better to take an average over all of the Support Vectors in S
  b <- numeric(nSV)
  ayK <- numeric(nSV)
  for (i in 1:nSV){
    for (m in 1:nSV){
      ayK[m] <- alpha[m]%*%rbf_kernel(Xv[m,],Xv[i,],gamma)
    }
    b[i]<-Yv[i]-sum(ayK)
  }
  w0 <- mean(b)
  #list(alpha=alpha, wstar=w, b=w0, nSV=nSV, Xv=Xv, Yv=Yv, gamma=gamma)
  list(alpha=alpha, b=w0, nSV=nSV, Xv=Xv, Yv=Yv, gamma=gamma)
}
model <-svrtrain(X,Y,C=12,gamma=24,esp=1e-10)


### Predict the Home value of each sample
svrpredict <- function(x,model){
  alpha<-model$alpha
  b<-model$b
  Yv<-model$Yv
  Xv<-model$Xv
  nSV<-model$nSV
  gamma<-model$gamma

  N<-nrow(x)
  ayK <- numeric(nSV)
  result<-numeric(N)
  for(k in 1:N){
    for(i in 1:nSV){
      ayK[i]<-alpha[i]*rbf_kernel(Xv[i,],x[k,],gamma)
    }
    result[k] <- sum(ayK)+b
  }
  return(result)
}

### Predict Y for testing dataset
test<-Data2[-train.rows,]
Yt<-test[,14]
Xt<- as.matrix(test[,1:13],ncol=13)
Predict<-svrpredict(Xt,model)

### Evaluate the performance
N<-length(Yt)
Tot_MAE<-0
Tot_RMSE<-0
for (i in 1:N){
  Tot_MAE<-Tot_MAE+abs(Yt[i]-Predict[i])
  Tot_RMSE<-Tot_RMSE+(Yt[i]-Predict[i])^2
}
MAE<-Tot_MAE/N
RMSE<-sqrt(Tot_RMSE/N)

MAE
RMSE

