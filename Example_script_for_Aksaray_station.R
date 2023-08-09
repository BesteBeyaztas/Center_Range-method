# Install following packages
library(MASS)
library(VAR.etp)
library(pracma)
library(psych)
library(tseries)
library(vars)
library(gdata)

# Number of forecast horizon and number of bootstrap
h <- 12
B <- 2000

# Load the data
max_data <- unmatrix(read.table("aksaray_max.txt", header=TRUE),byrow=TRUE)
min_data <- unmatrix(read.table("aksaray_min.txt", header=TRUE), byrow=TRUE)
mean_data <- unmatrix(read.table("aksaray_mean.txt", header=FALSE), byrow=TRUE)

# Data preprocessing
max.data <- c(max_data[1:(length(max_data)-h)])
min.data <- c(min_data[1:(length(min_data)-h)])

temp <- cbind(min.data, max.data)
nd <- nrow(temp)
temp_diff <- diff(temp,lag=12)

center_data <- (temp_diff[,1]+temp_diff[,2])/2
range_data <- (temp_diff[,2]-temp_diff[,1])/2

data_all <- cbind(center_data, range_data)
colnames(data_all) = c("Center", "Renge")

# Future values to be predicted
yh <- cbind(c(min_data[(length(min_data)-h+1):length(min_data)]),
c(max_data[(length(max_data)-h+1):length(max_data)]))
rownames(yh) <- 1:12

# Future center data
yh.c <- cbind(c(mean_data[(length(mean_data)-h+1):length(mean_data)]))

# ADF test
adf.test(min.data, alternative="stationary")
adf.test(max.data, alternative="stationary")

# Ljung-Box test and sample statistics
mardia(data_all)
describe(data_all)

est <- VAR(data_all, lag.max=14, type="const")
normality.test(est)
jarque.bera.test(min.data)
jarque.bera.test(max.data)


t = nrow(data_all)
N = 2
p = est$p
Nc <- N*p

# Function to calculate the roots for checking stationary condition
rootsPope <- function(Popecoef){
  cmpnion <- matrix(0, nrow=N*p, ncol=N*p)
  cmpnion[1:N, 1:(N*p)] <- Popecoef
  if (p > 1) {
    s <- 0
    for (l in (N + 1):(N * p)) {
      s <- s + 1
      cmpnion[l, s] <- 1
    }
  }
  rts <- eigen(cmpnion)$values
  rts <- Mod(rts)
  return(rts)
}

var.fitPope <- VAR.Pope(data_all, p, type="const")

phihatPope <- list()
k <- 1
for(i in 1:p){
  phihatPope[[i]] <- matrix(c(var.fitPope$coef[1,1:Nc][k], var.fitPope$coef[2,1:Nc][k], 
                              var.fitPope$coef[1,1:Nc][k+1], var.fitPope$coef[2,1:Nc][k+1]),2,2)
  k <- k+2
}


Amat <- c(do.call(cbind, phihatPope))
roots_res <- rootsPope(Amat)

#Residuals and re-scaled residuals
res1 <- var.fitPope$resid
res <- (res1-colMeans(res1))/apply(res1,2,sd) #forward residuals

boot.h_res <- list()
for(i in 1:h)
  boot.h_res[[i]] <- matrix(NA, ncol=N, nrow=B)

boot.c <- matrix(NA, ncol=h, nrow=B)


for(boot in 1:B){
  try({
    repeat{
      ###FRP Bootstrap Estimation(Begin)
      boot.data <- matrix(NA, ncol=N, nrow=(t+p))
      boot.data[1:p,] <- 0
      boot.err <- res[sample((t-p), (t+p), replace=TRUE),]
      
      for(i in (p+1):(t+p)){
        boot.data[i,] <- matrix(c(var.fitPope$coef[1,(Nc+1)], var.fitPope$coef[2,(Nc+1)]), ncol=1) +
          phihatPope[[1]] %*% boot.data[(i-1),] + phihatPope[[2]] %*% boot.data[(i-2),] + phihatPope[[3]] %*% boot.data[(i-3),]+
          phihatPope[[4]] %*% boot.data[(i-4),] + phihatPope[[5]] %*% boot.data[(i-5),] + phihatPope[[6]] %*% boot.data[(i-6),]+
          phihatPope[[7]] %*% boot.data[(i-7),] + phihatPope[[8]] %*% boot.data[(i-8),] + phihatPope[[9]] %*% boot.data[(i-9),]+
          phihatPope[[10]] %*% boot.data[(i-10),] + phihatPope[[11]] %*% boot.data[(i-11),] + phihatPope[[12]] %*% boot.data[(i-12),]+
          phihatPope[[13]] %*% boot.data[(i-13),] + boot.err[(i-p),]        
      }
      
      boot.data <- boot.data[(p+1):(t+p),]
      
      
      boot.fit <- VAR.Pope(boot.data, p, type="const")
      
      boot.phihat <- list()
      l <- 1
      for(j in 1:p){
        boot.phihat[[j]] <- matrix(c(boot.fit$coef[1,1:Nc][l], boot.fit$coef[2,1:Nc][l], 
                                     boot.fit$coef[1,1:Nc][l+1], boot.fit$coef[2,1:Nc][l+1]),2,2)
        l <- l+2
      }
      
      
      Bmat <- c(do.call(cbind, boot.phihat))
      bootroots <- rootsPope(Bmat)
      
      if(max(abs(bootroots)) < 1) break # Statioanry Condition for Bootstrap Step
    }
  }, silent=TRUE)
  ###FRP Bootstrap (End)
  
  
  #Step3 Y_{T+h} FRP Bootstrap (Begin)
  
  boot.h1 <- matrix(NA, ncol=N, nrow=(h+p))
  boot.h1[1:p,] <- data_all[(t-p+1):t,]
  h.err <- res[sample((t-p), (h+p), replace=TRUE),]
  
  for(j in (p+1):(h+p)){
    boot.h1[j,] <- matrix(c(boot.fit$coef[1,Nc], boot.fit$coef[2,Nc]), ncol=1) +
      boot.phihat[[1]] %*% boot.data[(j-1),] + boot.phihat[[2]] %*% boot.data[(j-2),] + boot.phihat[[3]] %*% boot.data[(j-3),]+
      boot.phihat[[4]] %*% boot.data[(j-4),] + boot.phihat[[5]] %*% boot.data[(j-5),] + boot.phihat[[6]] %*% boot.data[(j-6),]+
      boot.phihat[[7]] %*% boot.data[(j-7),] + boot.phihat[[8]] %*% boot.data[(j-8),] + boot.phihat[[9]] %*% boot.data[(j-9),]+
      boot.phihat[[10]] %*% boot.data[(j-10),] + boot.phihat[[11]] %*% boot.data[(j-11),] + boot.phihat[[12]] %*% boot.data[(j-12),]+
      boot.phihat[[13]] %*% boot.data[(j-13),] + h.err[(j-p),]
  }
  
  boot.h1 <- boot.h1[(p+1):(h+p),] 
  boot.h <- cbind(c(boot.h1[,1] - boot.h1[,2]), c(boot.h1[,1] + boot.h1[,2])) 
 
boot.min <- boot.data[,1]
boot.max <- boot.data[,2]

orig.boot.min <- diffinv(x=boot.min,lag=12,xi=min.data[1:12])
orig.boot.max <- diffinv(x=boot.max,lag=12,xi=max.data[1:12])

min_vals <- orig.boot.min[(nd-12+1):nd]
max_vals <- orig.boot.max[(nd-12+1):nd]

orig.min=diffinv(x=boot.h[,1],lag=12,xi=min_vals) 
orig.max=diffinv(x=boot.h[,2],lag=12,xi=max_vals) 

results <- cbind(orig.min[13:24],orig.max[13:24])
boot.c[boot,] <- (results[,1] + results[,2])/2

for(ij in 1:h){
  boot.h_res[[ij]][boot,] <- results[ij,]
}
  # Y_{T+h} FRP Bootstrap (End)
}

### FRP Bootstrap Quantiles
q11 <- quantile(boot.h_res[[1]][,1], c(0.025,0.975), na.rm=T)
q12 <- quantile(boot.h_res[[1]][,2], c(0.025,0.975), na.rm=T)

q21 <- quantile(boot.h_res[[2]][,1], c(0.025,0.975), na.rm=T)
q22 <- quantile(boot.h_res[[2]][,2], c(0.025,0.975), na.rm=T)

q31 <- quantile(boot.h_res[[3]][,1], c(0.025,0.975), na.rm=T)
q32 <- quantile(boot.h_res[[3]][,2], c(0.025,0.975), na.rm=T)

q41 <- quantile(boot.h_res[[4]][,1], c(0.025,0.975), na.rm=T)
q42 <- quantile(boot.h_res[[4]][,2], c(0.025,0.975), na.rm=T)

q51 <- quantile(boot.h_res[[5]][,1], c(0.025,0.975), na.rm=T)
q52 <- quantile(boot.h_res[[5]][,2], c(0.025,0.975), na.rm=T)

q61 <- quantile(boot.h_res[[6]][,1], c(0.025,0.975), na.rm=T)
q62 <- quantile(boot.h_res[[6]][,2], c(0.025,0.975), na.rm=T)

q71 <- quantile(boot.h_res[[7]][,1], c(0.025,0.975), na.rm=T)
q72 <- quantile(boot.h_res[[7]][,2], c(0.025,0.975), na.rm=T)

q81 <- quantile(boot.h_res[[8]][,1], c(0.025,0.975), na.rm=T)
q82 <- quantile(boot.h_res[[8]][,2], c(0.025,0.975), na.rm=T)

q91 <- quantile(boot.h_res[[9]][,1], c(0.025,0.975), na.rm=T)
q92 <- quantile(boot.h_res[[9]][,2], c(0.025,0.975), na.rm=T)

q101 <- quantile(boot.h_res[[10]][,1], c(0.025,0.975), na.rm=T)
q102 <- quantile(boot.h_res[[10]][,2], c(0.025,0.975), na.rm=T)

q111 <- quantile(boot.h_res[[11]][,1], c(0.025,0.975), na.rm=T)
q112 <- quantile(boot.h_res[[11]][,2], c(0.025,0.975), na.rm=T)

q121 <- quantile(boot.h_res[[12]][,1], c(0.025,0.975), na.rm=T)
q122 <- quantile(boot.h_res[[12]][,2], c(0.025,0.975), na.rm=T)
