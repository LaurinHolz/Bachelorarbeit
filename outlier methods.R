statistical_test <- function(x){
  
  n <- length(x)
  mu <- mean(x)
  sigma <- sd(x)
  
  z_value <-  (x-mu)/sigma
  outlier_position <- which(abs(z_value)>3)
    
  if(length(outlier_position)==0){
    return("Kein Outlier gefunden")
  }else{
    return(outlier_position)
  }  
}

Mahalanobis_distance <- function(data,alpha=0.01){
  
  mu <- colMeans(data)
  sigma <- matrix(unlist(cov(data)),ncol(data),ncol(data))
  mahalanobis_dist <- numeric(nrow(data))
  
  for(i in 1:nrow(data)){
    mahalanobis_dist[i] <- as.numeric((data[i,]-mu)) %*% solve(sigma) %*% as.numeric(t(data[i,]-mu))
  }
  
  p_value <- pchisq(mahalanobis_dist,df=ncol(data), lower.tail=FALSE)
  outlier_position <- which(p_value<alpha)
  
  if(length(outlier_position)==0){
    return("Kein Outlier gefunden")
  }else{
    return(outlier_position)
  }
} 

MCD <- function(data,alpha=0.01,h=0.75){
  
  cutoff <- qchisq(p = 1-alpha, df = ncol(data)*h)
  
  val <- covMcd(data)
  mu <- as.numeric(val$center)
  sigma <- matrix(val$cov,ncol(data),ncol(data))
  
  dist <- numeric(nrow(data))
  for (i in 1:nrow(data)) {
    dist[i] <- as.numeric((data[i,]-mu)) %*% solve(sigma)%*% as.numeric(t(data[i,]-mu))
  }
  outlier_position <- which(dist>cutoff)
  
  if(length(outlier_position)==0){
    return("Kein Outlier gefunden")
  }else{
    return(outlier_position)
  }
}
