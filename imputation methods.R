#Achtung: Lade und installiere zuerst Packages in "Numerische Experimente"

#Mittelwertimputation
mean_imp <- function(x){
  
  p <- ncol(x)
  n <- nrow(x)
  m <- colMeans(x,na.rm = TRUE)
  
  for (i in 1:p) {
    position <- which(is.na(x[,i]))
    x[position,i] <- m[i]
  }
  return(x)
}

#konditionale Mittelwertimputation
conditionalMean_imp <- function(x){
  
  p <- ncol(x)
  n <- nrow(x)
  m <- colMeans(x,na.rm = TRUE)
  
  #Mittelwertimputation für fehlenden ersten/letzten Eintrag
  for (i in 1:p) {
    
    if(is.na(x[1,i]))
      x[1,i] <- m[i]
    
    if(is.na(x[n,i]))
      x[n,i] <- m[i]
  }
  
  for(i in 1:p){
    col <- as.matrix(x[,i])
    bol <- FALSE
    md <- 0
    for (j in 2:(length(col)-1)) {
      #Fall 1: Eintrag ist beobachtet
      if(!is.na(col[j])){
        next
      #Fall 2: Eintrag ist unbeobachtet
      }else{
        md <- md+1
        #Fahre solange fort, bis nächster beobachteter Wert erreicht wird 
        while (!bol) {
          if(!is.na(col[j+md])){
            bol <-  TRUE
          }else{
            md <- md+1
          }
        }
        #Bilde den konditionalen Mittelwert
        col[j:(j+md-1)] <- (col[j-1]+col[j+md])/2
      }
      
      #Fahre in nächster Spalte fort, falls alle Einträge vollständig sind
      if(j+md < length(col)){
        j <- j+md-1
      }else{
        break
      }
    md <- 0
    bol <- FALSE
    }
    x[,i] <- col
  }
  return(x)
}

#Regressionsimputation
regression_imp <- function(x){
  
  if (is.data.frame(x) == TRUE) {
    x <- as.matrix(x)
  }
  
  p <- ncol(x)
  n <- nrow(x)
  
  for (i in 2:p) {
    position <- which(is.na(x[,i]))#Position der fehlenden Beobachtungen
    age <- which(is.na(x[,i]))+x[1,1]-1#Alter der fehlenden Beobachtungen
    fit <- lm(x[,i]~x[,1])#Lineare Regression mit Alter und aktueller Variable
    tmp <- summary(fit)
    coff1 <- tmp$coefficients[[1]]#Erhalte die Koeffizienten aus Regression
    coff2 <- tmp$coefficients[[2]]
    x[position,i] <- coff2*age+coff1#Stelle lineare Gleichung zur Berechnung auf
  }
  x <- as.data.frame(x)
  return(x)
}

#EM Algorithmus für multivariate Normalverteilungen
#Erklärung der verwendeten Funktionen in Kapitel 3.5: Implementierung
em_algo <- function(data,imputations=1,maxIt=1000,criterion=0.0001
                    ,rseed,start.mean,start.var){
  
  if(!is.matrix(data)){
    if(is.data.frame(data)){
      data <- as.matrix(data)
    }else{stop("data must be a data frame or matrix")
      
    }
  }
  
  norm::rngseed(rseed)
  data_imp <- vector("list",imputations)
  s <- norm::prelim.norm(data)
  start_list <- list(start.mean,start.var)
  start <- norm::makeparam.norm(s,start_list)
  theta <- norm::em.norm(s,start,maxIt,showits = FALSE)
  
  for (i in 1:imputations) {
    data_imp[[i]] <- as.data.frame(norm::imp.norm(s,theta,data))
  }
  
  return(data_imp)
}

#Die folgenden Funktionen werden für den Nachweis der Eigenschaften des 
#EM-Algorithmus benötigt. Die folgende Funktion "em_algorithm" stellt den
#normalen EM-Algorithmus dar. In der darin aufgerufenen Funktion "em" werden
#zusätzliche Informationen, wie die Likelihood-Funktion ausgewertet
em_algorithm <- function(data,imputations=10,maxIt=1000,criterion=0.0001,rseed,start.mean,start.var){
  
  if (is.data.frame(data) == FALSE) {
    stop("obsData argumment must be a data frame.")
  }
  
  norm::rngseed(rseed)
  s <- norm::prelim.norm(as.matrix(data))
  thetalist <- list(start.mean,start.var)
  start <- norm::makeparam.norm(s,thetalist)
  imps <- vector("list", imputations)

  #Berechnung der Thetamatrix, Speichere alle Iterationen von Theta
  #zur Berechnung der Konvergenzodnung
  thetahat <- em(s,start,showits = TRUE,maxIt,criterion)
  theta_star <- thetahat[[1]]
  theta_convergence <- thetahat[[2]]
  showtheta <- norm::getparam.norm(s,theta_star)
  
  for (i in 1:imputations) {
    imps[[i]] <- as.data.frame(norm::imp.norm(s, theta_star, 
                                              as.matrix(data)))
  }
  
  lambda <- numeric(length(theta_convergence))
  likelihood <- numeric(length(theta_convergence))
  
  #Berechnung der Konvergenzordnung
  for (t in 2:(length(lambda)-1)) {
    
    lambda[t] <- max(abs(theta_convergence[[t+1]]-theta_convergence[[t]]))/
               max(abs(theta_convergence[[t]]-theta_convergence[[t-1]]))
    
  }
  
  #Auswertung der Likelihood-Funktion
  for (t in 1:length(likelihood)){
    
    likelihood[t] <- norm::loglik.norm(s,theta_convergence[[t]])
  }
  
  iterations <- thetahat[[3]]
  
  #Rückgabe der imputierten Datensätze, Konvergenzordnung, Likelihood 
  #Auswertung, Anzahl Iterationen, Optimales Theta
  output <- list(imps,lambda,likelihood,iterations,showtheta)
  return(output)
}

em <- function (s, start, showits = TRUE, maxits = 1000, criterion = 1e-04, 
          prior){
  
  theta_convergence <- vector("list",maxits)
  
  s$x <- .na.to.snglcode(s$x, 999)
  if (missing(start)) {
    start <- .Fortran("stvaln", s$d, numeric(s$d), s$p, 
                      s$psi, PACKAGE = "norm")[[2]]
  }
  if (missing(prior)) {
    mle <- as.integer(1)
    tau <- numeric(1)
    m <- numeric(1)
    mu0 <- numeric(s$p)
    lambdainv <- matrix(0, s$p, s$p)
  }
  if (!(missing(prior))) {
    mle <- as.integer(0)
    tau <- as.numeric(prior[[1]])
    m <- as.numeric(prior[[2]])
    mu0 <- as.numeric(prior[[3]])
    lambdainv <- as.numeric(prior[[4]])
  }
  tmp <- as.integer(numeric(s$p))
  tobs <- .Fortran("tobsn", s$d, numeric(s$d), s$p, s$psi, 
                   s$n, s$x, s$npatt, s$r, s$mdpst, s$nmdp, tmp, PACKAGE = "norm")[[2]]
  it <- 0
  converged <- FALSE
  if (showits) 
    cat(paste("Iterations of EM:", "\n"))
  while ((!converged) & (it < maxits)) {
    old <- start
    start <- .Fortran("emn", s$d, old, start, tobs, s$p, 
                      s$psi, s$n, s$x, s$npatt, s$r, s$mdpst, s$nmdp, 
                      tmp, tmp, numeric(s$p), mle, tau, m, mu0, lambdainv, 
                      PACKAGE = "norm")[[3]]
    it <- it + 1
    theta_convergence[[it]] <- start
    
    if (showits) 
      cat(paste(format(it), "...", sep = ""))
    #Abbruchkriterium mit Thetamatrix (Alternativ Zeile 234: Likelihoodfunktion)
      converged <- max(abs(start-old)) <= criterion
      #norm::loglik.norm(s,start)-norm::loglik.norm(s,old)
  }
  if (showits) 
    cat("\n")
  
  conv <-  theta_convergence[-which(sapply(theta_convergence, is.null))]
  output <- list(start,conv,it)
  output
}
  
  
