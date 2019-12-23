library(coxme)
library(survival)

coxlmm = function(y, x = NULL, g){
	y = y
	x = as.matrix(x)
	g = as.matrix(g)
	colnames(y) = c("time", "status")
	time = y$time
	status = y$status
	y_surv = Surv(time,status)
	K = g%*%t(g)
	eig = eigen(K)
	evals <- eig$values
	evecs <- eig$vectors
	K12 = evecs%*%diag(sqrt(evals))%*%t(evecs)
	if(is.null(x)){
		fit <- coxme(y_surv~(K12|1))
		u=K12%*%unlist(fit$frail)
		beta_g = t(g)%*%ginv(K)%*%u
		lp     = g%*%beta_g
		return(list(fit = fit, beta_g = beta_g))
	}else {
		x = as.matrix(x)
		fit <- coxme(y_surv~x+(K12|1))
		u=K12%*%unlist(fit$frail)
		beta_x = fit$coefficients
		beta_g = t(g)%*%ginv(K)%*%u
		lp     = x%*%beta_x + g%*%beta_g
		return(list(fit = fit, beta_g = beta_g, beta_x = beta_x, lp = lp))
	}
}


coxlmm_pred = function(y, x = NULL, g, coxlmm_fit){
	y = y
	g = as.matrix(g)
	time = y$time
	status = y$status
	y_surv = Surv(time,status)
	beta_g = coxlmm_fit$beta_g

	Cindex = function (y, lp){
		    time <- y[, 1]
		    status <- y[, 2]
		    x <- lp
		    n <- length(time)
		    ord <- order(time, -status)
		    time <- time[ord]
		    status <- status[ord]
		    x <- x[ord]
		    wh <- which(status == 1)
		    total <- concordant <- 0
		    for (i in wh) {
			for (j in ((i + 1):n)) {
			    tt <- (time[j] > time[i])
			    if (is.na(tt)) 
				tt <- FALSE
			    if (tt) {
				total <- total + 1
				if (x[j] < x[i]) 
				  concordant <- concordant + 1
				if (x[j] == x[i]) 
				  concordant <- concordant + 0.5
			    }
			}
		    }
		   return(list(cindex = concordant/total))
		}

	if(is.null(x)){
		lpnew = g%*%beta_g
		cindex = Cindex(y_surv,lpnew)$cindex
		return(list(lpnew = lpnew, cindex = cindex))
	}else {
		x = as.matrix(x)
		beta_x = coxlmm_fit$beta_x
		lpnew = x%*%beta_x + g%*%beta_g
		cindex = Cindex(y_surv,lpnew)$cindex
		return(list(lpnew = lpnew, cindex = cindex))
	}
}




PVE = function(y, x, g, coxlmm_fit){
	y = y
	g = as.matrix(g)
	x = as.matrix(x)
	time = y$time
	status = y$status
	y_surv = Surv(time,status)
	beta_g = coxlmm_fit$beta_g
	beta_x = coxlmm_fit$beta_x
	G1 = var(x%*%beta_x); G2 = var(g%*%beta_g); e = pi*pi/6
	PCE = G1/(G1+G2+e); PGE = G2/(G1+G2+e)
	return(list(G1 = G1, G2 = G2, PCE = PCE, PGE = PGE))
	}


