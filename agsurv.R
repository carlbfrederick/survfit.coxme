agsurv <- function (y, x, wt, risk, survtype, vartype) 
{
#y = response (Surv object)
#x = model matrix
#wt = weights
#risk = exp(linear prediction)
#survtype = method for handling exact ties (1=Kaplan Meier, 2=Breslow, **3 = Efron**)
#vartype = method for dealing with variance in coxph (1="greenwood", 2="aalen", **3="efron"**)
    nvar <- ncol(as.matrix(x))
    status <- y[, ncol(y)]
    dtime <- y[, ncol(y) - 1]
    death <- (status == 1)
    time <- sort(unique(dtime))
    nevent <- as.vector(rowsum(wt * death, dtime))
    ncens <- as.vector(rowsum(wt * (!death), dtime))
    wrisk <- wt * risk
    rcumsum <- function(x) rev(cumsum(rev(x)))
    nrisk <- rcumsum(rowsum(wrisk, dtime))
    irisk <- rcumsum(rowsum(wt, dtime))
    if (ncol(y) == 2) {
        temp2 <- rowsum(wrisk * x, dtime)
        xsum <- apply(temp2, 2, rcumsum)
    }
#COUNTING PROCESS SURV OBJECT    
    else {
        delta <- min(diff(time))/2
        etime <- c(sort(unique(y[, 1])), max(y[, 1]) + delta)
        indx <- approx(etime, 1:length(etime), time, method = "constant", 
            rule = 2, f = 1)$y
        esum <- rcumsum(rowsum(wrisk, y[, 1]))
        nrisk <- nrisk - c(esum, 0)[indx]
        irisk <- irisk - c(rcumsum(rowsum(wt, y[, 1])), 0)[indx]
        xout <- apply(rowsum(wrisk * x, y[, 1]), 2, rcumsum)
        xin <- apply(rowsum(wrisk * x, dtime), 2, rcumsum)
        xsum <- xin - (rbind(xout, 0))[indx, , drop = F]
    }
    ndeath <- rowsum(status, dtime)
    ntime <- length(time)
    if (survtype == 1) {
        indx <- (which(status == 1))[order(dtime[status == 1])]
        km <- .C(Cagsurv4, as.integer(ndeath), as.double(risk[indx]), 
            as.double(wt[indx]), as.integer(ntime), as.double(nrisk), 
            inc = double(ntime))
    }
#DEFAULT CASE FOR COXME    
    if (survtype == 3 || vartype == 3) {
        xsum2 <- rowsum((wrisk * death) * x, dtime)
        erisk <- rowsum(wrisk * death, dtime)
#HERE IS THE CALL TO C functions        
        tsum <- .C(Cagsurv5, as.integer(length(nevent)), as.integer(nvar), 
            as.integer(ndeath), as.double(nrisk), as.double(erisk), 
            as.double(xsum), as.double(xsum2), sum1 = double(length(nevent)), 
            sum2 = double(length(nevent)), xbar = matrix(0, 
                length(nevent), nvar))
    }
    haz <- switch(survtype, nevent/nrisk, nevent/nrisk, nevent * 
        tsum$sum1)
    varhaz <- switch(vartype, nevent/(nrisk * ifelse(nevent >= 
        nrisk, nrisk, nrisk - nevent)), nevent/nrisk^2, nevent * 
        tsum$sum2)
    xbar <- switch(vartype, (xsum/nrisk) * haz, (xsum/nrisk) * 
        haz, nevent * tsum$xbar)
    result <- list(n = nrow(y), time = time, n.event = nevent, 
        n.risk = irisk, n.censor = ncens, hazard = haz, cumhaz = cumsum(haz), 
        varhaz = varhaz, ndeath = ndeath, xbar = apply(matrix(xbar, 
            ncol = nvar), 2, cumsum))
    if (survtype == 1) 
        result$surv <- km$inc
    result
}
