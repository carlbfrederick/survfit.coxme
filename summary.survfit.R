summary.survfit <- function (object, times, censored = FALSE, scale = 1, extend = FALSE, 
    rmean = getOption("survfit.rmean"), ...) 
{
    fit <- object
    if (!inherits(fit, "survfit")) 
        stop("summary.survfit can only be used for survfit objects")
    if (is.null(rmean)) 
        rmean <- "none"
    temp <- survmean(fit, scale = scale, rmean)
    table <- temp$matrix
    rmean.endtime <- temp$end.time
    if (!missing(times)) {
        if (!is.numeric(times)) 
            stop("times must be numeric")
        times <- sort(times)
    }
    nsurv <- if (is.matrix(fit$surv)) 
        nrow(fit$surv)
    else length(fit$surv)
    if (is.null(fit$strata)) {
        nstrat <- 1
        stemp <- rep(1L, nsurv)
        strata.names <- ""
    }
    else {
        nstrat <- length(fit$strata)
        stemp <- rep(1:nstrat, fit$strata)
        strata.names <- names(fit$strata)
    }
    if (missing(times)) {
        if (censored) 
            indx1 <- seq(along = fit$time)
        else indx1 <- which(fit$n.event > 0)
    }
    else {
        cfun <- function(x, data) diff(c(0, cumsum(c(0, data))[x]))
        indx1 <- n.risk <- n.event <- newtimes <- vector("list", 
            nstrat)
        n.enter <- vector("list", nstrat)
        n.censor <- vector("list", nstrat)
        n <- length(stemp)
        for (i in 1:nstrat) {
            who <- (1:n)[stemp == i]
            stime <- fit$time[who]
            if (is.null(fit$start.time)) 
                mintime <- min(stime, 0)
            else mintime <- fit$start.time
            ptimes <- times[times >= mintime]
            if (!extend) {
                maxtime <- max(stime)
                ptimes <- ptimes[ptimes <= maxtime]
            }
            newtimes[[i]] <- ptimes
            ntime <- length(stime)
            temp1 <- approx(c(mintime - 1, stime), 0:ntime, 
                xout = ptimes, method = "constant", f = 0, rule = 2)$y
            indx1[[i]] <- ifelse(temp1 == 0, 1, 1 + who[pmax(1, 
                temp1)])
            n.event[[i]] <- cfun(temp1 + 1, fit$n.event[who])
            if (!is.null(fit$n.censor)) {
                n.censor[[i]] <- cfun(temp1 + 1, fit$n.censor[who])
                j <- who[ntime]
                last.n <- fit$n.risk[j] - (fit$n.event[j] + 
                  fit$n.censor[j])
            }
            else {
                last.n <- 0
            }
            if (ntime == 1) 
                temp1 <- rep(1, length(ptimes))
            else temp1 <- approx(stime, 1:ntime, xout = ptimes, 
                method = "constant", f = 1, rule = 2)$y
            n.risk[[i]] <- ifelse(ptimes > max(stime), last.n, 
                fit$n.risk[who[temp1]])
        }
        times <- unlist(newtimes)
        n.risk <- unlist(n.risk)
        n.event <- unlist(n.event)
        n.enter <- unlist(n.enter)
        n.censor <- unlist(n.censor)
        indx1 <- unlist(indx1)
    }
    if (length(indx1) == length(fit$time) && all(indx1 == seq(along = fit$time))) {
        temp <- object
        temp$time <- temp$time/scale
        temp$table <- table
        if (!is.null(temp$strata)) 
            temp$strata <- factor(stemp, labels = strata.names)
    }
    else if (missing(times)) {
        temp <- object
        temp$time <- temp$time[indx1]/scale
        temp$table <- table
        for (j in c("n.risk", "n.event", "n.censor", "n.enter", 
            "surv", "std.err", "cumhaz", "lower", "upper")) {
            zed <- temp[[j]]
            if (!is.null(zed)) {
                if (is.matrix(zed)) 
                  temp[[j]] <- zed[indx1, , drop = FALSE]
                else temp[[j]] <- zed[indx1]
            }
        }
        if (!is.null(temp$strata)) 
            temp$strata <- factor(stemp[indx1], levels = 1:nstrat, 
                labels = strata.names)
    }
    else {
        temp <- list(n = object$n, time = times/scale, n.risk = n.risk, 
            n.event = n.event, conf.int = fit$conf.int, type = fit$type, 
            table = table)
        if (!is.null(n.censor)) 
            temp$n.censor <- n.censor
        if (!is.null(n.enter)) 
            temp$n.enter <- n.enter
        if (!is.null(fit$start.time)) 
            temp$start.time <- fit$start.time
        if (is.matrix(fit$surv)) {
            temp$surv <- rbind(1, fit$surv)[indx1, , drop = FALSE]
            if (!is.null(fit$std.err)) 
                temp$std.err <- rbind(0, fit$std.err)[indx1, 
                  , drop = FALSE]
            if (!is.null(fit$lower)) {
                temp$lower <- rbind(1, fit$lower)[indx1, , drop = FALSE]
                temp$upper <- rbind(1, fit$upper)[indx1, , drop = FALSE]
            }
            if (!is.null(fit$cumhaz)) 
                temp$cumhaz <- rbind(0, fit$cumhaz)[indx1, , 
                  drop = FALSE]
        }
        else {
            temp$surv <- c(1, fit$surv)[indx1]
            if (!is.null(fit$std.err)) 
                temp$std.err <- c(0, fit$std.err)[indx1]
            if (!is.null(fit$lower)) {
                temp$lower <- c(1, fit$lower)[indx1]
                temp$upper <- c(1, fit$upper)[indx1]
            }
            if (!is.null(fit$cumhaz)) 
                temp$cumhaz <- c(0, fit$cumhaz)[indx1]
        }
        if (!is.null(fit$strata)) {
            scount <- unlist(lapply(newtimes, length))
            temp$strata <- factor(rep(1:nstrat, scount), levels = 1:nstrat, 
                labels = strata.names)
        }
        if (length(rmean.endtime) > 0 && !is.na(rmean.endtime)) 
            temp$rmean.endtime <- rmean.endtime
        temp$call <- fit$call
        if (!is.null(fit$na.action)) 
            temp$na.action <- fit$na.action
    }
    if (!is.null(temp$std.err)) 
        temp$std.err <- temp$std.err * temp$surv
    class(temp) <- "summary.survfit"
    temp
}
