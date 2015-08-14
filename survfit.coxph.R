function (formula, newdata, se.fit = TRUE, conf.int = 0.95, 
    individual = FALSE, type, vartype, conf.type = c("log", 
        "log-log", "plain", "none"), censor = TRUE, id, na.action = na.pass, 
    ...) 
{
##PREP FUNCTION  
    Call <- match.call()
    Call[[1]] <- as.name("survfit")
    object <- formula
    if (!is.null(attr(object$terms, "specials")$tt)) 
        stop("The survfit function can not yet process coxph models with a tt term")
####setting type    
    if (missing(type)) {
        temp1 <- c("exact", "breslow", "efron")
        survtype <- match(object$method, temp1)
    }
    else {
        temp1 <- c("kalbfleisch-prentice", "aalen", "efron", 
            "kaplan-meier", "breslow", "fleming-harrington", 
            "greenwood", "tsiatis", "exact")
        survtype <- match(match.arg(type, temp1), temp1)
        survtype <- c(1, 2, 3, 1, 2, 3, 1, 2, 3)[survtype]
    }
####setting vartype
    if (missing(vartype)) {
        vartype <- survtype
    }
    else {
        temp2 <- c("greenwood", "aalen", "efron", "tsiatis")
        vartype <- match(match.arg(vartype, temp2), temp2)
        if (vartype == 4) 
            vartype <- 2
    }
####do confidence intervalse?
    if (!se.fit) 
        conf.type <- "none"
    else conf.type <- match.arg(conf.type)
####Are there strata
    has.strata <- !is.null(attr(object$terms, "specials")$strata)
####GET the model frame    
    if (is.null(object$y) || is.null(object[["x"]]) || !is.null(object$call$weights) || 
        (has.strata && is.null(object$strata)) || !is.null(attr(object$terms, 
        "offset"))) {
        mf <- model.frame(object)
    }
    else mf <- NULL
    if (is.null(mf)) 
        y <- object[["y"]]
    else {
        y <- model.response(mf)
        y2 <- object[["y"]]
        if (!is.null(y2) && any(as.matrix(y2) != as.matrix(y))) 
            stop("Could not reconstruct the y vector")
    }
    if (is.null(object[["x"]])) 
        x <- model.matrix.coxph(object, data = mf)
    else x <- object[["x"]]
    n <- nrow(y)
    if (n != object$n[1] || nrow(x) != n) 
        stop("Failed to reconstruct the original data set")
    if (is.null(mf)) 
        wt <- rep(1, n)
    else {
        wt <- model.weights(mf)
        if (is.null(wt)) 
            wt <- rep(1, n)
    }
    type <- attr(y, "type")
    if (type != "right" && type != "counting") 
        stop("Cannot handle \"", type, "\" type survival data")
    missid <- missing(id)
    if (!missid) 
        individual <- TRUE
    else if (missid && individual) 
        id <- rep(0, n)
    else id <- NULL
    if (individual && missing(newdata)) {
        stop("the id and/or individual options only make sense with new data")
    }
    if (individual && type != "counting") 
        stop("The individual option is  only valid for start-stop data")
    if (is.null(mf)) 
        offset <- 0
    else {
        offset <- model.offset(mf)
        if (is.null(offset)) 
            offset <- 0
    }
    Terms <- object$terms
    if (!has.strata) 
        strata <- rep(0L, n)
    else {
        stangle <- untangle.specials(Terms, "strata")
        strata <- object$strata
        if (is.null(strata)) {
            if (length(stangle$vars) == 1) 
                strata <- mf[[stangle$vars]]
            else strata <- strata(mf[, stangle$vars], shortlabel = TRUE)
        }
    }
    if (has.strata) {
        temp <- attr(Terms, "specials")$strata
        factors <- attr(Terms, "factors")[temp, ]
        strata.interaction <- any(t(factors) * attr(Terms, "order") > 
            1)
    }
    if (is.null(x) || ncol(x) == 0) {
        x <- matrix(0, nrow = n)
        coef <- 0
        varmat <- matrix(0, 1, 1)
        risk <- rep(exp(offset - mean(offset)), length = n)
    }
    else {
#HAVE CHANGED THE FOLLOWING LINE
        varmat <- vcov(object)
        coef <- ifelse(is.na(object$coefficients), 0, object$coefficients)
        xcenter <- object$means
        if (is.null(object$frail)) {
            x <- scale(x, center = xcenter, scale = FALSE)
            risk <- c(exp(x %*% coef + offset - mean(offset)))
        }
#PROBABLY NEED TO CHANGE THESE LINES ----        
        else {
            keep <- !is.na(match(dimnames(x)[[2]], names(coef)))
            x <- x[, keep, drop = F]
            risk <- exp(object$linear.predictor)
            x <- scale(x, center = xcenter, scale = FALSE)
        }
#---- To these lines        
    }
    subterms <- function(tt, i) {
        dataClasses <- attr(tt, "dataClasses")
        predvars <- attr(tt, "predvars")
        oldnames <- dimnames(attr(tt, "factors"))[[1]]
        tt <- tt[i]
        index <- match(dimnames(attr(tt, "factors"))[[1]], oldnames)
        if (length(index) > 0) {
            if (!is.null(predvars)) 
                attr(tt, "predvars") <- predvars[c(1, index + 
                  1)]
            if (!is.null(dataClasses)) 
                attr(tt, "dataClasses") <- dataClasses[index]
        }
        tt
    }
    temp <- untangle.specials(Terms, "cluster")
    if (length(temp$vars)) 
        Terms <- subterms(Terms, -temp$terms)
    
    
    if (missing(newdata)) {
        mf2 <- as.list(object$means)
        names(mf2) <- names(object$coefficients)
        mf2 <- as.data.frame(mf2)
        found.strata <- FALSE
    }
    else {
#CHANGE HERE ... what to do when !is.null(newdata)      
        if (!is.null(object$frail)) 
            stop("Newdata cannot be used when a model has frailty terms")
        Terms2 <- Terms
        if (!individual) 
            Terms2 <- delete.response(Terms)
        if (is.vector(newdata, "numeric")) {
            if (individual) 
                stop("newdata must be a data frame")
            if (is.null(names(newdata))) {
                stop("Newdata argument must be a data frame")
            }
            newdata <- data.frame(as.list(newdata))
        }
        if (missid) {
            if (has.strata && !strata.interaction) {
                found.strata <- TRUE
                tempenv <- new.env(, parent = emptyenv())
                assign("strata", function(..., na.group, shortlabel, 
                  sep) list(...), envir = tempenv)
                assign("list", list, envir = tempenv)
                for (svar in stangle$vars) {
                  temp <- try(eval(parse(text = svar), newdata, 
                    tempenv), silent = TRUE)
                  if (!is.list(temp) || any(unlist(lapply(temp, 
                    class)) == "function")) 
                    found.strata <- FALSE
                }
                if (found.strata) 
                  mf2 <- model.frame(Terms2, data = newdata, 
                    na.action = na.action, xlev = object$xlevels)
                else {
                  Terms2 <- subterms(Terms2, -attr(Terms2, "specials")$strata)
                  if (!is.null(object$xlevels)) {
                    myxlev <- object$xlevels[match(attr(Terms2, 
                      "term.labels"), names(object$xlevels), 
                      nomatch = 0)]
                    if (length(myxlev) == 0) 
                      myxlev <- NULL
                  }
                  else myxlev <- NULL
                  mf2 <- model.frame(Terms2, data = newdata, 
                    na.action = na.action, xlev = myxlev)
                }
            }
            else {
                mf2 <- model.frame(Terms2, data = newdata, na.action = na.action, 
                  xlev = object$xlevels)
                found.strata <- has.strata
            }
        }
        else {
            tcall <- Call[c(1, match(c("id", "na.action"), names(Call), 
                nomatch = 0))]
            tcall$data <- newdata
            tcall$formula <- Terms2
#HERE create xlevels            
            tcall$xlev <- object$xlevels
            tcall[[1]] <- as.name("model.frame")
            mf2 <- eval(tcall)
            found.strata <- has.strata
        }
    }
    if (has.strata && found.strata) {
        temp <- untangle.specials(Terms2, "strata")
        strata2 <- strata(mf2[temp$vars], shortlabel = TRUE)
        strata2 <- factor(strata2, levels = levels(strata))
        if (any(is.na(strata2))) 
            stop("New data set has strata levels not found in the original")
        Terms2 <- Terms2[-temp$terms]
    }
    else strata2 <- factor(rep(0, nrow(mf2)))
    if (individual) {
        if (missing(newdata)) 
            stop("The newdata argument must be present when individual=TRUE")
        if (!missid) {
            id <- model.extract(mf2, "id")
            if (is.null(id)) 
                stop("id=NULL is an invalid argument")
        }
        else id <- rep(1, nrow(mf2))
        x2 <- model.matrix(Terms2, mf2)[, -1, drop = FALSE]
        if (length(x2) == 0) 
            stop("Individual survival but no variables")
        x2 <- scale(x2, center = xcenter, scale = FALSE)
        offset2 <- model.offset(mf2)
        if (length(offset2) > 0) 
            offset2 <- offset2 - mean(offset)
        else offset2 <- 0
        y2 <- model.extract(mf2, "response")
        if (attr(y2, "type") != type) 
            stop("Survival type of newdata does not match the fitted model")
        if (attr(y2, "type") != "counting") 
            stop("Individual=TRUE is only valid for counting process data")
        y2 <- y2[, 1:2, drop = F]
        newrisk <- exp(c(x2 %*% coef) + offset2)
#HERE        
        result <- survfitcoxph.fit(y, x, wt, x2, risk, newrisk, 
            strata, se.fit, survtype, vartype, varmat, id, y2, 
            strata2)
    }
    else {
        if (missing(newdata)) {
            if (has.strata && strata.interaction) 
                stop("Models with strata by covariate interaction terms require newdata")
            x2 <- matrix(0, nrow = 1, ncol = ncol(x))
            offset2 <- 0
        }
        else {
            offset2 <- model.offset(mf2)
            if (length(offset2) > 0) 
                offset2 <- offset2 - mean(offset)
            else offset2 <- 0
            x2 <- model.matrix(Terms2, mf2)[, -1, drop = FALSE]
            x2 <- scale(x2, center = xcenter, scale = FALSE)
        }
        newrisk <- exp(c(x2 %*% coef) + offset2)
        result <- survfitcoxph.fit(y, x, wt, x2, risk, newrisk, 
            strata, se.fit, survtype, vartype, varmat)
        if (has.strata && found.strata) {
            if (is.matrix(result$surv)) {
                nr <- nrow(result$surv)
                indx1 <- split(1:nr, rep(1:length(result$strata), 
                  result$strata))
                rows <- indx1[as.numeric(strata2)]
                indx2 <- unlist(rows)
                indx3 <- as.integer(strata2)
                for (i in 2:length(rows)) rows[[i]] <- rows[[i]] + 
                  (i - 1) * nr
                indx4 <- unlist(rows)
                temp <- result$strata[indx3]
                names(temp) <- row.names(mf2)
                new <- list(n = result$n[indx3], time = result$time[indx2], 
                  n.risk = result$n.risk[indx2], n.event = result$n.event[indx2], 
                  n.censor = result$n.censor[indx2], strata = temp, 
                  surv = result$surv[indx4], cumhaz = result$cumhaz[indx4])
                if (se.fit) 
                  new$std.err <- result$std.err[indx4]
                result <- new
            }
        }
    }
    if (!censor) {
        kfun <- function(x, keep) {
            if (is.matrix(x)) 
                x[keep, , drop = F]
            else if (length(x) == length(keep)) 
                x[keep]
            else x
        }
        keep <- (result$n.event > 0)
        if (!is.null(result$strata)) {
            temp <- factor(rep(names(result$strata), result$strata), 
                levels = names(result$strata))
            result$strata <- c(table(temp[keep]))
        }
        result <- lapply(result, kfun, keep)
    }
    if (se.fit) {
        zval <- qnorm(1 - (1 - conf.int)/2, 0, 1)
        if (conf.type == "plain") {
            temp1 <- result$surv + zval * result$std.err * result$surv
            temp2 <- result$surv - zval * result$std.err * result$surv
            result <- c(result, list(upper = pmin(temp1, 1), 
                lower = pmax(temp2, 0), conf.type = "plain", 
                conf.int = conf.int))
        }
        if (conf.type == "log") {
            xx <- ifelse(result$surv == 0, 1, result$surv)
            temp1 <- ifelse(result$surv == 0, 0 * result$std.err, 
                exp(log(xx) + zval * result$std.err))
            temp2 <- ifelse(result$surv == 0, 0 * result$std.err, 
                exp(log(xx) - zval * result$std.err))
            result <- c(result, list(upper = pmin(temp1, 1), 
                lower = temp2, conf.type = "log", conf.int = conf.int))
        }
        if (conf.type == "log-log") {
            who <- (result$surv == 0 | result$surv == 1)
            xx <- ifelse(who, 0.1, result$surv)
            temp1 <- exp(-exp(log(-log(xx)) + zval * result$std.err/log(xx)))
            temp1 <- ifelse(who, result$surv + 0 * result$std.err, 
                temp1)
            temp2 <- exp(-exp(log(-log(xx)) - zval * result$std.err/log(xx)))
            temp2 <- ifelse(who, result$surv + 0 * result$std.err, 
                temp2)
            result <- c(result, list(upper = temp1, lower = temp2, 
                conf.type = "log-log", conf.int = conf.int))
        }
    }
    result$call <- Call
    indx <- match("surv", names(result))
    result <- c(result[1:indx], type = attr(y, "type"), result[-(1:indx)])
    if (is.R()) 
        class(result) <- c("survfit.cox", "survfit")
    else oldClass(result) <- "survfit.cox"
    result
}
