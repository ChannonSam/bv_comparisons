# FS-PLS 2019 code:

# can be used to refit a new offset for a set of vars
# return coefficients
refit_model <- function(data,inds){
        if (length(inds) > 0) {  
                x = as.matrix(data$dat[, inds])
                lm <- glm(data$y ~ x, family = getOption("family", "binomial"))
                
                coeff = summary(lm)$coeff[2 : (length(inds) + 1)]
                #print(summary(lm))
        }
        else{
                coeff = c()
        }
        coeff
}

# Aim is to project out all variations corresponding to vars variables
# (this is just indices) from the matrix Dall.
projOut <- function(Dall, vars) {
        P = Dall[, vars, drop = F]
        ncol = dim(Dall)[[2]]
        nrow = dim(Dall)[[1]]
        svd = svd(P)
        U = svd$u
        V = t(svd$v)
        if (length(vars) == 1) {
                D = as.matrix(svd$d)
        }else {
                D = diag(svd$d)
        }
        Vinv = solve(V)
        Dinv = solve(D)
        colU = dim(U)[[2]]
        proj = matrix(0, nrow = dim(U)[[2]], ncol = ncol)
        for (j in 1:ncol) { 	
                for (k in 1:colU) {
                        proj[k, j] = U[, k] %*% Dall[, j]
                }	
        }  
        Vinv %*% Dinv %*% proj
}

orthogAll <- function(data, variables){
        d = orthog(data, variables);
        d[, variables[1]] = data[, variables[1]]
        for (i in 2:length(variables)) {
                d1 = orthog(data, variables[1:(i - 1)])
                d[,variables[i]] = d1[,variables[i]]
                
        }
        d
}

orthog <- function(data, variables){
        
        Dall = as.matrix(data)
        W = matrix(0, nrow = length(variables), ncol = dim(Dall)[[2]])
        W = projOut(Dall, variables)
        R = Dall[, variables, drop = F] %*% W
        Dall - R 
}

# This finds the next best variable to select.
# data is an object with data$y and data$data. Rows are different individuals, and columns diff variables 
# variables is a list of variables already selected
# beta are the weights for those selected variables
# returns updated beta
# can optionally tell it which variable via best_i
find_best_glm <- function(data, variables, beta, var_thresh = getOption("var_thresh", 0.001), project = TRUE, best_i1 = NULL, excl = c()){
        family = getOption("family", "binomial")
        yTr = data$y
        #print(data$d[1:5,1:5])
        #print(yTr)
        #print(paste("besti", best_i1));
        Dall = as.matrix(data$dat)
        W = matrix(0, nrow = length(variables), ncol = dim(Dall)[[2]])
        offset = rep(0, length(yTr))
        if (length(variables >= 1)) {
                if (project) { 
                        print(project)
                        W = projOut(Dall, variables)
                }
                if (getOption("offset", FALSE)) {
                        print("offset");
                        offset = Dall[, variables, drop = F] %*% beta
                }
        }
        print(beta)
        R = Dall[, variables, drop = F] %*% W
        spval = 1.0
        if (!is.null(best_i1) && !is.na(best_i1)) {
                best_i = best_i1;
                x = Dall[, best_i] - R[, best_i]
                #print(yTr)
                #print(paste("besti",best_i, best_i1));
                #print(Dall[,best_i])
                #print(offset)
                #print(family)
                lm <- glm(yTr ~ x, family = family, offset = offset)
                
                sum <- summary(lm)
                spval <- sum$coefficients[8] 
        }else{
                
                todo = 1:dim(Dall)[[2]];
                overl = which(todo %in% excl)
                if (length(overl) > 0) todo = todo[-overl];
                best_i = 0
                #print(todo);
                for (i in todo) {
                        
                        x = Dall[, i] - R[, i]
                        vari = var(x, na.rm = TRUE);
                        #print(paste(vari, var_thresh))
                        if (!is.na(vari) && vari > var_thresh) {
                                
                                lm <- glm(yTr ~ x, family = family, offset = offset)
                                sum <- summary(lm)
                                pval <- sum$coefficients[8] ## t-test p-value
                                #print(spval)
                                if (pval >= 0 && pval < spval ) {
                                        spval = pval
                                        best_i = i
                                        
                                }
                                
                        }
                }
        }
        if (best_i == 0) print(paste('problem, best_i should not be 0.', best_i1, spval, dim(Dall)[[2]])); 
        x = Dall[, best_i] - R[, best_i]
        if (FALSE) {
                print('no shrinking');
                model <- glm(yTr ~ x, family = family, offset = offset) 
                coeffNew = summary(model)$coeff[2]
        }else{
                #coeff = beta - coeffNew*W[,best_i]
                lambda = getOption("lambda");
                if (!is.null(getOption("lambda2"))) lambda = getOption("lambda2");
                coeffNew = fitModelShrink(yTr, as.matrix(x), offset, lambda = lambda)
                
                #  const_term = const_term+coeffNew[1,]
                #print(coeffNew)
                coeffNew = coeffNew[-1,, drop = T]
        }
        if (coeffNew == 0) print("is zero!!")
        coeff = beta - coeffNew*W[, best_i]
        print(paste(getOption("lambda") , coeffNew))
        #print(coeff)
        #if(!pls) coeff=coeffNew
        #else coeff=rbind(beta,coeffNew[1,])
        
        list(pval = spval, variable = c(variables, best_i), coeffNew = coeffNew, coeff = c(coeff, coeffNew), var = var(Dall[, best_i], na.rm = T))
}

# this just returns the auc
# vars are indices referring to data$data
calcAUC <- function(coeff, vars, data) {
        auc = NaN
        family = getOption("family", "binomial")
        if (dim(data$dat)[[1]] > 1) {
                testD = data$dat
                yTs = data$y
                xnew = testD[, vars[length(vars)]];
                xM = as.matrix(testD[, vars]) %*% coeff
                x = xM[, 1]
                nonNA = which(!is.na(x))
                model <- glm(yTs[nonNA] ~ x[nonNA], family = family) 
                model1 <- glm(yTs[nonNA] ~ xnew[nonNA], family = family) 
                predictions = predict(model, x = x[nonNA], type = "response")
                #print(predictions)
                #print(yTs)
                #pred <- prediction(predictions, yTs[nonNA])
                #perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
                auc <- ci(yTs[nonNA], predictions)
                
                #auc<-performance(pred, measure = "auc")@y.values[[1]] 
                #plot(perf, col=rainbow(10))
        }
        #print(summary(model)$coeff)
        #print(summary(model)$coeff[8])
        c(auc, summary(model1)$coeff[8])
}

# following functions not currently in use
updown<-function(beta){
        r = rep(0,length(beta))
        for( i in 1:length(beta)){
                if(beta[i]>0) r[i] = 1 else r[i]= -1
        }
        r
}

updownModel<-function(model){
        model$beta = updown(model$beta);
        for(i in 1:length(model$betaG)){
                model$betaG[[i]] = updown(model$betaG[[i]]);
                model$betaGlobal[[i]] = updown(model$betaGlobal[[i]]);
                model$coeffOrth[[i]] = updown(model$coeffOrth[[i]]);
        }
        model
}

fitModelShrink <- function(yTr, x, offset, family = getOption("family", "binomial"),  useBF = getOption("useBF", FALSE), lambda = 0){
        #  debugprint('shrinking0')
        #  if((family=="binomial"||family=="multinomial")) alpha = 0
        #  else alpha = 1
        #print(x)
        alpha = getOption("alpha", 0)
        #alpha = .Machine$double.eps;
        ones = rep(1, length(yTr))
        x2 = cbind(x, ones)
        y2 = yTr
        #if(family=="binomial"&&getOption("method","fspls")=="fs"&&abs(offset[1])>1) offset=rep(0,length(offset))
        if (is.null(lambda)) {
                print("CROSS_VALIDATION");
                aa = cv.glmnet(x2, y2, family = family, alpha = alpha, offset = offset, lambda = getOption("testlambda", NULL))
                print(paste(min(aa$lambda), max(aa$lambda)))
                print(paste(aa$cvm[1], min(aa$cvm), aa$cvm[length(aa$cvm)]));
                
                #print(aa$cvlo)
                #	  print(aa$cvup)
                aamin0 = which(aa$cvm == min(aa$cvm))
                tol = getOption("tolerance", 1.0);
                aamin = which(aa$cvm <= min(aa$cvm) * tol)
                #print(aa$lambda[aamin])
                #        indsa = which(aa$cvm>=aa$cvlo[aamin] & aa$cvm <=aa$cvup[aamin])
                #       print(indsa)
                #	print(aamin)
                lambda =  max(aa$lambda[aamin])
                lambda0 =  max(aa$lambda[aamin0])
                options("lambda2" = lambda)
                print(paste("lambda", lambda, lambda0))
        }
        if (family != "multinomial") {
                #debugprint(paste('shrinking1', lambda, alpha, family))
                #	debugprint(x2)
                #	debugprint(y2)
                #	debugprint(offset)
                #print(y2)
                ridge = glmnet(x2, y2, family = family, lambda = lambda, alpha = alpha, offset = offset, intercept = F)
                #debugprint('shrinking2')
                rbeta <- coef(ridge, s = "lambda.1se")
                coeffNew = as.matrix(rbeta[-length(rbeta)])
        }else{
                ridge = glmnet(x2,y2,family = family,lambda = lambda, alpha = alpha, offset = offset, intercept = F)
                rbeta <- coef(ridge, s = "lambda.1se")
                rbeta = matrix(unlist(lapply(rbeta, as.matrix)), nrow = dim(x)[2] + 2)
                rbeta = rbeta - rbeta[,1]
                coeffNew = as.matrix(rbeta[ -(dim(x)[2] + 2), -1])
        }
        
        coeffNew
}



# train is an object with train$data and train$y  where y is outcome and data is data matrix (cols are variables)
# var_thresh ignores variables with var < var_thresh
# project - do we project out variables we already selected for search
# refit - do we keep the coeff trained for each variable separately, or refit as we go
# you can already specificy variables and beta selected
trainModel <- function(train, pv_thresh = 0.01, max=min(15, dim(train$data)[[2]]), project = TRUE, refit = FALSE, 
                       variables = c(), beta = c(), info = c(), auc = c(), auc_l = c(), auc_u = c(), sel_vars = c(), excl = c()) {
        spval = 0.0
        
        options("lambda2" = NULL);
        coeffOrth = list()
        coeffo = c();
        toremove = which(sel_vars %in% excl)
        if (length(toremove) > 0) sel_vars = sel_vars[-toremove];
        pv_thresh1 = 1.0;
        contains = c()
        global_index = length(variables) + 1
        variablesGlobal = list() # list of the subsets of variables selected
        betaGlobal = list() ## list of the subsets of betas selected at each stage 
        if (length(variables) > 0) {
                for (k in 1:length(variables)) {
                        variablesGlobal[[k]] = variables[1:k]
                        betaGlobal[[k]] = beta[1:k]
                }
        }
        while (global_index <= max && (spval <= pv_thresh1 ) && length(contains) <= 1) {
                best_i = NULL;
                if (global_index <= length(sel_vars)) best_i = sel_vars[global_index]
                if (!is.null(best_i)) pv_thresh1 = 1.0 else pv_thresh1 = pv_thresh
                
                print(global_index)
                model_fit = find_best_glm(train, variables, beta, project = project, best_i = best_i, excl = excl)
                print(global_index)
                coeffo = c(coeffo, model_fit$coeffNew);
                varia = model_fit$variable
                lastVar = varia[length(varia)]
                contains = which(varia == lastVar)
                spval = model_fit$pval # pval of selected model
                if (length(contains) <= 1 && spval <= pv_thresh1) {
                        variables = varia
                        if (!refit) {
                                beta = model_fit$coeff
                        } else {
                                beta = refit_model(train, variables)
                        }
                        
                        coeffOrth[[global_index]] = coeffo
                        betaGlobal[[global_index]] = beta
                        variablesGlobal[[global_index]] = variables
                        info[global_index] = spval; #paste(model_fit$var, spval,sep=" ")
                        auc1 = calcAUC(beta, variables, train)
                        auc = c(auc, auc1[2])
                        auc_l = c(auc_l,auc1[1])
                        auc_u = c(auc_u, auc1[3])
                        toprint = c(variables[[global_index]], auc1[2],auc1[4],info[[global_index]]);
                        names(toprint) = c("variable","auc","pv", "spval");
                        print(toprint);
                        #		print(paste,sep=" "))
                        global_index = global_index + 1
                }
                if (length(contains) > 1) print("STOPPED DUE TO REPEAT")
        }
        list(variablesGlobal = variablesGlobal, coeffOrth = coeffOrth, betaGlobal = betaGlobal, info = info, auc = auc, auc_l = auc_l, auc_u = auc_u, variables = variables, 
             beta = beta)
}

# Used to get a data frame with all products of columns to capture interactions
getCrossTerms <- function(train1, vars){
        data1 = centralise(train1$data)
        n = dim(data1)[[2]]
        m = length(vars)
        len = n*m
        dataRes = matrix(0,nrow = length(train1$y), ncol = len )
        k = 1
        nme = matrix(0, nrow = 2, ncol = len)
        for (i in 1:n) {
                for (j in 1:m) {
                        dataRes[,k] = data1[,i]*data1[,vars[j]]
                        nme[,k] = c(i,vars[j])
                        k = k + 1
                }
        }
        nme1 = rbind(1:n, rep(NA,n))
        list(data = cbind(train1$data,dataRes),y = train1$y, nme = cbind(nme1,nme))
}

getCrossTermsAll <- function(train1) {
        data1 = centralise(train1$data)
        n = dim(data1)[[2]]
        len = (n*(n - 1))/2
        dataRes = matrix(0, nrow = length(train1$y), ncol = len )
        k = 1
        nme = matrix(0, nrow = 2, ncol = len)
        for (i in 1:(n - 1)) {
                for (j in (i + 1):n) {
                        dataRes[,k] = data1[,i]*data1[,j]
                        nme[,k] = c(i,j)
                        k = k + 1
                }
        }
        nme1 = rbind(1:n, rep(NA,n))
        list(data = cbind(train1$data,dataRes), y = train1$y, nme = cbind(nme1,nme))
}

centralise <- function(data) {
        dimd = dim(data)
        print(dimd)
        mean = apply(data,2,mean)
        sd = apply(data,2,sd)
        for (i in 1:dimd[[2]]) {
                data[,i] = (data[,i] - rep(mean[i],dimd[[1]]))/sd[i]
        }
        data
}

# test is as for trainModel
# results is results from train model
testModels <- function(all, results, orth = FALSE){
        betaG = results$betaGlobal
        if (orth) betaG = results$coeffOrth
        variables = results$variables
        #info = results$info
        print(paste("varid","aucTest", "aucTrain", "variance", "spval"))
        print(variables)
        res  = matrix(NA, nrow = length(betaG), ncol = 12)
        
        for (k in 1:length(betaG)) {
                aucTest = calcAUC(betaG[[k]], variables[1:k], all$test)
                aucTrain = calcAUC(betaG[[k]], variables[1:k], all$train)
                aucComb = calcAUC(betaG[[k]], variables[1:k], all$comb)
                res[k,] = c(aucTrain[1:4], aucTest[1:4], aucComb[1:4])
                #print(paste(variables[k],aucTest[1],aucTest[2],sep="  "))
        }
        res1 = data.frame(res)
        names(res1) = c("AUC train lower","AUC train", "AUC train upper", "pv train", "AUC test lower","AUC test", "AUC test upper",  "pv test","AUC comb lower","AUC comb", "AUC comb upper","pv comb");
        res1
}

