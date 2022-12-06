library(mgcv)

train_gam <-
function(X,y,pars = list(numBasisFcts = 10))
{  
    X_matrix <- as.matrix(X)
    p <- dim(X_matrix)    

    if(!("numBasisFcts" %in% names(pars) ))
    { 
        pars$numBasisFcts = 10
    }   
    if(p[1]/p[2] < 3*pars$numBasisFcts)
    {
        pars$numBasisFcts <- ceiling(p[1]/(3*p[2]))
        cat("changed number of basis functions to    ", pars$numBasisFcts, "    in order to have enough samples per basis function\n")
    }
    coln <- rep("null",p[2]+1)
    Ks <- rep(pars$numBasisFcts, p[2]) # DYNAMIC K VALUE
    for (i in 1:p[2])
    {
        nunique <- length(unique(X_matrix[,i]))
        if (nunique < 100)
        {
            X_matrix[,i] <- as.factor(X_matrix[,i])
        }
        if (Ks[i] > nunique)
        {
            Ks[i] <- nunique - 1
        }
    }
    y_matrix <- as.matrix(y)
    p_y <- dim(y_matrix) 
    for (i in 1:p_y[2])
    {
        nunique <- length(unique(y_matrix[,i]))
        if (nunique < 100)
        {
            y_matrix[,i] <- as.factor(y_matrix[,i])
        }
    }   
    
    dat <- data.frame(y_matrix, X_matrix)
    for(i in 1:(p[2]+1))
    {        
        coln[i] <- paste("var",i,sep="")
    }
    colnames(dat) <- coln
    labs<-"var1 ~ "
    if(p[2] > 1)
    {
        for(i in 2:p[2])
        {
            labs<-paste(labs,"s(var",i,",k = ",Ks[i - 1],") + ",sep="")
        }
    }
    labs<-paste(labs,"s(var",p[2]+1,",k = ",Ks[p[2]],")",sep="")
    print(labs)
    mod_gam <- FALSE
    # THE PROBLEM IS HERE GAM 
    #  "Error in smooth.construct.tp.smooth.spec(object, data, knots) :
    # > A term has fewer unique covariate combinations than specified maximum
    # > degrees of freedom"
    # TRIED LOWERING THE K value but didnt work.
    try(mod_gam <- gam(formula=formula(labs), data=dat),silent = TRUE)
    if(typeof(mod_gam) == "logical")
    {
        cat("There was some error with gam. The smoothing parameter is set to zero.\n")
        labs<-"var1 ~ "
        if(p[2] > 1)
        {
            for(i in 2:p[2])
            {
                labs<-paste(labs,"s(var",i,",k = ",Ks[i - 1],",sp=0) + ",sep="")
            }
        }        
        labs<-paste(labs,"s(var",p[2]+1,",k = ",Ks[p[2]],",sp=0)",sep="")
        print(labs)
        #mod_gam <- gam(formula=formula(labs), data=dat)
        try(mod_gam <- gam(formula=formula(labs), data=dat), silent = TRUE) 
        if(typeof(mod_gam) == "logical") 
        {
            cat("New error with gam. Returning empty list.\n")
            return(list())
        }
        
    }
    result <- list()
    result$Yfit <- as.matrix(mod_gam$fitted.values)
    result$residuals <- as.matrix(mod_gam$residuals)
    result$model <- mod_gam 
    result$df <- mod_gam$df.residual     
    result$edf <- mod_gam$edf     
    result$edf1 <- mod_gam$edf1     
    
    # for degree of freedom see mod_gam$df.residual
    # for aic see mod_gam$aic
    return(result)
}


selGam <-
function(X,pars = list(cutOffPVal = 0.001, numBasisFcts = 10),output = TRUE,k)
{
    result <- list()
    p <- dim(as.matrix(X))
    if(p[2] > 1)
    {
        selVec <- rep(FALSE, p[2])
        mod_gam <- train_gam(X[,-k],as.matrix(X[,k]),pars)
        pValVec <- 0.0
        if (length(mod_gam) == 0) # This is to handle the error with GAM. Returning p-value vector of zeros. 
        {
            pValVec <- rep(0.0, length(selVec[-k]))
        }
        else{
          pValVec <- summary.gam(mod_gam$model)$s.pv
        }
        if(output)
        {
            cat("vector of p-values:", pValVec, "\n")
        }
        if(length(pValVec) != length(selVec[-k]))
        {
            show("This should never happen (function selGam).")
        }
        selVec[-k] <- (pValVec < pars$cutOffPVal)
    } else
    {
        selVec <- list()
    }
    return(selVec)
}

pruning <-
function(X, G, output = FALSE, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = 0.001, numBasisFcts = 10)) 
{
    p <- dim(G)[1]
    finalG <- matrix(0,p,p)
    for(i in 1:p)
    {
        parents <- which(G[,i]==1)
        lenpa <- length(parents)

        if(output)
        {
            cat("pruning variable:", i, "\n")
            cat("considered parents:", parents, "\n")
        }
        
        if(lenpa>0)
        {
            Xtmp <- cbind(X[,parents],X[,i])
            selectedPar <- pruneMethod(Xtmp, k = lenpa + 1, pars = pruneMethodPars, output = output)
            finalParents <- parents[selectedPar]
            finalG[finalParents,i] <- 1
        }
    }
    
    return(finalG)
}

print('SCRIPT R')
dataset <- read.csv(file='{PATH_DATA}', header=FALSE, sep=",")
dag <- read.csv(file='{PATH_DAG}', header=FALSE, sep=",")
set.seed(42)
pruned_dag <- pruning(dataset, dag, pruneMethod = selGam, pruneMethodPars = list(cutOffPVal = '{CUTOFF}', numBasisFcts = 10), output='{VERBOSE}')

write.csv(as.matrix(pruned_dag), row.names = FALSE, file = '{PATH_RESULTS}')
