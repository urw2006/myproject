# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
#  Collection of functions       #
#  Author: Udo Wittmann          #
#  Version: April 2015           #
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
library(ResourceSelection)      # blablabla  hjhjh
library(pROC)      

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                    FUNCTIONS                                        #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
# ----  infectionIntensity ---- #
# converts absolute epg counts into factors uninfected, light, moderate, heavy

infectionAnaemia  <-  function(hb, sex, age) { 
  
  # arguments: hb = haemoglobin level in g/dl
  #            sex = "M" or "F"
  #            age = numerical, age of individual
  # Returns:  the intensity of infection as
  #           'none': 'uninfected' , 'light', 'moderate', 'heavy' 
  #           as given in the document  WHO-2011 - helminth control in School age children.pdf p.ix
  #anaemia
  # Clinical condition reflecting an inadequate number of red
  # blood cells or an insuffi cient amount of haemoglobin (Hb)
  # in the blood. 
  # Children aged 5-11 years are considered to
  # be anaemic when the Hb concentration is below 115 g/litre
  # and severely anaemic when it is below 80 g/litre. Children
  # aged 12-14 years are considered to be anaemic when the Hb
  # concentration is below 120 g/litre and severely anaemic when
  # it is below 80 g/litre.
  # From Anna:
  # Anaemia was defined according to WHO guidelines (51) as a Hb level <11.5 g/dl for children 
  # 5-11 years old, <12.0 g/dl for children 12-14 years old, <12.0 g/dl for 
  # females ???15 years old and <13.0 g/dl for males ???15 years old.
  
 
  if(!(  (length(hb) == length(age))  &  (length(age) == length(sex))   )) stop("Arguments must have the same lengths.")
  
  N <- length(hb)
  out <- numeric(N)
  
  for(i in 1:N) {
    if( is.na(hb[i]))  out[i] <- NA
    else { out[i] <- 0     
    
    if(age[i] <= 11 & hb[i] < 11.5) out[i] <- 1
    if((age[i] >= 12 & age[i] <= 14) & hb[i] < 12) out[i] <- 1
    if((age[i] >= 15 & sex[i] == "F") & hb[i] < 12) out[i] <- 1
    if((age[i] >= 15 & sex[i] == "M") & hb[i] < 13) out[i] <- 1
     
          }
  }
 
  return(out)
} # end infectionAnaemia

#%%%%%%%%%%%%%%%%%%%%%%%

infectionIntensity  <-  function(ec, type) { 
  
  # arguments: ec = egg count per gram of faeces or 10ml of urine respectively. 
  #            type = c('asc' , 'trich' , 'hk' , 'sm' , 'sh' )   
  # Returns:  the intensity of infection as
  #           'none': 'uninfected' , 'light', 'moderate', 'heavy' 
  #           as given in the document  WHO-2011 - helminth control in School age children.pdf p.54
  #  sub function for one data observation: 
  ii <- function(ec, type) 
  { 
    if(is.na(ec)) out <- NA
    else {
      
      if(ec == 0) {out <- 'uninfected'}
      else {  
        out <- 
          switch(type , 
                 'asc' = {  ifelse(ec > 0 & ec < 5000, "light", 
                                   ifelse(ec >= 5000 & ec < 50000, "moderate", 
                                          ifelse(ec >= 50000, "heavy", NA)))
                   
                   
                 },
                 'trich' = {  ifelse(ec > 0 & ec < 1000, "light", 
                                     ifelse(ec >= 1000 & ec < 10000, "moderate", 
                                            ifelse(ec >= 10000, "heavy", NA)))
                   
                   
                 },
                 'hk' = {  ifelse(ec > 0 & ec < 2000, "light", 
                                  ifelse(ec >= 2000 & ec < 4000, "moderate", 
                                         ifelse(ec >= 4000, "heavy", NA)))
                   
                   
                 },
                 
                 'sm' = {  ifelse(ec > 0 & ec < 100, "light", 
                                  ifelse(ec >= 100 & ec < 400, "moderate", 
                                         ifelse(ec >= 400, "heavy", NA)))
                   
                   
                 },
                 'sh' = {  ifelse(ec > 0 & ec < 51, "light", 
                                  ifelse(ec >= 51, "heavy", NA))
                   
                   
                 }
                 
          ) # end switch
      }
    } #end of outer else
    return(out)
  } # end of subfunction ii
  
  out <- apply(   as.matrix(ec)  ,  1   ,  function(x) ii(x, type = type) )
  out <- factor(out, ordered = TRUE , levels = c("uninfected" , "light", "moderate", "heavy"))
  return(out)
} # end infectionIntensity

# %%%%%%%%%%%%
#### ggQQ ####
# %%%%%%%%%%%%

ggQQ <- function(vec, title) { 
# Plots a qq-normal-plot  
# Argument: vector of numbers
# Returns: qqnorm plot as ggplot-object (qqnorm doesn't do this)

  # following four lines from base R's qqline()
  y <- quantile(vec[!is.na(vec)], c(0.25, 0.75))
  x <- qnorm(c(0.25, 0.75))
  slope <- diff(y)/diff(x)
  int <- y[1L] - slope * x[1L]
  
  d <- data.frame(resids = vec)
  
  out <- ggplot(d, aes(sample = resids)) + stat_qq() + geom_abline(slope = slope, intercept = int) + ggtitle(title)
  
  return(out)
}

## Plotting histogram of random effects with ggPlot: 

ranef.histogram <- function(data.vector, binwidth, title)
{     
  d.plot <- data.frame(re = data.vector)
  names(d.plot) <- c("r")
  
  p  <- ggplot(d.plot, mapping=aes(x=r)) +
    geom_histogram( binwidth=binwidth, alpha=0.8, position="identity") +
    labs(x="Residuals", y="Counts", title=title) + 
    theme(title = element_text(size=8))
  
  return(p) 
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####            extract_legend             #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
extract_legend <- function(a.gplot){
  # returns the extracted legend from a ggplot plot a.plot
  # from http://stackoverflow.com/questions/13649473/add-a-common-legend-for-combined-ggplots  
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####                plot.ranres             #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot.ranef <- function(m, binwidths = c(0.4, 0.003)){
  # plots qq-norm plot and histogram of all random effects of a model
  # Argument: m is object of class glmmadmb
  #           binwidths is a vector of class numeric. It gives the binwidths for the histograms for each group. 
  #           length of binwidths must be number of groups, otherwise abort
  
  re <- ranef(m)
  
  nGroups <- length(re)
  namesGroups <- names(re)
  
  if (length(binwidths) != nGroups) stop("Lengths of binwidths must equal number of grouping variables.")
  
  plots <- list()
  nPlots <- 0
  
  for(i in 1:nGroups) { 
    
    nCol <- dim(re[[i]])[2]
    for(j in 1:nCol) {  nPlots <- nPlots + 1
                        plotTitle <- paste(namesGroups[i] , "-", colnames(re[[i]])[j])
                        plots[[nPlots]] <- ggQQ(re[[i]][, j], title = plotTitle)
                        
    }  # end inner for
  }  # end outer for  
  
  histos <- list()
  nPlots <- 0
  
  for(i in 1:nGroups) { 
    
    nCol <- dim(re[[i]])[2]
    for(j in 1:nCol) {  nPlots <- nPlots + 1
                        plotTitle <- paste(namesGroups[i] , "-", colnames(re[[i]])[j])
                        histos[[nPlots]] <- ranef.histogram(re[[i]][, j], binwidth = binwidths[i],  title = plotTitle)
                        
    }  # end inner for
  }  # end outer for  
  
  
  
  return( list(qqPlots = plots, histograms = histos) )
} # end plot.ranef



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#====                          HURDLE MODEL METHODS                        ======#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

humm <- function(x, ...) UseMethod("humm")

overdisp <- function(x, ...) UseMethod("overdisp") 

diagnostics <- function(x) UseMethod("diagnostics")

writeModel <- function(x) UseMethod("writeModel")



as.humm <- function(p, c) {
# takes two glmmADMB models c = count-model (zero-truncated model), p = prevalence-model (logit model)
# returns a humm - model object ("hurdle model")
  if (p$family %in% c("binom")  &  c$family %in% c("truncnbinom", "truncnbinom1", "truncpoiss" )) out <- list(prev = p, count = c, msg = NA)                                    
  else stop("The first model must have binom, the second a truncated family as distribution family,")
  
  class(out) <- "humm"
  
  return(out)
}


 
humm.default  <- function(form, data, family = c("binomial", "truncnbinom1") , iter = 40000, maxph = 5 , verbose = FALSE )
# fits a mixed hurdle model by fitting a count-model and a prevalence-model separately with glmmADMB. 
# Both models are fitted with the same model specification defined by formula. 
# Arguments:  form   formula as for a glmmADMB model
#             data   data set to be used, for the zero-truncated model the subset without zero-responses is used
#             family a character vector consisting of the two family specifications for the prevalence and count -model. 
#             iter   number of iterations, this value is passed to admb via  extra.args="-ndi <iter-value>" 
#             maxph  passed on with  admb.opts = admbControl( maxph = <maxph-value> )  
#
# Return-Value: returns a humm-object  
# to do: two separate output directories, organise warnings,   
{
  
  myWarnings <- NULL     # list to collect warnings
  myErrors <- NULL    # ditto
  prevalence <- NULL
  
  mHandler <- function(m)
  {
    message(m)
  }  
  
  wHandler <- function(w)  # handler for warnings  
  { 
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  } #end wHandler
  
  eHandler <- function(e)  # handler for errors
  {
    myErrors <<- conditionMessage(e)
    # invokeRestart("abort")
    NULL
  }  
  
   prevalence <- withCallingHandlers( tryCatch(
                  glmmadmb( form, 
                               data = data,
                               family = family[1], 
                               extra.args="-ndi 1000" ,
                               admb.opts = admbControl(noinit = FALSE , maxph = maxph ),
                               # save.dir = "S:/SCI - post 3 June 2011/Current programmes/Private Donor - Burundi/Impact analysis/Models and plots in R/Development Intensity/Comparison Baseline - last year/Pilot/glmmADMB_ouput" , 
                               verbose = verbose     
                               ),
                  message = mHandler,
                  warning = wHandler,
                  error = eHandler) #end tryCatch
                  ) # end withCallingHandler
  
  count <- withCallingHandlers( tryCatch(
    glmmadmb( form, 
              data = subset(data, eval(form[[2]]) > 0 ),   
              family = family[2], 
              extra.args="-ndi 1000" ,
              admb.opts = admbControl(noinit = FALSE , maxph = maxph ),
              # save.dir = "S:/SCI - post 3 June 2011/Current programmes/Private Donor - Burundi/Impact analysis/Models and plots in R/Development Intensity/Comparison Baseline - last year/Pilot/glmmADMB_ouput" , 
              verbose = verbose     
    ),
    warning = wHandler,
    error = eHandler) #end tryCatch
  ) # end withCallingHandler
  

  return(list(value = prevalence , warnings = myWarnings, error = myErrors ))
  
} # end humm.default



#%%%%%%%%%%%%%%%%%%%%%
####  print.humm  ####
#%%%%%%%%%%%%%%%%%%%%%

print.humm <- function(x, ...){
  
  ps <- summary(x[[1]])      # prevalence summary object
  cs <- summary(x[[2]])      # count model summary 
  
  
  cat("Calls:\n Prevalence\n")
  print(x[[1]]$call)
  cat(" Count\n")
  print(x[[2]]$call)
  
  cat("Coefficients:\n Prevalence\n")
  print(ps$coefficients)
  cat("Coefficients:\n Count\n")
  print(cs$coefficients)
  
} #  end print.humm

#%%%%%%%%%%%%%%%%%%%%%
#### summary.humm ####
#%%%%%%%%%%%%%%%%%%%%%

summary.humm <- function(model){
# Argument: model of class humm
  
  ps <- summary(model[[1]])      # prevalence summary object
  cs <- summary(model[[2]])      # count model summary 
  
  pcoefficients <- ps$coefficients
  ccoefficients <- cs$coefficients
  
  pS <- ps$S
  cS <- cs$S
  
  alpha <- cs$alpha
  sd_alpha <- cs$sd_alpha
  
  out <- list(pCall = model[[1]]$call, cCall = model[[2]]$call, pCoefficients = pcoefficients, cCoefficients = ccoefficients, pS = pS, cS = cS, alpha = alpha, sd_alpha = sd_alpha)
  
  class(out) <- "summary.humm"
  return(out)
  
} # end summary.humm


#%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### print.summary.humm ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

print.summary.humm <- function(x, ...) {
  
  cat("Calls:\n Prevalence\n")  
  print(x$pCall)
  cat(" Count\n")
  print(x$cCall)
  cat("\n")
  
  cat("Fixed Effects:\n Prevalence\n")
  printCoefmat(x$pCoefficients)
  cat("\n Count\n")
  printCoefmat(x$cCoefficients)
  cat("\n")
  cat("Negative binomial dispersion parameter: ", x$alpha, " (std.err.: ", x$sd_alpha, ")\n\n", sep="" )
  
  pGroupNames <- attributes(x$pS)$names
  cGroupNames <- attributes(x$cS)$names
  pNoGroups <- length(pGroupNames) 
  cNoGroups <- length(cGroupNames)
  
  
  cat("Random Effects Covariance-Matrices:\n Prevalence\n")
  for(i in 1:pNoGroups) {
    cat("Group=", pGroupNames[i], "\n", sep="")
    print(x$pS[[i]])
  }
   
  cat("\n Count\n")
  for(i in 1:cNoGroups) {
    cat("Group=", cGroupNames[i], "\n", sep="")
    print(x$cS[[i]])
  }
  cat("\n")
  
  
} # print.summary.humm
  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%
###### writeModel.humm ######
#%%%%%%%%%%%%%%%%%%%%%%%%%%%  

writeModel.humm <- function(m, wb, sheet, xy ){
# writes the summary of an object of class humm to an excel file
# Argument:
  
#  wb <- workbook
  
  sm <- summary(m)
  
  openxlsx::writeData(wb = wb, sheet = sheet, xy = xy,
                      x = "Fixed Effects:\n")
  
  openxlsx::writeData(wb = wb, sheet = sheet, xy = xy + c(1,7),
                      x = "Prevalence:")
  
  openxlsx::writeData(wb = wb, sheet = sheet, xy = xy + c(1,2),
                      x = round( sm$pCoefficients, 3) )
  
  openxlsx::writeData(wb = wb, sheet = sheet, xy = xy + c(1,4 + dim(sm$pCoefficients)[1] ),  x = "Count:")
  
  openxlsx::writeData(wb = wb, sheet = sheet, xy = xy + c(1,5 + dim(sm$pCoefficients)[1]),
                      x = round(sm$cCoefficients, 3))

  XY <- xy + c(1,4 + dim(sm$pCoefficients)[1] +  dim(sm$cCoefficients)[1])
# Random Effects
  openxlsx::writeData(wb = wb, sheet = sheet, xy = XY,
                    x = "Random Effects:\n")

  openxlsx::writeData(wb = wb, sheet = sheet, xy = XY + c(1,1),
                    x = "Prevalence:")

  openxlsx::writeData(wb = wb, sheet = sheet, xy = XY + c(1,2),
                    x = sm$pS  )
                    
  openxlsx::writeData(wb = wb, sheet = sheet, xy =  c(1,2 + 5 ),
                    x = "Count:")
                    
  openxlsx::writeData(wb = wb, sheet = sheet, xy =   c(1,3 +10),
                    x = sm$cS )
                      
                      
                      
} # end writeData.humm  

#%%%%%%%%%%%%%%%%%%%%%%%%%%%
######  logLik.humm   ######
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

logLik.humm <- function(model) {
# Argument: model of class humm
# Returns:  log likelihood  of this model
  out <- double()
  
  out[1] <- model[[1]]$loglik + model[[2]]$loglik
  out[2] <- model[[1]]$np + model[[2]]$np
  class(out) <- "logLik.humm" 
  
  attr(out, "pLogLik") <- model[[1]]$loglik
  attr(out, "cLogLik") <- model[[2]]$loglik
  attr(out, "pNp") <- model[[1]]$np
  attr(out, "cNp") <- model[[2]]$np
  
  out
}

print.logLik.humm <- function(x, ...){
  
  cat("logLik Hurdle Model:\n")
  cat(x[1], " (df=", x[2], ")\n" , sep="")
  
  cat("logLik Prevalence: ", attr(x, "pLogLik"), " (df=", attr(x, "pNp"), ")\n" , sep="")
  
  cat("logLik Count: ", attr(x, "cLogLik"), " (df=", attr(x, "cNp"), ")\n" , sep="")
} # end print.logLik.humm



#%%%%%%%%%%%%%%%%%%%%%%%%%%%
######     AIC.humm   ######
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

AIC.humm <- function(model) {
# Argument: model of class humm
# returns: AIC of model
  
  AIC(model[[1]]) + AIC(model[[2]])
  
} # end AIC.humm

#%%%%%%%%%%%%%%%%%%%%%%%%%%%
######  anova.humm    ######
#%%%%%%%%%%%%%%%%%%%%%%%%%%%

anova.humm <- function(x, ...) {
# likelihood ratio test for humm models
# Arguments: one or more objects of class humm
# Returns: Results from LRT test
  
#  match.call(expand.dots = TRUE) 
  list(...) -> args 
  
  if (length(args) > 1) stop("Not more than 2 models can be compared.") 
  
  model1 <- x
  model2 <- list(...)

 
  model1$np <- x[[1]]$np + x[[2]]$np
  model2$np <- args[[1]][[1]]$np + args[[1]][[2]]$np 
  model1$loglik <- x[[1]]$loglik + x[[2]]$loglik
  model2$loglik <- args[[1]][[1]]$loglik + args[[1]][[2]]$loglik
  
  D <- 2*abs(model1$loglik - model2$loglik)
  d.np <- abs(model1$np - model2$np)
  
  p <- pchisq(q = D, df = d.np, lower.tail = FALSE) 
  
  pD <- 2*(abs(x[[1]]$loglik - model2[[1]][[1]]$loglik))
  pd.np <-  abs(x[[1]]$np - model2[[1]][[1]]$np)
  
  pp <- pchisq(q = pD , df = pd.np, lower.tail = FALSE )
  
  cD <- 2*(abs(x[[2]]$loglik - model2[[1]][[2]]$loglik))
  cd.np <- abs(x[[2]]$np - model2[[1]][[2]]$np)
  cp <- pchisq(q = cD, df = cd.np , lower.tail = FALSE )
  
  return(list(chiSq = D, df = d.np, p.value = p, prev.p = pp, prev.df = pd.np , count.p = cp, count.df = cd.np))
  
}

# How to work with unspecified list of arguments:


f <- function(...) {
  list(...) -> args
  args[[1]] + args[[2]] + args[[3]]
}





#diagLogit <- function(m){
# Diagnostics for logit model of class glmmadmb
# Argument: object m of class glmmadmb with family = binom
  
#if ( (class(m) != "glmmadmb") | (m$family != "binom") ) stop("Object is not of correct type (glmmadmb with family = binom)") 
#}




#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####              hosmerLem                 #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hosmerLem <- function(m, g = 10){
# Hosmer - Lemeshow goodness of fit test for glmmadmb object with family binom
# Arguments: m is a model of class glmmadmb with family "binom" (i.e. a logistic model)
#            g is the number of groups the responses are split into for the Hosmer-Lemeshow test. 10 seems to be standard
# Returns:  an object of class htest (which is basically a list). The elements of this list are 
#           "statistic" "parameter" "p.value"   "method"    "data.name" "observed"  "expected" 
  
  if (!(class(m) %in% c("glmmadmb", "glmerMod"))) { stop("The argument is not of the correct type (class = glmmadmb or glmerMod, family = binom)" ) }
 # else {if( (m$family != "binom") ) stop("The argument is not of the correct type (class = glmmadmb, family = binom)") }

  if (class(m) == "glmmadmb") mat <- cbind( m$frame[,1] , prob = fitted(m)) 
  if (class(m) == "glmerMod") mat <- cbind( attr(m, "frame")[,1], prob = fitted(m))
  
  return( hoslem.test(mat[,1], mat[, "prob"], g = g) )
  
} # end hosmerLem

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####       Test for Overdispersion          #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#plogisCorrected <- function( beta , sigma2 , type = "mcculloch" ) {
  # beta = fixed effects as returned by glmer (as returned buy fixef((model)
  # sigma02 = variance of random effect as returned by VarCorr(model)
  # type: either "mcculloch" or "diggle" for the respective formulae
  
#  if(type == "mcculloch") { 
#    c2 <- ( (16*sqrt(3))/(15*pi)  )^2 
#    out <- plogis( beta/sqrt(1+c2*sigma2)  )   
#  } 
#  else { 
#    c <- ( 1 + 2*exp( -0.5*sigma2 )) / 6
#    out <- plogis( beta - 0.5*sigma2 * tanh( beta* c  )
#    )
#  }
#} # end plogisCorrected     

overdisp.glmmadmb <- function(model, dfJ=TRUE) { 
  ##  model: model object of class glmmadmb:  
  ##  dfJ: Boolean
  ##        TRUE: J < n  # see Hosmer, Lemeshow: "Applied Logistic Reg", p138
  ##        FALSE: J = n 
  ##  this function was basically taken from http://glmm.wikidot.com/faq
  vpars <- function(m) {
  # returns sum(1+2+...+ nrow(m)) , i.e. the number of distinct elements in a var-covar-matrix if none of the elements is zero.
    nrow(m)*(nrow(m)+1)/2
  } # end vpars
  
  is.diag <- function(m) {
  # tests if a matrix is a diagonal matrix, returns TRUE is yes, False otherwise
    
    return( all(m[!diag(nrow(m))] == 0))
    
  } # end is.diag
  
  count.vpars <- function(m) {
    
    if (is.diag(m)) return(nrow(m))
    else return(vpars(m))
    
  } # end count.vpars
  #---- end of local functions   
  
  model.df <- sum(sapply(VarCorr(model),count.vpars)) + length(fixef(model))     
  
  if(dfJ) {   
    rdf <- nrow( unique(model.frame(model)[,-1])  ) - model.df     # see Hosmer, Lemeshow: "Applied Logistic Reg", p138
  }
  else { 
    rdf <- nrow(model.frame(model))- model.df    # df = n - (p+1)
  }
  
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  
  prat <- Pearson.chisq/rdf
  Pearson.pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  
  data.frame(Pearson.chisq=Pearson.chisq,ratio=prat,Pearson.pvalue = Pearson.pval, model.df = model.df, r.df=rdf)
} # end overdisp.glmmadmb


plot.humm <- function(m) {
# creates residual plots of a hurdle model
# Argument: model m of class humm
  
  plot(fitted(m[[1]]), residuals(m[[1]]), main = "Prevalence", xlab = "Fitted" , ylab = "Pearson Residuals", pch=20)
  plot(sort( residuals(m[[1]])), main = "Prevalence", xlab = "Index", ylab = "Pearson Residuals sorted by size", pch=20)
  
#  plot(fitted(m[[2]]), m[[2]]$residuals, main = "Count", xlab = "Fitted" , ylab = "Pearson Residuals", pch=20)

# !! residuals(hm.trc[[2]]) does not work -> I need to check in Hilbe book the diagnostics for count data. Pearson resid for count models?
 
#  plot(sort( residuals(m[[2]])), main = "Count", xlab = "Index", ylab = "Pearson Residuals sorted by size", pch=20)
}  # end plot.humm
  


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####              diagnostics.humm         #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diagnostics.humm <- function(m, imodel = NULL){
# diagnostics of object of class humm
# Argument: object m of class humm
#           imodel is a count model of classglmmadmb. It should be the same model as the count- model of m but only with the intercept as "covariate". (for pseudo-R^2)
# returns:
  
  out <- list()
  
  out$pDisp <- overdisp(m[[1]])
  out$hosmerLem  <- hosmerLem(m[[1]]) 
  
  # ROC-curve: 
  out$ROC <- roc(m[[1]]$frame[, 1], fitted(m[[1]]) , auc = TRUE, ci = TRUE, plot = T ) 
  l <- out$ROC
  
  dfROC <- data.frame(specificity = l$specificities, sensitifity = l$sensitivities)
  out$ROC.plot <- ggplot(dfROC, aes(x = specificity, y =  sensitifity )) + geom_point() + geom_abline(slope = -1, intercept = 1) + ggtitle("ROC-Curve for Prevalence Model") 
  #plot(l$specificities, l$sensitivities, type = "l", xlim=c(1,0), ylim = c(0,1), main = "ROC Curve for Prevalence Model", xlab = "Specificity", ylab = "Sensitivity")
  #abline(1, -1)
  
  #### zero truncated count model ######
  df <- data.frame(obs = m[[2]]$frame[,1], expect = fitted(m[[2]]) )
  out$observed_expected.plot <- ggplot(df, aes(x = obs, y = expect ) ) + geom_point()+ ggtitle("Observed - Expected of Count Model")
  
  if (!(is.null(imodel))  ) out$Rp2 <- 1 - (exp( m[[2]]$loglik ) / (exp( imodel$loglik) ))
  
  
  class(out) <- "diagnostics.humm"
  out
  

} # end diagnostics

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####        print.diagnostics.humm          #####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print.diagnostics.humm <- function(x, ...) {
# print method for objects of class "diagnostics.humm"
  
  cat("Prevalence Model:\n") 
  cat("\n Chi-Sqare Test for Overdispersion: \n")
  print(x$pDisp)
  cat("\nHosmer - Lemeshow GOF Test: p-value ", x$hosmerLem$p.value, "\n", sep="")
  
  cat("Receiver Operating Characteristic (ROC), AUC: ", x$ROC$auc, ", CI(", x$ROC$ci[1], ", " , x$ROC$ci[3], ")", sep="")
  
  grid.arrange(x$observed_expected.plot, x$ROC.plot, ncol = 2, nrow = 1)
  
}

hl_grouped <- function(observed, predicted, group, g = 10){
  
  dfr <- data.frame(o = observed, p = predicted, g = group)
  
  no.groups <- length(unique(group))
  
  myTable <- dfr %>% group_by(g) %>% summarise(Predicted1 = sum(p), Predicted0 = sum(1 - p), 
                                               Observed1 = sum(o), Observed0 = sum(o == 0),
                                              N = n(), preProb = mean(p) )

  
  C <-  with(myTable, sum( ((Observed0 - Predicted0)^2) / (N * (1 - preProb)* preProb))  
                      + sum(  ((Observed1 - Predicted1)^2) / (N * preProb*(1 - preProb))  )   )
  C.p.value <- 1 - pchisq(C, df = no.groups - 2)
  
  # check: 
  CC <- with(myTable, sum( (Observed0 - Predicted0)^2  / Predicted0)  +  
                      sum( (Observed1 - Predicted1)^2  / Predicted1) )
        
  
  ordered.myTable <- myTable[order(myTable$preProb),]
  # ordered.myTable$gg <- rep(1:10, times = round(length(observed)))
  
  return(list( C = C, C.p.value = C.p.value, df = no.groups - 2, table = myTable ))
}
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
####         pseudo.R2                                ##########
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pseudo.R2 <- function(m) {
# computes the pseudo R2 for GLMM 
# as defined in  Nakagawa, Schielzeth (2013) "A general and simple method for obtaining R2 from generalized linear mixed-effect models"  
# Arguments: m model of class "glmerMod" or "glmmadmb"
# Returns:   R2 (marginal)
#            R2 (conditional)   
#            PCVs (not yet done)
  
extractSigmas <- function(m) {
# Returns the variance components of the random factors of model m    
    vm <- VarCorr(m)
    no.groups <- length(vm) 
    sigmas2 <- list()
    for(i in 1:no.groups) sigmas2[[i]] <- diag(vm[[i]])
    sigmas2[[no.groups + 1]] <- no.groups
    attributes(sigmas2)$names[no.groups + 1] <- "no.groups"
  
    
    return(sigmas2)
    
  } # end extractSigmas
# """"""""""""""""""""""""""""""""""  

varianceElements <- function(m) {
# returns the variance elements of a model   
  
  if (class(m) ==  "glmerMod" ) { 
          vc <- extractSigmas(m)
          
          # sum up the intercept variances of the different group variables:
          sigma2_l <- 0 
          for (j in 1:vc$no.groups) sigma2_l <- sigma2_l + vc[[j]][1]
          
          X <- model.matrix(m)  # design matrix
          b <- fixef(m)  # estimates of fixed effects coefficients
          
          sigma2_f <- as.numeric(var(X[,-1] %*% b[-1]))   
          
          fm <- family(m)
          if (fm[[1]] == "binomial" & fm[[2]] == "logit") { sigma2_d <- pi^2/3 
                                                            sigma2_e <- 0 } # assuming Boolean response   
          else  { sigma2_d <- "NA" 
                  sigma2_e <- "NA" }   
          
          
          return(list(s2_f = sigma2_f, s2_l = sigma2_l, s2_e = sigma2_e, s2_d = sigma2_d))
   
  } # end if (class == glmerMod)
  
#  case glmmadmb:  
  if (class(m) ==  "glmmadmb" ) { 
    vc <- extractSigmas(m)
    
    # sum up the intercept variances of the different group variables:
    sigma2_l <- 0 
    for (j in 1:vc$no.groups) sigma2_l <- sigma2_l + vc[[j]][1]
    
    X <- model.matrix(m)  # design matrix
    b <- fixef(m)  # estimates of fixed effects coefficients
    
    sigma2_f <- as.numeric(var(X[,-1] %*% b[-1]))  
    
    sigma2_d <- "NA" 
    sigma2_e <- "NA"
 
    if (m$family == "binom" & m$link == "logit") { sigma2_d <- pi^2/3 
                                                      sigma2_e <- 0 } # assuming Boolean response   
    
    if (m$family == "truncnbinom" & m$link == "log") { b0 <- m$b[1]
                                                       sigma2_d <-  log(1/exp(b0) + 1)    # ??????????????????????
                                                       epsilon  <- log(model.frame(m)[,1]) - log(fitted(m))   # ?? finding epsilon by hand ?
                                                       sigma2_e <-  var(epsilon)  } 
                                                             
    
    
    return(list(s2_f = sigma2_f, s2_l = sigma2_l, s2_e = sigma2_e, s2_d = sigma2_d))
    
  } # end if (class == glmerMod)
  
}  # end of varianceElements  
# """"""""""""""""""""""""""""""""""
  ve <- varianceElements(m)

  R2.marginal <- ve$s2_f / (ve$s2_f + ve$s2_l + ve$s2_e + ve$s2_d)
  R2.conditional <- (ve$s2_f + ve$s2_l )  / (ve$s2_f + ve$s2_l + ve$s2_e + ve$s2_d)  
  
return(list(variance_elements = ve, R2.marginal = R2.marginal, R2.conditional = R2.conditional ))
  
}  # end pseudo.R2


# logLik.humm <- function(hm)
# anova.humm <- function(hm)
# AIC.humm <- function(hm)
# diagMNB <- function(NBmodel)
# diagMLogit <- function(logitm)
# diag.humm 
# function like logisCorrect for marginal models or approximation of them.



