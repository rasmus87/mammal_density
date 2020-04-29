MCMC.predict.v2 <- function(object, newdata) {
  object2 <- MCMCglmm(fixed=object$Fixed$formula, 
                      random=object$Random$formula, 
                      rcov=object$Residual$formula, 
                      family=object$Residual$original.family,
                      data=newdata, 
                      nitt=1, 
                      thin=1,
                      burnin=0, 
                      ginverse=object$ginverse, 
                      verbose=FALSE, 
                      pr=TRUE)
  
  W <- cbind(object2$X, object2$Z)
  post.pred <- t(apply(object$Sol, 1, function(x){(W %*% x)@x}))[, 1:n.mam]
  
  se <- apply(post.pred, 2, sd)
  pred <- matrix(colMeans(post.pred), dim(post.pred)[2], 1)
  pred <- cbind(pred, coda::HPDinterval(mcmc(post.pred), prob=0.95), se)
  colnames(pred) <- c("fit", "lwrCI", "uprCI", "se")
  rownames(pred) <- 1:dim(pred)[1]
  
  colnames(post.pred) <- mam$Binomial.1.2
  return(list(post.pred, pred))
}