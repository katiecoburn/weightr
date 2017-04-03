#' Print Model Results
#'
#' This function allows you to print the model results.
#' @keywords weightr
#' @export
#' @importFrom stats model.matrix optim pchisq pnorm qnorm
#' @examples
#' \dontrun{
#' print.weightfunct(weightfunct(d,v))
#' }
print.weightfunct <- function(x){
  if (!inherits(x, "weightfunct"))
    stop("Argument 'x' must be an object of class \"weightfunct\".")
  
  ####### Unadjusted model ########
  
    digits <- 4
    cat("\n")
    cat("Unadjusted Model (k = ", length(x$effect), "):", sep="")
    cat("\n\n")
    if(x$fe == FALSE){
      cat("tau^2 (estimated amount of total heterogeneity): ", formatC(round(x[[1]]$par[1], digits = digits), digits = digits, format = "f"), " (SE = ", formatC(round(sqrt(diag(solve(x[[1]]$hessian)))[1], digits = digits),digits = digits, format = "f"), ")", sep="")
      cat("\n")
      cat("tau (square root of estimated tau^2 value): ", formatC(round(sqrt(x[[1]]$par[1]), digits = digits),digits = digits, format = "f"))
      cat("\n\n")
    }
    cat("Model Results:")
    cat("\n\n")
    if(x$fe == FALSE){
      unadj_est <- cbind(x[[1]]$par[2:(x$npred+2)])
      unadj_se <- cbind(sqrt(diag(solve(x[[1]]$hessian)))[2:(x$npred+2)])
      z_stat <- unadj_est/unadj_se
      p_val <- (2*pnorm(-abs(z_stat)))
      ci.lb <- unadj_est - qnorm(0.975) * unadj_se
      ci.ub <- unadj_est + qnorm(0.975) * unadj_se
    }
    if(x$fe == TRUE){
      unadj_est <- cbind(x[[1]]$par[1:(x$npred+1)])
      unadj_se <- cbind(sqrt(diag(solve(x[[1]]$hessian)))[1:(x$npred+1)])
      z_stat <- unadj_est/unadj_se
      p_val <- (2*pnorm(-abs(z_stat)))
      ci.lb <- unadj_est - qnorm(0.975) * unadj_se
      ci.ub <- unadj_est + qnorm(0.975) * unadj_se
    }
    res.table <- data.frame(matrix(c(unadj_est, unadj_se, z_stat, p_val, ci.lb, ci.ub), nrow=(x$npred+1), byrow=F),stringsAsFactors=FALSE)
    rowlabels <- rep(0, (x$npred+1))
    rowlabels[1] <- "Intercept"
    if(x$npred > 0){
      for(i in 2:(x$npred+1)){
        rowlabels[i] <- paste(c(colnames(x$XX)[i]))
      }
    }
    row.names(res.table) <- c(rowlabels)
    colnames(res.table) <- c("estimate","std.error","z-stat","p-val","ci.lb","ci.ub")
    res.table[,4] <- format.pval(res.table[,4])
    res.table[,c(1,2,3,5,6)] <- format(res.table[,c(1,2,3,5,6)], digits=4)
    print(res.table)
    
    ####### Adjusted model ########
    
      cat("\n")
      cat("Adjusted Model (k = ", length(x$effect), "):", sep="")
      cat("\n\n")
      if(x$fe == FALSE){
        if(is.null(x$weights)){
          if(is.nan(suppressWarnings(sqrt(diag(solve(x[[2]]$hessian)))[1]))){
            warning('The adjusted variance component is so close to zero that a border condition prevents a meaningful iterative solution. As long as other model estimates are still \nreasonable, the results are identical to those from a fixed-effects analysis.')
          }
          suppressWarnings(
          cat("tau^2 (estimated amount of total heterogeneity): ", formatC(round(x[[2]]$par[1], digits = digits),digits = digits, format = "f"), " (SE = ", formatC(round(sqrt(diag(solve(x[[2]]$hessian)))[1], digits = digits),digits = digits, format = "f"), ")", sep="")
          )
        }
        if(is.null(x$weights) == FALSE){
          cat("tau^2 (estimated amount of total heterogeneity): ", formatC(round(x[[2]]$par[1], digits = digits),digits = digits, format = "f"), " (SE = ", "---", ")", sep="")
        }
        cat("\n")
        cat("tau (square root of estimated tau^2 value): ", formatC(round(sqrt(x[[2]]$par[1]), digits = digits),digits = digits, format = "f"))
        cat("\n\n")
      }
      cat("Model Results:")
      cat("\n\n")
      
      if(x$fe == FALSE){
        if(is.null(x$weights)){
          adj_int_est <- cbind(x[[2]]$par[2:( (x$nsteps - 1) + (x$npred+2) )])
          suppressWarnings(
          adj_int_se <- cbind(sqrt(diag(solve(x[[2]]$hessian)))[2:( (x$nsteps - 1) + (x$npred+2) )]))
        }
        else{
          adj_int_est <- cbind(c(x[[2]]$par[2:( (x$npred+2) )], x$weights[2:length(x$weights)]))
          adj_int_se <- cbind(rep("---", length(x[[2]]$par[2:length(x[[2]]$par)])))
        }
      }
      if(x$fe == TRUE){
        if(is.null(x$weights)){
          adj_int_est <- cbind(x[[2]]$par[1:( (x$nsteps - 1) + (x$npred+1) )])
          suppressWarnings(
          adj_int_se <- cbind(sqrt(diag(solve(x[[2]]$hessian)))[1:( (x$nsteps - 1) + (x$npred+1) )]))
        }
        else{
          adj_int_est <- cbind(c(x[[2]]$par[1:( (x$npred+1) )], x$weights[2:length(x$weights)]))
          adj_int_se <- cbind(rep("---", length(x[[2]]$par[1:length(x[[2]]$par)])))
        }
      }
      
      if(is.null(x$weights)){
        z_stat_int <- adj_int_est/adj_int_se
        p_val_int <- (2*pnorm(-abs(z_stat_int)))
        ci.lb_int <- adj_int_est - qnorm(0.975) * adj_int_se
        ci.ub_int <- adj_int_est + qnorm(0.975) * adj_int_se
      }else{
        z_stat_int <- rep("---", length(x[[2]]$par[2:length(x[[2]]$par)]))
        p_val_int <- rep("---", length(x[[2]]$par[2:length(x[[2]]$par)]))
        ci.lb_int <- rep("---", length(x[[2]]$par[2:length(x[[2]]$par)]))
        ci.ub_int <- rep("---", length(x[[2]]$par[2:length(x[[2]]$par)]))
      }
      res.table <- data.frame(matrix(c(adj_int_est, adj_int_se, z_stat_int, p_val_int, ci.lb_int, ci.ub_int), nrow=(x$npred+1+(x$nsteps-1)), byrow=F),stringsAsFactors=FALSE)
      
      rowlabels1 <- rep(0, (x$npred+1))
      rowlabels1[1] <- "Intercept"
      if(x$npred > 0){
        for(i in 2:length(rowlabels1)){
          rowlabels1[i] <- paste(c(colnames(x$XX)[i]))
        }
      }
      rowlabels2 <- rep(0, (x$nsteps-1))
      for(i in 1:(length(rowlabels2))){
        rowlabels2[i] <- paste(c(x$steps[i], "< p <", x$steps[i + 1]), collapse=" ")
      }
      row.names(res.table) <- c(rowlabels1,rowlabels2)
      colnames(res.table) <- c("estimate","std.error","z-stat","p-val","ci.lb","ci.ub")
      if(is.null(x$weights)){
        res.table[,"p-val"] <- format.pval(res.table[,"p-val"])
      }
      res.table[,c(1,2,3,5,6)] <- format(res.table[,c(1,2,3,5,6)], digits=4)
      print(res.table)
      
      ####### LRT ########
      
        if(is.null(x$weights)){
          cat("\n")
          cat("Likelihood Ratio Test:")
          cat("\n")
          df <- length(x[[2]]$par) - length(x[[1]]$par)
          lrchisq <- 2*(abs(x[[1]]$value - x[[2]]$value))
          pvalue <- 1-pchisq(lrchisq,df)
          cat("X^2(df = ", df, ") = ", lrchisq, ", p-val = ", format.pval(pvalue), sep="")
        }else{
          cat("\n")
          cat("Note: The symbol --- appears because the user has specified weights,\nchoosing to use the Vevea and Woods model, which does not estimate \nweights for p-value intervals and therefore cannot produce meaningful \nstandard errors. The likelihood ratio test is also not interpretable.")
        }
      
      ####### Interval table ########
      
        
        if(x$table == TRUE){
          intervaltally <- function(p, steps) {
            p1 <- cut(p, breaks=c(-Inf,steps), labels=steps)
            return(p1) }
          pvalues <- as.numeric(table(intervaltally(x$p, x$steps)))
          
          sampletable <- function(p, pvalues, steps){
            nsteps <- length(steps)
            results <- matrix(nrow=length(pvalues),ncol=1)
            results[,1] <- pvalues
            rowlabels <- c(0, length(results[,1]))
            rowlabels[1] <- paste(c("p-values <", steps[1]), collapse="")
            for(i in 2:nsteps){
              rowlabels[i] <- paste(c(steps[i - 1], "< p-values <", steps[i]), collapse=" ")
            }
            resultsb <- data.frame(results, row.names=c(rowlabels))
            colnames(resultsb) <- c("Frequency")
            return(resultsb)
          }
          cat("\n")
          cat("\n")
          cat("Number of Effect Sizes per Interval:")
          cat("\n")
          cat("\n")
          format(print(sampletable(x$p, pvalues, x$steps)))
        }
      
      if(is.null(x$removed)==FALSE){
        cat("\n")
        cat("There were ", length(x$removed), "cases removed from your dataset due to the presence of missing data. To view the row numbers of these cases, use the attribute '$removed'.")
      }
      
      
      
    invisible()
}