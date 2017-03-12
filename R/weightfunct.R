#' Estimate the Vevea and Hedges (1995) Weight-Function Model
#'
#' This function allows the user to estimate the Vevea and Hedges (1995) weight-function model for publication bias.
#' @param effect a vector of meta-analytic effect sizes.
#' @param v a vector of meta-analytic sampling variances; needs to match up with the vector of effects, such that the first element in the vector of effect sizes goes with the first element in the vector of sampling variances, and so on.
#' @param pval defaults to \code{NULL}. A vector containing observed p-values for the corresponding effect sizes. If not provided, p-values are calculated.
#' @param steps a vector of p-value cutpoints. The default only distinguishes between significant and non-significant effects (p < 0.05).
#' @param mods defaults to \code{NULL}. A formula specifying the linear model.
#' @param weights defaults to \code{FALSE}. A vector of prespecified weights for p-value cutpoints to estimate the Vevea and Woods (2005) model.
#' @param fe defaults to \code{FALSE}. Indicates whether to estimate a fixed-effects model.
#' @param table defaults to \code{FALSE}. Indicates whether to print a table of the p-value intervals specified and the number of effect sizes per interval.
#' @importFrom stats model.matrix optim pchisq pnorm qnorm model.frame na.action na.omit
#' @keywords weightr
#' @details This function allows meta-analysts to estimate both the
#' weight-function model for publication bias that was originally published in
#' Vevea and Hedges (1995) and the modified version presented in Vevea and Woods
#' (2005). Users can estimate both of these models with and without predictors and
#' in random-effects or fixed-effects situations.
#'
#' The Vevea and Hedges (1995) weight-function model is a tool for modeling publication
#' bias using weighted distribution theory. The model first estimates an unadjusted
#' fixed-, random-, or mixed-effects model, where the observed effect sizes are
#' assumed to be normally distributed as a function of predictors. This unadjusted
#' model is no different from the traditional meta-analytic model. Next, the Vevea
#' and Hedges (1995) weight-function model estimates an adjusted model that includes
#' not only the original mean model, fixed-, random-, or mixed-effects, but a series
#' of weights for any pre-specified p-value intervals of interest. This produces mean,
#' variance component, and covariate estimates adjusted for publication bias, as well
#' as weights that reflect the likelihood of observing effect sizes in each specified
#' interval.
#'
#' It is important to remember that the weight for each
#' estimated p-value interval must be interpreted relative to the first interval,
#' the weight for which is fixed to 1 so that the model is identified. In other
#' words, a weight of 2 for an interval indicates that effect sizes in that p-value
#' interval are about twice as likely to be observed as those in the first interval.
#' Finally, it is also important to remember that the model uses p-value cutpoints
#' corresponding to one-tailed p-values. This allows flexibility in the selection
#' function, which does not have to be symmetric for effects in the opposite direction;
#' a two-tailed p-value of 0.05 can therefore be represented as p < .025 or p > .975.
#'
#' After both the unadjusted and adjusted meta-analytic models are estimated, a
#' likelihood-ratio test compares the two. The degrees of freedom for this test are
#' equal to the number of weights being estimated. If the likelihood-ratio test is
#' significant, this indicates that the adjusted model is a better fit for the data,
#' and that publication bias may be a concern.
#'
#' To estimate a large number of weights for p-value intervals, the Vevea and Hedges
#' (1995) model works best with large meta-analytic datasets. It may have trouble
#' converging and yield unreliable parameter estimates if researchers, for instance,
#' specify a p-value interval that contains no observed effect sizes. However,
#' meta-analysts with small datasets are still likely to be interested in assessing
#' publication bias, and need tools for doing so. Vevea and Woods (2005)
#' attempted to solve this problem by adapting the Vevea and Hedges (1995) model to
#' estimate fewer parameters. The meta-analyst can specify p-value cutpoints,
#' as before, and specify corresponding fixed weights for those cutpoints. Then the
#' model is estimated. For the adjusted model, only the variance component and mean
#' model parameters are estimated, and they are adjusted relative to the fixed weights.
#' For example, weights of 1 for each p-value interval specified describes a situation
#' where there is absolutely no publication bias, in which the adjusted estimates are
#' identical to the unadjusted estimates. By specifying weights that depart from 1 over various p-value intervals, meta-analysts can
#' examine how various one-tailed or two-tailed selection patterns would alter their
#' effect size estimates. If changing the pattern of weights drastically changes
#' the estimated mean, this is evidence that the data may be vulnerable to
#' publication bias.
#'
#' For more information, consult the papers listed in the References section here.
#' Also, feel free to email the maintainer of \code{weightr} at kcoburn@ucmerced.edu.
#' The authors are currently at work on a detailed package tutorial, which we
#' hope to publish soon.
#' @export
#' @return The function returns a list containing the following components: \code{output_unadj}, \code{output_adj}, \code{steps}, \code{mods}, \code{weights}, \code{fe}, \code{table}, \code{effect}, \code{v}, \code{npred}, \code{nsteps}, \code{p}, \code{XX}, \code{removed}.
#'
#' The results of the unadjusted and adjusted models are returned by selecting the first (\code{[[1]]}) and second (\code{[[2]]}) elements of the list, respectively. The parameters can be obtained by \code{[[1]]$par} or \code{[[2]]$par}. The order of parameters is as follows: variance component, mean or linear coefficients, and weights. (Note that if \code{weights}
#' are specified using the Vevea and Woods (2005) model, no standard errors, p-values, z-values, or confidence intervals
#' are provided for the adjusted model, as these are no longer meaningful. Also note that the variance component is not reported for fixed-effects models.)
#' \describe{
#'    \item{\code{unadj_est}}{the unadjusted model estimates}
#'    \item{\code{adj_est}}{the adjusted model estimates}
#'    \item{\code{steps}}{the specified p-value cutpoints}
#'    \item{\code{mods}}{the linear model formula, if one is specified}
#'    \item{\code{weights}}{the vector of weights for the Vevea and Woods (2005) model, if specified}
#'    \item{\code{fe}}{indicates whether or not a fixed-effects model was estimated}
#'    \item{\code{table}}{indicates whether a sample size table was produced}
#'    \item{\code{effect}}{the vector of effect sizes}
#'    \item{\code{v}}{the vector of sampling variances}
#'    \item{\code{npred}}{the number of predictors included; 0 represents an intercept-only model}
#'    \item{\code{nsteps}}{the number of p-value cutpoints}
#'    \item{\code{p}}{a vector of p-values for the observed effect sizes}
#'    \item{\code{XX}}{the model matrix; the first column of ones represents the intercept, and any other columns correspond to moderators}
#'    \item{\code{removed}}{effect sizes with missing data are removed by listwise deletion; any removed are provided here. Defaults to \code{NULL}}
#'    }
#'
#' @references Coburn, K. M. & Vevea, J. L. (2015). Publication bias as a function
#' of study characteristics. Psychological Methods, 20(3), 310.
#'
#' Vevea, J. L. & Hedges, L. V. (1995). A general linear model for
#' estimating effect size in the presence of publication bias. Psychometrika, 60(3),
#' 419-435.
#'
#' Vevea, J. L. & Woods, C. M. (2005). Publication bias in research synthesis:
#' Sensitivity analysis using a priori weight functions. Psychological Methods, 10(4),
#' 428-443.
#' @examples
#' \dontrun{
#' # Uses the default p-value cutpoints of 0.05 and 1:
#'
#' weightfunct(effect, v)
#'
#' # Estimating a fixed-effects model, again with the default cutpoints:
#'
#' weightfunct(effect, v, fe=TRUE)
#'
#' # Specifying cutpoints:
#'
#' weightfunct(effect, v, steps=c(0.01, 0.025, 0.05, 0.10, 0.20, 0.30, 0.50, 1.00))
#'
#' # Including a linear model, where moderators are denoted as 'mod1' and mod2':
#'
#' weightfunct(effect, v, mods=~mod1+mod2)
#'
#' # Specifying cutpoints and weights to estimate Vevea and Woods (2005):
#'
#' weightfunct(effect, v, steps=c(0.01, 0.05, 0.50, 1.00), weights=c(1, .9, .7, .5))
#'
#' # Specifying cutpoints and weights while including a linear model:
#'
#' weightfunct(effect, v, mods=~mod1+mod2, steps=c(0.05, 0.10, 0.50, 1.00), weights=c(1, .9, .8, .5))
#' }

weightfunct <- function(effect, v, steps=c(0.025,1.00), mods=NULL,
                        weights=NULL, fe=FALSE, table=FALSE, pval=NULL){

  neglike_unadj <- function(pars) {
    if(fe == FALSE){
      vc = pars[1]
      beta = pars[2:(npred+2)]

      mn = XX%*%beta
      eta = sqrt(v + vc)
      b = 1/2 * sum(((effect-mn)/eta)^2)
      c = sum(log(eta))
    }
    else{
      beta = pars[1:(npred+1)]
      mn = XX%*%beta
      eta = sqrt(v)
      b = 1/2 * sum(((effect-mn)/eta)^2)
      c = sum(log(eta))
    }
    return(b+c)
  }

  neglike_adj <- function(pars) {
    if(fe == FALSE){
      vc = pars[1]
      beta = pars[2:(npred+2)]
      if(is.null(weights) == FALSE){
        w = weights
      }
      else{
        w = c(1,pars[(npred+3):( (nsteps - 2) + (npred+3) )])
      }
      contrib = log(w[wt])

      mn = XX%*%beta
      a = sum(contrib)
      eta = sqrt(v + vc)
    }
    else{
      beta = pars[1:(npred+1)]
      if(is.null(weights) == FALSE){
        w = weights
      }
      else{
        w = c(1,pars[(npred+2):( (nsteps - 2) + (npred+2) )])
      }
      contrib = log(w[wt])

      mn = XX%*%beta
      a = sum(contrib)
      eta = sqrt(v)
    }
    b = 1/2 * sum(((effect-mn)/eta)^2)
    c = sum(log(eta))
    Bij <- matrix(rep(0,number*nsteps),nrow=number,ncol=nsteps)
    bi = -si * qnorm(steps[1])
    Bij[,1] = 1-pnorm((bi-mn)/eta)
    if(nsteps > 2){
      for(j in 2:(length(steps)-1)) {
        bi = -si * qnorm(steps[j])
        bilast = -si * qnorm(steps[j-1])
        Bij[,j] = pnorm((bilast-mn)/eta) - pnorm((bi-mn)/eta)
      }
    }
    bilast = -si * qnorm(steps[length(steps)-1])
    Bij[,length(steps)] = pnorm((bilast-mn)/eta)

    swbij = 0
    for(j in 1:length(steps)) swbij = swbij + w[j]*Bij[,j]
    d = sum(log(swbij))

    return(-a + b + c + d)
  }
  
  if(is.null(mods)){
    npred <- 0
    data <- data.frame(effect, v)
  }
  else{
    if(typeof(mods)=="language"){
        XX <- model.matrix(mods, model.frame(mods, na.action='na.pass'))
        npred <- dim(XX)[2]-1
        data <- data.frame(effect, v, XX)
    }
  }
  if(any(is.na(data))){
    data <- na.omit(data)
    removed <- as.numeric(na.action(data))
  }
  else{
    removed <- NULL
  }
  effect <- data[,1]
  v <- data[,2]
  if(npred == 0){
    XX <- cbind(rep(1,length(effect)))
  }
  else{
    XX <- as.matrix(data[,(3:(npred+3))])
  }

  if(length(effect)!=length(v)){
    stop('Your vector of effect sizes and your vector of sampling variances are not the same length. Please check your data.')
  }

  if(identical(effect,v)){
    stop('Your vector of effect sizes is exactly the same as your vector of sampling variances. Please check your data.')
  }

  if(min(v) < 0){
    stop('Sampling variances cannot be negative. Please check your data.')
  }

  si <- sqrt(v)
  if(is.null(pval)){
    p <- 1-pnorm(effect/sqrt(v))
  }
  else{
    p <- pval
  }
  if(max(steps)!=1){
    steps <- c(steps,1)
  }
  if(max(steps) > 1){
    stop('p-value cutpoints cannot be greater than 1.')
  }
  if(min(steps) < 0){
    stop('p-value cutpoints cannot be negative.')
  }
  if(length(unique(steps)) != length(steps)){
    stop('Two or more p-value cutpoints are identical.')
  }

  if(is.null(weights)){
    steps <- sort(steps)
  }
  if(is.null(weights) == FALSE){
    if(min(weights) < 0){
      stop('Weights for p-value intervals cannot be negative.')
    }
    if(length(weights)!=length(steps)){
      stop('The number of weights does not match the number of p-value intervals created.')
    }
    new <- cbind(steps, weights)
    steps <- new[order(steps),1]
    weights <- new[order(steps),2]
  }

  number <- length(effect)
  nsteps <- length(steps)


  wt <- rep(1,number)
  for(i in 1:number) {
    if(p[i] <= steps[1]) wt[i] = 1
    for(j in 2:nsteps) {
      if (steps[j-1] < p[i] && p[i] <= steps[j]) wt[i] = j
    }
    if(  p[i] > steps[nsteps-1] ) wt[i] = nsteps
  }

  intervaltally <- function(p, steps) {
    p1 <- cut(p, breaks=c(-Inf,steps), labels=steps)
    return(p1) }
  pvalues <- as.numeric(table(intervaltally(p, steps)))

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

  if(sum(table(intervaltally(p,steps)) == 0) >= 1){
    warning('At least one of the p-value intervals contains no effect sizes, leading to estimation problems. Consider re-specifying the cutpoints.')
  }

  if(sum( table(intervaltally(p,steps)) > 0 & table(intervaltally(p, steps)) <= 3) >= 1){
    warning('At least one of the p-value intervals contains three or fewer effect sizes, which may lead to estimation problems. Consider re-specifying the cutpoints.')
  }

  if(is.null(mods)){

    if(fe == FALSE){
  
      pars <- c(mean(v)/4, mean(effect), rep(0,(nsteps-1)))

      output_unadj <- optim(par=pars[1:2],fn=neglike_unadj,lower=c(0,-Inf),method="L-BFGS-B",hessian=TRUE)

      output_adj <- optim(par=pars,fn=neglike_adj,lower=c(0,-Inf, rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)
#       output_adj <- optimx(par=pars,fn=neglike_adj,lower=c(0,-Inf, rep(0.01,(nsteps-1))),hessian=TRUE)

      results <- list(output_unadj,output_adj, steps=steps, mods=mods, weights=weights, fe=fe, table=table, effect=effect, v=v, npred=npred, nsteps=nsteps, p=p, XX=XX, removed=removed)
#       print(output_adj)

#       model.print(results)

      if(is.null(weights)){
        suppressWarnings(
        results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                         adj_est=c(output_adj$par),
                         adj_se=c(sqrt(diag(solve(output_adj$hessian)))),
                         z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         z_adj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         p_adj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.lb_adj=c(output_adj$par - qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))),
                         ci.ub_adj=c(output_adj$par + qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))))
        )
      }
      if(is.null(weights) == FALSE){
        results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                         adj_est=c(output_adj$par),
                         z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))))
      }

    }

    if(fe == TRUE){

      pars <- c(mean(effect), rep(1,(nsteps-1)))

      output_unadj <- optim(par=pars[1],fn=neglike_unadj,lower=c(-Inf),method="L-BFGS-B",hessian=TRUE)

      output_adj <- optim(par=pars,fn=neglike_adj,lower=c(-Inf, rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

      results <- list(output_unadj,output_adj, steps=steps, mods=mods, weights=weights, fe=fe, table=table, effect=effect, v=v, npred=npred, nsteps=nsteps, p=p, removed=removed)

#       model.print(results)

      if(is.null(weights)){
        suppressWarnings(
        results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                         adj_est=c(output_adj$par),adj_se=c(sqrt(diag(solve(output_adj$hessian)))),
                         z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         z_adj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         p_adj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.lb_adj=c(output_adj$par - qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))),
                         ci.ub_adj=c(output_adj$par + qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))))
        )
      }
      if(is.null(weights) == FALSE){
        results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                         adj_est=c(output_adj$par),
                         z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                         p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                         ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                         ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))))
      }

    }


  }

  else{

    if(typeof(mods)=="language"){

      if(fe == FALSE){

        pars <- c(mean(v)/4, mean(effect), rep(0, npred), rep(1, (nsteps - 1)))

        output_unadj <- optim(par=pars[1:(npred+2)],fn=neglike_unadj,lower=c(0,rep(-Inf, (npred+1))),method="L-BFGS-B",hessian=TRUE)

        output_adj <- optim(par=pars,fn=neglike_adj,lower=c(0,rep(-Inf, (npred+1)),rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

        results <- list(output_unadj,output_adj, steps=steps, mods=mods, weights=weights, fe=fe, table=table, effect=effect, v=v, npred=npred, nsteps=nsteps, p=p, XX=XX, removed=removed)

#         model.print(results)

        if(is.null(weights)){
          suppressWarnings(
          results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                           adj_est=c(output_adj$par),adj_se=c(sqrt(diag(solve(output_adj$hessian)))),
                           z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           z_adj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           p_adj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.lb_adj=c(output_adj$par - qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))),
                           ci.ub_adj=c(output_adj$par + qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))))
          )
        }
        if(is.null(weights) == FALSE){
          suppressWarnings(
          results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                           adj_est=c(output_adj$par),
                           z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))))
          )
        }

      }

      if(fe == TRUE){

        pars <- c(mean(effect), rep(0, npred), rep(1, (nsteps - 1)))

        output_unadj <- optim(par=pars[1:(npred+1)],fn=neglike_unadj,lower=c(rep(-Inf, (npred+1))),method="L-BFGS-B",hessian=TRUE)

        output_adj <- optim(par=pars,fn=neglike_adj,lower=c(rep(-Inf, (npred+1)),rep(0.01,(nsteps-1))),method="L-BFGS-B",hessian=TRUE)

        results <- list(output_unadj,output_adj, steps=steps, mods=mods, weights=weights, fe=fe, table=table, effect=effect, v=v, npred=npred, nsteps=nsteps, p=p, XX=XX, removed=removed)

        # model.print(results)

        if(is.null(weights)){
          suppressWarnings(
          results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                           adj_est=c(output_adj$par),adj_se=c(sqrt(diag(solve(output_adj$hessian)))),
                           z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           z_adj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           p_adj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.lb_adj=c(output_adj$par - qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))),
                           ci.ub_adj=c(output_adj$par + qnorm(0.975) * sqrt(diag(solve(output_adj$hessian)))))
          )
        }
        if(is.null(weights) == FALSE){
          suppressWarnings(
          results2 <- list(unadj_est=c(output_unadj$par),unadj_se=c(sqrt(diag(solve(output_unadj$hessian)))),
                           adj_est=c(output_adj$par),
                           z_unadj=c(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))),
                           p_unadj=c(2*pnorm(-abs(output_unadj$par/sqrt(diag(solve(output_unadj$hessian)))))),
                           ci.lb_unadj=c(output_unadj$par - qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))),
                           ci.ub_unadj=c(output_unadj$par + qnorm(0.975) * sqrt(diag(solve(output_unadj$hessian)))))
          )
        }

      }

    }

    else{

      stop('Moderators must be entered as "mods= ~ x1 + x2"')

    }

  }

  # attr(results, 'Test') <- results2
  
  
  
  # This should be uncommented:
  # invisible(results2)
  
  class(results) <- c("weightfunct")
  return(results)
  
  
  
  #model.print(results)
  # return(list(results2,model.print(results)))
  
#   return(list(c(results2,model.print(results))))

}
