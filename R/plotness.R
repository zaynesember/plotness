#' Poisonness, Binomialness, and Negative Binomialness Plots
#' @author Zayne Sember <zsember@ucsd.edu>
#' @description  Creates a diagnostic distribution plot or dataframe. The function is adapted
#' from the vcd library's distplot function to produce ggplots and to allow the
#' option to return only the data needed for the plot so the user's preferred
#' graphics package can be used.
#' @importFrom vcd goodfit
#' @import ggplot2
#' @import stats
#' @param x either a vector of counts, a 1-way table of frequencies of counts or
#'  a data frame or matrix with frequencies in the first column and the
#'  corresponding counts in the second column.
#' @param type a character string indicating the distribution.
#' @param size the size argument for the binomial and negative binomial
#'  distribution. If set to NULL and type is "binomial", then size is taken to
#'  be the maximum count. If set to NULL and type is "nbinomial", then size is
#'  estimated from the data.
#' @param lambda parameter of the poisson distribution. If type is "poisson"
#'  and lambda is specified a leveled Poissonness plot is produced.
#' @param plot logical. When `TRUE` a ggplot is returned; otherwise a data frame
#'  containing data necessary for a plot is returned.
#' @param title a character string to be used as the title of the plot. If set
#'  to `NULL` a generic title is produced based on the distribution type.
#' @param legend logical. Should a legend be plotted?
#' @param label logical. Should fit data be printed on the plot?
#' @param conf_int logical. Should confidence intervals be plotted?
#' @param conf_level confidence level for confidence intervals.
#' @param xlim limits for the x-axis.
#' @param xlab a label for the x-axis.
#' @param ylab a label for the y-axis.
#' @return If `plot`=`TRUE`, returns a ggplot object; otherwise returns a
#'  data frame containing the data necessary to produce a plot.
#' @examples
#' plotness(rpois(15,10),type="poisson");
#' plotness(rbinom(15,10,prob=0.5),type="binomial");
#' plotness(rnbinom(15,10,prob=0.5),type="binomial");
#' plotness(rpois(15,10),type="poisson", plot=FALSE);
#' @export
plotness <-
  function(x, type = c("poisson", "binomial", "nbinomial"),
           size = NULL, lambda = NULL, plot=TRUE, title=NULL, legend = TRUE,
           label = TRUE, conf_int = TRUE, conf_level = 0.95, xlim = NULL,
           xlab = "Number of occurrences", ylab = "Distribution metameter")
  {
    if(is.vector(x)) {
      x <- table(x)
    }
    if(is.table(x)) {
      if(length(dim(x)) > 1) stop ("x must be a 1-way table")
      freq <- as.vector(x)
      count <- as.numeric(names(x))
    } else {
      if(!(!is.null(ncol(x)) && ncol(x) == 2))
        stop("x must be a 2-column matrix or data.frame")
      freq <- as.vector(x[,1])
      count <- as.vector(x[,2])
    }

    myindex <- (1:length(freq))[freq > 0]
    mycount <- count[myindex]
    myfreq <- freq[myindex]

    switch(match.arg(type),

           "poisson" = {
             par.ml <- suppressWarnings(goodfit(x, type = type)$par$lambda)

             phi <- function(nk, k, N, size = NULL)
               ifelse(nk > 0, lgamma(k + 1) + log(nk/N), NA)
             y <- phi(myfreq, mycount, sum(freq))
             if(!is.null(lambda)) y <- y + lambda - mycount * log(lambda)
             fm <- lm(y ~ mycount)
             par.estim <- exp(coef(fm)[2])
             names(par.estim) <- "lambda"
             txt <- "exp(slope)"
             if(!is.null(lambda)) {
               par.estim <- par.estim * lambda
               txt <- paste(txt, "x lambda")
             }
             legend.text <- paste(txt, "=", round(par.estim, digits = 3))
           },

           "binomial" = {
             if(is.null(size)) {
               size <- max(count)
               warning("size was not given, taken as maximum count")
             }
             par.ml <- suppressWarnings(goodfit(x, type = type, par = list(size = size))$par$prob)

             phi <- function(nk, k, N, size = NULL)
               log(nk) - log(N * choose(size, k))
             y <- phi(myfreq, mycount, sum(freq), size = size)
             fm <- lm(y ~ mycount)
             par.estim <- exp(coef(fm)[2])
             par.estim <- par.estim / (1 + par.estim)
             names(par.estim) <- "prob"
             legend.text <- paste("inv.logit(slope) =", round(par.estim, digits = 3))
           },

           "nbinomial" = {
             if(is.null(size)) {
               par.ml <- suppressWarnings(goodfit(x, type = type)$par)
               size <- par.ml$size
               par.ml <- par.ml$prob
             }else{
               xbar <- weighted.mean(mycount, myfreq)
               par.ml <- size / (size+xbar)
             }
             phi <- function(nk, k, N, size = NULL)
               log(nk) - log(N * choose(size + k - 1, k))
             y <- phi(myfreq, mycount, sum(freq), size = size)
             fm <- lm(y ~ mycount)
             par.estim <- 1 - exp(coef(fm)[2])
             names(par.estim) <- "prob"
             legend.text <- paste("1-exp(slope) =", round(par.estim, digits = 3))
           })

    yhat <- ifelse(myfreq > 1.5, myfreq - 0.67, 1/exp(1))
    yhat <- phi(yhat, mycount, sum(freq), size = size)
    if(!is.null(lambda)) yhat <- yhat + lambda - mycount * log(lambda)

    phat <- myfreq / sum(myfreq)
    ci.width <- qnorm(1-(1 - conf_level)/2) *
      sqrt(1-phat)/sqrt(myfreq - (0.25 * phat + 0.47)*sqrt(myfreq))
    RVAL <- cbind(count, freq, NA, NA, NA, NA, NA)
    RVAL[myindex,3:7] <- cbind(y,yhat,ci.width, yhat-ci.width, yhat + ci.width)
    RVAL <- as.data.frame(RVAL)
    names(RVAL) <- c("Counts", "Freq", "Metameter", "CI.center",
                     "CI.width", "CI.lower", "CI.upper")
    x_lim <- range(RVAL[,1])
    y_line <- predict(fm, newdata = data.frame(mycount = xlim))

    RVAL$y_line <- y_line

    if(plot){
      switch(match.arg(type),
             "poisson" = {
               if(is.null(title)){title="Poisonness Plot"}
               label_x = 0.9*quantile(RVAL$Count)[4]
               label_y = quantile(RVAL$Metameter)[2]
             },
             "binomial" = {
               if(is.null(title)){title="Binomialness Plot"}
               label_x = 0.9*quantile(RVAL$Count)[4]
               label_y = quantile(RVAL$Metameter)[2]
             },
             "nbinomial" = {
               if(is.null(title)){title="Negative Binomialness Plot"}
               label_x = 0.8*quantile(RVAL$Count)[2]
               label_y = quantile(RVAL$Metameter)[2]
             }
      )

      legend.text <- c(paste("slope =", round(coef(fm)[2], digits = 3)),
                       paste("intercept =", round(coef(fm)[1], digits = 3)),
                       "", paste(names(par.estim),": ML =", round(par.ml, digits=3)),
                       legend.text)
      legend.text <- paste(legend.text, collapse = "\n")

      p <- ggplot(RVAL,aes()) +
        geom_line(aes(x=Counts,y=y_line, color="Perfect distribution"),
                  size=0.75, key_glyph="point") +
        geom_point(aes(x=Counts, y=Metameter, color="Observed distribution"),
                   key_glyph="point") +
        xlab(xlab) +
        ylab(ylab) +
        ggtitle(title) +
        scale_colour_manual(values=c("red", "blue")) +
        theme_bw() +
        theme(legend.position="none")

      if(conf_int){ p <- p + geom_errorbar(aes(ymin=CI.lower,
                                         ymax=CI.upper,
                                         x=Counts,
                                         y=Metameter), color="red")}
      if(label){p <- p + annotate(geom="text", x=label_x, y=label_y,
                                   label=legend.text,hjust=0)}
      if(legend){p <- p + labs(color="", position="bottom") +
                          theme(legend.position="top")}
      return(p)
    }
    return(RVAL)
  }
