

##' Help calibrate survey weights to population totals
##'
##' given the popn dataset, the survey dataset, and a formula describing
##' the calibration equation, produce calibration weights
##'
##' Note that if you get warnings when running this, that may be a sign that the
##' survey data have missingness on some of the calibration variables. Missingness
##' will also cause total calibrated weights to differ from total un-calibrated weights.
##'
##' @param fml formula describing the calibration design (NB: this relies on
##'   varnames and levels being the same in pdat and sdat)
##' @param pdat population-level data
##' @param pweight weights for adding up rows of population-level data
##' @param sdesign survey design object for survey dataset (see survey package)
##' @param calfun optional arg for survey package's calibrate fn
##'
##' @return a vector of calibrated weights
calibratr <- function(fml, pdat, pweight, sdesign, calfun="linear") {

  pop.tots <- model.matrix(fml, data=pdat)

  pop.tots <- plyr::aaply(pop.tots,
                    2,
                    function(col, weights=pdat[,pweight]) {
                      return(sum(col*weights))
                    })

  res <- survey::calibrate(design = sdesign,
                   formula = fml,
                   calfun = calfun,
                   population = pop.tots)

  return(weights(res))

}

##' get weighted covar cell totals for the given dataset,
##' covars, and weight var
agg_covars <- function(dat, covars, weight) {
  dat %>%
    group_by_(.dots=covars) %>%
    summarise_(count = lazyeval::interp(~sum(weight), weight=as.name(weight)))
}

##' compare several different weights
##'
##' Return a dataframe that has comparisons for different weights;
##' useful for checking different calibration specifications
##'
##' @param datalist a list of dataframes
##' @param weights  a vector w/ names of weights,
##' @param names    a vector w/ description of each data frame
##' @param covars   to covariates to aggregate by
##' @param title
##'
##' @return A dataframe with different weight comparisons
compare_covar_margins <- function(datalist,
                                  weights,
                                  names,
                                  covars,
                                  title) {

  agg_data <- plyr::llply(1:length(datalist),
                    function(idx) {
                      this_agg <- agg_covars(dat=datalist[[idx]],
                                             covars=covars,
                                             weight=weights[idx])
                      this_agg$source <- names[idx]
                      return(this_agg)
                    })

  res <- do.call(plyr::rbind.fill, agg_data)

  res <- res %>%
    group_by(source) %>%
    mutate(frac=count/sum(count)) %>%
    ungroup()

  return(res)

}

plot_margin_comparison <- function(margin.dat, covars, title='') {

  ggplot(margin.dat) +
    geom_bar(aes_string(x=covars, y='frac', fill='source'),
             position='dodge', stat='identity') +
    #scale_x_discrete(drop=FALSE) +
    ggtitle(title)

}

## wrapper for comparing unweighted survey respondent distn
## and population
plot_weightcomp_margins <- function(data,
                                    weights,   # vector like ('poststrat'='pswt', 'base'='basewt')
                                    pdat,
                                    pweight,
                                    covar,
                                    title="") {

  #cmp <- compare_covar_margins(datalist=list(data, data, pdat),
  datalist <- rep(list(survey), length(weights))
  datalist <- c(datalist, list(pop))

  if(is.null(names(weights))) {
    names(weights) <- paste(weights)
  }

  if(is.null(names(pweight))) {
    names(pweight) <- paste(pweight)
  }

  cmp <- compare_covar_margins(datalist=datalist,
                               weights=c(weights, pweight),
                               names=c(names(weights), names(pweight)),
                               covars=covar,
                               title=title)

  plot_margin_comparison(cmp, covar, title=title)

}
