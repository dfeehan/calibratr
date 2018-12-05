

## helper functions for calibration

coef.var <- function(x) { return(sd(x) / mean(x)) }


## wrapper around the survey package implementation of calibration
##
## this assumes vars have the same name in the census and in
## the survey (ie, recodings have happened the same way)
##
## fml - the formula describing how to calibrate
## calfun - see survey package ?calibrate for more info
## base - the base survey design
## cdat - the census dataset with microdata
## caggdat - the census dataset with aggregate data
## ... - other args to pass to the calibrate function. see ?calibrate from survey package
##
## returns: the vector of calibrated weights
calib_weights <- function(fml, fml.agg=NULL, calfun="linear", base.design, cdat, caggdat=NULL, ...) {

  pop.tots.mat <- model.matrix(fml, data=cdat)

  pop.tots <- aaply(pop.tots.mat,
                    2,
                    function(col, weights=cdat$ibge_weight) {
                      return(sum(col*weights))
                    })

  ## the user can optionally specify aggregate-level quantites
  if (! is.null(fml.agg)) {
    pop.tots.agg.mat <- model.matrix(fml.agg, data=caggdat)

    # calculate totals for aggregate quantities
    pop.agg.tots <- aaply(pop.tots.agg.mat,
                          2,
                          # use the number of adults in the census block as the weight
                          function(col, weights=caggdat$tot_adult_all) {
                            return(sum(col*weights))
                          })

    all.tots <- c(pop.tots, pop.agg.tots[-1])

    calib.fml <- formula(as.Formula(terms(fml), formula(Formula(terms(fml.agg)), lhs=0)), collapse=TRUE)

    thisorder <- colnames(model.matrix(calib.fml, base.design$variables))

    # be sure that the order of the totals matches what is expected by calibrate()
    all.tots <- all.tots[thisorder]

  } else {
    all.tots <- pop.tots
    calib.fml <- fml
  }

  res <- calibrate(design = base.design,
                   formula = calib.fml,
                   calfun = calfun,
                   population = all.tots,
                   ...)

  cat("calibrated using : ", paste(calib.fml), "\n")
  cat("population totals: ", "\n")
  for( i in names(all.tots)) {
    cat(i, " : ", all.tots[i], "\n")
  }

  return(weights(res))
}

weighted_estimate <- function(df, fml, weight) {

  tots.mat <- model.matrix(fml, data=df)

  weighted.tots <- aaply(tots.mat,
                         2,
                         function(col, weights=df[,weight]) {
                           return(sum(col*weights))
                         })
  return(weighted.tots)

}

## data - the survey dataset
## fml - formula describing estimates
## weight.vec - vector of weight names, eg: c("weight1", "weight2")
## cdat - the census microdataset
compare_estimates <- function(df, fml, weight.vec, cdat) {

  sres <- ldply(weight.vec,
                function(thisw) {
                  res <- weighted_estimate(df, fml, thisw)
                  res <- data.frame(value=names(res),
                                    tot=res,
                                    wname=thisw,
                                    stringsAsFactors=FALSE)
                  return(res)
                })

  ## also calculate the total for each type of weights,
  ## so we can get the proportional cell sizes
  tots <- ldply(weight.vec,
                function(thisw) {
                  res <- data.frame(wgt_tot=sum(df[,thisw]),
                                    wname=thisw,
                                    stringsAsFactors=FALSE)
                })

  tots <- rbind(tots,
                data.frame(wgt_tot=sum(cdat[,'ibge_weight']),
                           wname='ibge_weight',
                           stringsAsFactors=FALSE))

  cres <- weighted_estimate(cdat, fml, 'ibge_weight')
  cres <- data.frame(value=names(cres), tot=cres, wname='ibge_weight')

  all.res <- rbind(sres, cres)

  all.res <- all.res %>% left_join(tots, by='wname') %>%
    mutate(frac = tot / wgt_tot)

  return(all.res)

}

plot_compare_estimates <- function(data, title="", xlab="", intercept=FALSE) {

  if(! intercept) {
    data <- data %>% filter(value != "(Intercept)")
  }

  ggplot(data) +
    geom_bar(aes(x=value,
                 y=tot,
                 fill=wname),
             stat='identity',
             position='dodge') +
    ggtitle(title) +
    xlab(xlab)

}

plot_compare_estimates_prop <- function(data, title="", xlab="", intercept=FALSE) {

  if(! intercept) {
    data <- data %>% filter(value != "(Intercept)")
  }

  ggplot(data) +
    geom_bar(aes(x=value,
                 y=frac,
                 fill=wname),
             stat='identity',
             position='dodge') +
    ggtitle(title) +
    xlab(xlab)

}

plot_weight_distn <- function(data, weight.vec, binwidth=NULL) {

  wdat <- ldply(weight.vec,
                function(thisw) {
                  return(setNames(data.frame(name=thisw, value=data[,thisw]),
                                  c('name', 'value')))

                })

  if(is.null(binwidth)) {
    res <- ggplot(wdat) +
      geom_vline(xintercept=0, lty=2) +
      geom_density(aes(x=value, color=name))
  } else {
    res <- ggplot(wdat) +
      geom_vline(xintercept=0, lty=2) +
      geom_histogram(aes(x=value, color=name), binwidth=binwidth)
  }

  return(res)
}

do_calibrations <- function(states,
                            base.dat,
                            census.microdat,
                            census.agg,
                            uber.calfun='linear',
                            ## hack / special case: these two states
                            ## can't be calibrated on irregular ('subnormal')
                            ## census blocks b/c the cities have none or we
                            ## have none in our sample
                            states.without.cb_subnorm = c('GO', 'TO'),
                            uber.bounds=c(-Inf, Inf)) {

  res <- ldply(states,
               function(this.state) {

                 cat("================================================\n")
                 cat("starting ", this.state, "\n")
                 cat("================================================\n")

                 state.base.dat <- base.dat %>%
                   filter(state_abbrev == this.state)

                 state.census.dat <- census.microdat %>%
                   filter(state_abbrev == this.state)

                 state.aggcensus.dat <- census.agg %>%
                   filter(state_abbrev == this.state)

                 ##########################
                 ## base weights
                 basedesign <- svydesign(id = ~ cod_setor_censitario + id_endereco,
                                         weights = ~ base_weight,
                                         data=state.base.dat)

                 worig <- weights(basedesign)

                 ##########################
                 ## tst_weight0
                 ## scale base weights up to city popn total to acount for nonresponse
                 state.base.dat$tst_weight00 <- calib_weights(~ 1,
                                                              base=basedesign,
                                                              calfun=uber.calfun,
                                                              bounds=uber.bounds,
                                                              cdat=state.census.dat)

                 ##########################
                 ## tst_weight1
                 ## age X sex
                 state.base.dat$tst_weight01 <- calib_weights(~ male*factor(agegp10v2),
                                                              base=basedesign,
                                                              calfun=uber.calfun,
                                                              bounds=uber.bounds,
                                                              cdat=state.census.dat)

                 ##########################
                 ## tst_weight2
                 ## age X sex + cb_subnorm
                 if (this.state %in% states.without.cb_subnorm) {

                   ## without cb_subnorm, this is the same as tst_weight1
                   state.base.dat$tst_weight02 <- state.base.dat$tst_weight01

                 } else {

                   state.base.dat$tst_weight02 <- calib_weights(fml = ~ male*factor(agegp10v2),
                                                                fml.agg = ~ factor(cb_subnorm),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)
                 }

                 ##########################
                 ## tst_weight3
                 ## age X sex + educ
                 state.base.dat$tst_weight03 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                factor(educ_recode2),
                                                              #fml.agg = ~ factor(cb_subnorm),
                                                              base=basedesign,
                                                              calfun=uber.calfun,
                                                              bounds=uber.bounds,
                                                              cdat=state.census.dat,
                                                              caggdat=state.aggcensus.dat)

                 ##########################
                 ## tst_weight4
                 ## age X sex + hhsize
                 state.base.dat$tst_weight04 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                factor(hh_num_elig_gp3),
                                                              base=basedesign,
                                                              calfun=uber.calfun,
                                                              bounds=uber.bounds,
                                                              cdat=state.census.dat,
                                                              caggdat=state.aggcensus.dat)

                 ##########################
                 ## tst_weight5
                 ## age X sex + marital
                 state.base.dat$tst_weight05 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                factor(marital2),
                                                              #fml.agg = ~ factor(cb_subnorm),
                                                              base=basedesign,
                                                              calfun=uber.calfun,
                                                              bounds=uber.bounds,
                                                              cdat=state.census.dat,
                                                              caggdat=state.aggcensus.dat)

                 ##########################
                 ## tst_weight6
                 ## age X sex + educ + favela
                 if (this.state %in% states.without.cb_subnorm) {

                   ## without cb_subnorm, this is the same as tst_weight3
                   state.base.dat$tst_weight06 <- state.base.dat$tst_weight03

                 } else {
                   state.base.dat$tst_weight06 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(educ_recode2),
                                                                fml.agg = ~ factor(cb_subnorm),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)
                 }

                 ##########################
                 ## tst_weight7
                 ## age X sex + hhsize + favela
                 if (this.state %in% states.without.cb_subnorm) {

                   ## without cb_subnorm, this is the same as tst_weight4
                   state.base.dat$tst_weight07 <- state.base.dat$tst_weight04

                 } else {
                   state.base.dat$tst_weight07 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(hh_num_elig_gp3),
                                                                fml.agg = ~ factor(cb_subnorm),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)
                 }

                 ##########################
                 ## tst_weight8
                 ## age X sex + marital + favela
                 if (this.state %in% states.without.cb_subnorm) {

                   ## without cb_subnorm, this is the same as tst_weight5
                   state.base.dat$tst_weight08 <- state.base.dat$tst_weight05

                 } else {
                   state.base.dat$tst_weight08 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(marital2),
                                                                fml.agg = ~ factor(cb_subnorm),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)
                 }

                 ##########################
                 ## tst_weight9
                 ## age X sex + marital + educ + favela
                 if (this.state %in% states.without.cb_subnorm) {

                   ## without cb_subnorm, this is ageXsex + marital + educ
                   state.base.dat$tst_weight09 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(marital2) +
                                                                  factor(educ_recode2),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)

                 } else {
                   state.base.dat$tst_weight09 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(marital2) +
                                                                  factor(educ_recode2),
                                                                fml.agg = ~ factor(cb_subnorm),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)
                 }

                 ##########################
                 ## tst_weight10
                 ## age X sex + hhsize(orig) + favela
                 ##
                 ## ie, use original hh size categories instead of the coarser ones
                 if (this.state %in% states.without.cb_subnorm) {

                   state.base.dat$tst_weight10 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(hh_num_elig_gp2),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)

                 } else {
                   state.base.dat$tst_weight10 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(hh_num_elig_gp2),
                                                                fml.agg = ~ factor(cb_subnorm),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)
                 }

                 ##########################
                 ## tst_weight11
                 ## age X sex + hhsize(new) + favela
                 ##
                 ## ie, use new hh size categories instead of the coarser ones
                 if (this.state %in% states.without.cb_subnorm) {

                   state.base.dat$tst_weight11 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(hh_num_elig_gp4),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)

                 } else {
                   state.base.dat$tst_weight11 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(hh_num_elig_gp4),
                                                                fml.agg = ~ factor(cb_subnorm),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)
                 }

                 ##########################
                 ## tst_weight11
                 ## age X sex + hhsize(new) + favela
                 ##
                 ## ie, use new hh size categories instead of the coarser ones
                 if (this.state %in% states.without.cb_subnorm) {

                   state.base.dat$tst_weight12 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(hh_single),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)

                 } else {
                   state.base.dat$tst_weight12 <- calib_weights(fml = ~ male*factor(agegp10v2) +
                                                                  factor(hh_single),
                                                                fml.agg = ~ factor(cb_subnorm),
                                                                base=basedesign,
                                                                calfun=uber.calfun,
                                                                bounds=uber.bounds,
                                                                cdat=state.census.dat,
                                                                caggdat=state.aggcensus.dat)
                 }



                 return(state.base.dat)
               })






  return(res)

}

summarize_calibrations <- function(calib.res) {

  all.states <- unique(calib.res$state_abbrev)
  all.weights <- c('base_weight',
                   colnames(calib.res %>% select(starts_with('tst'))))

  all.state.wgts <- expand.grid(state_abbrev=all.states,
                                weight = all.weights)

  wgt.summ <- ddply(all.state.wgts,
                    .(state_abbrev, weight),
                    #function(this.wgt, df=state.base.dat) {
                    function(row) {

                      this.wgt <- paste(row$weight)
                      this.state <- paste(row$state_abbrev)

                      df <- calib.res %>% filter(state_abbrev == this.state)

                      this.w <- as.data.frame(df)[,this.wgt]

                      res <- data.frame(state_abbrev = this.state,
                                        weight = this.wgt,
                                        mean = mean(this.w),
                                        min = min(this.w),
                                        max = max(this.w),
                                        median = median(this.w),
                                        coef.var = coef.var(this.w))
                      return(res)
                    })

  wgt.summ  <- wgt.summ %>%
    arrange(weight) %>%
    group_by(state_abbrev) %>%
    mutate(cv.base = first(coef.var))

  wgt.summ.summ <- wgt.summ %>%
    group_by(weight) %>%
    summarise(frac_nonpos = mean(min <= 0),
              mean.cv = mean(coef.var),
              mean.cv.ratio = mean(coef.var / cv.base))

  return(list(wgt.distn=wgt.summ, wgt.summ=wgt.summ.summ))

}

