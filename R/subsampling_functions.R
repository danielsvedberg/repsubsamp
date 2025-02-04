library(tidyverse) #library important for data-wrangling
library(lubridate) #library important for managing date-time data
library(stringr)
library(fitdistrplus)
library(philentropy)
#library(bettermc)
library(MASS)
library(minpack.lm)
library(parallel)


#' get multinom params
#'
#' This function gets the multinomial distribution parameters for each unique
#' category in col
#' @param col is a dataframe column you want to get parameters for
#' @return params, a dataframe with a column containing each category in col, and the rate of each category in the rate column
#' @keywords params
#' @export
#' @examples
#' get_multinom_params()
get_multinom_params = function(col){
  len = length(col)
  df = data.frame(x = col)
  df = df %>%
    group_by(x) %>%
    summarise(n = n()) %>%
    mutate(rate = n/len)
  params = setNames(df$rate,df$x)
  return(params)
}

#' get_grouped_norm_params
#'
#' This function returns the normal distribution for each unique

get_grouped_norm_params = function(col){
  len = length(col)
  d = data.frame(x = col)
  d = d %>%
    group_by(x) %>%
    summarise(n = n()) %>%
    mutate(prop = n/len)
  avg = mean(d$prop, na.rm = TRUE)
  sd = sd(d$prop, na.rm = TRUE)
  params = c('mean' = avg, 'sd' = sd)
  return(params)
}

get_binom_params = function(col){
  k = sum(col==TRUE, na.rm=TRUE)
  n = length(col)
  p = k/n
  params = c('k' = k,'n' = n,'p'= p)
  return(params)
}

get_exp_params = function(col){
  lam = 1/mean(col)
  params = c('lambda'=lam)
  return(params)
}

get_norm_params = function(col){
  mean = mean(col, na.rm=TRUE)
  sd = sd(col, na.rm=TRUE)
  params = c('mean' = mean,'sd' = sd)
  return(params)
}

get_gamma_params = function(col){
  gamfit = try(fitdist(col, distr = "gamma", method = 'mle'))
  if(!inherits(gamfit,"try-error")){
    params = gamfit$estimate
  } else {
    params = c(NA,NA)
  }
  return(params)
}

get_poiss_params = function(col){
  lambda = mean(col, na.rm = TRUE)
  params = c('lambda' = lambda)
  return(params)
}

#get LR returns Likelihood Ratio of L(fullcol)/l(fullcol), where L is the PDF built on fullcol data, and l is built on subcol data
getLikeParams= function(col, type){
  prmsfuns = list(poiss = get_poiss_params,
                  groupednorm = get_grouped_norm_params,
                  gamma = get_gamma_params,
                  binom = get_binom_params,
                  norm = get_norm_params,
                  multinom = get_multinom_params,
                  exp = get_exp_params)
  fun = prmsfuns[[type]]
  prms = fun(col)
  return(prms)
}

#dmultinom(c(1,0,0), prob = c(0,0.25,0.25,0.5)) you can get the likelihood of each individual event like this
dmultinomLogLike = function(col, params){
  idx = names(unlist(params))
  df = tibble(idx = idx, n = unlist(map(idx, function(x) sum(col %in% x))))
  like = dmultinom(df$n, prob = params,log=TRUE)
  return(like)
}

dgammaLogLike = function(col, params){
  like = dgamma(col,params[1],params[2], log = TRUE)
  return(like)
}

dbinomLogLike = function(col, params){
  n_succ = unlist(sum(col))
  n = length(col)
  like = dbinom(n_succ, n, prob = params[3], log = TRUE)
  return(like)
}

dgroupnormLogLike = function(col, params){
  len = length(col)
  d = data.frame(x = col)
  d = d %>%
    group_by(x) %>%
    summarise(prop = n()/len)
  like = dnorm(d$prop, params[1],params[2], log = TRUE)
  return(like)
}

dnormLogLike = function(col,params){
  like = dnorm(col,params[1],params[2], log = TRUE)
  return(like)
}

dexpLogLike = function(col, params){
  like = dexp(col, params[1], log = TRUE)
  return(liike)
}

dpoissLogLike = function(col, params){
  like = sum(log((params[1]^col * exp(-params[1])) / factorial(col)))
  return(like)
}

getLogLike = function(col,params,type){
  if(all(!is.na(params))){
    dnormlog = function() dnormLogLike(col, params)
    dmultinomlog = function() dmultinomLogLike(col,params)
    dgamlog = function() dgammaLogLike(col, params)
    dbinomlog = function() dbinomLogLike(col, params)
    dpoisslog = function() dpoissLogLike(col, params)
    dgroupnormlog = function() dgroupnormLogLike(col, params)
    dexplog = function() dexpLogLike(col, params)
    lffuns = list(exp = dexplog, gamma = dgamlog, binom = dbinomlog, norm = dnormlog, multinom = dmultinomlog, groupednorm = dgroupnormlog, poiss = dpoisslog)
    like = lffuns[[type]]()
    return(like)
  } else{
    return(NA)
  }
}

normalize = function(vec){
  tot = sum(vec, na.rm=TRUE)
  normed = vec/tot
  return(normed)
}

dgroupedNormlogProb = function(col,params){
  len = length(col)
  d = data.frame(x = col)
  d = d %>%
    group_by(x) %>%
    summarise(prop = n()/len)
  like = dnorm(unique(d$prop), params[1],params[2])
  prob = log(normalize(like))
  return(prob)
}

dnormlogProb = function(col,params){
  like = dnorm(unique(col), params[1], params[2])
  #if(!inherits(like,"try-error")){
  prob = log(normalize(like))
  return(prob)
}

dgammalogProb = function(col, params){
  like = dgamma(unique(col),params[1],params[2])
  # if(!inherits(like,"try-error")){
  prob = log(normalize(like))
  return(prob)
  # } else {
  #   return(NA)
  # }
}

dbinomlogProb = function(params){
  #if(all(!is.na(params))){
    prob0 = log(1-params[3])
    prob1 = log(params[3])
    prob = c(prob1, prob0)
    return(prob)
  # } else {
  #   return(c(NA,NA))
  # }
}

dexplogProb = function(col, params){
  prob = dexp(col, rate = params[1])
  return(prob)
}

dpoisslogProb = function(col, params){
  prob = dpois(col, lambda = params[1])
  return(prob)
}

getProb = function(col,params,type){
  if(all(!is.na(params))){
    normProb = function() dnormlogProb(col, params)
    groupnormProb = function() dgroupedNormlogProb(col, params)
    gamProb = function() dgammalogProb(col, params)
    multinomProb = function() return(log(params))
    binomProb = function() dbinomlogProb(params)
    poissProb = function() dpoisslogProb(col, params)
    expProb = function() dexplogProb(col, params)
    lffuns = list(exp = expProb,
                  norm = normProb,
                  groupednorm = groupnormProb,
                  gamma = gamProb,
                  multinom = multinomProb,
                  binom = binomProb,
                  poiss = poissProb
                  )
    prob = lffuns[[type]]()
    return(prob)
  } else {
    return(NA)
  }
}

getPrmsDf = function(df,colbinds,flag = NA){
  cols = names(colbinds)
  distnms = unname(colbinds)
  metadf = data.frame(cols = cols, distnms = distnms)
  metadf = metadf %>%
    group_by(cols) %>%
    mutate(prms = list(getLikeParams(unlist(df[cols]),distnms)))
  if(!is.na(flag)){
    newname1 = paste(flag,'prms',sep="")
    metadf = metadf %>%
      rename(!!newname1 := prms)
  }
  return(metadf)
}

logkld = function(p,q){
  kld = try(sum((exp(p)*(p-q)),na.rm=TRUE))
  if(!inherits(kld, 'try-error')){
    return(kld)
  } else {
    return(NA)
  }
}

logLLR = function(a,b){
  LLR = try(sum(a,na.rm=TRUE)-sum(b,na.rm=TRUE))
  if(!inherits(LLR,'try-error')){
    return(LLR)
  } else {
    return(NA)
  }
}

bootStrapRate = function(idcol, nIter, intervals = c(0.975,0.5,0.025)){
  nSubs = length(idcol)
  as = function() sum(sample(idcol, replace = TRUE))/nSubs
  rates = replicate(nIter, as())
  rates = quantile(rates, intervals)
  df = data_frame(tiles = c("upper","median","lower"), values = rates)
  return(df)
}

permTest = function(col1,col2, nIter){
  nS1 = length(col1)
  nS2 = length(col2)
  all = c(col1,col2)
  rtdiff = function() {
    samps = sample(all, replace = FALSE)
    r1 = sum(samps[1:nS1])/nS1
    r2 = sum(samps[-nS2:-1])/nS2
    #r1 = sum(sample(all, nS1))/nS1
    #r2 = sum(sample(all, nS2))/nS2
    diff = r1-r2
    return(diff)
  }
  diffs = replicate(nIter, rtdiff())
  #diffs = quantile(diffs, c(0.95,0.5,0.05))
  return(diffs)
}

ORpermTest = function(col1, col2, nIter){
  nS1 = length(col1)
  nS2 = length(col2)
  all = c(col1,col2)
  permOR = function() {
    samps = sample(all, replace = FALSE)
    SS1 = samps[1:nS1]
    SS2 = samps[-nS1:-1]
    pos1 = sum(SS1)
    neg1 = length(SS1) - pos1
    pos2 = sum(SS2)
    neg2 = length(SS2) - pos2
    r1 = pos1/neg1
    r2 = pos2/neg2
    OR = r1/r2
    return(OR)
  }
  ORs = replicate(nIter, permOR())
  return(ORs)
}

RRpermTest = function(col1, col2, nIter){
  nS1 = length(col1)
  nS2 = length(col2)
  all = c(col1,col2)
  permRR = function() {
    samps = sample(all, replace = FALSE)
    r1 = sum(samps[1:nS1])/nS1
    r2 = sum(samps[-nS1:-1])/nS2
    RR = r1/r2
    return(RR)
  }
  RRs = replicate(nIter, permRR())
  return(RRs)
}

calculateRelativeRisk <- function(df, eventColumn, conditioningColumn){
  # Calculate the contingency table
  contingencyTable <- table(df[[eventColumn]], df[[conditioningColumn]])
  # Calculate the relative risk
  RR <- (contingencyTable[2, 2] / sum(contingencyTable[2, ])) / (contingencyTable[1, 2] / sum(contingencyTable[1, ]))
  return(RR)
}
relativeRiskPermutationTest <- function(df, eventColumn, conditioningColumn, numPermutations = 1000){
  # Calculate the observed relative risk
  observedRR <- calculateRelativeRisk(df, eventColumn, conditioningColumn)
  # Initialize a variable to hold the number of permutations with a higher relative risk
  permRR = function() {
    # Perform a permutation on the event column
    permutedEvents <- sample(df[[eventColumn]])
    # Create a temporary dataframe with the permuted events
    permutedDf <- df
    permutedDf[[eventColumn]] <- permutedEvents
    # Calculate the relative risk for the permuted data
    permutedRR <- calculateRelativeRisk(permutedDf, eventColumn, conditioningColumn)
    return(permutedRR)
  }
  RRs = replicate(permRR, numPermutations)

  return(RRs)
}


generateCandidateSubSamples = function(idcol, nSubs, nIter){
  ss = function() sample(idcol, nSubs)
  candidates = replicate(nIter, ss())
  return(candidates)
}

generateAutoSamples = function(idcol, nIter){
  nSubs = length(idcol)
  as = function() sample(idcol, replace = TRUE)
  candidates = replicate(nIter, as())
  return(candidates)
}

jointEval = function(df, joinmetadf){
  resdf = joinmetadf %>%
    group_by(cols) %>%
    mutate(fGf_LL = list(getLogLike(unlist(df[cols]),unlist(prms),distnms)),
           fGs_LL = list(getLogLike(unlist(df[cols]),unlist(subprms),distnms)),
           fGf_prob = list(getProb(unlist(df[cols]),unlist(prms),distnms)),
           fGs_prob = list(getProb(unlist(df[cols]),unlist(subprms),distnms)),
           f_LLR = logLLR(unlist(fGs_LL),unlist(fGf_LL)), #https://books.google.com/books?id=cligOwrd7XoC&pg=PA84#v=onepage&q&f=false
           f_KLD = logkld(unlist(fGs_prob),unlist(fGf_prob))) %>%
    ungroup() %>%
    bind_rows(summarise(.,
                        across(where(is.numeric), sum),
                        across(where(is.character), ~"Total")))

  return(resdf)
}

evalSingleSub = function(df,subIds,idCol,colbinds){
  subdf = df %>%
    group_by(!!sym(idCol)) %>%
    filter(!!sym(idCol) %in% subIds)
  metadf = getPrmsDf(df, colbinds)
  submetadf = getPrmsDf(subdf,colbinds,flag = 'sub')
  metadf = merge(metadf,submetadf,by = c("cols","distnms"))
  res = jointEval(df,metadf)
  return(res)
}

#' evalsub evaluates the similarity between a subsample and the population
#'
#' @param df is the population dataframe
#' @param subIds is a list of the IDs in idCol of df that compose the subsample
#' @param idCol is a string of the column name for the id column
#' @param metadf is the metadata dataframe, with columns:
#' cols: a column of strings naming the columns used to model the population
#' distnms: a column of strings naming the distributions used to model the variable specified by the corresponding entry in cols
#' prms: a column of lists, where each list contains the parameters for the population data for the variable specified by the row

evalSub = function(df,subIds,idCol,metadf){
  subIds = unlist(subIds)
  sub = df %>%
    group_by(!!sym(idCol)) %>%
    filter(!!sym(idCol) %in% subIds)
  print(sub)
  colbinds = setNames(metadf$distnms,metadf$cols)
  submetadf = getPrmsDf(sub,colbinds,flag = 'sub') #creation of a metadf-equivalent for the subsample, where column 'prms' becomes 'subprms'
  joinmetadf = merge(metadf,submetadf, by = c("cols","distnms"))
  joinmetadf = jointEval(df,joinmetadf)
  KLD = tail(joinmetadf$f_KLD,1)
  LLR = tail(joinmetadf$f_LLR,1)
  res = c('LLR' = LLR, 'KLD' = KLD)
  return(res)
}


#' multiEval: evaluate many samples against a population
#'
#' This function evaluates the similarity between a population of data and many candidate sub/re-samples of rows in the population
#' Similarity between the two is measured by the log-likelihood that the sample resembles the population data,
#' according to the variables contained in columns of df, specified in colbinds
#'
#' @param df is the population dataframe. df needs to have a column matching idCol, as well as a column for each column named in colbinds
#' @param candidates is a matrix of ids for the candidate subsamples. Each column is an iteration, each row is an ID in the subsample
#' @param idCol is the string for the ID column that identifies each row
#' @param colbinds is a named list, where each item is a the name of a column in df,
#' and each corresponding name is the string-identifier for the distribution used to model the variable in the column
#' distributions and their identifiers:
#' -multinomial: 'multinom'
#' -binomial: 'binom'
#' -normal: 'norm'
#' -gamma: 'gamma'
#' @param mp is a boolean argument to specify if you wnat to parallelize execution across subsamples
#' @returns outdf a dataframe where each row is an iteration in candidates,
#' with the following columns:
#' -iterno: the iteration id
#' -ids: a list of the row ids in the sample
#' -LLR: the log-likelihood ratio; the (log) likelihood of the population (LoP) data given the subsample's paramter estimates, divided by LoP given the population's parameter estimates
#' -KLD: Kullback-leibler divergence between the sample's distributions and the population's distributions
#' @seealso [generateAutoSamples()], [generateCandidateSubSamples()], which return matrices that can be input into candidates
#' @examples multiEval(population_df,candidate_matrtix,'IDs', colbinds)
multiEval = function(df, candidates, idCol, colbinds, mp = TRUE){
  fullmetadf = getPrmsDf(df,colbinds)
  outdf = as.data.frame(t(candidates)) %>%
    mutate(iterno = 1) %>%
    mutate(iterno = cumsum(iterno)) %>%
    pivot_longer(!iterno, names_to = NULL, values_to = 'id') %>%
    group_by(iterno) %>%
    summarise(id = list(id)) %>%
    rename(!!idCol := id)

  subIds = outdf[[idCol]]
  print(subIds)
  fun = function(x) evalSub(df,x,idCol,fullmetadf)
  if(mp == TRUE){
    ncores = parallel::detectCores() - 2
    print(ncores)
    res = parallel::mclapply(subIds,fun, mc.cores=ncores)#, mc.progress = TRUE)
  } else {
    res = lapply(subIds, fun)
  }
  res = bind_rows(res)
  outdf = cbind(outdf, res)
  return(outdf)
}

###code for subsample size optimization
itersampsize <- function(df, idCol, colbinds, nIter, min_samps, max_samps){
  ncores = parallel::detectCores()/2
  nSamps = seq(min_samps, max_samps, 10) #make a vector spanning min and max samps, skipping every x
  fun = function(nS){
    candidates = generateCandidateSubSamples(df[[idCol]],nS,nIter)
    res = multiEval(df, candidates, idCol, colbinds, mp=FALSE)
    res$nSamps = nS
    message(paste("testing ", nS, " samples"))
    return(res)
  }
  fullres = parallel::mclapply(nSamps, fun, mc.cores = ncores)#, mc.progress = TRUE)
  fullres = bind_rows(fullres)
  return(fullres)
}

################################################################################
#OLD SHIT BELOW#################################################################
################################################################################

get_burg_stats = function(data) {
  dataStat = data %>%
    ungroup() %>%
    summarise(nCases = n(),
              nDocsLam = 1/mean(n_docs),
              kUnanimous = sum(unanimous==TRUE, na.rm = TRUE),
              pUnanimous = kUnanimous/nCases,
              knytSal = sum(nytSalience==TRUE, na.rm=TRUE),
              pnytSal = knytSal/nCases,
              kcqSal = sum(cqSalience==TRUE, na.rm=TRUE),
              pcqSal = kcqSal/nCases
    )
  return(dataStat)
}

#get_old_case_nos returns list of unique case numbers (usCite format) from csv where there exists a usCite column
get_old_case_nos = function(fn){
  old_subset = read_csv(fn)
  old_subset = old_subset %>%
    filter(!is.na(usCite)) %>%
    dplyr::select(usCite) %>%
    distinct()
  old_cases = old_subset$usCite
  return(old_cases)
}

#returns the best subsample for cases containing a specific action code (subsetCode)
#with size of nSamps from nIter iterations
makeSubsetFrame = function(burgerDocuments, subcases, ac, label){
  subFrame = burgerDocuments %>%
    group_by(id) %>%
    filter(ac %in% unlist(action_code)) %>%
    filter(ids %in% unlist(subcases)) %>%
    distinct() %>%
    mutate(justice_name = paste(unlist(justice_name)),
           caseName = paste(unlist(caseName),collapse = ";"),
           action_code = paste(action_code)) %>%
    distinct()

  filename = paste("action",ac,label,"subset.csv",sep="")
  write_csv(subFrame, filename)
  return(subFrame)
}

get_best_kls = function(kls){
  bestkls = kls %>%
    filter(!is.na(KLD)) %>%
    group_by(n_samps) %>%
    summarise(bestKLD = min(KLD))
  return(bestkls)
}

model_samp_size = function(df, col){
  #fit exponential decay curve to the data for best model, then perform Kneedle algorithm to find point of maximum curvature
  #THANK YOU MENG XU on RPUBS
  col = as.symbol(col)

  model0 <- eval(bquote(lm(log(.(col)) ~ nSamps, data = df))) #lm(log(bestKLD - theta0) ~ n_samps, data=bestkls)
  alpha0 = exp(coef(model0)[1])
  beta0 = coef(model0)[2]
  start = list(alpha = alpha0, beta = beta0)#, theta = theta0)
  model = eval(bquote(nlsLM(.(col) ~ alpha * exp(beta * nSamps), data = df, start = start)))
  return(model)
}

expcurvature = function(x,model){
  coeffs = coef(model)
  alpha = coeffs[[1]]
  beta = coeffs[[2]]
  numerator = exp(x*beta)
  denomnum = (alpha^2*beta^2*exp(2*beta*x))+1
  denomdenomdenom = denomnum^3
  denomdenom = ((alpha^4*beta^6*exp(2*x*beta))/denomdenomdenom) + ((alpha^2*beta^4)/denomdenomdenom)
  denom = sqrt(denomnum/denomdenom)
  curve = numerator/denom
  return(curve)
}

get_max_curve = function(model){
  x  = 0:1000
  curve = expcurvature(x,model)
  max_curve = x[which.max(curve)]
  return(max_curve)
}

plot_kls_model = function(bestkls, model, ac, tag = NULL){
  n_samps = bestkls$n_samps
  y_pred = predict(model, new_data = bestkls)
  bestkls$y_pred = y_pred

  scale = max(bestkls$bestKLD) / max(bestkls$curvature)
  kldVsNbest = ggplot(data=bestkls)+
    geom_point(aes(x=n_samps, y=bestKLD))+
    geom_line(aes(x= n_samps, y = y_pred, color = "decay model of KLD"))+
    geom_line(aes(x = n_samps, y = curvature*scale , color = "curvature"))+
    scale_y_continuous(sec.axis = sec_axis(~./scale, name = "curvature"))+
    labs(y = "Kullback-Leibler Divergence",
         x = "number of subsamples") +
    scale_color_manual(values = c("dodgerblue1", "deeppink4"), name = "")+
    theme(legend.position = "bottom") +
    ggtitle(paste("KLD vs # of subsamples for action code",ac, tag))
  name = paste('KLDvsNbest',ac, tag, '.png', sep = "")
  ggsave(name, kldVsNbest,height = 5, width = 5)
  return(kldVsNbest)
}

