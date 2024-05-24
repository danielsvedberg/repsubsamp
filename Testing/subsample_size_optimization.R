library(tidyverse) #library important for data-wrangling
library(lubridate) #library important for managing date-time data
library(stringr)
library(fitdistrplus)
source("R/subsampling_functions.R")
source("Testing/example_preprocessing.R")

BCOWD400 = BCOWDCases %>%
  group_by(caseNo) %>%
  filter(400 %in% unlist(action_code))
cols = c('n_docs','term','nytSalience','cqSalience','unanimous')
distnms = c('poiss','multinom','binom','binom','binom')
colbinds = setNames(distnms, cols)

auto400cands = generateAutoSamples(BCOWD400$caseNo, 10000)
auto400 = multiEval(BCOWD400,auto400cands,'caseNo',colbinds)
auto400 = auto400 %>%
  mutate(LambdaLR = -2*LLR,
         LR = exp(LLR))

percentiles400LambdaLR = quantile(auto400$LambdaLR, c(0.05,0.5,0.95))
percentiles400KLD = quantile(auto400$KLD, c(0.05,0.5,0.95))
percentiles400LR = quantile(auto400$LR, c(0.05,0.5,0.95))

opti400 = itersampsize(BCOWD400, 'caseNo', colbinds, 1000, 20,600)
testcases = opti400$caseNo[1]
testres = evalSingleSub(BCOWD400,unlist(testcases),'caseNo',colbinds)
bestopti400 = opti400 %>%
  mutate(LR = exp(LLR),
         LambdaLR = -2*LLR,
         score = LR+KLD) %>%
  group_by(nSamps) %>%
  filter(LLR == max(LLR))
write_csv(bestopti400, 'bestopti400.csv')

opti400test = opti400 %>%
  mutate(LambdaLR = LLR*-2)
model400LambdaLR = model_samp_size(bestopti400, 'LambdaLR')
curve = get_max_curve(model400LambdaLR)

model400KLD= model_samp_size(bestopti400, 'KLD')
curve = get_max_curve(model400KLD)

BCOWD400sal = BCOWD400 %>%
  ungroup() %>%
  filter(nytSalience == TRUE)

auto40salcands = generateAutoSamples(BCOWD400sal$caseNo, 20000)
auto400sal = multiEval(BCOWD400sal, auto400salcands, 'caseNo', colbinds)
auto400sal = auto400sal %>%
  mutate(LambdaLR = -2*LLR,
         LR = exp(LLR))
pctiles400salKLD = quantile(auto400sal$KLD, c(0.05,0.5,0.95))
pctiles400salLambdaLR = quantile(auto400sal$LambdaLR, c(0.05,0.5,0.95))


BCOWD410 = BCOWDCases %>%
  filter(410 %in% unlist(action_code))
auto410cands = generateAutoSamples(BCOWD410$caseNo, 20000)
auto410 = multiEval(BCOWD410, auto410cands, 'caseNo', colbinds)
auto410 = auto410 %>%
  mutate(LambdaLR = -2*LLR,
         LR = exp(LLR))
pctiles410KLD = quantile(auto410$KLD, c(0.05,0.5,0.95))
pctiles410LambdaLR = quantile(auto410$LambdaLR, c(0.05,0.5,0.95))




BCOWD410 = BCOWDCases %>%
  group_by(caseNo) %>%
  filter(410 %in% unlist(action_code))
cols = c('n_docs','term','nytSalience','cqSalience','unanimous')
distnms = c('gamma','groupednorm','binom','binom','binom')
colbinds = setNames(distnms, cols)
opti410 = itersampsize(BCOWD410, 'caseNo', colbinds, 1000, 20,350)
testcases = opti410$caseNo[1]
testres = evalSingleSub(BCOWD410,unlist(testcases),'caseNo',colbinds)
bestopti410 = opti410 %>%
  mutate(LR = exp(LLR),
         LambdaLR = -2*LLR,
         score = LLRT+KLD) %>%
  group_by(nSamps) %>%
  filter(score == max(LLR))
write_csv(bestopti410, 'bestopti410.csv')
