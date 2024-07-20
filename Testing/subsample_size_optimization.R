library(tidyverse) #library important for data-wrangling
library(lubridate) #library important for managing date-time data
library(stringr)
library(fitdistrplus)

source("R/subsampling_functions.R") #load the sub sampling functions
source("Testing/example_preprocessing.R") #load and preprocess the test data
#critically, the pre-processing makes it so each row represents a unique case

folder = paste(getwd(), "Testing", sep="/")

################################################################################
#example 1: get the cases with code 400 and subsample them
################################################################################
#extract cases from BCOWD with code 400
BCOWD400 = BCOWDCases %>%
  group_by(caseNo) %>%
  filter(400 %in% unlist(action_code))

#set the column bindings by creating named list linking column names to distribution names
cols = c('n_docs','term','nytSalience','cqSalience','unanimous') #list of columns to be modeled
distnms = c('poiss','multinom','binom','binom','binom') #matching list of distribution labels for each column
colbinds = setNames(distnms, cols) #make a named list where each distribution label can be called by the column name

###generate & characterize the null distribution:
#generateAutoSamples resamples identifiers with replacement ("autosampling")
#takes in arguments of ids (list) and #iterations (int)
#returns a matrix of ids [#ids, #iterations]
auto400cands = generateAutoSamples(BCOWD400$caseNo, 1000)
#multiEval calculates the distance between the population and each sample in the matrix, according the columns specified in colbinds
auto400 = multiEval(BCOWD400,auto400cands,'caseNo',colbinds) #characterize the distance of each autosample to the population

###characterize null distribution:
#get the 95% confidence interval of the distances for:
percentiles400LLR = quantile(auto400$LLR, c(0.025,0.5,0.975)) #Log-likelihood ratio
percentiles400KLD = quantile(auto400$KLD, c(0.025,0.5,0.975)) #Kullback-Leibler divergence

###iteratively optimize sample size
#itersampsize takes the dataframe, id column name (str), column bindings (named list)
#number of repeats per iteration (int), minimum sub-sample size (int), and maximum sub-sample size (int)
opti400 = itersampsize(BCOWD400, 'caseNo', colbinds, 80, 100, 700)
thresh = percentiles400LLR['50%']

#plot a histogram of the LLRs in auto400
ggplot() +
  geom_histogram(data = auto400, aes(x=LLR)) +
  geom_vline(xintercept = thresh) +
  geom_label(aes(x=thresh, y=0, hjust='left', vjust='bottom', label=paste("50th %-tile:", round(thresh,2)))) +
  ylab("count") +
  xlab("Log-Likelihood Ratio") +
  expand_limits(x=0)+
  coord_flip() +
  ggtitle("null-distribution of LLR")
fnsvg = paste(folder, "LLR_null_dist.svg", sep="/")
fnpng = paste(folder, "LLR_null_dist.png", sep="/")
ggsave(fnsvg, width=6, height=4)
ggsave(fnpng, width=6, height=4)

#plot a histogram of the LLRs in auto400
ggplot() +
  geom_histogram(data = auto400, aes(x=exp(LLR))) +
  geom_vline(xintercept =exp(thresh)) +
  geom_label(aes(x=exp(thresh), y=0, hjust='left', vjust='bottom', label=paste("50th %-tile:", exp(thresh)))) +
  ylab("count") +
  xlab("Log-Likelihood Ratio") +
  expand_limits(x=0)+
  coord_flip() +
  ggtitle("null-distribution of LLR")

#plot a histogram
ggplot() +
  geom_histogram(data = auto400, aes(x=log(-LLR))) +
  geom_vline(xintercept = log(-thresh)) +
  geom_label(aes(x=log(-thresh), y=0, hjust='left', vjust='bottom', label=paste("50th %-tile:", round(log(-thresh),2)))) +
  ylab("count") +
  coord_flip() +
  ggtitle("null-distribution is log-normally distributed")
fnsvg = paste(folder, "LogLLR_null_dist.svg", sep="/")
fnpng = paste(folder, "LogLLR_null_dist.png", sep="/")
ggsave(fnsvg, width=6, height=4)
ggsave(fnpng, width=6, height=4)

#plot opti vs LLR
ggplot(data = opti400, aes(x=nSamps, y=LLR))+
  geom_point(alpha=0.1)+
  stat_summary(fun=mean, geom="line", color="blue") +
  geom_hline(yintercept=thresh) +
  geom_label(aes(y=-100, hjust='right', vjust='top', x=700, label=paste("Null-dist. 50th %-tile:", round(thresh,2)))) +
  ylab("log-likelihood ratio") +
  xlab("Sub-sample size") +
  ggtitle("Sub-sample size vs distance to population after 100 iterations")
fnsvg = paste(folder, "LLRvsNSamps.svg", sep="/")
fnpng = paste(folder, "LLRvsNSamps.png", sep="/")
ggsave(fnsvg, width=6, height=4)
ggsave(fnpng, width=6, height=4)


logthresh = round(log(-thresh),2)
ggplot(data = opti400, aes(x=nSamps, y=log(-LLR)))+
  geom_point(alpha=0.1)+
  stat_summary(fun=mean, geom="line", color="blue") +
  geom_smooth(method="lm") +
  geom_hline(yintercept=logthresh) +
  geom_label(aes(y=logthresh, hjust='left', vjust='bottom', x=100, label=paste("Null-dist. 50th %-tile:", logthresh))) +
  ylab("Log negative log-likelihood ratio") +
  xlab("Sub-sample size") +
  ggtitle("Sub-sample size vs distance to population after 100 iterations")
fnsvg = paste(folder, "LogLLRvsNSamps.svg", sep="/")
fnpng = paste(folder, "LogLLRvsNSamps.png", sep="/")
ggsave(fnsvg, width=6, height=4)
ggsave(fnpng, width=6, height=4)

expthresh = exp(thresh)
ggplot(data = opti400, aes(x=nSamps, y=exp(LLR)))+
  geom_point(alpha=0.1)+
  stat_summary(fun=mean, geom="line", color="blue") +
  geom_hline(yintercept=expthresh) +
  geom_label(aes(y=expthresh, hjust='left', vjust='top', x=100, label=paste("Null-dist. 50th %-tile:\n", round(expthresh,10)))) +
  ylim(0, expthresh)+
  ylab("likelihood ratio") +
  xlab("Sub-sample size") +
  ggtitle("Sub-sample size vs distance to population after 100 iterations")
fnsvg = paste(folder, "LRvsNSamps.svg", sep="/")
fnpng = paste(folder, "LRvsNSamps.png", sep="/")
ggsave(fnsvg, width=6, height=4)
ggsave(fnpng, width=6, height=4)


modelopti = opti400 %>%
  group_by(nSamps) %>%
  summarise(lmeanLLR = mean(log(-LLR), na.rm=TRUE),
            lsdLLR = sd(log(-LLR), na.rm=TRUE)) %>%
  mutate(prob_thresh = pnorm(log(-thresh), mean=lmeanLLR, sd=lsdLLR, lower.tail=TRUE),
         est_n_iter = 1/prob_thresh,
         LLR1k =
         LLR10k = qnorm(1/10000, mean=lmeanLLR, sd=lsdLLR),
         LLR100k = qnorm(1/100000, mean=lmeanLLR, sd=lsdLLR))


ggplot() +
  geom_line(data = modelopti, aes(x=nSamps, y=est_n_iter)) +
  ylim(0, 1e+04) +
  xlab("number of subsamples") +
  ylab("est. # iterations to achieve threshold")

#idea number 1: change itersampsize to iterate in a half-halving algorithm instead of in order, going up a half if the iteration is successful, down if it fails
#idea number 2: model each iteration using a gaussian, which can then be used to model the probability of a successful iteration in n-iterations


testcases = opti400$caseNo[1]
testres = evalSingleSub(BCOWD400,unlist(testcases),'caseNo',colbinds)

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
