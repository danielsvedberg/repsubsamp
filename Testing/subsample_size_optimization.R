library(tidyverse) #library important for data-wrangling
library(lubridate) #library important for managing date-time data
library(stringr)
library(fitdistrplus)

source("R/subsampling_functions.R") #load the sub sampling functions
folder = paste(getwd(), "Testing", sep="/")
BCOWDCases = read_csv("data_raw/BCOWDCases.csv")
#critically, the pre-processing makes it so each row represents a unique case

################################################################################
#example 1: get the cases with code 400 and subsample them
################################################################################
#extract cases from BCOWD with code 400
BCOWD400 = BCOWDCases %>%
  filter(str_detect(action_code, '400'))

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
auto400 = lists_to_text(auto400)

###characterize null distribution:
#get the 95% confidence interval of the distances for:
percentiles400LLR = quantile(auto400$LLR, c(0.025,0.5,0.975)) #Log-likelihood ratio
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

#plot a histogram of the log-normal distribution of the LLRs
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

###iteratively optimize sample size
#itersampsize takes the dataframe, id column name (str), column bindings (named list)
#number of repeats per iteration (int), minimum sub-sample size (int), and maximum sub-sample size (int)
opti400 = itersampsize(BCOWD400, 'caseNo', colbinds, 80, 100, 700)


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


modelopti = opti400 %>%
  group_by(nSamps) %>%
  summarise(lmeanLLR = mean(log(-LLR), na.rm=TRUE),
            lsdLLR = sd(log(-LLR), na.rm=TRUE)) %>%
  mutate(prob_thresh = pnorm(log(-thresh), mean=lmeanLLR, sd=lsdLLR, lower.tail=TRUE),
         est_n_iter = 1/prob_thresh,
         LLR10k = qnorm(1/10000, mean=lmeanLLR, sd=lsdLLR),
         LLR100k = qnorm(1/100000, mean=lmeanLLR, sd=lsdLLR))

opti400 = opti400 %>%
  mutate(LLLR = log(-LLR))
linear_model <- lm(LLLR ~ nSamps, data = opti400)
predictions <- predict(linear_model, newdata = data.frame(nSamps = opti400$nSamps), interval = "prediction", level = 1 - 1/100000)
# Create a new dataframe with the predictions
prediction_df <- data.frame(nSamps = opti400$nSamps, fit = predictions[, "fit"], lwr = predictions[, "lwr"], upr = predictions[, "upr"])

mean_df <- opti400 %>%
  group_by(nSamps) %>%
  summarise(mean_LLLR = mean(LLLR))

logthresh = round(log(-thresh),2)
test = prediction_df %>%
  filter(lwr < logthresh)
minsamps = min(test$nSamps)

# Plot the data, mean line, and lower bound quantile line
ggplot(opti400, aes(x = nSamps, y = LLLR)) +
  geom_point(alpha = 0.3) +  # Scatter plot of all observations
  geom_line(data = mean_df, aes(x = nSamps, y = mean_LLLR, color='data average')) +  # Mean line
  geom_hline(aes(yintercept=logthresh, color='Null-dist. 50th %-tile')) +
  geom_line(data = prediction_df, aes(x = nSamps, y = lwr, color='est. best LLLR after\n100k iterations')) +  # Lower bound quantile line
  geom_vline(aes(xintercept=minsamps, color='est. minimum\nsub-sample size')) +
  labs(title = "Distance to population vs sample-size") +
  xlab("sub-sample size") +
  ylab("log(-LLR)") +
  theme(legend.position="bottom",legend.title=element_blank())
fnsvg = paste(folder, "LLLRvsNSamps.svg", sep="/")
fnpng = paste(folder, "LLLRvsNSamps.png", sep="/")
ggsave(fnsvg, width=6, height=4)
ggsave(fnpng, width=6, height=4)
#conclusion: we need 300 auto-samples across 100k iterations
#side idea number 1: change itersampsize to iterate in a half-halving algorithm instead of in order, going up a half if the iteration is successful, down if it fails
#side idea number 2: [DONE] model each iteration using a gaussian, which can then be used to model the probability of a successful iteration in n-iterations


candidates = generateCandidateSubSamples(BCOWD400$caseNo, 300, 100000)
cand_evals = multiEval(BCOWD400, candidates, 'caseNo', colbinds)

best_subsample = cand_evals %>%
  filter(LLR == max(LLR))
best_LLR = best_subsample$LLR

#plot a histogram of the LLRs in auto400
ggplot() +
  geom_histogram(data = auto400, aes(x=LLR)) +
  geom_vline(aes(xintercept = thresh, color = 'null-dist. 50th %-tile')) +
  #geom_label(aes(x=thresh, y=0, hjust='left', vjust='bottom', label=paste("50th %-tile:", round(thresh,2)))) +
  geom_vline(aes(xintercept = best_LLR, color='sub-sample LLR'))+
  ylab("count") +
  xlab("Log-Likelihood Ratio") +
  theme(legend.position="bottom",legend.title=element_blank()) +
  expand_limits(x=0)+
  coord_flip() +
  ggtitle("null-distribution of LLR")
fnsvg = paste(folder, "LLR_null_dist_w_subsample.svg", sep="/")
fnpng = paste(folder, "LLR_null_dist_w_subsample.png", sep="/")
ggsave(fnsvg, width=6, height=4)
ggsave(fnpng, width=6, height=4)
