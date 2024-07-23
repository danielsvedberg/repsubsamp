# Repsubsamp: extraction and validation of "representative" subsamples
Daniel A. Svedberg (2024)
## Introduction
Data is often delivered in a dataframe format, 
but sometimes we need to extract a subsample from these data, 
while ensuring that it is "representative" of the original data. 

In our specific case, we made use of the
[Burger Court Opinion Writing Database](http://supremecourtopinions.wustl.edu/?rt=dataset/archive)
(BCOWD) [1] to study the effect of bargaining strategies used by Supreme Court justices,
through analysis of the memos they wrote to each other in the majority-opinion writing process.
Because our analysis required a legal expert to manually read and score each bargaining exchange, 
we sought to extract the smallest-necessary subsample of the data that was still representative of the original data. 

More specifically, we needed a subsample that contained cases with similar proportions of various variables, 
like: the year the case was decided, the margin of the vote, and the salience of the issue. 

To solve this problem, the approach employed by repsubsamp is to characterize the 
distribution of each *norming variable*, in the population as well as in the subsample,
and quantify the similarity between the subsample and the population data by calculating the 
likelihood-ratio that the distributions of the norming variables in the subsample 
are the same as in the population data.

Given the tools to quantify the similarity of a subsample and the population data, 
we also needed a way to determine how similar is similar enough. 

To solve this problem, the approach employed by repsubsamp is to perform a permutation test 
of the likelihood that the subsample is *not* significanly different from the population data, 
by resampling the population data with replacement to calculate the null-distribution of the likelihood-ratio. 
We then consider any sub-sample that is at least as similar to the population data as 50% 
of resamples in the null-distribution to be "representative" of the population data.

Simply put, repsubsamp integrates tools for:
1) creating sub-samples from a dataframe
2) calculating the similarity of the sub-sample to the population data
3) calculating the null-distribution of resamples' similarity to population data, and using it to benchmark the "representativeness" of a sub-sample 


## Similarity metric

The Likelihood-ratio test statistic used to quantify similarity of a subsample to the population ($\lambda_{LR}$), 
which equals the log of the product of the likelihood ratios 
across each norming variable *i* of *n* possible norming variables, 
where a variable’s likelihood ratio equals the likelihood 
$(\mathcal{L})$ of the population data $(x_0)$ given the parameter 
estimates generated from the entire population $(\theta_0)$, 
divided by the $\mathcal{L}$ of $x_0$ given the parameter estimates 
generated from the sample $(\theta_s)$:

$$ 
\lambda_{LLR}=\ln{\left[\prod_{i=1}^{n}\frac{\mathcal{L}_i\left(\theta_0|x_0\right)}{\mathcal{L}_i\left(\theta_s|x_0\right)}\right]}\ 
$$

The closer the $\lambda_{LLR}$ is to zero, the more similar the subsample is to the population data.

## Example Pipeline

The repsubsamp example pipeline can be found in ` repsubsamp/Testing/subsample_size_optimization.R `.

The core functions of repsubsamp are found in ` repsubsamp/R/subsampling_functions.R `.
```R
source("R/subsampling_functions.R") #load the sub sampling functions
```

This example also makes use of several libraries:
```R
library(tidyverse) #used to wrangle data
library(lubridate) #library important for managing date-time data
library(stringr) #used to process some strings
library(fitdistrplus) #used to fit distributions to data
```

The example data can be found in ` repsubsamp/data_raw/BCOWDCases.csv `. 
This dataframe includes the BCOWD, as well as data from the Supreme Court Database (SCBD) [2] 
and the Issue Salience Database [3]--data not found in BCOWD 
that were eventually used to validate the subsample.

```R
BCOWDCases <- read_csv("repsubsamp/data_raw/BCOWDCases.csv")
> head(BCOWDCases)
# A tibble: 6 × 14
  caseNo usVol  usPg n_docs action_code                term caseId caseName nytSalience cqSalience voteMargin     n unanimous n_cases
  <chr>  <dbl> <dbl>  <dbl> <chr>                     <dbl> <chr>  <chr>    <lgl>       <lgl>           <dbl> <dbl> <lgl>       <dbl>
1 396_13   396    13      8 300;130;101;101;300;160;…  1969 1969-… SIMPSON… FALSE       FALSE               5     1 FALSE           1
2 396_19   396    19     43 103;810;100;810;103;103;…  1969 1969-… ALEXAND… TRUE        FALSE               8     1 FALSE           1
3 396_28   396    28     13 300;300;120;120;120;120;…  1969 1969-… DEBACKE… FALSE       FALSE               4     1 FALSE           1
4 396_41   396    41      9 300;320;400;300;300;101;…  1969 1969-… BROCKIN… FALSE       FALSE               8     1 FALSE           1
5 396_45   396    45     11 300;410;300;203;120;120;…  1969 1969-… HALL et… FALSE       FALSE               4     1 FALSE           1
6 396_57   396    57     13 610;300;300;100;100;100;…  1969 1969-… ANDERSO… FALSE       FALSE               7     1 FALSE           1
```

In this specific example, we want to analyze a specific subset of the data--every 
case in which a "suggestion to majority opinion" memo was circulated,
which is recorded by action code #400:

```R
BCOWD400 = BCOWDCases %>% #extract cases from BCOWD with code 400
  filter(str_detect(action_code, '400'))
```

Next, we want to set up the norming variables that we want to use to characterize the population data. 
This is done by setting "column bindings" that link selected column names in the dataframe to the distribution names

```R
#set the column bindings by creating named list linking column names to distribution names
cols = c('n_docs','term','nytSalience','cqSalience','unanimous') #list of columns to be modeled
distnms = c('poiss','multinom','binom','binom','binom') #matching list of distribution labels for each column
colbinds = setNames(distnms, cols) #make a named list where each distribution label can be called by the column name
```
Whenever you do this step for your own data, you will need to
determine the appropriate distribution for each norming variable you choose. 
Repsubsamp currently supports the following distributions:
* 'norm' for Normal
* 'exp' for Exponential
* 'binom' for Binomial
* 'multinom' for Multinomial
* 'poiss' for Poisson
* 'gamma' for Gamma (warning, this is slow)

Next, we generate a null-distribution of the likelihood-ratio test statistic. 

First, we generate 1000 resamples of the IDs in the population data, with replacement:
```R
auto400cands = generateAutoSamples(BCOWD400$caseNo, 1000) #arguments: IDs (list) and # of iterations (int)
returns a matrix of resampled IDs [#ids, #iterations]
```
Each resample contains the same number of IDs as the population data, 
but the IDs are drawn with replacement to create some variation. 

Then, we calculate the likelihood-ratio test statistic for each resample:
```R
#multiEval calculates the distance between the population and each sample in the matrix, according the columns specified in colbinds
#arguments: population dataframe, matrix of resampled IDs, ID column name, and column bindings
auto400 = multiEval(BCOWD400,auto400cands,'caseNo',colbinds) 
#returns a dataframe with each resample's similarity to the population
```
The 'LLR' column in the dataframe contains the likelihood-ratio test statistic for each resample.

Next, we want to get the 50th percentile of the null-distribution of the likelihood-ratio test statistic,
which is the similarity threshold of our permutation test:
```R
#get the 95% confidence interval of the distances for:
percentiles400LLR = quantile(auto400$LLR, c(0.025,0.5,0.975)) #Log-likelihood ratio
thresh = percentiles400LLR['50%']
```

Let's visualize the null distribution by plotting a histogram with the similarity threshold:
```R
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
```
![image of the null distribution](https://github.com/danielsvedberg/repsubsamp/blob/9c38b7f96547f2b506a066133d6eeeaf21185719/Testing/LLR_null_dist.png)

Our next task is to identify how large of a subsample is needed to achieve 
similarity that is greater than the threshold. 
The approach to this problem is to iterate over many different subsample sizes, 
and compare each iteration to the similarity threshold. For this analysis, 
we will use only 100 iterations per subsample size, 
to characterize the relationship between subsample size and similarity to the population data.

```R
#itersampsize takes the dataframe, id column name (str), column bindings (named list)
#number of repeats per iteration (int), minimum sub-sample size (int), and maximum sub-sample size (int)
opti400 = itersampsize(BCOWD400, 'caseNo', colbinds, 100, 100, 700)
```

Let's plot the results of this comparison:
```R
#plot opti vs LLR
ggplot(data = opti400, aes(x=nSamps, y=LLR))+
  geom_point(alpha=0.1)+
  stat_summary(fun=mean, geom="line", color="blue") +
  geom_hline(yintercept=thresh) +
  geom_label(aes(y=-100, hjust='right', vjust='top', x=700, label=paste("Null-dist. 50th %-tile:", round(thresh,2)))) +
  ylab("log-likelihood ratio") +
  xlab("Sub-sample size") +
  ggtitle("Sub-sample size vs distance to population after 100 iterations")
```
![subsample similarities over many different subsample sizes](https://github.com/danielsvedberg/repsubsamp/blob/9c38b7f96547f2b506a066133d6eeeaf21185719/Testing/LLRvsNSamps.png)

Here we can see, that given around 100 iterations, a subsample size of around 400 cases (of 1298) 
will achieve the similarity threshold. 

Consider that the similarity of a subsample to the population data is stochastic;
therefore, the likelihood that a subsample achieves the similarity threshold
is a function of both the subsample size *and the number of iterations*. Ergo, increasing the number of 
also increases the likelihood that a smaller subsample can reach the similarity threshold. 

The distribution of subsample similarities across iterations enables us to forecast the highest
similarity achievable for every subsample size, given a chosen number of iterations.

The largest number of iterations I want to run on my laptop is 100,000,
so I will forecast the smallest subsample size likely to achieve the similarity threshold
given 100,000 iterations. You might choose a larger number if your computer is faster than mine. 

To do this, I will use the distribution of subsample similarities across iterations,
and to create a linear model of the 1/100,000-th percentile for each subsample size. 

In order to create a linear model, the data must first be linearized. 
If you look examine the data above, the distribution of LLRs is roughly log-normal. 
Therefore, I will normalize the distribution of LLRs by taking its log, before creating the model:

```R
#plot a histogram of the log-normal distribution of the LLRs
ggplot() +
  geom_histogram(data = auto400, aes(x=log(-LLR))) +
  geom_vline(xintercept = log(-thresh)) +
  geom_label(aes(x=log(-thresh), y=0, hjust='left', vjust='bottom', label=paste("50th %-tile:", round(log(-thresh),2)))) +
  ylab("count") +
  coord_flip() +
  ggtitle("null-distribution is log-normally distributed")
```
![log of the null distribution](https://github.com/danielsvedberg/repsubsamp/blob/9c38b7f96547f2b506a066133d6eeeaf21185719/Testing/LogLR_null_dist.png)

While this transformation is not perfect, it is close enough to create a linear model:

![log of the subsample size iterations](https://github.com/danielsvedberg/repsubsamp/blob/9c38b7f96547f2b506a066133d6eeeaf21185719/Testing/LogLLRvsNSamps.png)

Next, we create use the linear model of the subsample size similarity vs iterations, 
to find the 1/100,000-th percentile for each subsample size, and plot the results:
```R
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
  ```

![log of the subsample size iterations with model](https://github.com/danielsvedberg/repsubsamp/blob/9c38b7f96547f2b506a066133d6eeeaf21185719/Testing/LLLRvsNSamps.png)

This model forecasts that a subsample size of at least 300 cases will achieve the 
similarity threshold within 100,000 iterations. Therefore, 
we will generate a final sub-sample with a size of 300 cases over 100,000 iterations:

```R
candidates = generateCandidateSubSamples(BCOWD400$caseNo, 300, 100000) #generate 100,000 subsamples of 300 cases
cand_evals = multiEval(BCOWD400, candidates, 'caseNo', colbinds) #evaluate the similarity of each subsample to the population

best_subsample = cand_evals %>% #find the best subsample in the candidates
  filter(LLR == max(LLR))
best_LLR = best_subsample$LLR
```
As forecasted, the best subsample from 100,0000 iterations achieved the similarity threshold:

![null distribution with best subsample LLR](https://github.com/danielsvedberg/repsubsamp/blob/9c38b7f96547f2b506a066133d6eeeaf21185719/Testing/LLR_null_dist_w_subsample.png)

## Conclusion
Here I demonstrated how to use repsubsamp to:

1) establish and characterize norming variables in a dataframe
2) generate a null-distribution of the likelihood-ratio statistic to establish a threshold similarity target
3) iterate over a range of sub-sample sizes to characterize how similarity varies with subsample size
4) forecast the smallest subsample size likely to achieve the similarity threshold within a chosen number of iterations
5) generate a smallest-possible subsample that achieves the similarity threshold within the chosen number of iterations

While this example used the BCOWD in a specific use-case, repsubsamp can be used with any dataframe that contains a variety of columns that can be modeled using different distributions.

## References 
[1] Wahlbeck, Paul J., James F. Spriggs II, and Forrest Maltzman. 2011. The Burger Court Opinion-Writing Database. http://supremecourtopinions.wustl.edu

[2] Harold J. Spaeth, Lee Epstein, et al. 2023 Supreme Court Database, Version 2021 Release 1. URL: http://supremecourtdatabase.org

[3] Lee Epstein and Jeffrey A. Segal. 2000. “Measuring Issue Salience.” American
Journal of Political Science. 44(1): 66-83.