library(tidyverse) #library important for data-wrangling
library(lubridate) #library important for managing date-time data
library(stringr)

###READNIG IN DATA AND PREPROCESSING############################################

#read the BCOWD court opinion writing database pre-processing
BCOWD = read_csv('data_raw/Burger Court Opinion Writing Database.csv') #load dataset from csv file.
BCOWD$date = as.Date(as.character(BCOWD$date),"%Y%m%d")

#preprocessing BCOWD to be case-centered data
BCOWDCases = BCOWD %>%
  dplyr::select(c(id, term, us, action1:action5)) %>% #grab columns of interest: case ID, term, US index, and the action codes
  separate(us, c("usVol", "usPg")) %>% #separate the US case ID into volume and page, for later merging with scbd
  pivot_longer(cols = action1:action5, names_to = "action_number", values_to = "action_code") %>% #longform the action codes
  dplyr::select(-action_number) %>% #eliminate the incidentally created column
  filter(!is.na(action_code)) %>% #remove all rows with NA in action code
  ungroup() %>%
  distinct() %>%
  group_by(usVol, usPg) %>% #now grouping by case
  mutate(action_code = list(action_code)) %>% #put each action code in each case into a singular list
  summarise(n_docs = n(), #consolidate rows so each case just one row
            action_code = action_code) %>%
  distinct() %>%
  group_by(usVol,usPg) %>%
  mutate_at(c('usVol', 'usPg'), as.numeric) #make the indexing columns numeric format for merging with scbd

#scbd preprocessing
scbd = read_csv('data_raw/SCDB_2021_01_caseCentered_Docket.csv') #load scbd database
salience = read_csv('data_raw/WashU Salience Data.csv') #load salience data
scbd = left_join(scbd, salience, by = c("caseId")) #merge scbd and salience along caseId
#preprocessing scbd data
scbd = scbd %>%
  separate(docketId, c(NA,NA,"subCasePriority")) %>%
  mutate(subCasePriority = as.integer(subCasePriority)) %>% #use docket ids to find lead case
  filter(subCasePriority == 1) %>% #in consolidated cases, filter out everything but the lead case
  separate(dateArgument, c("mon", "day", NA)) %>% #split date into separate columns, month and day
  mutate(year = term) %>% #rename the term column year
  relocate(year, .after = day) %>%
  unite("date", mon:year, sep = "/") %>% #create a new date column with matching format
  mutate(dateArgument = as.Date(date, "%m/%d/%Y")) %>%
  separate(usCite, c("usVol", NA, NA, "usPg")) %>%
  mutate(usVol = as.integer(usVol),
         usPg = as.integer(usPg)) %>%
  dplyr::select(usVol, usPg, term, caseId, caseName, nytSalience, cqSalience, majVotes, minVotes) %>%
  distinct() %>%
  group_by(caseId) %>%
  mutate(caseName = list(caseName)) %>%
  mutate(voteMargin = majVotes-minVotes) %>%
  dplyr::select( -majVotes, -minVotes) %>%
  group_by(caseId) %>%
  mutate(n = n()) %>%
  distinct()

#the SCBD case IDs in blacklist have database issues which make them incompatible with our analysis
blacklist = list("1974-123", "1978-101","1978-145","1979-154","1980-078","1971-078","1973-059")
BCOWDCases = left_join(BCOWDCases, scbd, by = c("usVol", "usPg")) #merge salience into BCOWD data, only keeping cases that were in BCOWD to begin with
BCOWDCases = BCOWDCases %>%
  filter(!(caseId %in% blacklist),
         !is.na(caseId),
         !is.na(voteMargin),
         !is.na(nytSalience),
         !is.na(cqSalience)) %>%
  mutate(nytSalience = (nytSalience==1),
         cqSalience = (cqSalience==1),
         unanimous = (voteMargin==9)) %>%
  unite(caseNo, usVol:usPg, remove = FALSE) %>%
  mutate(n_cases = n())



