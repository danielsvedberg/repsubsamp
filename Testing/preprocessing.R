library(tidyverse) #library important for data-wrangling
library(lubridate) #library important for managing date-time data
library(stringr)
#read the burger court opinion writing database pre-processing
justices = c("WEBurger", "HLBlack", "WODouglas", "JHarlan2", "WJBrennan",
             "PStewart", "BRWhite", "TMarshall", "HABlackmun", "LFPowell",
             "WHRehnquist","JPStevens", "SDOConnor", "Conference")
recipients = list("CJ" = "WEBurger", "BL" = "HLBlack", "WD" = "WODouglas", "JH" = "JHarlan2", "WB" = "WJBrennan",
                  "PS" = "PStewart", "BW" = "BRWhite", "TM" = "TMarshall", "HB" = "HABlackmun", "LP" = "LFPowell",
                  "WR" = "WHRehnquist", "JS" = "JPStevens", "SO" = "SDOConnor", "Conf" = "Conference", "NULL" = "NULL")

justiceNums = c("WEBurger"=1, "HLBlack"=2, "WODouglas"=3, "JHarlan2"=4, "WJBrennan"=5,
                "PStewart"=6, "BRWhite"=7, "TMarshall"=8, "HABlackmun"=9, "LFPowell"=10,
                "WHRehnquist"=11, "JPStevens"=12, "SDOConnor"=13, "Conference"=14, "NULL"= NA)

getBurgerDocuments = function(burger){
  burger = read_csv('2009-08-06.csv') #load dataset from csv file.
  burger$date = as.Date(as.character(burger$date),"%Y%m%d")

  burgerDocuments = burger %>%
    dplyr::select(c(id, us, term, action1:action5, justice, date, recipient, draft, docket1)) %>%
    #get case ID into universal form
    separate(us, c("usVol", "usPg"), remove = FALSE) %>%
    filter(!is.na(usVol), !is.na(usPg)) %>%
    mutate(justice = as.numeric(justice),
           usVol = as.numeric(usVol),
           usPg = as.numeric(usPg)) %>%
    unite(caseNo, usVol:usPg, remove = FALSE) %>% #create caseNo variable
    separate(docket1, c("docket1", NA), sep = " ") %>% #get primary docket
    #coalesce action code columns into one column
    group_by(id) %>%
    pivot_longer(cols = action1:action5, names_to = "action_number", values_to = "action_code") %>%
    dplyr::select(-action_number) %>%
    filter(!is.na(action_code)) %>%
    mutate(action_code = list(action_code)) %>%
    distinct() %>%
    mutate(nullcol = "NULL",
           justice_name = justices[justice],
           recipient = coalesce(recipient,nullcol),
           recipient = list(str_split(recipient,",")),
           recipient_name = list(as.character(recipients[unlist(recipient)])),
           recipient = list(unname(justiceNums[unlist(recipient_name)]))) %>% #,
           #joined_name = list(as.character(recipients[unlist(justices_joined)]))) %>%
    mutate_at(c('usVol', 'usPg'), as.numeric) %>%
    distinct() %>%
    group_by(id) %>%
    mutate(threat = 410 %in% unlist(action_code),
           suggest = 400 %in% unlist(action_code),
           barg = any(c(400,410) %in% unlist(action_code)),
           threatOnly = threat & !(suggest),
           suggestOnly = suggest & !(threat),
           threat_suggest = all(threat,suggest)) %>%
    filter(caseNo != "424_1") %>% #remove buckley bc that shit sucks and violates all conventions
    ungroup() %>%
    dplyr::select(-us, -nullcol)
}

burgerDocuments1 = getBurgerDocuments(burger)

addSCBD = function(burgerDocuments){
  scbd_justice = read_csv('SCDB_2022_01_justiceCentered_Citation.csv')
  scbd_jk = scbd_justice %>% #justice key, this will get the author name from the ID numbers in the majOpinWriter column
    dplyr::select(justice, justiceName) %>%
    ungroup() %>%
    distinct() %>%
    rename(majOpinWriter = justice, opAuth = justiceName) #change names in the key to match the majOpinWriter column

  scbd_justice = scbd_justice %>%
    mutate(dateArgument = as.Date(dateArgument, "%m/%d/%Y"),
           majority = majority - 1,
           voteMajOp = (vote == 1)) %>%
    group_by(caseId) %>%
    mutate(nMajOpVotes = sum(voteMajOp, na.rm = TRUE)) %>%
    group_by(docketId) %>%
    separate(docketId, c(NA,NA,"subCasePriority")) %>%
    separate(usCite, c("usVol", NA, NA, "usPg")) %>%
    mutate(usVol = as.integer(usVol),
           usPg = as.integer(usPg)) %>%
    group_by(caseIssuesId) %>%
    mutate(subCasePriority = as.integer(subCasePriority)) %>% #use docket ids to find lead case
    filter(subCasePriority == 1) %>% #in consolidated cases, filter out everything but the lead case
    mutate(usVol = as.integer(usVol),
           usPg = as.integer(usPg)) %>%
    ungroup() %>%
    distinct() %>%
    group_by(justiceName) %>%
    arrange(justiceName, term, caseId) %>%
    mutate(first_term = min(term, na.rm=TRUE),
           terms_served = term-first_term,
           cases_served = 1,
           cases_served = cumsum(cases_served)) %>%
    rename(justice_name = justiceName) %>%
    ungroup() %>%
    dplyr::select(caseId, justice_name, dateArgument, cases_served, terms_served, majority, voteMajOp, nMajOpVotes) %>%
    distinct()

  #scbd preprocessing
  scbd = read_csv('SCDB_2021_01_caseCentered_Docket.csv') #load scbd
  salience = read_csv('WashU Salience Data.csv') #load salience
  scbd = left_join(scbd, salience, by = c("caseId")) #merge scbd and salience
  scbd = left_join(scbd, scbd_jk, by = c("majOpinWriter")) %>%
    filter(caseId != "1973-060") %>%
    group_by(docketId) %>%
    separate(docketId, c(NA,NA,"subCasePriority")) %>%
    group_by(caseIssuesId) %>%
    mutate(subCasePriority = as.integer(subCasePriority)) %>% #use docket ids to find lead case
    filter(subCasePriority == 1) %>% #in consolidated cases, filter out everything but the lead case
    separate(usCite, c("usVol", NA, NA, "usPg")) %>%
    mutate(usVol = as.integer(usVol),
           usPg = as.integer(usPg),
           totVotes = sum(majVotes,minVotes, na.rm=TRUE)) %>%
    group_by(usVol, usPg) %>%
    rename(scbd_term = term) %>%
    dplyr::select(usVol, usPg, scbd_term, caseId, caseName, nytSalience, cqSalience, majVotes, minVotes, opAuth, chief) %>% #, votesForMaj) %>%#, dateArgument) %>%
    group_by(caseId) %>%
    distinct() %>%
    mutate(caseName = list(caseName)) %>%
    mutate(voteMargin = majVotes-minVotes) %>%
    dplyr::select(-minVotes) %>%
    ungroup()

  scbd = left_join(scbd_justice,scbd, by = c("caseId"))

  burgerDocuments = left_join(burgerDocuments, scbd, by = c("usVol", "usPg","justice_name")) %>%
    mutate(term = coalesce(term, scbd_term))%>%
    dplyr::select(-scbd_term) %>%
    relocate(date, .after= usPg) %>%
    relocate(docket1, .after = usPg) %>%
    relocate(action_code, .after = docket1) %>%
    relocate(justice_name, .after = justice) %>%
    group_by(id) %>%
    mutate(year = as.integer(format(date, "%Y")),
           term = as.integer(term),
           year = coalesce(year, term)) %>%
    mutate(oyez_link = paste("oyez.org/cases", term, docket1, sep = "/")) %>%
    group_by(justice,term) %>%
    mutate(yrCases = length(unique(caseNo)),
           yrBargs = sum(barg, na.rm = TRUE),
           yrBargRate = yrBargs/yrCases) %>%
    group_by(justice) %>%
    mutate(avgBargRate = mean(yrBargRate, na.rm = TRUE)) %>%
    ungroup()

  burgerIDs = burgerDocuments %>%
    ungroup() %>%
    dplyr::select(justice, justice_name) %>%
    distinct() %>%
    filter(!is.na(justice))

  majCoals = left_join(scbd_justice, burgerIDs, by = "justice_name") %>%
    filter(!is.na(justice),
           majority == 1) %>%
    group_by(caseId) %>%
    summarise(maj_vote = list(justice),
              maj_vote_names = list(justice_name))

  burgerDocuments = left_join(burgerDocuments, majCoals, by = c("caseId")) %>%
    dplyr::select(-caseId, -usVol, -usPg)


  return(burgerDocuments)
}
burgerDocuments2 = addSCBD(burgerDocuments1)

parse_joins = function(burgerDocuments){
  burgerDocuments = burgerDocuments %>%
    group_by(id) %>%
    arrange(caseNo, date, id) %>%
    mutate(dissOpCirc = any(c(120,121,130,160) %in% unlist(action_code)),
           altOpCirc = any(c(110,111,112,113,140,150,151) %in% unlist(action_code)),
           majOpCirc = any(c(100,101,104) %in% unlist(action_code)),
           anyOpCirc = any(dissOpCirc,altOpCirc,majOpCirc),
           majOpMemo = any(c(200,230,800) %in% unlist(action_code)),
           majOpJoin = any(c(300,301) %in% unlist(action_code)),
           join_memo = any(majOpCirc,majOpMemo,majOpJoin),
           join_alt = any(c(311, 313, 320) %in% unlist(action_code)),
           withdraw_memo = any(c(700, 701) %in% unlist(action_code)),
           non_join = any(dissOpCirc,withdraw_memo),
           join_value = case_when(join_memo == TRUE ~ 1,
                                  non_join == TRUE ~ 0)) %>%
    mutate(opAuthConf = case_when(any(c(majOpCirc,majOpMemo)) & recipient_name == 'Conference' ~ justice_name),
           opAuthRec = case_when((any(c(join_memo, barg, withdraw_memo)) & (recipient_name != 'Conference') & (length(unlist(recipient)) == 1)) ~ as.character(recipient_name))) %>%
    group_by(caseNo, date) %>%
    fill(opAuthConf, .direction = "downup") %>%
    group_by(caseNo) %>%
    fill(opAuthConf, .direction = "downup") %>%
    group_by(id) %>%
    mutate(opAuthConf = coalesce(opAuthRec,opAuthConf),
           opAuth = coalesce(opAuthConf, opAuth),
           opAuthID = justiceNums[opAuth],) %>%
    group_by(caseNo) %>%
    fill(c(opAuth), .direction = "downup") %>%
    group_by(caseNo,opAuth,justice,date) %>%
    fill(join_value) %>%
    group_by(caseNo,opAuth,justice) %>%
    arrange(caseNo,opAuth,justice,date) %>%
    fill(join_value, .direction = "down") %>%
    ungroup() %>%
    replace_na(list(join_value = 0)) %>%
    arrange(caseNo, date, id) %>%
    group_by(caseNo) %>%
    filter(!is.na(chief)) %>%
    mutate(bargStart = min(date[((anyOpCirc == TRUE) | (join_memo == TRUE))], na.rm=TRUE),#ifelse(majOpCirc == TRUE, min(date[majOpCirc == TRUE]), NA),
           n_case_documents = n(),
           nMajOpCirc = cumsum(majOpCirc),
           nAnyOpCirc = cumsum(anyOpCirc)) %>%
    group_by(id) %>%
    mutate(exposure_time = date - bargStart) %>%
    ungroup()
return(burgerDocuments)
}
burgerDocuments3 = parse_joins(burgerDocuments2)

addMQS = function(data){
  MQS = read_csv("justicesMQS.csv")
  justIni = list("WEB",  "HLB",  "WOD",  "JMH",   "WJB",  "PS",   "BRW",  "TM",   "HAB",  "LFP",  "WHR",  "JPS",  "SDO", "")
  justNms = list("WEBurger","HLBlack","WODouglas","JHarlan2","WJBrennan","PStewart","BRWhite","TMarshall","HABlackmun","LFPowell","WHRehnquist","JPStevens","SDOConnor","unsigned")
  names(justIni) = justNms

  MQS_ = MQS %>%
    ungroup() %>%
    dplyr::select('term', 'justiceName', 'post_med') %>%
    filter(justiceName %in% justNms)

  MQS_A = MQS_ %>%
    rename(recipient_name = justiceName, RecipientMQS = post_med)

  MQS_B = MQS_ %>%
    rename(justice_name = justiceName, BargainerMQS = post_med)

  MQS_C = MQS_ %>%
    rename(opAuth = justiceName, opAuthMQS = post_med)

  MQS_D = MQS_ %>%
    rename(maj_vote_names = justiceName, majVoteMQS = post_med)


  data = left_join(data,MQS_B, by = c("justice_name","term")) #note, SDO does not have 1980 MQS, trying year instead of term, switch back if problem
  data = left_join(data,MQS_C, by = c("opAuth","term"))

  data = data %>%
    unnest(c(recipient,recipient_name))
  data = left_join(data,MQS_A, by = c("recipient_name","term"))
  data = data %>%
    group_by(id) %>%
    mutate(recipient = list(recipient),
           recipient_name = list(recipient_name),
           RecipientMQS = list(RecipientMQS)) %>%
    distinct() %>%
    unnest(c(maj_vote_names))

  data = left_join(data,MQS_D, by = c("maj_vote_names", "term")) %>%
    group_by(id) %>%
    mutate(maj_vote_names = list(maj_vote_names),
           majVoteMQS = list(majVoteMQS),
           majMQSAvg = mean(unlist(majVoteMQS)),
           majMQSSD = sd(unlist(majVoteMQS))) %>%
    distinct()

  return(data)
}

burgerDocuments4 = addMQS(burgerDocuments3)

addCoalition = function(burgerDocuments){
  mergeMQS = burgerDocuments %>%
    ungroup() %>%
    dplyr::select(caseNo,justice,BargainerMQS, term) %>%
    distinct()

  bD = burgerDocuments %>%
    ungroup() %>%
    dplyr::select(caseNo, date, term, justice, join_value, opAuth, opAuthID, opAuthMQS) %>%
    filter(!is.na(date)) %>%
    group_by(caseNo,justice, date, term, opAuth, opAuthID, opAuthMQS) %>%
    arrange(justice,caseNo,opAuth,date) %>%
    mutate(join_value = tail(join_value,1)) %>%
    distinct() %>%
    pivot_wider(names_from = justice, values_from = join_value) %>%
    group_by(caseNo, opAuth, opAuthID, opAuthMQS) %>%
    arrange(caseNo, opAuth, date) %>%
    fill(everything(), .direction = "down") %>%
    group_by(caseNo,opAuth,opAuthID, opAuthMQS, date, term) %>%
    mutate_at(vars(c(-caseNo,-date,-term, -opAuth, -opAuthID, -opAuthMQS)), ~ sum(unlist(.x),na.rm=TRUE)) %>%
    ungroup() %>%
    pivot_longer(!c(caseNo,date,term, opAuth, opAuthID,opAuthMQS), names_to = "justice", values_to = "joined") %>%
    mutate(justice = as.numeric(justice),
           joined = ifelse(justice == opAuthID, 1, joined)) %>%
    left_join(mergeMQS, by = c("caseNo","term","justice")) %>%
    dplyr::select(-term) %>%
    filter(joined == 1) %>%
    distinct() %>%
    group_by(caseNo, date, opAuth, opAuthID, opAuthMQS) %>%
    arrange(caseNo,opAuth, date) %>%
    summarise(coalition = list(unique(justice)),
              coalMQS = list(BargainerMQS),
              coalSize = length(unlist(coalition))) %>%
    group_by(caseNo,opAuth) %>%
    mutate(prev_coalition = lag(coalition, default = list(unique(unname(opAuthID)))),
           prev_coalMQS = lag(coalMQS, default = list(unique(unname(opAuthMQS)))),
           prev_coalSize = lag(coalSize, default = 1),
           max_coalSize = max(coalSize, na.rm=TRUE)) %>%
    group_by(caseNo,date) %>%
    mutate(prev_coalAvg = mean(unlist(prev_coalMQS),na.rm=TRUE),
           prev_coalSD = sd(unlist(prev_coalMQS)),
           prev_coalSD = replace_na(prev_coalSD, 0)) %>%
    group_by(caseNo) %>%
    dplyr::select(-opAuthID,-opAuthMQS) %>%
    ungroup() %>%
    distinct()

  burgerDocuments_ = left_join(burgerDocuments, bD, by = c("caseNo", "date", "opAuth")) %>%
    group_by(caseNo, opAuth, date) %>%
    group_by(id) %>%
    replace_na(list(nytSalience = 0)) %>%
    mutate(RecipientMQS = ifelse(recipient_name == "Conference",majMQSAvg, mean(unlist(RecipientMQS))),
           MQSdiff2rec = BargainerMQS - mean(unlist(RecipientMQS),na.rm=TRUE),
           MQSdiff2Auth = BargainerMQS - opAuthMQS,
           MQSdiff2coal = BargainerMQS - prev_coalAvg,
           MQSdiff2maj = BargainerMQS - mean(unlist(majVoteMQS), na.rm=TRUE),
           absMQSdiff2rec = abs(MQSdiff2rec),
           absMQSdiff2Auth = abs(MQSdiff2Auth),
           absMQSdiff2coal = abs(MQSdiff2coal),
           absMQSdiff2maj = abs(MQSdiff2maj),
           absRecipientMQS = abs(mean(unlist(RecipientMQS),na.rm=TRUE)),
           absBargMQS = abs(BargainerMQS),
           prevCoalMarg2MajVote = majVotes - prev_coalSize)
  return(burgerDocuments_)
}
burgerDocuments5 = addCoalition(burgerDocuments4)
burgerDocuments = burgerDocuments5

cn = colnames(burgerDocuments)
for (i in cn){
  if (is.logical(burgerDocuments[[i]])){
    burgerDocuments[[i]] = as.integer(burgerDocuments[[i]])
  }
}

lists_to_text <- function(df, delimiter = ",") {
  df %>%
    mutate(across(where(is.list), ~map_chr(., paste, collapse = delimiter)))
}

burgerDocSave = lists_to_text(burgerDocuments)

write_csv(burgerDocSave, file = "preprocessedBurgerDocuments.csv")
