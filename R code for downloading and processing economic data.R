library(censusapi)
library(tidyverse)
library(tigris)
library(sf)
mykey <- 'baaaffb5ed3accd8dfa53c6f827659d43fcdfa21'
options(tigris_use_cache = TRUE)


excluded.states <- c("AS","GU","MP","PR","UM","VI","AK", "HI","DC")
states <- unique(subset(fips_codes, !(fips_codes$state %in% excluded.states)))
st <- paste0("state:",unique(states$state_code))

timber.sic <- c("0811",	"0831","0851","2411","2421","2426","2429","2431","2435","2436","2439","2441",
            "2448","2449","2451","2452","2491","2493","2499","2611","2621","2631","3553","5031")

tot.emp <- map(st, function(x) getCensus(name="cbp",
                                             vintage = "1987",
                                             key = mykey,
                                             vars = c("SIC_TTL","EMP", "EMP_F"),
                                             region = "county:*",
                                             regionin = x,
                                             sic = "00")) %>% 
  do.call(rbind, .) %>% 
  filter(., county != "999") %>% 
  select(., state, county, TOTALEMP=EMP) 

timber.state.1987 <- lapply(st, function(x) map(timber, possibly(function(y) getCensus(name = "cbp",
                                             vintage = "1987",
                                             key = mykey,
                                             vars = c("SIC_TTL","EMP", "EMP_F"),
                                             region = "county:*",
                                             regionin = x,
                                             sic = y), NA))) 
  


timber.1987.df <- lapply(seq_along(timber.state.1987), function(x) do.call(rbind,timber.state.1987[[x]][1:24])) %>% 
  do.call(rbind, .) %>% 
  filter(., county != "999") %>% 
  mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
         EMP = ifelse(EMP_F == "a", 1, EMP),
         EMP = ifelse(EMP_F == "b", 20, EMP),
         EMP = ifelse(EMP_F == "c", 100, EMP),
         EMP = ifelse(EMP_F == "e", 250, EMP),
         EMP = ifelse(EMP_F == "f", 500, EMP),
         EMP = ifelse(EMP_F == "g", 1000, EMP),
         EMP = ifelse(EMP_F == "h", 2500, EMP),
         EMP = ifelse(EMP_F == "i", 5000, EMP),
         EMP = ifelse(EMP_F == "j", 10000, EMP),
         EMP = ifelse(EMP_F == "k", 25000, EMP),
         EMP = ifelse(EMP_F == "l", 50000, EMP),
         EMP = ifelse(EMP_F == "m", 100000, EMP)) %>% 
  group_by(state, county) %>% 
  summarise(AllTimber = sum(as.numeric(EMP), na.rm=TRUE)) %>% 
  right_join(., tot.emp) %>% 
  mutate(AllTimber = ifelse(is.na(AllTimber), 0, AllTimber)) %>% 
  mutate(propTimber = AllTimber/as.numeric(TOTALEMP))


# Total Employment functions ----------------------------------------------


get_total_sic <- function(states, years){
  tot.emp <-lapply(seq_along(years), function(yr) 
    map(states, possibly(function(st) getCensus(name = "cbp",
                                                vintage = years[yr],
                                                key = mykey,
                                                vars = c("SIC_TTL","EMP", "EMP_F"),
                                                region = "county:*",
                                                regionin = st,
                                                sic = "00"), NA)))
  tot.emp.df <- lapply(seq_along(tot.emp), function(x) do.call(rbind,tot.emp[[x]][1:length(states)]) %>% 
                         mutate(year = years[x])) %>% 
    do.call(rbind, .) %>% 
    filter(., county != "999") %>% 
    mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
                                          EMP = ifelse(EMP_F == "a", 19, EMP),
                                          EMP = ifelse(EMP_F == "b", 99, EMP),
                                          EMP = ifelse(EMP_F == "c", 249, EMP),
                                          EMP = ifelse(EMP_F == "e", 499, EMP),
                                          EMP = ifelse(EMP_F == "f", 999, EMP),
                                          EMP = ifelse(EMP_F == "g", 2499, EMP),
                                          EMP = ifelse(EMP_F == "h", 4999, EMP),
                                          EMP = ifelse(EMP_F == "i", 9999, EMP),
                                          EMP = ifelse(EMP_F == "j", 24999, EMP),
                                          EMP = ifelse(EMP_F == "k", 49999, EMP),
                                          EMP = ifelse(EMP_F == "l", 99999, EMP),
                                          EMP = ifelse(EMP_F == "m", 100000, EMP)) %>% 
    select(., state, county, year, TOTALEMP=EMP)
}


get_total_naics97 <- function(states, years){
  tot.emp <-lapply(seq_along(years), function(yr) map(states, possibly(function(st) getCensus(name = "cbp",
                                                                                              vintage = years[yr],
                                                                                              key = mykey,
                                                                                              vars = c("NAICS1997_TTL","EMP", "EMP_F"),
                                                                                              region = "county:*",
                                                                                              regionin = st,
                                                                                              naics1997 = "00"), NA)))
  
  tot.emp.df <- lapply(seq_along(tot.emp), function(x) do.call(rbind,tot.emp[[x]][1:length(states)]) %>% 
                         mutate(year = years[x]))%>% 
    do.call(rbind, .) %>% 
    filter(., county != "999") %>% 
    mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
           EMP = ifelse(EMP_F == "a", 19, EMP),
           EMP = ifelse(EMP_F == "b", 99, EMP),
           EMP = ifelse(EMP_F == "c", 249, EMP),
           EMP = ifelse(EMP_F == "e", 499, EMP),
           EMP = ifelse(EMP_F == "f", 999, EMP),
           EMP = ifelse(EMP_F == "g", 2499, EMP),
           EMP = ifelse(EMP_F == "h", 4999, EMP),
           EMP = ifelse(EMP_F == "i", 9999, EMP),
           EMP = ifelse(EMP_F == "j", 24999, EMP),
           EMP = ifelse(EMP_F == "k", 49999, EMP),
           EMP = ifelse(EMP_F == "l", 99999, EMP),
           EMP = ifelse(EMP_F == "m", 100000, EMP)) %>% 
    select(., state, county, year, TOTALEMP=EMP)
}

get_total_naics02 <- function(states, years){
  tot.emp <-lapply(seq_along(years), function(yr) 
    map(states, possibly(function(st) getCensus(name = "cbp",
                                                vintage = years[yr],
                                                key = mykey,
                                                vars = c("NAICS2002_TTL","EMP", "EMP_F"),
                                                region = "county:*",
                                                regionin = st,
                                                naics2002 = "00"), NA)))
  
  tot.emp.df <- lapply(seq_along(tot.emp), function(x) do.call(rbind,tot.emp[[x]][1:length(states)]) %>% 
                         mutate(year = years[x]))%>% 
    do.call(rbind, .) %>% 
    filter(., county != "999") %>% 
    mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
           EMP = ifelse(EMP_F == "a", 19, EMP),
           EMP = ifelse(EMP_F == "b", 99, EMP),
           EMP = ifelse(EMP_F == "c", 249, EMP),
           EMP = ifelse(EMP_F == "e", 499, EMP),
           EMP = ifelse(EMP_F == "f", 999, EMP),
           EMP = ifelse(EMP_F == "g", 2499, EMP),
           EMP = ifelse(EMP_F == "h", 4999, EMP),
           EMP = ifelse(EMP_F == "i", 9999, EMP),
           EMP = ifelse(EMP_F == "j", 24999, EMP),
           EMP = ifelse(EMP_F == "k", 49999, EMP),
           EMP = ifelse(EMP_F == "l", 99999, EMP),
           EMP = ifelse(EMP_F == "m", 100000, EMP)) %>% 
    select(., state, county, year, TOTALEMP=EMP)
}



get_total_naics07 <- function(states, years){
  tot.emp <-lapply(seq_along(years), function(yr) 
    map(states, possibly(function(st) getCensus(name = "cbp",
                                                vintage = years[yr],
                                                key = mykey,
                                                vars = c("NAICS2007_TTL","EMP", "EMP_F"),
                                                region = "county:*",
                                                regionin = st,
                                                naics2007 = "00"), NA)))
  
  tot.emp.df <- lapply(seq_along(tot.emp), function(x) do.call(rbind,tot.emp[[x]][1:length(states)]) %>% 
                         mutate(year = years[x]))%>% 
    do.call(rbind, .) %>% 
    filter(., county != "999") %>% 
    mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
           EMP = ifelse(EMP_F == "a", 19, EMP),
           EMP = ifelse(EMP_F == "b", 99, EMP),
           EMP = ifelse(EMP_F == "c", 249, EMP),
           EMP = ifelse(EMP_F == "e", 499, EMP),
           EMP = ifelse(EMP_F == "f", 999, EMP),
           EMP = ifelse(EMP_F == "g", 2499, EMP),
           EMP = ifelse(EMP_F == "h", 4999, EMP),
           EMP = ifelse(EMP_F == "i", 9999, EMP),
           EMP = ifelse(EMP_F == "j", 24999, EMP),
           EMP = ifelse(EMP_F == "k", 49999, EMP),
           EMP = ifelse(EMP_F == "l", 99999, EMP),
           EMP = ifelse(EMP_F == "m", 100000, EMP)) %>% 
    select(., state, county, year, TOTALEMP=EMP)
}

get_total_naics12 <- function(states, years){
  tot.emp <-lapply(seq_along(years), function(yr) 
    map(states, possibly(function(st) getCensus(name = "cbp",
                                                vintage = years[yr],
                                                key = mykey,
                                                vars = c("NAICS2012_TTL","EMP", "EMP_F"),
                                                region = "county:*",
                                                regionin = st,
                                                naics2012 = "00"), NA)))
  
  tot.emp.df <- lapply(seq_along(tot.emp), function(x) do.call(rbind,tot.emp[[x]][1:length(states)]) %>% 
                         mutate(year = years[x]))%>% 
    do.call(rbind, .) %>% 
    filter(., county != "999") %>% 
    mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
           EMP = ifelse(EMP_F == "a", 19, EMP),
           EMP = ifelse(EMP_F == "b", 99, EMP),
           EMP = ifelse(EMP_F == "c", 249, EMP),
           EMP = ifelse(EMP_F == "e", 499, EMP),
           EMP = ifelse(EMP_F == "f", 999, EMP),
           EMP = ifelse(EMP_F == "g", 2499, EMP),
           EMP = ifelse(EMP_F == "h", 4999, EMP),
           EMP = ifelse(EMP_F == "i", 9999, EMP),
           EMP = ifelse(EMP_F == "j", 24999, EMP),
           EMP = ifelse(EMP_F == "k", 49999, EMP),
           EMP = ifelse(EMP_F == "l", 99999, EMP),
           EMP = ifelse(EMP_F == "m", 100000, EMP)) %>% 
    select(., state, county, year, TOTALEMP=EMP)
}
# Sector functions --------------------------------------------------------


get_sect_sic <- function(sector, states, years){
  sector.state.year <- lapply(seq_along(years), function(yr)lapply(seq_along(states), function(st) map(sector, possibly(function(y) getCensus(name = "cbp",
                                                                                         vintage = years[yr],
                                                                                         key = mykey,
                                                                                         vars = c("SIC_TTL","EMP", "EMP_F"),
                                                                                         region = "county:*",
                                                                                         regionin = states[st],
                                                                                         sic = y), NA))))
  sector.state <- lapply(seq_along(sector.state.year), 
                          function(x) map(seq_along(states), 
                          function(y) do.call(rbind,sector.state.year[[x]][[y]][1:length(sector)])))
  sector.yr <- lapply(seq_along(sector.state),
                      function(x) do.call(rbind, sector.state[[x]][1:length(states)]))
  sector.df <- lapply(seq_along(sector.yr), function(x) bind_rows(sector.yr[[x]]) %>% 
                        mutate(year = years[x])) %>% 
    do.call(rbind, .) %>% 
    filter(., county != "999") %>% 
    mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
           EMP = ifelse(EMP_F == "a", 19, EMP),
           EMP = ifelse(EMP_F == "b", 99, EMP),
           EMP = ifelse(EMP_F == "c", 249, EMP),
           EMP = ifelse(EMP_F == "e", 499, EMP),
           EMP = ifelse(EMP_F == "f", 999, EMP),
           EMP = ifelse(EMP_F == "g", 2499, EMP),
           EMP = ifelse(EMP_F == "h", 4999, EMP),
           EMP = ifelse(EMP_F == "i", 9999, EMP),
           EMP = ifelse(EMP_F == "j", 24999, EMP),
           EMP = ifelse(EMP_F == "k", 49999, EMP),
           EMP = ifelse(EMP_F == "l", 99999, EMP),
           EMP = ifelse(EMP_F == "m", 100000, EMP)) 
}

get_sect_naics97 <- function(sector, states, years){
  sector.state.year <- lapply(seq_along(years), function(yr)lapply(seq_along(states), function(st) map(sector, possibly(function(y) getCensus(name = "cbp",
                                                                                                                                              vintage = years[yr],
                                                                                                                                              key = mykey,
                                                                                                                                              vars = c("NAICS1997_TTL","EMP", "EMP_F"),
                                                                                                                                              region = "county:*",
                                                                                                                                              regionin = states[st],
                                                                                                                                              naics1997 = y), NA))))
  sector.state <- lapply(seq_along(sector.state.year), 
                         function(x) map(seq_along(states), 
                                         function(y) do.call(rbind,sector.state.year[[x]][[y]][1:length(sector)])))
  sector.yr <- lapply(seq_along(sector.state),
                      function(x) do.call(rbind, sector.state[[x]][1:length(states)]))
  sector.df <- lapply(seq_along(sector.yr), function(x) bind_rows(sector.yr[[x]]) %>% 
                        mutate(year = years[x])) %>% 
    do.call(rbind, .) %>% 
    filter(., county != "999") %>% 
    mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
           EMP = ifelse(EMP_F == "a", 19, EMP),
           EMP = ifelse(EMP_F == "b", 99, EMP),
           EMP = ifelse(EMP_F == "c", 249, EMP),
           EMP = ifelse(EMP_F == "e", 499, EMP),
           EMP = ifelse(EMP_F == "f", 999, EMP),
           EMP = ifelse(EMP_F == "g", 2499, EMP),
           EMP = ifelse(EMP_F == "h", 4999, EMP),
           EMP = ifelse(EMP_F == "i", 9999, EMP),
           EMP = ifelse(EMP_F == "j", 24999, EMP),
           EMP = ifelse(EMP_F == "k", 49999, EMP),
           EMP = ifelse(EMP_F == "l", 99999, EMP),
           EMP = ifelse(EMP_F == "m", 100000, EMP))  
}


get_sect_naics02 <- function(sector, states, years){
  sector.state.year <- lapply(seq_along(years), function(yr)lapply(seq_along(states), function(st) map(sector, possibly(function(y) getCensus(name = "cbp",
                                                                                                                                              vintage = years[yr],
                                                                                                                                              key = mykey,
                                                                                                                                              vars = c("NAICS2002_TTL","EMP", "EMP_F"),
                                                                                                                                              region = "county:*",
                                                                                                                                              regionin = states[st],
                                                                                                                                              naics2002 = y), NA))))
  sector.state <- lapply(seq_along(sector.state.year), 
                         function(x) map(seq_along(states), 
                                         function(y) do.call(rbind,sector.state.year[[x]][[y]][1:length(sector)])))
  sector.yr <- lapply(seq_along(sector.state),
                      function(x) do.call(rbind, sector.state[[x]][1:length(states)]))
  sector.df <- lapply(seq_along(sector.yr), function(x) bind_rows(sector.yr[[x]]) %>% 
                        mutate(year = years[x])) %>% 
    do.call(rbind, .) %>% 
    filter(., county != "999") %>% 
    mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
           EMP = ifelse(EMP_F == "a", 19, EMP),
           EMP = ifelse(EMP_F == "b", 99, EMP),
           EMP = ifelse(EMP_F == "c", 249, EMP),
           EMP = ifelse(EMP_F == "e", 499, EMP),
           EMP = ifelse(EMP_F == "f", 999, EMP),
           EMP = ifelse(EMP_F == "g", 2499, EMP),
           EMP = ifelse(EMP_F == "h", 4999, EMP),
           EMP = ifelse(EMP_F == "i", 9999, EMP),
           EMP = ifelse(EMP_F == "j", 24999, EMP),
           EMP = ifelse(EMP_F == "k", 49999, EMP),
           EMP = ifelse(EMP_F == "l", 99999, EMP),
           EMP = ifelse(EMP_F == "m", 100000, EMP)) 
}

get_sect_naics07 <- function(sector, states, years){
  sector.state.year <- lapply(seq_along(years), function(yr)lapply(seq_along(states), function(st) 
    map(sector, possibly(function(y) getCensus(name = "cbp",
                                              vintage = years[yr],
                                              key = mykey,
                                              vars = c("NAICS2007_TTL","EMP", "EMP_F"),
                                              region = "county:*",
                                              regionin = states[st],
                                              naics2007 = y), NA))))
  sector.state <- lapply(seq_along(sector.state.year), 
                         function(x) map(seq_along(states), 
                                         function(y) do.call(rbind,sector.state.year[[x]][[y]][1:length(sector)])))
  sector.yr <- lapply(seq_along(sector.state),
                      function(x) do.call(rbind, sector.state[[x]][1:length(states)]))
  sector.df <- lapply(seq_along(sector.yr), function(x) bind_rows(sector.yr[[x]]) %>% 
                        mutate(year = years[x])) %>% 
    do.call(rbind, .) %>% 
    filter(., county != "999") %>% 
    mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
           EMP = ifelse(EMP_F == "a", 19, EMP),
           EMP = ifelse(EMP_F == "b", 99, EMP),
           EMP = ifelse(EMP_F == "c", 249, EMP),
           EMP = ifelse(EMP_F == "e", 499, EMP),
           EMP = ifelse(EMP_F == "f", 999, EMP),
           EMP = ifelse(EMP_F == "g", 2499, EMP),
           EMP = ifelse(EMP_F == "h", 4999, EMP),
           EMP = ifelse(EMP_F == "i", 9999, EMP),
           EMP = ifelse(EMP_F == "j", 24999, EMP),
           EMP = ifelse(EMP_F == "k", 49999, EMP),
           EMP = ifelse(EMP_F == "l", 99999, EMP),
           EMP = ifelse(EMP_F == "m", 100000, EMP)) 
}

get_sect_naics12 <- function(sector, states, years){
  sector.state.year <- lapply(seq_along(years), function(yr)lapply(seq_along(states), function(st) 
    map(sector, possibly(function(y) getCensus(name = "cbp",
                                               vintage = years[yr],
                                               key = mykey,
                                               vars = c("NAICS2012_TTL","EMP", "EMP_F"),
                                               region = "county:*",
                                               regionin = states[st],
                                               naics2012 = y), NA))))
  sector.state <- lapply(seq_along(sector.state.year), 
                         function(x) map(seq_along(states), 
                                         function(y) do.call(rbind,sector.state.year[[x]][[y]][1:length(sector)])))
  sector.yr <- lapply(seq_along(sector.state),
                      function(x) do.call(rbind, sector.state[[x]][1:length(states)]))
  sector.df <- lapply(seq_along(sector.yr), function(x) bind_rows(sector.yr[[x]]) %>% 
                        mutate(year = years[x])) %>% 
    do.call(rbind, .) %>% 
    filter(., county != "999") %>% 
    mutate(., EMP_F = ifelse(is.na(EMP_F), "Z", EMP_F),
           EMP = ifelse(EMP_F == "a", 19, EMP),
           EMP = ifelse(EMP_F == "b", 99, EMP),
           EMP = ifelse(EMP_F == "c", 249, EMP),
           EMP = ifelse(EMP_F == "e", 499, EMP),
           EMP = ifelse(EMP_F == "f", 999, EMP),
           EMP = ifelse(EMP_F == "g", 2499, EMP),
           EMP = ifelse(EMP_F == "h", 4999, EMP),
           EMP = ifelse(EMP_F == "i", 9999, EMP),
           EMP = ifelse(EMP_F == "j", 24999, EMP),
           EMP = ifelse(EMP_F == "k", 49999, EMP),
           EMP = ifelse(EMP_F == "l", 99999, EMP),
           EMP = ifelse(EMP_F == "m", 100000, EMP))
}
# Batch downloads ---------------------------------------------------------


states <- st

years.sic <- as.character(seq(from=1987, to = 1997, by=1))
years.naics97 <- as.character(seq(from=1998, to =2002, by =1))
years.naics02 <- as.character(seq(from=2003, to=2007, by=1))
years.naics07 <- as.character(seq(from=2008, to=2011, by=1))
years.naics12 <- as.character(seq(from=2012, to=2015, by=1))


timber.sic <- c("0811",	"0831","0851","2411","2421","2426","2429","2431","2435","2436","2439","2441",
                  "2448","2449","2451","2452","2491","2493","2499","2611","2621","2631","3553","5031")
timber.naics97 <- c("1131","1132","1133","1153","3211","3212","3219",	"3221","33321","42131")
timber.naics02 <- c("1131","1132","1133","1153","3211","3212","3219",	"3221", "33321","42331")
timber.naics07 <- c("1131","1132","1133","1153","3211","3212","3219",	"3221", "33321","42331")
timber.naics12 <- c("1131","1132","1133","1153","3211","3212","3219",	"3221", "33321","42331")


mine.sic <- c("1011","1021","1031","1041","1044","1061","1081","1094","1099",
              "1221","1222","1231","1241","1311","1321","1381","1382","1389",
              "1411","1422","1423","1429","1442","1446","1455","1459","1474",
              "1475","1479","1481","1499","2911","2951","2952","2992","2999",
              "3211","3221","3229","3231","3241","3251","3253","3255","3259",
              "3261","3262","3263","3264","3269","3271","3272","3273","3274",
              "3275","3281","3291","3292","3295","3296","3297","3299","3312",
              "3313","3315","3316","3317","3321","3322","3324","3325","3331",
              "3334","3339","3341","3351","3353","3354","3355","3356","3357",
              "3363","3364","3365","3366","3369","3399","3462","3463","3466",
              "3532","3533","4612","4922","4923","4924","4925","5051","5052",
              "5171","5172","5983","5984","5989")

mine.naics97 <- c("2111","2121","2122","2123","2131","2212","3241","3271","3272",
                  "3273","3274","3279","3311","3312","3313","3314","3315","3321",
                  "33313","4215","4227","45431","4861","4862","54136")

sic.total.emp <- get_total_sic(states=states, years=years.sic)
naics97.total.emp <- get_total_naics97(states=states, years = years.naics97)
naics02.total.emp <- get_total_naics02(states=states, years=years.naics02)
naics07.total.emp <- get_total_naics07(states=states, years=years.naics07)
naics12.total.emp <- get_total_naics12(states=states, years = years.naics12)

sic.timber.emp <- get_sect_sic(sector=timber.sic,states=states, years=years.sic) %>% 
  group_by(year,state, county) %>% 
  summarise(AllTimber = sum(as.numeric(EMP), na.rm=TRUE))
naics97.timber.emp <- get_sect_naics97(sector=timber.naics97, years = years.naics97, states=states)%>% 
  group_by(year,state, county) %>% 
  summarise(AllTimber = sum(as.numeric(EMP), na.rm=TRUE)) 
naics02.timber.emp <- get_sect_naics02(sector=timber.naics02, years = years.naics02, states=states)%>% 
  group_by(year,state, county) %>% 
  summarise(AllTimber = sum(as.numeric(EMP), na.rm=TRUE))
naics07.timber.emp <- get_sect_naics07(sector=timber.naics07, states=states, years=years.naics07)%>% 
  group_by(year,state, county) %>% 
  summarise(AllTimber = sum(as.numeric(EMP), na.rm=TRUE))
naics12.timber.emp <- get_sect_naics12(sector=timber.naics12, states=states, years=years.naics12) %>% 
  group_by(year,state, county) %>% 
  summarise(AllTimber = sum(as.numeric(EMP), na.rm=TRUE))

sic.mine.emp <- get_sect_sic(sector=mine.sic,states=states, years=years.sic) %>% 
  group_by(year,state, county) %>% 
  summarise(AllMine = sum(as.numeric(EMP), na.rm=TRUE))
naics97.mine.emp <- get_sect_naics97(sector=mine.naics97, years = years.naics97, states=states)%>% 
  group_by(year,state, county) %>% 
  summarise(AllMine = sum(as.numeric(EMP), na.rm=TRUE)) 
naics02.mine.emp <- get_sect_naics02(sector=mine.naics97, years = years.naics02, states=states)%>% 
  group_by(year,state, county) %>% 
  summarise(AllMine = sum(as.numeric(EMP), na.rm=TRUE)) 
naics07.mine.emp <- get_sect_naics07(sector=mine.naics97, years = years.naics07, states=states)%>% 
  group_by(year,state, county) %>% 
  summarise(AllMine = sum(as.numeric(EMP), na.rm=TRUE)) 
naics12.mine.emp <- get_sect_naics12(sector=mine.naics97, years = years.naics12, states=states)%>% 
  group_by(year,state, county) %>% 
  summarise(AllMine = sum(as.numeric(EMP), na.rm=TRUE)) 


total.emp <- rbind(sic.total.emp, naics97.total.emp, naics02.total.emp, naics07.total.emp, naics12.total.emp)
timber.emp <- rbind(sic.timber.emp, naics97.timber.emp, naics02.timber.emp, naics07.timber.emp, naics12.timber.emp)
mine.emp <- rbind(sic.mine.emp, naics97.mine.emp, naics02.mine.emp, naics07.mine.emp, naics12.mine.emp)

timber.join <- left_join(total.emp, timber.emp)
timber.join$AllTimber <- ifelse(is.na(timber.join$AllTimber), 0, timber.join$AllTimber) 
timber.join$propTimber <- timber.join$AllTimber/as.numeric(timber.join$TOTALEMP)
timber.join$propTimber <- round(ifelse(timber.join$propTimber > 1, 1, timber.join$propTimber), digits=4)

mine.join <- left_join(total.emp, mine.emp)
mine.join$AllMine <- ifelse(is.na(mine.join$AllMine), 0, mine.join$AllMine) 
mine.join$propMine <- mine.join$AllMine/as.numeric(mine.join$TOTALEMP)
mine.join$propMine <- round(ifelse(mine.join$propMine > 1, 1, mine.join$propMine), digits=4)



obj <- ls()
rem.obj <- obj[obj!="mine.join"]
rm(list=rem.obj)


prj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

all.years <- seq(from=1987, to=2015, by=1)
yr <- all.years
library(raster) #masks select so don't load when running the download scripts
get_geog_rst <- function(yr,fips, emp){
  if(yr < 2000){
    cty.shp <- counties(fips,year=1990, cb = TRUE)%>% as(., "sf") %>% st_transform(., crs=prj)  
    yr.data <- emp %>% 
      filter(., year == as.character(yr)) %>% 
      left_join(., cty.shp, by=c("state"="ST", "county"="CO")) 
    }
  if(yr >= 2000 & yr <2010){
    cty.shp <- counties(fips,year=2000)%>% as(., "sf") %>% st_transform(., crs=prj)
    yr.data <- emp %>% 
      filter(., year == as.character(yr)) %>% 
      left_join(., cty.shp, by=c("state"="STATEFP00", "county"="COUNTYFP00"))
    }
  cty.shp <- counties(fips, year=2010)%>% as(., "sf") %>% st_transform(., crs=prj)
  yr.data <- emp %>% 
    filter(., year == as.character(yr)) %>% 
    left_join(., cty.shp, by=c("state"="STATEFP10", "county"="COUNTYFP10")) %>% 
    st_sf(.)
  b <- fasterize::raster(extent(yr.data), res=270)
  val.name <-colnames(emp)[6]
  rst.out <- fasterize::fasterize(yr.data, b, field=val.name, fun="max")
}

excluded.states <- c("AS","GU","MP","PR","UM","VI","AK", "HI","DC")
states <- unique(subset(fips_codes, !(fips_codes$state %in% excluded.states)))
st.fips <- unique(states$state_code)
geo.yr <- sample(all.years,7)
timber.rst1 <- map(all.years[1:7], function(x) get_geog_rst(yr=x, fips=st.fips, emp=timber.join))
timber.rst2 <- map(all.years[8:14], function(x) get_geog_rst(yr=x, fips=st.fips, emp=timber.join))
timber.rst3 <- map(all.years[15:21], function(x) get_geog_rst(yr=x, fips=st.fips, emp=timber.join))
timber.rst4 <- map(all.years[22:29], function(x) get_geog_rst(yr=x, fips=st.fips, emp=timber.join))

tim.brick1 <- do.call(brick, timber.rst1)
rm(timber.rst1)
tim.brick2 <- do.call(brick, timber.rst2)
rm(timber.rst2)
tim.brick3 <- do.call(brick, timber.rst3)
rm(timber.rst3)
tim.brick4 <- do.call(brick, timber.rst4)

tim.brick <- stack(tim.brick1, tim.brick2, tim.brick3, tim.brick4)
names(tim.brick) <- all.years
writeRaster(tim.brick, filename = names(tim.brick),bylayer=TRUE, format = "GTiff")

mine.rst1 <- map(all.years[1:7], function(x) get_geog_rst(yr=x, fips=st.fips, emp=mine.join))
mine.brick1 <- do.call(brick, mine.rst1)
rm(mine.rst1)

mine.rst2 <- map(all.years[8:14], function(x) get_geog_rst(yr=x, fips=st.fips, emp=mine.join))
mine.brick2 <- do.call(brick, mine.rst2)
rm(mine.rst2)

mine.rst3 <- map(all.years[15:21], function(x) get_geog_rst(yr=x, fips=st.fips, emp=mine.join))
mine.brick3 <- do.call(brick, mine.rst3)
rm(mine.rst3)

mine.rst4 <- map(all.years[22:29], function(x) get_geog_rst(yr=x, fips=st.fips, emp=mine.join))
mine.brick4 <- do.call(brick, mine.rst4)

mine.brick <- stack(mine.brick1, mine.brick2, mine.brick3, mine.brick4)
names(mine.brick) <- all.years
writeRaster(mine.brick, filename = names(mine.brick),bylayer=TRUE, format = "GTiff")


