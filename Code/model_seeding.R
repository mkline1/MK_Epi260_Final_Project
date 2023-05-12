library(tidyverse)
library(lubridate)
library(MMWRweek)

#figuring out infected vs recovered:
cdcdat <- read_csv("/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/data_table_for_weekly_case_trends__the_united_states.csv",
                   skip = 2)
#fix date column using lubridate
cdcdat <- cdcdat |> mutate(Date= mdy(Date))

#using latent period of 3 days, recovery period of 5 days, and 120 days for the time people spend in the R compartment

#to get infected and exposed, go one week back and imagine starting there
# the infected as the people infected 1 week ago 
#the people exposed during that time are the people infected in the current week (5+3 > 7)


###infected
#the number of people who will be in the infected compartment are the people within the last 5 days
#weekly data, so take the last week
currently_infected <- cdcdat$`Weekly Cases`[2]

#the number exposed is the people 
currently_exposed <- cdcdat$`Weekly Cases`[1]



#figure out when 120 days ago was (those are the people who will be non-susceptible)
startdate <- cdcdat$Date[2] - 120
#figure out which date in the dataframe is closest to this date
i_from <- which(abs(cdcdat$Date-startdate) == min(abs(cdcdat$Date - startdate)))


#sum cases starting with the start date
recovered_cases <- cdcdat |> 
  filter(Date >= cdcdat$Date[i_from], Date < cdcdat$Date[2]) |> 
  summarize(total_cases  = sum(`Weekly Cases`)) |> pull()

#these are absolute numbers, figure out proportions assuming a population size of 330x10^6
total_pop <- 3.3e8
prop_i <- currently_infected/total_pop
prop_e <- currently_exposed/total_pop
prop_R <- recovered_cases/total_pop

#save these to multiply by population size for model seeding
# saveRDS(prop_i, "prop_i_start.RDS")
# saveRDS(prop_e, "prop_e_start.RDS")
# saveRDS(prop_R, "prop_R_start.RDS")


#now figure out proportion vaccinated 
#read in CDC vaccination data
vaxdat <- read_csv("/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/covid19_vaccinations_in_the_united_states.csv",
                   skip = 4)
usvax <- vaxdat[1,]


#we want total residents with updated bivalent booster in the US
total_boost <- usvax$`Residents with an updated (bivalent) booster dose`
boost_prop <- total_boost/total_pop
#saveRDS(boost_prop, "boost_prop_start.RDS")



#attempting with age structure
#https://data.cdc.gov/Public-Health-Surveillance/Rates-of-COVID-19-Cases-or-Deaths-by-Age-Group-and/3rge-nu2a
ac_dat <- read_csv("/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/Rates_of_COVID-19_Cases_or_Deaths_by_Age_Group_and_Vaccination_Status.csv")

#this data goes through the end of september 2022
#use MMRW week 37 for infected, 38 for exposed, and then use all the way back through week 21 for recovered
ac_dat_filt <- ac_dat |> filter(`MMWR week` >= 202221)

infected_ages <- ac_dat_filt |> filter(`MMWR week` == 202237) |>
  filter(`Age group` != "all_ages_adj") |>
  select(month, `MMWR week`, `Age group`, `Vaccinated with outcome`, `Fully vaccinated population`, `Unvaccinated with outcome`, `Unvaccinated population`) |>
  mutate(i_vax_prop = `Vaccinated with outcome`/`Fully vaccinated population`,
         i_unvax_prop = `Unvaccinated with outcome` / `Unvaccinated population`)
#need to make it match my 10 age groups; 
#<10 is 5-11
#10-20 is 12-17
#20-30 is 18-29
#30-40 is 30-49
#40-50 is 30-49
#50-60 is 50-60
#60-70 is 65-79
#70-80 is 65-79
#80+ is 80+

i_vax_ages <- c(infected_ages$i_vax_prop[1],infected_ages$i_vax_prop[2], infected_ages$i_vax_prop[3], infected_ages$i_vax_prop[4],
                infected_ages$i_vax_prop[4], infected_ages$i_vax_prop[5], infected_ages$i_vax_prop[6], infected_ages$i_vax_prop[6],
                infected_ages$i_vax_prop[7])
#saveRDS(i_vax_ages, "i_vax_ages_prop.RDS")
#multiply these proportions by the number of people in each age group that are vaccinated
i_unvax_ages <- c(infected_ages$i_unvax_prop[1],infected_ages$i_unvax_prop[2], infected_ages$i_unvax_prop[3], infected_ages$i_unvax_prop[4],
                  infected_ages$i_unvax_prop[4], infected_ages$i_unvax_prop[5], infected_ages$i_unvax_prop[6], infected_ages$i_unvax_prop[6],
                  infected_ages$i_unvax_prop[7])
#saveRDS(i_unvax_ages, "i_unvax_ages_prop.RDS")

exposed_ages <- ac_dat_filt |> filter(`MMWR week` == 202238) |>
  filter(`Age group` != "all_ages_adj") |>
  select(month, `MMWR week`, `Age group`, `Vaccinated with outcome`, `Fully vaccinated population`, `Unvaccinated with outcome`, `Unvaccinated population`) |>
  mutate(e_vax_prop = `Vaccinated with outcome`/`Fully vaccinated population`,
         e_unvax_prop = `Unvaccinated with outcome` / `Unvaccinated population`)

e_vax_ages <- c(exposed_ages$e_vax_prop[1],exposed_ages$e_vax_prop[2], exposed_ages$e_vax_prop[3], exposed_ages$e_vax_prop[4],
                exposed_ages$e_vax_prop[4], exposed_ages$e_vax_prop[5], exposed_ages$e_vax_prop[6], exposed_ages$e_vax_prop[6],
                exposed_ages$e_vax_prop[7])
#saveRDS(e_vax_ages, "e_vax_ages_prop.RDS")

e_unvax_ages <- c(exposed_ages$e_unvax_prop[1],exposed_ages$e_unvax_prop[2], exposed_ages$e_unvax_prop[3], exposed_ages$e_unvax_prop[4],
                exposed_ages$e_unvax_prop[4], exposed_ages$e_unvax_prop[5], exposed_ages$e_unvax_prop[6], exposed_ages$e_unvax_prop[6],
                exposed_ages$e_unvax_prop[7])
#saveRDS(e_unvax_ages, "e_unvax_ages_prop.RDS")

#those are the proportions of vaccinated and unvaccinated people in each age group that are infected and exposed
#now do recovered
recovered_ages <- ac_dat_filt |>
  filter(`MMWR week` < 202237) |>
  filter(`Age group` != "all_ages_adj") |>
  select(month, `MMWR week`, `Age group`, `Vaccinated with outcome`, `Fully vaccinated population`, `Unvaccinated with outcome`, `Unvaccinated population`) |>
  group_by(`Age group`) |>
  summarize(Vaxxed_outcome = sum(`Vaccinated with outcome`), Vaxxed_pop = mean(`Fully vaccinated population`),
            Unvaxxed_outcome = sum(`Unvaccinated with outcome`), Unvaxxed_pop = mean(`Unvaccinated population`)) |>
  mutate(r_vax_prop = Vaxxed_outcome/Vaxxed_pop,
         r_unvax_prop = Unvaxxed_outcome/Unvaxxed_pop)
r_vax_ages <- c(recovered_ages$r_vax_prop[1],recovered_ages$r_vax_prop[2], recovered_ages$r_vax_prop[3], recovered_ages$r_vax_prop[4],
                recovered_ages$r_vax_prop[4], recovered_ages$r_vax_prop[5], recovered_ages$r_vax_prop[6], recovered_ages$r_vax_prop[6],
                recovered_ages$r_vax_prop[7])
r_unvax_ages <- c(recovered_ages$r_unvax_prop[1],recovered_ages$r_unvax_prop[2], recovered_ages$r_unvax_prop[3], recovered_ages$r_unvax_prop[4],
                recovered_ages$r_unvax_prop[4], recovered_ages$r_unvax_prop[5], recovered_ages$r_unvax_prop[6], recovered_ages$r_unvax_prop[6],
                recovered_ages$r_unvax_prop[7])

#saveRDS(r_vax_ages, "r_vax_ages_prop.RDS")
#saveRDS(r_unvax_ages, "r_unvax_ages_prop_take2.RDS")

#now I just need the proportion of each age group that is vaccinated
vaxdatages <- read_csv("/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/Archive__COVID-19_Vaccination_and_Case_Trends_by_Age_Group__United_States.csv")
#as of october 2022, so aligns with the other dataset
#filter to latest date data was collected which was 10/10/2022
vaxdatages <- vaxdatages |> filter(`Date Administered`  == '10/10/2022 12:00:00 AM')
#0-10 average the <2, 2-4, and 5-11 (1-3)
#10-20 use 12-17 (4)
#20-30 use 18-24 (5)
#30-40 use 25-49 (6)
#40-50 use 25-49 (6)
#50-60 use 50-64 (7)
#60-70 use 65+ (8)
#70-80 use 65+ (8)
#80+ use use 65+ (8)
vaxpropages <- c(mean(c(vaxdatages$Series_Complete_Pop_pct_agegroup[1],vaxdatages$Series_Complete_Pop_pct_agegroup[2],vaxdatages$Series_Complete_Pop_pct_agegroup[3])),
                 vaxdatages$Series_Complete_Pop_pct_agegroup[4], vaxdatages$Series_Complete_Pop_pct_agegroup[5], vaxdatages$Series_Complete_Pop_pct_agegroup[6],
                 vaxdatages$Series_Complete_Pop_pct_agegroup[6], vaxdatages$Series_Complete_Pop_pct_agegroup[7], vaxdatages$Series_Complete_Pop_pct_agegroup[8],
                 vaxdatages$Series_Complete_Pop_pct_agegroup[8], vaxdatages$Series_Complete_Pop_pct_agegroup[8]
                 
)
#saveRDS(vaxpropages, "prop_ages_vax.RDS")
                 



















  
 