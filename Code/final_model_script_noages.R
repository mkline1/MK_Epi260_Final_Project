library(tidyverse)
library(ggplot2)
library(odin)
library(gridExtra)

no_age_SEIR <- odin({
  #Q1 model script
  
  #Differential equations of states in 'deriv' form
  deriv(S_v) <- p_b*V_b*S + p_a*alpha*Sp*S + w_vi*R_v - lambda_v*S_v - w_v*S_v
  deriv(S) <- w_v*S_v + w_i*R - p_b*V_b*S - p_a*alpha*Sp*S - lambda*S
  deriv(E_v) <- lambda_v*S_v + p_b*V_b*E + p_a*alpha*Sp*E - w_v*E_v - sigma*E_v 
  deriv(E) <- lambda*S + w_v*E_v - p_b*V_b*E - p_a*alpha*Sp*E - sigma*E
  deriv(I_v) <- sigma*E_v - gamma*IFR*(1-Ve_d)*I_v - gamma*(1-(IFR*(1-Ve_d)))*I_v
  deriv(I) <- sigma*E - gamma*IFR*I - gamma*(1-IFR)*I
  deriv(D_v) <- gamma*IFR*(1-Ve_d)*I_v
  deriv(D) <- gamma*IFR*I
  deriv(R_v) <- gamma*(1-(IFR*(1-Ve_d)))*I_v + p_b*V_b*R + p_a*alpha*(1-Se)*R - w_vi*R_v - w_v*R_v
  deriv(R) <- gamma*(1-IFR)*I + w_v*R_v - p_b*V_b*R - p_a*alpha*(1-Se)*R - w_i*R
   
  #Total population
  N <- S_v + S + E_v + E + I_v + I + R_v + R #should this have D_v and D?
  lambda <- mu*c*(I+I_v)/N
  lambda_v <- (1-Ve_i)*mu*c*(I+I_v)/N
  
  
  ## Initial conditions
  initial(S_v) <- S_v_ini
  initial(S) <- S_ini
  initial(E_v) <- E_v_ini
  initial(E) <- E_ini
  initial(I_v) <- I_v_ini
  initial(I) <- I_ini
  initial(D_v) <- D_v_ini
  initial(D) <- D_ini
  initial(R_v) <- R_v_ini
  initial(R) <- R_ini
  
  ## parameters
  ###initial states (to be user defined in each run)
  S_v_ini <- user()
  S_ini <- user()
  E_v_ini <- user()
  E_ini <- user()
  I_v_ini <- user()
  I_ini <- user()
  D_v_ini <- user()
  D_ini <- user()
  R_v_ini <- user()
  R_ini <- user()
  
  ###parameter values (to be user defined in each run)
  V_b <- user()
  p_b <- user()
  p_a <- user()
  alpha <- user()
  Sp <- user()
  w_vi <- user()
  w_i <- user()
  w_v <- user()
  sigma <- user()
  gamma <- user()
  IFR <- user()
  Se <- user()
  mu <- user()
  Ve_d <- user()
  Ve_i <- user()
  c <- user()
  
})

#model seeding
#load in info from data calculations that i will need to for seeding with no age structure
boosted_prop <- readRDS('/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/no_age_structure/boost_prop_start.RDS')
prop_e_start <- readRDS('/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/no_age_structure/prop_e_start.RDS')
prop_i_start <- readRDS('/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/no_age_structure/prop_i_start.RDS')
prop_r_start <- readRDS('/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/no_age_structure/prop_r_start.RDS')

start_N <- 10^6
n_vax <- boosted_prop*start_N
E_v_ini <- n_vax * prop_e_start
I_v_ini <- n_vax*prop_i_start
R_v_ini <- n_vax*prop_r_start
S_v_ini <- n_vax - E_v_ini - I_v_ini - R_v_ini
novax <- start_N - n_vax
E_ini <- novax*prop_e_start
I_ini <- novax*prop_i_start 
R_ini <- novax*prop_r_start
S_ini <- novax - E_ini - I_ini - R_ini







scen1 <- no_age_SEIR$new(S_v_ini = S_v_ini,
                        S_ini = S_ini,
                        E_v_ini = E_v_ini,
                        E_ini = E_ini,
                        I_v_ini = I_v_ini,
                        I_ini = I_ini,
                        D_v_ini = 0,
                        D_ini = 0,
                        R_v_ini = R_v_ini,
                        R_ini = R_ini,
                        V_b = 0,
                        p_b = 0,
                        p_a = 0,
                        alpha = 0,
                        Sp = 0,
                        w_vi = 0.008333333,
                        w_i = 0.008333333,
                        w_v = 0.004166667,
                        sigma = 0.3333333,
                        gamma = 0.2,
                        IFR = 0.025,
                        Se = 0,
                        mu = 0.7,
                        Ve_d = 0.87,
                        Ve_i = 0.3,
                        c = 13)


t = seq(0,1825)
waning_noab_noannualvax <- as_tibble(scen1$run(t)) |>
  pivot_longer(cols = 2:11) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + 
  geom_line() +
  xlab("Time Step") + 
  ggtitle("No Antibody Testing, No Additional Vaccination") +
  theme_minimal()
scen1_total_deaths <- as_tibble(scen1$run(t)) |> 
  mutate(total_deaths = D + D_v) |>
  pull(total_deaths)
scen1_endpoint_deaths <- scen1_total_deaths[length(scen1_total_deaths)]


yearlyab <- no_age_SEIR$new(S_v_ini = S_v_ini,
                            S_ini = S_ini,
                            E_v_ini = E_v_ini,
                            E_ini = E_ini,
                            I_v_ini = I_v_ini,
                            I_ini = I_ini,
                            D_v_ini = 0,
                            D_ini = 0,
                            R_v_ini = R_v_ini,
                            R_ini = R_ini,
                            V_b = 0,
                            p_b = 0,
                            p_a = 1,
                            alpha = 0.002739726,
                            Sp = 0.9,
                            w_vi = 0.008333333,
                            w_i = 0.008333333,
                            w_v = 0.004166667,
                            sigma = 0.3333333,
                            gamma = 0.2,
                            IFR = 0.025,
                            Se = 0.9,
                            mu = 0.7,
                            Ve_d = 0.87,
                            Ve_i = 0.3,
                            c = 13)

t = seq(0,1825)
yearlyab_testing <- as_tibble(yearlyab$run(t)) |>
  pivot_longer(cols = 2:11) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + 
  geom_line() +
  xlab("Time Step") + 
  labs(title = "Yearly Home Antibody Testing Only",
       subtitle = "90% Sensitivity, 90% Specificity, 100% Participation") +
  theme_minimal()

yearlyab_total_deaths <- as_tibble(yearlyab$run(t)) |> 
  mutate(total_deaths = D + D_v) |>
  pull(total_deaths)
yearlyab_endpoint_deaths <- yearlyab_total_deaths[length(yearlyab_total_deaths)]


##scenario 2, waning but additional yearly vaccination 
scen2 <- no_age_SEIR$new(S_v_ini = S_v_ini,
                         S_ini = S_ini,
                         E_v_ini = E_v_ini,
                         E_ini = E_ini,
                         I_v_ini = I_v_ini,
                         I_ini = I_ini,
                         D_v_ini = 0,
                         D_ini = 0,
                         R_v_ini = R_v_ini,
                         R_ini = R_ini,
                       V_b = 0.002739726,
                       p_b = 1,
                       p_a = 0,
                       alpha = 0,
                       Sp = 0,
                       w_vi = 0.008333333,
                       w_i = 0.008333333,
                       w_v = 0.004166667,
                       sigma = 0.3333333,
                       gamma = 0.2,
                       IFR = 0.025,
                       Se = 0,
                       mu = 0.7,
                       Ve_d = 0.87,
                       Ve_i = 0.3,
                       c = 13)

t = seq(0,1825)
waning_noab <- as_tibble(scen2$run(t)) |>
  pivot_longer(cols = 2:11) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + 
  geom_line() +
  xlab("Time Step") + 
  ggtitle("No Antibody Testing, Yearly Vaccination") +
  theme_minimal()

scen2_total_deaths <- as_tibble(scen2$run(t)) |> 
  mutate(total_deaths = D + D_v) |>
  pull(total_deaths)
scen2_endpoint_deaths <- scen2_total_deaths[length(scen2_total_deaths)]

##scenario 3: no yearly vaccination, but use of a home antibody test (90% sensitivity/specificity, 100% uptake)
scen3 <- no_age_SEIR$new(S_v_ini = S_v_ini,
                         S_ini = S_ini,
                         E_v_ini = E_v_ini,
                         E_ini = E_ini,
                         I_v_ini = I_v_ini,
                         I_ini = I_ini,
                         D_v_ini = 0,
                         D_ini = 0,
                         R_v_ini = R_v_ini,
                         R_ini = R_ini,
                         V_b = 0,
                         p_b = 0,
                         p_a = 1,
                         alpha = 0.03333333,
                         Sp = 0.9,
                         w_vi = 0.008333333,
                         w_i = 0.008333333,
                         w_v = 0.004166667,
                         sigma = 0.3333333,
                         gamma = 0.2,
                         IFR = 0.025,
                         Se = 0.9,
                         mu = 0.7,
                         Ve_d = 0.87,
                         Ve_i = 0.3,
                         c = 13)


t = seq(0,1825)
waning_ab_noyearly <- as_tibble(scen3$run(t)) |>
  pivot_longer(cols = 2:11) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + 
  geom_line() +
  xlab("Time Step") + 
  labs(title = "Monthly Home Antibody Testing, No Yearly Vaccination",
       subtitle = "90% Sensitivity, 90% Specificity, 100% Participation") +
  theme_minimal()

scen3_total_deaths <- as_tibble(scen3$run(t)) |> 
  mutate(total_deaths = D + D_v) |>
  pull(total_deaths)
scen3_endpoint_deaths <- scen3_total_deaths[length(scen3_total_deaths)]


## Scenario 4, both yearly vaccination and antibody testing
scen4 <- no_age_SEIR$new(S_v_ini = S_v_ini,
                         S_ini = S_ini,
                         E_v_ini = E_v_ini,
                         E_ini = E_ini,
                         I_v_ini = I_v_ini,
                         I_ini = I_ini,
                         D_v_ini = 0,
                         D_ini = 0,
                         R_v_ini = R_v_ini,
                         R_ini = R_ini,
                         V_b = 0.002739726,
                         p_b = 1,
                         p_a = 1,
                         alpha = 0.03333333,
                         Sp = 0.9,
                         w_vi = 0.008333333,
                         w_i = 0.008333333,
                         w_v = 0.004166667,
                         sigma = 0.3333333,
                         gamma = 0.2,
                         IFR = 0.025,
                         Se = 0.9,
                         mu = 0.7,
                         Ve_d = 0.87,
                         Ve_i = 0.3,
                         c = 13)
t = seq(0,1825)
waning_ab_yearly <- as_tibble(scen4$run(t)) |>
  pivot_longer(cols = 2:11) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + 
  geom_line() +
  xlab("Time Step") + 
  labs(title = "Monthly Home Antibody Testing, Additionally Yearly Vaccination",
       subtitle = "90% Sensitivity, 90% Specificity, 100% Participation") +
  theme_minimal()

scen4_total_deaths <- as_tibble(scen4$run(t)) |> 
  mutate(total_deaths = D + D_v) |>
  pull(total_deaths)
scen4_endpoint_deaths <- scen4_total_deaths[length(scen4_total_deaths)]


allscensplot <- grid.arrange(waning_noab_noannualvax, yearlyab_testing, waning_noab,
             waning_ab_noyearly, waning_ab_yearly, ncol = 3)
#ggsave("noage_scenarios.png", allscensplot, device = 'png', width = 20, height = 10)

total_deaths <- data.frame(deaths = rbind(scen1_endpoint_deaths, yearlyab_endpoint_deaths,scen2_endpoint_deaths,
      scen3_endpoint_deaths, scen4_endpoint_deaths), scenario = c(1:5))
final_death_plot <- total_deaths |> 
  ggplot(aes(x = scenario, y = deaths, fill = as.factor(scenario))) + geom_bar(stat = "identity") +
  theme_minimal() + scale_fill_discrete(name = "Scenario",
                                        labels = c("No additional vaccination", "Yearly Antibody Test Only",
                                                   "Yearly Vaccination",
                                                   "Monthly use of home antibody test only",
                                                   "Monthly use of home antibody test + yearly vaccination")) +
  labs(title = "Total Deaths after 5 years")
#ggsave("noagefinal_death_plot.png", final_death_plot, device = "png", width = 10, height = 7)

deathsscenset1 <- grid.arrange(deaths_over_time,final_death_plot, ncol = 2)
#ggsave("noages_scen_set1_deaths.png", deathsscenset1, width =15)




#Plot total deaths over time in all 4 scenarios
deaths_df <- data.frame("t" = as_tibble(scen1$run(t))[,1],
                        "scenario1" = rowSums(as_tibble(scen1$run(t))[,8:9]), 
                        "scenario2" = rowSums(as_tibble(yearlyab$run(t))[,8:9]),
                        "scenario3" = rowSums(as_tibble(scen2$run(t))[,8:9]),
                        "scenario4" = rowSums(as_tibble(scen3$run(t))[,8:9]),
                        "scenario5" = rowSums(as_tibble(scen4$run(t))[,8:9]))
deaths_over_time <- deaths_df |> 
  pivot_longer(cols = 2:6) |>
  ggplot(aes(x = t, y = value, color = name)) + 
  geom_line() + 
  theme_minimal() +
  ylab("Deaths") + 
  scale_color_discrete(name = "Scenario",
                      labels = c("No additional vaccination", "Yearly Antibody Test Only",
                                 "Yearly Vaccination",
                                 "Monthly use of home antibody test only",
                                 "Monthly use of home antibody test + yearly vaccination")) +
  ggtitle("Deaths Over Time")

#ggsave("noages_deaths_over_time.png",deaths_over_time, device = 'png')


#second set of scenarios varying the proportion of people participating in the antibody testing
propquarter <- no_age_SEIR$new(S_v_ini = S_v_ini,
                               S_ini = S_ini,
                               E_v_ini = E_v_ini,
                               E_ini = E_ini,
                               I_v_ini = I_v_ini,
                               I_ini = I_ini,
                               D_v_ini = 0,
                               D_ini = 0,
                               R_v_ini = R_v_ini,
                               R_ini = R_ini,
                         V_b = 0.002739726,
                         p_b = 1,
                         p_a = 0.25,
                         alpha = 0.03333333,
                         Sp = 0.9,
                         w_vi = 0.008333333,
                         w_i = 0.008333333,
                         w_v = 0.004166667,
                         sigma = 0.3333333,
                         gamma = 0.2,
                         IFR = 0.025,
                         Se = 0.9,
                         mu = 0.7,
                         Ve_d = 0.87,
                         Ve_i = 0.3,
                         c = 13)

t = seq(0,1825)
quarter_testing <- as_tibble(propquarter$run(t)) |>
  pivot_longer(cols = 2:11) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + 
  geom_line() +
  xlab("Time Step") + 
  labs(title = "Monthly Home Antibody Testing, Additionally Yearly Vaccination",
       subtitle = "90% Sensitivity, 90% Specificity, 25% Participation") +
  theme_minimal()

propquarter_total_deaths <- as_tibble(propquarter$run(t)) |> 
  mutate(total_deaths = D + D_v) |>
  pull(total_deaths)
propquarter_endpoint_deaths <- propquarter_total_deaths[length(propquarter_total_deaths)]


prophalf <- no_age_SEIR$new(S_v_ini = S_v_ini,
                            S_ini = S_ini,
                            E_v_ini = E_v_ini,
                            E_ini = E_ini,
                            I_v_ini = I_v_ini,
                            I_ini = I_ini,
                            D_v_ini = 0,
                            D_ini = 0,
                            R_v_ini = R_v_ini,
                            R_ini = R_ini,
                               V_b = 0.002739726,
                            p_b = 1,
                               p_a = 0.50,
                               alpha = 0.03333333,
                               Sp = 0.9,
                               w_vi = 0.008333333,
                               w_i = 0.008333333,
                               w_v = 0.004166667,
                               sigma = 0.3333333,
                               gamma = 0.2,
                               IFR = 0.025,
                               Se = 0.9,
                               mu = 0.7,
                               Ve_d = 0.87,
                               Ve_i = 0.3,
                               c = 13)

t = seq(0,1825)
half_testing <- as_tibble(prophalf$run(t)) |>
  pivot_longer(cols = 2:11) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + 
  geom_line() +
  xlab("Time Step") + 
  labs(title = "Monthly Home Antibody Testing, Additionally Yearly Vaccination",
       subtitle = "90% Sensitivity, 90% Specificity, 50% Participation") +
  theme_minimal()

prophalf_total_deaths <- as_tibble(prophalf$run(t)) |> 
  mutate(total_deaths = D + D_v) |>
  pull(total_deaths)
prophalf_endpoint_deaths <- prophalf_total_deaths[length(prophalf_total_deaths)]


propthreequart <- no_age_SEIR$new(S_v_ini = S_v_ini,
                                  S_ini = S_ini,
                                  E_v_ini = E_v_ini,
                                  E_ini = E_ini,
                                  I_v_ini = I_v_ini,
                                  I_ini = I_ini,
                                  D_v_ini = 0,
                                  D_ini = 0,
                                  R_v_ini = R_v_ini,
                                  R_ini = R_ini,
                            V_b = 0.002739726,
                            p_b = 1,
                            p_a = 0.75,
                            alpha = 0.03333333,
                            Sp = 0.9,
                            w_vi = 0.008333333,
                            w_i = 0.008333333,
                            w_v = 0.004166667,
                            sigma = 0.3333333,
                            gamma = 0.2,
                            IFR = 0.025,
                            Se = 0.9,
                            mu = 0.7,
                            Ve_d = 0.87,
                            Ve_i = 0.3,
                            c = 13)

t = seq(0,1825)
threequart_testing <- as_tibble(propthreequart$run(t)) |>
  pivot_longer(cols = 2:11) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + 
  geom_line() +
  xlab("Time Step") + 
  labs(title = "Monthly Home Antibody Testing, Additionally Yearly Vaccination",
       subtitle = "90% Sensitivity, 90% Specificity, 75% Participation") +
  theme_minimal()

propthreequart_total_deaths <- as_tibble(propthreequart$run(t)) |> 
  mutate(total_deaths = D + D_v) |>
  pull(total_deaths)
propthreequart_endpoint_deaths <- propthreequart_total_deaths[length(propthreequart_total_deaths)]



diff_partic <- grid.arrange(yearlyab_testing,waning_noab, quarter_testing, 
                            half_testing, threequart_testing, waning_ab_yearly, ncol =2)

ggsave("noage_ab_participation.png", diff_partic, device = "png", width = 13, height = 10)

#Plot total deaths over time in all 4 scenarios
deaths_df2 <- data.frame("t" = as_tibble(scen2$run(t))[,1],
                        "No_Ab_Yearly_Vax" = rowSums(as_tibble(scen2$run(t))[,8:9]), 
                        "Yearly_Ab_only"  = rowSums(as_tibble(yearlyab$run(t))[,8:9]),
                        "25perc_Ab_and_Yearly_Vax" = rowSums(as_tibble(propquarter$run(t))[,8:9]),
                        "50perc_Ab_and_Yearly_Vax" = rowSums(as_tibble(prophalf$run(t))[,8:9]),
                        "75perc_Ab_and_Yearly_Vax" = rowSums(as_tibble(propthreequart$run(t))[,8:9]),
                        "100perc_Ab_and_Yearly_Vax" = rowSums(as_tibble(scen4$run(t))[,8:9]))
deaths_over_time2 <- deaths_df2 |> 
  pivot_longer(cols = 2:7) |>
  ggplot(aes(x = t, y = value, color = name)) + 
  geom_line() + 
  theme_minimal() +
  xlab("Deaths") +
  scale_color_discrete(name = "Scenario",
                       labels = c("Yearly Vaccination Only", "Monthly Antibody Testing, 100% Participation ",
                                  "Monthly Antibody Testing, 25% Participation",
                                  "Monthly Antibody Testing, 50% Participation",
                                  "Monthly Antibody Testing, 75% Participation",
                                  "Yearly Antibody Testing, 100% Participation"
                                  )) +
  ggtitle("Deaths Over Time")


#ggsave("noage_deaths_over_time_participation.png", deaths_over_time2, device = "png")
final_deaths <- deaths_df2[nrow(deaths_df2),] |>
  pivot_longer(cols = 2:7) |>
  ggplot(aes(x = reorder(name, -value), y = value, fill = name)) +
  geom_bar(stat = "identity") + 
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  xlab("Scenario") +
  ylab("Deaths") +
  labs(title = "Final Deaths after 5 Years", subtitle = "Differing Antibody Test Participation",) +
  scale_fill_discrete(name = "Scenario",
                       labels = c("Yearly Vaccination Only", "Monthly Antibody Testing, 100% Participation ",
                                  "Monthly Antibody Testing, 25% Participation",
                                  "Monthly Antibody Testing, 50% Participation",
                                  "Monthly Antibody Testing, 75% Participation",
                                  "Yearly Antibody Testing, 100% Participation"))
  
#ggsave("noage_finaldeaths_participation.png",final_deaths)









