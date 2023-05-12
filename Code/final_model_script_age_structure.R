
library(tidyverse)
library(ggplot2)
library(odin)
library(gridExtra)

age_structured <- odin({
  #Differential equations of states in 'deriv' form
  deriv(S_v[]) <- p_b[i]*V_b[i]*S[i] + p_a[i]*alpha[i]*Sp*S[i] + w_vi*R_v[i] - lambda_v[i]*S_v[i] - w_v*S_v[i]
  deriv(S[]) <- w_v*S_v[i] + w_i*R[i] - p_b[i]*V_b[i]*S[i] - p_a[i]*alpha[i]*Sp*S[i] - lambda[i]*S[i]
  deriv(E_v[]) <- lambda_v[i]*S_v[i] + p_b[i]*V_b[i]*E[i] + p_a[i]*alpha[i]*Sp*E[i] - w_v*E_v[i] - sigma*E_v[i] 
  deriv(E[]) <- lambda[i]*S[i] + w_v*E_v[i] - p_b[i]*V_b[i]*E[i] - p_a[i]*alpha[i]*Sp*E[i] - sigma*E[i]
  deriv(I_v[]) <- sigma*E_v[i] - gamma*IFR[i]*(1-Ve_d)*I_v[i] - gamma*(1-(IFR[i]*(1-Ve_d)))*I_v[i]
  deriv(I[]) <- sigma*E[i] - gamma*IFR[i]*I[i] - gamma*(1-IFR[i])*I[i]
  deriv(D_v[]) <- gamma*IFR[i]*(1-Ve_d)*I_v[i]
  deriv(D[]) <- gamma*IFR[i]*I[i]
  deriv(R_v[]) <- gamma*(1-(IFR[i]*(1-Ve_d)))*I_v[i] + p_b[i]*V_b[i]*R[i] + p_a[i]*alpha[i]*(1-Se)*R[i] - w_vi*R_v[i] - w_v*R_v[i]
  deriv(R[]) <- gamma*(1-IFR[i])*I[i] + w_v*R_v[i] - p_b[i]*V_b[i]*R[i] - p_a[i]*alpha[i]*(1-Se)*R[i] - w_i*R[i]
  
  #age structured mixing
  c[,] <- user() #this is the age-structured contact matrix
  N[] <- S_v[i] + S[i] + E_v[i] + E[i] + I_v[i] + I[i] + R_v[i] + R[i] #idea here is to get total alive in each age group #not sure this is coded correctly
  s_ij[,] <- c[i,j] * (I[i] + I_v[i])/N[i]
  lambda[] <- mu[i]*sum(s_ij[i,])
  lambda_v[] <- (1-Ve_i)*mu[i]*sum(s_ij[i,])
  
  
  ## Initial conditions
  initial(S_v[]) <- S_v_ini[i]
  initial(S[]) <- S_ini[i]
  initial(E_v[]) <- E_v_ini[i]
  initial(E[]) <- E_ini[i]
  initial(I_v[]) <- I_v_ini[i]
  initial(I[]) <- I_ini[i]
  initial(D_v[]) <- D_v_ini[i]
  initial(D[]) <- D_ini[i]
  initial(R_v[]) <- R_v_ini[i]
  initial(R[]) <- R_ini[i]
  
  ## parameters
  ###initial states (to be user defined in each run)
  S_v_ini[] <- user()
  S_ini[] <- user()
  E_v_ini[] <- user()
  E_ini[] <- user()
  I_v_ini[] <- user()
  I_ini[] <- user()
  D_v_ini[] <- user()
  D_ini[] <- user()
  R_v_ini[] <- user()
  R_ini[] <- user()
  
  ###parameter values (to be user defined in each run)

  V_b[] <- user()
  p_a[] <- user()
  p_b[] <- user()
  alpha[] <- user()
  Sp <- user()
  w_vi <- user()
  w_i <- user()
  w_v <- user()
  sigma <- user()
  gamma <- user()
  IFR[] <- user()
  Se <- user()
  mu[] <- user()
  Ve_d <- user()
  Ve_i <- user()
  
  ###define number of age groups and dimensions of compartments and parameters
  N_age <- user() #this is a number of age groups
  dim(S_v) = N_age
  dim(S) = N_age
  dim(E_v) = N_age
  dim(E) = N_age
  dim(I_v) = N_age
  dim(I) = N_age
  dim(R_v) = N_age
  dim(R) = N_age
  dim(D_v) = N_age
  dim(D) = N_age
  
  dim(S_v_ini) = N_age
  dim(S_ini) = N_age
  dim(E_v_ini) = N_age
  dim(E_ini) = N_age
  dim(I_v_ini) = N_age
  dim(I_ini) = N_age
  dim(R_v_ini) = N_age
  dim(R_ini) = N_age
  dim(D_v_ini) = N_age
  dim(D_ini) = N_age
  
  dim(c) <- c(N_age, N_age)
  dim(s_ij) <- c(N_age, N_age)
  
  dim(lambda) <- N_age
  dim(lambda_v) <- N_age
  dim(p_a) <- N_age
  dim(p_b) <- N_age
  dim(alpha) <- N_age
  dim(IFR) <- N_age
  dim(mu) <- N_age
  dim(N) <- N_age
  dim(V_b) <- N_age
  

  
})

cmat <- readRDS("C_USA_bytens_overall.RData") #contact matrix for US in 10 year age bins from Bubar
ages <- readRDS("age_demographics_USA.RData")  #vector of proportion of pop in each age bin from Bubar
N_total <- ages[length(ages)]
avec <- ages[1:9]
#say that the total population is 10^6, figure out how many people will be in each age group in total
age_pops <- 10^6 *avec

#determine the proportion of each age group that is vaccinated
prop_ages_vax <- readRDS('/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/age_structure/prop_ages_vax.RDS')
age_pops_vax <- age_pops*prop_ages_vax
age_pops_novax <- age_pops - age_pops_vax

#determine what fraction of the vaccinated go in infected, exposed, recovered, and susceptible:
i_vax_ages_prop <- readRDS('/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/age_structure/i_vax_ages_prop.RDS')
I_v_ini <- age_pops_vax*i_vax_ages_prop
e_vax_ages_prop <- readRDS('/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/age_structure/e_vax_ages_prop.RDS')
E_v_ini <- age_pops_vax*e_vax_ages_prop
r_vax_ages_prop <- readRDS('/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/age_structure/r_vax_ages_prop.RDS')
R_v_ini <- age_pops_vax*r_vax_ages_prop
S_v_ini <- age_pops_vax - I_v_ini - E_v_ini - R_v_ini

#determine what fraction of the unvaccinated go in infeccted, exposed, recovered and susceptible:
i_unvax_ages_prop <- readRDS('/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/age_structure/i_unvax_ages_prop.RDS')
I_ini <- age_pops_novax * i_unvax_ages_prop
e_unvax_ages_prop <- readRDS('/Users/madeleinekline/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/age_structure/e_unvax_ages_prop.RDS')
E_ini <- age_pops_novax * e_unvax_ages_prop
r_unvax_ages_prop <- readRDS("~/Dropbox (Harvard University)/G1/2023_Spring_Semester_Classes/Epi260/Final_project/seeding_values/age_structure/r_unvax_ages_prop.RDS")
R_ini <- age_pops_novax * r_unvax_ages_prop
S_ini <- age_pops_novax - I_ini - E_ini - R_ini 

#test that all these together equal 10^6
sum(S_ini + S_v_ini + I_ini+ I_v_ini + E_ini + E_v_ini + R_ini + R_v_ini)

#for test have no one antibody testing, no yearly vaccination
test <- age_structured$new(S_v_ini = S_v_ini,
                           S_ini = S_ini,
                           E_v_ini = E_v_ini,
                           E_ini = E_ini,
                           I_v_ini = I_v_ini,
                           I_ini = I_ini,
                           R_v_ini = R_v_ini,
                           R_ini = R_ini,
                           D_v_ini = c(0,0,0,
                                        0,0,0,
                                        0,0,0),
                           D_ini = c(0,0,0,
                                      0,0,0,
                                      0,0,0),
                           N_age = 9,
                           c = cmat,
                           p_a = c(0,0,0,
                                 0,0,0,
                                 0,0,0),
                           p_b = c(0,0,0,
                                   0,0,0,
                                   0,0,0),
                           alpha = c(0,0,0,
                                     0,0,0,
                                     0,0,0),
                           mu = c(0.4, 0.38, 0.79,
                                   0.86, 0.8, 0.82,
                                   0.88, 0.74, 0.74),
                           V_b = c(0,0,0,
                                   0,0,0,
                                   0,0,0),
                           Sp = 0,
                           w_vi = 0.008333333,
                           w_i = 0.008333333,
                           w_v = 0.004166667,
                           sigma = 0.3333333,
                           gamma = 0.2,
                           IFR = c(0.00001, 0.00003, 0.0001,
                                    0.0004, 0.0012, 0.0040,
                                    0.00136, 0.0455, 0.1524),
                           Se = 0,
                           Ve_d = 0.87,
                           Ve_i = 0.3)

t = seq(0,1825)
y <- test$run(t)
test_plot <- as_tibble(y) |>
  pivot_longer(cols = 2:91) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + 
  geom_line() +
  xlab("Time Step") + 
  ggtitle("Test") +
  theme_minimal()

test_deaths <- as_tibble(y) |>
  pivot_longer(cols = 56:73) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + 
  geom_line() +
  xlab("Time Step") + 
  ggtitle("Test Just Deaths") +
  theme_minimal()





#sum the deaths by age group; 
tmp <- as_tibble(y) |>
  mutate(d_under10 = `D_v[1]` + `D[1]`,
         d_10_20 = `D_v[2]` + `D[2]`,
         d_20_30 = `D_v[3]` + `D[3]`,
         d_30_40 = `D_v[4]` + `D[4]`,
         d_40_50 = `D_v[5]` + `D[5]`,
         d_50_60 = `D_v[6]` + `D[6]`,
         d_60_70 = `D_v[7]` + `D[7]`,
         d_70_80 = `D_v[8]` + `D[8]`,
         d_80plus = `D_v[9]` + `D[9]`) |>
  select(t, d_under10, d_10_20, d_20_30, d_30_40,
         d_40_50,  d_50_60, d_60_70, d_70_80, d_80plus) |>
  pivot_longer(cols = 2:10) |>
  ggplot(aes(x = t, y = value, group = name, color = name)) + geom_line() +
  xlab("Time Step") +
  ggtitle("Test Deaths by Age Group") +
  theme_minimal()

total_deaths_over_time <- as_tibble(y) |>
  mutate(total_deaths = `D_v[1]` + `D[1]` + `D_v[2]` + `D[2]` + `D_v[3]` + `D[3]` + 
           `D_v[4]` + `D[4]` + `D_v[5]` + `D[5]` + `D_v[6]` + `D[6]` + `D_v[7]` + `D[7]` + `D_v[8]` + `D[8]` + `D_v[9]` + `D[9]`) |>
  select(t, total_deaths) 

total_deaths_plot <- total_deaths_over_time |>
  ggplot(aes(x = t,  y= total_deaths)) + geom_line() +
  xlab("Time Step") +
  ggtitle("Test Total Deaths Over Time") +
  theme_minimal()

#write a function that produces a bunch of these plots for each scenario 
as_runmodel <- function(title, t, S_v_ini,S_ini,E_v_ini,E_ini,I_v_ini,I_ini,R_v_ini,
                        R_ini,D_v_ini,D_ini,N_age,c,p_a,p_b,alpha,mu,V_b,Sp,w_vi,
                        w_i,w_v,sigma,gamma,IFR,Se,Ve_d,Ve_i){
  mod <- age_structured$new(S_v_ini = S_v_ini,
                            S_ini = S_ini,
                            E_v_ini = E_v_ini,
                            E_ini = E_ini,
                            I_v_ini = I_v_ini,
                            I_ini = I_ini,
                            R_v_ini = R_v_ini,
                            R_ini = R_ini,
                            D_v_ini = D_v_ini,
                            D_ini = D_ini,
                            N_age = N_age,
                            c = c,
                            p_a = p_a,
                            p_b = p_b,
                            alpha = alpha,
                            mu = mu,
                            V_b = V_b,
                            Sp = Sp,
                            w_vi = w_vi,
                            w_i = w_i,
                            w_v = w_v,
                            sigma = sigma,
                            gamma = gamma,
                            IFR = IFR,
                            Se = Se,
                            Ve_d =Ve_d,
                            Ve_i = Ve_i)
  y <- mod$run(t)
  plot1 <- as_tibble(y) |>
    pivot_longer(cols = 2:91) |>
    ggplot(aes(x = t, y = value, group = name, color = name)) + 
    geom_line() +
    xlab("Time Step") + 
    labs(title = "All States Over Time",
         subtitle = title) +
    theme_minimal()
  plot2 <-  as_tibble(y) |>
    pivot_longer(cols = 56:73) |>
    ggplot(aes(x = t, y = value, group = name, color = name)) + 
    geom_line() +
    xlab("Time Step") + 
    labs(title = "Deaths in All Ages Groups Vaxxed vs Unvaxxed",
         subtitle  = title) +
    theme_minimal()
  plot3 <- as_tibble(y) |>
    mutate(d_under10 = `D_v[1]` + `D[1]`,
           d_10_20 = `D_v[2]` + `D[2]`,
           d_20_30 = `D_v[3]` + `D[3]`,
           d_30_40 = `D_v[4]` + `D[4]`,
           d_40_50 = `D_v[5]` + `D[5]`,
           d_50_60 = `D_v[6]` + `D[6]`,
           d_60_70 = `D_v[7]` + `D[7]`,
           d_70_80 = `D_v[8]` + `D[8]`,
           d_80plus = `D_v[9]` + `D[9]`) |>
    select(t, d_under10, d_10_20, d_20_30, d_30_40,
           d_40_50,  d_50_60, d_60_70, d_70_80, d_80plus) |>
    pivot_longer(cols = 2:10) |>
    ggplot(aes(x = t, y = value, group = name, color = name)) + geom_line() +
    xlab("Time Step") +
    ylab("Deaths") +
    scale_color_discrete(name = "Age Group",
                         labels = c("10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+", "<10")) +
    labs(title = "Deaths by Age Group",
         subtitle = title) +
    theme_minimal()
  total_deaths_over_time <- as_tibble(y) |>
    mutate(t = as.numeric(t),
           total_deaths = as.numeric(`D_v[1]` + `D[1]` + `D_v[2]` + `D[2]` + `D_v[3]` + `D[3]` + 
             `D_v[4]` + `D[4]` + `D_v[5]` + `D[5]` + `D_v[6]` + `D[6]` + `D_v[7]` + `D[7]` + `D_v[8]` + `D[8]` + `D_v[9]` + `D[9]`) )|>
    select(t, total_deaths )
  
  return(list(y, plot1, plot2, plot3, total_deaths_over_time))
}


scen0_age <- as_runmodel(title = "No Ab or Yearly Vax", t = seq(0,1825),
                         S_v_ini = S_v_ini,
                         S_ini = S_ini,
                         E_v_ini = E_v_ini,
                         E_ini = E_ini,
                         I_v_ini = I_v_ini,
                         I_ini = I_ini,
                         R_v_ini = R_v_ini,
                         R_ini = R_ini,
            D_v_ini = c(0,0,0,
                        0,0,0,
                        0,0,0),
            D_ini = c(0,0,0,
                      0,0,0,
                      0,0,0),
            N_age = 9,
            c = cmat,
            p_a = c(0,0,0,
                  0,0,0,
                  0,0,0),
            p_b = c(0,0,0,
                    0,0,0,
                    0,0,0),
            alpha = c(0,0,0,
                      0,0,0,
                      0,0,0),
            mu = c(0.4, 0.38, 0.79,
                   0.86, 0.8, 0.82,
                   0.88, 0.74, 0.74),
            V_b = c(0,0,0,
                    0,0,0,
                    0,0,0),
            Sp = 0,
            w_vi = 0.008333333,
            w_i = 0.008333333,
            w_v = 0.004166667,
            sigma = 0.3333333,
            gamma = 0.2,
            IFR = c(0.00001, 0.00003, 0.0001,
                    0.0004, 0.0012, 0.0040,
                    0.00136, 0.0455, 0.1524),
            Se = 0,
            Ve_d = 0.87,
            Ve_i = 0.3 )

scen1_age <- as_runmodel(title = "Only Yearly Antibody", t = seq(0,1825),
                         S_v_ini = S_v_ini,
                         S_ini = S_ini,
                         E_v_ini = E_v_ini,
                         E_ini = E_ini,
                         I_v_ini = I_v_ini,
                         I_ini = I_ini,
                         R_v_ini = R_v_ini,
                         R_ini = R_ini,
                         D_v_ini = c(0,0,0,
                                     0,0,0,
                                     0,0,0),
                         D_ini = c(0,0,0,
                                   0,0,0,
                                   0,0,0),
                         N_age = 9,
                         c = cmat,
                         p_a = c(1,1,1,
                                 1,1,1,
                                 1,1,1),
                         p_b = c(0,0,0,
                                 0,0,0,
                                 0,0,0),
                         alpha = c(0.002739726,0.002739726,0.002739726,
                                   0.002739726,0.002739726,0.002739726,
                                   0.002739726,0.002739726,0.002739726),
                         mu = c(0.4, 0.38, 0.79,
                                0.86, 0.8, 0.82,
                                0.88, 0.74, 0.74),
                         V_b = c(0,0,0,
                                 0,0,0,
                                 0,0,0),
                         Sp = 0,
                         w_vi = 0.008333333,
                         w_i = 0.008333333,
                         w_v = 0.004166667,
                         sigma = 0.3333333,
                         gamma = 0.2,
                         IFR = c(0.00001, 0.00003, 0.0001,
                                 0.0004, 0.0012, 0.0040,
                                 0.00136, 0.0455, 0.1524),
                         Se = 0,
                         Ve_d = 0.87,
                         Ve_i = 0.3 )


scen2_age <- as_runmodel(title = "No Ab, Yearly Vax", t = seq(0,1825), 
                         S_v_ini = S_v_ini,
                         S_ini = S_ini,
                         E_v_ini = E_v_ini,
                         E_ini = E_ini,
                         I_v_ini = I_v_ini,
                         I_ini = I_ini,
                         R_v_ini = R_v_ini,
                         R_ini = R_ini,
                         D_v_ini = c(0,0,0,
                                     0,0,0,
                                     0,0,0),
                         D_ini = c(0,0,0,
                                   0,0,0,
                                   0,0,0),
                         N_age = 9,
                         c = cmat,
                         p_a = c(0,0,0,
                               0,0,0,
                               0,0,0),
                         p_b = c(1,1,1,
                                 1,1,1,
                                 1,1,1),
                         alpha = c(0,0,0,
                                   0,0,0,
                                   0,0,0),
                         mu = c(0.4, 0.38, 0.79,
                                0.86, 0.8, 0.82,
                                0.88, 0.74, 0.74),
                         V_b = c(0.002739726,0.002739726, 0.002739726,
                                 0.002739726,0.002739726,0.002739726,
                                 0.002739726,0.002739726,0.002739726),
                         Sp = 0,
                         w_vi = 0.008333333,
                         w_i = 0.008333333,
                         w_v = 0.004166667,
                         sigma = 0.3333333,
                         gamma = 0.2,
                         IFR = c(0.00001, 0.00003, 0.0001,
                                 0.0004, 0.0012, 0.0040,
                                 0.00136, 0.0455, 0.1524),
                         Se = 0,
                         Ve_d = 0.87,
                         Ve_i = 0.3 )

scen3_age <- as_runmodel(title = "Monthly Ab testing in all, No Yearly Vax", t = seq(0,1825),
                         S_v_ini = S_v_ini,
                         S_ini = S_ini,
                         E_v_ini = E_v_ini,
                         E_ini = E_ini,
                         I_v_ini = I_v_ini,
                         I_ini = I_ini,
                         R_v_ini = R_v_ini,
                         R_ini = R_ini,
                                      D_v_ini = c(0,0,0,
                                                  0,0,0,
                                                  0,0,0),
                                      D_ini = c(0,0,0,
                                                0,0,0,
                                                0,0,0),
                                      N_age = 9,
                                      c = cmat,
                                      p_a = c(1,1,1,
                                            1,1,1,
                                            1,1,1),
                                      p_b = c(0,0,0,
                                              0,0,0,
                                              0,0,0),
                                      alpha = c(0.3333333,0.3333333,0.3333333,
                                                0.3333333,0.3333333,0.3333333,
                                                0.3333333,0.3333333,0.3333333),
                                      mu = c(0.4, 0.38, 0.79,
                                             0.86, 0.8, 0.82,
                                             0.88, 0.74, 0.74),
                                      V_b = c(0,0,0,
                                              0,0,0,
                                              0,0,0),
                                      Sp = 0,
                                      w_vi = 0.008333333,
                                      w_i = 0.008333333,
                                      w_v = 0.004166667,
                                      sigma = 0.3333333,
                                      gamma = 0.2,
                                      IFR = c(0.00001, 0.00003, 0.0001,
                                              0.0004, 0.0012, 0.0040,
                                              0.00136, 0.0455, 0.1524),
                                      Se = 0,
                                      Ve_d = 0.87,
                                      Ve_i = 0.3 )
scen4_age <- as_runmodel(title = "Monthly Ab testing in all,Yearly Vax", t = seq(0,1825),S_v_ini = S_v_ini,
                         S_ini = S_ini,
                         E_v_ini = E_v_ini,
                         E_ini = E_ini,
                         I_v_ini = I_v_ini,
                         I_ini = I_ini,
                         R_v_ini = R_v_ini,
                         R_ini = R_ini,
                         D_v_ini = c(0,0,0,
                                     0,0,0,
                                     0,0,0),
                         D_ini = c(0,0,0,
                                   0,0,0,
                                   0,0,0),
                         N_age = 9,
                         c = cmat,
                         p_a = c(1,1,1,
                               1,1,1,
                               1,1,1),
                         p_b = c(1,1,1,
                                 1,1,1,
                                 1,1,1),
                         alpha = c(0.3333333,0.3333333,0.3333333,
                                   0.3333333,0.3333333,0.3333333,
                                   0.3333333,0.3333333,0.3333333),
                         mu = c(0.4, 0.38, 0.79,
                                0.86, 0.8, 0.82,
                                0.88, 0.74, 0.74),
                         V_b = c(0.002739726, 0.002739726,0.002739726,
                                 0.002739726,0.002739726,0.002739726,
                                 0.002739726,0.002739726,0.002739726),
                         Sp = 0,
                         w_vi = 0.008333333,
                         w_i = 0.008333333,
                         w_v = 0.004166667,
                         sigma = 0.3333333,
                         gamma = 0.2,
                         IFR = c(0.00001, 0.00003, 0.0001,
                                 0.0004, 0.0012, 0.0040,
                                 0.00136, 0.0455, 0.1524),
                         Se = 0,
                         Ve_d = 0.87,
                         Ve_i = 0.3 )
colnames(scen0_age[[5]]) = c("t", "scen0")
colnames(scen1_age[[5]]) = c("t", "scen1")
colnames(scen2_age[[5]]) = c("t", "scen2")
colnames(scen3_age[[5]]) = c("t", "scen3")
colnames(scen4_age[[5]]) = c("t", "scen4")
alldeaths_df <- left_join(scen0_age[[5]], scen1_age[[5]]) |> left_join(scen2_age[[5]]) |>
  left_join(scen3_age[[5]]) |> left_join(scen4_age[[5]])
age_deaths_fourscens <- alldeaths_df |>
  pivot_longer(cols = 2:6) |>
  ggplot(aes(x = t, y = value, color = name)) + 
  geom_line() + 
  scale_color_discrete(name = "Scenarios",
                       labels = c("No Antibody Testing, No Yearly Vaccination",
                                  "Yearly Antibody Testing in Everyone",
                                  "No Antibody Testing, Yearly Vaccination in Everyone",
                                  "Monthly Antibody Testing in Everyone, No Yearly Vaccination",
                                  "Monthly Antibody Testing and Yearly Vaccination in Everyone")) +
  theme_minimal() +
  ggtitle("Deaths Over Time")

#ggsave("age_deaths_fourscens.png", age_deaths_fourscens)

deaths_by_age_groups_scen <- grid.arrange(scen0_age[[4]], scen1_age[[4]], scen2_age[[4]], scen3_age[[4]] + ylim(0,30000), scen4_age[[4]] + ylim(0,30000))
#ggsave("age_deaths_by_age_group.png", deaths_by_age_groups_scen)



#now replicate the same 4 scenarios as with no age structure; 25%, 50% and 75% participation in all ages
#with additional yearly vaccination in all
#scen2age is comparison, with only yearly vax
#assuming everyone gets vaccinated yearly, but not everyone participates in monthly antibody testing
scen2_part <- as_runmodel(title = "Monthly Ab testing in all, Yearly Vax \n 90% Sensitivty/Specificity, 25% Participation", t = seq(0,1825), S_v_ini = S_v_ini,
                          S_ini = S_ini,
                          E_v_ini = E_v_ini,
                          E_ini = E_ini,
                          I_v_ini = I_v_ini,
                          I_ini = I_ini,
                          R_v_ini = R_v_ini,
                          R_ini = R_ini,
                          D_v_ini = c(0,0,0,
                                      0,0,0,
                                      0,0,0),
                          D_ini = c(0,0,0,
                                    0,0,0,
                                    0,0,0),
                          N_age = 9,
                          c = cmat,
                          p_a = c(0.25,0.25,0.25,
                                0.25,0.25,0.25,
                                0.25,0.25,0.25),
                          p_b = c(1,1,1,1,1,1,1,1,1),
                          alpha = c(0.3333333,0.3333333,0.3333333,
                                    0.3333333,0.3333333,0.3333333,
                                    0.3333333,0.3333333,0.3333333),
                          mu = c(0.4, 0.38, 0.79,
                                 0.86, 0.8, 0.82,
                                 0.88, 0.74, 0.74),
                          V_b = c(0.002739726, 0.002739726,0.002739726,
                                  0.002739726,0.002739726,0.002739726,
                                  0.002739726,0.002739726,0.002739726),
                          Sp = 0,
                          w_vi = 0.008333333,
                          w_i = 0.008333333,
                          w_v = 0.004166667,
                          sigma = 0.3333333,
                          gamma = 0.2,
                          IFR = c(0.00001, 0.00003, 0.0001,
                                  0.0004, 0.0012, 0.0040,
                                  0.00136, 0.0455, 0.1524),
                          Se = 0,
                          Ve_d = 0.87,
                          Ve_i = 0.3 )

scen3_part <- as_runmodel(title = "Monthly Ab testing in all, Yearly Vax \n 90% Sensitivty/Specificity, 50% Participation", t = seq(0,1825), S_v_ini = S_v_ini,
                          S_ini = S_ini,
                          E_v_ini = E_v_ini,
                          E_ini = E_ini,
                          I_v_ini = I_v_ini,
                          I_ini = I_ini,
                          R_v_ini = R_v_ini,
                          R_ini = R_ini,
                          D_v_ini = c(0,0,0,
                                      0,0,0,
                                      0,0,0),
                          D_ini = c(0,0,0,
                                    0,0,0,
                                    0,0,0),
                          N_age = 9,
                          c = cmat,
                          p_a = c(0.5,0.5,0.5,
                                0.5,0.5,0.5,
                                0.5,0.5,0.5),
                          alpha = c(0.3333333,0.3333333,0.3333333,
                                    0.3333333,0.3333333,0.3333333,
                                    0.3333333,0.3333333,0.3333333),
                          mu = c(0.4, 0.38, 0.79,
                                 0.86, 0.8, 0.82,
                                 0.88, 0.74, 0.74),
                          p_b = c(1,1,1,1,1,1,1,1,1),
                          V_b = c(0.002739726, 0.002739726,0.002739726,
                                  0.002739726,0.002739726,0.002739726,
                                  0.002739726,0.002739726,0.002739726),
                          Sp = 0,
                          w_vi = 0.008333333,
                          w_i = 0.008333333,
                          w_v = 0.004166667,
                          sigma = 0.3333333,
                          gamma = 0.2,
                          IFR = c(0.00001, 0.00003, 0.0001,
                                  0.0004, 0.0012, 0.0040,
                                  0.00136, 0.0455, 0.1524),
                          Se = 0,
                          Ve_d = 0.87,
                          Ve_i = 0.3 )

scen4_part <- as_runmodel(title = "Monthly Ab testing in all, Yearly Vax \n 90% Sensitivty/Specificity, 75% Participation", t = seq(0,1825), S_v_ini = S_v_ini,
                          S_ini = S_ini,
                          E_v_ini = E_v_ini,
                          E_ini = E_ini,
                          I_v_ini = I_v_ini,
                          I_ini = I_ini,
                          R_v_ini = R_v_ini,
                          R_ini = R_ini,
                          D_v_ini = c(0,0,0,
                                      0,0,0,
                                      0,0,0),
                          D_ini = c(0,0,0,
                                    0,0,0,
                                    0,0,0),
                          N_age = 9,
                          c = cmat,
                          p_a = c(0.75,0.75,0.75,
                                0.75,0.75,0.75,
                                0.75,0.75,0.75),
                          p_b = c(1,1,1,1,1,1,1,1,1),
                          alpha = c(0.3333333,0.3333333,0.3333333,
                                    0.3333333,0.3333333,0.3333333,
                                    0.3333333,0.3333333,0.3333333),
                          mu = c(0.4, 0.38, 0.79,
                                 0.86, 0.8, 0.82,
                                 0.88, 0.74, 0.74),
                          V_b = c(0.002739726, 0.002739726,0.002739726,
                                  0.002739726,0.002739726,0.002739726,
                                  0.002739726,0.002739726,0.002739726),
                          Sp = 0,
                          w_vi = 0.008333333,
                          w_i = 0.008333333,
                          w_v = 0.004166667,
                          sigma = 0.3333333,
                          gamma = 0.2,
                          IFR = c(0.00001, 0.00003, 0.0001,
                                  0.0004, 0.0012, 0.0040,
                                  0.00136, 0.0455, 0.1524),
                          Se = 0,
                          Ve_d = 0.87,
                          Ve_i = 0.3 )

participation_age <- grid.arrange(scen2_age[[4]], scen2_part[[4]] + ylim(0,30000), scen3_part[[4]]+ ylim(0,30000), scen4_part[[4]]+ ylim(0,30000), scen4_age[[4]]+ ylim(0,30000))
#ggsave("participation_scenarios_age.png", participation_age)


#try one last set of scenarios:
#different strategies by age group. 
#scen1: yearly vaccines in everyone, no ab testing -> already done, this is scen2_age
#scen2: yearly vaccines only in those older than 70, no vaccines in anybody else
#scen3: no yearly vaccines, antibody testing in everyone -> already done, this is scen3_age
#scen4: no yearly vaccines, antibody testing only in those older than 70
#scen5: yearly vaccines in everyone, antibody testing only in those older 70
#scen6: yearly vaccines and antibody testing in everyone -> already done, scen4_age


scen2_agediff <- as_runmodel(title = "Yearly Vax >70, no other additional Vax", t = seq(0,1825), S_v_ini = S_v_ini,
                             S_ini = S_ini,
                             E_v_ini = E_v_ini,
                             E_ini = E_ini,
                             I_v_ini = I_v_ini,
                             I_ini = I_ini,
                             R_v_ini = R_v_ini,
                             R_ini = R_ini,
                             D_v_ini = c(0,0,0,
                                         0,0,0,
                                         0,0,0),
                             D_ini = c(0,0,0,
                                       0,0,0,
                                       0,0,0),
                             N_age = 9,
                             c = cmat,
                             p_a = c(0,0,0,
                                   0,0,0,
                                   0,0,0),
                             alpha = c(0,0,0,
                                       0,0,0,
                                       0,0,0),
                             mu = c(0.4, 0.38, 0.79,
                                    0.86, 0.8, 0.82,
                                    0.88, 0.74, 0.74),
                             V_b = c(0,0,0,
                                     0,0,0,
                                     0,0.002739726,0.002739726),
                             p_b = c(0,0,0,
                                     0,0,0,
                                     0,1,1),
                             Sp = 0,
                             w_vi = 0.008333333,
                             w_i = 0.008333333,
                             w_v = 0.004166667,
                             sigma = 0.3333333,
                             gamma = 0.2,
                             IFR = c(0.00001, 0.00003, 0.0001,
                                     0.0004, 0.0012, 0.0040,
                                     0.00136, 0.0455, 0.1524),
                             Se = 0,
                             Ve_d = 0.87,
                             Ve_i = 0.3 )


scen4_agediff <- as_runmodel(title = "No Yearly Vax, Monthly Ab testing in >70", t = seq(0,1825), S_v_ini = S_v_ini,
                             S_ini = S_ini,
                             E_v_ini = E_v_ini,
                             E_ini = E_ini,
                             I_v_ini = I_v_ini,
                             I_ini = I_ini,
                             R_v_ini = R_v_ini,
                             R_ini = R_ini,
                             D_v_ini = c(0,0,0,
                                         0,0,0,
                                         0,0,0),
                             D_ini = c(0,0,0,
                                       0,0,0,
                                       0,0,0),
                             N_age = 9,
                             c = cmat,
                             p_a = c(0,0,0,
                                   0,0,0,
                                   0,1,1),
                             p_b = c(0,0,0,
                                     0,0,0,
                                     0,0,0),
                             alpha = c(0,0,0,
                                       0,0,0,
                                       0,0.3333333,0.3333333),
                             mu = c(0.4, 0.38, 0.79,
                                    0.86, 0.8, 0.82,
                                    0.88, 0.74, 0.74),
                             V_b = c(0, 0,0,
                                     0,0,0,
                                     0,0,0),
                             Sp = 0,
                             w_vi = 0.008333333,
                             w_i = 0.008333333,
                             w_v = 0.004166667,
                             sigma = 0.3333333,
                             gamma = 0.2,
                             IFR = c(0.00001, 0.00003, 0.0001,
                                     0.0004, 0.0012, 0.0040,
                                     0.00136, 0.0455, 0.1524),
                             Se = 0,
                             Ve_d = 0.87,
                             Ve_i = 0.3 )



#scen5: yearly vaccines in everyone, antibody testing only in those older 70
scen5_agediff <- as_runmodel(title = "Yearly Vax in all, Monthly Ab testing in >70", t = seq(0,1825), S_v_ini = S_v_ini,
                             S_ini = S_ini,
                             E_v_ini = E_v_ini,
                             E_ini = E_ini,
                             I_v_ini = I_v_ini,
                             I_ini = I_ini,
                             R_v_ini = R_v_ini,
                             R_ini = R_ini,
                             D_v_ini = c(0,0,0,
                                         0,0,0,
                                         0,0,0),
                             D_ini = c(0,0,0,
                                       0,0,0,
                                       0,0,0),
                             N_age = 9,
                             c = cmat,
                             p_a = c(0,0,0,
                                   0,0,0,
                                   0,1,1),
                             p_b = c(1,1,1,1,1,1,1,1,1),
                             alpha = c(0,0,0,
                                       0,0,0,
                                       0,0.3333333,0.3333333),
                             mu = c(0.4, 0.38, 0.79,
                                    0.86, 0.8, 0.82,
                                    0.88, 0.74, 0.74),
                             V_b = c(0.002739726, 0.002739726,0.002739726,
                                     0.002739726,0.002739726,0.002739726,
                                     0.002739726,0.002739726,0.002739726),
                             Sp = 0,
                             w_vi = 0.008333333,
                             w_i = 0.008333333,
                             w_v = 0.004166667,
                             sigma = 0.3333333,
                             gamma = 0.2,
                             IFR = c(0.00001, 0.00003, 0.0001,
                                     0.0004, 0.0012, 0.0040,
                                     0.00136, 0.0455, 0.1524),
                             Se = 0,
                             Ve_d = 0.87,
                             Ve_i = 0.3 )




#scen1: yearly vaccines in everyone, no ab testing -> already done, this is scen2_age
#scen2: yearly vaccines only in those older than 70, no vaccines in anybody else
#scen3: no yearly vaccines, antibody testing in everyone -> already done, this is scen3_age
#scen4: no yearly vaccines, antibody testing only in those older than 70
#scen5: yearly vaccines in everyone, antibody testing only in those older 70
#scen6: yearly vaccines and antibody testing in everyone -> already done, scen4_age
diff_age_strats <- grid.arrange(scen2_age[[4]], scen2_agediff[[4]], scen3_age[[4]] + ylim(0,30000), 
             scen4_agediff[[4]]+ ylim(0,30000), scen5_agediff[[4]]+ ylim(0,30000), scen4_age[[4]]+ ylim(0,30000))
#ggsave('diff_age_strats.png', diff_age_strats)

colnames(scen2_age[[5]]) = c("t", "scen1_agediff")
colnames(scen2_agediff[[5]]) = c("t", "scen2_agediff")
colnames(scen3_age[[5]]) = c("t", "scen3_agediff")
colnames(scen4_agediff[[5]]) = c("t", "scen4_agediff")
colnames(scen5_agediff[[5]]) = c("t", "scen5_agediff")
colnames(scen4_age[[5]]) = c("t", "scen6_diff")

alldeaths_df_agediff <- left_join(scen2_age[[5]], scen2_agediff[[5]]) |>
  left_join(scen3_age[[5]]) |>
  left_join(scen4_agediff[[5]]) |>
  left_join(scen5_agediff[[5]]) |>
  left_join(scen4_age[[5]])
  
  
agediff_deaths_scens <- alldeaths_df_agediff |>
  pivot_longer(cols = 2:7) |>
  ggplot(aes(x = t, y = value, color = name)) + 
  geom_line() + 
  scale_color_discrete(name = "Scenario",
                       labels = c("Yearly Vaccines in Everyone, No Antibody Testing",
                                  "Yearly Vaccines in those >70, No Vaccines Otherwise",
                                  "No Yearly Vaccines, Monthly Antibody Testing in Everyone",
                                  "No Yearly Vaccines, Monthly Antibody Testing in those >70 Only",
                                  "Yearly Vaccines in Everyone, Monthly Antibody Testing in those >70 Only",
                                  "Yearly Vaccines in Everyone, Monthly Antibody Testing in Everyone")) +
  theme_minimal() +
  ggtitle("Deaths Over Time")

#ggsave("agediff_deaths_scens.png",agediff_deaths_scens )

make_final_deaths_plot <- function(list){#takes a list of model run list object 5s as an output
  for(i in 1:length(list)){
    colnames(list[[i]]) = c("t", paste("scenario", i))
  }
  
  alldeath_df <- left_join(list[[1]], list[[2]])
  for(i in 3:length(list)){
    alldeath_df <- left_join(alldeath_df, list[[i]])
  }
  
  deaths_over_time <- alldeath_df |>
    pivot_longer(cols = 2:ncol(alldeath_df)) |>
    ggplot(aes(x = t, y = value, color = name)) + 
    geom_line() + 
    theme_minimal() +
    ggtitle("Deaths Over Time")
  final_deaths <- alldeath_df[nrow(alldeath_df),] |>
    pivot_longer(cols = 2:ncol(alldeath_df)) |>
    ggplot(aes(x = reorder(name, -value), y = value, fill = name)) +
    geom_bar(stat = "identity") + 
    theme_minimal() +
    theme(axis.text.x=element_blank()) +
    xlab("Scenario") +
    ylab("Deaths")
  
  return(list(deaths_over_time, final_deaths))
}

#first set of scenarios
deaths_scenarios1 <- make_final_deaths_plot(list(scen0_age[[5]], scen1_age[[5]],scen2_age[[5]], scen3_age[[5]], scen4_age[[5]]))[[2]] + 
  labs(title = "Total Deaths after 5 Years", subtitle = "Different Yearly Vaccination vs Monthly Antibody Strategies in All Ages") +
  scale_fill_discrete(name = "Scenario",
                       labels = c("No Antibody Testing, No Yearly Vaccination",
                                  "Yearly Antibody Testing in Everyone",
                                  "No Antibody Testing, Yearly Vaccination in Everyone",
                                  "Monthly Antibody Testing in Everyone, No Yearly Vaccination",
                                  "Monthly Antibody Testing and Yearly Vaccination in Everyone"))
#ggsave("ages_scenarioset1_finaldeaths.png", deaths_scenarios1)

#second set of scenarios (participation)
deaths_scenarios2 <- make_final_deaths_plot(list(scen2_age[[5]], scen2_part[[5]], scen3_part[[5]], scen4_part[[5]], scen4_age[[5]]))
scens2_dot <- deaths_scenarios2[[1]] + labs(title = "Deaths Over Time", subtitle = "Differing Antibody Testing Participation") + ylab("Deaths") + 
  scale_color_discrete(name = "Scenario", labels = c("Yearly Vaccination in Everyone, No Antibody Testing",
                               "Yearly Vaccination in Everyone, Monthly Antibody Testing in 25%",
                               "Yearly Vaccination in Everyone, Monthly Antibody Testing in 50%",
                               "Yearly Vaccination in Everyone, Monthly Antibody testing in 75%",
                               "Yaerly Vaccination in Everyone, Monthly Antibody Testing in Everyone"))
scens2_fd <- deaths_scenarios2[[2]] + labs(title = "Total Deaths after 5 Years", subtitle = "Differing Antibody Testing Participation") +  
  scale_fill_discrete(name = "Scenario", labels = c("Yearly Vaccination in Everyone, No Antibody Testing",
                                                     "Yearly Vaccination in Everyone, Monthly Antibody Testing in 25%",
                                                     "Yearly Vaccination in Everyone, Monthly Antibody Testing in 50%",
                                                     "Yearly Vaccination in Everyone, Monthly Antibody testing in 75%",
                                                     "Yaerly Vaccination in Everyone, Monthly Antibody Testing in Everyone"))
scens2_deaths <- grid.arrange(scens2_dot, scens2_fd, ncol = 2)
#ggsave("ages_scenarioset2_deaths.png", scens2_deaths, width =17, height = 10)


diff_age_strats_deaths <- make_final_deaths_plot(list(scen1_age[[5]], scen2_agediff[[5]], scen3_age[[5]], 
                                scen4_agediff[[5]], scen5_agediff[[5]], scen4_age[[5]]))

scens3_dot <- diff_age_strats_deaths[[1]] + scale_color_discrete(name = "Scenario",
                                                   labels = c("Yearly Vaccines in Everyone, No Antibody Testing",
                                                              "Yearly Vaccines in those >70, No Vaccines Otherwise",
                                                              "No Yearly Vaccines, Monthly Antibody Testing in Everyone",
                                                              "No Yearly Vaccines, Monthly Antibody Testing in those >70 Only",
                                                              "Yearly Vaccines in Everyone, Monthly Antibody Testing in those >70 Only",
                                                              "Yearly Vaccines in Everyone, Monthly Antibody Testing in Everyone")) +
  ylab("Deaths")
scens3_fd <- diff_age_strats_deaths[[2]] + scale_fill_discrete(name = "Scenario",
                                                     labels = c("Yearly Vaccines in Everyone, No Antibody Testing",
                                                                "Yearly Vaccines in those >70, No Vaccines Otherwise",
                                                                "No Yearly Vaccines, Monthly Antibody Testing in Everyone",
                                                                "No Yearly Vaccines, Monthly Antibody Testing in those >70 Only",
                                                                "Yearly Vaccines in Everyone, Monthly Antibody Testing in those >70 Only",
                                                                "Yearly Vaccines in Everyone, Monthly Antibody Testing in Everyone")) + 
  labs(title = "Total Deaths After 5 Years", subtitle = "Different Age-Specific Vaccination Strategies")
scens3_deaths <- grid.arrange(scens3_dot, scens3_fd, ncol = 2)

#ggsave("diff_age_strats_deaths.png", scens3_deaths, width = 20)


