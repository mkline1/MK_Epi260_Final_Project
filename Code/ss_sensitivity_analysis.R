#goal: determine for what regions of sensitivity / specificity combinations having everyone antibody test monthly outperforms yearly vaccination
library(tidyverse)
library(gridExtra)
source("final_model_script_age_structure.R")


#write a function that runs the model, calculates the final number of deaths, and compares it to final deaths in the comparator scenario
#comparator scenario is yearly vaccination. Run this and calculate total deaths:
comp <- as_runmodel(title = "No Ab, Yearly Vax", t = seq(0,1825), 
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
#extract final deaths
comp_deaths <- comp[[5]]$total_deaths[nrow(comp[[5]])]


compare_deaths_rat <- function(sensitivity, specificity){
  mod <- as_runmodel(title = "Monthly Ab testing in all, No Yearly Vax", t = seq(0,1825),
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
                     Sp = specificity,
                     w_vi = 0.008333333,
                     w_i = 0.008333333,
                     w_v = 0.004166667,
                     sigma = 0.3333333,
                     gamma = 0.2,
                     IFR = c(0.00001, 0.00003, 0.0001,
                             0.0004, 0.0012, 0.0040,
                             0.00136, 0.0455, 0.1524),
                     Se = sensitivity,
                     Ve_d = 0.87,
                     Ve_i = 0.3 )
  mod_final_deaths <- mod[[5]]$total_deaths[nrow(mod[[5]])]
  return(mod_final_deaths / comp_deaths)
}

sensitivity_range <- seq(0,1,0.1)
specificity_range <- seq(0,1,0.1)
df <- data.frame(expand.grid("Sensitivity" = sensitivity_range, "Specificity" = specificity_range)) |> mutate("Ratio" = NA)
for(i in 1:nrow(df)){
  df$Ratio[i] = compare_deaths_rat(df$Sensitivity[i], df$Specificity[i])
}

#pivot wider this dataframe
abss_heatmap <- df |> ggplot(aes(Sensitivity, Specificity, fill = Ratio)) + geom_tile() +
  scale_fill_gradient2(name = "Ratio of Deaths \n Compared to Yearly Vaccination", low = "blue", mid ="white", high = "red",  midpoint = 1) +theme_minimal() +
  labs(title = "Sensitivity and Specificty of Home Antibody Test",
       subtitle = "Antibody test taken by everyone in all ages monthly")
#ggsave("AbSS_heatmap.png",abss_heatmap)
sens1_plot <- df |> filter(Sensitivity == 1) |> mutate(col = ifelse(Ratio <= 1, "blue", "red")) |>
  ggplot(aes(Specificity, Ratio)) + geom_line(color = "blue") + geom_line(aes(color = col)) +
  scale_color_manual(values = c("blue", "red"), 
                     name = "Color",
                     labels = c("Better than yearly vaccination",
                                "Worse than yearly vaccination")) + theme_minimal() +
  geom_hline(yintercept = 1, color = "grey", linetype = 3) +
  geom_vline(xintercept = 0.1, color = "grey", linetype = 2) +
  ylab("Ratio of deaths compared to yearly vaccination") +
  labs(title = "Home Antibody Test Deaths with Varying Specificity",
       subtitle = "Sensitivity of 1")
ab_ss_plots <- grid.arrange(abss_heatmap, sens1_plot, ncol = 2)
#ggsave("ab_ss_plots.png", ab_ss_plots, width =15)


#undertake similar analyses looking at: 1) participation in baseline vaccination and antibody testing
compare_deaths_rat_part <- function(p_a, p_b){
  mod <- as_runmodel(title = "Monthly Ab testing in all, No Yearly Vax", t = seq(0,1825),
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
                     p_a = c(rep(p_a,9)),
                     p_b = c(rep(p_b,9)),
                     alpha = c(0.3333333,0.3333333,0.3333333,
                               0.3333333,0.3333333,0.3333333,
                               0.3333333,0.3333333,0.3333333),
                     mu = c(0.4, 0.38, 0.79,
                            0.86, 0.8, 0.82,
                            0.88, 0.74, 0.74),
                     V_b = c(0.002739726,0.002739726,0.002739726,
                             0.002739726,0.002739726,0.002739726,
                             0.002739726,0.002739726,0.002739726),
                     Sp = 0.9,
                     w_vi = 0.008333333,
                     w_i = 0.008333333,
                     w_v = 0.004166667,
                     sigma = 0.3333333,
                     gamma = 0.2,
                     IFR = c(0.00001, 0.00003, 0.0001,
                             0.0004, 0.0012, 0.0040,
                             0.00136, 0.0455, 0.1524),
                     Se = 0.9,
                     Ve_d = 0.87,
                     Ve_i = 0.3 )
  mod_final_deaths <- mod[[5]]$total_deaths[nrow(mod[[5]])]
  return(mod_final_deaths / comp_deaths)
}
  
part_range <- seq(0,1,0.1)
part_df <- data.frame(expand.grid("P_a" = part_range, "P_b" = part_range)) |> mutate("Ratio" = NA)
for(i in 1:nrow(part_df)){
  part_df$Ratio[i] = compare_deaths_rat_part(part_df$P_a[i], part_df$P_b[i])
}
#pivot wider this dataframe to make heat map
part_heatmap <- part_df |> ggplot(aes(P_a, P_b, fill = Ratio)) + geom_tile() +
  scale_fill_gradient2(name = "Ratio of Deaths \n Compared to Yearly Vaccination 100% Participation", low = "blue", mid ="white", high = "red",  midpoint = 1) +theme_minimal() +
  labs(title = "Participation Rate in Yearly Vaccination and Monthly Home Antibody Tests",
       subtitle = "Antibody Test Sensitivity/Specificity = 0.9")
#ggsave("Participation_heatmap.png", part_heatmap)


#need to redo this with higher resolution 
compare_deaths_rat_part_p_a<- function(p_a){
  mod <- as_runmodel(title = "Monthly Ab testing in all, No Yearly Vax", t = seq(0,1825),
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
                     p_a = c(rep(p_a,9)),
                     p_b = c(rep(0,9)),
                     alpha = c(0.3333333,0.3333333,0.3333333,
                               0.3333333,0.3333333,0.3333333,
                               0.3333333,0.3333333,0.3333333),
                     mu = c(0.4, 0.38, 0.79,
                            0.86, 0.8, 0.82,
                            0.88, 0.74, 0.74),
                     V_b = c(0.002739726,0.002739726,0.002739726,
                             0.002739726,0.002739726,0.002739726,
                             0.002739726,0.002739726,0.002739726),
                     Sp = 0.9,
                     w_vi = 0.008333333,
                     w_i = 0.008333333,
                     w_v = 0.004166667,
                     sigma = 0.3333333,
                     gamma = 0.2,
                     IFR = c(0.00001, 0.00003, 0.0001,
                             0.0004, 0.0012, 0.0040,
                             0.00136, 0.0455, 0.1524),
                     Se = 0.9,
                     Ve_d = 0.87,
                     Ve_i = 0.3 )
  mod_final_deaths <- mod[[5]]$total_deaths[nrow(mod[[5]])]
  return(mod_final_deaths / comp_deaths)
}

p_a_range <- seq(0,1,0.01)
p_a_rats <- unlist(lapply(p_a_range, compare_deaths_rat_part_p_a))
part_plot <-data.frame(p_a = p_a_range, ratio = p_a_rats) |> mutate(col = ifelse(ratio <= 1, "blue", "red")) |>
  ggplot(aes(p_a, ratio)) + geom_line(aes(color = "red")) + geom_line(aes(color = col)) + scale_color_manual(values = c("blue", "red"), 
                                                                                                             name = "Color",
                                                                                                             labels = c("Better than yearly vaccination",
                                                                                                                        "Worse than yearly vaccination")) +
  theme_minimal() + geom_hline(yintercept = 1, linetype = 3, color = "grey") +
  geom_vline(xintercept = 0.05, linetype = 2, color = "grey") +
  labs(title = "Monthly Antibody Testing Participation Resulting Deaths Compared to Yearly Vaccination in All",
       subtitle = "P_b = 0, Sensitivity/Specificity = 0.9")
part_plots <- grid.arrange(part_heatmap, part_plot, ncol = 2)
#ggsave("participation_plots.png", part_plots, width = 20)

#now look at antibody testing frequency over time
compare_deaths_rat_freq <- function(days){
  mod <- as_runmodel(title = "Monthly Ab testing in all, No Yearly Vax", t = seq(0,1825),
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
                     alpha = c(rep(1/days,9)),
                     mu = c(0.4, 0.38, 0.79,
                            0.86, 0.8, 0.82,
                            0.88, 0.74, 0.74),
                     V_b = c(0,0,0,
                             0,0,0,
                             0,0,0),
                     Sp = 0.9,
                     w_vi = 0.008333333,
                     w_i = 0.008333333,
                     w_v = 0.004166667,
                     sigma = 0.3333333,
                     gamma = 0.2,
                     IFR = c(0.00001, 0.00003, 0.0001,
                             0.0004, 0.0012, 0.0040,
                             0.00136, 0.0455, 0.1524),
                     Se = 0.9,
                     Ve_d = 0.87,
                     Ve_i = 0.3 )
  mod_final_deaths <- mod[[5]]$total_deaths[nrow(mod[[5]])]
  return(mod_final_deaths / comp_deaths)
}

freq_opts <- c(1:365) #this ranges from daily to yearly
frequency_deaths <- lapply(freq_opts, compare_deaths_rat_freq)
head(unlist(frequency_deaths))

frequency_plot <- data.frame("freq" = freq_opts, "deaths" = unlist(frequency_deaths)) |>
  mutate(col = ifelse(deaths <= 1, 'blue', "red")) |>
  ggplot(aes(freq_opts, deaths, color = col)) + geom_line() + scale_color_manual(values = c("blue", "red"), 
                                                                                 name = "Color",
                                                                                 labels = c("Better than yearly vaccination",
                                                                                            "Worse than yearly vaccination")) +
  xlab("Testing Frequency in Days") +
  geom_hline(yintercept = 1, linetype = 3, color = "grey") +
  geom_vline(xintercept = 62, linetype =2, color = "black") + ylab("Ratio of deaths compared to yearly vaccination") +
  theme_minimal() +
  labs(title = "Home Antibody Testing Frequency")

#ggsave("testing_frequency_analysis.png", frequency_plot)

#ggsave("ab_part_plot.png", grid.arrange(abss_heatmap, sens1_plot, part_heatmap, part_plot, ncol = 2), width =13)




