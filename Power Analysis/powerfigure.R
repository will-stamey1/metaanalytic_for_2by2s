library(tidyr)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())

a <- readxl::read_xlsx("Power Analysis/power_results_8-14.xlsx")

# a %>% rename("Assurance" = "assurrance",
#              "Power" = "power",
#              "Power 2" = "power2") %>% 
a %>%  
  pivot_longer(cols = c("power", "power2", "assurrance"), values_to = "rate") %>% 
  ggplot(aes(x = n, y = rate)) + 
  geom_point(aes(color = name)) +
  geom_smooth(method = "gam", formula = y ~ s(x), mapping = aes(color = name), se = FALSE) +
  scale_color_manual(values = c("pink", "lightblue", "forestgreen"), 
                     labels = c("Assurance", "Power", "Power w/Nuisance Prior")) + 
  labs(color = "",
       x = "Sample Size",
       y = "Rejection Rate")
