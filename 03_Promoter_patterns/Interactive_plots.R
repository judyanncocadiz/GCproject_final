#interactive plot
library(ggplot2)
library(plotly)

AP <- ggplot(prmtr.means, aes(position, ratio, color = sample)) +
  geom_line() +
  facet_wrap(~ group) +
  theme_bw() +
  theme(legend.position = "none")

ggplotly(AP, tooltip = "sample")

HS <- ggplot(prmtr.means, aes(position, ratio, color = sample)) +
  geom_line() +
  facet_wrap(~ group) +
  theme_bw() +
  theme(legend.position = "none")

ggplotly(HS, tooltip = "sample")

HS2 <- ggplot(prmtr.means, aes(position, ratio, color = sample)) +
  geom_line() +
  facet_wrap(~ group) +
  theme_bw() +
  theme(legend.position = "none")

ggplotly(HS2, tooltip = "sample")

HB <- ggplot(prmtr.means, aes(position, ratio, color = sample)) +
  geom_line() +
  facet_wrap(~ group) +
  theme_bw() +
  theme(legend.position = "none")

ggplotly(HB, tooltip = "sample")

LA <- ggplot(prmtr.means, aes(position, ratio, color = sample)) +
  geom_line() +
  facet_wrap(~ group) +
  theme_bw() +
  theme(legend.position = "none")

ggplotly(LA, tooltip = "sample")

#save to plotly website
Sys.setenv("plotly_username" = "judyanncocadiz")
Sys.setenv("plotly_api_key" = "Op5slVPOQ7w2opv9D1ov")

api_create(AP, "All_Promoters")
api_create(HS, "High_Stad")
api_create(HS2, "High_Stad2")
api_create(HB, "High_Blood")
api_create(LA, "Low_All")
