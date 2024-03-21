library(ggplot2)
library(plotly)

#EXONS
exon_percents <- readRDS(paste0(output, "/exon2.RDS"))

exon_percents <- exon_percents %>% pivot_longer(cols = c(percent_shorts_in_gf,
                                                         percent_longs_in_gf),
                                                names_to = "short.or.long",
                                                values_to = "percent_of_fragments_overlapping_exons")
exon <- exon_percents %>%
  mutate(group = case_when(substr(sample,nchar(sample),nchar(sample)) == "e" ~ "high_gc",
                           substr(sample,nchar(sample),nchar(sample)) == "t" ~ "low_gc",
                           substr(sample,1,1) == "P" ~ "delfi_gc",
                           substr(sample,1,1) == "E" ~ "healthy"))

e <- ggplot(exon, aes(x=short.or.long,y=percent_of_fragments_overlapping_exons)) +
  geom_point(aes(color = sample)) +
  geom_line(aes(group = sample, color = sample)) +
  facet_wrap(~group, nrow = 1) +
  scale_x_discrete(name = "Short or Long", labels=c("percent_shorts_in_gf" = "% Shorts", "percent_longs_in_gf" = "% Longs")) +
  scale_y_continuous(name = "% of fragments") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("% of fragments overlapping exons")

ggplotly(e, tooltip = "sample")

#PROMOTERS
promoter_percents <- readRDS(paste0(output, "/promoter2.RDS"))

promoter_percents <- promoter_percents %>%
  pivot_longer(cols = c(percent_shorts_in_gf,
                        percent_longs_in_gf),
               names_to = "short.or.long",
               values_to = "percent_of_fragments_overlapping_promoters")
promoter <- promoter_percents %>%
  mutate(group = case_when(substr(sample,nchar(sample),nchar(sample)) == "e" ~ "high_gc",
                           substr(sample,nchar(sample),nchar(sample)) == "t" ~ "low_gc",
                           substr(sample,1,1) == "P" ~ "delfi_gc",
                           substr(sample,1,1) == "E" ~ "healthy"))

p <- ggplot(promoter, aes(x=short.or.long,y=percent_of_fragments_overlapping_promoters)) +
  geom_point(aes(color = sample)) +
  geom_line(aes(group = sample, color = sample)) +
  facet_wrap(~group, nrow = 1) +
  scale_x_discrete(name = "Short or Long", labels=c("percent_shorts_in_gf" = "% Shorts", "percent_longs_in_gf" = "% Longs")) +
  scale_y_continuous(name = "% of fragments") +
  theme_bw() +
  theme(legend.position = "none")+
  ggtitle("% of fragments overlapping promoters")

ggplotly(p, tooltip = "sample")

#5UTRs
fiveUTR_percents <- readRDS(paste0(output, "/fiveUTR_2.RDS"))
fiveUTR_percents <- fiveUTR_percents %>% pivot_longer(cols = c(percent_shorts_in_gf,
                                                               percent_longs_in_gf),
                                                      names_to = "short.or.long",
                                                      values_to = "percent_of_fragments_overlapping_fiveUTRs")
fiveUTR <- fiveUTR_percents %>%
  mutate(group = case_when(substr(sample,nchar(sample),nchar(sample)) == "e" ~ "high_gc",
                           substr(sample,nchar(sample),nchar(sample)) == "t" ~ "low_gc",
                           substr(sample,1,1) == "P" ~ "delfi_gc",
                           substr(sample,1,1) == "E" ~ "healthy"))

f <- ggplot(fiveUTR, aes(x=short.or.long,y=percent_of_fragments_overlapping_fiveUTRs)) +
  geom_point(aes(color = sample)) +
  geom_line(aes(group = sample, color = sample)) +
  facet_wrap(~group, nrow = 1) +
  scale_x_discrete(name = "Short or Long", labels=c("percent_shorts_in_gf" = "% Shorts", "percent_longs_in_gf" = "% Longs")) +
  scale_y_continuous(name = "% of fragments") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("% of fragments overlapping 5UTRs")

ggplotly(f, tooltip = "sample")

#G4s
G4_percents <- readRDS(paste0(output, "/G4_2.RDS"))
G4_percents <- G4_percents %>%
  pivot_longer(cols = c(percent_shorts_in_gf,
                        percent_longs_in_gf),
               names_to = "short.or.long",
               values_to = "percent_of_fragments_overlapping_G4s")
G4 <- G4_percents %>%
  mutate(group = case_when(substr(sample,nchar(sample),nchar(sample)) == "e" ~ "high_gc",
                           substr(sample,nchar(sample),nchar(sample)) == "t" ~ "low_gc", 
                           substr(sample,1,1) == "P" ~ "delfi_gc",
                           substr(sample,1,1) == "E" ~ "healthy"))

g <- ggplot(G4, aes(x=short.or.long,y=percent_of_fragments_overlapping_G4s)) +
  geom_point(aes(color = sample)) +
  geom_line(aes(group = sample, color = sample)) +
  facet_wrap(~group, nrow = 1) +
  scale_x_discrete(name = "Short or Long", labels=c("percent_shorts_in_gf" = "% Shorts", "percent_longs_in_gf" = "% Longs")) +
  scale_y_continuous(name = "% of fragments") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("% of fragments overlapping G4s")

ggplotly(g, tooltip = "sample")

#HIGH IN STAD
high_stad_percents <- readRDS(paste0(output, "/high_stad_pros_2.RDS"))
high_stad_percents <- high_stad_percents %>% pivot_longer(cols = c(percent_shorts_in_gf,
                                                                   percent_longs_in_gf),
                                                          names_to = "short.or.long",
                                                          values_to = "percent_of_fragments_STAD_highly_expressed_promoters")

high_stad <- high_stad_percents %>% mutate(group = case_when(substr(sample,nchar(sample),nchar(sample)) == "e" ~ "high_gc",
                                                             substr(sample,nchar(sample),nchar(sample)) == "t" ~ "low_gc",
                                                             substr(sample,1,1) == "P" ~ "delfi_gc",
                                                             substr(sample,1,1) == "E" ~ "healthy"))


s <- ggplot(high_stad, aes(x = short.or.long, y= percent_of_fragments_STAD_highly_expressed_promoters)) +
  geom_point(aes(color = sample)) +
  geom_line(aes(group = sample, color = sample)) +
  facet_wrap(~group, nrow = 1) +
  scale_x_discrete(name = "Short or Long", labels=c("percent_shorts_in_gf" = "% Shorts", "percent_longs_in_gf" = "% Longs")) +
  scale_y_continuous(name = "% of fragments") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("% of fragments overlapping highly expressed promoters in STAD")

ggplotly(s, tooltip = "sample")

#HIGH IN BLOOD
high_blood_percents <- readRDS(paste0(output, "/short_vs_long_percents/high_blood_pros_2.RDS"))
high_blood_percents <- high_blood_percents %>%
  pivot_longer(cols = c(percent_shorts_in_gf, percent_longs_in_gf),
               names_to = "short.or.long",
               values_to = "percent_of_fragments_blood_highly_expressed_promoters")

high_blood <- high_blood_percents %>%
  mutate(group = case_when(substr(sample,nchar(sample),nchar(sample)) == "e" ~ "high_gc",
                           substr(sample,nchar(sample),nchar(sample)) == "t" ~ "low_gc",
                           substr(sample,1,1) == "P" ~ "delfi_gc",
                           substr(sample,1,1) == "E" ~ "healthy"))

b <- ggplot(high_blood, aes(x=short.or.long,y=percent_of_fragments_blood_highly_expressed_promoters)) +
  geom_point(aes(color = sample)) +
  geom_line(aes(group = sample, color = sample)) +
  facet_wrap(~group, nrow = 1) +
  scale_x_discrete(name = "Short or Long", labels=c("percent_shorts_in_gf" = "% Shorts", "percent_longs_in_gf" = "% Longs")) +
  scale_y_continuous(name = "% of fragments") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("% of fragments overlapping highly expressed promoters in blood")

ggplotly(b, tooltip = "sample")

#LOW IN BOTH
low_in_both_percents <- readRDS(paste0(output, "/short_vs_long_percents/low_in_both_pros_2.RDS"))
low_in_both_percents <- low_in_both_percents %>%
  pivot_longer(cols = c(percent_shorts_in_gf, percent_longs_in_gf),
               names_to = "short.or.long",
               values_to = "percent_of_fragments_low_promoters")

low <- low_in_both_percents %>%
  mutate(group = case_when(substr(sample,nchar(sample),nchar(sample)) == "e" ~ "high_gc",
                           substr(sample,nchar(sample),nchar(sample)) == "t" ~ "low_gc", 
                           substr(sample,1,1) == "P" ~ "delfi_gc",
                           substr(sample,1,1) == "E" ~ "healthy"))

l <- ggplot(low, aes(x=short.or.long,y=percent_of_fragments_low_promoters)) +
  geom_point(aes(color = sample)) +
  geom_line(aes(group = sample, color = sample)) +
  facet_wrap(~group, nrow = 1) +
  scale_x_discrete(name = "Short or Long", labels=c("percent_shorts_in_gf" = "% Shorts", "percent_longs_in_gf" = "% Longs")) +
  scale_y_continuous(name = "% of fragments") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("% of fragments overlapping promoters low in STAD & blood")

ggplotly(l, tooltip = "sample")

#HIGH IN STAD 2

high_stad2_percents <- readRDS(paste0(output, "/high_stad2_pros.RDS"))
high_stad2_percents <- high_stad2_percents %>% pivot_longer(cols = c(percent_shorts_in_gf,
                                                                     percent_longs_in_gf),
                                                            names_to = "short.or.long",
                                                            values_to = "percent_of_fragments_STAD2_highly_expressed_promoters")

high_stad2 <- high_stad2_percents %>% mutate(group = case_when(substr(sample,nchar(sample),nchar(sample)) == "e" ~ "high_gc",
                                                               substr(sample,nchar(sample),nchar(sample)) == "t" ~ "low_gc",
                                                               substr(sample,1,1) == "P" ~ "delfi_gc",
                                                               substr(sample,1,1) == "E" ~ "healthy"))


s2 <- ggplot(high_stad2, aes(x = short.or.long, y= percent_of_fragments_STAD2_highly_expressed_promoters)) +
  geom_point(aes(color = sample)) +
  geom_line(aes(group = sample, color = sample)) +
  facet_wrap(~group, nrow = 1) +
  scale_x_discrete(name = "Short or Long", labels=c("percent_shorts_in_gf" = "% Shorts", "percent_longs_in_gf" = "% Longs")) +
  scale_y_continuous(name = "% of fragments") +
  theme_bw() +
  theme(legend.position = "none") +
  ggtitle("% of fragments overlapping highly expressed promoters in STAD 2")

ggplotly(s2, tooltip = "sample")

#save to plotly website
Sys.setenv("plotly_username" = "judyanncocadiz")
Sys.setenv("plotly_api_key" = "Op5slVPOQ7w2opv9D1ov")

api_create(p, "bw_Promoters")
api_create(e, "bw_Exons")
api_create(f, "bw_5UTRs")
api_create(g, "bw_G4s")
api_create(s, "bw_High in STAD")
api_create(s2, "bw_High in STAD 2")
api_create(b, "bw_High in blood")
api_create(l, "bw_Low in STAD & blood")

