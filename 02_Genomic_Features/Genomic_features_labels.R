library(ggrepel)

low_in_both_percents <- readRDS(paste0(output, "/short_vs_long_percents/low_in_both_pros_2.RDS"))

low_in_both_percents <- low_in_both_percents %>% pivot_longer(cols = c(percent_shorts_in_gf,
                                                         percent_longs_in_gf),
                                                names_to = "short.or.long",
                                                values_to = "percent_of_fragments_low_promoters")

low <- low_in_both_percents %>%
  mutate(group = case_when(substr(sample,nchar(sample),nchar(sample)) == "e" ~ "high_gc",
                           substr(sample,nchar(sample),nchar(sample)) == "t" ~ "low_gc",
                           substr(sample,1,1) == "P" ~ "delfi_gc",
                           substr(sample,1,1) == "E" ~ "healthy"))
low


lp <- ggplot(low, aes(x = short.or.long, y = percent_of_fragments_low_promoters)) + 
  geom_point(aes(color = sample)) +
  geom_line(aes(group = sample, color = sample)) +
  facet_wrap(~group, nrow = 1) +
  scale_x_discrete(name = "Short or Long", labels=c("percent_shorts_in_gf" = "% Shorts", "percent_longs_in_gf" = "% Longs")) +
  scale_y_continuous(name = "% of fragments") +
  theme(legend.position = "none") +
  ggtitle("% of fragments overlapping promoters low in STAD & blood")
lp

data_ends <- low %>% filter(short.or.long == "percent_shorts_in_gf")
data_ends
data_start <- low %>% filter(short.or.long == "percent_longs_in_gf")
data_start

lp + 
  geom_label_repel(aes(label = sample, color = sample), data = data_ends, size = 1.5, hjust = "outward", 
                  min.segment.length = 0, nudge_x = 2.5, direction = "y") +
  geom_label_repel(aes(label = sample, color = sample), data = data_start, size = 1.5, hjust = "outward", 
                  min.segment.length = 0, nudge_x = -2.5, direction = "y")

ggsave(paste0(output, "/plots3/low_split.png"), last_plot())