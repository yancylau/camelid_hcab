#!/usr/bin/env Rscript
#


v_usage <- v_cdna %>% 
  select(v_call) %>% distinct() %>%
  separate(v_call, sep = "_", into = c(NA, "type", NA), remove = F) %>%
  group_by(type) %>% distinct() %>% summarise(count = n()) %>%
  mutate(total = ifelse(type == "vh", 115, 55)) %>%
  mutate(used = count / total) %>% 
  mutate(unused = 1 - used) %>%
  gather(key = "usage", value = "prop", -c(type, count, total))


p_v_usage <- ggplot(v_usage, aes(x = "", y = prop, fill = usage)) +
  geom_bar(stat = "identity") + 
  geom_text(aes(y = prop, label =  round(prop, 4)), position = position_fill(reverse = F, vjust = 0.5)) +
  coord_polar(theta = "y") +
  theme_void() +
  facet_grid(. ~ type)


pdf("results/v_cdna/v_gene_usage.pdf")
p_v_usage
dev.off()
