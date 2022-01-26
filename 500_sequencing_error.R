######################################################################################
#
# Calculate the sequencing error on forward/reverse primers
# Date: 2021-01-04
#


forward <- read_tsv("../imported/v_rearranged/sequencing_error/forward.tab", col_names = FALSE) %>%
  select(1:3) %>% set_names(c("sample", "count", "seq"))
reverse <- read_tsv("../imported/v_rearranged/sequencing_error/reverse.tab", col_names = FALSE) %>%
  select(1:3) %>% set_names(c("sample", "count", "seq"))



fwd_seq <- forward %>% mutate(correct = grepl("GAGGACA(C|T)(G|A|T)GCC(G|A|C)TGTA(C|T)(A|T)(A|C)CTGT", seq)) 
rev_seq <- reverse %>% mutate(correct = grepl("GGAACAGTTCAACAGCACGTACC", seq)) 

fwd_error_rate <- fwd_seq %>%
  group_by(sample, correct) %>% summarise(count = sum(count)) %>% 
  group_by(sample) %>% add_tally(count, name = "total_count") %>%
  subset(!correct) %>%
  mutate(proportion = count / total_count) %>%
  mutate(base_error_rate = (1 - proportion) ^ (1/23)) %>%
  mutate(primer = "forward")
rev_error_rate <- rev_seq %>%
  group_by(sample, correct) %>% summarise(count = sum(count)) %>% 
  group_by(sample) %>% add_tally(count, name = "total_count") %>%
  subset(!correct) %>%
  mutate(proportion = count / total_count) %>%
  mutate(base_error_rate = (1 - proportion) ^ (1/24)) %>%
  mutate(primer = "reverse")

error_rate <- bind_rows(fwd_error_rate, rev_error_rate)
save(error_rate, file = "data/sequencing_error.Rda")


p_erate_1 <- ggplot(error_rate, aes(x = sample, y = 1 - base_error_rate, color = primer)) +
  geom_line(aes(group = primer)) +
  # scale_y_continuous(limits = c(0, 0.005)) +
  theme_classic() +
  theme(legend.position = "top") +
  labs(x = "", y = "Average base error rate", title = "", color = "Primer")
p_erate_1

p_erate_2 <- ggplot(error_rate, aes(x= primer, y = base_error_rate)) +
  geom_violin() +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, aes(fill = primer)) +
  # scale_y_continuous(limits = c(0.998, 1)) +
  theme_classic() +
  # theme(axis.ticks.x = element_blank()) +
  labs(x = "", y = "Average base error rate", title = "", fill = "Primer")
p_erate_2



pdf("results/figure_C4_error_rate.pdf", width = 4, height = 4)
p_erate_1
p_erate_2
dev.off()



sites <- 0:70
prob <- 1 - min(error_rate$base_error_rate)
prob <- data.frame(site = sites, prob = dbinom(sites, size = 70, prob))

ggplot(prob, aes(x = site, y = prob)) +
  geom_bar(stat = "identity", width = 0.1) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 70, 10)) +
  scale_y_continuous(breaks = seq(0, 0.8, 0.1)) +
  labs(x = "site", y = "Probability")

