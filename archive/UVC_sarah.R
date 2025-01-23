library(tidyverse)
library(readxl)


ctrl <- read_csv("cumulativeResults-nbl01.csv")
ctrl$type <- "control"
trt <- read_csv("cumulativeResults-nbl07.csv")
trt$type <- "300J_m2_254"
trt2 <- read_csv("cumulativeResults-nbl09.csv")
trt2$type <- "2000J_m2_254"

dat <- rbind(ctrl[,c(7,11)], trt[,c(7, 10)])
dat <- rbind(dat, trt2[,c(7,10)])
dat$type <- factor(dat$type)

p1 <- dat %>%
  ggplot(aes(x = `Mean-G`, fill = type)) +
  geom_density(alpha = 0.2) +
  # geom_histogram(aes(fill = type), alpha = 0.5, bins = 50, color = "black") +
  ggsci::scale_fill_jco() +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,2)) + 
  geom_vline(xintercept = quantile(ctrl$`Mean-G`, probs = 0.999), lty = 2)

p2 <- dat %>%
  ggplot(aes(x = `Mean-G`, fill = type)) +
  geom_density(alpha = 0.2) +
  # geom_histogram(aes(fill = type), alpha = 0.5, bins = 50, color = "black") +
  ggsci::scale_fill_jco() +
  scale_x_continuous(limits = c(10,100), breaks = seq(10,100,10))

ggpubr::ggarrange(p1, p2, nrow = 1)



