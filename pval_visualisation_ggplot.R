#<<<<<<<<<<<<<<<<<HEAD 

# load packages
library(ggpubr)
library(rstatix)

head(ToothGrowth)

# calculate adjusted pvalues
stat.test <- ToothGrowth %>%
  group_by(dose) %>%
  t_test(len ~ supp) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

# Create a box plot
bxp <- ggboxplot(
  ToothGrowth, x = "supp", y = "len", facet.by = "dose",
  fill = "supp", palette = c("#00AFBB", "#E7B800")
)

# Add adjusted p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "supp", dodge = 0.8)
bxp + stat_pvalue_manual(
  stat.test,  label = "p.adj", tip.length = 0
)

# Add 10% spaces between the adjusted p-value labels and the plot border
bxp + stat_pvalue_manual(
  stat.test,  label = "p.adj", tip.length = 0
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Only show stars for adjusted pvals
bxp + stat_pvalue_manual(
  stat.test,  label = "p.adj.signif", tip.length = 0
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) 

