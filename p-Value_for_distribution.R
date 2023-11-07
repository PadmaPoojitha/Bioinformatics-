

library(ggpubr)
library(rstatix)
library(ggplot2)
# data = read.table("Hek_distribution.txt", sep="\t", header=TRUE)
data = read.table("3_prime_PSU_distribution.txt", sep="\t", header=TRUE)
#View(data)


ggplot(data, aes(x = PSU_distribution, fill = Cell_line)) +
  geom_density(alpha = .3)+ 
  xlab("Psi distribution across normalized gene length") + ylab("Density")


ggplot(data) +
  aes(x = Cell_line, y = PSU_distribution) +
  theme_minimal() + stat_compare_means(method = "wilcox.test")


bxp <- ggboxplot(data, x = "Cell_line", y = "PSU_distribution", fill = "#00AFBB")
bxp + stat_compare_means()

res <- wilcox.test(PSU_distribution ~ Cell_line, data = data, exact = FALSE)
res

kw = kruskal.test(PSU_distribution ~ Cell_line, data = data) 
kw

test <- chisq.test(table(data$Cell_line, data$PSU_distribution))
test

#################
ggplot(data, aes(x = PSU_distribution, fill = Cell_line)) +
  geom_density(stat = "count", position = "stack")+ 
  xlab(expression(psi~"distribution across normalized gene length")) + ylab("Density")


ggplot(data, aes(x = PSU_distribution, fill = Cell_line)) +
  geom_density(alpha = .3)+ 
  xlab(expression(psi~"distribution across normalized gene length")) + ylab("Density")