library(ggplot2)
# data = read.table("Hek_distribution.txt", sep="\t", header=TRUE)
data = read.table("combined_PSU_freq_per_gene.txt", sep="\t", header=TRUE)
View(data)
ggplot(data, aes(x = PSU_frequency_per_gene, fill = Cell_line)) +
  geom_density(alpha = .3)+ 
  xlab("Psi frequency per gene") + ylab("Density")

ggplot(data, aes(x = PSU_frequency_per_gene, fill = Cell_line)) +
  geom_density(stat = "count", position = "stack")+ 
  xlab("Number of genes") + ylab(expression(psi~"frequency distribution"))

ggplot(data, aes(x = PSU_frequency_per_gene, fill = Cell_line)) +
  geom_density(stat = "count", position = "stack", adjust= 300)+ 
  xlab("Number of genes") + ylab(expression(psi~"frequency distribution"))
