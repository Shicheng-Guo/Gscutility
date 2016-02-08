#######################################################################################################################
###   Title : ggplot2 for barplot with error-bar and apply Arial font in PDF output
###   Lab:  Kun's Lab (http://zhang.openwetware.org/)   
###   Time :  Feb/29/2016
#######################################################################################################################

n.per.group<-10
alpha<-0.05 # for a (1.00-alpha)=95% confidence interval
# Simulate raw data for an experiment or observational study.
data.raw <- data.frame(
  treatment=rep(c('A','B'), each=n.per.group),
  value=c(rnorm(n.per.group, 2), rnorm(n.per.group, 3))   
)
# This data frame calculates statistics for each treatment.
data.summary <- data.frame(
  treatment=levels(data.raw$treatment),
  mean=tapply(data.raw$value, data.raw$treatment, mean),
  n=tapply(data.raw$value, data.raw$treatment, length),
  sd=tapply(data.raw$value, data.raw$treatment, sd)
)
# Precalculate standard error of the mean (SEM)
data.summary$sem <- data.summary$sd/sqrt(data.summary$n)
# Precalculate margin of error for confidence interval
data.summary$me <- qt(1-alpha/2, df=data.summary$n)*data.summary$sem

# install and load Arial for PDF
install.packages("extrafont")
library("extrafont")
font_import()
fonts()
loadfonts() # loadfonts(device = "postscript")

# Load ggplot2 library
require(ggplot2)
# Use ggplot to draw the bar plot using the precalculated 95% CI.
pdf("barplot.ggplot2.pdf", height = 3, width = 3, family = "Arial", paper = "special", onefile = FALSE)
ggplot(data.summary, aes(x = treatment, y = mean)) +  
  geom_bar(position = position_dodge(), stat="identity", fill="blue") + 
  geom_errorbar(aes(ymin=mean-me, ymax=mean+me)) +
  ggtitle("Bar plot with 95% confidence intervals") + # plot title
  theme_bw() + # remove grey background (because Tufte said so)
  theme(panel.grid.major = element_blank()) # remove x and y major grid lines (because Tufte said so)
dev.off() # Close PNG

# Plot one standard error (standard error of the mean/SEM)
pdf("barplot.ggplot2.pdf", height = 3, width = 3, family = "Arial", paper = "special", onefile = FALSE)
ggplot(data.summary, aes(x = treatment, y = mean)) +  
  geom_bar(position = position_dodge(), stat="identity", fill="blue") + 
  geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem)) +
  ggtitle("Bar plot with standard error as error bars") + 
  theme_bw() +
  theme(panel.grid.major = element_blank())
dev.off()

# Plot one standard deviation
pdf("barplot.ggplot2.pdf", height = 3, width = 3, family = "Arial", paper = "special", onefile = FALSE)
ggplot(data.summary, aes(x = treatment, y = mean)) +  
  geom_bar(position = position_dodge(), stat="identity", fill="blue") + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd)) +
  ggtitle("Bar plot with standard deviation as error bars") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank())
dev.off()
