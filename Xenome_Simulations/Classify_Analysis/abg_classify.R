abg_classify_graft_1 <- read.csv("~/Xenograft_RStudio/Xenome_Simulations/Classify_Analysis/abg_classify_graft_1.csv")
abg_classify_graft_2 <- read.csv("~/Xenograft_RStudio/Xenome_Simulations/Classify_Analysis/abg_classify_graft_2.csv")
abg_classify_total_1 <- read.csv("~/Xenograft_RStudio/Xenome_Simulations/Classify_Analysis/abg_classify_total_1.csv")
abg_classify_total_2 <- read.csv("~/Xenograft_RStudio/Xenome_Simulations/Classify_Analysis/abg_classify_total_2.csv")

total_graft <- sum(abg_classify_graft_1$total_graft_reads + abg_classify_graft_2$total_graft_reads)
total_g_in_g <- sum(abg_classify_graft_1$graft_in_graft + abg_classify_graft_2$graft_in_graft)
total_g_in_b <- sum(abg_classify_graft_1$graft_in_both +  abg_classify_graft_2$graft_in_both)
total_g_in_a <- sum(abg_classify_graft_1$graft_in_ambiguous +  abg_classify_graft_2$graft_in_ambiguous)
total_g_in_n <- sum(abg_classify_graft_1$graft_in_neither +  abg_classify_graft_2$graft_in_neither)
total_g_in_h <- sum(abg_classify_graft_1$graft_in_host +  abg_classify_graft_2$graft_in_host)

total_graft==total_g_in_a+total_g_in_b+total_g_in_g+total_g_in_n+total_g_in_h

graft_destination <- data.frame(
    group = c("Graft", "Both", "Ambiguous", "Neither", "Host"),
    reads = c(total_g_in_g,total_g_in_b,total_g_in_a,total_g_in_n,total_g_in_h)
)
head(graft_destination)

library(ggplot2)
# Barplot 
graft_div_bp <- ggplot(graft_destination, aes(x="", y=reads, fill = group)) + geom_bar(width = 1, stat = "identity")
graft_div_bp

#PieChart
graft_div_pie <- graft_div_bp + coord_polar("y", start=0)
graft_div_pie


total_by_cg<-c(total_g_in_g,total_g_in_b,total_g_in_a,total_g_in_n,total_g_in_h)
total_by_cg_perc<- 100 *total_by_cg/total_graft
total_by_cg_perc

graft_by_perc <- data.frame(
    sim_no = abg_classify_graft_1$sim.no.,
    total_graft = c(rep(100, 18)), 
    graft_fraction = c(0.20, 0.45, 0.70, 0.25, 0.50, 0.75,0.30,0.55,0.80,0.35,0.50,0.60,0.40,0.85,0.65,0.90,0.95,1.00),
    g_in_g = abg_classify_graft_1$graft_in_graft / abg_classify_graft_1$total_graft_reads, 
    g_in_b = abg_classify_graft_1$graft_in_both / abg_classify_graft_1$total_graft_reads, 
    g_in_a = abg_classify_graft_1$graft_in_ambiguous / abg_classify_graft_1$total_graft_reads, 
    g_in_n = abg_classify_graft_1$graft_in_neither / abg_classify_graft_1$total_graft_reads, 
    g_in_h = abg_classify_graft_1$graft_in_host / abg_classify_graft_1$total_graft_reads 
)

#line plots

p1<- ggplot(graft_by_perc, aes(x=graft_fraction, y=g_in_g)) + geom_point() + geom_smooth(method=lm) + ggtitle("Percent of Graft Reads classified in GRAFT bucket")
p2<- ggplot(graft_by_perc, aes(x=graft_fraction, y=g_in_b)) + geom_point() + geom_smooth(method=lm) + ggtitle("Percent of Graft Reads classified in BOTH bucket")
p3<- ggplot(graft_by_perc, aes(x=graft_fraction, y=g_in_a)) + geom_point() + geom_smooth(method=lm) + ggtitle("Percent of Graft Reads classified in AMBIGUOUS bucket")
p4<- ggplot(graft_by_perc, aes(x=graft_fraction, y=g_in_n)) + geom_point() + geom_smooth(method=lm) + ggtitle("Percent of Graft Reads classified in NEITHER bucket")
p5<- ggplot(graft_by_perc, aes(x=graft_fraction, y=g_in_h)) + geom_point() + geom_smooth(method=lm) + ggtitle("Percent of Graft Reads classified in HOST bucket")

# Trial number Graft Fraction 0.20 - 0.90 - read length 100pm. total read count 60 mil


library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, nrow = 3, ncol = 2)

# by perc of each bucket
graft_by_perc_cg <- data.frame(
    sim_no = abg_classify_graft_1$sim.no.,
    total_graft = c(rep(100, 18)), 
    graft_fraction = c(0.20, 0.45, 0.70, 0.25, 0.50, 0.75,0.30,0.55,0.80,0.35,0.50,0.60,0.40,0.85,0.65,0.90,0.95,1.00),
    g_in_g = abg_classify_graft_1$graft_in_graft / abg_classify_total_1$total_in_graft, 
    g_in_b = abg_classify_graft_1$graft_in_both / abg_classify_total_1$total_in_both, 
    g_in_a = abg_classify_graft_1$graft_in_ambiguous / abg_classify_total_1$total_in_ambiguous, 
    g_in_n = abg_classify_graft_1$graft_in_neither / abg_classify_total_1$total_in_neither, 
    g_in_h = abg_classify_graft_1$graft_in_host / abg_classify_total_1$total_in_host
)

p1<- ggplot(graft_by_perc_cg, aes(x=graft_fraction, y=g_in_g)) + geom_point() + geom_smooth(method=lm) + ggtitle("Percent of Graft Reads in GRAFT bucket vs. Graft Fraction")
p2<- ggplot(graft_by_perc_cg, aes(x=graft_fraction, y=g_in_b)) + geom_point() + geom_smooth(method=lm) + ggtitle("Percent of Graft Reads in BOTH bucket vs. Graft Fraction")
p3<- ggplot(graft_by_perc_cg, aes(x=graft_fraction, y=g_in_a)) + geom_point() + geom_smooth(method=lm) + ggtitle("Percent of Graft Reads in AMBIGUOUS bucket vs. Graft Fraction")
p4<- ggplot(graft_by_perc_cg, aes(x=graft_fraction, y=g_in_n)) + geom_point() + geom_smooth(method=lm) + ggtitle("Percent of Graft Reads in NEITHER bucket vs. Graft Fraction")
p5<- ggplot(graft_by_perc_cg, aes(x=graft_fraction, y=g_in_h)) + geom_point() + geom_smooth(method=lm) + ggtitle("Percent of Graft Reads in HOST bucket vs. Graft Fraction")

library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, nrow = 3, ncol = 2)

