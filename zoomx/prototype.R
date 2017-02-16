Somatic=matrix(c(0, 40, 3, 53),
       nrow = 2,
       dimnames = list(Tumor = c("ft", "bg"),
                       Normal = c("ft", "bg")))
fisher.test(Somatic, alternative = "greater")

binom.test(3,53,(1/41),alternative="greater")
binom.test(49,190,(1/201),alternative="greater")
binom.test(7950,144004,(7836/143309),alternative="greater")

# Statistic Looks OK

# Let's 

#  G\T M  T
#  M   3  1
#  T   1  3

# P(T|N)