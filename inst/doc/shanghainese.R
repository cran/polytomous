### R code from vignette source 'shanghainese.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: data
###################################################
library(polytomous)
data(shanghainese)
head(shanghainese)

summary(shanghainese)


###################################################
### code chunk number 2: FUNCTION.INTR crosstabulation
###################################################
table(shanghainese$FUNCTION=="INTR", shanghainese$TOPIC_MARKER)


###################################################
### code chunk number 3: FUNCTION.INTR topic-marker-wise proportions
###################################################
round(table(shanghainese$FUNCTION=="INTR", shanghainese$TOPIC_MARKER)[2,]*100/100)


###################################################
### code chunk number 4: FUNCTION.INTR feature-wise proportions
###################################################
round(table(shanghainese$FUNCTION=="INTR", shanghainese$TOPIC_MARKER)[2,]*100/197)


###################################################
### code chunk number 5: FUNCTION.INTR chisq.test (eval = FALSE)
###################################################
## chisq.test(table(shanghainese$FUNCTION=="INTR", shanghainese$TOPIC_MARKER))[c("statistic",
## "parameter","p.value")]


###################################################
### code chunk number 6: FUNCTION.INTR associations
###################################################
associations(table(shanghainese$FUNCTION=="INTR",
shanghainese$TOPIC_MARKER))[c("uc.RC","uc.CR")]


###################################################
### code chunk number 7: FUNCTION.INTR chisq.posthoc std.pearson.residuals
###################################################
chisq.posthoc(table(shanghainese$FUNCTION=="INTR",
shanghainese$TOPIC_MARKER))$cells$std.pearson.residuals


###################################################
### code chunk number 8: FUNCTION.INTR alternative
###################################################
chisq.test(table(shanghainese$FUNCTION=="INTR", shanghainese$TOPIC_MARKER))$stdres


###################################################
### code chunk number 9: FUNCTION.INTR chisq.posthoc signs
###################################################
chisq.posthoc(table(shanghainese$FUNCTION=="INTR", shanghainese$TOPIC_MARKER),
std.pearson.residual.min=2)$cells$std.pearson.residuals.sign


###################################################
### code chunk number 10: FUNCTION/TOPIC_MARKER crosstabulation
###################################################
table(shanghainese$FUNCTION, shanghainese$TOPIC_MARKER)


###################################################
### code chunk number 11: FUNCTION chisq.test
###################################################
chisq.test(table(shanghainese$FUNCTION, shanghainese$TOPIC_MARKER))$p.value


###################################################
### code chunk number 12: FUNCTION associations
###################################################
associations(table(shanghainese$FUNCTION,
shanghainese$TOPIC_MARKER))[c("uc.RC","uc.CR")]


###################################################
### code chunk number 13: FUNCTION univariate summary results
###################################################
chisq.posthoc(table(shanghainese$FUNCTION,
shanghainese$TOPIC_MARKER))$cells$std.pearson.residuals.sign


###################################################
### code chunk number 14: univariate table creation
###################################################

shanghainese$TOPIC_LENGTH <- factor(shanghainese$TOPIC_LENGTH)

shanghainese.logical <- multinomial2logical(shanghainese, outcome="TOPIC_MARKER",
variables=c("TOPIC_LENGTH", "TOPIC_POS", "FUNCTION", "COMMENT_TYPE", "GENRE"),
variable.value.separator=".")

shanghainese.univariate <- nominal(TOPIC_MARKER ~ ., data=shanghainese.logical)


###################################################
### code chunk number 15: univariate complete results
###################################################
print(summary(shanghainese.univariate), max.print=NA)


###################################################
### code chunk number 16: ANOVA
###################################################
shanghainese$TOPIC_LENGTH <- as.numeric(as.character(shanghainese$TOPIC_LENGTH))

summary(aov(TOPIC_LENGTH ~ TOPIC_MARKER, data=shanghainese))


###################################################
### code chunk number 17: Kruskal-Wallis
###################################################
kruskal.test(TOPIC_LENGTH ~ TOPIC_MARKER, data=shanghainese)


###################################################
### code chunk number 18: TukeyHSD
###################################################
TukeyHSD(aov(TOPIC_LENGTH ~ TOPIC_MARKER, data=shanghainese))


###################################################
### code chunk number 19: topic-length-means
###################################################
sapply(levels(shanghainese$TOPIC_MARKER),
   function(i) mean(shanghainese$TOPIC_LENGTH[shanghainese$TOPIC_MARKER==i]))


###################################################
### code chunk number 20: Wilcox-tests
###################################################
wilcox.test(shanghainese$TOPIC_LENGTH[shanghainese$TOPIC_MARKER=="mo"],
shanghainese$TOPIC_LENGTH[shanghainese$TOPIC_MARKER=="a"])

wilcox.test(shanghainese$TOPIC_LENGTH[shanghainese$TOPIC_MARKER=="a"],
shanghainese$TOPIC_LENGTH[shanghainese$TOPIC_MARKER=="ne"])

wilcox.test(shanghainese$TOPIC_LENGTH[shanghainese$TOPIC_MARKER=="ne"],
shanghainese$TOPIC_LENGTH[shanghainese$TOPIC_MARKER=="ma"])

wilcox.test(shanghainese$TOPIC_LENGTH[shanghainese$TOPIC_MARKER=="ma"],
shanghainese$TOPIC_LENGTH[shanghainese$TOPIC_MARKER=="zi"])


###################################################
### code chunk number 21: bivariate basic crosstabulation
###################################################
table(shanghainese.logical[["TOPIC_LENGTH.1"]], shanghainese.logical[["FUNCTION.INTR"]])


###################################################
### code chunk number 22: FUNCTION associations
###################################################
associations(table(shanghainese.logical[["TOPIC_LENGTH.1"]],
shanghainese.logical[["FUNCTION.INTR"]]))[c("uc.RC","uc.CR")]


###################################################
### code chunk number 23: bivariate table creation
###################################################
shanghainese.bivariate <- nominal(. ~ ., data=shanghainese.logical[-1])
summary(shanghainese.bivariate)


###################################################
### code chunk number 24: bivariate uc greater than .3
###################################################
subset(summary(shanghainese.bivariate)$sumry.table, uc.12>.3 | uc.21>.3)


###################################################
### code chunk number 25: polytomous model fit
###################################################
shanghainese$TOPIC_LENGTH <- as.numeric(shanghainese$TOPIC_LENGTH)

shanghainese.polytomous <- polytomous(TOPIC_MARKER ~ TOPIC_LENGTH + TOPIC_POS +
FUNCTION + COMMENT_TYPE + GENRE, data=shanghainese)


###################################################
### code chunk number 26: polytomous model summary
###################################################
print(summary(shanghainese.polytomous), max.print=NA)


###################################################
### code chunk number 27: shanghainese factor relevel
###################################################
levels(shanghainese$FUNCTION)
shanghainese$FUNCTION <- relevel(shanghainese$FUNCTION,"EMPH")
levels(shanghainese$FUNCTION)


###################################################
### code chunk number 28: polytomous model performance
###################################################
shanghainese.polytomous$statistics$R2.likelihood

shanghainese.polytomous$statistics$accuracy


###################################################
### code chunk number 29: polytomous recall
###################################################
shanghainese.polytomous$statistics$recall.predicted


###################################################
### code chunk number 30: polytomous model prediction crosstabulation
###################################################
shanghainese.polytomous$statistics$crosstable


###################################################
### code chunk number 31: polytomous model predicted probabilities
###################################################
head(shanghainese.polytomous$fitted)


###################################################
### code chunk number 32: probs figure (eval = FALSE)
###################################################
## plot(shanghainese.polytomous, values="probabilities", panes="multiple")


###################################################
### code chunk number 33: polytomous model sorted probabilities
###################################################
apply(shanghainese.polytomous$fitted,1,sort)[,1:5]


###################################################
### code chunk number 34: polytomous model sorted probabilities two highest
###################################################
apply(shanghainese.polytomous$fitted,1,sort)[4:5,1:5]


###################################################
### code chunk number 35: number of two-ways interchangeable contexts
###################################################
length(which(apply(apply(shanghainese.polytomous$fitted,1,sort)[4:5,],2,function(x) all(x>0.3))))


###################################################
### code chunk number 36: indices of best exemplars
###################################################
apply(shanghainese.polytomous$fitted, 2, function(x) order(x, decreasing=T))[1:5,]


###################################################
### code chunk number 37: exemplar probabilities
###################################################
round(shanghainese.polytomous$fitted[222,],3)

round(shanghainese.polytomous$fitted[185,],3)

round(shanghainese.polytomous$fitted[40,],3)


###################################################
### code chunk number 38: most interchangeable cases
###################################################
order(apply(shanghainese.polytomous$fitted,1,sd))[1:5]


###################################################
### code chunk number 39: probabilities for interchageable case
###################################################
round(shanghainese.polytomous$fitted[94,],3)


