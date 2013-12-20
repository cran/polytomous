### R code from vignette source 'exemplars2prototypes.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: data
###################################################
library(polytomous)
data(think)


###################################################
### code chunk number 2: logical
###################################################
names(think)

think.logical <- multinomial2logical(data=think, outcome="Lexeme", variables=names(think)[2:24])


###################################################
### code chunk number 3: formula
###################################################
names(think.logical)[-1]


###################################################
### code chunk number 4: formula2
###################################################
grep("(Other)|(None)|(FiniteVerbChain)|(Overt)",
   names(think.logical)[-1],value=T,invert=T)


###################################################
### code chunk number 5: formula3
###################################################
paste(grep("(Other)|(None)|(FiniteVerbChain)|(Overt)",
   names(think.logical)[-1],value=T,invert=T),collapse="")
paste("Lexeme",paste(grep("(Other)|(None)|(FiniteVerbChain)|(Overt)",
   names(think.logical)[-1],value=T,invert=T),collapse=" + "),sep=" ~ ")
as.formula(paste("Lexeme",
   paste(grep("(Other)|(None)|(FiniteVerbChain)|(Overt)",
   names(think.logical)[-1],value=T,invert=T),collapse=" + "),sep=" ~ "))
think.formula <- as.formula(paste("Lexeme", 
   paste(grep("(Other)|(None)|(FiniteVerbChain)|(Overt)",
   names(think.logical)[-1],value=T,invert=T),collapse=" + "),sep=" ~ "))


###################################################
### code chunk number 6: polytomous
###################################################
think.polytomous <- polytomous(think.formula, data=think.logical)

print(summary(think.polytomous),max.print=NA)


###################################################
### code chunk number 7: exemplars
###################################################
extract.exemplars(think.polytomous, n.clusters=100)


###################################################
### code chunk number 8: prototypes
###################################################
extract.prototypes(think.polytomous)


