install.packages("mlbench")
library(mlbench)
?BreastCancer
data(BreastCancer)
dim(BreastCancer)
sum(is.na(BreastCancer))
head(BreastCancer)

BreastCancer<-na.omit(BreastCancer)
dim(BreastCancer)
head(BreastCancer)
sum(is.na(BreastCancer))


t1 = table(BreastCancer$Class)
names(dimnames(t1)) = "Class"
print(t1)

pairdata<-data.frame(BreastCancer$Cl.thickness, BreastCancer$Cell.size, BreastCancer$Cell.shape,
          BreastCancer$Marg.adhesion, BreastCancer$Epith.c.size, BreastCancer$Bare.nuclei, 
          BreastCancer$Bl.cromatin, BreastCancer$Normal.nucleoli, BreastCancer$Mitoses)
          

pairs(pairdata)
plot(BreastCancer$Cell.shape, BreastCancer$Cell.size)

#some correlation between
#Cell.size and Cell.shape

#linear regression on all cont. variables

lm1 = lm(as.numeric(Class == "malignant") ~ Cl.thickness+Cell.size+Cell.shape+
         Marg.adhesion+Epith.c.size+Bare.nuclei+Bl.cromatin+Normal.nucleoli+Mitoses,
         data=BreastCancer)
summary(lm1)

#most signif factors are thickness, bare nuclei, bl.chromatin and normal nucleoli

#lm on these four

lm2 = lm(as.numeric(Class == "malignant") ~ Cl.thickness+Bare.nuclei+Bl.cromatin+Normal.nucleoli,
         data = BreastCancer)
summary(lm2)

#not terrible, collinearity?

library(tidyverse)
library(carData)
library(caret)

car::vif(lm2)

#no major issues with these variables

#logistic regression on one variable

BreastCancer$numClass = as.numeric(BreastCancer$Class == "malignant")
glm1 = glm(numClass ~ Cl.thickness, data = BreastCancer , family = "binomial")

summary(glm1)
summary(glm1$fitted.values)

plot(BreastCancer$Cl.thickness, as.numeric(BreastCancer$Class == "malignant"),col="red")
points(glm1$data$Cl.thickness, glm1$fitted.values, col = "black", pch = 4)
curve(predict(glm1,data.frame(Cl.thickness = x),type="resp"),col="blue",lwd=2,add=TRUE)

#curve - not sure what is happening

BreastCancer$numClass = as.numeric(BreastCancer$Class == "malignant")
glm2 = glm(numClass ~ Bare.nuclei, data = BreastCancer , family = "binomial")

summary(glm2)
summary(glm2$fitted.values)

plot(BreastCancer$Bare.nuclei, as.numeric(BreastCancer$Class == "malignant"),col="red")
points(glm2$data$Bare.nuclei, glm2$fitted.values, col = "black", pch = 4)
curve(predict(glm1,data.frame(Bare.nuclei = x),type="resp"),col="blue",lwd=2,add=TRUE)

#this curve stuff isn't working but that's ok as i'm going to move on to LDN

#logistic regression with multiple variables (4 most significant)

glm3<-(glm(numClass ~ Cl.thickness+Bare.nuclei+Bl.cromatin+Normal.nucleoli, data = BreastCancer,
           family = "binomial"))


#now getting an error 'glm.fit: fitted probabilities numerically 0 or 1 occurred'
#perfect selection? this probably signals instability
#also talk about low variation (some classes only have 5 or 6 possible values?)

summary(glm3)

#which is a GOOD REASON TO TRY linear discriminant analysis
#LDA
#four variables

library(MASS)

ld1<-lda(Class~ Cl.thickness+Bare.nuclei+Bl.cromatin+Normal.nucleoli,
         data = BreastCancer, na.action = na.omit)

ld1

ld1_pred = predict(ld1)
ld1_pred_class = ld1_pred$class
table(BreastCancer$Class, ld1_pred_class)

mean(ld1_pred_class != BreastCancer$Class)

#perfectly good model, ostensibly
#let's validate it
#sample size is quite small so i will start with LOOCV and then I might bootstrap it


#neither of these attempts at LOOCV is working

#K-fold CV on the glm

BreastCancer$numClass<- as.numeric(BreastCancer$Class == "malignant")

K = 5
rnd = sample(1:nrow(BreastCancer),nrow(BreastCancer))
N = nrow(BreastCancer)
nk = rep(floor(N/K),K)
leftover = N - K*floor(N/K)
if (leftover > 0){
  nk[1:leftover] = nk[1:leftover]+1
}
shuffleBC = BreastCancer[rnd,]
shuffleBC$Kfold = rep(1:K,nk)

g_vec = list()
pred_vec = list()
mse_vec = rep(NA, K)
for (i in 1:length(mse_vec)){
  training = subset(shuffleBC, Kfold != i)
  testing = subset(shuffleBC, Kfold == i)
  g1 = glm(numClass ~ Cl.thickness+Bare.nuclei+Bl.cromatin+Normal.nucleoli, data = training,
           family = "binomial")
  pred = predict(g1, testing, type="response")
  g_vec[[i]] = g1
  pred_vec[[i]] = pred
  mse_vec[i] = mean((testing$numClass - pred)^2)
  }
sum(mse_vec*nk/N)

#fitted probabilities numerically 0 or 1 occurred

table(ifelse(pred_vec[[1]]>0.5,1,0),testing$numClass)

#now it doesn't want to make the table because the lengths of numClass and pred_vec differ by 1

length(testing$numClass)
length(pred_vec)
length(pred_vec[2])
length(ifelse(pred_vec[[1]]>0.5,1,0))

#i don't know how to fix this short of appending an extra value to numClass

numClassNew<-append(testing$numClass, 0, after=length(testing$numClass))
length(numClassNew)
numClassNew       

table(ifelse(pred_vec[[1]]>0.5,1,0),numClassNew)

#not a good outcome, and obviously the data is biased too now

#let's try this again with the lda


#it's still giving the non-numeric argument to binary operator error
#i don't know how to fix this










