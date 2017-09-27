#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     02.07.2017
#modificado: 02.07.2017

library("e1071")

trainset <- read.csv("../data/golf.csv",header=TRUE)
testset <- read.csv("../data/golf-testset.csv",header=TRUE)

#Recoding dummy variables
trainset[,"Outlook"] <- factor(trainset[,"Outlook"])
trainset[,"Wind"] <- factor(trainset[,"Wind"])
tmp1 = model.matrix(~trainset[,"Outlook"]-1)
tmp2 = model.matrix(~trainset[,"Wind"]-1)
colnames(tmp1) = paste("Outlook",levels(trainset[,"Outlook"]),sep="_")
colnames(tmp2) = paste("Wind",levels(trainset[,"Wind"]),sep="_")
trainset = cbind(tmp1,tmp2,trainset[,c("Temperature","Humidity","Play")])

testset[,"Outlook"] <- factor(testset[,"Outlook"])
testset[,"Wind"] <- factor(testset[,"Wind"])
tmp1 = model.matrix(~testset[,"Outlook"]-1)
tmp2 = model.matrix(~testset[,"Wind"]-1)
colnames(tmp1) = paste("Outlook",levels(testset[,"Outlook"]),sep="_")
colnames(tmp2) = paste("Wind",levels(testset[,"Wind"]),sep="_")
testset = cbind(tmp1,tmp2,testset[,c("Temperature","Humidity","Play")])

#Perform SVM
mod <- svm(Play ~ .,                       #model design
                  trainset,                #training set
                  scale=TRUE,              #scaling stored for testing
                  type="C-classification", #SVM algorithm
                  kernel="radial",         #kernel type
                  gamma=1.0,               #kernel gamma
                  cost=1.0,                #C-constant (cost)
                  tolerance=0.001,         #convergence epsilon
                  epsilon=0.0)             #loss-function epsilon
res <- predict(mod,                        #model
               testset)                    #testing set

#Confusion Matrix
conf_mat <- table(testset[,c("Play")],res)
conf_mat

#Error rate
err_rt <- 100*(sum(conf_mat) - sum(diag(conf_mat)))/sum(conf_mat)
err_rt
