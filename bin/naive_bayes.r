#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     02.07.2017
#modificado: 02.07.2017

library("e1071")

trainset <- read.csv("../data/golf.csv",header=TRUE)
testset <- read.csv("../data/golf-testset.csv",header=TRUE)

mod <- naiveBayes(Play ~ Outlook + Wind,      #model design
                  trainset,                   #training set    
                  laplace=0)                  #Laplace correction
res <- predict(mod,                           #model
               testset[,c("Outlook","Wind")]) #testing set

#Confusion Matrix
conf_mat <- table(testset[,c("Play")],res)
conf_mat

#Error rate
err_rt <- 100*(sum(conf_mat) - sum(diag(conf_mat)))/sum(conf_mat)
err_rt
