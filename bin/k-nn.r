#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     02.07.2017
#modificado: 02.07.2017

library("class")

trainset <- read.csv("../data/golf.csv",header=TRUE)
testset <- read.csv("../data/golf-testset.csv",header=TRUE)

res <- knn(trainset[,c("Temperature","Humidity")], #attributes of training set
           testset[,c("Temperature","Humidity")],  #testing set
           trainset[,"Play"],                      #classes of trainset
           k=1)                                    #k-nearest neighbour

#Confusion Matrix
conf_mat <- table(testset[,c("Play")],res)
conf_mat

#Error rate
err_rt <- 100*(sum(conf_mat) - sum(diag(conf_mat)))/sum(conf_mat)
err_rt
