#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     02.07.2017
#modificado: 02.07.2017

library("nnet")
library("caret")
library("NeuralNetTools")

trainset <- read.csv("../data/ripley-set.csv",header=TRUE)
trainset[,"label"] = factor(trainset[,"label"],levels=c("0","1"))

#Create testing set by sampling with replacement
n <- nrow(trainset)
bootindex <- sample(1:n,size=n,replace=TRUE)
testset <- trainset[bootindex,]

#Select settings for artificial neural network
nn_decay <- 0.3  #weight decay (improved version of learning rate)
nn_maxit <- 500  #maximum number of iterations
grid1 <- expand.grid(size=1:5,decay=nn_decay)
train1 <- train(label ~.,              #model design
                trainset,              #training set
                maxit=nn_maxit,        #maximum number of iterations
                tuneGrid=grid1,        #parameter space to explore
                metric="Accuracy",     #metric to evaluate models
                method="nnet",trace=FALSE)
nn_size <- train1$bestTune[1,1] #number of nodes in hidden layer

#Perform Neural Network fit
mod <- nnet(label ~.,       #model design
            trainset,       #training set
            size=nn_size,decay=nn_decay,maxit=nn_maxit)
print(mod)
summary(mod)
plotnet(mod,cex_val=0.7,circle_cex=3)
olden(mod,bar_plot=TRUE)
res <- predict(mod,          #model
               testset,      #testing set
               type="class") #returns classification/probabilities
               
#Confusion Matrix
conf_mat <- table(testset[,c("label")],res)
conf_mat

#Error rate
err_rt <- 100*(sum(conf_mat) - sum(diag(conf_mat)))/sum(conf_mat)
err_rt
