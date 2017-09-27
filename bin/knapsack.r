#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     02.07.2017
#modificado: 02.07.2017

library("GA")

#Define data 
n <- c("pocketknife","beans","potatoes","unions","sleepingbag","rope","compass")                # Items 
p <- c(10,20,15,2,30,10,30) # Profits 
w <- c(1,5,10,1,7,5,1)      # Weights
W <- 20                     # Knapsack ’s capacity 

#Define fitness function 
knapsack <- function(x){ 
    if(sum(x*w) > W)
	    return(0)
    else
		return(sum(x*p))
}

# Run SGA
SGA <- ga(type="binary", 
  fitness=knapsack,         #fitness function
  nBits=length(n),          #chromosome length
  popSize=100,              #population size
  pcrossover=0.8,           #crossover rate
  pmutation=0.1,            #mutation rate
  elitism=5,                #number of best individuals sure to be selected
  maxiter=100,              #number of generations
  names=n,                  #name of "genes"
  seed=101)

res <- SGA@solution
sum(res)                    #Total number of selected items
sum(res*p)                  #Total profit of selected items
sum(res*w)                  #Total weight of selected items