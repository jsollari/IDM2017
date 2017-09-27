#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     27.06.2017
#modificado: 28.06.2017

#1. FUNCTIONS
source("misc_v2.2.r")

#2. READ DATA
f1 = "../data/lmktattract.csv"
f2 = "../data/lmktmobil.csv"
x_all = read.table(file=f1,header=TRUE,sep=",",dec=".",row.names=1)
y_all = read.table(file=f2,header=TRUE,sep=",",dec=".",row.names=1)
dim(x_all)
dim(y_all)
#str(x_all)
#str(y_all)
#summary(x_all)
#summary(y_all)

# 3. ANALYSE DATA

#remove attributes with missing data
ina = apply(x_all,2,function(x){any(is.na(x))})
x = x_all[,!ina]
dim(x)

#remove data entries with missing data
y = y_all[,"lmktm_Total",drop=FALSE]
ina = is.na(y[,1])
x = x[!ina,]
y = y[!ina,,drop=FALSE]
dim(x)
dim(y)

#reduce number of attributes
f1 = "../results/lmktattract_2/datred.log"
x = reduce_predictors_v2(x,y[,1],
  thr1=0.90,  #upper threshold for correlation between predictors 
  method="pearson",
  thr2=NULL,  #lower threshold for correlation between predictors and response variables [if thr2==NULL then reduces till maxsize]
  thr3=0.00,  #lower threshold for coefficient of variation [if thr3==NULL then reduces till maxsize]
  thr4=Inf,   #upper threshold for Variance Inflation Factor (VIF)
  maxsize=30, #maximum number of predictors [if maxsize==NULL then reduces till nsubj - 1]
  f1=f1)
dim(x)

#save final data
f1 = "../results/lmktattract_2/datsel.csv"
write.table(t(apply(x,2,summary)),file=f1,sep=",",dec=".",
  row.names=TRUE,col.names=NA)

#perform GA
f1 = "../results/lmktattract_2/modsel"
f2 = paste(f1,".RDS",sep="")
f3 = paste(f1,".log",sep="")
f4 = paste(f1,"_1.csv",sep="")
f5 = paste(f1,"_2.csv",sep="")
if(!file.exists(f3) | !file.exists(f4) | !file.exists(f5)){
    mod1 = mod_select_lm_v2(x,y,
      maxs1=10,        #Maximum number of attributes
      popsize = 100,   #Population size
      mutrate = 0.001, #Per locus (i.e. per term) mutation rate, [0,1]
      sexrate = 0.1,   #Sexual reproduction rate, [0,1]
      imm = 0.0,       #Immigration rate, [0,1]
      deltaM=1e-6,     #Stop Rule: change in mean IC
      deltaB=1e-6,     #Stop Rule: change in best IC
      conseq = 5,      #Stop Rule: times with no improvement
      nreps = 4,       #Number of repeats
      f1=f1)
    save(mod1,file=f2)
}else{
    load(file=f2)
}

#fit best model
f1 = "../results/lmktattract_2/fit"
f2 = paste(f1,".log",sep="")
f3 = paste(f1,".csv",sep="")
if(!file.exists(f2) | !file.exists(f3)){
    fit_lm(mod1$formula,mod1$data,f1=f1,main="LMkt Mobility")
}
