#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     28.06.2017
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
main1 = "Labour market attractiveness"  #title for plots
col1 = c(rbind(hsv(h=seq(0,1,len=8),s=1/3,v=5/5)[1:7],  #very light
               hsv(h=seq(0,1,len=8),s=2/3,v=5/5)[1:7],  #light
               hsv(h=seq(0,1,len=8),s=3/3,v=5/5)[1:7],  #normal
               hsv(h=seq(0,1,len=8),s=3/3,v=4/5)[1:7])) #dark
names(col1) = rownames(x_all)           #colours for EU28

#normalize and scale data
w_all = c(1,rep(1/5,5),1,1,1,1,1,rep(1/6,6),rep(1/2,2),rep(1/3,3),rep(1/20,20),
  rep(1/2,2),1,1,1,1,rep(1/27,27),1)
x_norm = t(apply(scale(x_all),1,function(x){x*w_all}))

#calculate distances
#from help(dist): "Missing values are allowed, and are excluded from all 
#  computations involving the rows within which they occur. If some columns are
#  excluded in calculating a [...] distance, the sum is scaled up proportionally
#  to the number of columns used."
x_dist = dist(x_norm,method="euclidean")

#perform Social Network Analysis
f1 = "../results/lmktattract_1/net"
plot_sna_v2(x_dist,col1=col1,edge_thr=0.65,f1=f1,main=main1)
analyse_net(x_dist,f1=f1)

#perform Partition Around Medoids analysis
f1 = "../results/lmktattract_1/pam"
plot_pam_v2(x_dist,col1=col1,f1=f1,main=main1)
kbest = 10
f1 = paste("../results/lmktattract_1/pam_grps",kbest,sep="")
plot_pam_v2(x_dist,kbest=kbest,col1=col1,f1=f1,main=main1)
