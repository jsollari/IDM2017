#autor:      Joao Sollari Lopes
#local:      INE, Lisboa
#criado:     27.06.2017
#modificado: 27.06.2017

## 1. FUNCTIONS
{
# Key:
# a - auxiliar
# d - dataframe
# f - file name
# i - iterator
# l - list
# n - vector size
# s - selected samples
# t - title
# v - vector
# w - weights

# Perform Partition around medoids clustering (wide format) [allow for specific number of partitions]
# dst1 - distance matrix of variables to be analysed
# kbest - specify number of partitions [if is.null(kbest) then calculate best number of partitions]
# col1 - colors for points
# border - add border to columns
# f1 - output filename
# main - title for plots
plot_pam_v2 = function(dst1,kbest=NULL,col1=NULL,border=1,f1,main=NULL){
    library("cluster")
    kmax <- dim(as.matrix(dst1))[1]-1
    asw <- rep(0,kmax)
    for(k in 2:kmax)
        asw[k] = pam(dst1,k)$silinfo$avg.width
    if(is.null(kbest)){
        kbest = which.max(asw)
    }
    pres1 = pam(dst1,kbest)
    d1 = cbind(1:kmax,asw)
    colnames(d1) = c("K","Average silhouette width")
    d2 = silhouette(pres1)
    f2 = paste(f1,"_1.csv",sep="")
    f3 = paste(f1,"_2.csv",sep="")
    f4 = paste(f1,"_3.csv",sep="")
    write.table(d1,f2,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",
      col.names=TRUE,row.names=FALSE)
    write.table(pres1$clusinfo,f3,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",
      col.names=TRUE,row.names=FALSE)
    write.table(d2,f4,quote=TRUE,sep=",",eol="\n",na="NA",dec=".",col.names=NA)
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7.0,height=3.5,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7.0,height=3.5,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfcol=c(1,2),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        plot(d1,type="h",main=main)
        points(kbest,asw[kbest],col="red",type= "h")
        plot(d2,main=main,sub="",xlab="Silhouette width",do.n.k=FALSE,
          border=border,do.clus.stat=FALSE,mgp=c(1.0,0.2,0.0),
          col=col1[rownames(d2)])
        par(oldpar)
        invisible(dev.off())
    }
}

# min-max transformation
# y - vector
# miny - minimum value to scale y
# maxy - maximum value to scale y
minmax = function(y,miny=NULL,maxy=NULL){
    if(is.null(miny)){miny = min(y)}
    if(is.null(maxy)){maxy = max(y)}
    return((y - miny)/(maxy - miny))
}

# Inverse min-max transformation
# y - vector
# miny - minimum value to scale y
# maxy - maximum value to scale y
minmax_inv = function(y,miny,maxy){
    return(y*(maxy - miny) + miny)
}

# Construct network (wide format) [allow for very many nodes]
# dst1 - distance matrix of variables to be analysed
# col1 - colors for points
# vmode - network layout
# edge_thr - threshold for edges to plot
# edge_lwd - scale for edges
# label_cex - size of node labels
# f1 - output filename
# main - title for plots
plot_sna_v2 = function(dst1,col1=NULL,vmode="fruchtermanreingold",edge_thr=0.80,
  edge_lwd=10,label_cex=0.8,f1,main=NULL){
    library("sna")
    m1 = as.matrix(dst1)
    m1 = -m1                            #invert distances
    mind = min(m1)
    maxd = max(m1)
    m1 = minmax(m1,mind,maxd)           #min-max transformation
    m1[m1 < edge_thr] = NA              #remove dist < thresh
    m1 = ((m1 - edge_thr))/(1 - edge_thr) #rescale dist
    edge_col = "lightblue"
    vertex_cex=1.5
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7,height=10,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7,height=10,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mar=c(1.1,1.1,2.1,1.1),mgp=c(0,0,0),xpd=TRUE,cex.main=0.9)
        gplot(m1,gmode="graph",displaylabels=TRUE,label.cex=label_cex,
          label.pos=5,mode=vmode,jitter=FALSE,main=paste(main,
          " [minmax(-dist) > ",edge_thr,"]",sep=""),vertex.cex=vertex_cex,
          vertex.col=col1[colnames(m1)],vertex.border=1,vertex.lty=1,
          vertices.last=TRUE,edge.lwd=edge_lwd,edge.lty=1,edge.lty.neg=1,
          edge.col=edge_col)
        par(oldpar)
        invisible(dev.off())
    }
}

# Summary statistics for network
# dst1 - distance matrix of variables to be analysed
# f1 - output filename
analyse_net = function(dst1,f1){
    library("sna")
    m1 = as.matrix(dst1)
    m1 = -m1                            #invert distances
    mind = min(m1)
    maxd = max(m1)
    m1 = (m1 - mind)/(maxd - mind)      #min-max transformation
    f2 = paste(f1,"_degree.csv",sep="")
    f3 = paste(f1,"_eigen.csv",sep="")
    f4 = paste(f1,"_dens.csv",sep="")
    #Degree - number of links incident upon a node
    net_deg = degree(m1,gmode="graph")
    s1 = order(net_deg,decreasing=TRUE)
    d1 = cbind(rownames(m1)[s1],net_deg[s1])
    colnames(d1) = c("Node","Degree")
    write.table(d1,f2,row.names=FALSE,col.names=TRUE,sep=",",dec=".",eol="\n",
      na="NA")
    #Eigenvector centrality - Influence of a node in a network
    net_evc = evcent(m1,gmode="graph",use.eigen=TRUE)
    s1 = order(abs(net_evc),decreasing=TRUE)
    d1 = cbind(rownames(m1)[s1],net_evc[s1])
    colnames(d1) = c("Node","Eigenvector")
    write.table(d1,f3,row.names=FALSE,col.names=TRUE,sep=",",dec=".",eol="\n",
      na="NA")
    #Density - the proportion of direct ties in a network relative to the total number possible
    nvert = nrow(m1)
    m2 = m1; diag(m2) = NA; m2[upper.tri(m2)] = NA
    nedge = sum(m2 != 0 & !is.na(m2))
    tedge = nties(m1,mode="graph")
    dens = gden(m1,mode="graph")
    d1 = data.frame(nvert,nedge,tedge,dens)
    write.table(d1,f4,row.names=FALSE,col.names=TRUE,sep=",",dec=".",eol="\n",
      na="NA")
}

# Reduce number of predictors
# dat1 - data frame of variables to be analysed (only predictors)
# y - response variable
# thr1 - upper threshold for correlation between predictors 
# thr2 - lower threshold for correlation between predictors and response variables [if thr2==NULL then reduces till maxsize]
# thr3 - lower threshold for coefficient of variation [if thr3==NULL then reduces till maxsize]
# thr4 - upper threshold for Variance Inflation Factor (VIF)
# maxsize - maximum number of predictors [if maxsize==NULL then reduces till nsubj - 1]
# f1 - output filename
reduce_predictors_v2 = function(dat1,y,thr1=0.95,method="spearman",thr2=NULL,
  thr3=NULL,thr4=10,maxsize=NULL,f1){
# Note: Sqrt(R^2) for continuous variables is the correlation between variables. Sqrt(R^2)
# for continuous and categorical variables is the correlation between observed and predicted
# values.
    if(is.null(maxsize)){
        maxsize = nrow(dat1) - 1
    }
    if(is.null(thr2)){
        doStep2 = TRUE
    }else if(thr2 > 0){
        doStep2 = TRUE
    }else{
        doStep2 = FALSE
    }
    if(is.null(thr3)){
        doStep3 = TRUE
    }else if(thr3 > 0){
        doStep3 = TRUE
    }else{
        doStep3 = FALSE
    }
    write("Step1: Highly correlated between themselves",file=f1,append=FALSE)
    if(thr1 < 1){                       #remove variables highly correlated between themselves
        vars = colnames(dat1)
        rm_vars = c()
        for(i1 in vars){
            for(i2 in vars){
                if(i1 < i2 & !i1 %in% rm_vars){
                    rres = cor(dat1[,i1],dat1[,i2],method=method,use="complete")
                    if(abs(rres) > thr1){
                        write(paste("  Remove ",i2," (r_",i1," = ",round(rres,2),
                          ")",sep=""),file=f1,append=TRUE)
                        rm_vars = c(rm_vars,i2)
                    }
                }
            }
        }
        dat1 = dat1[,!vars %in% rm_vars]
    }else{
        write("<Skip step>",file=f1,append=TRUE)
    }
    write("Step2: Low correlated with response variable",file=f1,append=TRUE)
    if(doStep2){                        #remove variables with low correlation with dependent variable
        vars = colnames(dat1)
        if(is.factor(y)){
            y = as.character(y)
        }
        cors = apply(dat1,2,function(x){sqrt(summary(lm(y~x))$r.squared)})
        if(is.null(thr2)){
            n_rm = ncol(dat1) - maxsize #n_rm, so that, nvars = nsubj - 1
            if(n_rm <= 0){
                thr2 = 0.00
            }else{
                thr2 = sort(cors)[n_rm]
            }
        }
        rm_vars = cors <= thr2
        if(sum(rm_vars) == 0){
            write("<None removed>",file=f1,append=TRUE)
        }else{
            dat1 = dat1[,!rm_vars]
            tmp = cbind(vars[rm_vars],cors[rm_vars])
            tmp = apply(tmp,1,function(x){paste("  Remove ",x[1]," (r = ",
              round(as.numeric(x[2]),2),")",sep="")})
            write(tmp,file=f1,append=TRUE)
        }
    }else{
        write("<Skip step>",file=f1,append=TRUE)
    }
    write("Step3: Variables with low variation",file=f1,append=TRUE)
    if(doStep3){                        #remove variables with low variation
        vars = colnames(dat1)
        sds = apply(dat1,2,sd,na.rm=TRUE)
        means = apply(dat1,2,mean,na.rm=TRUE)
        if(is.null(thr3)){
            n_rm = ncol(dat1) - maxsize #n_rm, so that, nvars = nsubj - 1
            if(n_rm <= 0){
                thr3 = 0.00
            }else{
                thr3 = sort(sds/means)[n_rm]
            }
        }
        rm_vars = sds/means <= thr3
        if(sum(rm_vars) == 0){
            write("<None removed>",file=f1,append=TRUE)
        }else{
            dat1 = dat1[,!rm_vars]
            tmp = cbind(vars[rm_vars],sds/means[rm_vars])
            tmp = apply(tmp,1,function(x){paste("  Remove ",x[1]," (r = ",
              round(as.numeric(x[2]),2),")",sep="")})
            write(tmp,file=f1,append=TRUE)
        }
    }else{
        write("<Skip step>",file=f1,append=TRUE)
    }
    write("Step4: Variables with low VIF",file=f1,append=TRUE)
    if(thr4 < Inf){                     #remove variables using VIF (reduce multicollinearity)
        vars = colnames(dat1)
        keep_vars = stepwise_vif(dat1,thresh=thr4,trace=FALSE)
        if(sum(!keep_vars) == 0){
            write("<None removed>",file=f1,append=TRUE)
        }else{
            dat1 = dat1[,keep_vars]
            tmp = paste("  Remove",vars[!vars %in% keep_vars])
            write(tmp,file=f1,append=TRUE)
        }
    }else{
        write("<Skip step>",file=f1,append=TRUE)
    }
    return(dat1)
}
 
# Stepwise Variance Inflation Factors (VIF) selection to reduce collinearity
# [From "Marcus W Beck" in https://gist.github.com/fawda123/4717702#file-vif_fun-r]
# in_frame - data frame of variables to be analyzed
# thresh - threshold for VIF
# trace - print output of each iteration
stepwise_vif = function(in_frame,thresh=10,trace=TRUE){
    if(class(in_frame) != 'data.frame'){
        in_frame = data.frame(in_frame)
    }
    #get initial vif value for all comparisons of variables
    var_names = names(in_frame)
    vif_init = vif(in_frame)
    vif_max = max(as.numeric(vif_init),na.rm=TRUE)
    if(vif_max < thresh){
        if(trace==TRUE){
            print(data.frame(vif_init))
            cat('\n')
            cat(paste('All variables have VIF < ', thresh,', max VIF ',
              round(vif_max,2), sep=''),'\n\n')
        }
        return(var_names)
    }else{
        in_dat = in_frame
        #backwards selection of explanatory variables, stops when all VIF values are below 'thresh'
        while(vif_max >= thresh){
            var_names = names(in_dat)
            vif_vals = vif(in_dat)
            imax = which(vif_vals == max(as.numeric(vif_vals),na.rm=TRUE))[1]
            vif_max = as.numeric(vif_vals[imax])
            if(vif_max < thresh){
                break
            }
            if(trace==TRUE){
                print(data.frame(vif_vals))
                cat('\n')
                cat('removed: ',names(vif_vals)[imax],vif_max,'\n\n')
                flush.console()
            }
            in_dat = in_dat[,!names(in_dat) %in% names(vif_vals)[imax]]
        }
        return(names(in_dat))
    }
}

# Calculate Variance Inflation Factors (VIF) at a model-level
# [based on HH::vif()]
# dat1 - data frame of variables to be analysed (only predictors)
vif = function(dat1){
#see https://en.wikipedia.org/wiki/Variance_inflation_factor
#see car::vif() for calculation of VIF at variable-level
#see HH::vif() for calculation of VIF at model-level
#see faraway::vif() for calculation of VIF at model-level
#see fmsb::vif() for calculation of VIF at model-level
# Note: VIF can be defined at variable-level (using the partial R^2) and at the model
# level (using R^2, i.e. coefficient of determination).
    vars = colnames(dat1)
    nvars = length(vars)
    r2 = vector("numeric",length=nvars)
    names(r2) = vars
    for(i1 in 1:nvars){
        tmp = lm(dat1[,i1] ~ data.matrix(dat1[,-i1]),na.action="na.omit")
        r2[i1] = 1/(1 - summary(tmp)$r.squared)
    }
    return(r2)
}

# Automatic model selection for multiple linear regression using General LM [simplified version]
# dat1 - data frame of independet variables (predictors)
# resp - dependent variable (response variable)
# maxs1 -  maximum number of terms in candidate models {"hard"}
# crit - information criteria to use {"aic","aicc","bic","qaic","qaicc"}
# thrs  - threshold for nmods to perform exhaustive search
# popsize - Population size
# mutrate - Per locus (i.e. per term) mutation rate, [0,1]
# sexrate - Sexual reproduction rate, [0,1]
# imm - Immigration rate, [0,1]
# deltaM - Stop Rule: change in mean IC
# deltaB - Stop Rule: change in best IC
# conseq - Stop Rule: times with no improvement
# f1 - output filename
mod_select_lm_v2 = function(dat1,resp,maxs1="hard",crit="aicc",thrs=100000,
  popsize=100,mutrate=0.001,sexrate=0.1,imm=0.3,deltaM=0.05,deltaB=0.05,
  conseq=5,nreps=2,f1){
    library("glmulti")
    d1 = data.frame(resp,dat1)
    colnames(d1) = sapply(colnames(d1),function(x){gsub("[-|.]","_",x)}) #remove symbols used in "formulas"
    nvars = ncol(d1)
    vars = colnames(d1)
    hmax = (nrow(d1) - 1) - 3           #hard max (i.e. nparam = nsubj - 3)
    maxs1 = min(maxs1,hmax)             #prevent maxs1 > hard max
    nmods = sum(sapply(1:maxs1,function(x){choose(nvars-1,x)})) + 1
    if(nmods < thrs){                   #exhaustive search
        res2 = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,maxsize=maxs1,
          method="h",crit=crit,plotty=FALSE,report=FALSE,includeobjects=TRUE,
          marginality=FALSE,model=TRUE,fitfunction="lm")
    }else{                              #genetic algorithm
        res1 = vector("list",length=nreps)
        for(i1 in 1:nreps){
            f2 = paste(f1,"_rep",i1,"_1",sep="")
            res1[[i1]] = glmulti(y=vars[1],xr=vars[2:nvars],data=d1,level=1,
              maxsize=maxs1,method="g",crit=crit,confsetsize=100,plotty=FALSE,
              report=FALSE,name=f2,includeobjects=TRUE,marginality=FALSE,
              model=TRUE,popsize=popsize,mutrate=mutrate,sexrate=sexrate,
              imm=imm,deltaM=deltaM,deltaB=deltaB,conseq=conseq,
              fitfunction="lm")
        }
        if(nreps > 1){
            res2 = glmulti::consensus(xs=res1,confsetsize=100)
        }else{
            res2 = res1
        }
    }
    t1 = weightable(res2)
    t2 = t1[t1[,crit] <= min(t1[,crit]) + 2,]
    t3 = coef(res2)
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,"_1.csv",sep="")
    f4 = paste(f1,"_2.csv",sep="")
    capture.output(print(res2),file=f2,append=FALSE,type="output")
    write.table(t2,file=f3,sep=",",row.names=TRUE,col.names=NA)
    write.table(t3,file=f4,sep=",",row.names=TRUE,col.names=NA)
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7.0,height=10,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7.0,height=10,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfcol=c(3,1),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          oma=c(0,0,0,0),tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        par(mar=c(2.2,10.2,1.2,1.2))
        plot(res2,type="s",cex.names=0.6)
        par(mar=c(2.2,2.2,1.2,1.2))
        plot(res2,type="p")
        plot(res2,type="w")
        par(oldpar)
        invisible(dev.off())
    }
    a1 = formula(res2@objects[[1]]$terms) #formula
    a2 = res2@objects[[1]]$model        #data
    l1 = list(a1,a2)
    names(l1) = c("formula","data")
    return(l1)
}

# Fit multiple linear regression using General LS
# form1 - formula of the model to fit
# dat1 - data frame of variables to fit 
# f1 - output filename
# main - title for plots
fit_lm = function(form1,dat1,f1,main){
#see https://en.wikipedia.org/wiki/General_linear_model
# Note: Multiple linear regression assumptions:
# 1. The error terms need to be independent.
# 2. The error variance is homoskedastic.
# 3. The errors are normally distributed with mean = 0
# 4. The independent variables should be independent from each other (i.e. little or no multicollinearity).
# 5. Linearity between independent variables and dependent variable.
# 6. Large sample sizes (i.e. 10-30 cases per parameter).
    library("car")
    fit1 = lm(form1,data=dat1)
    sfit1 = summary(fit1)
    t1 = Anova(fit1,type="II")
    f2 = paste(f1,".log",sep="")
    f3 = paste(f1,".csv",sep="")
    capture.output(print(sfit1),file=f2,append=FALSE,type="output")
    write.table(t1,file=f3,sep=",",row.names=TRUE,col.names=NA)
    vars = all.vars(form1)
    y = dat1[,vars[1]]
    f = fit1$fitted.values
    r = fit1$residuals
    for(i1 in c("eps","tiff")){
        if(i1 == "eps"){
            postscript(file=paste(f1,".eps",sep=""),width=7.0,height=7.0,
              colormodel="rgb",horizontal=FALSE,onefile=FALSE,paper="special",
              family="Helvetica")
        }else{
            tiff(file=paste(f1,".tif",sep=""),units="in",width=7.0,height=7.0,
              res=300,compression="lzw",family="sans")
        }
        oldpar = par(mfcol=c(2,2),mar=c(2.2,2.2,1.2,1.2),mgp=c(1.0,0.2,0.0),
          tcl=-0.2,cex=0.8,cex.lab=0.8,cex.axis=0.8,cex.main=0.9)
        qqnorm(r,col=4); qqline(r)                        #wg-err_norm: qqplot 
        plot(f,sqrt(abs(r)),xlab="Fitted",ylab="sqrt(abs(Residuals))",col=4)
        abline(h=0,lty=2)                                 #wg-err_var: scatterplot
        hist(r,freq=FALSE,main="",col=4,xlab="Residuals") #wg-err_norm: histogram
        plot(f,y,xlab="Fitted values",ylab="Values",main=main,col=4)
        abline(c(0,1),lty=2)                              #goodness-of-fit: scatterplot
        par(oldpar)
        invisible(dev.off())
    }
}

}

