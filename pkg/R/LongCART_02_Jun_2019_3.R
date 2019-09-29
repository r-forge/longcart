#-------------------------------------------------------------------------
#  LongCART algorithm to construct regression tree with Longitudinal data
#  Date: 28-Sep-2016
#  Developed by: Madan G Kundu (Email: madan_g.kundu@novartis.com)
#------------------------------------------------------------------------- 
##### Following checks need to be done
#1. Check: No missing value for patid and response variable
#2. Check: Grouping varaible should have only observations for each patid
#3. Check: Length of gvars and tgvars should be equal
#4. Check: The variables should exist
############## End of Checks


#--------------------------------------------------------#
#   Stability test for categorical grouping variable     #
#--------------------------------------------------------#

StabCat<- function(data, patid, fixed, splitvar)
{
  #------------Checks
  if(!exists(as.character(substitute(data)), envir=sys.frame(-1L))) stop("Dataset does not exist\n")
  if(!is.data.frame(data)) stop("Dataset does not exist\n")
  
  #-------- Checks for patid variable
  if(!patid %in% colnames(data))  stop("The column ", patid, " containing subjects id is missing in dataset.\n")
  data$patid<- data[[patid]]
  #-------------Drop observations with missing patid
  data<- data[!is.na(data[[patid]]),]
  
  #-------- Checks for Y variable
  Y.name<- as.character(attr(as.Formula(fixed), "lhs"))
  if(!Y.name %in% colnames(data))  
    stop("The column ", Y.name, " containing subjects id is missing in dataset.\n")
  #-------------Drop observations with missing Y and patid
  data<- data[!is.na(data[[Y.name]]),]
  
  #--------- Check whether splitvar exist or not
  if(!splitvar %in% colnames(data))  
    stop("The column ", splitvar, " is missing in dataset.\n")
  

sig.stab<- NA
message("Stability Test for Categorical grouping variable \n", appendLF=F)
data.t<- data[!is.na(splitvar),]
splitvar<- data[[splitvar]]
splitvar<- splitvar[!is.na(splitvar)]

temp<- gsub(pattern=" ", replacement="", fixed, fixed=T)[3]
noint=grepl(pattern="-1", temp, fixed=T)

form<-Formula(fixed)
Y.name<- as.character(attr(form, "lhs"))
X.name1<- as.character(attr(form, "rhs"))
X.name1<- gsub(pattern=" ", replacement="", x=X.name1)
X.name1<- gsub(pattern="-1", replacement="", x=X.name1)
X.name<- unlist(strsplit(X.name1, "+", fixed = TRUE))

Y<- data.t[[Y.name]]
p<- length(X.name)

if(X.name1=='1') X<- rep(1, length(Y))
if(X.name1!='1'){
  for(p.i in 1:p){
		if(p.i==1) X<-data.t[[X.name[p.i]]]
		if(p.i>1)  X<-cbind(X, data.t[[X.name[p.i]]])
	}
}

id<- data.t$patid
Group<- tapply(splitvar, id, unique)
Group<- Group[!is.na(Group)]
G<- length(unique(Group))

## Fit mixed model
fit<- NULL
if (X.name1=='1') {
  p<- 1
  fit <- try(fit <- lme(fixed=Y~1, random=~1|id) , silent=T)
} else if (noint) {
  p<- p
  fit <- try(fit <- lme(fixed=Y~X-1, random=~1|id) , silent=T)
} else if (!noint) {
  p<- p+1
  fit <- try(fit <- lme(fixed=Y~X, random=~1|id) , silent=T)
}

if(!is.null(fit)){
#Get sd_int and sigma.e
sd_int<- as.numeric(VarCorr(fit)[1,2])
sd_err<- as.numeric(VarCorr(fit)[2,2])

#Get beta.hat
beta <- matrix(fit$coeff$fixed,ncol=1)

#Get J
idvec<- unique(id)
N<- length(idvec)
J=solve(vcov(fit))/N

u<- matrix(0, N, p)

for(i in 1:N)	{
	sub=idvec[i]
	Y.i<- subset(Y, id==sub)
	X.i<- subset(X, id==sub)
	if (noint|X.name1=='1') X.i<- model.matrix(~X.i-1)
	if (!noint) X.i<- model.matrix(~X.i)
	n.i<- length(Y.i)
	if(p==1) X.i<- matrix(X.i, nrow=n.i, ncol=1)
	if(n.i==1) V.i<- matrix(sd_int^2+sd_err^2, nrow=n.i, ncol=n.i) 
		else V.i<- matrix(sd_int^2, nrow=n.i, ncol=n.i) + diag(rep(sd_err^2, n.i))
	u[i,]<- t(t(X.i)%*%solve(V.i)%*%(Y.i-X.i%*%beta))
	}

u.grp<- W<- W.thetahat<- matrix(NA, nrow=G, ncol=p)
n.grp<- table(splitvar)

for(pi in 1:p) u.grp[,pi]<- tapply(u[,pi], Group, sum)

#Test statistic
Test.stat=0
for(Gi in 1:G){
	txx<- matrix(u.grp[Gi,], ncol=1)
	Test.stat<- Test.stat+t(txx)%*%solve(n.grp[Gi]*J)%*%txx
	}

#p value
sig.stab<- 1- pchisq(Test.stat, (G-1)*p)

message('Test.statistic=', c(Test.stat), ', p-value=', c(sig.stab), ' \n', appendLF=F)
}
ret<- list(pval=c(sig.stab))
ret
}

#--------------------------------------------------------#
#   Stability test for continuous grouping variable      #
#--------------------------------------------------------#


StabCont<- function(data, patid, fixed, splitvar)
{
CDF.D <- function(x)
	{
	k <- seq(1:20)
	1 + 2 * sum((-1)^k * exp(- 2 * k^2 * x^2))
	}

#------------Checks
if(!exists(as.character(substitute(data)), envir=sys.frame(-1L))) stop("Dataset does not exist\n")
if(!is.data.frame(data)) stop("Dataset does not exist\n")

#-------- Checks for patid variable
if(!patid %in% colnames(data))  stop("The column ", patid, " containing subjects id is missing in dataset.\n")
data$patid<- data[[patid]]
#-------------Drop observations with missing patid
data<- data[!is.na(data[[patid]]),]

#-------- Checks for Y variable
Y.name<- as.character(attr(as.Formula(fixed), "lhs"))
if(!Y.name %in% colnames(data))  
  stop("The column ", Y.name, " containing subjects id is missing in dataset.\n")
#-------------Drop observations with missing Y and patid
data<- data[!is.na(data[[Y.name]]),]

#--------- Check whether splitvar exist or not
if(!splitvar %in% colnames(data))  
  stop("The column ", splitvar, " is missing in dataset.\n")

sig.stab<- NA
message("Stability Test for Continuous grouping variable \n", appendLF=F)
data.t<- data[!is.na(splitvar),]
splitvar<- data[[splitvar]]
splitvar<- splitvar[!is.na(splitvar)]

temp<- gsub(pattern=" ", replacement="", fixed, fixed=T)[3]
noint=grepl(pattern="-1", temp, fixed=T)

form<-Formula(fixed)
Y.name<- as.character(attr(form, "lhs"))
X.name1<- as.character(attr(form, "rhs"))
X.name1<- gsub(pattern=" ", replacement="", x=X.name1)
X.name1<- gsub(pattern="-1", replacement="", x=X.name1)
X.name<- unlist(strsplit(X.name1, "+", fixed = TRUE))

Y<- data.t[[Y.name]]
p<- length(X.name)

if(X.name1=='1') X<- rep(1, length(Y))
if(X.name1!='1'){
  for(p.i in 1:p){
		if(p.i==1) X<-data.t[[X.name[p.i]]]
		if(p.i>1)  X<-cbind(X, data.t[[X.name[p.i]]])
	}
}

id<- data.t$patid
Group<- tapply(splitvar, id, unique)
Group<- Group[!is.na(Group)]
G<- length(unique(Group))

## Fit mixed model
fit<- NULL
if (X.name1=='1') {
  p<- 1
  fit <- try(fit <- lme(fixed=Y~1, random=~1|id) , silent=T)
} else if (noint) {
  p<- p
  fit <- try(fit <- lme(fixed=Y~X-1, random=~1|id) , silent=T)
} else if (!noint) {
  p<- p+1
  fit <- try(fit <- lme(fixed=Y~X, random=~1|id) , silent=T)
}

if(!is.null(fit)){
#Get sd_int and sigma.e
sd_int<- as.numeric(VarCorr(fit)[1,2])
sd_err<- as.numeric(VarCorr(fit)[2,2])

#Get beta.hat
beta <- matrix(fit$coeff$fixed,ncol=1)

#Get J
idvec<- unique(id)
N<- length(idvec)
J=solve(vcov(fit))/N

if(p==1){
	J.sqrt<- sqrt(J)
	} else {
	J.eig <- eigen(J)
	J.sqrt <- J.eig$vectors %*% diag(sqrt(J.eig$values)) %*% solve(J.eig$vectors)
	}

u<- matrix(0, N, p)

for(i in 1:N)	{
	sub=idvec[i]
	Y.i<- subset(Y, id==sub)
	n.i<- length(Y.i)
	X.i<- subset(X, id==sub)
	if (noint|X.name1=='1') X.i<- model.matrix(~X.i-1)
  if (!noint) X.i<- model.matrix(~X.i)
	if(p==1) X.i<- matrix(X.i, nrow=n.i, ncol=1)
	if(n.i==1) V.i<- matrix(sd_int^2+sd_err^2, nrow=n.i, ncol=n.i) 
		else V.i<- matrix(sd_int^2, nrow=n.i, ncol=n.i) + diag(rep(sd_err^2, n.i))
	u[i,]<- t(t(X.i)%*%solve(V.i)%*%(Y.i-X.i%*%beta))
	}

u.grp<- W<- W.thetahat<- matrix(NA, nrow=G, ncol=p)
n.grp<- table(splitvar)
t.grp<- cumsum(n.grp)/N

for(pi in 1:p)	{
	u.grp[,pi]<- tapply(u[,pi], Group, sum)
	W.thetahat[,pi]<- cumsum(u.grp[,pi])/sqrt(N)
	}

M.thetahat<- solve(J.sqrt)%*%t(W.thetahat)

#Test statistic & p-value
Test.stat<- apply(abs(M.thetahat),1,max)
out.p<- 1- unlist(lapply(Test.stat, FUN=CDF.D))
out.p.adj<- p.adjust(out.p, method = "hochberg")
message('Test.statistic=', Test.stat, ', Adj. p-value=', out.p.adj, ' \n', appendLF=F)
sig.stab<- min(out.p.adj)
}
ret<- list(pval=c(sig.stab))
ret
}



#--------------------------------------------------------#
#   Final function						   #
#--------------------------------------------------------#
LongCART<- function(data, patid, fixed, gvars, tgvars, minsplit=40, minbucket=20, alpha=0.05, coef.digits=2, print.lme=F)
	{

#--------------------------------------------------------#
#   Obtaining log-likelihood from a given model		   #
#--------------------------------------------------------#
llike<- function(data, fixed)
{
  if(exists('txx.out')) rm(txx.out)
  loglik<- NA
  try(txx.out<- lme(fixed, data=data, method="ML", random=~1|patid,
                na.action=na.omit), silent=T)
  if(exists('txx.out')) loglik<- logLik(txx.out)[1]
  loglik
}

#--------------------------------------------------------#
#   LRT at a given split					   #
#--------------------------------------------------------#
test4split<- function(cutoff, data, fixed, splitvar)
{
  LRTS<-NA
  message(cutoff, appendLF=F)
  left.data<- data[data[[splitvar]]<cutoff, ]     #Updated on 02-Jun-2019
  right.data<- data[data[[splitvar]]>=cutoff, ]  #Updated on 02-Jun-2019
  
  rn.ll<- llike(data=data, fixed=fixed)
  left.ll<- llike(data=left.data, fixed=fixed)
  right.ll<- llike(data=right.data, fixed=fixed)
  
  LRTS<- 2*(left.ll+right.ll-rn.ll)
  message(" ", appendLF=F)
  LRTS
}

#--------------------------------------------------------#
#   Choosing best split						   #
#--------------------------------------------------------#
single.group<- function(data, fixed, splitvar, minbucket)
{
  xcut<- NA
  improve<- NA
  
  Y.name<- as.character(attr(as.Formula(fixed), "lhs"))
  data<- data[!is.na(data[[Y.name]]),]
  data<- data[!is.na(data[[splitvar]]),] 
  
  #Unique patid and splitvar combination
  temp<- unique(data[,c("patid", splitvar)])
  N<- tapply(temp$patid, temp[[splitvar]], length)
  temp1<- as.data.frame(list(cuts=as.numeric(names(N)), N=N))
  temp1$left<- cumsum(temp1$N) - temp1$N
  temp1$right<- sum(temp1$N) - temp1$left
  temp2<- temp1[temp1$left>=minbucket & temp1$right>=minbucket,]
  vals<- temp2$cuts
  
  message("Evaluations of cutoffs for maximum improvements (Maximum cutoff value = ", max(vals), ")\n", appendLF=F)
  LRT.vec<- lapply(vals, FUN=test4split, data=data, fixed=fixed, 
                   splitvar=splitvar)
  message(".\n", appendLF=F)
  LRT.vec<- unlist(LRT.vec)
  
  if(any(!is.na(LRT.vec))){
    index<-which.max(LRT.vec) 
    xcut=vals[index]
    improve=LRT.vec[index]
  } 
  ret<- list(xcut=xcut, improve=improve)
  ret
}


#--------------------------------------------------------#
#   Choosing best partitioning variable			   #
#--------------------------------------------------------#
bestsplit<- function (data, fixed, gvars, tgvars, node.name, minbucket, alpha)
	{
	ngvars<- length(gvars)
	best.gvar<- NA
	best.cutoff<- NA
	improve<- NA
	min.pval.adj<- NA
	stab.pval<- numeric(length=ngvars)
	for (v in 1:ngvars)
		{
		stab.pval[v]<- NA
		splitvar<-gvars[v]
		message("\nSplitting variable: ", splitvar, '\n', appendLF=F)
		data1<- data[!is.na(splitvar),]
		G<- length(unique(data1[[splitvar]]))
		message('G=', G, '\n', appendLF=F)
		if (G>1){
			if(tgvars[v]==0) stab.pval[v]<- StabCat(data=data1, patid="patid", fixed=fixed, splitvar=splitvar)$pval
			if(tgvars[v]==1) stab.pval[v]<- StabCont(data=data1, patid="patid", fixed=fixed, splitvar=splitvar)$pval
			}
		else message("Instability test was NOT performed. \n", appendLF=F)	
		}
	message('\n stab.pval=', stab.pval, '\n', appendLF=F)
	if(any(!is.na(stab.pval))){
		stab.pval.adj<- p.adjust(stab.pval, method = "hochberg")
		sel.v<- which.min(stab.pval)
		min.pval.adj<- stab.pval.adj[sel.v]
		best.gvar<- gvars[sel.v]
		message('\n stab.pval.adj=', stab.pval.adj, '\n', appendLF=F)
		message('\n alpha=', alpha, '\n', appendLF=F)
		if (min.pval.adj<alpha){
		  bestcut<- single.group(data, fixed, best.gvar, minbucket)
		  best.cutoff<- bestcut$xcut
		  improve<- bestcut$improve
		  }
		}
	return(list(node=node.name, gvar=best.gvar, cutoff=best.cutoff, improve=improve, pval=min.pval.adj))
	}

#--------------------------------------------------------#
#   Coefficient Estimate							   #
#--------------------------------------------------------#
coeff.Estimate<- function(LMEobj, coef.digits=2)
	{
  ct<- LMEobj$tTable
  p<- nrow(ct)
  ct.names<- row.names(ct)
  yval<- " "
  for (ct.i in 1:p) yval<- paste0(yval, '+', round(ct[ct.i,1],coef.digits),ct.names[ct.i])
  yval<- gsub(pattern="(Intercept)", replacement="", x=yval, fixed=T)
  yval<- gsub(pattern="+-", replacement="-", yval, fixed=T)
  yval<- gsub(pattern=" +", replacement="", yval, fixed=T)
  yval<- gsub(pattern=" ", replacement="", yval, fixed=T)
  
  yval
	}


#--------------------------------------------------------#
#   Iterative splitting						   #
#--------------------------------------------------------#

rsplit<- function(data, fixed, gvars, tgvars, id, split, alpha, minsplit, minbucket, Rate, loglik, env=parent.frame(), coef.digits=2)
	{
	s.var<- unlist(split[2])
	s.cut<- unlist(split[3])
  s.improve<- unlist(split[4])
	s.pval<- unlist(split[5])

	N<- length(unique(data$patid))	

	if (id==1) {
		env$idlist<- list(unique(data[['patid']]))
		env$Treeout<- data.frame(id=id, N=N, yval=Rate, splitvar=s.var, 
                          cutoff=s.cut, pstab=s.pval, loglik=loglik, improve=s.improve, stringsAsFactors = FALSE)
		}
	else {
		env$idlist<- c(env$idlist, list(unique(data[['patid']])))
		env$Treeout<- rbind(env$Treeout, c(id, N, Rate, s.var, s.cut, s.pval, loglik, s.improve))
		}

	if (!is.na(s.var) && !is.na(s.cut))
		{
		data<- data[!is.na(data[[s.var]]),]
		##Left Node
		id_l=id*2
		message("---------------------------------------- \n", appendLF=F)
		message("NODE ", id_l, "- Rule:", s.var, " <", s.cut, '\n', appendLF=F) #Updated on 02-Jun-2019
		message("---------------------------------------- \n", appendLF=F)
		data_l<- subset(data, data[[s.var]]<s.cut)              #Updated on 02-Jun-2019
		left.subj<- length(unique(data_l$patid))	
		message("No. of individual in left node: ", left.subj, ' \n', appendLF=F)
		
	  	tjj<- NULL
		tjj<- summary(lme(fixed, data=data_l, method="REML", random=~1|patid))
                if(print.lme) cat("NODE ", id_l, "- Rule:", s.var, " <", s.cut, "\n")
		if(print.lme) print(tjj$tTable)
		Rate<- coeff.Estimate(tjj, coef.digits=coef.digits)
		loglik<- logLik(lme(fixed, data=data_l, method="REML", random=~1|patid))

		if (left.subj>=minsplit)
			{
			message("\nDECISION: Go to the next level \n", appendLF=F)
			split_l<- bestsplit(data=data_l, fixed=fixed, gvars=gvars, tgvars=tgvars, node.name=id_l, minbucket=minbucket, alpha=alpha)
			rsplit(data=data_l, fixed=fixed, gvars=gvars, tgvars=tgvars, id=id_l, split=split_l, alpha=alpha, minsplit=minsplit, minbucket=minbucket, Rate=Rate, loglik=loglik, env=env, coef.digits=coef.digits)
			}
		else {
			env$idlist<- c(env$idlist, list(unique(data_l[['patid']])))
			env$Treeout<- rbind(env$Treeout, c(id_l, left.subj, Rate, NA, NA, NA, loglik, NA))
			}
	
		##Right Node
		id_r=id*2+1
		message("---------------------------------------- \n", appendLF=F)
		message("NODE ", id_r, "- Rule:", s.var, " >=", s.cut, "\n", appendLF=F)   #Updated on 02-Jun-2019
		message("---------------------------------------- \n", appendLF=F)
		data_r<- subset(data, data[[s.var]]>=s.cut)                #Updated on 02-Jun-2019
		right.subj<- length(unique(data_r$patid))	
		message("No. of individual in right node: ", right.subj, " \n", appendLF=F)

	  	tjj<- NULL
		tjj<- summary(lme(fixed, data=data_r, method="REML", random=~1|patid))
                if(print.lme) cat("NODE ", id_r, "- Rule:", s.var, " >=", s.cut, "\n")
		if(print.lme) print(tjj$tTable)
		Rate<- coeff.Estimate(tjj, coef.digits=coef.digits)
		loglik<- logLik(lme(fixed, data=data_r, method="REML", random=~1|patid))
		
		if (right.subj>=minsplit)
			{
			message("\nDECISION: Go to the next level \n", appendLF=F)
			split_r<- bestsplit(data=data_r, fixed=fixed, gvars=gvars, tgvars=tgvars, node.name=id_r, minbucket=minbucket, alpha=alpha)
			rsplit(data=data_r, fixed=fixed, gvars=gvars, tgvars=tgvars, id=id_r, split=split_r, alpha=alpha, minsplit=minsplit, minbucket=minbucket, Rate=Rate, loglik=loglik, env=env, coef.digits=coef.digits)
			}
		else {
			env$idlist<- c(env$idlist, list(unique(data_r[['patid']])))
			env$Treeout<- rbind(env$Treeout, c(id_r, right.subj, Rate, NA, NA, NA, loglik, NA))
			}
		}
	else 	message("\nDECISION: NO more splitting required \n", appendLF=F)			
	}



  #------------Checks
  if(!exists(as.character(substitute(data)))) stop("Dataset does not exist\n")
  if(!is.data.frame(data)) stop("Dataset does not exist\n")
  
  #-------- Checks for patid variable
  if(!patid %in% colnames(data))  stop("The column ", patid, " containing subjects id is missing in dataset.\n")
  data$patid<- data[[patid]]
  #-------------Drop observations with missing patid
  data<- data[!is.na(data[["patid"]]),]
  
  #-------- Checks for Y variable
  Y.name<- as.character(attr(as.Formula(fixed), "lhs"))
  if(!Y.name %in% colnames(data))  
    stop("The column ", Y.name, " containing subjects id is missing in dataset.\n")
  #-------------Drop observations with missing Y and patid
  data<- data[!is.na(data[[Y.name]]),]
  
  #--------- Check whether length of gvars matches with tgvars
  if (length(gvars)!=length(tgvars))
    stop("gvars and tgvars are not of equal length. \n")
  
  #--------- Check tgvars does not include NA values
  if (any(is.na(tgvars)))
    stop("tgvars cannot have NA value. \n")

  #--------- Check whether all variables listed in gvars exist or not
  for(var.i in 1:length(gvars)){
      if(!gvars[var.i] %in% colnames(data))  stop("The column ", gvars[var.i], " is missing in dataset.\n")
    }

  #--------- Check whether one subject is assoicated with more than one partitioning variable.
  for(var.i in 1:length(gvars)){
    temp<- as.data.frame(list(patid=cbind(data["patid"], gvar=data[gvars[var.i]])))
    colnames(temp)=c("patid", "gvar")
    temp<- unique(temp)
    temp<- temp[!is.na("gvar"),]
    n.gvar.pat<- tapply(temp$gvar, temp$patid, function(x) length(unique(x)))
    if(any(n.gvar.pat!=1)) stop("One subject is associated with more than one value of ", gvars[var.i],".\n")
    rm(temp, n.gvar.pat)
  }
  
	LongCART.env<- new.env() #Define a new environment
	message("------------------------------------------ \n", appendLF=F)
	message("           ROOT NODE: NODE 1               \n", appendLF=F)
	message("------------------------------------------ \n", appendLF=F)
  
	tjj<- NULL
	tjj<- summary(lme(fixed, data=data, method="REML", random=~1|patid))
        if(print.lme) cat("ROOT NODE: NODE 1\n")
	if(print.lme) print(tjj$tTable)
	Rate<- coeff.Estimate(tjj, coef.digits=coef.digits)
	loglik<- logLik(lme(fixed, data=data, method="REML", random=~1|patid))

	split<- bestsplit(data=data, fixed=fixed, gvars=gvars, tgvars=tgvars, node.name=1, minbucket=minbucket, alpha=alpha)
	rsplit(data=data, fixed=fixed, gvars=gvars, tgvars=tgvars, id=1, split=split, alpha=alpha, minsplit=minsplit, minbucket=minbucket, Rate=Rate, loglik=loglik, env=LongCART.env, coef.digits=coef.digits)

	Treeout<- LongCART.env$Treeout
	colnames(Treeout)<- c("ID", "n", "yval", "var", 
                  	    "index", "p (Instability)", 'loglik', 'improve')
	row.names(Treeout)<- NULL

	Treeout[,1]<- as.numeric(Treeout[,1])
	Treeout[,2]<- as.numeric(Treeout[,2])
	Treeout[,5]<- round(as.numeric(Treeout[,5]), digits=2)
	Treeout[,6]<- round(as.numeric(Treeout[,6]), digits=3)
	Treeout[,7]<- round(as.numeric(Treeout[,7]), digits=0)
	Treeout[,8]<- round(as.numeric(Treeout[,8]), digits=0)
	Treeout$Terminal<- ifelse(is.na(Treeout[,5]), TRUE, FALSE)
	print(Treeout)
	
	#Statistics for measuring improvment in TREE
	##determine p
	form<-Formula(fixed)
	Y.name<- as.character(attr(form, "lhs"))
	X.name<- as.character(attr(form, "rhs"))
	X.name<- gsub(pattern=" ", replacement="", x=X.name)
	X.name<- unlist(strsplit(X.name, "+", fixed = TRUE))
	p<- length(X.name)+1

	AIC.tree<- 2*sum(Treeout[,7]*Treeout[,9])-2*(p+2)*sum(Treeout[,9])
	AIC.root<- 2*Treeout[1,7]-2*(p+2)
	improve.AIC<- AIC.tree-AIC.root
	logLik.tree<- sum(Treeout[,7]*Treeout[,9])
	logLik.root<- Treeout[1,7]
	Deviance<- 2*(logLik.tree-logLik.root)
	LRT.df<- p*sum(Treeout[,9])-p
	LRT.p<- 1-pchisq(Deviance, LRT.df)

	message('AIC(tree)=', AIC.tree, '   AIC(Root)=', AIC.root, '\n', appendLF=F)
	message('logLikelihood (tree)=', logLik.tree, '   logLikelihood (Root)=', logLik.root, '\n', appendLF=F)
	message('Deviance=', Deviance, ' (df=', LRT.df, ', p-val=', LRT.p, ') \n', appendLF=F)

	#Assigning information to the subject
	sel<- as.numeric(rownames(Treeout[Treeout$Terminal,]))
	idlist<- LongCART.env$idlist
	subj.class<- NULL
	for(i in 1:length(sel)) {
		temp<-Treeout[sel[i], c(1,3,4)]
		rownames(temp)<- NULL
		subj.class<- rbind(subj.class, cbind(idlist[[sel[i]]], temp))
		}
	
	subj.class<- as.data.frame(subj.class)
	names(subj.class)<- c('patid', 'node', 'n', 'yval')

	#--- Create FRAME, SPLITS, CPTABLE and FUNCTIONS objects for using PLOT.RPART
	Treeout1<- as.data.frame(Treeout[,c("ID", "var", "n", "yval", 
	                                    "loglik", "Terminal", "index", "improve")])
	row.names(Treeout1)<- Treeout1$ID
	Treeout1$var=ifelse(Treeout1$Terminal, "<leaf>",Treeout1$var)
	Treeout1$dev=-Treeout1$loglik
	Treeout1$wt<- Treeout1$count<- Treeout1$n
	Treeout1$ncompete<- Treeout1$nsurrogate<- 0 
	Treeout1$ncat=-1
	
	frame<- Treeout1[c("var", "n", "wt", "dev", "yval", "ncompete", "nsurrogate")]
	
	splits<- Treeout1[!Treeout1$Terminal,]
	
	splits<- splits[c("count", "ncat", "improve", "index")]
	splits<- as.matrix(splits)

	temp<- frame[frame$var=="<leaf>",]
	cptable<- 0:(dim(temp)[1]-1)
	
	functions<- NULL
	functions$text<- function (yval, dev, wt, ylevel, digits, n, use.n) {
	  if (use.n) 
	    paste0(yval, "\nn=", n)
	  else yval
	}
	
	ret<- list(Treeout, frame, splits, cptable, functions,
	           p, AIC.tree, AIC.root, improve.AIC, logLik.tree, 
             logLik.root, Deviance, LRT.df, LRT.p, subj.class)
	names(ret)<- c('Treeout', 'frame', 'splits', 'cptable', 'functions', 
                'p', 'AIC.tree', 'AIC.root', 'improve.AIC', 'logLik.tree', 
                     'logLik.root', 'Deviance', 'LRT.df', 'LRT.pval', 'subj.class')
	class(ret)<- c("rpart", "LongCART")
	ret
	}

#---End of code



