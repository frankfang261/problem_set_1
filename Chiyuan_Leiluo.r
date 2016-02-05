library(foreign)
library(sandwich)
library(stargazer)
library(lmtest)
library(multiwayvcov)
library(AER)
library(Matching)
############################
set.seed(1234)
Dxc=rnorm(100)   #Dc star
Dc=vector(length=100)
nc=rnorm(100)  
for (i in 1:100)
{
  Dc[i]=0
  if(Dxc[i]>0)
    Dc[i]=1
  
}
Y=vector(length=10000)
for (i in 1:100)
{
  for(j in 1:100)
    Y[100*(i-1)+j]=Dc[i]+nc[i]+Dc[i]*nc[i]+rnorm(1)
}
cDc=vector(length=10000)
cnc=vector(length=10000)
for (i in 1:100)
{
  for(j in 1:100)
  {
    cDc[100*(i-1)+j]=Dc[i]
    cnc[100*(i-1)+j]=nc[i]
  }
}
state=vector(length=10000)
for (i in 1:100)
{
  for(j in 1:100)
  {
    state[100*(i-1)+j]=i
    
  }
}
myData=data.frame(cbind(Y,cDc,cnc,state))





###a
rega=lm(Y~cDc,data=myData)
summary(rega)




###b
#regb=rlm(Y~cDc,data=myData) different approach from Stata's robust
#summary(regb)

coeftest(rega,vcov=vcovHC(rega,"HC1")) # reproduces STATA robust standard error


###c
rega.vcov=cluster.vcov(rega,myData$state)  #clustered covariance matrix
coeftest(rega,rega.vcov)                    #coeefficient t test


###e
collapsedY=vector(length=100)
for (i in 1:100)
{
  s=0
  for(j in 1:100)
  {
    s=s+Y[100*(i-1)+j]
    
  }
  collapsedY[i]=s/100
}
state2=(1:100)

myData2=data.frame(cbind(collapsedY,Dc,state2))
names(myData2)=c('collapsedY','cDc','state')
rege=lm(collapsedY~cDc,data=myData2)
summary(rege)

###f

coeftest(rege,vcov=vcovHC(rege,"HC1"))   #STATA robust standard error
coeftest(rege,vcov=vcovHC(rege,"HC2"))   #HC1 standard error
coeftest(rege,vcov=vcovHC(rege,"HC3"))   #HC3 standard error

#g monte carlo simlation
testa=0    #num of rejections of Beta_Dc=1 for each case
testb=0
testc=0
testd=0
teste=0
for (ll in 1:1000)
{
  #data generation
  Dxc=rnorm(100)   #Dc star
  Dc=vector(length=100)
  nc=rnorm(100)  
  for (i in 1:100)
  {
    Dc[i]=0
    if(Dxc[i]>0)
      Dc[i]=1
    
  }
  Y=vector(length=10000)
  for (i in 1:100)
  {
    for(j in 1:100)
      Y[100*(i-1)+j]=Dc[i]+nc[i]+Dc[i]*nc[i]+rnorm(1)
  }
  cDc=vector(length=10000)
  cnc=vector(length=10000)
  for (i in 1:100)
  {
    for(j in 1:100)
    {
      cDc[100*(i-1)+j]=Dc[i]
      cnc[100*(i-1)+j]=nc[i]
    }
  }
  state=vector(length=10000)
  for (i in 1:100)
  {
    for(j in 1:100)
    {
      state[100*(i-1)+j]=i
      
    }
  }
  myData=data.frame(cbind(Y,cDc,cnc,state))
  collapsedY=vector(length=100)
  for (i in 1:100)
  {
    s=0
    for(j in 1:100)
    {
      s=s+Y[100*(i-1)+j]
      
    }
    collapsedY[i]=s/100
  }
  state2=(1:100)
  myData2=data.frame(cbind(collapsedY,Dc,state2))
  names(myData2)=c('collapsedY','cDc','state')
  
  
  #t test 1
  rega=lm(Y~cDc,data=myData)
  rega1=summary(rega)
  
  if(abs(rega1$coefficients[2,1]-1)/rega1$coefficients[2,2]>1.96)
  {   testa=testa+1} 
  
  
  #t test 2
  rega.vcov=cluster.vcov(rega,myData$state)  #clustered covariance matrix
  regb1=coeftest(rega,rega.vcov) 
  if(abs(regb1[2,1]-1)/regb1[2,2]>1.96)
  {   testb=testb+1} 
  
  #t test 3
  rege=lm(collapsedY~cDc,data=myData2)
  rege1=summary(rege)
  rege1
  if(abs(rege1$coefficients[2,1]-1)/rege1$coefficients[2,2]>1.96)
  {   testc=testc+1} 
  
  #t test 4
  regf1 = coeftest(rege,vcov=vcovHC(rege,"HC1"))
  if(abs(regf1[2,1]-1)/regf1[2,2]>1.96)
  {   testd=testd+1} 
  
  #t test 5
  rege1=coeftest(rege,vcov=vcovHC(rege,"HC2"))
  
  if(abs(rege1[2,1]-1)/rege1[2,2]>1.96)
  {   teste=teste+1} 
  
  #  if (ll %% 10 == 0) print(ll/10)  #counter percentage has been done
}
testa/1000
testb/1000
testc/1000
testd/1000
teste/1000



#####################
####Question 2#######
#####################

DataQ2=read.table("D:/R Workspace/dataset/PSET1-2016.csv",sep=",",header=T)
summary(DataQ2)

####Question 2a)
    ####summarizing the data

stargazer(DataQ2,type="latex",out="D:/R Workspace/PS1/data2.tex")
#####Question 2bi) reduced form

reg2b.i=lm(vote02~treat_real+age+female+vote00+vote98+newreg+persons,data=DataQ2)
summary(reg2b.i)

stargazer(reg2b.i,type="latex",out = "D:/R Workspace/PS1/reduced form.tex")

coeftest(reg2b.i,vcov=vcovHC(reg2b.i,"HC1")) #STATA robust standard error
#####Question 2bii)  first stage
reg2b.ii=lm(contact~treat_real+age+female+vote00+vote98+newreg+persons,data=DataQ2)
stargazer(reg2b.ii,type="latex",out = "D:/R Workspace/PS1/first stage.tex")
#####Wald Estimator

DataQ2.wald.Tre=na.omit(subset(DataQ2,DataQ2$treat_real>0))
DataQ2.wald.NTre=na.omit(subset(DataQ2,DataQ2$treat_real==0))

wald.estimator=(mean(predict(reg2b.i,DataQ2.wald.Tre))-mean(predict(reg2b.i,DataQ2.wald.NTre)))/(mean(predict(reg2b.ii,DataQ2.wald.Tre)-mean(predict(reg2b.ii,DataQ2.wald.NTre))))

##### control function approach

reg2b.iv=lm(vote02~reg2b.ii$residuals+age+female+vote00+vote98+newreg+persons,data=DataQ2.matching)


####Question 2biii)  2SLS
reg2b.iii=ivreg(vote02~contact+age+female+vote00+vote98+newreg+persons|treat_real+age+female+vote00+vote98+newreg+persons,data=DataQ2)

summary(reg2b.iii)
stargazer(reg2b.iii,type="latex",out = "D:/R Workspace/PS1/ivreg.tex")
####Question 2c)
reg2c=lm(vote02~contact,data=DataQ2)
reg2c.controls=lm(vote02~contact+age+female+vote00+vote98+newreg+persons,data=DataQ2)
reg2c.controls2=lm(vote02~contact+age+female+vote00+vote98+newreg+persons
                   +contact*age+contact*female+contact*vote00+contact*vote98
                   +contact*newreg+contact*persons,data=DataQ2)

stargazer(reg2c,reg2c.controls,reg2c.controls2,type = "latex",out = "D:/R Workspace/PS1/ols.tex")

####Question 2di)
  # generating propensity score
ps.2di=glm(contact~age+female+vote00+vote98+newreg+persons,data=DataQ2,family =binomial(link = "logit") )
DataQ2.matching=DataQ2
DataQ2.matching$pscore=predict(ps.2di,DataQ2.matching,type="response")
DataQ2.matching=na.omit(DataQ2.matching)


####Kernal density
density1=density(DataQ2.matching$pscore,kernel = "epanechnikov")
density2=density(DataQ2.matching$pscore,bw=0.05,kernel = "epanechnikov")
density3=density(DataQ2.matching$pscore,kernel = "triangular")

DataQ2.contact=subset(DataQ2.matching,contact==1)
DataQ2.nocontact=subset(DataQ2.matching,contact==0)
DataQ2.matching$weights=DataQ2.matching$contact+(1-DataQ2.matching$contact)*DataQ2.matching$pscore/(1-DataQ2.matching$pscore)
#####weighted least square
reg2.wls=lm(vote02~contact+age+female+vote00+vote98+newreg+persons,data=DataQ2.matching,weights = DataQ2.matching$weights)
summary(reg2.wls)
stargazer(reg2.wls,type="latex",out = "D:/R Workspace/PS1/wls.tex")
a=mean(DataQ2.contact$vote02)
b=sum(DataQ2.nocontact$vote02*DataQ2.nocontact$pscore/(1-DataQ2.nocontact$pscore))
beta=a-b/sum(DataQ2.nocontact$pscore/(1-DataQ2.nocontact$pscore))
##hence beta is 0.04257
density4=density(DataQ2.contact$pscore,kernel = "epanechnikov")
density5=density(DataQ2.nocontact$pscore,kernel = "epanechnikov")
plot(density5,col=1)           #propensity score kernel density of no contact
lines(density4,col=3)          #propensity score kernel density of contact
####question 2dv)
match.out1=matchit(contact~age+female+vote00+vote98+newreg+persons,data=DataQ2.matching,
                  method = "nearest")
summary(match.out1)
match.out1.ps=matchit(contact~pscore,data=DataQ2.matching,method="nearest")
summary(match.out1[match.out1$nn])


#########################################################################
###not sure how three other match methods are done in R##################
#########################################################################

match.result.ATT=Match(DataQ2.matching$vote02,DataQ2.matching$contact,DataQ2.matching$pscore,estimand = "ATT")
match.result.ATE=Match(DataQ2.matching$vote02,DataQ2.matching$contact,DataQ2.matching$pscore,estimand = "ATE")
summary(match.result.ATT)
plot(density1,main="ep kernel plot")
plot(density2,main="ep kernel bandwith 0.05 plot")
plot(density3,main="ep kernel tri plot")


plot(density4,main="control and treated plot")
lines(density5,type="h")

