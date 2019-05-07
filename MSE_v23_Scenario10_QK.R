# Operating model for LA, MS, and AL Spotted Seatrout

# Nate Jermain
# June 29th, 2018

####### packages needed######### 
library(r4ss)




############## Load recruitment time series from control model ######################
setwd("C:/Users/w10007346/Documents/MSE")
load("Control.master.output.v23.RData")


La.Rzeros<-master.output$La.recruits
Ms.Rzeros<-master.output$Ms.recruits
Al.Rzeros<-master.output$Al.recruits 


############ List output for performance indices #########################
master.output<-list(list(Al.recruits=list(), Ms.recruits=list(), La.recruits=list(),
                         Al.harvest=list(), Ms.harvest=list(), La.harvest=list(),
                         Al.SSB=list(), Ms.SSB=list(), La.SSB=list(), Al.RSD=list(),
                         Ms.RSD=list(), La.RSD=list()), Al.regulations=list(), 
                    Ms.regulations=list(), La.regulations=list(), Al.gradient=list(), 
                    Ms.gradient=list(), La.gradient=list(), Al.SPR=list(), Ms.SPR=list(), 
                    La.SPR=list(), Al.SPRTAIL=list(), Ms.SPRTAIL=list(), La.SPRTAIL=list())



####### Arguments for MSE Run ###########################################
# hatchery release number
relquant<-5000000
# SPR deviations
SPRTargets<-c(18,20,30)
startMLL<-c(12,15,14)
# time series
projection<-60

# Discard mortality rate 
DiscMort<-.1

# Maturity
maturity<-c(0,.8,1,1,1,1,1)


# Minimum length limit to selectivity harvest
MLL.11<-c(0,1,1,1,1,1,1,1)
MLL.12<-c(0,.9,1,1,1,1,1,1)
MLL.13<-c(0,.8,1,1,1,1,1,1)
MLL.14<-c(0,.6,1,1,1,1,1,1)
MLL.15<-c(0,.3,.9,1,1,1,1,1)
MLL.16<-c(0,.1,.7,.9,1,1,1,1)

# Release selectivity
R.MLL.11<-1-MLL.11
R.MLL.12<-1-MLL.12
R.MLL.13<-1-MLL.13
R.MLL.14<-1-MLL.14
R.MLL.15<-1-MLL.15
R.MLL.16<-1-MLL.16

# MLL as lists 
Lengthlimit<-list(MLL.11,MLL.12, MLL.13, MLL.14, MLL.15, MLL.16)
ReleaseLlimit<-list(R.MLL.11,R.MLL.12, R.MLL.13, R.MLL.14, R.MLL.15, R.MLL.16)


# Assign the starting minimum length limit to carrier that will change after assessments
la.new.MLL<-MLL.12
ms.new.MLL<-MLL.15
al.new.MLL<-MLL.14

# empty vectors for tracking regulation changes
al.regulations<-0
ms.regulations<-0
la.regulations<-0

al.Rel.regulations<-0
ms.Rel.regulations<-0
la.Rel.regulations<-0

# empty vectors for tracking gradients through the assessment model
al.gradient<-0
ms.gradient<-0
la.gradient<-0

# empty vectors for tracking terminal SPR values
Al.T.SPR<-0
Ms.T.SPR<-0
La.T.SPR<-0

# Natural mortality
M<-c(.875,.495,.385,.33,.30,.285,.275,.27)

##################### Hatchery Input #################################################
# Length at release 30 cm
# 1 million released by each state
# how many die before becoming age0?
# Use length-age relationship from Powell et al. 2004

# July 1st birthday (La Assessment 2014), grow linearly for 6 months
wkbins<-seq(38,184,7) # only interested in lengths after 30mm so start with 38 days old, determined from L-A formula

Jlen<-0
for(i in 1:length(wkbins)){
  Jlen[i]<--10.56+.8834*wkbins[i]
}

# convert standard length to total length using Hein et al. 1980

Jtlen<-0
for(i in 1:length(wkbins)){
  Jtlen[i]<-(Jlen[i]+2.0520)/.8369
}
Jtlen

# grow based on vonbert from age .5 years to 1


L<-0
for (i in 1:26 ){
  linf<-54.7   # Parameter values from Alabama FI sampling 
  k<-.33
  to<-1.17
  L[i]<-linf*(1-exp(-k*(seq(.5,1,length.out=26)[i]+to)))
}
L


Jtlen<-c(Jtlen,L/.1)

# Natural mortality using Lorenzen 2005 estimate where average M for 1cm fish is 15 yr
Jtlen<-Jtlen*.1 # mm to cm
MJ<-0
M1<-15/52.143 # divide 15/yr to weeks
for (i in 1:length(Jtlen)){
  MJ[i]<-M1*(1/Jtlen[i])
}


# survival in the wild based on natural mortality MJ
NJt<-relquant
for (i in 1:length(Jtlen)){
  NJt[i+1]<-NJt[i]*exp(-MJ[i])
}

hatcheryA1<-(NJt[length(NJt)])/1000 #number (thousands) of stocked age1 recruits per state

########### START SIMULATIONS FOR LOOP ##############
# start timer
library(tictoc)
tic("timer")

for (s in 1:1000){

################# Calculate first year's numbers at age ################
#Louisiana
la.Agebegin<-La.Rzeros[[s]][1]
la.Select<-MLL.12
la.R.Select<-R.MLL.12
la.FHmult<-1.69
la.FRmult<-.03
la.FHmort<-0
la.FRmort<-0
la.Catchbegin<-0
la.DDiscardsbegin<-0

for (i in 1:7){ # 6 year plus group 
  # build vectors 
  la.FHmort[i]<-la.FHmult*la.Select[i]
  la.FRmort[i]<-la.FRmult*la.R.Select[i]
  if (i==2){la.Agebegin[2]<-la.Agebegin[2]+hatcheryA1}
  la.Agebegin[i+1]<-la.Agebegin[i]*exp(-(M[i]+la.FRmort[i]+la.FHmort[i]))
  la.Catchbegin[i]<-(la.Agebegin[i]-la.Agebegin[i+1])*(la.FHmort[i]/(la.FHmort[i]+la.FRmort[i]+M[i]))
  la.DDiscardsbegin[i]<-(la.Agebegin[i]-la.Agebegin[i+1])*(la.FRmort[i]/(la.FHmort[i]+la.FRmort[i]+M[i]))
}


# Beginning values to be discarded after 1st run
# put in a master list of output
la.out<-list(NumbersatAge=list(), 
             Catch=list(), 
             DDiscards=list())

la.out$NumbersatAge[[1]]<-la.Agebegin
la.out$Catch[[1]]<-la.Catchbegin
la.out$DDiscards[[1]]<-la.DDiscardsbegin

# Mississippi

ms.Agebegin<-Ms.Rzeros[[s]][1]
ms.Select<-MLL.15
ms.R.Select<-R.MLL.15
ms.FHmult<-1.45
ms.FRmult<-.04
ms.FHmort<-0
ms.FRmort<-0
ms.Catchbegin<-0
ms.DDiscardsbegin<-0

for (i in 1:7){ # 6 year plus group
  # build vectors 
  ms.FHmort[i]<-ms.FHmult*ms.Select[i]
  ms.FRmort[i]<-ms.FRmult*ms.R.Select[i]
  if (i==2){ms.Agebegin[2]<-ms.Agebegin[2]+hatcheryA1}
  ms.Agebegin[i+1]<-ms.Agebegin[i]*exp(-(M[i]+ms.FRmort[i]+ms.FHmort[i]))
  ms.Catchbegin[i]<-(ms.Agebegin[i]-ms.Agebegin[i+1])*(ms.FHmort[i]/(ms.FHmort[i]+ms.FRmort[i]+M[i]))
  ms.DDiscardsbegin[i]<-(ms.Agebegin[i]-ms.Agebegin[i+1])*(ms.FRmort[i]/(ms.FHmort[i]+ms.FRmort[i]+M[i]))
}


# put in a master list of output
ms.out<-list(NumbersatAge=list(), 
             Catch=list(), 
             DDiscards=list())

ms.out$NumbersatAge[[1]]<-ms.Agebegin
ms.out$Catch[[1]]<-ms.Catchbegin
ms.out$DDiscards[[1]]<-ms.DDiscardsbegin

# Alabama

al.Agebegin<-Al.Rzeros[[s]][1]
al.Select<-MLL.14
al.R.Select<-R.MLL.14
al.FHmult<-.99
al.FRmult<-.05
al.FHmort<-0
al.FRmort<-0
al.Catchbegin<-0
al.DDiscardsbegin<-0

for (i in 1:7){ # 6 year plus group
  # build vectors 
  al.FHmort[i]<-al.FHmult*al.Select[i]
  al.FRmort[i]<-al.FRmult*al.R.Select[i]
  if (i==2){al.Agebegin[2]<-al.Agebegin[2]+hatcheryA1}
  al.Agebegin[i+1]<-al.Agebegin[i]*exp(-(M[i]+al.FRmort[i]+al.FHmort[i]))
  al.Catchbegin[i]<-(al.Agebegin[i]-al.Agebegin[i+1])*(al.FHmort[i]/(al.FHmort[i]+al.FRmort[i]+M[i]))
  al.DDiscardsbegin[i]<-(al.Agebegin[i]-al.Agebegin[i+1])*(al.FRmort[i]/(al.FHmort[i]+al.FRmort[i]+M[i]))
}


# put in a master list of output
al.out<-list(NumbersatAge=list(), 
             Catch=list(), 
             DDiscards=list())

al.out$NumbersatAge[[1]]<-al.Agebegin
al.out$Catch[[1]]<-al.Catchbegin
al.out$DDiscards[[1]]<-al.DDiscardsbegin

########## Reset MLLs to original starting values for the first year  #############

la.new.MLL<-MLL.12
ms.new.MLL<-MLL.15
al.new.MLL<-MLL.14

######## Copy uncorrupted .dat template files to the assessment model folder #################


current.folder<-"C:/Users/w10007346/Documents/Assesment Model/Uncorrupted data files"
ms.folder<-"C:/Users/w10007346/Documents/Assesment Model/Mississippi"
la.folder<-"C:/Users/w10007346/Documents/Assesment Model/Louisiana"
al.folder<-"C:/Users/w10007346/Documents/Assesment Model/Alabama"

ms.file<-list.files(current.folder, "MSss3.dat")
file.copy("C:/Users/w10007346/Documents/Assesment Model/Uncorrupted data files/MSss3.dat", ms.folder, overwrite = T)
file.copy("C:/Users/w10007346/Documents/Assesment Model/Uncorrupted data files/LAss3.dat", la.folder, overwrite = T)
file.copy("C:/Users/w10007346/Documents/Assesment Model/Uncorrupted data files/ALss3.dat", al.folder, overwrite = T)

########## FOR LOOP BEGINS HERE #################################################

for (y in 1:projection){
  
  
  
  ########## Louisiana ####################################################
  
  
  # if not the initial year, build the population using last years Num-at-age 
  # and new recruits now named la.pop
  if(y>1){
    # Natural mortality using only female specific right now
    la.Select<-as.vector(la.new.MLL)
    la.R.Select<-as.numeric(as.character(unlist(ReleaseLlimit[which(sapply(Lengthlimit, function(x) all(la.new.MLL %in% x))==TRUE)])))
    la.FHmort<-0
    la.FRmort<-0
    la.Catch<-0
    la.DDiscards<-0
    
    for (i in 1:7){ # 6 year plus group, all dead after inclusion into that group
      # build vectors 
      la.FHmort[i]<-la.FHmult*la.Select[i]
      la.FRmort[i]<-la.FRmult*la.R.Select[i]
      if (i==2){la.pop[2]<-la.pop[2]+hatcheryA1}
      la.pop[i+1]<-la.out$NumbersatAge[[y-1]][i]*exp(-(M[i]+la.FRmort[i]+la.FHmort[i]))
      if (i==7){la.pop[7]<-la.pop[7]+la.pop[8]}
      la.Catch[i]<-(la.out$NumbersatAge[[y-1]][i]-la.pop[i+1])*(la.FHmort[i]/(la.FHmort[i]+la.FRmort[i]+M[i]))
      la.DDiscards[i]<-(la.out$NumbersatAge[[y-1]][i]-la.pop[i+1])*(la.FRmort[i]/(la.FHmort[i]+la.FRmort[i]+M[i]))
    }
    la.LDiscards<-la.DDiscards/.1  #get number of live discards
    la.NewDDiscards<-la.LDiscards*DiscMort # the new dead discards 
    # update population with the difference between old and new dead discards
    la.pop<-(la.DDiscards-la.NewDDiscards)+la.pop[1:7]
    
    
    # transfer this years information to the list
    la.out$NumbersatAge[[y]]<-la.pop
    la.out$Catch[[y]]<-la.Catch
    la.out$DDiscards[[y]]<-la.NewDDiscards
    
  }
  
  
  ########## Mississippi ####################################################
  
  
  # if not the initial year, build the population using last years Num-at-age 
  # and new recruits now named la.pop
  if(y>1){
    # Natural mortality using only female specific right now
    ms.Select<-as.vector(ms.new.MLL)
    ms.R.Select<-as.numeric(as.character(unlist(ReleaseLlimit[which(sapply(Lengthlimit, function(x) all(ms.new.MLL %in% x))==TRUE)])))
    ms.FHmort<-0
    ms.FRmort<-0
    ms.Catch<-0
    ms.DDiscards<-0
    
    for (i in 1:7){ # 6 year plus group, all dead after inclusion into that group
      # build vectors 
      ms.FHmort[i]<-ms.FHmult*ms.Select[i]
      ms.FRmort[i]<-ms.FRmult*ms.R.Select[i]
      if (i==2){ms.pop[2]<-ms.pop[2]+hatcheryA1}
      ms.pop[i+1]<-ms.out$NumbersatAge[[y-1]][i]*exp(-(M[i]+ms.FRmort[i]+ms.FHmort[i]))
      if (i==7){ms.pop[7]<-ms.pop[7]+ms.pop[8]}
      ms.Catch[i]<-(ms.out$NumbersatAge[[y-1]][i]-ms.pop[i+1])*(ms.FHmort[i]/(ms.FHmort[i]+ms.FRmort[i]+M[i]))
      ms.DDiscards[i]<-(ms.out$NumbersatAge[[y-1]][i]-ms.pop[i+1])*(ms.FRmort[i]/(ms.FHmort[i]+ms.FRmort[i]+M[i]))
    }
    ms.LDiscards<-ms.DDiscards/.1  #get number of live discards
    ms.NewDDiscards<-ms.LDiscards*DiscMort # the new dead discards 
    # update population with the difference between old and new dead discards
    ms.pop<-(ms.DDiscards-ms.NewDDiscards)+ms.pop[1:7]
    
    
    # transfer this years information to the list
    ms.out$NumbersatAge[[y]]<-ms.pop
    ms.out$Catch[[y]]<-ms.Catch
    ms.out$DDiscards[[y]]<-ms.NewDDiscards
    
  }
  
  
  ########## Alabama ####################################################
  
  
  # if not the initial year, build the population using last years Num-at-age 
  # and new recruits now named la.pop
  if(y>1){
    # Natural mortality using only female specific right now
    al.Select<-as.vector(al.new.MLL)
    al.R.Select<-as.numeric(as.character(unlist(ReleaseLlimit[which(sapply(Lengthlimit, function(x) all(al.new.MLL %in% x))==TRUE)])))
    al.FHmort<-0
    al.FRmort<-0
    al.Catch<-0
    al.DDiscards<-0
    
    for (i in 1:7){ # 6 year plus group, all dead after inclusion into that group
      # build vectors 
      al.FHmort[i]<-al.FHmult*al.Select[i]
      al.FRmort[i]<-al.FRmult*al.R.Select[i]
      if (i==2){al.pop[2]<-al.pop[2]+hatcheryA1}
      al.pop[i+1]<-al.out$NumbersatAge[[y-1]][i]*exp(-(M[i]+al.FRmort[i]+al.FHmort[i]))
      if (i==7){al.pop[7]<-al.pop[7]+al.pop[8]}
      al.Catch[i]<-(al.out$NumbersatAge[[y-1]][i]-al.pop[i+1])*(al.FHmort[i]/(al.FHmort[i]+al.FRmort[i]+M[i]))
      al.DDiscards[i]<-(al.out$NumbersatAge[[y-1]][i]-al.pop[i+1])*(al.FRmort[i]/(al.FHmort[i]+al.FRmort[i]+M[i]))
    }
    al.LDiscards<-al.DDiscards/.1  #get number of live discards
    al.NewDDiscards<-al.LDiscards*DiscMort # the new dead discards 
    # update population with the difference between old and new dead discards
    al.pop<-(al.DDiscards-al.NewDDiscards)+al.pop[1:7]
    
    
    # transfer this years information to the list
    al.out$NumbersatAge[[y]]<-al.pop
    al.out$Catch[[y]]<-al.Catch
    al.out$DDiscards[[y]]<-al.NewDDiscards
    
  }
  
  
  
  
  
  ###### Stock Recruitment ##################################################
# get next years recruitment based on the control run time series 
  
  if (y<projection){
    la.pop<-La.Rzeros[[s]][y+1]
    ms.pop<-Ms.Rzeros[[s]][y+1]
    al.pop<-Al.Rzeros[[s]][y+1]
  }
  

############### Estimation component #########################################################


  
  
### Mississippi
  setwd("C:/Users/w10007346/Documents/Assesment model/Mississippi")
if (y==1){MStp<-SS_readdat("MSss3.dat", verbose=T, echoall=T) }

  # start at year 1
  MStp$styr<-1
  
  # now update list with replacement for the first year, and append each successive year
  
  MStp$endyr<-y
  MStp$N_catch<-y
  #Catch
  if (y==1){
    MStp$init_equil<-sum(ms.out$Catch[[1]])
    firstcatch<-c(sum(ms.out$Catch[[1]]),1,1)
    MStp$catch<-as.vector(firstcatch)
  } else{
    newcatch<-c(sum(ms.out$Catch[[y]]),y,1) # y is the yth year in the projection
    MStp$catch<-as.data.frame(rbind(MStp$catch,as.vector(newcatch)))
  }
  MStp$N_cpue<-y
  # INDICES
  if (y==1){
    firstindex<-c(1,1,2,sum(ms.out$NumbersatAge[[1]]), .25)
    MStp$CPUE<-as.vector(firstindex)
  } else{
    newindices<-c(y,1,2,sum(ms.out$NumbersatAge[[y]]), .25)
    MStp$CPUE<-as.data.frame(rbind(MStp$CPUE,as.vector(newindices)))
  }
  # Age composition
  MStp$N_agecomp<-y*2 #number of observations
  ### fishing fleet
  if (y==1){
    firstagecomp<-c(1,1,1,0,0,1,-1,-1,sum(ms.out$Catch[[1]]), ms.out$Catch[[1]], rep(999,7))
    MStp$agecomp<-as.vector(firstagecomp)
  } else{
    newagecomp<-c(y,1,1,0,0,1,-1,-1,sum(ms.out$Catch[[y]]), ms.out$Catch[[y]], rep(999,7))
    MStp$agecomp<-as.data.frame(rbind(MStp$agecomp,as.vector(newagecomp)))
  }
  ### Survey
  if (y==1){
    firstsurveyagecomp<-c(1,1,2,0,0,1,-1,-1,sum(ms.out$NumbersatAge[[1]][2:7]), c(0,ms.out$NumbersatAge[[1]][2:7]), rep(999,7))
    MStp$agecomp<-as.data.frame(rbind(MStp$agecomp,as.vector(firstsurveyagecomp)))
  } else{
    newsurveyagecomp<-c(y,1,2,0,0,1,-1,-1,sum(ms.out$NumbersatAge[[y]][2:7]), c(0,ms.out$NumbersatAge[[y]][2:7]), rep(999,7))
    MStp$agecomp<-as.data.frame(rbind(MStp$agecomp,as.vector(newsurveyagecomp)))
  }
  


# only write data file for years after burn in 
# write list to .dat for specific assessment years
if (y==30|y==33|y==36|y==39|y==42|y==45|y==48|y==51|y==54|y==57){SS_writedat(MStp, "C:/Users/w10007346/Documents/Assesment model/Mississippi/MSss3.dat", overwrite=T)
  }




#### Louisiana

setwd("C:/Users/w10007346/Documents/Assesment model/Louisiana")
if (y==1){LAtp<-SS_readdat("LAss3.dat", verbose=T, echoall=T)} 


  # start at year 1
  LAtp$styr<-1
  
  # now update list with replacement for the first year, and append each successive year
  
  LAtp$endyr<-y
  LAtp$N_catch<-y
  #Catch
  if (y==1){
    LAtp$init_equil<-sum(la.out$Catch[[1]])
    firstcatch<-c(sum(la.out$Catch[[1]]),1,1)
    LAtp$catch<-as.vector(firstcatch)
  } else{
    newcatch<-c(sum(la.out$Catch[[y]]),y,1) # y is the yth year in the projection
    LAtp$catch<-as.data.frame(rbind(LAtp$catch,as.vector(newcatch)))
  }
  LAtp$N_cpue<-y
  # INDICES
  if (y==1){
    firstindex<-c(1,1,2,sum(la.out$NumbersatAge[[1]]), .25)
    LAtp$CPUE<-as.vector(firstindex)
  } else{
    newindices<-c(y,1,2,sum(la.out$NumbersatAge[[y]]), .25)
    LAtp$CPUE<-as.data.frame(rbind(LAtp$CPUE,as.vector(newindices)))
  }
  # Age composition
  LAtp$N_agecomp<-y*2 #number of observations
  ### fishing fleet
  if (y==1){
    firstagecomp<-c(1,1,1,0,0,1,-1,-1,sum(la.out$Catch[[1]]), la.out$Catch[[1]], rep(999,7))
    LAtp$agecomp<-as.vector(firstagecomp)
  } else{
    newagecomp<-c(y,1,1,0,0,1,-1,-1,sum(la.out$Catch[[y]]), la.out$Catch[[y]], rep(999,7))
    LAtp$agecomp<-as.data.frame(rbind(LAtp$agecomp,as.vector(newagecomp)))
  }
  ### Survey
  if (y==1){
    firstsurveyagecomp<-c(1,1,2,0,0,1,-1,-1,sum(la.out$NumbersatAge[[1]][2:7]), c(0,la.out$NumbersatAge[[1]][2:7]), rep(999,7))
    LAtp$agecomp<-as.data.frame(rbind(LAtp$agecomp,as.vector(firstsurveyagecomp)))
  } else{
    newsurveyagecomp<-c(y,1,2,0,0,1,-1,-1,sum(la.out$NumbersatAge[[y]][2:7]), c(0,la.out$NumbersatAge[[y]][2:7]), rep(999,7))
    LAtp$agecomp<-as.data.frame(rbind(LAtp$agecomp,as.vector(newsurveyagecomp)))
  }
  
  # only write data file for years after burn in 
if (y==30|y==33|y==36|y==39|y==42|y==45|y==48|y==51|y==54|y==57){
# write list to .dat
SS_writedat(LAtp, "C:/Users/w10007346/Documents/Assesment model/Louisiana/LAss3.dat", overwrite=T)
 
  }


###### Alabama

setwd("C:/Users/w10007346/Documents/Assesment model/Alabama")
if (y==1){ALtp<-SS_readdat("ALss3.dat", verbose=T, echoall=T)} 



  # start at year 1
  ALtp$styr<-1
  
  # now update list with replacement for the first year, and append each successive year
  
  ALtp$endyr<-y
  ALtp$N_catch<-y
  #Catch
  if (y==1){
    ALtp$init_equil<-sum(al.out$Catch[[1]])
    firstcatch<-c(sum(al.out$Catch[[1]]),1,1)
    ALtp$catch<-as.vector(firstcatch)
  } else{
    newcatch<-c(sum(al.out$Catch[[y]]),y,1) # y is the yth year in the projection
    ALtp$catch<-as.data.frame(rbind(ALtp$catch,as.vector(newcatch)))
  }
  ALtp$N_cpue<-y
  # INDICES
  if (y==1){
    firstindex<-c(1,1,2,sum(al.out$NumbersatAge[[1]]), .25)
    ALtp$CPUE<-as.vector(firstindex)
  } else{
    newindices<-c(y,1,2,sum(al.out$NumbersatAge[[y]]), .25)
    ALtp$CPUE<-as.data.frame(rbind(ALtp$CPUE,as.vector(newindices)))
  }
  # Age composition
  ALtp$N_agecomp<-y*2 #number of observations
  ### fishing fleet
  if (y==1){
    firstagecomp<-c(1,1,1,0,0,1,-1,-1,sum(al.out$Catch[[1]]), al.out$Catch[[1]], rep(999,7))
    ALtp$agecomp<-as.vector(firstagecomp)
  } else{
    newagecomp<-c(y,1,1,0,0,1,-1,-1,sum(al.out$Catch[[y]]), al.out$Catch[[y]], rep(999,7))
    ALtp$agecomp<-as.data.frame(rbind(ALtp$agecomp,as.vector(newagecomp)))
  }
  ### Survey
  if (y==1){
    firstsurveyagecomp<-c(1,1,2,0,0,1,-1,-1,sum(al.out$NumbersatAge[[1]][2:7]), c(0,al.out$NumbersatAge[[1]][2:7]), rep(999,7))
    ALtp$agecomp<-as.data.frame(rbind(ALtp$agecomp,as.vector(firstsurveyagecomp)))
  } else{
    newsurveyagecomp<-c(y,1,2,0,0,1,-1,-1,sum(al.out$NumbersatAge[[y]][2:7]), c(0,al.out$NumbersatAge[[y]][2:7]), rep(999,7))
    ALtp$agecomp<-as.data.frame(rbind(ALtp$agecomp,as.vector(newsurveyagecomp)))
  }
  

  # only write data file for years after burn in 
if (y==30|y==33|y==36|y==39|y==42|y==45|y==48|y==51|y==54|y==57){
# write list to .dat
SS_writedat(ALtp, "C:/Users/w10007346/Documents/Assesment model/Alabama/ALss3.dat", overwrite=T)
  }






#################### Assessment model ##############################


if (y==30|y==33|y==36|y==39|y==42|y==45|y==48|y==51|y==54|y==57){   # only years after 25, read in all the data and conduct the assessment 

  
  # Mississippi
  setwd("C:/Users/w10007346/Documents/Assesment model/Mississippi")
  cmd_name<- "MSss3.exe"
  system2(cmd_name)
  
  MSassess<-SS_output(dir="C:/Users/w10007346/Documents/Assesment model/Mississippi", forecast=F,covar=F)
  MSSPR<-MSassess$sprseries$SPR
  
  file.remove("MSss3.dat")  #remove the data file so it can be rewritten next year without errors 
  
  # Louisiana
  # Run SS3 on this years data
  setwd("C:/Users/w10007346/Documents/Assesment model/Louisiana")
  cmd_name<- "LAss3.exe"
  system2(cmd_name)
  
  LAassess<-SS_output(dir="C:/Users/w10007346/Documents/Assesment model/Louisiana", forecast=F,covar=F)
  LASPR<-LAassess$sprseries$SPR
  
  file.remove("LAss3.dat")  #remove the data file so it can be rewritten next year without errors 
  
  #Alabama
  setwd("C:/Users/w10007346/Documents/Assesment model/Alabama")
  cmd_name<- "ALss3.exe"
  system2(cmd_name)
  
  ALassess<-SS_output(dir="C:/Users/w10007346/Documents/Assesment model/Alabama", forecast=F,covar=F)
  ALSPR<-ALassess$sprseries$SPR
  
  file.remove("ALss3.dat")  #remove the data file so it can be rewritten next year without errors 

  al.gradient[y]<-ALassess$maximum_gradient_component
  ms.gradient[y]<-MSassess$maximum_gradient_component
  la.gradient[y]<-LAassess$maximum_gradient_component
  
  }  # end of if statement 
  


#################### Regulation Determination ##########################################
  # only for years after assessments
if (y==31|y==34|y==37|y==40|y==43|y==46|y==49|y==52|y==55|y==58){
ALT5<-mean(tail(ALSPR,1))
LAT5<-mean(tail(LASPR,1))  
MST5<-mean(tail(MSSPR,1))

# set up to determine if SPR is increasing as secondary qualification of regulation change
al.tails<-tail(ALSPR,5)
ms.tails<-tail(MSSPR,5)
la.tails<-tail(LASPR,5)
index<-seq(1,5,1)
al.inc<-summary(lm(al.tails~index))$coefficients[2] # coefficient of increase/decrease of spr over 5 years
ms.inc<-summary(lm(ms.tails~index))$coefficients[2]
la.inc<-summary(lm(la.tails~index))$coefficients[2]

# based on the target SPR plus or minus 5 percent, change MLL up or down one
# second if statement keeps regulations restricted to the provided options, so the index doesnt go to negative 
# soft code the SPR thresholds later
# if SPR is decreasing by less than a substantial amount (.005) or increasing, change regulations
if (ALT5>.35){if(al.inc>-.005){ if (which(sapply(Lengthlimit, function(x) all(al.new.MLL %in% x))==TRUE)[1]!=1) {al.new.MLL<-unlist(Lengthlimit[which(sapply(Lengthlimit, function(x) all(al.new.MLL %in% x))==TRUE)[1]-1])} else {al.new.MLL=MLL.11}}}
# if SPR isnt increasing by a substantial amount (.005) change the regulations
if (ALT5<.25){if(al.inc<.005){ if (which(sapply(Lengthlimit, function(x) all(al.new.MLL %in% x))==TRUE)[1]!=6){al.new.MLL<-unlist(Lengthlimit[which(sapply(Lengthlimit, function(x) all(al.new.MLL %in% x))==TRUE)[1]+1])} else{al.new.MLL=MLL.16}}}

if (MST5>.25){if(ms.inc>-.005){  if (which(sapply(Lengthlimit, function(x) all(ms.new.MLL %in% x))==TRUE)[1]!=1){ms.new.MLL<-unlist(Lengthlimit[which(sapply(Lengthlimit, function(x) all(ms.new.MLL %in% x))==TRUE)[1]-1])} else{ms.new.MLL=MLL.11}}}
if (MST5<.15){ if(ms.inc<.005){ if (which(sapply(Lengthlimit, function(x) all(ms.new.MLL %in% x))==TRUE)[1]!=6){ms.new.MLL<-unlist(Lengthlimit[which(sapply(Lengthlimit, function(x) all(ms.new.MLL %in% x))==TRUE)[1]+1])} else{ms.new.MLL=MLL.16}}}

if (LAT5>.23){ if(la.inc>-.005){if (which(sapply(Lengthlimit, function(x) all(la.new.MLL %in% x))==TRUE)[1]!=1){la.new.MLL<-unlist(Lengthlimit[which(sapply(Lengthlimit, function(x) all(la.new.MLL %in% x))==TRUE)[1]-1])} else{la.new.MLL=MLL.11}}}
if (LAT5<.13){ if(la.inc<.005){ if (which(sapply(Lengthlimit, function(x) all(la.new.MLL %in% x))==TRUE)[1]!=6){la.new.MLL<-unlist(Lengthlimit[which(sapply(Lengthlimit, function(x) all(la.new.MLL %in% x))==TRUE)[1]+1])} else{la.new.MLL=MLL.16}}}


}
  # track regulation changes through the projection window
  al.regulations[y]<-which(sapply(Lengthlimit, function(x) all(al.new.MLL %in% x))==TRUE)
  ms.regulations[y]<-which(sapply(Lengthlimit, function(x) all(ms.new.MLL %in% x))==TRUE)
  la.regulations[y]<-which(sapply(Lengthlimit, function(x) all(la.new.MLL %in% x))==TRUE)[1]
  
  } # end of for loop for projection



############ Assessment model diagnostics ###############################
# ms.regulations
# al.regulations
# la.regulations
# 
# al.Rel.regulations
# ms.Rel.regulations
# la.Rel.regulations
# 
# ms.gradient
# al.gradient
# la.gradient



# # Alabama
# SSplotCatch(ALassess, subplots=1:15, print=T)
# SSplotIndices(ALassess, subplots=1:9, print=T)
# SSplotComps(ALassess, kind="AGE", print=T)
# # 
# # # Mississippi
# SSplotCatch(MSassess, subplots=1:15, print=T)
# SSplotIndices(MSassess, subplots=1:9, print=T)
# SSplotComps(MSassess, kind="AGE", print=T)
# # 
# # # Louisiana
# SSplotCatch(LAassess, subplots=1:15, print=T)
# # SSplotIndices(LAassess, subplots=1:9, print=T)
# # SSplotComps(LAassess, kind="AGE", print=T)

################ MSE output ###############################
la.annual.harvest<-lapply(la.out$Catch, sum)
ms.annual.harvest<-lapply(ms.out$Catch, sum)
al.annual.harvest<-lapply(al.out$Catch, sum)

master.output$La.harvest[[s]]<-unlist(la.annual.harvest)
master.output$Al.harvest[[s]]<-unlist(al.annual.harvest)
master.output$Ms.harvest[[s]]<-unlist(ms.annual.harvest)

# Recruitment of each year for each simulation, to be pulled in for subsequent treatment runs
master.output$`La.recruits`[[s]]<-unlist(lapply(la.out$NumbersatAge, `[[`,1))
master.output$`Ms.recruits`[[s]]<-unlist(lapply(ms.out$NumbersatAge, `[[`,1))
master.output$`Al.recruits`[[s]]<-unlist(lapply(al.out$NumbersatAge, `[[`,1))


# Spawning biomass 
maturity<-c(0,.8,1,1,1,1,1)

## SPB multiplier
# Median length for each age in cm
L<-0
for (i in 1:6){
  linf<-604.52   # Parameter values from Dippold thesis 
  a<-1.74
  b<-.54
  L[i]<-linf/(1+a*(exp(-b*i)))
}

Eggsperindivid<--3500000+147000*L
Eggsperindivid<-c(0,Eggsperindivid)


la.num.age<-lapply(la.out$NumbersatAge, `[`,1:7) # get only ages 0-6
Spawners.la<-lapply(la.num.age, "*", maturity) # how many in each age are mature
la.SPBy<-lapply(Spawners.la, "*", Eggsperindivid)
la.SPBy.sum<-lapply(la.SPBy, sum)
master.output$`La.SSB`[[s]]<-unlist(la.SPBy.sum)

ms.num.age<-lapply(ms.out$NumbersatAge, `[`,1:7) # get only ages 0-6
Spawners.ms<-lapply(ms.num.age, "*", maturity) # how many in each age are mature
ms.SPBy<-lapply(Spawners.ms, "*", Eggsperindivid)
ms.SPBy.sum<-lapply(ms.SPBy, sum)
master.output$`Ms.SSB`[[s]]<-unlist(ms.SPBy.sum)

al.num.age<-lapply(al.out$NumbersatAge, `[`,1:7) # get only ages 0-6
Spawners.al<-lapply(al.num.age, "*", maturity) # how many in each age are mature
al.SPBy<-lapply(Spawners.al, "*", Eggsperindivid)
al.SPBy.sum<-lapply(al.SPBy, sum)
master.output$`Al.SSB`[[s]]<-unlist(al.SPBy.sum)

# plot(unlist(la.SPBy.sum), ylim=c(100000,1700000000000))
# points(unlist(al.SPBy.sum))
# points(unlist(ms.SPBy.sum))

# Relative stock density 
L # I dont have a age length key that includes age zeros, and no one cares too much 
# about specifically RSD, so I will measure the proportion of the population older than 3
la.natage.sum<-lapply(la.num.age, sum) # total number of fish per year
la.old.fish<-lapply(la.out$NumbersatAge, `[`,4:7)
la.old.fish.sum<-lapply(la.old.fish, sum) # total number of age 3 and older individuals
la.rsd<-unlist(la.old.fish.sum)/unlist(la.natage.sum) #percent age 3 or older in the population
master.output$`La.RSD`[[s]]<-la.rsd

ms.natage.sum<-lapply(ms.num.age, sum) # total number of fish per year
ms.old.fish<-lapply(ms.out$NumbersatAge, `[`,4:7)
ms.old.fish.sum<-lapply(ms.old.fish, sum) # total number of age 3 and older individuals
ms.rsd<-unlist(ms.old.fish.sum)/unlist(ms.natage.sum) #percent age 3 or older in the population
master.output$`Ms.RSD`[[s]]<-ms.rsd

al.natage.sum<-lapply(al.num.age, sum) # total number of fish per year
al.old.fish<-lapply(al.out$NumbersatAge, `[`,4:7)
al.old.fish.sum<-lapply(al.old.fish, sum) # total number of age 3 and older individuals
al.rsd<-unlist(al.old.fish.sum)/unlist(al.natage.sum) #percent age 3 or older in the population
master.output$`Al.RSD`[[s]]<-al.rsd

# store regulation changes for each simulation
master.output$Al.regulations[[s]]<-al.regulations
master.output$Ms.regulations[[s]]<-ms.regulations
master.output$La.regulations[[s]]<-la.regulations

# store gradients of each stock assessment
master.output$Al.gradient[[s]]<-al.gradient
master.output$Ms.gradient[[s]]<-ms.gradient
master.output$La.gradient[[s]]<-la.gradient

# store SPR estimates for each state 
master.output$Al.SPR[[s]]<-ALSPR
master.output$Ms.SPR[[s]]<-MSSPR
master.output$La.SPR[[s]]<-LASPR

# store tail estimates of SPR for checking regulation changes
master.output$Al.SPRTAIL[[s]]<-ALT5
master.output$Ms.SPRTAIL[[s]]<-MST5
master.output$La.SPRTAIL[[s]]<-LAT5

} # end of for loop for each simulation

# timer end
toc(quiet=F)



# Save output
getwd()
setwd("C:/Users/w10007346/Documents/MSE")
# annotate datafile for 
master.output$comments<-"assessment every 3 years, regulation changes if terminal year SPR is diff by .05, 
initial F values are fixed at correct values, 30 year burn in"
master.output$enhance.details<-"zero stockers, 9% discard, base SPR"

##
# hold the phone
##
save(master.output, file="Scenario10.master.output.v24.RData")


# diagnose regulation determination, how often does it change, does it meet spr target values
listss<-list(list())
for ( i in 1:100){
  listss[[i]]<-master.output$La.regulations[[i]][30:60]
}
listss


master.output$Ms.SPR[[5]]
master.output$Al.SPR[[1]]
master.output$Al.regulations[[100]]
master.output$Ms.SPR[[100]]

ms.means<-lapply(master.output$Ms.SPR, mean)[1:10]
length(which(ms.means>.25))/100
length(which(ms.means<.15))/100 # small percentage of runs had a mean spr more than .05 diff than target value

la.means<-lapply(master.output$La.SPR, mean)[1:10]
length(which(la.means>.23))/100
length(which(la.means<.13))/100  # louisiana very stable

al.means<-lapply(master.output$Al.SPR, mean)[1:10]
length(which(al.means>.35))/100
length(which(al.means<.25))/100

