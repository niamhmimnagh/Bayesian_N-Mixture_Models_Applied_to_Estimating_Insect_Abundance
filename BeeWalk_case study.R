## Load Required Packages
library(tidyr)
library(dplyr)
library(R2jags)
library(stringr)


#########################################################################################################
# Data Preparation
data <- readxl::read_excel("bees.xlsx")
data$date<-substr(data$StartDate, 1, 10)      #remove characters from end of date string
data<-data[which(!is.na(data$date)),]         # remove NA values
data<-data[data$Year %in% c(2016,2019),]      # Looking at 2016 and 2019

data<-data[,c(4,5,11,12,15,30, 32, 34, 35, 36, 37,43,44,50)] # Removing columns we are not interested in
data$date<-as.factor(data$date)

data$'Site section'<-str_trim(sub('.*(?=.{3}$)', '', data$`Site section`, perl=T))  # Removing site name at beginning of section name
data$abundance<-rowSums(data[,c(8:11)])                                             # Total abundance rather than abundance divided into queen, worker etc.
data<-data[,-c(8:11)]


species<-c("a white-tailed bumblebee", "Buff-tailed Bumblebee", "Garden Bumblebee", "Tree Bumblebee", "Early Bumblebee",  "Common Carder", "Red-tailed Bumblebee", "Honeybee")
data<-data[(data$species %in% species),]
colnames(data)<-c("site", "section","lat","long","length", "temperature", "species", "Year","Month","date", "abundance")

data<-data[complete.cases(data),]



# Transects have varying numbers of sections. We need all transects examined to have the same number of sections.
# Let T=2 sections for each transect
data$section<-replace(data$section, data$section %in% c("S1","S3", "S5", "S7", "S9", "S11", "S13", "S15"), 1)
data$section<-replace(data$section, data$section %in% c("S2","S4", "S6", "S8", "S10", "S12", "S14", "nds"), 2)
data$section<-as.numeric(data$section)
data<-data[data$Month == 6,]                # Just looking at data collected in June


x<-table(data$site, data$Year)
x<-x[rowSums(x==0)%in% c(0),]
set.seed(20)
sites<-sample(rownames(x),60)
data<-data[data$site %in% sites,] # Taking random sample of 60 sites


# Getting each section ready separately - 
# Section 1
data1<-data[data$section==1,]
data1<-right_join(data1, x, by=c("date", "site"))
data1<-aggregate(abundance~species+site+length+lat+long+temperature+Year, data=data1, FUN=sum)
data1$temperature<-as.numeric(data1$temperature)
data1.1<-aggregate(temperature~site+Year, data=data1, FUN=max)
data2<-left_join(data1, data1.1, by=c("site",  "Year"))
data2.1<-aggregate(abundance~species+site+length+lat+long+temperature.y+Year, data=data2, FUN=sum)


# Section 2
data1<-data[data$section==2,]
data1<-aggregate(abundance~species+site+length+lat+long+temperature+Year, data=data1, FUN=sum)
data1$temperature<-as.numeric(data1$temperature)
data1.1<-aggregate(temperature~site+Year, data=data1, FUN=max)
data1.1<-left_join(data1, data1.1, by=c("site", "Year"))
data2.2<-aggregate(abundance~species+site+length+lat+long+temperature.y+Year, data=data1.1, FUN=sum)

# Putting sections back together
data2<-rbind(data2.1, data2.2)
data2$section<-c(rep(1, nrow(data2.1)), c(rep(2,nrow(data2.2))))


# Adding new rows with 0 abundance where there are missing combinations of species, visit, site and section
a1<-data.frame(table(data2$species, data2$Year, data2$site, data2$section))
colnames(a1)<-c("species",  "Year", "site", "section", "Frequency")

newdata<-as_tibble(data2)
newdata$Year<-as.factor(newdata$Year)
newdata$section<-as.factor(newdata$section)
for(i in 1:nrow(a1)){
  if(a1[i,5]==0){
    newdata<- add_row(newdata,
                      species=a1[i, "species"],
                      site=a1[i,"site"],
                      section=a1[i, "section"],
                      Year=a1[i,"Year"])
  }
}

newdata[,8][is.na(newdata[,8]),] <- 0 # new rows have 0 abundance


# Assigning right values for lat, long, length, date, temp to these new rows
newdata<-newdata[order(newdata$site, newdata$Year),]
a<-newdata %>% fill(Year)

a<-a[order(a$site),]
a<-a %>% fill(lat, .direction="downup")
a<-a %>% fill(long, .direction="downup")
a<-a %>% fill(length, .direction="downup")

a<-a[order(a$site, a$Year),]
a<-a %>% fill(temperature.y, .direction="downup")



# We now have the dataset we require
# We just need to reformat into an array of dimension R,T,S,K
R<-length(unique(a$site)) # number of transects
T<-2                      # number of sections in each transect
S<-length(species)        # number of species
K<-2                      # number of years

a1 <- reshape2::dcast(a, site +Year+species+lat+long+length ~ section, value.var="abundance") # Long to wide format



# Creating an intermediary list "d" so that I can fill in Y[i,t,s,k] more easily
y<-split(a1, a1$Year)
d<-vector(mode="list", length=length(y))

for(i in 1:length(y)){
  d[[i]]<-split(y[[i]], y[[i]]["species"])
}

names(d)<-names(y) # assigning names
Year<-temperature<-list(list(), list(), list(), list()) # initialising dates and temperature variables
lat<-d[[1]][[1]][["lat"]] # latitude of each site
long<-d[[1]][[1]][["long"]] # longitude of each site
length<-d[[1]][[1]][["length"]] # length of each site

for(k in 1:K){  
  Year[[k]]<-d[[k]][[1]][["Year"]] # date and temperature each vary with timepoint
  
  for(s in 1:S){
    d[[k]][[s]]<-as.data.frame(d[[k]][[s]])
    rownames(d[[k]][[s]])<-d[[k]][[s]][,"site"] # assigning site names
    d[[k]][[s]]<-d[[k]][[s]][-c(1:6)] # removing all variables except abundances
  }
}



# Creating Y[i,t,s,k]
# This can be used in the MNM model, or we can select one species and timepoint for the original N-mix model
Y<-array(dim=c(R, T,S,K))

for(i in 1:R){
  for(t in 1:T){
    for(s in 1:S){
      for(k in 1:K){
        Y[i,t,s,k]<-d[[k]][[s]][[t]][i]
      }
    }
  }
}

dimnames(Y)[[1]]<-rownames(d[[1]][[1]])
dimnames(Y)[[2]]<-names(d[[1]][[1]])
dimnames(Y)[[4]]<-names(d)
dimnames(Y)[[3]]<-names(d[[1]])


day<-matrix(ncol=K, nrow=R)
for(k in 1:K){
  day[,k]<-Year[[k]]
}

################################################################################
## JAGS code - MNM model
model_code = '
model {
# Poisson Hurdle Abundance
  for(s in 1:S){
    for (i in 1:R) { 
      for(k in 1:K){
        x[i,s,k] ~ dbern(1-theta)                       # Occupancy component
        log(lambda[i,s,k]) <- a[i,s,k] 
        count[i,s,k] ~ dpois(lambda[i,s,k])T(1,)        # Count component       
        N[i,s,k]<-ifelse(x[i,s,k]==0, 0, count[i,s,k])
      }
    }
  }
  
  for(i in 1:R){
    for(s in 1:S){
      for(k in 1:K){           
        mu[i,s,k]<-beta[1,s]*lat[i]+beta[2,s]*long[i]+beta[3,s]*lat[i]*long[i]+beta[4,s]*day[i,k] 
        for(t in 1:T){   
          Y[i,t,s,k] ~ dbin(probability[i,t,s,k], N[i,s,k])
          logit(probability[i,t,s,k]) <- gamma[i,t,k,s]
          gamma[i,t,k,s]~dunif(-4,-2)
        }
      }
    }
  }
  
  
  # Random effects a with mean vector mu, and variance-covariance matrix cov
  for(i in 1:R){
    for(k in 1:K){
      a[i,1:S,k] ~ dmnorm(mu[i,,k], precision[,])
    }
  }
  
  
  # Wishart prior on precision with df=S+1 and diagonal matrix Omega
  df<-S+1
  precision[1:S,1:S] ~ dwish(Omega[,], df)
  covariance[1:S,1:S]<-inverse(precision[,])
  
  
  # Correlations and Standard deviations
  for (s in 1:S){    
    for (s1 in 1:S){
      cor[s,s1]<-covariance[s,s1]/sqrt(covariance[s,s]*covariance[s1, s1])
    }
  }
  
  for(i in 1:4){
    for(s in 1:S){
      beta[i,s]~dnorm(0,0.01)
    }
  }
  theta~dbeta(1,1)
}
'


# Initial values -
x.inits<-apply(Y, c(1,3,4),function(z) ifelse(any(z>0), 1, 0))
jags.inits <- list(count = apply(Y, c(1,3,4), max)+1, x=x.inits) 

# Parameters to monitor
model_parameters =  c( "probability", "N", "beta", "cor", "theta")

model_data = list("R" = R,  "Y"=Y, "T"=T, "S"=S,"day"=day, "Omega"=diag(S), "lat"=lat, "long"=long, "K"=K)
# Run the model
model_run2 = jags(data = model_data,
                  inits=rep(list(jags.inits), 4),
                  parameters.to.save = model_parameters,
                  model.file = textConnection(model_code),
                  n.chains = 4,
                  n.iter = 20000,
                  n.burnin = 10000,
                  n.thin = 10)
