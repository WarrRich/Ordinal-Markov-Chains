####################
# R Code to accompany the paper entitled "Quantitative Modeling for Threat Analysis using Ordinal Markov Chains"
####################

# Setting a seed for reproducibility
set.seed(1)

####################
# Function definitions

# Function to get a sample from the Dirichlet(alpha) distribution
rdirichlet <- function(alpha) {
  k <- length(alpha)
  samples <- rgamma(k,alpha,1)
  samples/sum(samples)
}

# Function to normalize the row of a matrix (so that they sum to 1)
normalize <- function(m) {
  for (i in 1:nrow(m)) {
    m[i,] <- m[i,]/sum(m[i,])
  }
  m
}

# Function to sample a proability matrix given a mean of "mat" and concentration of "c"
sample.p <- function(mat,c=100) {
  output <- mat
  for (i in 1:nrow(mat)) {
     output[i,] <- rdirichlet(mat[i,]*c)
  }
  output
}

# Function to compute quantities from Section 4.2
# The fundamental matrix (Equation 12), the average time to absorption, 
# and the probability of attaining each absorbing state (Equation 13) are returned in a list
get.process.quants <- function(P) {
  trans.states <- which(diag(P)<1)
  Q <- P[trans.states,trans.states]; R <- P[trans.states,setdiff(1:nrow(P),trans.states)]
  N <- solve(diag(length(trans.states))-Q)
  tmp.names <- rownames(N)
  rownames(N) <- colnames(N)
  colnames(N) <- tmp.names
  T <- apply(N,1,sum)
  B <- N%*%R
  list(Avg.Num.Visits=N,Avg.Time.2.Absorb=T,Prob.Absorb=B)
}



####################
# Computing for Section 3.2

####
# Definitions

# Defining the P_O matrix
P <- matrix(c(0,0.333,0.667,0,0,1,0,0,0.333,0.167,0,0.5,0,0,0,1),ncol=4,byrow=T)

####
# Results Computation

# Obtaining Equation 9
round(P%*%P,3)

# Applying Algorithm 1 (using the infinity norm) to find the Equilibrium Probabilities (P^*_O)
P_0 <- P
while (max(apply(abs(P%*%P_0-P_0),1,sum)) > 10e-6) {
  P_0 <- P%*%P_0
}
round(P_0,3)

####################
# Computing for Section 3.3.1

####
# Definitions

# Defining the P_O matrix
P_N <- matrix(c(0,0.01,0.99,0,0,1,0,0,0.1,0.05,0,0.85,0,0,0,1),ncol=4,byrow=T)

# Algorithm 1 (using the infinity norm) to find the Equilibrium Probabilities (P^*_N)
algo1 <- function(P,e=10e-6) {
  P0 <- P
  while (max(apply(abs(P%*%P0-P0),1,sum)) > e) {
    P0 <- P%*%P0
  }
  P0
}

####
# Results Computation

round(algo1(P_N),3)

####################
# Computing for Section 3.4.1

####
# Definitions
algo3 <- function(n) {
  condition <- FALSE
  while (!condition) {
    x <- rdirichlet(rep(1,n))
    condition <- prod(order(x)==(1:n))
  }
  x
}

####
# Results Computation

# First row of proability transition matrix
row1 <- round(algo3(2),2)
# Third row of proability transition matrix
row3 <- round(algo3(3),2)

P_E <- matrix(NA,nrow=4,ncol=4) 
#Finding the rows of the P_E matrix
P_E[1,] <- c(0,row1,0) # First row
P_E[2,] <- c(0,1,0,0) # Second row
P_E[3,] <- c(row3[2:1],0,row3[3]) # Third row
P_E[4,] <- c(0,0,0,1) # Fourth row

# equilibrium probability matrix
round(algo1(P_E),2)

####################
# Computing for Section 4.1

round(rdirichlet(rep(1/3,3)*10),4)
round(rdirichlet(rep(1/3,3)*10),4)
round(rdirichlet(rep(1/3,3)*10),4)
round(rdirichlet(rep(1/3,3)*10),4)

round(rdirichlet(rep(1/3,3)*100),4)
round(rdirichlet(rep(1/3,3)*100),4)
round(rdirichlet(rep(1/3,3)*100),4)
round(rdirichlet(rep(1/3,3)*100),4)


####################
# Computing for Section 5

####
# Definitions

# Defining the P_D matrix (Equation 14)
P_D <- matrix(c(0,1,0,0,0,0,0,0,0,0,0,0,0.5,0.5,0,0,0,0,0,0,0,0,0,0.675,0.325,0,0,0,0,0,
0.03,0,0.325,0,0.5,0.125,0,0,0,0.03,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0.97,0,0,0.03,
0,0,0,0,0,0,0,0.97,0.03,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1),
  ncol=10,byrow=T)
P_D <- normalize(P_D)

## Concentration of Inf gives mean values of:
get.process.quants(P_D)

# Defining the P_F matrix (Equation 15)
P_F <- matrix(c(0,1,0,0,0,0,0,0,0,0,0,0,0.5,0.5,0,0,0,0,0,0,0,0,0,0.667,0.333,0,0,0,0,0,
0.1,0,0.267,0,0.333,0.2,0,0,0,0.1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0.667,0,0,0.333,
0,0,0,0,0,0,0,0.667,0.333,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1),
ncol=10,byrow=T)
P_F <- normalize(P_F)

## Concentration of Inf gives mean values of:
get.process.quants(P_F)

# Defining function to calculate the absorption into State 8 ("Attack Succeeds") given starting in State 1
prob1to8.samps <- function(c,p,n=10000) {
  output <- rep(NA,n)
  for (i in 1:n) {
    P1 <- sample.p(p,c=c)
    temp.out <- get.process.quants(P1)
    output[i] <- temp.out$Prob.Absorb[1,1]
  }
  output
}

####
# Results Computation

## Reproducing Figure 9a
### Obtaining Monte Carlo samples of the probability of starting in State 1 and being Absorbed in State 8
data1to8 <- matrix(NA,nrow=100,ncol=4); data1to8f <- matrix(NA,nrow=100,ncol=4)
for (i in 1:100) {
  temp <- prob1to8.samps(i*10,P_D,100000)
  data1to8[i,] <- c(quantile(temp,c(.025,.5,.975)),mean(temp))
  temp <- prob1to8.samps(i*10,P_F,100000)
  data1to8f[i,] <- c(quantile(temp,c(.025,.5,.975)),mean(temp))
}
#### Plot 
#pdf("attackSucceedsBW-Both.pdf",width=7,height=6)
plot((1:100)*10,data1to8f[,1],type="n",ylim=c(0,1),xlim=c(0,1000),
     xlab="Concentration",ylab="Probability",main="Probability Attack Succeeds")
lines((1:100)*10,data1to8[,1],col="black",lwd=2,lty=1)
lines((1:100)*10,data1to8[,3],col="black",lwd=2,lty=1);lines((1:100)*10,data1to8[,4],col="gray",lwd=2)
lines((1:100)*10,data1to8f[,1],col="black",lwd=2,lty=2)
lines((1:100)*10,data1to8f[,3],col="black",lwd=2,lty=2);lines((1:100)*10,data1to8f[,4],col="gray",lwd=2,lty=2)
lines(c(500,500),c(0,1),col="grey30",lty=3)
legend("bottomright",legend=c(expression(italic("P")[D] ~ " 95% Prob Interval"), expression(italic("P")[D] ~ 
     " Average"),expression(italic("P")[F] ~ " 95% Prob Interval"), expression(italic("P")[F] ~ " Average")), 
     col=c("black","gray","black","gray"),lwd=c(2,2,2,2),lty=c(1,1,2,2))
#dev.off()

## Reproducing Figure 9b
### Obtaining Monte Carlo samples of the probability of starting in State 1 and being absorbed in State 8
### with a concentration of 500
Our.Method <- prob1to8.samps(500,P_D,100000)
Fish.Method <- prob1to8.samps(500,P_F,100000)
#### Plot 
#pdf("attackSucceedsComp.pdf",width=7,height=4.5)
plot(density(Fish.Method),xlim=c(.2,.9),xlab='Probability Attack Succeeds',ylab='Density',ylim=c(0,20),main='',lwd=2,lty=2,type = "n")
abline(h=0,col="white")
par(new=T)
hist(Fish.Method,xlim=c(.2,.9),ylim=c(0,20),freq=F,ylab='',xlab='',main='')
par(new=T)
plot(density(Fish.Method),xlim=c(.2,.9),xlab='Probability Attack Succeeds',ylab='Density',ylim=c(0,20),main='',lwd=2,lty=2)
par(new=T)
hist(Our.Method,xlim=c(.2,.9),ylim=c(0,20),freq=F,ylab='',xlab='',main='')
par(new=T)
plot(density(Our.Method),xlim=c(.2,.9),xlab='',ylab='',ylim=c(0,20),axes=F,main='',lwd=2,lty=1)
legend("topright",legend=c(expression(italic("P")[D]),expression(italic("P")[F])),lwd=c(2,2),lty=c(1,2))
abline(h=0,col="grey30",lty=3)
#dev.off()


