#
require(deSolve)

# This defines the differential equations
rhs <- function(t,x,parms){
# t is time, x is the vector of state variables, parms is an annoying
# vector of parameter that I ignore
  S1 <- x["S1"]  
  I1 <- x["I1"] 
  R1 <- x["R1"] 
  S2 <- x["S2"]  
  I2 <- x["I2"] 
  R2 <- x["R2"]
  S3 <- x["S3"]
  I3 <- x["I3"]
  R3 <- x["R3"]
  F1 <- v1^2*theta1*zeta1*I1 +  v1*v2*theta2*zeta1*I2 + v1*v3*theta3*zeta1*I3
  F2 <- v1*v2*theta1*zeta2*I1 + v2^2*theta2*zeta2*I2 + v2*v3*theta3*zeta2*I3
  F3 <- v1*v3*theta1*zeta3*I1 + v2*v3+theta2*zeta3*I2 + v3^2*theta3*zeta3*I3
  dS1 <- -F1*S1
  dI1 <- F1*S1-gamma*I1
  dR1 <- gamma*I1
  dS2 <- -F2*S2
  dI2 <- F2*S2-gamma*I2
  dR2 <- gamma*I2
  dS3 <- -F3*S3-gamma*I3
  dI3 <- F3*S3-gamma*I3
  dR3 <- gamma*I3
  return(list(c(dS1,dI1,dR1,dS2,dI2,dR2,dS3,dI3,dR3)))  
}

# Define summary file
summ <- data.frame(case=c("baseline"),
   R0=1,lambda=1,finalS1=0,maxI1=0,maxI1time=0,
                 finalS2=0,maxI2=0,maxI2time=0,
                 finalS3=0,maxI3=0,maxI3time=0)
summnms <- c("R0","lambda","finalS1","maxI1","maxI1time","finalS2","maxI2","maxI2time","finalS3","maxI3","maxI3time")

# Define times to solve
tmax <- 365
times <- seq(from=0,to=tmax,by=1)



# Parameters used by all cases
#N <- 3.27e8 #population of the US
N1 <- 500000
N2 <- 250000
N3 <- 150000
gamma <- 1/14 # Recovery rate

# Define initial condition, one infected person in each group
init <- c(N1-1,1,0,N2-1,1,0,N3-1,1,0)
names(init) <- c("S1","I1","R1","S2","I2","R2","S3","I3","R3")

# Baseline case: Everyone the same
theta1 <- 1
zeta1 <- 1
v1 <- 2.34e-5
theta2 <- 1
zeta2 <- 1
v2 <- 2.34e-5
theta3 <- 1
zeta3 <- 1
v3 <- 2.34e-5

R0 <- v1^2*theta1*zeta1*N1/gamma + v2^2*theta2*zeta2*N2/gamma + v3^2*theta3*zeta3*N3/gamma
lambda <- v1^2*theta1*zeta1*N1 + v2^2*theta2*zeta2*N2-gamma + v3^2*theta3*zeta3*N3/gamma
SIRout <- as.data.frame(ode(y = init,times = times, func = rhs))
finalS1 <- SIRout$S1[nrow(SIRout)]/N1
maxI1 <- max(SIRout$I1)/N1
maxI1time <- SIRout$time[which.max(SIRout$I1)]
finalS2 <- SIRout$S2[nrow(SIRout)]/N2
maxI2 <- max(SIRout$I2)/N2
maxI2time <- SIRout$time[which.max(SIRout$I2)]
finalS3 <- SIRout$S3[nrow(SIRout)]/N3
maxI3 <- max(SIRout$I3)/N3
maxI3time <- SIRout$time[which.max(SIRout$I3)]
summ[summ$case=="baseline",summnms] <- 
    c(R0,lambda,finalS1,maxI1,maxI1time,finalS2,maxI2,maxI2time,finalS3,maxI3,maxI3time)
baseline <- SIRout


plot(I1 ~ time,baseline)
plot(I2 ~ time,baseline)





print(summ)

# Graph it
pdf("I1.pdf")
par(mfrow=c(2,2),mar=c(5,6,4,2))
# Baseline graph
plot(I1 ~ time,baseline)
plot(I2 ~ time,baseline)

# highv graph
plot(I(I1/N1) ~ time,highv,type="l",lwd=2,ylim=c(0,0.4),
    ylab="Fraction infected",
     cex.axis=1.7,cex.lab=1.7,col="black")
lines(I(I2/N2) ~ time,highv,lwd=2,col="gray")
lines(I(I1/N1) ~ time,baseline,lwd=1,lty=2,col="blue")
legend("topright",c("Category 1","Category 2","baseline"),cex=0.9,
      col=c("black","gray","blue"),lty=c(1,1,2),lwd=c(2,2,1))
title(main="Increased velocity")

# hightheta graph
plot(I(I1/N1) ~ time,hightheta,type="l",lwd=2,ylim=c(0,0.4),
    ylab="Fraction infected",
     cex.axis=1.7,cex.lab=1.7,col="black")
lines(I(I2/N2) ~ time,hightheta,lwd=2,col="gray")
lines(I(I1/N1) ~ time,baseline,lwd=1,lty=2,col="blue")
title(main="Increased infectiousness")

#  highzeta graph
plot(I(I1/N1) ~ time,highzeta,type="l",lwd=2,ylim=c(0,0.4),
    ylab="Fraction infected",
     cex.axis=1.7,cex.lab=1.7,col="black")
lines(I(I2/N2) ~ time,highzeta,lwd=2,col="gray")
lines(I(I1/N1) ~ time,baseline,lwd=1,lty=2,col="blue")
title(main="Increased susceptibility")
dev.off()

pdf("S1.pdf")
par(mfrow=c(2,2),mar=c(5,6,4,2))
# Baseline graph
plot(I(S1/N1) ~ time,baseline,type="l",lwd=2,ylim=c(0,1),
    ylab="Fraction susceptible",
     cex.axis=1.7,cex.lab=1.7,col="blue")
title(main="Baseline case")

# highv graph
plot(I(S1/N1) ~ time,highv,type="l",lwd=2,ylim=c(0,1),
    ylab="Fraction susceptible",
     cex.axis=1.7,cex.lab=1.7,col="black")
lines(I(S2/N2) ~ time,highv,lwd=2,col="gray")
lines(I(S1/N1) ~ time,baseline,lwd=1,lty=2,col="blue")
legend("topright",c("Category 1","Category 2","baseline"),cex=0.9,
      col=c("black","gray","blue"),lty=c(1,1,2),lwd=c(2,2,1))
title(main="Increased velocity")

# hightheta graph
plot(I(S1/N1) ~ time,hightheta,type="l",lwd=2,ylim=c(0,1),
    ylab="Fraction susceptible",
     cex.axis=1.7,cex.lab=1.7,col="black")
lines(I(S2/N2) ~ time,hightheta,lwd=2,col="gray")
lines(I(S1/N1) ~ time,baseline,lwd=1,lty=2,col="blue")
title(main="Increased infectiousness")

#  highzeta graph
plot(I(S1/N1) ~ time,highzeta,type="l",lwd=2,ylim=c(0,1),
    ylab="Fraction susceptible",
     cex.axis=1.7,cex.lab=1.7,col="black")
lines(I(S2/N2) ~ time,highzeta,lwd=2,col="gray")
lines(I(S1/N1) ~ time,baseline,lwd=1,lty=2,col="blue")
title(main="Increased susceptibility")
dev.off()

