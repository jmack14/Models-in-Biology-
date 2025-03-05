require(deSolve)

#Define the function corresponding to the differential equation to be solved
exp<-function(t,y,params){
  S<-y[1]
  E<-y[2]
  Q<-y[3]
  I<-y[4]
  J<-y[5]
  R<-y[6]
  with(as.list(params),{
    dS<-Pi-(S[beta*I+eJ*beta*J]/N)-mu*S
    dE<-(S[beta*I+eJ*beta*J]/N)-(y1+k1+mu)*E
    dQ<-y1*E-(k2+mu)*Q
    dI<-k1*E-(y2+d1+sigma1+mu)*I
    dJ<-y2*I+k2*Q-(sigma2+d2+mu)*J
    dR<-sigma1*I+sigma2*J-mu*R
    res<-c(dS,dE,dQ,dI,dJ,dR)
    list(res)
  })
}

# The initial conditions
S0<- 4000000
E0<-6
I0<-1
Q0=0
R0=0
J0=0


# The parameters for the model
Pi<-136
beta<-0.2
mu<-0.000034
eJ<-0.36
y1<-0.1
y2<-0.5
k1<-0.1
k2<-0.125
d1<-0.0079
d2<-0.0068
sigma1<-0.0377
sigma2<-0.0386
params<-c(Pi=Pi)
N=4000007



# The time values that you want to solve the differential equation for
times <- seq(0,10,0.01)

# The lsoda command solves the differential equation and saves the results in a data frame
nsol <- as.data.frame(lsoda(y=c(S=S0,E=E0,Q=Q0,I=I0,J=J0,R=R0), times, exp, params))



