# Wiener Process
set.seed (123)
N <- 100 # number of end-points of the grid including T
T <- 1 # length of the interval [0,T] in time units
delta <- T/N # time increment
W <- numeric(N+1) # initialization of the vector W
t <- seq (0,T, length =N+1)
 for (i in 2:(N+1))
   W[i] <- W[i-1] + rnorm (1) * sqrt (delta)
   plot (t,W, type ="l", main ="Wiener process",
         ylim =c(-1,1))
         
# Brownian Motion as Random Walk
set.seed(123)
n <- 10 # far from the CLT, number of random variables
T <- 1
t <- seq (0,T,length =100)
#runif generates n random numbers from the uniform
#distribution in (0, 1) and it transforms these into a
#sequence of zeros and ones, then TRUE or FALSE
S <- cumsum (2*( runif (n ) >0.5) -1) #maps 0 to -1 and 1 to 1,calculate S_n W <- sapply (t, function (x) ifelse (n*x >0,S[n*x],0)) #select from S
W <- as.numeric(W)/ sqrt (n)
plot (t,W, type ="l",ylim =c( -1 ,1))
n <- 100 # closer to the CLT
S <- cumsum (2*( runif(n) >0.5) -1)
W <- sapply (t, function (x) ifelse (n*x >0,S[n*x],0))
W <- as.numeric(W)/ sqrt (n)
lines (t,W, lty =2)
n <- 1000 # quite close to the limit
S <- cumsum (2*( runif(n) >0.5) -1)
W <- sapply (t, function (x) ifelse (n*x >0,S[n*x],0))
W <- as.numeric(W)/ sqrt(n)
lines (t,W, lty =3)

# Brownian Motion as L2[0, T] Expansion
set.seed (123)
phi <- function (i,t,T){
(2*sqrt(2*T))/((2 *i +1)*pi)*sin(((2*i+1)*pi*t)/(2*T)) }
T <-
N <-
t <-
W <-
n <-
Z <-
for (i in 2:(N+1))
1 # length of time interval
100 # number of points
seq (0,T,length=N+1) # time increment vector numeric(N+1) # initialization of vector W
10 # number of random variables
rnorm(n) # random numbers distributed standard normal
W[i] <- sum (Z*sapply (1:n, function(x) phi(x,t[i],T))) plot (t,W, type ="l",ylim =c( -1 ,1))
n <- 50
Z <- rnorm (n)
for (i in 2:(N+1))
W[i] <- sum (Z*sapply (1:n, function(x) phi(x,t[i],T))) lines (t,W,lty =2)
n <- 100
Z <- rnorm (n)
for (i in 2:(N+1))
W[i] <- sum (Z*sapply (1:n, function(x) phi(x,t[i],T))) lines (t,W,lty =3)

# Geometric Brownian Motion
set.seed(123)
r <- 1
sigma <-0.5
x <- 10 # start point
N <- 100 # number of end points of the grid including T
T <- 1 # length of the interval [0,T] in time units
Delta <- T/N # time increment
W <- numeric(N+1) # initialization of the vector W
t <- seq(0,T,length=N+1) # time vector
for (i in 2:(N+1))
  W[i] <- W[i-1]+rnorm(1)*sqrt(Delta) # Wiener process
S <- x*exp((r-sigma^2/2)*t+sigma*W)
plot(t,S,type="l",main="Geometric Brownian Motion")

# Brownian Bridge
set.seed (123)
N <- 100 # number of end points of the grid including T
T <- 1 # length of the interval [0 ,T] in time units
Delta <- T/N # time increment
W <- numeric(N+1) # initialization of the vector W
t <- seq(0,T,length=N+1) # time vector
for (i in 2:(N+1))
  W[i] <- W[i-1]+rnorm(1)*sqrt(Delta)
a <- 0 # starting point
b <- -1 #end point
BB <- a+W-t/T*(W[N+1]-b+a) #Brownian bridge
plot (t,BB , type ="l")
abline (h=-1, lty =3)

# SDE Package-BM Function
BM <- function (x=0, t0 =0, T=1, N =100){
  if(T <= t0) stop ("wrong times")
  dt <- (T-t0)/N #time increment
  t <- seq(t0 ,T, length=N+1) #time vector
  X <- ts(cumsum(c(x,rnorm(N)*sqrt(dt))),
          start=t0, deltat=dt)
  return (invisible(X))
}

# SDE Package-BBridge Function
BBridge <- function (x=0, y=0, t0 =0, T=1, N =100){
  if(T <= t0) stop ("wrong times")
  dt <- (T-t0)/N
  t <- seq(t0, T, length=N+1)
  X <- c(0,cumsum(rnorm(N)*sqrt(dt)))
  BB <- x+X-(t-t0)/(T-t0)*(X[N+1]-y+x)
  X <- ts(BB, start=t0 , deltat=dt)
  return(invisible(X))
}

# SDE Package-GBM Function
# x = starting point at time 0
# r = interest rate
# sigma = square root of the volatility
GBM <- function (x, r=0, sigma , T=1, N=100){
  tmp <- BM(T=T,N=N)
  S <- x*exp((r-sigma^2/2)*time(tmp)+sigma*as.numeric(tmp))
  X <- ts(S, start=0, deltat=1/N)
  return(invisible(X))
}
