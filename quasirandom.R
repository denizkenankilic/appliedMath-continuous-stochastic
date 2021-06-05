#Randtoolbox
install.packages("randtoolbox")
library(randtoolbox)
#Halton sequences
# uniform variate
unihalt<-halton(n = 10, dim = 5) #n is number of observations 
#and dim is dimension
# normal variate
normalhalt<-halton(n = 10, dim = 5, normal = TRUE)
#some plots
par(mfrow = c(2, 2), cex = 0.75)
hist(halton(n = 5000, dim = 1), main = "Uniform Halton",
     xlab = "x", col = "steelblue3", border = "white")
hist(halton(n = 5000, dim = 1, norm = TRUE), main = "Normal Halton",
     xlab = "x", col = "steelblue3", border = "white")
#halton plots
halton1<-halton(n=500, dim=2)
halton2<-halton(n=500, dim=3)
halton3<-halton(n=500,dim=101)
par(mfrow = c(2, 2))
plot(halton3[,1],halton3[,2])
plot(halton3[,2],halton3[,3])
plot(halton3[,21],halton3[,101])
plot(halton3[,51],halton3[,101])

#Sobol sequences
sobol1<-sobol(n = 500, dim = 101, scrambling = 3)
par(mfrow = c(2, 2))
plot(sobol1[,1],sobol1[,2])
plot(sobol1[,2],sobol1[,3])
plot(sobol1[,21],sobol1[,101])
plot(sobol1[,51],sobol1[,101])