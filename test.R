# density estimate with default parameters for 
# different sample sizes

x = seq(-3,3,length.out = 1000)
for(n in c(10,100,1000,10000)){
  plot(density(rnorm(n)),
       main=paste0("n=",n),
       xlim=c(-3,3),
       ylim=c(0,0.5))
  lines(x,dnorm(x),col="red")
}
