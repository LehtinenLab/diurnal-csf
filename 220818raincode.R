require(rain)

ldttr <- read.csv("/Users/pkalugin/Desktop/full_runs/220820ldrainttr.csv",header=FALSE)
ldwt <- read.csv("/Users/pkalugin/Desktop/full_runs/220820ldrainwt.csv",header=FALSE)
ddttr  <- read.csv("/Users/pkalugin/Desktop/full_runs/220820ddrainttr.csv",header=FALSE)
ldttr3000 <- read.csv("/Users/pkalugin/Desktop/full_runs/220820ldrainttr3000.csv",header=FALSE)
ldwt3000 <- read.csv("/Users/pkalugin/Desktop/full_runs/220820ldrainwt3000.csv",header=FALSE)
ddttr3000  <- read.csv("/Users/pkalugin/Desktop/full_runs/220820ddrainttr3000.csv",header=FALSE)
#lddatads30 <- read.csv("/Users/pkalugin/Desktop/full_runs/220818ldrainds30.csv",header=FALSE)
#dddatads30 <- read.csv("/Users/pkalugin/Desktop/full_runs/220818ddrainds30.csv",header=FALSE)
#lddatads300 <- read.csv("/Users/pkalugin/Desktop/full_runs/220818ldrainds300.csv",header=FALSE)
#dddatads300 <- read.csv("/Users/pkalugin/Desktop/full_runs/220818ddrainds300.csv",header=FALSE)
#lddatads1500 <- read.csv("/Users/pkalugin/Desktop/full_runs/220818ldrainds1500.csv",header=FALSE)
#dddatads1500 <- read.csv("/Users/pkalugin/Desktop/full_runs/220818ddrainds1500.csv",header=FALSE)

ldttrd <- data.matrix(ldttr, rownames.force = NA)
ldwtd <- data.matrix(ldwt, rownames.force = NA)
ddttrd <- data.matrix(ddttr, rownames.force = NA)
ldttr3000d <- data.matrix(ldttr3000, rownames.force = NA)
ldwt3000d <- data.matrix(ldwt3000, rownames.force = NA)
ddttr3000d <- data.matrix(ddttr3000, rownames.force = NA)

#ldrainds1500 <- rain(lddatads1500d[,2:4],period=24,deltat=1500*lddatad[2,1], peak.border=c(0.1,0.9), verbose=TRUE, method="longitudinal")
#ddrainds1500 <- rain(dddatads1500d[,2],period=24,deltat=1500*lddatad[2,1], peak.border=c(0.1,0.9), verbose=TRUE, method="longitudinal")
ldttr3000drain <- rain(ldttr3000d[,2:11],period=24,deltat=ldttr3000d[2,1], peak.border=c(0.1,0.9), verbose=TRUE, method="longitudinal")
ldwt3000drain <- rain(ldwt3000d[,2:6],period=24,deltat=ldwt3000d[2,1], peak.border=c(0.1,0.9), verbose=TRUE, method="longitudinal")
ddttr3000drain <- rain(ddttr3000d[,2:10],period=24,deltat=ddttr3000d[2,1], peak.border=c(0.1,0.9), verbose=TRUE, method="longitudinal")

ldttr3000drain
ldwt3000drain
ddttr3000drain



set.seed(123)
times <- c(1: 24) * 2
sin <- 1 + 0.5 * sin(times / 24 * 2 * pi) + rnorm(24, 0, 0.3)
saw <- rep(13:24 / 18 , 2) + rnorm(24, 0, 0.3)
measure <- cbind(sin, saw)
require('lattice')
xyplot(t(measure)~rep(times, each=2) | c('sin', 'saw'), layout = c(1, 2), type = 'o', xlab = 'time', ylab = 'value', cex.lab = 0.6)

rainresult <- rain(measure, period=24, deltat=2, peak.border=c(0.1,0.9), verbose=FALSE)