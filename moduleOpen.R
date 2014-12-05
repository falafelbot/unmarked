# -----------------------------------------------------------------------
# --------------- Fit ppen population models to real data ---------------
# -----------------------------------------------------------------------








# ------------------------ MacKenzie et al. (2003) ----------------------


# PART I. Import data and create unmarkedFrame

alfl.data <- read.csv("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/alfl0506.csv?attredirects=0&d=1", row.names=1)
str(alfl.data)

library(unmarked)
alfl.y <- alfl.data[,1:6]
alfl.y[alfl.y>1] <- 1
alfl.umf <- unmarkedMultFrame(y=alfl.y,
    siteCovs=alfl.data[,c("woody", "struct")],
    obsCovs=list(time=alfl.data[,c("time1_05", "time2_05", "time3_05",
                                   "time2_06", "time3_06", "time3_06")],
                 date=alfl.data[,c("date1_05", "date2_05", "date3_05",
                                   "date2_06", "date3_06", "date3_06")]),
    numPrimary=2)
summary(alfl.umf)

# Standardize covariates
siteCovs(alfl.umf) <- scale(siteCovs(alfl.umf))
obsCovs(alfl.umf) <- scale(obsCovs(alfl.umf))
summary(alfl.umf)



# PART II. Fit models

(fm1 <- colext(~1, ~1, ~1, ~1, alfl.umf))
(fm2 <- colext(~woody, ~1, ~1, ~time+date, alfl.umf))
(fm3 <- colext(~woody+struct, ~1, ~1, ~time+date+struct, alfl.umf))



# PART III. Do stuff with the fitted models.

# Model selection
fms <- fitList("phi(.)col(.)ext(.)p(.)"                 = fm1,
               "phi(woody)col(.)ext(.)p(time+date)"     = fm2,
               "phi(woody+struct)col(.)ext(.)p(time+date+struct)" = fm3)
modSel(fms)

# Back-transform colonization and extinction probabilities
backTransform(fm2, type="col")
backTransform(fm2, type="ext")







# ------------------- Open N-mixture model -------------------------------





# PART I. Import data and create unmarkedFrame


alfl.data <- read.csv("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/alfl0506.csv?attredirects=0&d=1", row.names=1)
str(alfl.data)


library(unmarked)
alfl.umf <- unmarkedFramePCO(y=alfl.data[,1:6],
    siteCovs=alfl.data[,c("woody", "struct")],
# we could have yearlySiteCovs here.
    obsCovs=list(time=alfl.data[,c("time1_05", "time2_05", "time3_05",
                                   "time2_06", "time3_06", "time3_06")],
                 date=alfl.data[,c("date1_05", "date2_05", "date3_05",
                                   "date2_06", "date3_06", "date3_06")]),
    numPrimary=2)

# Standardize covariates after making the UMF
siteCovs(alfl.umf) <- scale(siteCovs(alfl.umf))
obsCovs(alfl.umf) <- scale(obsCovs(alfl.umf))
summary(alfl.umf)



# PART II. Fit models

# It would be better to model omega/gamma too. Not done here b/c it's slow
(fm1 <-  pcountOpen(~1, ~1, ~1, ~1, alfl.umf, K=30))
(fm2 <- pcountOpen(~woody, ~1, ~1, ~date+time, alfl.umf, K=30))
(fm3 <- pcountOpen(~woody+struct, ~1, ~1, ~date+time+struct, alfl.umf,
                   K=30))
(fm4 <- pcountOpen(~woody, ~1, ~1, ~date+time, alfl.umf, K=30,
                   dynamics="autoreg"))


# PART III. Do stuff with fitted models.

# Model selection
fms <- fitList("lam(.)gam(.)om(.)p(.)"                 = fm1,
               "lam(woody)gam(.)om(.)p(date+time)"     = fm2,
               "lam(woody+struct)gam(.)om(.)p(date+time+struct)" = fm3,
               "lam(woody)gam(.)om(.)p(date+time)AR"   = fm4)
modSel(fms)


# Recruitment rate and apparent survival
backTransform(fm2, type="gamma")
backTransform(fm2, type="omega")




# Estimate population size for both years
names(fm2)
N.hat <- function(fm) {
    lambda <- predict(fm, type="lambda")[,1] # abundance at each site
    gamma <- predict(fm, type="gamma")[,1]   # recruitment rate
    omega <- predict(fm, type="omega")[,1]   # apparent survival rate
    N <- matrix(NA, length(lambda), 2)       # expected population size
    N[,1] <- lambda
    N[,2] <- lambda*omega + gamma
    N.total <- colSums(N)
    return(N.total)
}

N.hat(fm2)


# Get confidence intervals for our estimates of these derived parameters
N.hats <- parboot(fm2, N.hat, nsim=20, report=1) # nsim should be >>100
N.hats












# ------------------- Temporary emigration N-mixture model ---------------

# PART I. Import data, create unmarkedFrame

# This is removal sampling data, also known as time to detection
# Each point count was divided into 2 5-min periods.
# Only new detections are recorded
alfl05.remData <- read.csv("http://sites.google.com/site/unmarkedinfo/home/webinars/2012-january/data/alfl05.remData.csv?attredirects=0&d=1",
                           row.names=1)
head(alfl05.remData)

# Notice how the counts diminish during each primary period
colSums(alfl05.remData[,1:6])

# Create a "visit" covariate
visit <- matrix(as.character(1:3), nrow(alfl05.remData), 3, byrow=TRUE)

# Create unmarkedFrame. NOTE: type="removal"
umf <- unmarkedFrameGMM(y=alfl05.remData[,1:6],
    siteCovs=alfl05.remData[,c("woody", "struct")],
    yearlySiteCovs=list(
        visit=visit,
        date=alfl05.remData[,c("date.1","date.2","date.3")],
        time=alfl05.remData[,c("time.1","time.2","time.3")]),
    type="removal", numPrimary=3)

# Standardize covariates (only the continuous ones)
siteCovs(umf) <- scale(siteCovs(umf))
ysc <- yearlySiteCovs(umf)
yearlySiteCovs(umf) <- cbind(visit=ysc[,1], scale(ysc[,2:3]))
summary(umf)

# Fit models
(fm1 <- gmultmix(~1, ~1, ~1, umf, K=40))
(fm2 <- gmultmix(~1, ~1, ~date+time, umf, K=40))
(fm3 <- gmultmix(~woody, ~1, ~date+time, umf, K=40))

# Use the parametric bootstrap to assess model fit
set.seed(3453)
pb.fit <- parboot(fm3, statistic=SSE, nsim=50, report=10)
pb.fit # fail to reject null hypothesis (no overdispersion)

# Compute density in sampled area
D.fn <- function(fm) {
    E.M <- predict(fm, type="lambda")[,1]
    E.phi <- predict(fm, type="phi")[,1]
    plot.area <- pi*50^2 / 1e4 # hectares
    nSites <- length(E.M)
    E.D <- sum(E.M * E.phi) / (plot.area*nSites)
    return(E.D)
}
D.fn(fm3) # 3.59 ALFL/hectare

# Use the parametric bootstrap to compute confidence interval for density
pb.D <- parboot(fm3, statistic=D.fn, nsim=50, report=5)
pb.D # 95% CI = (2.4, 4.7)
plot(pb.D)
