#Figure 2 : Dissimilarity-Overlap Curve (DOC)
#Note: Below functions are adapted from the DOC package to fit linear model rather than linear mixed-effects models

#The general function
my.DOC<-function (otu, R = 100, subr = NULL, pair = NULL, mov.avg = 5, 
                  ci = c(0.025, 0.5, 0.975), span = 0.2, degree = 1, family = "symmetric", 
                  iterations = 4, surface = "interpolate", cores = 1) 
{
  if (!is.null(pair)) {
    if (length(pair) != 2) 
      stop("There should only be two names in pair!")
  }
  otun <- apply(otu, 2, function(x) as.numeric(x)/sum(as.numeric(x)))
  if (!is.null(pair)) {
    Dis.Over <- DOC.do(otun, pair = pair)
  }
  else {
    Dis.Over <- DOC.do(otun)
  }
  Bootstrap <- my_DOC_boot(Dis.Over, R = R, subr = subr, pair = pair, 
                           mov.avg = mov.avg, span = span, degree = degree, family = family, 
                           iterations = iterations, surface = surface, cores = cores)
  LCIS <- DOC.ci(Bootstrap, ci = ci)
  LOWESS <- DOC.loess(Dis.Over, pair = pair, span = span, degree = degree, 
                      family = family, iterations = iterations, surface = surface)
  LME <- as.data.frame(Bootstrap[[2]])
  colnames(LME) <- "Slope"
  NEG <- as.data.frame(Bootstrap[[3]])
  colnames(NEG) <- "Neg.Slope"
  FNS <- as.data.frame(Bootstrap[[4]])
  colnames(FNS) <- "Fns"
  Final <- list(DO = Dis.Over[[3]], LME = LME, LOWESS = LOWESS, 
                NEG = NEG, FNS = FNS, BOOT = Bootstrap[[1]], CI = LCIS)
  Final <- lapply(Final, function(k) data.frame(k, Data = "Real"))
  class(Final) <- "DOC"
  return(Final)
}

#The bootstrap function
my_DOC_boot<- function (do, R = 100, subr = NULL, pair = NULL, mov.avg = 5, 
                        span = 0.2, degree = 1, family = "symmetric", iterations = 4, 
                        surface = "interpolate", cores = 1) 
{
  #由于不在函数空间内,一些包需要提前加载
  require(foreach)
  require(lme4)
  require(snow)
  require(doSNOW)
  OL <- do[[1]]
  DIS <- do[[2]]
  rownames(OL) <- 1:nrow(OL)
  rownames(DIS) <- 1:nrow(DIS)
  colnames(OL) <- 1:ncol(OL)
  colnames(DIS) <- 1:ncol(DIS)
  xs <- seq(0, 1, by = 0.001)
  if (cores == 1) {
    registerDoSEQ()
  }
  else {
    cl <- makeCluster(cores)
    registerDoSNOW(cl)
  }
  message("Running bootstraps")
  pb <- txtProgressBar(max = R, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  i <- NULL
  llboot <- foreach(i = 1:R, .options.snow = opts, .export = "mov.avg") %dopar% 
    {
      if (is.null(pair)) {
        if (is.null(subr)) {
          Samp <- sample(rownames(OL), replace = TRUE)
        }
        else {
          Samp <- sample(rownames(OL), subr, replace = FALSE)
        }
        OL.sub <- OL[Samp, Samp]
        DIS.sub <- DIS[Samp, Samp]
        OL.tri <- OL.sub[upper.tri(OL.sub)]
        DIS.tri <- DIS.sub[upper.tri(DIS.sub)]
      }
      else {
        if (is.null(subr)) {
          Sampr <- sample(rownames(OL), replace = TRUE)
          Sampc <- sample(colnames(OL), replace = TRUE)
        }
        else {
          Sampr <- sample(rownames(OL), subr, replace = FALSE)
          Sampc <- sample(colnames(OL), subr, replace = FALSE)
        }
        OL.sub <- OL[Sampr, Sampc]
        DIS.sub <- DIS[Sampr, Sampc]
        OL.tri <- as.numeric(OL.sub)
        DIS.tri <- as.numeric(DIS.sub)
      }
      DF.l <- data.frame(y = c(DIS.tri), x = c(OL.tri))
      LOW <- loess(y ~ x, data = DF.l, span = span, degree = degree, 
                   family = family, iterations = iterations, surface = surface)
      LOW.P <- data.frame(predict(LOW, xs))
      colnames(LOW.P) <- paste("rJSD Boot", i)
      ####################
      rowCol <- expand.grid(rownames(OL.sub), colnames(OL.sub))
      if(is.null(pair)){
        labs <- rowCol[as.vector(upper.tri(OL.sub, diag = F)), 
        ]
        tris <- cbind(labs, OL.sub[upper.tri(OL.sub, diag = F)], 
                      DIS.sub[upper.tri(DIS.sub, diag = F)])
        colnames(tris) <- c("Row", "Col", "OL", "DIS")
        Tris <- as.data.frame(na.omit(tris))
        Tris.sub <- Tris[Tris$OL >= median(OL.sub, na.rm = T),]
        fit <- lm(DIS ~ OL, data = Tris.sub) # Using lm fit rather than the original lmer fit
        Est <- coef(summary(fit))[, "Estimate"][2]
      }
      else{
        tris <- cbind(rowCol,OL.tri,DIS.tri)
        colnames(tris) <- c('Row','Col','OL','DIS')
        Tris <- as.data.frame(na.omit(tris))
        Tris.sub <- Tris[Tris$OL >= median(OL.sub, na.rm =T),]
        fit <- lm(DIS ~ OL, data = Tris.sub)
        Est <- coef(summary(fit))[,'Estimate'][2]
      }
      ######################
      ma <- function(x, n = mov.avg) {
        filter(x, rep(1/n, n), sides = 2)
      }
      low.ma <- ma(LOW.P[, 1])
      slope <- diff(low.ma)/diff(xs)
      point <- which(slope > 0)
      neg.slope <- xs[point[length(point)]]
      Fns <- sum(OL.tri > neg.slope, na.rm = T)/length(OL.tri)
      Final <- list(LOW.P, Est, neg.slope, Fns)
      return(Final)
    }
  if (cores > 1) {
    stopCluster(cl)
    registerDoSEQ()
  }
  LOWES <- do.call(cbind, lapply(llboot, function(x) x[[1]]))
  LME <- do.call(cbind, lapply(llboot, function(x) x[[2]]))
  NEG <- do.call(cbind, lapply(llboot, function(x) x[[3]]))
  FNS <- do.call(cbind, lapply(llboot, function(x) x[[4]]))
  LOWES <- cbind(xs, LOWES)
  colnames(LOWES)[1] <- "Overlap"
  Final <- list(LOWES, as.numeric(LME), as.numeric(NEG), as.numeric(FNS))
  return(Final)
}

#DOC analysis
doc_obj <- my.DOC(f_norm %>% t,R=1000,cores = 7) #f_norm->The mycobiome profile with sample as rows and OTU as columns

#Calculate the Nso
nso <- coef(summary(lm(rJSD~Overlap,data=doc_obj$DO[f_plant_doc$DO$Overlap>=median(f_plant_doc$DO$Overlap),])))[,'Estimate'][2]

#Calculate the P value of Nso based on 1000 times bootstrap
sum(doc_obj$LME$Slope>=0)/1000