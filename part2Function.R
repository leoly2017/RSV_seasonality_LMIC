library(ggthemes)
# Prepare MAP data
WorldData <- map_data('world')
WorldData %>% filter(region != "Antarctica") -> WorldData
WorldData <- fortify(WorldData)
efficacyGenerater <- function(x, rate, n) {
  x0 <- 3 * x /(rate*rate +rate +1)
  return(x0 * c(rate^(0:(n-1)), rep(0,12-n)))
}
getNEpi <- function(dataset = Seasonality, month.window = 5, cut.off = 2) {
  epiMatrix <- as.matrix(dataset[paste(month.abb, "epi2", sep = "")])
  epiMatrix <- cbind(epiMatrix, epiMatrix)
  resDta <- data.frame(matrix(nrow = nrow(dataset), ncol = 12))
  names(resDta) <- paste(month.abb, month.window, cut.off, sep = ".")
  for(i in 1:12){
    resDta[,i] <-rowSums( epiMatrix[,(1:month.window) + i -1])
  }
  resDta <- resDta >= cut.off
  return(resDta)
}
getBirthMonth <- function(rowOfCountry) {
  resDf <- expand.grid(Month = 1:12,
                       Age = factor(c("<1m", "1-<3m", "<3m","3-<6m", "6-<9m", "9-<12m", "<6m", "<12m")),
                       Outcome = c("RSV-ALRI", "RSV-ALRIcwi", "RSV-ALRIHos"))
  seasonality <- as.numeric(rowOfCountry[paste("p", month.abb, sep = "")])
  resDf$Calendar <- seasonality
  resDf$Country <- rowOfCountry$Country[1]
  resDf$Birth <- NA
  resDf$DE <- rowOfCountry$epi2[1]
  ALRI <- Relative_IR$ALRI[1:12]
  sALRI <- Relative_IR$`ALRI-CWI`[1:12]
  Hos <- Relative_IR$Hos[1:12]
  matrix.ALRI <- matrix(nrow = 12, ncol = 12, data = rep(seasonality, each = 12) * ALRI,
                        byrow = FALSE)
  matrix.sALRI <- matrix(nrow = 12, ncol = 12, data = rep(seasonality, each = 12) * sALRI,
                         byrow = FALSE)
  matrix.Hos <- matrix(nrow = 12, ncol = 12, data = rep(seasonality, each = 12) * Hos,
                       byrow = FALSE)
  matrix.ALRI2 <- matrix.ALRI[,c(1:12,1:12)]
  matrix.sALRI2 <- matrix.sALRI[,c(1:12,1:12)]
  matrix.Hos2 <- matrix.Hos[,c(1:12,1:12)]
  genBirth <- function(input.matrix) {
    resMatrix <- matrix(nrow = 12, ncol = 12, data = 0)
    for(i in 1:12){
      resMatrix[i,] <- input.matrix[i,(0:11)+i]
    }
    resMatrix <- as.data.frame(t(resMatrix))
    resMatrix$`<1m` <- with(resMatrix, V1/sum(V1))
    resMatrix$`1-<3m` <- with(resMatrix, (V2+V3)/sum(V2+V3))
    resMatrix$`<3m` <- with(resMatrix, (V1+V2+V3)/sum(V1+V2+V3))
    resMatrix$`3-<6m` <- with(resMatrix, (V4+V5+V6)/sum(V4+V5+V6))
    resMatrix$`6-<9m` <- with(resMatrix, (V7+V8+V9)/sum(V7+V8+V9))
    resMatrix$`9-<12m` <- with(resMatrix, (V10+V11+V12)/sum(V10+V11+V12))
    resMatrix$`<6m` <- with(resMatrix, (V1+V2+V3+V4+V5+V6)/sum(V1+V2+V3+V4+V5+V6))
    resMatrix$`<12m` <- with(resMatrix, (V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12)/sum(V1+V2+V3+V4+V5+V6+V7+V8+V9+V10+V11+V12))
    resMatrix <- resMatrix[-(1:12)]
    resMatrix <- gather(resMatrix, `<1m`, `1-<3m`,`<3m`,`3-<6m`,`6-<9m`,`9-<12m`,`<6m`,`<12m`, key = "Age", value = "IR")
    return(resMatrix)
  }
  B.ALRI <- genBirth(matrix.ALRI2)
  B.sALRI <- genBirth(matrix.sALRI2)
  B.Hos <- genBirth(matrix.Hos2)
  resDf$Birth <- rbind(rbind(B.ALRI, B.sALRI), B.Hos)$IR
  return(resDf)
}
getECON <- function(rowOfCountry,
                    mAb = efficacyGenerater(0.701, 1, 5),
                    mAbs = efficacyGenerater(0.701, 1, 5),
                    mAbh = efficacyGenerater(0.784, 1, 5),
                    mV = efficacyGenerater(0.394, 1, 3),
                    mVs = efficacyGenerater(0.394, 1, 3),
                    mVh = efficacyGenerater(0.444, 1, 3),
                    deno = 6, mV.coverage = "ANC4"
                    ) {
  resDf <-expand.grid(Approach = c("A", "B", "C", "D","YR"),
                      Immunisation = c("mAb", "mV"),
                      Outcome =c("RSV-ALRI", "RSV-ALRIcwi", "RSV-ALRIHos"))
  resDf <- resDf[!(resDf$Immunisation=="mV" & (resDf$Approach=="C"|resDf$Approach=="D")),]
  resDf$dose <- NA
  resDf$coverage <- NA
  resDf$ideal.prop <- NA
  resDf$real.prop <- NA
  resDf$prop.prop <- NA
  resDf$efficiency <- NA
  resDf$DE <- sum(as.logical(rowOfCountry[paste(month.abb, "epi2", sep = "")]))
  resDf$Approach2 <- rep(c("S1", "S2.mAb", "S3.mAb","S4.mAb","YR", "S1", "S2.mV", "YR"), 3)
  resDf$Country <- rowOfCountry$Country
  seasonality <- as.numeric(rowOfCountry[paste("p", month.abb, sep = "")])
  ALRI <- Relative_IR$ALRI[1:12]
  sALRI <- Relative_IR$`ALRI-CWI`[1:12]
  Hos <- Relative_IR$Hos[1:12]
  matrix.ALRI <- matrix(nrow = 12, ncol = 12, data = rep(seasonality, each = 12) * ALRI,
                        byrow = FALSE)
  matrix.sALRI <- matrix(nrow = 12, ncol = 12, data = rep(seasonality, each = 12) * sALRI,
                         byrow = FALSE)
  matrix.Hos <- matrix(nrow = 12, ncol = 12, data = rep(seasonality, each = 12) * Hos,
                       byrow = FALSE)
  S1 <- rep(as.logical(rowOfCountry[paste(month.abb, "epi2", sep = "")]),2)
  S2.mAb <- rep(as.logical(rowOfCountry[paste(month.abb, ".3.2", sep = "")]),2)
  S3.mAb <- rep(as.logical(rowOfCountry[paste(month.abb, ".4.2", sep = "")]),2)
  S4.mAb <- rep(as.logical(rowOfCountry[paste(month.abb, ".5.2", sep = "")]),2)
  S2.mV <- rep(as.logical(rowOfCountry[paste(month.abb, ".3.2", sep = "")]),2)
  YR <- rep(TRUE, 24)
  resDf$dose <- rep(
    c(sum(S1)/2, sum(S2.mAb)/2, sum(S3.mAb)/2, sum(S4.mAb)/2,12,
      sum(S1)/2, sum(S2.mV)/2, 12), 3
  )
  resDf$coverage <- c(rep(rowOfCountry$mAb, 5),
                      rep(rowOfCountry[[mV.coverage]], 3),
                      rep(rowOfCountry$mAb, 5),
                      rep(rowOfCountry[[mV.coverage]], 3),
                      rep(rowOfCountry$mAb, 5),
                      rep(rowOfCountry[[mV.coverage]], 3))
  getMatrixMask <- function(n.months = 5, logic24, efficacy){
    resMatrix <- matrix(nrow = 12, ncol = 12, data = 0)
    resMatrix[1,] <- logic24[14:25 - 1]
    resMatrix[2,] <- logic24[14:25 - 2]
    resMatrix[3,] <- logic24[14:25 - 3]
    if(n.months ==5) {
      resMatrix[4,] <- logic24[14:25 - 4]
      resMatrix[5,] <- logic24[14:25 - 5]
    }
    resMatrix <- resMatrix * efficacy
    return(resMatrix[1:deno,])
  }
  #ALRI
  resDf$ideal.prop[1] <-
    sum(matrix.ALRI[1:deno,] * getMatrixMask(n.months = 5, S1, mAb))/sum(matrix.ALRI[1:deno,]) * 100
  resDf$ideal.prop[2] <- 
    sum(matrix.ALRI[1:deno,] * getMatrixMask(n.months = 5, S2.mAb, mAb))/sum(matrix.ALRI[1:deno,]) *100
  resDf$ideal.prop[3] <- 
    sum(matrix.ALRI[1:deno,] * getMatrixMask(n.months = 5, S3.mAb, mAb))/sum(matrix.ALRI[1:deno,]) *100
  resDf$ideal.prop[4] <- 
    sum(matrix.ALRI[1:deno,] * getMatrixMask(n.months = 5, S4.mAb, mAb))/sum(matrix.ALRI[1:deno,]) *100
  resDf$ideal.prop[5] <- 
    sum(matrix.ALRI[1:deno,] * getMatrixMask(n.months = 5, YR, mAb))/sum(matrix.ALRI[1:deno,]) *100
  resDf$ideal.prop[6] <- 
    sum(matrix.ALRI[1:deno,] * getMatrixMask(n.months = 3, S1, mV))/sum(matrix.ALRI[1:deno,]) *100
  resDf$ideal.prop[7] <- 
    sum(matrix.ALRI[1:deno,] * getMatrixMask(n.months = 5, S2.mV, mV))/sum(matrix.ALRI[1:deno,]) *100
  resDf$ideal.prop[8] <- 
    sum(matrix.ALRI[1:deno,] * getMatrixMask(n.months = 5, YR, mV))/sum(matrix.ALRI[1:deno,]) *100
  resDf$efficiency[1] <- resDf$ideal.prop[1]/resDf$dose[1]/(resDf$ideal.prop[5]/12)
  resDf$efficiency[2] <- resDf$ideal.prop[2]/resDf$dose[2]/(resDf$ideal.prop[5]/12)
  resDf$efficiency[3] <- resDf$ideal.prop[3]/resDf$dose[3]/(resDf$ideal.prop[5]/12)
  resDf$efficiency[4] <- resDf$ideal.prop[4]/resDf$dose[4]/(resDf$ideal.prop[5]/12)
  resDf$efficiency[5] <- resDf$ideal.prop[5]/resDf$dose[5]/(resDf$ideal.prop[5]/12)
  resDf$efficiency[6] <- resDf$ideal.prop[6]/resDf$dose[6]/(resDf$ideal.prop[8]/12)
  resDf$efficiency[7] <- resDf$ideal.prop[7]/resDf$dose[7]/(resDf$ideal.prop[8]/12)
  resDf$efficiency[8] <- resDf$ideal.prop[8]/resDf$dose[8]/(resDf$ideal.prop[8]/12)
  resDf$prop.prop[1] <- resDf$ideal.prop[1]/resDf$ideal.prop[5]
  resDf$prop.prop[2] <- resDf$ideal.prop[2]/resDf$ideal.prop[5]
  resDf$prop.prop[3] <- resDf$ideal.prop[3]/resDf$ideal.prop[5]
  resDf$prop.prop[4] <- resDf$ideal.prop[4]/resDf$ideal.prop[5]
  resDf$prop.prop[5] <- resDf$ideal.prop[5]/resDf$ideal.prop[5]
  resDf$prop.prop[6] <- resDf$ideal.prop[6]/resDf$ideal.prop[8]
  resDf$prop.prop[7] <- resDf$ideal.prop[7]/resDf$ideal.prop[8]
  resDf$prop.prop[8] <- resDf$ideal.prop[8]/resDf$ideal.prop[8]
  #sALRI
  resDf$ideal.prop[9] <-
    sum(matrix.sALRI[1:deno,] * getMatrixMask(n.months = 5, S1, mAbs))/sum(matrix.sALRI[1:deno,]) * 100
  resDf$ideal.prop[10] <- 
    sum(matrix.sALRI[1:deno,] * getMatrixMask(n.months = 5, S2.mAb, mAbs))/sum(matrix.sALRI[1:deno,]) *100
  resDf$ideal.prop[11] <- 
    sum(matrix.sALRI[1:deno,] * getMatrixMask(n.months = 5, S3.mAb, mAbs))/sum(matrix.sALRI[1:deno,]) *100
  resDf$ideal.prop[12] <- 
    sum(matrix.sALRI[1:deno,] * getMatrixMask(n.months = 5, S4.mAb, mAbs))/sum(matrix.sALRI[1:deno,]) *100
  resDf$ideal.prop[13] <- 
    sum(matrix.sALRI[1:deno,] * getMatrixMask(n.months = 5, YR, mAbs))/sum(matrix.sALRI[1:deno,]) *100
  resDf$ideal.prop[14] <- 
    sum(matrix.sALRI[1:deno,] * getMatrixMask(n.months = 3, S1, mVs))/sum(matrix.sALRI[1:deno,]) *100
  resDf$ideal.prop[15] <- 
    sum(matrix.sALRI[1:deno,] * getMatrixMask(n.months = 5, S2.mV, mVs))/sum(matrix.sALRI[1:deno,]) *100
  resDf$ideal.prop[16] <- 
    sum(matrix.sALRI[1:deno,] * getMatrixMask(n.months = 5, YR, mVs))/sum(matrix.sALRI[1:deno,]) *100
  resDf$efficiency[9] <- resDf$ideal.prop[9]/resDf$dose[9]/(resDf$ideal.prop[13]/12)
  resDf$efficiency[10] <- resDf$ideal.prop[10]/resDf$dose[10]/(resDf$ideal.prop[13]/12)
  resDf$efficiency[11] <- resDf$ideal.prop[11]/resDf$dose[11]/(resDf$ideal.prop[13]/12)
  resDf$efficiency[12] <- resDf$ideal.prop[12]/resDf$dose[12]/(resDf$ideal.prop[13]/12)
  resDf$efficiency[13] <- resDf$ideal.prop[13]/resDf$dose[13]/(resDf$ideal.prop[13]/12)
  resDf$efficiency[14] <- resDf$ideal.prop[14]/resDf$dose[14]/(resDf$ideal.prop[16]/12)
  resDf$efficiency[15] <- resDf$ideal.prop[15]/resDf$dose[15]/(resDf$ideal.prop[16]/12)
  resDf$efficiency[16] <- resDf$ideal.prop[16]/resDf$dose[16]/(resDf$ideal.prop[16]/12)
  resDf$prop.prop[9] <- resDf$ideal.prop[9]/resDf$ideal.prop[13]
  resDf$prop.prop[10] <- resDf$ideal.prop[10]/resDf$ideal.prop[13]
  resDf$prop.prop[11] <- resDf$ideal.prop[11]/resDf$ideal.prop[13]
  resDf$prop.prop[12] <- resDf$ideal.prop[12]/resDf$ideal.prop[13]
  resDf$prop.prop[13] <- resDf$ideal.prop[13]/resDf$ideal.prop[13]
  resDf$prop.prop[14] <- resDf$ideal.prop[14]/resDf$ideal.prop[16]
  resDf$prop.prop[15] <- resDf$ideal.prop[15]/resDf$ideal.prop[16]
  resDf$prop.prop[16] <- resDf$ideal.prop[16]/resDf$ideal.prop[16]
  #Hos
  resDf$ideal.prop[17] <-
    sum(matrix.Hos[1:deno,] * getMatrixMask(n.months = 5, S1, mAbh))/sum(matrix.Hos[1:deno,]) * 100
  resDf$ideal.prop[18] <- 
    sum(matrix.Hos[1:deno,] * getMatrixMask(n.months = 5, S2.mAb, mAbh))/sum(matrix.Hos[1:deno,]) *100
  resDf$ideal.prop[19] <- 
    sum(matrix.Hos[1:deno,] * getMatrixMask(n.months = 5, S3.mAb, mAbh))/sum(matrix.Hos[1:deno,]) *100
  resDf$ideal.prop[20] <- 
    sum(matrix.Hos[1:deno,] * getMatrixMask(n.months = 5, S4.mAb, mAbh))/sum(matrix.Hos[1:deno,]) *100
  resDf$ideal.prop[21] <- 
    sum(matrix.Hos[1:deno,] * getMatrixMask(n.months = 5, YR, mAbh))/sum(matrix.Hos[1:deno,]) *100
  resDf$ideal.prop[22] <- 
    sum(matrix.Hos[1:deno,] * getMatrixMask(n.months = 3, S1, mVh))/sum(matrix.Hos[1:deno,]) *100
  resDf$ideal.prop[23] <- 
    sum(matrix.Hos[1:deno,] * getMatrixMask(n.months = 5, S2.mV, mVh))/sum(matrix.Hos[1:deno,]) *100
  resDf$ideal.prop[24] <- 
    sum(matrix.Hos[1:deno,] * getMatrixMask(n.months = 5, YR, mVh))/sum(matrix.Hos[1:deno,]) *100
  resDf$efficiency[17] <- resDf$ideal.prop[17]/resDf$dose[17]/(resDf$ideal.prop[21]/12)
  resDf$efficiency[18] <- resDf$ideal.prop[18]/resDf$dose[18]/(resDf$ideal.prop[21]/12)
  resDf$efficiency[19] <- resDf$ideal.prop[19]/resDf$dose[19]/(resDf$ideal.prop[21]/12)
  resDf$efficiency[20] <- resDf$ideal.prop[20]/resDf$dose[20]/(resDf$ideal.prop[21]/12)
  resDf$efficiency[21] <- resDf$ideal.prop[21]/resDf$dose[21]/(resDf$ideal.prop[21]/12)
  resDf$efficiency[22] <- resDf$ideal.prop[22]/resDf$dose[22]/(resDf$ideal.prop[24]/12)
  resDf$efficiency[23] <- resDf$ideal.prop[23]/resDf$dose[23]/(resDf$ideal.prop[24]/12)
  resDf$efficiency[24] <- resDf$ideal.prop[24]/resDf$dose[24]/(resDf$ideal.prop[24]/12)
  resDf$prop.prop[17] <- resDf$ideal.prop[17]/resDf$ideal.prop[21]
  resDf$prop.prop[18] <- resDf$ideal.prop[18]/resDf$ideal.prop[21]
  resDf$prop.prop[19] <- resDf$ideal.prop[19]/resDf$ideal.prop[21]
  resDf$prop.prop[20] <- resDf$ideal.prop[20]/resDf$ideal.prop[21]
  resDf$prop.prop[21] <- resDf$ideal.prop[21]/resDf$ideal.prop[21]
  resDf$prop.prop[22] <- resDf$ideal.prop[22]/resDf$ideal.prop[24]
  resDf$prop.prop[23] <- resDf$ideal.prop[23]/resDf$ideal.prop[24]
  resDf$prop.prop[24] <- resDf$ideal.prop[24]/resDf$ideal.prop[24]
  resDf$real.prop <- resDf$ideal.prop*resDf$coverage/100
  return(resDf)
}
getECONM <- function(rowsOfCountry,mAb = efficacyGenerater(0.701, 1, 5),
                     mAbs = efficacyGenerater(0.701, 1, 5),
                     mAbh = efficacyGenerater(0.784, 1, 5),
                     mV = efficacyGenerater(0.394, 1, 3),
                     mVs = efficacyGenerater(0.394, 1, 3),
                     mVh = efficacyGenerater(0.444, 1, 3),
                     deno = 6, mV.coverage = "ANC4"
                     ) {
  row.average <- rowsOfCountry[rowsOfCountry$YearIndex=="Average",]
  row.rest <- rowsOfCountry[rowsOfCountry$YearIndex!="Average",]
  n.year <- nrow(row.rest)
  resDf <-   expand.grid(Approach = c("A", "B", "C", "D","YR"),
                         Immunisation = c("mAb", "mV"),
                         Outcome = c("RSV-ALRI", "RSV-ALRIcwi", "RSV-ALRIHos"),
                         YearIndex = 1:n.year)
  resDf <- resDf[!(resDf$Immunisation=="mV" & (resDf$Approach=="C"|resDf$Approach=="D")),]
  resDf$dose <- NA
  resDf$coverage <- NA
  resDf$ideal.prop <- NA
  resDf$real.prop <- NA
  resDf$prop.prop <- NA
  resDf$efficiency <- NA
  resDf$DE <- sum(as.logical(row.average[paste(month.abb, "epi2", sep = "")]))
  resDf$Approach2 <- rep(c("S1", "S2.mAb", "S3.mAb", "S4.mAb","YR", "S1", "S2.mV", "YR"), 3)
  resDf$Country <- row.average$Country
  seasonality <- c(t(row.rest[paste("p", month.abb, sep = "")]))
  ALRI <- Relative_IR$ALRI[1:12]
  sALRI <- Relative_IR$`ALRI-CWI`[1:12]
  Hos <- Relative_IR$Hos[1:12]
  matrix.ALRI <- matrix(nrow = 12 , ncol = 12*n.year, 
                        data = rep(seasonality, each = 12) * ALRI,
                        byrow = FALSE)
  matrix.sALRI <- matrix(nrow = 12, ncol = 12*n.year,
                         data = rep(seasonality, each = 12) * sALRI,
                         byrow = FALSE)
  matrix.Hos <- matrix(nrow = 12, ncol = 12*n.year,
                       data = rep(seasonality, each = 12) * Hos,
                       byrow = FALSE)
  S1 <- rep(as.logical(row.average[paste(month.abb, "epi2", sep = "")]),2)
  S2.mAb <- rep(as.logical(row.average[paste(month.abb, ".3.2", sep = "")]),2)
  S3.mAb <- rep(as.logical(row.average[paste(month.abb, ".4.2", sep = "")]),2)
  S4.mAb <- rep(as.logical(row.average[paste(month.abb, ".5.2", sep = "")]),2)
  S2.mV <- rep(as.logical(row.average[paste(month.abb, ".3.2", sep = "")]),2)
  YR <- rep(TRUE, 24)
  resDf$dose <- rep(
    c(sum(S1)/2, sum(S2.mAb)/2,sum(S3.mAb)/2, sum(S4.mAb)/2, 12,
      sum(S1)/2, sum(S2.mV)/2, 12), 3
  )
  resDf$coverage <- c(rep(row.average$mAb, 5),
                      rep(row.average[[mV.coverage]], 3),
                      rep(row.average$mAb, 5),
                      rep(row.average[[mV.coverage]], 3),
                      rep(row.average$mAb, 5),
                      rep(row.average[[mV.coverage]], 3))
  getMatrixMask <- function(n.months = 5, logic24, efficacy){
    resMatrix <- matrix(nrow = 12, ncol = 12, data = 0)
    resMatrix[1,] <- logic24[14:25 - 1]
    resMatrix[2,] <- logic24[14:25 - 2]
    resMatrix[3,] <- logic24[14:25 - 3]
    if(n.months ==5) {
      resMatrix[4,] <- logic24[14:25 - 4]
      resMatrix[5,] <- logic24[14:25 - 5]
    }
    resMatrix <- resMatrix * efficacy
    return(resMatrix[1:deno,])
  }
  for(i in 1:n.year) {
    #ALRI
    resDf$ideal.prop[1+ (i-1)*24] <-
      sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S1, mAb))/sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)]) * 100
    resDf$ideal.prop[2+ (i-1)*24] <- 
      sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S2.mAb, mAb))/sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[3+ (i-1)*24] <- 
      sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S3.mAb, mAb))/sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[4+ (i-1)*24] <- 
      sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S4.mAb, mAb))/sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[5+ (i-1)*24] <- 
      sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, YR, mAb))/sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[6+ (i-1)*24] <- 
      sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 3, S1, mV))/sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[7+ (i-1)*24] <- 
      sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S2.mV, mV))/sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[8+ (i-1)*24] <- 
      sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, YR, mV))/sum(matrix.ALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$efficiency[1+ (i-1)*24] <- resDf$ideal.prop[1+ (i-1)*24]/resDf$dose[1+ (i-1)*24]/(resDf$ideal.prop[5+ (i-1)*24]/12)
    resDf$efficiency[2+ (i-1)*24] <- resDf$ideal.prop[2+ (i-1)*24]/resDf$dose[2+ (i-1)*24]/(resDf$ideal.prop[5+ (i-1)*24]/12)
    resDf$efficiency[3+ (i-1)*24] <- resDf$ideal.prop[3+ (i-1)*24]/resDf$dose[3+ (i-1)*24]/(resDf$ideal.prop[5+ (i-1)*24]/12)
    resDf$efficiency[4+ (i-1)*24] <- resDf$ideal.prop[4+ (i-1)*24]/resDf$dose[4+ (i-1)*24]/(resDf$ideal.prop[5+ (i-1)*24]/12)
    resDf$efficiency[5+ (i-1)*24] <- resDf$ideal.prop[5+ (i-1)*24]/resDf$dose[5+ (i-1)*24]/(resDf$ideal.prop[5+ (i-1)*24]/12)
    resDf$efficiency[6+ (i-1)*24] <- resDf$ideal.prop[6+ (i-1)*24]/resDf$dose[6+ (i-1)*24]/(resDf$ideal.prop[8+ (i-1)*24]/12)
    resDf$efficiency[7+ (i-1)*24] <- resDf$ideal.prop[7+ (i-1)*24]/resDf$dose[7+ (i-1)*24]/(resDf$ideal.prop[8+ (i-1)*24]/12)
    resDf$efficiency[8+ (i-1)*24] <- resDf$ideal.prop[8+ (i-1)*24]/resDf$dose[8+ (i-1)*24]/(resDf$ideal.prop[8+ (i-1)*24]/12)
    resDf$prop.prop[1+ (i-1)*24] <- resDf$ideal.prop[1+ (i-1)*24]/resDf$ideal.prop[5+ (i-1)*24]
    resDf$prop.prop[2+ (i-1)*24] <- resDf$ideal.prop[2+ (i-1)*24]/resDf$ideal.prop[5+ (i-1)*24]
    resDf$prop.prop[3+ (i-1)*24] <- resDf$ideal.prop[3+ (i-1)*24]/resDf$ideal.prop[5+ (i-1)*24]
    resDf$prop.prop[4+ (i-1)*24] <- resDf$ideal.prop[4+ (i-1)*24]/resDf$ideal.prop[5+ (i-1)*24]
    resDf$prop.prop[5+ (i-1)*24] <- resDf$ideal.prop[5+ (i-1)*24]/resDf$ideal.prop[5+ (i-1)*24]
    resDf$prop.prop[6+ (i-1)*24] <- resDf$ideal.prop[6+ (i-1)*24]/resDf$ideal.prop[8+ (i-1)*24]
    resDf$prop.prop[7+ (i-1)*24] <- resDf$ideal.prop[7+ (i-1)*24]/resDf$ideal.prop[8+ (i-1)*24]
    resDf$prop.prop[8+ (i-1)*24] <- resDf$ideal.prop[8+ (i-1)*24]/resDf$ideal.prop[8+ (i-1)*24]
    #sALRI
    resDf$ideal.prop[9+ (i-1)*24] <-
      sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S1, mAbs))/sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)]) * 100
    resDf$ideal.prop[10+ (i-1)*24] <- 
      sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S2.mAb, mAbs))/sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[11+ (i-1)*24] <- 
      sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S3.mAb, mAbs))/sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[12+ (i-1)*24] <- 
      sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S4.mAb, mAbs))/sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[13+ (i-1)*24] <- 
      sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, YR, mAbs))/sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[14+ (i-1)*24] <- 
      sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 3, S1, mVs))/sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[15+ (i-1)*24] <- 
      sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S2.mV, mVs))/sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[16+ (i-1)*24] <- 
      sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, YR, mVs))/sum(matrix.sALRI[1:deno,(1:12)+12*(i-1)]) *100
    resDf$efficiency[9+ (i-1)*24] <- resDf$ideal.prop[9+ (i-1)*24]/resDf$dose[9+ (i-1)*24]/(resDf$ideal.prop[13+ (i-1)*24]/12)
    resDf$efficiency[10+ (i-1)*24] <- resDf$ideal.prop[10+ (i-1)*24]/resDf$dose[10+ (i-1)*24]/(resDf$ideal.prop[13+ (i-1)*24]/12)
    resDf$efficiency[11+ (i-1)*24] <- resDf$ideal.prop[11+ (i-1)*24]/resDf$dose[11+ (i-1)*24]/(resDf$ideal.prop[13+ (i-1)*24]/12)
    resDf$efficiency[12+ (i-1)*24] <- resDf$ideal.prop[12+ (i-1)*24]/resDf$dose[12+ (i-1)*24]/(resDf$ideal.prop[13+ (i-1)*24]/12)
    resDf$efficiency[13+ (i-1)*24] <- resDf$ideal.prop[13+ (i-1)*24]/resDf$dose[13+ (i-1)*24]/(resDf$ideal.prop[13+ (i-1)*24]/12)
    resDf$efficiency[14+ (i-1)*24] <- resDf$ideal.prop[14+ (i-1)*24]/resDf$dose[14+ (i-1)*24]/(resDf$ideal.prop[16+ (i-1)*24]/12)
    resDf$efficiency[15+ (i-1)*24] <- resDf$ideal.prop[15+ (i-1)*24]/resDf$dose[15+ (i-1)*24]/(resDf$ideal.prop[16+ (i-1)*24]/12)
    resDf$efficiency[16+ (i-1)*24] <- resDf$ideal.prop[16+ (i-1)*24]/resDf$dose[16+ (i-1)*24]/(resDf$ideal.prop[16+ (i-1)*24]/12)
    resDf$prop.prop[9+ (i-1)*24] <- resDf$ideal.prop[9+ (i-1)*24]/resDf$ideal.prop[13+ (i-1)*24]
    resDf$prop.prop[10+ (i-1)*24] <- resDf$ideal.prop[10+ (i-1)*24]/resDf$ideal.prop[13+ (i-1)*24]
    resDf$prop.prop[11+ (i-1)*24] <- resDf$ideal.prop[11+ (i-1)*24]/resDf$ideal.prop[13+ (i-1)*24]
    resDf$prop.prop[12+ (i-1)*24] <- resDf$ideal.prop[12+ (i-1)*24]/resDf$ideal.prop[13+ (i-1)*24]
    resDf$prop.prop[13+ (i-1)*24] <- resDf$ideal.prop[13+ (i-1)*24]/resDf$ideal.prop[13+ (i-1)*24]
    resDf$prop.prop[14+ (i-1)*24] <- resDf$ideal.prop[14+ (i-1)*24]/resDf$ideal.prop[16+ (i-1)*24]
    resDf$prop.prop[15+ (i-1)*24] <- resDf$ideal.prop[15+ (i-1)*24]/resDf$ideal.prop[16+ (i-1)*24]
    resDf$prop.prop[16+ (i-1)*24] <- resDf$ideal.prop[16+ (i-1)*24]/resDf$ideal.prop[16+ (i-1)*24]
    #Hos
    resDf$ideal.prop[17+ (i-1)*24] <-
      sum(matrix.Hos[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S1, mAbh))/sum(matrix.Hos[1:deno,(1:12)+12*(i-1)]) * 100
    resDf$ideal.prop[18+ (i-1)*24] <- 
      sum(matrix.Hos[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S2.mAb, mAbh))/sum(matrix.Hos[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[19+ (i-1)*24] <- 
      sum(matrix.Hos[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S3.mAb, mAbh))/sum(matrix.Hos[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[20+ (i-1)*24] <- 
      sum(matrix.Hos[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S4.mAb, mAbh))/sum(matrix.Hos[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[21+ (i-1)*24] <- 
      sum(matrix.Hos[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, YR, mAbh))/sum(matrix.Hos[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[22+ (i-1)*24] <- 
      sum(matrix.Hos[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 3, S1, mVh))/sum(matrix.Hos[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[23+ (i-1)*24] <- 
      sum(matrix.Hos[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, S2.mV, mVh))/sum(matrix.Hos[1:deno,(1:12)+12*(i-1)]) *100
    resDf$ideal.prop[24+ (i-1)*24] <- 
      sum(matrix.Hos[1:deno,(1:12)+12*(i-1)] * getMatrixMask(n.months = 5, YR, mVh))/sum(matrix.Hos[1:deno,(1:12)+12*(i-1)]) *100
    resDf$efficiency[17+ (i-1)*24] <- resDf$ideal.prop[17+ (i-1)*24]/resDf$dose[17+ (i-1)*24]/(resDf$ideal.prop[21+ (i-1)*24]/12)
    resDf$efficiency[18+ (i-1)*24] <- resDf$ideal.prop[18+ (i-1)*24]/resDf$dose[18+ (i-1)*24]/(resDf$ideal.prop[21+ (i-1)*24]/12)
    resDf$efficiency[19+ (i-1)*24] <- resDf$ideal.prop[19+ (i-1)*24]/resDf$dose[19+ (i-1)*24]/(resDf$ideal.prop[21+ (i-1)*24]/12)
    resDf$efficiency[20+ (i-1)*24] <- resDf$ideal.prop[20+ (i-1)*24]/resDf$dose[20+ (i-1)*24]/(resDf$ideal.prop[21+ (i-1)*24]/12)
    resDf$efficiency[21+ (i-1)*24] <- resDf$ideal.prop[21+ (i-1)*24]/resDf$dose[21+ (i-1)*24]/(resDf$ideal.prop[21+ (i-1)*24]/12)
    resDf$efficiency[22+ (i-1)*24] <- resDf$ideal.prop[22+ (i-1)*24]/resDf$dose[22+ (i-1)*24]/(resDf$ideal.prop[24+ (i-1)*24]/12)
    resDf$efficiency[23+ (i-1)*24] <- resDf$ideal.prop[23+ (i-1)*24]/resDf$dose[23+ (i-1)*24]/(resDf$ideal.prop[24+ (i-1)*24]/12)
    resDf$efficiency[24+ (i-1)*24] <- resDf$ideal.prop[24+ (i-1)*24]/resDf$dose[24+ (i-1)*24]/(resDf$ideal.prop[24+ (i-1)*24]/12)
    resDf$prop.prop[17+ (i-1)*24] <- resDf$ideal.prop[17+ (i-1)*24]/resDf$ideal.prop[21+ (i-1)*24]
    resDf$prop.prop[18+ (i-1)*24] <- resDf$ideal.prop[18+ (i-1)*24]/resDf$ideal.prop[21+ (i-1)*24]
    resDf$prop.prop[19+ (i-1)*24] <- resDf$ideal.prop[19+ (i-1)*24]/resDf$ideal.prop[21+ (i-1)*24]
    resDf$prop.prop[20+ (i-1)*24] <- resDf$ideal.prop[20+ (i-1)*24]/resDf$ideal.prop[21+ (i-1)*24]
    resDf$prop.prop[21+ (i-1)*24] <- resDf$ideal.prop[21+ (i-1)*24]/resDf$ideal.prop[21+ (i-1)*24]
    resDf$prop.prop[22+ (i-1)*24] <- resDf$ideal.prop[22+ (i-1)*24]/resDf$ideal.prop[24+ (i-1)*24]
    resDf$prop.prop[23+ (i-1)*24] <- resDf$ideal.prop[23+ (i-1)*24]/resDf$ideal.prop[24+ (i-1)*24]
    resDf$prop.prop[24+ (i-1)*24] <- resDf$ideal.prop[24+ (i-1)*24]/resDf$ideal.prop[24+ (i-1)*24]
  }
  resDf$real.prop <- resDf$ideal.prop*resDf$coverage/100
  return(resDf)
}
getEachCountry <- function(rowOfCountry) {
  mAb = 0.701
  mAbs = 0.701
  mAbh = 0.784
  mV = 0.394
  mVs = 0.394
  mVh = 0.444
  resDf <-expand.grid(Approach = c("A", "B", "C", "D","YR"),
                      Immunisation = c("mAb", "mV"),
                      Outcome =c("RSV-ALRI", "RSV-ALRIcwi", "RSV-ALRIHos"))
  resDf <- resDf[!(resDf$Immunisation=="mV" & (resDf$Approach=="C"|resDf$Approach=="D")),]
  resDf$prop <- NA
  resDf$dose.prop <- NA
  resDf$dose <- NA
  resDf$coverage <- NA
  resDf$targeted.prop <- NA
  resDf$dose.targeted.prop <- NA
  resDf$prop.prop <- NA
  resDf$dose.prop.ratio <- NA
  resDf$DE <- sum(as.logical(rowOfCountry[paste(month.abb, "epi2", sep = "")]))
  resDf$Approach2 <- rep(c("S1", "S2.mAb", "S3.mAb","S4.mAb","YR", "S1", "S2.mV", "YR"), 3)
  resDf$Country <- rowOfCountry$Country
  # POP = 1000
  seasonality <- as.numeric(rowOfCountry[paste("p", month.abb, sep = "")])
  ALRI <- c(rep(rowOfCountry$est11, 3),
            rep(rowOfCountry$est12, 3),
            rep(rowOfCountry$est13, 3),
            rep(rowOfCountry$est14, 3))
  sALRI <- c(rep(rowOfCountry$est21, 3),
             rep(rowOfCountry$est22, 3),
             rep(rowOfCountry$est23, 3),
             rep(rowOfCountry$est24, 3))
  Hos <- c(rep(0.361, 3),
           rep(0.283, 3),
           rep(0.215, 3),
           rep(0.142, 3))
  matrix.ALRI <- matrix(nrow = 12, ncol = 12, data = rep(seasonality, each = 12) * ALRI,
                        byrow = FALSE)
  matrix.sALRI <- matrix(nrow = 12, ncol = 12, data = rep(seasonality, each = 12) * sALRI,
                         byrow = FALSE)
  matrix.Hos <- matrix(nrow = 12, ncol = 12, data = rep(seasonality, each = 12) * Hos,
                       byrow = FALSE)
  S1 <- rep(as.logical(rowOfCountry[paste(month.abb, "epi2", sep = "")]),2)
  S2.mAb <- rep(as.logical(rowOfCountry[paste(month.abb, ".3.2", sep = "")]),2)
  S3.mAb <- rep(as.logical(rowOfCountry[paste(month.abb, ".4.2", sep = "")]),2)
  S4.mAb <- rep(as.logical(rowOfCountry[paste(month.abb, ".5.2", sep = "")]),2)
  S2.mV <- rep(as.logical(rowOfCountry[paste(month.abb, ".3.2", sep = "")]),2)
  YR <- rep(TRUE, 24)
  resDf$dose <- rep(
    c(sum(S1)/2, sum(S2.mAb)/2, sum(S3.mAb)/2, sum(S4.mAb)/2,12,
      sum(S1)/2, sum(S2.mV)/2, 12), 3
  )
  resDf$coverage <- c(rep(rowOfCountry$mAb, 5),
                      rep(rowOfCountry$ANC4, 3),
                      rep(rowOfCountry$mAb, 5),
                      rep(rowOfCountry$ANC4, 3),
                      rep(rowOfCountry$mAb, 5),
                      rep(rowOfCountry$ANC4, 3))
  getMatrixMask <- function(n.months = 5, logic24){
    resMatrix <- matrix(nrow = 12, ncol = 12, data = 0)
    resMatrix[1,] <- logic24[14:25 - 1]
    resMatrix[2,] <- logic24[14:25 - 2]
    resMatrix[3,] <- logic24[14:25 - 3]
    if(n.months ==5) {
      resMatrix[4,] <- logic24[14:25 - 4]
      resMatrix[5,] <- logic24[14:25 - 5]
    }
    return(resMatrix)
  }
  resDf$targeted.prop[1] <- 
    sum(matrix.ALRI[getMatrixMask(n.months = 5, S1)])/sum(matrix.ALRI) *100
  resDf$targeted.prop[2] <- 
    sum(matrix.ALRI[getMatrixMask(n.months = 5, S2.mAb)])/sum(matrix.ALRI) *100
  resDf$targeted.prop[3] <- 
    sum(matrix.ALRI[getMatrixMask(n.months = 5, S3.mAb)])/sum(matrix.ALRI) *100
  resDf$targeted.prop[4] <- 
    sum(matrix.ALRI[getMatrixMask(n.months = 5, S4.mAb)])/sum(matrix.ALRI) *100
  resDf$targeted.prop[5] <- 
    sum(matrix.ALRI[getMatrixMask(n.months = 5, YR)])/sum(matrix.ALRI) *100
  resDf$targeted.prop[6] <- 
    sum(matrix.ALRI[getMatrixMask(n.months = 3, S1)])/sum(matrix.ALRI) *100
  resDf$targeted.prop[7] <- 
    sum(matrix.ALRI[getMatrixMask(n.months = 3, S2.mV)])/sum(matrix.ALRI) *100
  resDf$targeted.prop[8] <- 
    sum(matrix.ALRI[getMatrixMask(n.months = 3, YR)])/sum(matrix.ALRI) *100
  resDf$prop[1] <- resDf$targeted.prop[1] * rowOfCountry$mAb/100 * mAb
  resDf$prop[2] <- resDf$targeted.prop[2] * rowOfCountry$mAb/100 * mAb
  resDf$prop[3] <- resDf$targeted.prop[3] * rowOfCountry$mAb/100 * mAb
  resDf$prop[4] <- resDf$targeted.prop[4] * rowOfCountry$mAb/100 * mAb
  resDf$prop[5] <- resDf$targeted.prop[5] * rowOfCountry$mAb/100 * mAb
  resDf$prop[6] <- resDf$targeted.prop[6] * rowOfCountry$ANC4/100 * mV
  resDf$prop[7] <- resDf$targeted.prop[7] * rowOfCountry$ANC4/100 * mV
  resDf$prop[8] <- resDf$targeted.prop[8] * rowOfCountry$ANC4/100 * mV
  
  resDf$targeted.prop[9] <- 
    sum(matrix.sALRI[getMatrixMask(n.months = 5, S1)])/sum(matrix.sALRI) *100
  resDf$targeted.prop[10] <- 
    sum(matrix.sALRI[getMatrixMask(n.months = 5, S2.mAb)])/sum(matrix.sALRI) *100
  resDf$targeted.prop[11] <- 
    sum(matrix.sALRI[getMatrixMask(n.months = 5, S3.mAb)])/sum(matrix.sALRI) *100
  resDf$targeted.prop[12] <- 
    sum(matrix.sALRI[getMatrixMask(n.months = 5, S4.mAb)])/sum(matrix.sALRI) *100
  resDf$targeted.prop[13] <- 
    sum(matrix.sALRI[getMatrixMask(n.months = 5, YR)])/sum(matrix.sALRI) *100
  resDf$targeted.prop[14] <- 
    sum(matrix.sALRI[getMatrixMask(n.months = 3, S1)])/sum(matrix.sALRI) *100
  resDf$targeted.prop[15] <- 
    sum(matrix.sALRI[getMatrixMask(n.months = 3, S2.mV)])/sum(matrix.sALRI) *100
  resDf$targeted.prop[16] <- 
    sum(matrix.sALRI[getMatrixMask(n.months = 3, YR)])/sum(matrix.sALRI) *100
  resDf$prop[9] <- resDf$targeted.prop[9] * rowOfCountry$mAb/100 * mAbs
  resDf$prop[10] <- resDf$targeted.prop[10] * rowOfCountry$mAb/100 * mAbs
  resDf$prop[11] <- resDf$targeted.prop[11] * rowOfCountry$mAb/100 * mAbs
  resDf$prop[12] <- resDf$targeted.prop[12] * rowOfCountry$mAb/100 * mAbs
  resDf$prop[13] <- resDf$targeted.prop[13] * rowOfCountry$mAb/100 * mAbs
  resDf$prop[14] <- resDf$targeted.prop[14] * rowOfCountry$ANC4/100 * mVs
  resDf$prop[15] <- resDf$targeted.prop[15] * rowOfCountry$ANC4/100 * mVs
  resDf$prop[16] <- resDf$targeted.prop[16] * rowOfCountry$ANC4/100 * mVs
  
  
  resDf$targeted.prop[17] <- 
    sum(matrix.Hos[getMatrixMask(n.months = 5, S1)])/sum(matrix.Hos) *100
  resDf$targeted.prop[18] <- 
    sum(matrix.Hos[getMatrixMask(n.months = 5, S2.mAb)])/sum(matrix.Hos) *100
  resDf$targeted.prop[19] <- 
    sum(matrix.Hos[getMatrixMask(n.months = 5, S3.mAb)])/sum(matrix.Hos) *100
  resDf$targeted.prop[20] <- 
    sum(matrix.Hos[getMatrixMask(n.months = 5, S4.mAb)])/sum(matrix.Hos) *100
  resDf$targeted.prop[21] <- 
    sum(matrix.Hos[getMatrixMask(n.months = 5, YR)])/sum(matrix.Hos) *100
  resDf$targeted.prop[22] <- 
    sum(matrix.Hos[getMatrixMask(n.months = 3, S1)])/sum(matrix.Hos) *100
  resDf$targeted.prop[23] <- 
    sum(matrix.Hos[getMatrixMask(n.months = 3, S2.mV)])/sum(matrix.Hos) *100
  resDf$targeted.prop[24] <- 
    sum(matrix.Hos[getMatrixMask(n.months = 3, YR)])/sum(matrix.Hos) *100
  resDf$prop[17] <- resDf$targeted.prop[17] * rowOfCountry$mAb/100 * mAbh
  resDf$prop[18] <- resDf$targeted.prop[18] * rowOfCountry$mAb/100 * mAbh
  resDf$prop[19] <- resDf$targeted.prop[19] * rowOfCountry$mAb/100 * mAbh
  resDf$prop[20] <- resDf$targeted.prop[20] * rowOfCountry$mAb/100 * mAbh
  resDf$prop[21] <- resDf$targeted.prop[21] * rowOfCountry$mAb/100 * mAbh
  resDf$prop[22] <- resDf$targeted.prop[22] * rowOfCountry$ANC4/100 * mVh
  resDf$prop[23] <- resDf$targeted.prop[23] * rowOfCountry$ANC4/100 * mVh
  resDf$prop[24] <- resDf$targeted.prop[24] * rowOfCountry$ANC4/100 * mVh
  
  resDf$dose.prop <- resDf$prop/resDf$dose
  ## Add prop and ratio
  for(i in 0:2){
    resDf$prop.prop[1 + i*8] <- resDf$prop[1+ i*8]/resDf$prop[5+ i*8] *100
    resDf$prop.prop[2+ i*8] <- resDf$prop[2+ i*8]/resDf$prop[5+ i*8] *100
    resDf$prop.prop[3+ i*8] <- resDf$prop[3+ i*8]/resDf$prop[5+ i*8] *100
    resDf$prop.prop[4+ i*8] <- resDf$prop[4+ i*8]/resDf$prop[5+ i*8] *100
    resDf$prop.prop[5+ i*8] <- resDf$prop[5+ i*8]/resDf$prop[5+ i*8] *100
    resDf$prop.prop[6+ i*8] <- resDf$prop[6+ i*8]/resDf$prop[8+ i*8] *100
    resDf$prop.prop[7+ i*8] <- resDf$prop[7+ i*8]/resDf$prop[8+ i*8] *100
    resDf$prop.prop[8+ i*8] <- resDf$prop[8+ i*8]/resDf$prop[8+ i*8] *100
    resDf$dose.prop.ratio[1+ i*8] <- resDf$dose.prop[1+ i*8]/resDf$dose.prop[5+ i*8]
    resDf$dose.prop.ratio[2+ i*8] <- resDf$dose.prop[2+ i*8]/resDf$dose.prop[5+ i*8]
    resDf$dose.prop.ratio[3+ i*8] <- resDf$dose.prop[3+ i*8]/resDf$dose.prop[5+ i*8]
    resDf$dose.prop.ratio[4+ i*8] <- resDf$dose.prop[4+ i*8]/resDf$dose.prop[5+ i*8]
    resDf$dose.prop.ratio[5+ i*8] <- resDf$dose.prop[5+ i*8]/resDf$dose.prop[5+ i*8]
    resDf$dose.prop.ratio[6+ i*8] <- resDf$dose.prop[6+ i*8]/resDf$dose.prop[8+ i*8]
    resDf$dose.prop.ratio[7+ i*8] <- resDf$dose.prop[7+ i*8]/resDf$dose.prop[8+ i*8]
    resDf$dose.prop.ratio[8+ i*8] <- resDf$dose.prop[8+ i*8]/resDf$dose.prop[8+ i*8] 
  }
  
  
  resDf$dose.targeted.prop <- resDf$targeted.prop/resDf$dose
  return(resDf)
}
getEachCountryM <- function(rowsOfCountry) {
  mAb = 0.701
  mAbs = 0.701
  mAbh = 0.784
  mV = 0.394
  mVs = 0.394
  mVh = 0.444
  row.average <- rowsOfCountry[rowsOfCountry$YearIndex=="Average",]
  row.rest <- rowsOfCountry[rowsOfCountry$YearIndex!="Average",]
  n.year <- nrow(row.rest)
  resDf <-   expand.grid(Approach = c("A", "B", "C", "D","YR"),
                         Immunisation = c("mAb", "mV"),
                         Outcome = c("RSV-ALRI", "RSV-ALRIcwi", "RSV-ALRIHos"),
                         YearIndex = 1:n.year)
  resDf <- resDf[!(resDf$Immunisation=="mV" & (resDf$Approach=="C"|resDf$Approach=="D")),]
  resDf$prop <- NA
  resDf$dose.prop <- NA
  resDf$dose <- NA
  resDf$coverage <- NA
  resDf$targeted.prop <- NA
  resDf$dose.targeted.prop <- NA
  resDf$prop.prop <- NA
  resDf$dose.prop.ratio <- NA
  resDf$DE <- sum(as.logical(row.average[paste(month.abb, "epi2", sep = "")]))
  resDf$Approach2 <- rep(c("S1", "S2.mAb", "S3.mAb", "S4.mAb","YR", "S1", "S2.mV", "YR"), 3)
  resDf$Country <- row.average$Country
  # POP = 1000
  seasonality <- c(t(row.rest[paste("p", month.abb, sep = "")]))
  ALRI <- c(rep(row.average$est11, 3),
            rep(row.average$est12, 3),
            rep(row.average$est13, 3),
            rep(row.average$est14, 3))
  sALRI <- c(rep(row.average$est21, 3),
             rep(row.average$est22, 3),
             rep(row.average$est23, 3),
             rep(row.average$est24, 3))
  Hos <- c(rep(0.361, 3),
           rep(0.283, 3),
           rep(0.215, 3),
           rep(0.142, 3))
  matrix.ALRI <- matrix(nrow = 12 , ncol = 12*n.year, 
                        data = rep(seasonality, each = 12) * ALRI,
                        byrow = FALSE)
  matrix.sALRI <- matrix(nrow = 12, ncol = 12*n.year,
                         data = rep(seasonality, each = 12) * sALRI,
                         byrow = FALSE)
  matrix.Hos <- matrix(nrow = 12, ncol = 12*n.year,
                       data = rep(seasonality, each = 12) * Hos,
                       byrow = FALSE)
  S1 <- rep(as.logical(row.average[paste(month.abb, "epi2", sep = "")]),2)
  S2.mAb <- rep(as.logical(row.average[paste(month.abb, ".3.2", sep = "")]),2)
  S3.mAb <- rep(as.logical(row.average[paste(month.abb, ".4.2", sep = "")]),2)
  S4.mAb <- rep(as.logical(row.average[paste(month.abb, ".5.2", sep = "")]),2)
  S2.mV <- rep(as.logical(row.average[paste(month.abb, ".3.2", sep = "")]),2)
  YR <- rep(TRUE, 24)
  resDf$dose <- rep(
    c(sum(S1)/2, sum(S2.mAb)/2,sum(S3.mAb)/2, sum(S4.mAb)/2, 12,
      sum(S1)/2, sum(S2.mV)/2, 12), 3
  )
  resDf$coverage <- c(rep(row.average$mAb, 5),
                      rep(row.average$ANC4, 3),
                      rep(row.average$mAb, 5),
                      rep(row.average$ANC4, 3),
                      rep(row.average$mAb, 5),
                      rep(row.average$ANC4, 3))
  getMatrixMask <- function(n.months = 5, logic24){
    resMatrix <- matrix(nrow = 12, ncol = 12, data = FALSE)
    resMatrix[1,] <- logic24[14:25 - 1]
    resMatrix[2,] <- logic24[14:25 - 2]
    resMatrix[3,] <- logic24[14:25 - 3]
    if(n.months ==5) {
      resMatrix[4,] <- logic24[14:25 - 4]
      resMatrix[5,] <- logic24[14:25 - 5]
    }
    return(resMatrix)
  }
  for(i in 1:n.year){
    resDf$targeted.prop[1 + (i-1)*24] <- 
      sum(matrix.ALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S1)])/
      sum(matrix.ALRI/n.year) *100
    resDf$targeted.prop[2+ (i-1)*24] <- 
      sum(matrix.ALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S2.mAb)])/
      sum(matrix.ALRI/n.year) *100
    resDf$targeted.prop[3+ (i-1)*24] <- 
      sum(matrix.ALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S3.mAb)])/
      sum(matrix.ALRI/n.year) *100
    resDf$targeted.prop[4+ (i-1)*24] <- 
      sum(matrix.ALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S4.mAb)])/
      sum(matrix.ALRI/n.year) *100
    resDf$targeted.prop[5+ (i-1)*24] <- 
      sum(matrix.ALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, YR)])/
      sum(matrix.ALRI/n.year) *100
    resDf$targeted.prop[6+ (i-1)*24] <- 
      sum(matrix.ALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 3, S1)])/
      sum(matrix.ALRI/n.year) *100
    resDf$targeted.prop[7+ (i-1)*24] <- 
      sum(matrix.ALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 3, S2.mV)])/
      sum(matrix.ALRI/n.year) *100
    resDf$targeted.prop[8+ (i-1)*24] <- 
      sum(matrix.ALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 3, YR)])/
      sum(matrix.ALRI/n.year) *100
    resDf$prop[1+ (i-1)*24] <- 
      resDf$targeted.prop[1+ (i-1)*24] * row.average$mAb/100 * mAb
    resDf$prop[2+ (i-1)*24] <-
      resDf$targeted.prop[2+ (i-1)*24] * row.average$mAb/100 * mAb
    resDf$prop[3+ (i-1)*24] <- 
      resDf$targeted.prop[3+ (i-1)*24] * row.average$mAb/100 * mAb
    resDf$prop[4+ (i-1)*24] <- 
      resDf$targeted.prop[4+ (i-1)*24] * row.average$mAb/100 * mAb
    resDf$prop[5+ (i-1)*24] <- 
      resDf$targeted.prop[5+ (i-1)*24] * row.average$mAb/100 * mAb
    resDf$prop[6+ (i-1)*24] <- 
      resDf$targeted.prop[6+ (i-1)*24] * row.average$ANC4/100 * mV
    resDf$prop[7+ (i-1)*24] <- 
      resDf$targeted.prop[7+ (i-1)*24] * row.average$ANC4/100 * mV
    resDf$prop[8+ (i-1)*24] <- 
      resDf$targeted.prop[8+ (i-1)*24] * row.average$ANC4/100 * mV
    
    resDf$targeted.prop[9 + (i-1)*24] <- 
      sum(matrix.sALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S1)])/
      sum(matrix.sALRI/n.year) *100
    resDf$targeted.prop[10+ (i-1)*24] <- 
      sum(matrix.sALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S2.mAb)])/
      sum(matrix.sALRI/n.year) *100
    resDf$targeted.prop[11+ (i-1)*24] <- 
      sum(matrix.sALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S3.mAb)])/
      sum(matrix.sALRI/n.year) *100
    resDf$targeted.prop[12+ (i-1)*24] <- 
      sum(matrix.sALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S4.mAb)])/
      sum(matrix.sALRI/n.year) *100
    resDf$targeted.prop[13+ (i-1)*24] <- 
      sum(matrix.sALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, YR)])/
      sum(matrix.sALRI/n.year) *100
    resDf$targeted.prop[14+ (i-1)*24] <- 
      sum(matrix.sALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 3, S1)])/
      sum(matrix.sALRI/n.year) *100
    resDf$targeted.prop[15+ (i-1)*24] <- 
      sum(matrix.sALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 3, S2.mV)])/
      sum(matrix.sALRI/n.year) *100
    resDf$targeted.prop[16+ (i-1)*24] <- 
      sum(matrix.sALRI[,(1:12)+12*(i-1)][getMatrixMask(n.months = 3, YR)])/
      sum(matrix.sALRI/n.year) *100
    resDf$prop[9+ (i-1)*24] <- 
      resDf$targeted.prop[9+ (i-1)*24] * row.average$mAb/100 * mAbs
    resDf$prop[10+ (i-1)*24] <-
      resDf$targeted.prop[10+ (i-1)*24] * row.average$mAb/100 * mAbs
    resDf$prop[11+ (i-1)*24] <- 
      resDf$targeted.prop[11+ (i-1)*24] * row.average$mAb/100 * mAbs
    resDf$prop[12+ (i-1)*24] <-
      resDf$targeted.prop[12+ (i-1)*24] * row.average$mAb/100 * mAbs
    resDf$prop[13+ (i-1)*24] <- 
      resDf$targeted.prop[13+ (i-1)*24] * row.average$mAb/100 * mAbs
    resDf$prop[14+ (i-1)*24] <- 
      resDf$targeted.prop[14+ (i-1)*24] * row.average$ANC4/100 * mVs
    resDf$prop[15+ (i-1)*24] <- 
      resDf$targeted.prop[15+ (i-1)*24] * row.average$ANC4/100 * mVs
    resDf$prop[16+ (i-1)*24] <- 
      resDf$targeted.prop[16+ (i-1)*24] * row.average$ANC4/100 * mVs
    
    resDf$targeted.prop[17 + (i-1)*24] <- 
      sum(matrix.Hos[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S1)])/
      sum(matrix.Hos/n.year) *100
    resDf$targeted.prop[18+ (i-1)*24] <- 
      sum(matrix.Hos[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S2.mAb)])/
      sum(matrix.Hos/n.year) *100
    resDf$targeted.prop[19+ (i-1)*24] <- 
      sum(matrix.Hos[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S3.mAb)])/
      sum(matrix.Hos/n.year) *100
    resDf$targeted.prop[20+ (i-1)*24] <- 
      sum(matrix.Hos[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, S4.mAb)])/
      sum(matrix.Hos/n.year) *100
    resDf$targeted.prop[21+ (i-1)*24] <- 
      sum(matrix.Hos[,(1:12)+12*(i-1)][getMatrixMask(n.months = 5, YR)])/
      sum(matrix.Hos/n.year) *100
    resDf$targeted.prop[22+ (i-1)*24] <- 
      sum(matrix.Hos[,(1:12)+12*(i-1)][getMatrixMask(n.months = 3, S1)])/
      sum(matrix.Hos/n.year) *100
    resDf$targeted.prop[23+ (i-1)*24] <- 
      sum(matrix.Hos[,(1:12)+12*(i-1)][getMatrixMask(n.months = 3, S2.mV)])/
      sum(matrix.Hos/n.year) *100
    resDf$targeted.prop[24+ (i-1)*24] <- 
      sum(matrix.Hos[,(1:12)+12*(i-1)][getMatrixMask(n.months = 3, YR)])/
      sum(matrix.Hos/n.year) *100
    resDf$prop[17+ (i-1)*24] <- 
      resDf$targeted.prop[17+ (i-1)*24] * row.average$mAb/100 * mAbh
    resDf$prop[18+ (i-1)*24] <-
      resDf$targeted.prop[18+ (i-1)*24] * row.average$mAb/100 * mAbh
    resDf$prop[19+ (i-1)*24] <- 
      resDf$targeted.prop[19+ (i-1)*24] * row.average$mAb/100 * mAbh
    resDf$prop[20+ (i-1)*24] <-
      resDf$targeted.prop[20+ (i-1)*24] * row.average$mAb/100 * mAbh
    resDf$prop[21+ (i-1)*24] <- 
      resDf$targeted.prop[21+ (i-1)*24] * row.average$mAb/100 * mAbh
    resDf$prop[22+ (i-1)*24] <- 
      resDf$targeted.prop[22+ (i-1)*24] * row.average$ANC4/100 * mVh
    resDf$prop[23+ (i-1)*24] <- 
      resDf$targeted.prop[23+ (i-1)*24] * row.average$ANC4/100 * mVh
    resDf$prop[24+ (i-1)*24] <- 
      resDf$targeted.prop[24+ (i-1)*24] * row.average$ANC4/100 * mVh
  }
  
  
  resDf$dose.prop <- resDf$prop/resDf$dose
  for (j in 1:n.year) {
    for (i in 0:2) {
      resDf$prop.prop[1 + i*8+(j-1) *24] <- resDf$prop[1+ i*8+(j-1) *24]/resDf$prop[5+ i*8+(j-1) *24] *100
      resDf$prop.prop[2+ i*8+(j-1) *24] <- resDf$prop[2+ i*8+(j-1) *24]/resDf$prop[5+ i*8+(j-1) *24] *100
      resDf$prop.prop[3+ i*8+(j-1) *24] <- resDf$prop[3+ i*8+(j-1) *24]/resDf$prop[5+ i*8+(j-1) *24] *100
      resDf$prop.prop[4+ i*8+(j-1) *24] <- resDf$prop[4+ i*8+(j-1) *24]/resDf$prop[5+ i*8+(j-1) *24] *100
      resDf$prop.prop[5+ i*8+(j-1) *24] <- resDf$prop[5+ i*8+(j-1) *24]/resDf$prop[5+ i*8+(j-1) *24] *100
      resDf$prop.prop[6+ i*8+(j-1) *24] <- resDf$prop[6+ i*8+(j-1) *24]/resDf$prop[8+ i*8+(j-1) *24] *100
      resDf$prop.prop[7+ i*8+(j-1) *24] <- resDf$prop[7+ i*8+(j-1) *24]/resDf$prop[8+ i*8+(j-1) *24] *100
      resDf$prop.prop[8+ i*8+(j-1) *24] <- resDf$prop[8+ i*8+(j-1) *24]/resDf$prop[8+ i*8+(j-1) *24] *100
      resDf$dose.prop.ratio[1+ i*8+(j-1) *24] <- resDf$dose.prop[1+ i*8+(j-1) *24]/resDf$dose.prop[5+ i*8+(j-1) *24]
      resDf$dose.prop.ratio[2+ i*8+(j-1) *24] <- resDf$dose.prop[2+ i*8+(j-1) *24]/resDf$dose.prop[5+ i*8+(j-1) *24]
      resDf$dose.prop.ratio[3+ i*8+(j-1) *24] <- resDf$dose.prop[3+ i*8+(j-1) *24]/resDf$dose.prop[5+ i*8+(j-1) *24]
      resDf$dose.prop.ratio[4+ i*8+(j-1) *24] <- resDf$dose.prop[4+ i*8+(j-1) *24]/resDf$dose.prop[5+ i*8+(j-1) *24]
      resDf$dose.prop.ratio[5+ i*8+(j-1) *24] <- resDf$dose.prop[5+ i*8+(j-1) *24]/resDf$dose.prop[5+ i*8+(j-1) *24]
      resDf$dose.prop.ratio[6+ i*8+(j-1) *24] <- resDf$dose.prop[6+ i*8+(j-1) *24]/resDf$dose.prop[8+ i*8+(j-1) *24]
      resDf$dose.prop.ratio[7+ i*8+(j-1) *24] <- resDf$dose.prop[7+ i*8+(j-1) *24]/resDf$dose.prop[8+ i*8+(j-1) *24]
      resDf$dose.prop.ratio[8+ i*8+(j-1) *24] <- resDf$dose.prop[8+ i*8+(j-1) *24]/resDf$dose.prop[8+ i*8+(j-1) *24] 
    }
  }
  resDf$dose.targeted.prop <- resDf$targeted.prop/resDf$dose
  return(resDf)
}

getSummary <- function(Input, digit = 1) {
  return(
    paste(
      round(quantile(Input, probs = 0.5, na.rm = TRUE), digit),
      " (",
      round(quantile(Input, probs = 0.25, na.rm = TRUE), digit),
      "-",
      round(quantile(Input, probs = 0.75, na.rm = TRUE), digit),
      ")",
      sep = ""
    )
  )
}
compileRes <- function(dataset) {
  dataset <- do.call(rbind, dataset); row.names(dataset) <- NULL
  dataset$Country[dataset$Country == "Dominican Republic"] <- "Dominican Rep"
  dataset$Country_By_DE <- paste(dataset$Country, dataset$DE, sep = "_")
  dataset$Country_By_DE <- factor(dataset$Country_By_DE, 
                                  levels =unique(dataset$Country_By_DE[order(dataset$DE, 
                                                                             decreasing = FALSE)]))
  return(dataset)
}