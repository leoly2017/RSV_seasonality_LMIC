library(plyr)
library(dplyr)
library(ggplot2)
library(maps)
library(mapdata)
library(rworldmap)
library(rgeos)
library(geosphere)
drawHeat <- function(dataset, filename = "heatAll") {
  dataset <- dataset[c("Region", "Country", "Latitude", paste("p", month.abb, sep = ""))]
  len <- nrow(dataset)
  dataset$RC <- paste(dataset$Region, ", ", dataset$Country, sep = "")
  dataset <- dataset[order(dataset$Latitude),]
  dataset$RC <- factor(dataset$RC, levels = unique(dataset$RC))
  dataset2 <- gather(dataset, pJan, pFeb, pMar, pApr, pMay, pJun,
                     pJul, pAug, pSep, pOct, pNov, pDec, value = "AAP",
                     key = "month")
  dataset2$month <- factor(dataset2$month, levels = paste("p", month.abb, sep = ""),
                           labels = month.abb)
  dataset2$AAP <- dataset2$AAP * 100
  plotRes <-   ggplot(dataset2, aes(month, RC)) + 
    geom_tile(aes(fill = AAP), color = "white") + 
    scale_fill_gradient2(low = "deepskyblue", mid = "white", high = "red3",
                         limits = c(0, 40), midpoint = 10, na.value = "red3") +
    scale_x_discrete(labels = 1:12) + 
    geom_hline(yintercept = sum(dataset$Latitude <0) + 0.5)+
    geom_hline(yintercept = sum(dataset$Latitude < -23.5) + 0.5, linetype = "dashed")+
    geom_hline(yintercept = sum(dataset$Latitude < 23.5) + 0.5, linetype = "dashed")+
    theme(text = element_text(size = 10), axis.title.y = element_blank(),
          legend.key.size = unit(0.5, "inches")) + 
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20))
  ggsave(plotRes, filename = paste("plots/",filename,".jpg", sep = ""), width = 8, height = len/4)
}
drawHeat2 <- function(dataset, filename = "heatAll") {
  dataset <- dataset[c("Country", paste("p", month.abb, sep = ""))]
  centroids <- read_excel("dataset/Country_coordinates.xlsx")
  centroids$Country[centroids$Country=="Vietnam"] <- "Viet Nam"
  dataset <- left_join(dataset, centroids)
  len <- nrow(dataset)
  dataset$RC <- paste(dataset$Country, sep = "")
  dataset <- dataset[order(dataset$y),]
  dataset$RC <- factor(dataset$RC, levels = unique(dataset$RC))
  dataset2 <- gather(dataset, pJan, pFeb, pMar, pApr, pMay, pJun,
                     pJul, pAug, pSep, pOct, pNov, pDec, value = "AAP",
                     key = "month")
  dataset2$month <- factor(dataset2$month, levels = paste("p", month.abb, sep = ""),
                           labels = month.abb)
  dataset2$AAP <- dataset2$AAP * 100
  plotRes <-   ggplot(dataset2, aes(month, RC)) + 
    geom_tile(aes(fill = AAP), color = "white") + 
    scale_fill_gradient2(low = "white", mid = "white", high = "red3",
                         limits = c(0, 40), midpoint = 10, na.value = "red3") +
    scale_x_discrete(labels = 1:12) + 
    geom_hline(yintercept = sum(dataset$y <0) + 0.5)+
    geom_hline(yintercept = sum(dataset$y < -23.44) + 0.5, linetype = "dashed")+
    geom_hline(yintercept = sum(dataset$y < 23.44) + 0.5, linetype = "dashed")+
    theme(text = element_text(size = 10), axis.title.y = element_blank(),
          legend.key.size = unit(0.4, "inches")) + 
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15),
          axis.text.x = element_text(size = 15), axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 15))
  ggsave(plotRes, filename = paste("figures/Seasonality/",filename,".pdf", sep = ""), width = 8, height = len/4)
}
drawHeatMultiyear2 <- function(dataset) {
  dataset <- dataset[c("Region", "Country", "epi2", paste("p", month.abb, sep = ""), "YearIndex")]
  len <- nrow(dataset)
  dataset0 <- dataset[dataset$YearIndex == "Average", c("Country", "epi2")]
  names(dataset0)[2] <- "epi20"
  dataset <- left_join(dataset, dataset0)
  dataset$RC <- paste(dataset$Country, dataset$epi20,sep = "_")
  dataset <- dataset[order(dataset$epi20),]
  dataset$RC <- factor(dataset$RC, levels = unique(dataset$RC))
  dataset$YearIndex <- factor(dataset$YearIndex, levels = rev(c("Average", 1:10)))
  dataset2 <- gather(dataset, pJan, pFeb, pMar, pApr, pMay, pJun,
                     pJul, pAug, pSep, pOct, pNov, pDec, value = "AAP",
                     key = "month")
  dataset2$month <- factor(dataset2$month, levels = paste("p", month.abb, sep = ""),
                           labels = month.abb)
  dataset2$AAP <- dataset2$AAP * 100
  plotRes <-   ggplot(dataset2, aes(month, YearIndex)) + 
    geom_tile(aes(fill = AAP), color = "white") + 
    scale_fill_gradient2(low = "white", mid = "white", high = "red3",
                         limits = c(0, 40), midpoint = 10, na.value = "red3",
                         name = "AP") +
    scale_x_discrete(labels = 1:12) + 
    facet_wrap(~RC, ncol = 5) +
  theme(text = element_text(size = 15), axis.title.y = element_blank(),
        legend.key.size = unit(0.3, "inches")) + 
    theme(legend.text = element_text(size = 15), legend.title = element_text(size = 15),
          axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 15))
  ggsave(plotRes, filename = "figures/Seasonality/heatMult.pdf", width = 12, height = 9)
}
drawHeatMultiyear <- function(dataset) {
  dataset <- dataset[c("Region", "Country", paste("p", month.abb, sep = ""), "YearIndex")]
  len <- nrow(dataset)
  dataset$RC <- paste(dataset$Region, ", ", dataset$Country, sep = "")
  dataset$RC <- factor(dataset$RC, levels = rev(unique(dataset$RC)))
  dataset$YearIndex <- factor(dataset$YearIndex, levels = rev(c("Average", 1:10)))
  dataset2 <- gather(dataset, pJan, pFeb, pMar, pApr, pMay, pJun,
                     pJul, pAug, pSep, pOct, pNov, pDec, value = "AAP",
                     key = "month")
  dataset2$month <- factor(dataset2$month, levels = paste("p", month.abb, sep = ""),
                           labels = month.abb)
  dataset2$AAP <- dataset2$AAP * 100
  plotRes <-   ggplot(dataset2, aes(month, YearIndex)) + 
    geom_tile(aes(fill = AAP), color = "white") + 
    scale_fill_gradient2(low = "deepskyblue", mid = "white", high = "red3",
                         limits = c(0, 40), midpoint = 10, na.value = "red3",
                         name = "AP") +
    scale_x_discrete(labels = 1:12) + 
    facet_wrap(~RC, nrow = 9)
  theme(text = element_text(size = 15), axis.title.y = element_blank(),
        legend.key.size = unit(0.5, "inches")) + 
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 20),
          axis.text.x = element_text(size = 20), axis.title.x = element_text(size = 20))
  ggsave(plotRes, filename = "plots/heatMultiyear.jpg", width = 8, height = 16)
}

drawAnySite <- function(dataset, siteColor, filename, var) {
  tiff(paste("plots/",filename, ".tif", sep = ""), width = 1000, height = 500)
  map("worldHires", bg = "aliceblue", fill = TRUE, col = "white", ylim = c(-60,90))
  abline(h = 0, lty = 1)
  abline(h = 23.5, lty = 2)
  abline(h = -23.5, lty = 2)
  abline(h = 35, lty = 3)
  abline(h = -35, lty = 3)
  points(dataset$Longitude, dataset$Latitude, col = siteColor[var], 
         pch = 19)
  legend(-180,-10, 
         unique(var[order(var)]), 
         pch = 19, col = siteColor[unique(var[order(var)])], title = NULL,  
         bg = "white")
  text(x = rep(185,5),
       y = c(-32, -20.5, 3, 26.5, 38),
       c("-35", "-23.5", "0", "23.5", "35"),
       cex = 1)
  dev.off()
}

drawHist <- function(dataset) {
  mu <- ddply(dataset, "latgp2", summarise, grp.mean=mean(epi2))
  ggsave(
    ggplot(data = LMICAll, aes(x = epi2)) +
      geom_histogram(aes(y=..density..),binwidth = 1, colour = "black", fill = "white") +
      geom_density(alpha=.2,adjust = 1,fill="#FF6666") +
      geom_vline(data = mu, aes(xintercept = grp.mean, colour = "red"), linetype = "dashed",
                 show.legend = FALSE) + 
      scale_x_continuous(limits = c(1,9), breaks = c(1:9)) +
      facet_wrap(~latgp2, nrow = 2) + 
      theme_bw() +
      labs(x = "duration of epidemics"),
    filename = "plots/hist.jpg",
    height = 8, width = 4
  )
  
}

genMCI <- function(vector) {
  c(
    Mean = round(mean(vector),2),
    LCI = round(mean(vector)-1.96*sd(vector)/sqrt(length(vector)),2),
    UCI = round(mean(vector)+1.96*sd(vector)/sqrt(length(vector)),2)
  )
}

genEpiTable <- function(dataset) {
  res <- dataset[c("WHORegion", "Country", "Region", paste(month.abb, "epi2", sep = ""))]
  names(res) <- c(
    "WHORegion", "Country", "Region", month.abb
  )
  return(res[order(res["WHORegion"], res["Country"]),])
}

exportForPart2 <- function() {
  centroids <- read_excel("dataset/Country_coordinates.xlsx")
  LMICAllt <- left_join(LMICAll, centroids)
  LMICAllt$dist <- distHaversine(p1 = as.matrix(LMICAllt[c("Longitude", "Latitude")]),
                                 p2 = as.matrix(LMICAllt[c("x", "y")]))/1000
  distMin <- LMICAllt %>% 
    group_by(Country) %>%
    dplyr::summarise(dist = min(dist))
  LMICAllt <- left_join(distMin, LMICAllt)
  write.csv(LMICAll[c("Dataset ID","Region", "Country", "Latitude", "Longitude", paste(
    "p", month.abb, sep = ""
  ), paste(
    month.abb, "epi2",sep = ""
  ), "epi2")], file = "dataset/AllSites.csv", row.names = FALSE)
  write.csv(LMICAllt[c("Dataset ID","Region", "Country", "Latitude", "Longitude", paste(
    "p", month.abb, sep = ""
  ), paste(
    month.abb, "epi2",sep = ""
  ), "epi2", "nPeak")], file = "dataset/CountryLevel.csv", row.names = FALSE)
  
  LMICAllt <- left_join(LMICMultiyear, centroids)
  LMICAllt$dist <- distHaversine(p1 = as.matrix(LMICAllt[c("Longitude", "Latitude")]),
                                 p2 = as.matrix(LMICAllt[c("x", "y")]))/1000
  distMin <- LMICAllt %>% 
    group_by(Country) %>%
    dplyr::summarise(dist = min(dist))
  LMICAllt <- left_join(distMin, LMICAllt)
  write.csv(LMICAllt[c("Dataset ID","YearIndex", "Region", "Country", "Latitude", "Longitude", paste(
    "p", month.abb, sep = ""
  ), paste(
    month.abb, "epi2",sep = ""
  ), "epi2", "nPeak")], file = "dataset/CountryLevelMulti.csv", row.names = FALSE)
}

genDosingSchedule <- function() {
  dataset <- MainData[c("Country", paste(month.abb, "epi2", sep = ""),
                        paste(month.abb, ".3.2", sep = ""),
                        paste(month.abb, ".4.2", sep = ""),
                        paste(month.abb, ".5.2", sep = ""))]
  centroids <- read_excel("dataset/Country_coordinates.xlsx")
  centroids$Country[centroids$Country=="Vietnam"] <- "Viet Nam"
  dataset <- left_join(dataset, centroids)
  len <- nrow(dataset)
  dataset$RC <- paste(dataset$Country, sep = "")
  dataset <- dataset[order(dataset$y),]
  dataset$RC <- factor(dataset$RC, levels = unique(dataset$RC))
  datasetA <- gather(dataset[c("RC", paste(month.abb, "epi2", sep = ""))],
                     -RC,
                     value = "Admin",
                     key = "month")
  datasetA$month <- substr(datasetA$month, 1,3)
  datasetA$Approach <- "A" 
  datasetB <- gather(dataset[c("RC", paste(month.abb, ".3.2", sep = ""))],
                     -RC,
                     value = "Admin",
                     key = "month")
  datasetB$month <- substr(datasetB$month, 1,3)
  datasetB$Approach <- "B" 
  datasetC <- gather(dataset[c("RC", paste(month.abb, ".4.2", sep = ""))],
                     -RC,
                     value = "Admin",
                     key = "month")
  datasetC$month <- substr(datasetC$month, 1,3)
  datasetC$Approach <- "C" 
  datasetD <- gather(dataset[c("RC", paste(month.abb, ".5.2", sep = ""))],
                     -RC,
                     value = "Admin",
                     key = "month")
  datasetD$month <- substr(datasetD$month, 1,3)
  datasetD$Approach <- "D" 
  dataset2 <- Reduce(rbind, list(datasetA, datasetB, datasetC, datasetD))
  dataset2$month <- factor(dataset2$month, levels = month.abb)
  plotRes <- ggplot(dataset2, aes(month, RC)) + 
    geom_tile(aes(fill = Admin), color = "grey") + 
    scale_fill_manual(values = c("white","deepskyblue3")) +
    scale_x_discrete(labels = 1:12) +
    facet_wrap(~Approach, ncol = 4) +
    theme(text = element_text(size = 10), axis.title.y = element_blank(),
          legend.position = "none")+ 
    theme(axis.text.x = element_text(size = 10), 
          axis.title.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  ggsave(plotRes, filename = "figures/mAb_schedule.pdf", height = 8, width = 10)
  plotRes <- ggplot(dataset2[dataset2$Approach %in% c("A", "B"),], aes(month, RC)) + 
    geom_tile(aes(fill = Admin), color = "grey") + 
    scale_fill_manual(values = c("white","deepskyblue3")) +
    scale_x_discrete(labels = 1:12) +
    facet_wrap(~Approach, ncol = 4) +
    theme(text = element_text(size = 10), axis.title.y = element_blank(),
          legend.position = "none")+ 
    theme(axis.text.x = element_text(size = 10), 
          axis.title.x = element_text(size = 12),
          axis.text.y = element_text(size = 12))
  ggsave(plotRes, filename = "figures/mV_schedule.pdf", height = 8, width = 6)
}