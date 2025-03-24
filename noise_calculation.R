################################################################################
#############          Pelagic Synthesis           #############################
#############             MAR-2025                 #############################
#############          Double Integration          #############################
## by: Alexandra Cabanelas 
################################################################################

## ------------------------------------------ ##
#            Packages -----
## ------------------------------------------ ##
library(tidyverse)
library(astsa)

## ------------------------------------------ ##
#            Data -----
## ------------------------------------------ ##

# PDO
PDO <- read.csv(file.path("output",
                             "CCE",
                             "PDO_CCE_integrated.csv")) %>%
  select(time, pdo, pdoInt) %>%
  filter(!is.na(pdoInt))

# MEI
MEI <- read.csv(file.path("output",
                          "PAL",
                          "MEI_PAL_integrated.csv")) %>%
  select(time, mei, meiInt) %>%
  filter(!is.na(meiInt))

## ------------------------------------------ ##
#            TIDY -----
## ------------------------------------------ ##
#make ts object
PDOts <- ts(PDO$pdo, start = c(1941, 2), frequency = 12)
pdoIntts <- ts(PDO$pdoInt, start = c(1941, 2), frequency = 12)

MEIts <- ts(MEI$mei, start = c(1979, 2), frequency = 12)
meiIntts <- ts(MEI$meiInt, start = c(1979, 2), frequency = 12)

## ------------------------------------------ ##
#           PDO  -----
## ------------------------------------------ ##
specPDO <- spectrum(PDOts, log = "no")  
plot(log10(specPDO$freq), log10(specPDO$spec), type = "l",
     xlab = "log10(Frequency)", ylab = "log10(Variance)",
     main = "PDO Spectrum")

# get slope
fit <- lm(log10(specPDO$spec) ~ log10(specPDO$freq))
abline(fit, col = "blue", lwd = 2)
# slope
legend("bottomleft", legend = paste0("Slope: ", 
                                   round(coef(fit)[2], 2)), bty = "n")
specPDO_plot <- recordPlot()

specPDO1 <- spec.pgram(PDOts, taper = 0.1, log = "yes")
plot(specPDO1$freq, specPDO1$spec, type = "l", 
     xlab = "Frequency", 
     ylab = "Variance")

specPDO1 <- mvspec(PDOts, log = "no", plot = TRUE)


## PDO INT
specInt <- spectrum(pdoIntts, log = "no")  
plot(log10(specInt$freq), log10(specInt$spec), type = "l",
     xlab = "log10(Frequency)", ylab = "log10(Variance)",
     main = "Integrated PDO Spectrum")

# get slope
fitInt <- lm(log10(specInt$spec) ~ log10(specInt$freq))
abline(fitInt, col = "blue", lwd = 2)
# slope
legend("bottomleft", legend = paste0("Slope: ", 
                                   round(coef(fitInt)[2], 2)), 
       bty = "n")
specInt_plot <- recordPlot()


## ------------------------------------------ ##
#           MEI  -----
## ------------------------------------------ ##
specMEI <- spectrum(MEIts, log = "no")  
plot(log10(specMEI$freq), log10(specMEI$spec), type = "l",
     xlab = "log10(Frequency)", ylab = "log10(Variance)",
     main = "MEI Spectrum")

# get slope
fitMEI <- lm(log10(specMEI$spec) ~ log10(specMEI$freq))
abline(fitMEI, col = "blue", lwd = 2)
# slope
legend("bottomleft", legend = paste0("Slope: ", 
                                   round(coef(fitMEI)[2], 2)), 
       bty = "n")
specMEI_plot <- recordPlot()


## MEI INT
specMeiInt <- spectrum(meiIntts, log = "no")  
plot(log10(specMeiInt$freq), log10(specMeiInt$spec), type = "l",
     xlab = "log10(Frequency)", ylab = "log10(Variance)",
     main = "Integrated MEI Spectrum")

# get slope
fitMeiInt <- lm(log10(specMeiInt$spec) ~ log10(specMeiInt$freq))
abline(fitMeiInt, col = "blue", lwd = 2)
# slope
legend("bottomleft", legend = paste0("Slope: ", 
                                   round(coef(fitMeiInt)[2], 2)), 
       bty = "n")
specMeiInt_plot <- recordPlot()




specPDO_plot
specInt_plot
specMEI_plot
specMeiInt_plot

plots <- list(
  specPDO_plot = specPDO_plot,
  specInt_plot = specInt_plot,
  specMEI_plot = specMEI_plot,
  specMeiInt_plot = specMeiInt_plot
)

for (name in names(plots)) {
  png(filename = file.path("figures", paste0(name, ".png")),
      width = 6, height = 4.5, units = "in", res = 300)
  replayPlot(plots[[name]])
  dev.off()
}
