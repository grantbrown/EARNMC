makeKikwitGraphPDFs = function(fullSpecResultsFilename=NA,
                               underspecResultsFilename=NA){
    if (is.na(fullSpecResultsFilename)){
        print("No full specification results file name provided, using included result data.")
        data("KikwitFullSpecResults")
    }
    else{
        load(fullSpecResultsFilename)
        KikwitFullSpecResults = outData
    }
    if (is.na(underspecResultsFilename)){
        print("No underspecification results file name provided, using included result data.")
        data("KikwitUnderspecResults")
    }
    else{
        load(underspecResultsFilename)
        KikwitUnderspecResults = outData
    }

    CI_width = 2
    MainLineWidth = 3
    ThresholdWidth = 2


    r01 = KikwitUnderspecResults$R0Estimates[[1]][[1]]
    r02 = KikwitUnderspecResults$R0Estimates[[1]][[2]]

    min.r0 = min(min(KikwitUnderspecResults$R0Estimates[[1]][[1]]$LB), min(KikwitUnderspecResults$R0Estimates[[1]][[2]]$LB))
    max.r0 = max(max(KikwitUnderspecResults$R0Estimates[[1]][[1]]$UB), max(KikwitUnderspecResults$R0Estimates[[1]][[2]]$UB) + 0.1)
    r0.ylim = c(min.r0, max.r0)

    pdf(file="./EA_RO_Comparison.pdf", width = 12, height = 8)
        #par(mfrow = c(2,1))
        layout(matrix(c(1,2,3,2), 2, 2, byrow = TRUE))
        


        plot(r01$mean[1:(length(r01$mean) - 1)], ylim = r0.ylim, main = "Traditional R0: Underspecified Intensity\n Posterior Mode and 90% CI",
             xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
        abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
        abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
        lines(1:(length(r01$mean) - 1), r01$LB[1:(length(r01$LB)-1)], lty = 2, lwd = CI_width)
        lines(1:(length(r01$mean) - 1), r01$UB[1:(length(r01$UB)-1)], lty = 2, lwd = CI_width)

        barplot(t(KikwitUnderspecResults$simResults$I_star), main = "New Cases", xlab = "Day", ylab = "Cases")
        axis(side = 1, at = seq(0, (length(r01$mean)), 50))

        plot(r02$mean[1:(length(r02$mean) - 1)], ylim = r0.ylim, main = "EA-R: Underspecified Intensity\n Posterior Mode and 90% CI",
             xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
        abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
        abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
        lines(1:(length(r02$mean) - 1), r02$LB[1:(length(r02$LB)-1)], lty = 2, lwd = CI_width)
        lines(1:(length(r02$mean) - 1), r02$UB[1:(length(r02$UB)-1)], lty = 2, lwd = CI_width)
    dev.off()

    pdf(file="./EA_RO_Comparison_nodata.pdf", width = 12, height = 12)
        #layout(matrix(c(1,2,3,2), 2, 2, byrow = TRUE))
        par(mfrow = c(2,1))
       


        plot(r01$mean[1:(length(r01$mean) - 1)], ylim = r0.ylim, main = "Traditional R0: Underspecified Intensity\n Posterior Mode and 90% CI",
             xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
        abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
        abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
        lines(1:(length(r01$mean) - 1), r01$LB[1:(length(r01$LB)-1)], lty = 2, lwd = CI_width)
        lines(1:(length(r01$mean) - 1), r01$UB[1:(length(r01$UB)-1)], lty = 2, lwd = CI_width)


        plot(r02$mean[1:(length(r02$mean) - 1)], ylim = r0.ylim, main = "EA-R: Underspecified Intensity\n Posterior Mode and 90% CI",
             xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
        abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
        abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
        lines(1:(length(r02$mean) - 1), r02$LB[1:(length(r02$LB)-1)], lty = 2, lwd = CI_width)
        lines(1:(length(r02$mean) - 1), r02$UB[1:(length(r02$UB)-1)], lty = 2, lwd = CI_width)
    dev.off()

    png(filename="./EA_RO_Comparison_nodata.png", width = 800, height = 600)
        #layout(matrix(c(1,2,3,2), 2, 2, byrow = TRUE))
        par(mfrow = c(2,1))
       


        plot(r01$mean[1:(length(r01$mean) - 1)], ylim = r0.ylim, main = "Traditional R0: Underspecified Intensity\n Posterior Mode and 90% CI",
             xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
        abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
        abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
        lines(1:(length(r01$mean) - 1), r01$LB[1:(length(r01$LB)-1)], lty = 2, lwd = CI_width)
        lines(1:(length(r01$mean) - 1), r01$UB[1:(length(r01$UB)-1)], lty = 2, lwd = CI_width)


        plot(r02$mean[1:(length(r02$mean) - 1)], ylim = r0.ylim, main = "EA-R: Underspecified Intensity\n Posterior Mode and 90% CI",
             xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
        abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
        abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
        lines(1:(length(r02$mean) - 1), r02$LB[1:(length(r02$LB)-1)], lty = 2, lwd = CI_width)
        lines(1:(length(r02$mean) - 1), r02$UB[1:(length(r02$UB)-1)], lty = 2, lwd = CI_width)
    dev.off()


    r01 = KikwitFullSpecResults$R0Estimates[[1]][[1]]
    r02 = KikwitFullSpecResults$R0Estimates[[1]][[2]]

    min.r0 = min(min(KikwitFullSpecResults$R0Estimates[[1]][[1]]$LB), min(KikwitFullSpecResults$R0Estimates[[1]][[2]]$LB))
    max.r0 = max(max(KikwitFullSpecResults$R0Estimates[[1]][[1]]$UB), max(KikwitFullSpecResults$R0Estimates[[1]][[2]]$UB) + 0.1)
    r0.ylim = c(min.r0, max.r0)

    pdf(file="EA_R0_Comparison2.pdf", width = 12, height = 12)
        par(mfrow = c(2,1))
        plot(r01$mean[1:(length(r01$mean) - 1)], ylim = r0.ylim, main = "Traditional R0\n Posterior Mode and 90% CI",
             xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
        abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
        abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
        lines(1:(length(r01$mean) - 1), r01$LB[1:(length(r01$LB)-1)], lty = 2, lwd = CI_width)
        lines(1:(length(r01$mean) - 1), r01$UB[1:(length(r01$UB)-1)], lty = 2, lwd = CI_width)

        plot(r02$mean[1:(length(r02$mean) - 1)], ylim = r0.ylim, main = "EA-R\n Posterior Mode and 90% CI",
             xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
        abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
        abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
        lines(1:(length(r02$mean) - 1), r02$LB[1:(length(r02$LB)-1)], lty = 2, lwd = CI_width)
        lines(1:(length(r02$mean) - 1), r02$UB[1:(length(r02$UB)-1)], lty=2,lwd = CI_width)


    dev.off()

    png(filename="EA_R0_Comparison2.png", width = 800, height = 600)
        par(mfrow = c(2,1))
        plot(r01$mean[1:(length(r01$mean) - 1)], ylim = r0.ylim, main = "Traditional R0\n Posterior Mode and 90% CI",
             xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
        abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
        abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
        lines(1:(length(r01$mean) - 1), r01$LB[1:(length(r01$LB)-1)], lty = 2, lwd = CI_width)
        lines(1:(length(r01$mean) - 1), r01$UB[1:(length(r01$UB)-1)], lty = 2, lwd = CI_width)

        plot(r02$mean[1:(length(r02$mean) - 1)], ylim = r0.ylim, main = "EA-R\n Posterior Mode and 90% CI",
             xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
        abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
        abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
        lines(1:(length(r02$mean) - 1), r02$LB[1:(length(r02$LB)-1)], lty = 2, lwd = CI_width)
        lines(1:(length(r02$mean) - 1), r02$UB[1:(length(r02$UB)-1)], lty=2,lwd = CI_width)
    dev.off()

    pdf(file="EA_R0_Comparison_Combined.pdf", width=12, height=6)
    layout(matrix(c(1,2,3,
                    4,2,5), nrow = 2, ncol = 3, byrow = TRUE))
    
    # Plot 1
    r01 = KikwitUnderspecResults$R0Estimates[[1]][[1]]
    r02 = KikwitUnderspecResults$R0Estimates[[1]][[2]]

    min.r0 = min(min(KikwitUnderspecResults$R0Estimates[[1]][[1]]$LB), min(KikwitUnderspecResults$R0Estimates[[1]][[2]]$LB))
    max.r0 = max(max(KikwitUnderspecResults$R0Estimates[[1]][[1]]$UB), max(KikwitUnderspecResults$R0Estimates[[1]][[2]]$UB) + 0.1)
    r0.ylim = c(min.r0, max.r0)

    plot(r01$mean[1:(length(r01$mean) - 1)], ylim = r0.ylim, main = "Traditional R0: Underspecified Intensity\n Posterior Mode and 90% CI",
         xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
    abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
    abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
    lines(1:(length(r01$mean) - 1), r01$LB[1:(length(r01$LB)-1)], lty = 2, lwd = CI_width)
    lines(1:(length(r01$mean) - 1), r01$UB[1:(length(r01$UB)-1)], lty = 2, lwd = CI_width)

    # Plot 2
    barplot(t(KikwitUnderspecResults$simResults$I_star), main = "New Cases", xlab = "Day", ylab = "Cases")
    axis(side = 1, at = seq(0, (length(r01$mean)), 50))

    # Plot 3
    r01 = KikwitFullSpecResults$R0Estimates[[1]][[1]]
    r02 = KikwitFullSpecResults$R0Estimates[[1]][[2]]

    min.r0 = min(min(KikwitFullSpecResults$R0Estimates[[1]][[1]]$LB), min(KikwitFullSpecResults$R0Estimates[[1]][[2]]$LB))
    max.r0 = max(max(KikwitFullSpecResults$R0Estimates[[1]][[1]]$UB), max(KikwitFullSpecResults$R0Estimates[[1]][[2]]$UB) + 0.1)
    r0.ylim = c(min.r0, max.r0)

    plot(r01$mean[1:(length(r01$mean) - 1)], ylim = r0.ylim, main = "Traditional R0\n Posterior Mode and 90% CI",
         xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
    abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
    abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
    lines(1:(length(r01$mean) - 1), r01$LB[1:(length(r01$LB)-1)], lty = 2, lwd = CI_width)
    lines(1:(length(r01$mean) - 1), r01$UB[1:(length(r01$UB)-1)], lty = 2, lwd = CI_width)

    # Plot 4
    r01 = KikwitUnderspecResults$R0Estimates[[1]][[1]]
    r02 = KikwitUnderspecResults$R0Estimates[[1]][[2]]

    min.r0 = min(min(KikwitUnderspecResults$R0Estimates[[1]][[1]]$LB), min(KikwitUnderspecResults$R0Estimates[[1]][[2]]$LB))
    max.r0 = max(max(KikwitUnderspecResults$R0Estimates[[1]][[1]]$UB), max(KikwitUnderspecResults$R0Estimates[[1]][[2]]$UB) + 0.1)
    r0.ylim = c(min.r0, max.r0)

    plot(r02$mean[1:(length(r02$mean) - 1)], ylim = r0.ylim, main = "EA-R: Underspecified Intensity\n Posterior Mode and 90% CI",
         xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
    abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
    abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
    lines(1:(length(r02$mean) - 1), r02$LB[1:(length(r02$LB)-1)], lty = 2, lwd = CI_width)
    lines(1:(length(r02$mean) - 1), r02$UB[1:(length(r02$UB)-1)], lty = 2, lwd = CI_width)

    # Plot 5
    r01 = KikwitFullSpecResults$R0Estimates[[1]][[1]]
    r02 = KikwitFullSpecResults$R0Estimates[[1]][[2]]

    min.r0 = min(min(KikwitFullSpecResults$R0Estimates[[1]][[1]]$LB), min(KikwitFullSpecResults$R0Estimates[[1]][[2]]$LB))
    max.r0 = max(max(KikwitFullSpecResults$R0Estimates[[1]][[1]]$UB), max(KikwitFullSpecResults$R0Estimates[[1]][[2]]$UB) + 0.1)
    r0.ylim = c(min.r0, max.r0)

    plot(r02$mean[1:(length(r02$mean) - 1)], ylim = r0.ylim, main = "EA-R\n Posterior Mode and 90% CI",
         xlab = "Day", ylab = "Reproductive Number", type = "l", lwd = MainLineWidth)
    abline(h = seq(0,5,0.2), lty = 3,col="lightgrey")
    abline(h = 1.0, col = "blue", lwd = ThresholdWidth, lty = 2)
    lines(1:(length(r02$mean) - 1), r02$LB[1:(length(r02$LB)-1)], lty = 2, lwd = CI_width)
    lines(1:(length(r02$mean) - 1), r02$UB[1:(length(r02$UB)-1)], lty=2,lwd = CI_width)
    dev.off()

}



