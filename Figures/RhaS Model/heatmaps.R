library(ggplot2)
library(latex2exp) 
library(ggpubr)
library(scales)
setwd("~/Desktop/Bio_DTP_Oxford/Rotations/Antonis/PID_Rotation/Figures/RhaS_PI_ No Dilution of Zs")
data = read.csv("Eigenvalues_alphas.csv")
p <- ggplot() 
p <- p + geom_raster(data = data , aes(x = alpha_1, y = alpha_2, fill = Im.lambda_4.))
p <- p + coord_equal()
p <- p + scale_fill_gradientn(TeX("Im($\\lambda_4$)"), limits=c(-0.01, 8),
                              colours=c("green", "yellow", "orange", "red"),
                              values = rescale(c(0,
                                                 1,
                                                 4,
                                                 8)))
p <- p + ggtitle(TeX("$\\delta = 0.00039$")) +
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))


data1 = read.csv("Eigenvalues_alphas_0.1xDeg.csv")
p1 <- ggplot() 
p1 <- p1 + geom_raster(data = data1 , aes(x = alpha_1, y = alpha_2, fill = Im.lambda_4.))
p1 <- p1 + coord_equal()
p1 <- p1 + scale_fill_gradientn(TeX("Im($\\lambda_4$)"), limits=c(-0.01, 8),
                              colours=c("green", "yellow", "orange", "red"),
                              values = rescale(c(0,
                                                 1,
                                                 4,
                                                 8)))
p1 <- p1 + ggtitle(TeX("$\\delta = 0.000039$")) + 
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))
                              
data2 = read.csv("Eigenvalues_alphas_0.5xDeg.csv")
p2<- ggplot() 
p2 <- p2 + geom_raster(data = data2 , aes(x = alpha_1, y = alpha_2, fill = Im.lambda_4.))
p2<- p2 + coord_equal()
p2 <- p2 + scale_fill_gradientn(TeX("Im($\\lambda_4$)"), limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(0,
                                                   1,
                                                   4,
                                                   8)))
p2 <- p2 + ggtitle(TeX("$\\delta = 0.000195$")) + 
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))



data3 = read.csv("Eigenvalues_alphas_2xDeg.csv")
p3 <- ggplot() 
p3 <- p3 + geom_raster(data = data3 , aes(x = alpha_1, y = alpha_2, fill = Im.lambda_4.))
p3 <- p3 + coord_equal()
p3 <- p3 + scale_fill_gradientn(TeX("Im($\\lambda_4$)"), limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(0,
                                                   1,
                                                   4,
                                                   8)))
p3 <- p3 + ggtitle(TeX("$\\delta = 0.00078$")) + 
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))


data4 = read.csv("Eigenvalues_alphas_10xDeg.csv")
p4 <- ggplot() 
p4 <- p4 + geom_raster(data = data4 , aes(x = alpha_1, y = alpha_2, fill = Im.lambda_4.))
p4 <- p4 + coord_equal()
p4 <- p4 + scale_fill_gradientn(TeX("Im($\\lambda_4$)"), limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(0,
                                                   1,
                                                   4,
                                                   8)))
p4 <- p4 + ggtitle(TeX("$\\delta = 0.0039$")) + 
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))


data5 = read.csv("Eigenvalues_alphas_20xDeg.csv")
p5 <- ggplot() 
p5 <- p5 + geom_raster(data = data5 , aes(x = alpha_1, y = alpha_2, fill = Im.lambda_4.))
p5 <- p5 + coord_equal()
p5 <- p5 + scale_fill_gradientn(TeX("Im($\\lambda_4$)"), limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(0,
                                                   1,
                                                   4,
                                                   8)))
p5 <- p5 + ggtitle(TeX("$\\delta = 0.0078$")) + 
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))

F = ggarrange(p, p1, p2, p3, p4, p5, 
              ncol = 2, nrow = 3)                                
                                                                                            

annotate_figure(F, 
                TeX("Im($\\lambda_4$) for RhaS Model"),
                fig.lab.pos = "top",
                fig.lab.face = "bold",
                fig.lab.size = 18)
