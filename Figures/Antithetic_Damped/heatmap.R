library(ggplot2)
library(latex2exp) 
library(ggpubr)
library(scales)
data = read.csv("eigenvalues_alphas.csv")
data$alpha_1 = data$alpha_1/60/60
data$alpha_2 = data$alpha_2/60/60
p <- ggplot() 
p <- p + geom_raster(data = data , aes(x = alpha_1, y = alpha_2, fill = Imag))
p <- p + coord_equal()
p <- p + scale_fill_gradientn(TeX("Im($\\lambda_3$)"), limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(min(data$Imag),
                                                   1,
                                                   4,
                                                   max(data$Imag))))
p <- p + ggtitle(TeX("$\\delta_1 = 0.00039$ and   $\\delta_2 = 0.00039$")) +
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))


data1 = read.csv("eigenvalues_alphas_10xDegStar.csv")
data1$alpha_1 = data1$alpha_1/60/60
data1$alpha_2 = data1$alpha_2/60/60
p1 <- ggplot() 
p1 <- p1 + geom_raster(data = data1 , aes(x = alpha_1, y = alpha_2, fill = imaginary.part))
p1 <- p1 + coord_equal()
p1 <- p1 + scale_fill_gradientn(TeX("Im($\\lambda_3$)"),limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(min(data1$imaginary.part),
                                                   1,
                                                   4,
                                                   max(data1$imaginary.part))))
p1 <- p1 + ggtitle(TeX("$\\delta_1 = 0.00039$ and   $\\delta_2 = 0.0039$")) +
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))

data2 = read.csv("eigenvalues_alphas_0xDegStar.csv")
data2$alpha_1 = data2$alpha_1/60/60
data2$alpha_2 = data2$alpha_2/60/60
p2 <- ggplot() 
p2 <- p2 + geom_raster(data = data2 , aes(x = alpha_1, y = alpha_2, fill = imaginary.part))
p2 <- p2 + coord_equal()
p2 <- p2 + scale_fill_gradientn(TeX("Im($\\lambda_3$)"),limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(min(data2$imaginary.part),
                                                   1,
                                                   4,
                                                   max(data2$imaginary.part))))
p2 <- p2 + ggtitle(TeX("$\\delta_1 = 0.00039$ and   $\\delta_2 = 0$")) +
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))



data3 = read.csv("eigenvalues_alphas_0xDeg.csv")
data3$alpha_1 = data3$alpha_1/60/60
data3$alpha_2 = data3$alpha_2/60/60
p3 <- ggplot() 
p3 <- p3 + geom_raster(data = data3 , aes(x = alpha_1, y = alpha_2, fill = Imag))
p3 <- p3 + coord_equal()
p3 <- p3 + scale_fill_gradientn(TeX("Im($\\lambda_3$)"),limits=c(-0.01, 8),
                              colours=c("green", "yellow", "orange", "red"),
                              values = rescale(c(min(data3$Imag),
                                                 1,
                                                 4,
                                                 max(data3$Imag))))
p3 <- p3 + ggtitle(TeX("$\\delta_1 = 0$ and   $\\delta_2 = 0.00039$")) +
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))


data4 = read.csv("eigenvalues_alphas_0.1xDeg0xDegStar.csv")
data4$alpha_1 = data4$alpha_1/60/60
data4$alpha_2 = data4$alpha_2/60/60
p4 <- ggplot() 
p4 <- p4 + geom_raster(data = data4 , aes(x = alpha_1, y = alpha_2, fill = Imag))
p4 <- p4 + coord_equal()
p4 <- p4 + scale_fill_gradientn(TeX("Im($\\lambda_3$)"),limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(min(data4$Imag),
                                                   1,
                                                   4,
                                                   max(data4$Imag))))
p4 <- p4 + ggtitle(TeX("$\\delta_1 = 0.000039$ and   $\\delta_2 = 0$")) +
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))

data5 = read.csv("eigenvalues_alphas_10xDeg10xDegStar.csv")
data5$alpha_1 = data5$alpha_1/60/60
data5$alpha_2 = data5$alpha_2/60/60
p5 <- ggplot() 
p5 <- p5 + geom_raster(data = data5 , aes(x = alpha_1, y = alpha_2, fill = Imag))
p5 <- p5 + coord_equal()
p5 <- p5 + scale_fill_gradientn(TeX("Im($\\lambda_3$)"),limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(min(data5$Imag),
                                                   1,
                                                   4,
                                                   max(data5$Imag))))
p5 <- p5 + ggtitle(TeX("$\\delta_1 = 0.0039$ and   $\\delta_2 = 0.0039$")) +
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))


data6 = read.csv("eigenvalues_alphas_0.1xDeg0.1xDegStar.csv")
data6$alpha_1 = data6$alpha_1/60/60
data6$alpha_2 = data6$alpha_2/60/60
p6 <- ggplot() 
p6 <- p6 + geom_raster(data = data6 , aes(x = alpha_1, y = alpha_2, fill = Imag))
p6 <- p6 + coord_equal()
p6 <- p6 + scale_fill_gradientn(TeX("Im($\\lambda_3$)"),limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(min(data6$Imag),
                                                   1,
                                                   4,
                                                   max(data6$Imag))))
p6 <- p6 + ggtitle(TeX("$\\delta_1 = 0.000039$ and   $\\delta_2 = 0.000039$")) +
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))


data7 = read.csv("eigenvalues_alphas_1xDeg0.1xDegStar.csv")
data7$alpha_1 = data7$alpha_1/60/60
data7$alpha_2 = data7$alpha_2/60/60
p7 <- ggplot() 
p7 <- p7 + geom_raster(data = data7 , aes(x = alpha_1, y = alpha_2, fill = Imag))
p7 <- p7 + coord_equal()
p7 <- p7 + scale_fill_gradientn(TeX("Im($\\lambda_3$)"),limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(min(data7$Imag),
                                                   1,
                                                   4,
                                                   max(data7$Imag))))
p7 <- p7 + ggtitle(TeX("$\\delta_1 = 0.00039$ and   $\\delta_2 = 0.000039$")) +
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))

data8 = read.csv("eigenvalues_alphas_0.1xDeg1xDegStar.csv")
data8$alpha_1 = data8$alpha_1/60/60
data8$alpha_2 = data8$alpha_2/60/60
p8 <- ggplot() 
p8 <- p8 + geom_raster(data = data8 , aes(x = alpha_1, y = alpha_2, fill = Imag))
p8 <- p8 + coord_equal()
p8 <- p8 + scale_fill_gradientn(TeX("Im($\\lambda_3$)"),limits=c(-0.01, 8),
                                colours=c("green", "yellow", "orange", "red"),
                                values = rescale(c(min(data8$Imag),
                                                   1,
                                                   4,
                                                   max(data8$Imag))))
p8 <- p8 + ggtitle(TeX("$\\delta_1 = 0.000039$ and   $\\delta_2 = 0.00039$")) +
  xlab(TeX("$\\alpha_1$")) + ylab(TeX("$\\alpha_2$"))


F = ggarrange(p, p1, p2, p3, p4, p5, p6, p7, p8,
          ncol = 3, nrow = 3)

annotate_figure(F, 
                TeX("Im($\\lambda_3$) for antithetic model + damping"),
                fig.lab.pos = "top",
                fig.lab.face = "bold",
                fig.lab.size = 18)


