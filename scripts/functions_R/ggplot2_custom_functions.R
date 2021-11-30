# Done with the precious help of this webpage : https://eliocamp-github-io.translate.goog/codigo-r/2018/06/tu-propio-geom-smooth/?_x_tr_sl=auto&_x_tr_tl=en&_x_tr_hl=fr&_x_tr_pto=nui

#mean for "_smooth" ggplot2 functions
AngleSmooth <- function(formula, data, weights, n = 0.5) {
  ff <- data.frame()
  yi <- unique(round(data$x, 2))
  for (k in yi) {
    ww <- subset(data, round(x, 2) == k)
    f <- mean.angle.pi(data=ww$y)
    dd <- data.frame(matrix(ncol = 2, nrow = nrow(ww)))
    dd[1:nrow(ww),1] <- f
    dd[,2] <- ww$x
    setnames(dd, 1:2, c("y","x"))
    ff <- rbind(ff, dd)
    
  }
  model <- list(x = ff$x, pred = ff$y)
  class(model) <- "angle_smooth"
  return(model)
}

predictdf.angle_smooth <- function(model, xseq, se=FALSE, level) {
  data.frame(x = model$x, y = model$pred)
}

#mean angle for "summary" ggplot2 functions
Angle_Mean <- function(data) {
  data.frame(y = mean.angle.pi(data))
}

#standard deviation for error bars in stat_summary ggplot2 function 
st_dev_angle <- function(data){
  sd <- mean(sqrt( modulo.all((data - mean.angle.pi(data)))^2 ))
  data.frame( ymin=(mean.angle.pi(data))-sd,
             ymax=(mean.angle.pi(data))+sd
             )
}

#standard deviation value for stat_summary ggplot2 function 
st_dev_abs <- function(data){
  sd <- mean(sqrt( modulo.all((data - mean.angle.pi(data)))^2 ))
  data.frame( sda = sd
  )
}

