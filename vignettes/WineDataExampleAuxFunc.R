# function that helps in plotting the functional data
df_plot_fda <- function(S_data, time_df, s=1, a= 0.8, n_sample = 10){
  # to do
  # check that the length of the time data agree with the dim of the data
  # write an if statment that in case of empty time_df generate a one based on the dim S_data
  df_plot <- as.data.frame(t(S_data))
  # selecting samples
  df_plot_n <- df_plot[,1:n_sample]
  # add the time var
  df_plot_n$time <- time_df
  df_plot_melt_n <- melt(df_plot_n, id.vars=c("time"))
  # ploting
  plot_n <- ggplot(df_plot_melt_n, aes(x=time, y = value, color=variable)) + geom_line(size=s, alpha=a) + theme_minimal() + theme(legend.position = "none", axis.title = element_blank(), axis.text.y = element_blank())  
  res <- list(plot_n, df_plot_melt_n)
  return(res)
}
# prepare the data
data_prepare <- function(f_data, t_data){
  colnames(f_data) <- NULL
  colnames(t_data) <- NULL
  f_data <- t(f_data)
  ready_data <- cbind(t_data, f_data)
  ready_data <- as.matrix(ready_data)
  return(ready_data)
}
# Prepare knots
Knots_prepare <- function(selected_knots, Time){
  knots_normalized <- selected_knots / max(selected_knots)
  knots_normalized = knots_normalized*(max(Time)- min(Time)) + min(Time)
  # taking the first three decimal number
  knots_normalized = as.numeric(format(round(knots_normalized,4), nsmall = 1))
  return(knots_normalized)
}
#
GetProjCovEig <- function(f_ready_data, ready_knots){
  ProjObj <- project(f_ready_data,ready_knots)
  Sigma=cov(ProjObj$coeff)
  Spect=eigen(Sigma,symmetric = T)
  EigenSp=lincomb(ProjObj$basis,t(Spect$vec))
  Proj_C_EV <- list(ProjObj = ProjObj, Sigma=Sigma, Spect=Spect, EigenSp=EigenSp)
  return(Proj_C_EV)
}
#
compare_fit <- function(N_plots, SubsampleFractions, ProjCovEigObj){
  if(!is.vector(SubsampleFractions)){
    stop("SubsampleFractions has to be a vector of subsamples")
  }
  if(!(length(SubsampleFractions) == N_plots)){
    stop("length of SubsampleFractions has to be equal to the N_plots")
  }
  EgnFncts <- list()
  for(i in 1:N_plots){
    EgnFncts[i] = subsample(ProjCovEigObj$EigenSp,1:SubsampleFractions[i])
  }
  return(EgnFncts)
}
#
get_Compare_Plot <- function(ready_f_data, ProjCovEigObj, portions = c(), EgnFuncObj, sampleNumber){
  C_mat=ProjCovEigObj$ProjObj$coeff %*% ProjCovEigObj$Spect$vec
  return({matplot(ready_f_data[,1],ready_f_data[,sampleNumber],type='l',lty=1,xlab='',ylab='')
    lines(ProjCovEigObj$ProjObj$sp,sID=sampleNumber-1,col='red',lty=2,lwd=1)
    if(length(portions > 0)){
      lines(lincomb(EgnFuncObj[[1]],C_mat[1,1:portions[1],drop=F]),col='green')
      lines(lincomb(EgnFuncObj[[2]],C_mat[1,1:portions[2],drop=F]),col='brown')
      lines(lincomb(EgnFuncObj[[3]],C_mat[1,1:portions[3],drop=F]),col='blue',lwd=1)}
  })
}
#
get_amse_diff <- function(f_ready_data, ProjCovEigObj, t_data){
  if(is.data.frame(t_data)){
    f_hat <- evspline(object = ProjCovEigObj$ProjObj$sp, x = t_data[,1])
  }
  if(is.vector(t_data)){
    f_hat <- evspline(object = ProjCovEigObj$ProjObj$sp, x = t_data)
  }
  f_diff <- f_hat - f_ready_data
  AMSE_diff <- dim(f_ready_data)[1] * amse(f_diff[,-1])
  # AMSE_diff <- amse(f_diff[,-1])
  AMSE_diff_list = list(f_hat, f_diff, AMSE_diff)
  return(AMSE_diff_list)
}
#
plot_EigenfunctionEigenValueScaled <- function(ProjCovEigObj, EigenNumber, mrgn=2, type='l', bty="n",col='deepskyblue4',
                                               lty=1, lwd=2, xlim=NULL, ylim = NULL, xlab="", ylab = "", vknots=TRUE){
  if(!is.numeric(EigenNumber)){
    stop(" please insert the number of eigenfunctions as EigenMuber")
  }
  y = evspline(ProjCovEigObj$EigenSp, sID = 1:EigenNumber)
  Arg=y[,1]
  Val=y[,-1,drop=F]
  # if(is.null(xlim)){
  #   xlim = range(Arg)
  # }
  # if(is.null(ylim)){
  #   ylim = range(Val)
  # }
  plot(Arg,Val[,1]*sqrt(ProjCovEigObj$Spect$values[1]),type=type,bty=bty,col=col,xlim=xlim,ylim=ylim,
       xlab=xlab,ylab=ylab,lty=lty,lwd=lwd)
  ourcol=c( 'darkorange3', 'goldenrod', 'darkorchid4',
            'darkolivegreen4', 'deepskyblue', 'red4',
            'slateblue','deepskyblue4')
  for(i in 2:EigenNumber){
    lines(Arg,Val[,i]*sqrt(ProjCovEigObj$Spect$values[i]),col=ourcol[(i-2)%%8+1],lty=lty,lwd=lwd)
  }
  if(vknots){
    abline(v = ProjCovEigObj$EigenSp@knots, lty = 3, lwd = 0.5)
  }
  abline(h = 0, lwd = 0.5)
}