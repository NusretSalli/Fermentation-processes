# 3D phase plane analysis by the cube.R extension for grind.R 
# This script defines cube() and run3d() and uses the plot3D library.
# The viewpoint is defined by the defaults theta=40 and phi=40,
# where theta gives the azimuthal direction and phi the colatitude.
# One can use plotdev(theta=60,phi=40) to change the viewpoint afterwards.
# plotdev() "refreshes" the picture (which can be used to make a pdf).
# One can set par(mar=c(0,0,0,0)) to reduce the margins.
# Rob de Boer, Utrecht University

library(plot3D)

cube <- function(x=80, y=100, z=100, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1, log="", npixels=150, state=s, parms=p, odes=model, time=0, eps=0, show=NULL, zero=TRUE, lwd=1, ...) {
  # Make a 3D-phase space with nullclines
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(args_run,args_plot,names(formals(contour3D)),"ticktype","nticks")])
    if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
  }
  if (!is.null(dots)) dots_run <- dots[names(dots) %in% args_run]
  else dots_run <- NULL
  if (!is.numeric(x)) x <- index(x,names(state))
  if (!is.numeric(y)) y <- index(y,names(state))
  if (!is.numeric(z)) z <- index(z,names(state))
  if (!is.null(show)) ishows <- index(show, names(state))
  else ishows <- c(x, y, z)
  nvar <- length(state)
  if (zero) state[1:nvar] <- rep(0,nvar)
  logx <- ifelse(grepl('x',log), TRUE, FALSE)
  logy <- ifelse(grepl('y',log), TRUE, FALSE)
  logz <- ifelse(grepl('z',log), TRUE, FALSE)
  if (logx) xc <- 10^seq(log10(xmin),log10(xmax),length.out=npixels)
  else xc <- seq(xmin+eps,xmax,length.out=npixels)
  if (logy) yc <- 10^seq(log10(ymin),log10(ymax),length.out=npixels)
  else yc <- seq(ymin+eps,ymax,length.out=npixels)
  if (logz) zc <- 10^seq(log10(zmin),log10(zmax),length.out=npixels)
  else zc <- seq(zmin+eps,zmax,length.out=npixels)
  xvar <- names(state)[x]; yvar <- names(state)[y]; zvar <- names(state)[z]
  
  npixels2 <- npixels^2
  vparms <- as.list(parms)
  vparms <- lapply(vparms,rep,vparms,npixels2)
  vstate <- as.list(state)
  #vstate<-lapply(vstate,rep,vstate,npixels2);vstate[[x]]<-0;vstate[[y]]<-0
  for (j in seq(1,nvar)) if (j!=x & j!=y & j!=z) vstate[[j]]<-rep(vstate[[j]],npixels2);
  #xzy
  vstate[[x]] <- rep.int(xc, npixels)
  vstate[[z]] <- rep.int(zc, rep.int(npixels, npixels))  #outer(xc,zc)
  vstate[[y]] <- rep.int(yc[1], npixels2)
  dxzy1 <- odes(time,vstate,vparms)[[1]]
  dim(dxzy1) <- c(npixels,npixels,nvar)
  vstate[[y]] <- rep.int(yc[npixels], npixels2)
  dxzy2 <- odes(time,vstate,vparms)[[1]]
  dim(dxzy2) <- c(npixels,npixels,nvar)
  #xyz
  vstate[[y]] <- rep.int(yc, rep.int(npixels, npixels))
  vstate[[z]] <- rep.int(zc[1], npixels2)
  dxyz1 <- odes(time,vstate,vparms)[[1]]
  dim(dxyz1) <- c(npixels,npixels,nvar)
  vstate[[z]] <- rep.int(zc[npixels], npixels2)
  dxyz2 <- odes(time,vstate,vparms)[[1]]
  dim(dxyz2) <- c(npixels,npixels,nvar)
  #yzx
  vstate[[y]] <- rep.int(yc, npixels)
  vstate[[z]] <- rep.int(zc, rep.int(npixels, npixels))
  vstate[[x]] <- rep.int(xc[1], npixels2)
  dyzx1 <- odes(time,vstate,vparms)[[1]]
  dim(dyzx1) <- c(npixels,npixels,nvar)
  vstate[[x]] <- rep.int(xc[npixels], npixels2)
  dyzx2 <- odes(time,vstate,vparms)[[1]]
  dim(dyzx2) <- c(npixels,npixels,nvar)
  
  add <- FALSE
  for (i in ishows) {
    try(contour3D(x=xc,y=yc[1],z=zc,colvar=dxzy1[,,i],levels=0,col=colors[i],add=add,lwd=lwd,xlim=c(xmin,xmax),ylim=c(ymin,ymax),zlim=c(zmin,zmax),xlab=names(state[x]),ylab=names(state[y]),zlab=names(state[z]),...),silent=TRUE)
    add <- TRUE
    try(contour3D(x=xc,y=yc[npixels],z=zc,colvar=dxzy2[,,i],levels=0,col=colors[i],add=TRUE,lwd=lwd,...),silent=TRUE)
    
    try(contour3D(x=xc,y=yc,z=zc[1],colvar=dxyz1[,,i],levels=0,col=colors[i],add=TRUE,lwd=lwd,...),silent=TRUE)
    try(contour3D(x=xc,y=yc,z=zc[npixels],colvar=dxyz2[,,i],levels=0,col=colors[i],add=TRUE,lwd=lwd,...),silent=TRUE)
    
    try(contour3D(x=xc[1],y=yc,z=zc,colvar=dyzx1[,,i],levels=0,col=colors[i],add=TRUE,lwd=lwd,...),silent=TRUE)
    try(contour3D(x=xc[npixels],y=yc,z=zc,colvar=dyzx2[,,i],levels=0,col=colors[i],add=TRUE,lwd=lwd,...),silent=TRUE)
  }
}

run3d <- function(x=1, y=2, z=3, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1, log="", col=1, add=FALSE, state=s, ...) {
  dots <- list(...)
  
  if (!is.null(dots)) {
    
    args_run3d <- c(names(formals(scatter3D)),"ticktype","nticks")
    
    unknown <- names(dots[!names(dots) %in% c(args_run,args_run3d)])
    
    if (length(unknown)>0) warning(paste("Unknown argument(s):",unknown,sep=" "))
    
    dots_run <- dots[names(dots) %in% args_run]
    
    dots_l3d <- dots[names(dots) %in% args_run3d]
    
  } else dots_run <- NULL
  
  data <- do.call('run',c(list(state=state,timeplot=FALSE,table=TRUE),dots_run))
  
  if (!add) do.call('lines3D',c(list(x=data[,1+x],y=data[,1+y],z=data[,1+z],xlim=c(xmin,xmax),ylim=c(ymin,ymax),zlim=c(zmin,zmax),xlab=names(state[x]),ylab=names(state[y]),zlab=names(state[z]),col=col),dots_l3d))
  
  else do.call('lines3D',c(list(x=data[,1+x],y=data[,1+y],z=data[,1+z],col=col,add=TRUE),dots_l3d))
  
  f <- state
  
  f[1:length(f)] <- as.numeric(data[nrow(data),2:(length(s)+1)])
  
  return(f)
}



