geom_arch_flip <- function(data, ..., n = 25, max.height = 10, bottom = TRUE){

  args <- list(...)
  args.aes <- biovizBase::parseArgsForAes(args)
  args.non <- biovizBase::parseArgsForNonAes(args)
  if("y" %in% names(args.aes))
    y.name <- quo_name(args.aes$y)
  else
    y.name <- NULL

  ## check required argument
  if(!all(c("x", "xend") %in% names(args.aes)))
    stop("x, xend, are requried in aes(), need to be passed into geom_arch()")
  startY <- rlang::eval_tidy(args.aes$y, data)
  endY <- rlang::eval_tidy(args.aes$yend, data)

  if("height" %in% names(args.aes)){
    if(!is.numeric(args.aes$height)){
      h <- rlang::eval_tidy(args.aes$height, data)
    }else{
      if(length(args.aes$height) == 1)
        h <- rep(args.aes$height, length(startY))
      else
        stop("unequal length of heights specified")
    }}else{
      h <- rep(max.height, length(startY))
    }
  if("x" %in% names(args.aes))
    x <- rlang::eval_tidy(args.aes$x, data)
  else
    x <- rep(0, length(startY))
  args.aes2 <- args.aes[!(names(args.aes) %in% c("x", "y", "group",
                                                 "hjust", "xend", "yend"))]
  xx<-c()
  yy<-c()
  for(i in 1:n){
    ang<-i*pi/(2*n)
    xx[i]<-sin(ang)
    yy[i]<-cos(ang)
  }
  ##takes the quarter of the curve calculated, flips a copy over the y axis
  ##reduces time spent in for loop
  if(bottom){
    yy <- c(1,yy,rev(-yy),-1, 1)
    xx <- c(0,xx,rev(xx), 0, 0)
  }else{
    yy <- c(1,yy,rev(-yy),-1)
    xx <- c(0,xx,rev(xx), 0)
  }
  ##SETS UP DATAFRAME TO KEEP TRACK OF ALL POINTS TO DRAW ALL ARCHES
  junc <- rep(seq_along(startY), each = length(yy))
  startY <- rep(startY, each = length(yy))
  endY <- rep(endY, each = length(yy))
  h <- rep(h, each = length(yy))
  x <- rep(x, each = length(yy))
  jump <- abs(endY - startY)
  jumpAdj <- if (length(jump)) max(jump) / max(abs(h)) else NA
  apoint <- data.frame(yy = yy * (abs(startY - endY) / 2) + (startY + endY) / 2,
                       xx = xx * h + x, junc,
                       s = ((abs(h) - jump / jumpAdj)) /
                         if (length(jump)) max(jump) else NA)
  data$junc <- seq_len(nrow(data))
  apoint <- merge(apoint, data, by = "junc")
  args.aes <- list(x = as.name("xx"),
                   y = as.name("yy"),
                   group = as.name("junc"))

  aesres <- do.call(aes, c(args.aes, args.aes2))
  if(nrow(apoint)){
    reslst <- c(list(data = apoint), list(aesres),args.non)
    p <- do.ggcall(geom_polygon, reslst)
    if("ylab" %in% names(args.non)){
      ylab <- quo_name(args.non$ylab)
    }else if(length(y.name)){
      ylab <- y.name
    }else{
      ylab <- ""
    }
    p <- list(p, ggplot2::ylab(ylab))
  }else{
    p <- NULL
  }
  p
}

geom_arch_flip2 <- function(data, ..., n = 25, max.height = 10, bottom = FALSE){


  args <- list(...)
  args.aes <- biovizBase::parseArgsForAes(args)
  args.non <- biovizBase::parseArgsForNonAes(args)
  if("y" %in% names(args.aes))
    y.name <- quo_name(args.aes$y)
  else
    y.name <- NULL

  ## check required argument
  if(!all(c("x", "xend") %in% names(args.aes)))
    stop("x, xend, are requried in aes(), need to be passed into geom_arch()")
  startY <- rlang::eval_tidy(args.aes$y, data)
  endY <- rlang::eval_tidy(args.aes$yend, data)

  if("height" %in% names(args.aes)){
    if(!is.numeric(args.aes$height)){
      h <- rlang::eval_tidy(args.aes$height, data)
    }else{
      if(length(args.aes$height) == 1)
        h <- rep(args.aes$height, length(startY))
      else
        stop("unequal length of heights specified")
    }}else{
      h <- rep(max.height, length(startY))
    }
  if("x" %in% names(args.aes))
    x <- rlang::eval_tidy(args.aes$x, data)
  else
    x <- rep(0, length(startY))
  args.aes2 <- args.aes[!(names(args.aes) %in% c("x", "y", "group",
                                                 "hjust", "xend", "yend"))]
  xx<-c()
  yy<-c()
  for(i in 1:n){
    ang<-i*pi/(2*n)
    xx[i]<-sin(ang)
    yy[i]<-cos(ang)
  }
  ##takes the quarter of the curve calculated, flips a copy over the y axis
  ##reduces time spent in for loop
  if(bottom){
    yy <- c(1,yy,rev(-yy),-1, 1)
    xx <- c(0,xx,rev(xx), 0, 0)
  }else{
    yy <- c(1,yy,rev(-yy),-1)
    xx <- c(0,xx,rev(xx), 0)
  }
  ##SETS UP DATAFRAME TO KEEP TRACK OF ALL POINTS TO DRAW ALL ARCHES
  junc <- rep(seq_along(startY), each = length(yy))
  startY <- rep(startY, each = length(yy))
  endY <- rep(endY, each = length(yy))
  h <- rep(h, each = length(yy))
  x <- rep(x, each = length(yy))
  jump <- abs(endY - startY)
  jumpAdj <- if (length(jump)) max(jump) / max(abs(h)) else NA
  apoint <- data.frame(yy = yy * (abs(startY - endY) / 2) + (startY + endY) / 2,
                       xx = xx * h + x, junc,
                       s = ((abs(h) - jump / jumpAdj)) /
                         if (length(jump)) max(jump) else NA)
  data$junc <- seq_len(nrow(data))
  apoint <- merge(apoint, data, by = "junc")
  args.aes <- list(x = as.name("xx"),
                   y = as.name("yy"),
                   group = as.name("junc"))

  aesres <- do.call(aes, c(args.aes, args.aes2))
  if(nrow(apoint)){
    reslst <- c(list(data = apoint), list(aesres),args.non)
    p <- do.ggcall(geom_path, reslst)
    if("ylab" %in% names(args.non)){
      ylab <- quo_name(args.non$ylab)
    }else if(length(y.name)){
      ylab <- y.name
    }else{
      ylab <- ""
    }
    p <- list(p, ggplot2::ylab(ylab))
  }else{
    p <- NULL
  }
  p
}


do.ggcall <- function(fun, args) {
  do.call(fun, filterArgs(fun, args))
}

filterArgs <- function(fun, args,
                       layerArgs=args[names(args) %in% c("geom", "stat")])
{
  resolveGeneric <- function(fun, args) {
    if (is(fun, "genericFunction")) {
      method <- selectMethod(fun, class(args$data))
      if (method@defined == "ANY") {
        ggfun <- get0(fun@generic, getNamespace("ggplot2"),
                      mode="function")
        if (!is.null(ggfun)) { # a generic overriding a ggplot2 function
          fun <- ggfun
        }
      }
    }
    fun
  }
  fun <- resolveGeneric(fun, args)
  ggplot2 <- !is(fun, "genericFunction")
  if (ggplot2) {
    aes <- vapply(args, is, "uneval", FUN.VALUE=logical(1L))
    args[aes] <- lapply(args[aes], filterArgs, fun=fun, layerArgs=layerArgs)
    if (is.null(names(args))) {
      args <- args[aes]
    } else {
      args <- ggplot2:::rename_aes(args)
      layer <- do.call(fun, layerArgs)
      validArgs <- c(names(formals(fun)),
                     layer$geom$aesthetics(),
                     layer$stat$aesthetics(),
                     layer$geom$parameters(TRUE),
                     layer$stat$parameters(TRUE))
      args <- args[names(args) %in% validArgs | aes]
    }
  }
  args
}

combineAes2 <- function(keep, lose){

  keep.nms <- names(keep)
  lose.nms <- names(lose)
  if("ymin" %in% keep.nms && "y" %in% lose.nms){
    lose$y <- keep$ymin
  }
  if("ymax" %in% keep.nms && "yend" %in% lose.nms){
    lose$yend <- keep$ymax
  }
  nms <- intersect(lose.nms, keep.nms)

  if(length(nms)){
    return(c(keep, lose[setdiff(lose.nms, keep.nms)]))
  }else{
    return(c(keep, lose))
  }
}


combineAes <- function(keep, lose){

  keep.nms <- names(keep)
  lose.nms <- names(lose)

  nms <- intersect(lose.nms, keep.nms)

  if(length(nms)){
    return(c(keep, lose[setdiff(lose.nms, keep.nms)]))
  }else{
    return(c(keep, lose))
  }
}
