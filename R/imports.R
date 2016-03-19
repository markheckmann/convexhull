#//////////////////////////////////////////////////////////////////////////////
#
#       code copied from other packages to remove import overhead
#
#//////////////////////////////////////////////////////////////////////////////


# code from likert::recode

#' recode function from likert package
#' 
#' Allows vector input of recodes.
#' 
#' @param x The vector whose values will be recoded.
#' @param from  The old values in x to be recoded.
#' @param to	The new values.
#' @param to.class	an 'as.' function representing the desired vector type (i.e.
#'   as.character, as.numeric, as.logical, as.numeric).
#' @export
#' @keywords internal
#' 
recode3 <- function (x, from, to, to.class = NULL) 
{
  if (is.null(to.class)) {
    if (is.character(to)) {
      to.class <- as.character
    }
    else {
      to.class <- as.integer
    }
  }
  if (length(from) != length(to)) {
    stop("The length of from and to do not match")
  }
  r <- rep(to.class(NA), length(x))
  if (is.factor(x) & is.numeric(from)) {
    from = levels(x)[from]
  }
  for (i in seq_along(from)) {
    r[x == from[i]] = to[i]
  }
  return(to.class(r))
}



#### recode from car ####

# # recode function (J. Fox)
# # last modified 2014-08-04 by J. Fox
# 
# recode <- function(var, recodes, as.factor.result, as.numeric.result=TRUE, levels){
#   lo <- -Inf
#   hi <- Inf
#   recodes <- gsub("\n|\t", " ", recodes)
#   recode.list <- rev(strsplit(recodes, ";")[[1]])
#   is.fac <- is.factor(var)
#   if (missing(as.factor.result)) as.factor.result <- is.fac
#   if (is.fac) var <- as.character(var)
#   result <- var
#   for (term in recode.list){
#     if (0 < length(grep(":", term))) {
#       range <- strsplit(strsplit(term, "=")[[1]][1],":")
#       low <- try(eval(parse(text=range[[1]][1])), silent=TRUE)
#       if (class(low) == "try-error"){
#         stop("\n  in recode term: ", term, 
#              "\n  message: ", low)
#       }
#       high <- try(eval(parse(text=range[[1]][2])), silent=TRUE)
#       if (class(high) == "try-error"){
#         stop("\n  in recode term: ", term, 
#              "\n  message: ", high)
#       }
#       target <- try(eval(parse(text=strsplit(term, "=")[[1]][2])), silent=TRUE)
#       if (class(target) == "try-error"){
#         stop("\n  in recode term: ", term, 
#              "\n  message: ", target)
#       }
#       result[(var >= low) & (var <= high)] <- target
#     }
#     else if (0 < length(grep("^else=", squeezeBlanks(term)))) {
#       target <- try(eval(parse(text=strsplit(term, "=")[[1]][2])), silent=TRUE)
#       if (class(target) == "try-error"){
#         stop("\n  in recode term: ", term, 
#              "\n  message: ", target)
#       }
#       result[1:length(var)] <- target
#     }
#     else {
#       set <- try(eval(parse(text=strsplit(term, "=")[[1]][1])), silent=TRUE)
#       if (class(set) == "try-error"){
#         stop("\n  in recode term: ", term, 
#              "\n  message: ", set)
#       }
#       target <- try(eval(parse(text=strsplit(term, "=")[[1]][2])), silent=TRUE)
#       if (class(target) == "try-error"){
#         stop("\n  in recode term: ", term, 
#              "\n  message: ", target)
#       }
#       for (val in set){
#         if (is.na(val)) result[is.na(var)] <- target
#         else result[var == val] <- target
#       }
#     }
#   }
#   if (as.factor.result) {
#     result <- if (!missing(levels)) factor(result, levels=levels) 
#     else as.factor(result)
#   }
#   else if (as.numeric.result && (!is.numeric(result))) {
#     result.valid <- na.omit(result)
#     opt <- options("warn"=-1)
#     result.valid <- as.numeric(result.valid)
#     options(opt)
#     if (!any(is.na(result.valid))) result <- as.numeric(result)
#   }
#   result
# }
# 



#### draw.circle from plotrix ####

# from the plotrix package

getYmult <- function () 
{
  if (dev.cur() == 1) {
    warning("No graphics device open.")
    ymult <- 1
  }
  else {
    xyasp <- par("pin")
    xycr <- diff(par("usr"))[c(1, 3)]
    ymult <- xyasp[1]/xyasp[2] * xycr[2]/xycr[1]
  }
  return(ymult)
}


draw_circle <- function (x, y, radius, nv = 100, border = NULL, col = NA, lty = 1, 
          lwd = 1) 
{
  xylim <- par("usr")
  plotdim <- par("pin")
  ymult <- getYmult()
  angle.inc <- 2 * pi/nv
  angles <- seq(0, 2 * pi - angle.inc, by = angle.inc)
  if (length(radius) < length(x)) 
    radius <- rep(radius, length.out = length(x))
  if (length(col) < length(radius)) 
    col <- rep(col, length.out = length(radius))
  for (circle in 1:length(radius)) {
    xv <- cos(angles) * radius[circle] + x[circle]
    yv <- sin(angles) * radius[circle] * ymult + y[circle]
    polygon(xv, yv, border = border, col = col[circle], lty = lty, 
            lwd = lwd)
  }
  invisible(list(x = xv, y = yv))
}



#### row.match from prodlim ####

row_match <- function (x, table, nomatch = NA) 
{
  if (class(table) == "matrix") 
    table <- as.data.frame(table)
  if (is.null(dim(x))) 
    x <- as.data.frame(matrix(x, nrow = 1))
  cx <- do.call("paste", c(x[, , drop = FALSE], sep = "\r"))
  ct <- do.call("paste", c(table[, , drop = FALSE], sep = "\r"))
  match(cx, ct, nomatch = nomatch)
}
