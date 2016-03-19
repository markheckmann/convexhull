#//////////////////////////////////////////////////////////////////////////////
#
#                               UTILITY FUNCTIONS 
#
#//////////////////////////////////////////////////////////////////////////////



# #' Set alpha value to color
# #' 
# #' Given a color in any standard format, the color is 
# #' converted to hex and an alpha value is added.
# #' 
# #' @param col A color in a standard format 
# #'  (e.g. \code{"#000000"}, \code{"black"}, \code{1}).
# #' @param alpha Alpha value to add to color 
# #' (default \code{=1}, i.e. opaque).
# #' @return hex color value.
# #' @export
# #' @examples
# #' 
# #'  set_alpha_color_value("black", .5)
# #'  set_alpha_color_value(1:3, .5)
# #'  set_alpha_color_value("#000000", .5)
# #'  
# set_alpha_color_value <- function(col, alpha=1) 
# {
#   k <- col2rgb(col) / 255
#   rgb(red = k["red", ], 
#       green = k["green", ],
#       blue = k["blue", ],
#       alpha = alpha, maxColorValue=1)  
# }
# 
# 


# #### +------------- Misc  -------------- ####
# 
# 
# # Add a class to class attributes
# #
# add_class <- function(x, add="")
# {
#   if (add != "" & !inherits(x, add)) { 
#     class(x) <- c(add, class(x))
#   }
#   x
# }
# 
# 
# #' Argument verification using partial matching plus numerics
# #' 
# #' The function is almost the same as \code{base::match.arg} excepts 
# #' that also the numerical index of the \code{choices} vector can be 
# #' used for matching. Setting \code{numerics=FALSE} is
# #' the same as using plain \code{match.arg}.
# #' 
# #' @inheritParams base::match.arg
# #' @param numerics Are numeric indexes sllowed for matching?
# #' @keywords internal
# #' @examples
# #' 
# #' # use an index for matching
# #' match.arg2("r", c("rows", "columns"))
# #' match.arg2(1, c("rows", "columns"))
# #' 
# match.arg2 <- function(arg, choices, several.ok=FALSE, numerics=TRUE) 
# {
#   if (numerics & is.numeric(arg))          
#     arg <- choices[arg]
#   match.arg(arg, choices, several.ok=several.ok)
# }
# 
# 
# 
# #### +------------- Logical functions -------------- ####
# 
# # does the dataframe contain the variables x, y, and z?
# #
# has_xyz <- function(x)
# {
#   all(c("x", "y", "z") %in% colnames(x))
# }
# 
# 
# stop_if_not_has_xyz <- function(x)
# {
#   if (! has_xyz(x))
#     stop("'x' must contain the variables: ", paste(vars, ""), call. = FALSE)
# }
# 
# 
# is_list_of_matrices <- function(x)
# {
#   is.list(x) & all(sapply(x, is.matrix))
# }
# 
# 
# stop_if_not_is_list_of_matrices <- function(x)
# {
#   if (! is_list_of_matrices(x))
#     stop("'x' must be a list", call. = FALSE)
# }
# 
# 
# stop_if_not_matrix <- function(x)
# {
#   if (! is.matrix(x))
#     stop("'x' must be a matrix", call. = FALSE)
# }
# 
# 
# # Do all matrices in a list have the same size (n x p)
# # x:  list of matrices
# #
# all_matrices_same_size <- function(x)
# {
#   all_matrices_same_nrow(x) & all_matrices_same_ncol(x)
# }
# 
# 
# all_matrices_same_nrow <- function(x)
# {
#   stop_if_not_is_list_of_matrices(x)
#   nr <- sapply(x, nrow)
#   all_equal_vec(nr)       # all values identical?
# }
# 
# 
# all_matrices_same_ncol <- function(x)
# {
#   stop_if_not_is_list_of_matrices(x)
#   nc <- sapply(x, ncol)
#   all_equal_vec(nc)       # all values identical?
# }
# 
# 
# stop_if_not_all_matrices_same_size <- function(x)
# {
#   if (! all_matrices_same_size(x))
#     stop("All matrices must have the same size")  
# }
# 
# 
# stop_if_not_all_matrices_same_nrow <- function(x)
# {
#   if (! all_matrices_same_nrow(x))
#     stop("All matrices must have the same number of rows")  
# }
# 
# 
# stop_if_not_all_matrices_same_ncol <- function(x)
# {
#   if (! all_matrices_same_ncol(x))
#     stop("All matrices must have the same number of columns")  
# }
# 
# 
# 
# # does the dataframe contain all necessary columns
# # for a elementinfo object?
# #
# is_elementinfo <- function(x)
# {
#   required <- c("dataset_id", "element", "x", "y", "z")
#   is.data.frame(x) & all(required %in% names(x))
# }
# 
# 
# stop_if_not_is_elementinfo <- function(x)
# {
#   if (!is_elementinfo(x))
#     stop("'x' is no dataframe or lacks one of the following columns:\n", 
#          "\tdataset_id  element  x  y  z", call.=FALSE)
# }
# 
# 
# # does the dataframe contain all necessary columns
# # for a constructinfo object?
# #
# is_constructinfo <- function(x)
# {
#   required <- c("dataset_id", "pole_low", "pole_high", "x", "y", "z")
#   is.data.frame(x) & all(required %in% names(x))
# }
# 
# 
# stop_if_not_is_constructinfo <- function(x)
# {
#   if (!is_constructinfo(x))
#     stop("'x' is no dataframe or lacks one of the following columns:\n", 
#          "\tdataset_id  pole_low  pole_high  x  y  z", call.=FALSE)
# }
# 
# 
# # does the dataframe contain variables necessary to calculate themes, 
# # i.e. orientation and theme_id
# #
# constructinfo_has_themes <- function(x)
# {
#   all(c("orientation", "theme_id") %in% names(x))
# }
# 
# 
# 
# #### +------------- Vectors functions -------------- ####
# 
# # Note that this section may also contain matrix methods 
# # for some generic functions
# 
# 
# extend_vec_to_length <- function(v, len=1) 
# {
#   u = v / sum(v^2)^.5   # unit vector in axis direction
#   u * len               # stretch vector to multiple of u
# }
# 
# 
# #' Length of vectors (euclidean)
# #' 
# #' Generic function for length of vectors. Can be applied to numeric vectors and
# #' matrices (rowwise or columnwise).
# #'
# #' @param x A vector or a matrix.
# #' @param along If \code{x} is a matrix, calculate vector length for
# #'   \code{1=rows}, \code{2=columns}.
# #' @export
# #' @keywords internal
# #' @examples
# #' # vector
# #' vec_length(1:3)          
# #' 
# #' # matrix
# #' x <- matrix(1:4, 2)
# #' vec_length(x)            # length of row vectors
# #' vec_length(x, along=2)   # length of column vectors
# #'   
# vec_length <- function(x, along=1, ...)
# {
#   UseMethod("vec_length")
# }
# 
# 
# #' @export
# #' @keywords internal
# vec_length.numeric <- function(x)
# {
#   x <- unlist(x)
#   sum(x^2)^.5
# }
# 
# #' @export
# #' @keywords internal
# vec_length.matrix <- function(x, along=1)
# {
#   apply(x, along, vec_length.numeric)
# }
# 
# 
# # deprecated (use vec_length instead)
# row_vec_length <- function(x)
# {
#   apply(x, 1, vec_length)
# }
# 
# 
# 
# #' Maximal vector length (euclidean) in matrix or list of matrices
# #' 
# #' Get length of longest row or column vector in a matrix or a list
# #' of matrices. If applied to a vector the result is identical to the 
# #' length of the vector.
# #' 
# #' @param x A matriox or list of matrices.
# #' @param along Calculate vector lengths for \code{1=rows}, \code{2=columns} of
# #'   (each) matrix.
# #' @export
# #' @keywords internal
# #' @examples
# #' # matrix
# #' x <- matrix(1:4, 2)
# #' max_vec_length(x)             # along rows
# #' max_vec_length(x, along=2)    # along columns
# #' 
# #' # list of matrices
# #' l <- random_matrix_list(2,2,2)
# #' max_vec_length(l)
# #' max_vec_length(l, along=2)
# #' 
# max_vec_length <- function(x, along=1, ...)
# {
#   UseMethod("max_vec_length")
# }
# 
# #' @export
# #' @keywords internal
# max_vec_length.numeric <- function(x)
# {  
#   max(vec_length(x))     # get length of vector
# }
# 
# 
# #' @export
# #' @keywords internal
# max_vec_length.matrix <- function(x, along=1)
# {  
#   max(vec_length(x, along=along))     # get length of longest row vector (most outer point)
# }
# 
# #' @export
# #' @keywords internal
# max_vec_length.data.frame <- function(x, along=1)
# {  
#   m <- as.matrix(x)
#   max_vec_length(m, along=along)
# }
# 
# 
# #' @export
# #' @keywords internal
# max_vec_length.list <- function(x, along=1)
# {
#   if (! all(sapply(x, is.matrix)))
#     stop("'x' must be a list of matrices")
#   
#   max(sapply(x, vec_length, along=along))     # get length of longest row vector (most outer point)
# }
# 
# 
# 
# 
# 
# #' Scale length of vectors in matrix
# #' 
# #' Scale length of vector in matrix or list of matrices 
# #' so the longest row or column vector has a certain length.
# #' 
# #' @note Note that the maximal length across all matrices is used for scaling.
# #' To apply this approach to each matrix seperately use \code{lapply}.
# #' 
# #' @param x A list of matrices.
# #' @param len Length of longest vector after scaling.
# #' @inheritParams max_vec_length
# #' @export
# #' @keywords internal
# #' @examples
# #' # vector
# #' scale_vec_length(1:3)  
# #' scale_vec_length(1:3, len=2) 
# #' 
# #' # matrix
# #' x <- matrix(1:4, 2)
# #' scale_vec_length(x)
# #' scale_vec_length(x, along=2)
# #' 
# #' # list of matrices
# #' l <- random_matrix_list(2,3,4)
# #' scale_vec_length(l)           
# #' scale_vec_length(l, along=2)
# #'  
# scale_vec_length <- function(x, len=1, along=1, ...)
# {
#   UseMethod("scale_vec_length")
# }
# 
# 
# #' @export
# #' @keywords internal
# scale_vec_length.numeric <- function(x, len=1)
# {
#   mx <- max_vec_length(x)   # get length of row vector
#   u <- x / mx               # devide by vector length to get unit length         
#   u * len                   # multiply matrix by desired max length
# }
# 
# 
# #' @export
# #' @keywords internal
# scale_vec_length.matrix <- function(x, len=1, along=1)
# {
#   mx <- max_vec_length(x, along)    # get longest row vector in matrix
#   u <- x / mx                       # divide by longest row vector to get unit length         
#   u * len                           # multiply matrix by desired max length
# }
# 
# 
# #' @export
# #' @keywords internal
# scale_vec_length.data.frame <- function(x, len=1, along=1)
# {
#   x <- as.matrix(x)
#   scale_vec_length(x, len, along)
# }
# 
# 
# 
# #' @export
# #' @keywords internal
# scale_vec_length.list <- function(x, len=1, along=1)
# {
#   stop_if_not_is_list_of_matrices(x)
#   
#   mx <- max_vec_length(x, along)              # get longest row vector
#   u <- mapply("/", x, mx, SIMPLIFY = FALSE)   # devide by longest row vector to get unit length         
#   mapply("*", u, len, SIMPLIFY = FALSE)       # multiply each matrix by desired max length
# }
# 
# 
# 
# #' Norm length of vectors in matrix
# #' 
# #' Norm length of all vectors in matrix or list of matrices 
# #' so they have the same length.
# #' 
# #' @param x A list of matrices.
# #' @param len Length of vectors after scaling.
# #' @inheritParams max_vec_length
# #' @export
# #' @keywords internal
# #' @examples
# #' # vector
# #' norm_vec_length(1:3)          # has length 1 now
# #' norm_vec_length(1:3, len=2)  
# #' 
# #' # matrix
# #' x <- matrix(1:4, 2)
# #' norm_vec_length(x)
# #' 
# #' # list of matrices
# #' l <- random_matrix_list(2,3,4)
# #' norm_vec_length(l)           
# #' norm_vec_length(l, along=2)
# #'  
# norm_vec_length <- function(x, len=1, along=1, ...)
# {
#   UseMethod("norm_vec_length")
# }
# 
# 
# #' @export
# #' @keywords internal
# norm_vec_length.numeric <- function(x, len=1, ...)
# {
#   scale_vec_length(x, len=len)
# }
# 
# 
# #' @export
# #' @keywords internal
# norm_vec_length.matrix <- function(x, len=1, along=1)
# {
#   # plyr::aaply(x, .margins=along, scale_vec_length, len=len)  # does not work
#   m <- x
#   if (along == 1) {
#     for (i in seq_along_rows(m)) {
#       m[i, ] <- scale_vec_length(x[i, ], len = len)
#     }  
#   } else if (along == 2) {
#     for (i in seq_along_cols(m)) {
#       m[ , i] <- scale_vec_length(x[ , i], len = len)
#     }   
#   }
#   m
# }
# 
# 
# #' @export
# #' @keywords internal
# norm_vec_length.data.frame <- function(x, len=1, along=1)
# {
#   x <- as.matrix(x)
#   norm_vec_length(x, len, along)
# }
# 
# 
# #' @export
# #' @keywords internal
# norm_vec_length.list <- function(x, len=1, along=1)
# {
#   stop_if_not_is_list_of_matrices(x)
#   lapply(x, norm_vec_length, len=len, along=along)
# }
# 
# 
# 
# 
# # get maximum length
# extend_vec_to_length <- function(v, len=1) 
# {
#   u = v / sum(v^2)^.5   # unit vector in axis direction
#   u * len               # stretch vector to multiple of u
# }
# 
# 
# # unipolar
# norm_to_circle <- function(R, j=NULL, dim=1:2, len=1) 
# {
#   R.new <- R[, dim]
#   if (is.null(j))
#     jj <- 1L:nrow(R)
#   for (j in jj) {
#     v = R[j, dim]                           # get biplot axis vector
#     R.new[j, dim] = extend_vec_to_length(v, len = len)   # extend to length of axis plus overhead
#   }
#   is.na(R.new) <- is.na(R.new)
#   R.new
# }
# 
# 
# 
# 
# #### +------------- Matrix functions   --------------- ####
# 
# 
# angles_between_all_row_vectors <- function(x)
# {
#   res <- NA_real_
#   nn <- 1L:nrow(x)
#   k <- 0
#   for (i in nn) {
#     for (j in nn) {
#       if (i < j ) {
#         k <- k + 1 
#         res[k] <- vec_angle(x[i, ], x[j, ], rad = FALSE, mode = 1) 
#       }
#     }
#   }
#   res
# }
# # angles_between_all_row_vectors(diag(3))
# 
# 
# 
# #' @rdname name-matrix
# #' @export
# random_colnames <- function(x, w=1)
# {
#   p <- ncol(x)
#   colnames(x) <- words::sentences(w, s=p)
#   x
# }
# 
# 
# #' @rdname name-matrix
# #' @export
# random_rownames <- function(x, w=1, bipolar=TRUE, sep=":")
# {
#   n <- nrow(x)
#   if (bipolar) {
#     rownames(x) <- paste0(words::sentences(w, s=n), sep, 
#                           words::sentences(w, s=n))
#   } else {
#     rownames(x) <- words::sentences(w, s=n)  
#   }
#   x
# }
# 
# 
# #' Add random row and/or column names to matrix
# #' 
# #' @param x A matrix.
# #' @param w Number of words for names.
# #' @param rownames,colnames Logical. Add random row / columns names?
# #' @param bipolar Logical. Are rownames bipolar (default \code{TRUE}).
# #' @param sep Seperator symbol between poles in case bipolar is \code{TRUE}.
# #' 
# #' @rdname name-matrix
# #' @export
# #' @keywords internal
# #' 
# name_matrix <- function(x, w=1, rownames=TRUE, colnames=TRUE, bipolar=TRUE, sep=":")
# {
#   if (rownames)
#     x <- random_colnames(x, w)
#   if (colnames)
#     x <- random_rownames(x, w, bipolar, sep)
#   x
# }
# 
# 
# 
# #' @export
# #' @keywords internal
# get_rownames.matrix <- function(x, default="")
# {
#   rn <- rownames(x)
#   if ( is.null(rn) ) {
#     rn <- rep(default, nrow(x))
#   }
#   names(rn) <- NULL
#   rn
# }
# 
# 
# #' @export
# #' @keywords internal
# get_rownames.list <- function(x, default="")
# {
#   # are list elements matrices
#   if (!all(sapply(x, is.matrix)))
#     stop("All list elements must be matrices.")
#   
#   rns <- lapply(x, get_rownames.matrix, default=default)
#   as.vector(unlist(rns))
# }
# 
# 
# #' @export
# #' @keywords internal
# get_rownames.default <- function(x, ...)
# {
#   stop("get_rownames needs a matrix or a list of matrices as input.", call. = FALSE)
# }
# 
# 
# #' Get rownames from matrix or list of matrices.
# #' 
# #' The rownames of a matrix or a list of matrices are retrieved. If the matrix
# #' has no rownames, a vector containing the default value for each row of the
# #' matrix is returned.
# #' 
# #' @param x A matrix or a list of matrices.
# #' @param default The default value if no row name exists.
# #' @return A vector of rownames or default values.
# #' @export
# #' @keywords internal
# #' @examples
# #' 
# #' # matrix
# #' x <- as.matrix(mtcars)
# #' get_rownames(x)
# #' 
# #' # list of matrices
# #' l <- list(x,x)
# #' get_rownames(l)
# #' 
# get_rownames <- function(x, default="", ...)
# {
#   UseMethod("get_rownames")
# }
# 
# 
# 
# 
# 
# # consecutive indexes for list of matrices
# #
# consecutive_index <- function(l, fun="nrow", mode="entry")
# {
#   fun <- match.fun(fun)    # nrow or ncol
#   
#   if (! is.list(l))
#     stop("'l' must be a list.", call. = FALSE)
#   s <- sapply(l, function(x) is.matrix(x))
#   if (!all(s))
#     stop("All list entries must be matrices.")
#   
#   # get consecutive row index
#   if (mode == "entry") {
#     r <- rep(1:length(l), sapply(l, fun))
#   } else  {
#     r <- lapply(l, function(x) 1L:fun(x))
#     r <- as.vector(unlist(r))                    # unlist and remove names
#   }
#   r  
# }
# 
# 
# 
# #' consecutive index for matrix list entries
# #' 
# #' Produces two types of consecutive indexes.
# #' 
# #' @param l A list.
# #' @param mode Type of index. \code{1="entry"} (default): consecutive index for
# #'   each list entry repeated as often as there are matrix/dataframe rows in
# #'   this entry. \code{2="rows"}: consecutive index from 1 to number of rows for
# #'   each matrix/dataframe in list, \code{3="columns"}: consecutive index from 1 to number of columns for
# #'   each matrix/dataframe in list. 
# #' @return Numeric vector.
# #' @export
# #' @keywords internal
# #' @rdname consecutive-index
# #' @examples
# #' 
# #' library(procrustes)
# #' x <- random_matrix_list(3, 4, 2)
# #' 
# #' # row indexes
# #' consecutive_rowindex(x, mode="entry")
# #' consecutive_rowindex(x, mode="rows")
# #' 
# #' # column indexes
# #' consecutive_colindex(x, mode="entry")
# #' consecutive_colindex(x, mode="columns")
# #' 
# consecutive_rowindex <- function(l, mode="entry")
# {
#   mode <- match.arg2(mode, c("entry", "rows"))
#   consecutive_index(l, fun="nrow", mode=mode)
# }
# 
# 
# #' @export
# #' @keywords internal
# #' @rdname consecutive-index
# consecutive_colindex <- function(l, mode="entry")
# {
#   mode <- match.arg2(mode, c("entry", "columns"))
#   consecutive_index(l, fun="ncol", mode=mode)
# }
# 
# 
# #### +------------- Dataframe functions  --------------- ####
# 
# 
# # shorten character columns. Can be used fopr custom printing functions
# # for dataframes.
# shorten_catcols <- function(x, width=10)
# {
#   for (j in 1L:ncol(x)) {
#     v <- x[[j]]
#     if ( is.character(v) ) {
#       x[[j]] <- stringr::str_sub(v, start=1, end=width)
#     }
#   }
#   x
# }
# 
# 
# #' Generate sequence along rows or columns of matrix or dataframe
# #' 
# #' @param x A matrix or a dataframe.
# #' @return Nuneric vector. Sequence along rows or columns. Returns
# #'   \code{integer()} for empty dataframes.
# #' @rdname seq-along
# #' @export
# #' @keywords internal
# #' @examples
# #' # standard dataframe
# #' seq_along_rows(mtcars)
# #' seq_along_cols(mtcars)
# #' 
# #' # empty dataframe
# #' d <- data.frame()
# #' seq_along_rows(d)
# #' seq_along_cols(d)
# #' 
# seq_along_rows <- function(x)
# {
#   if (!is.data.frame(x) & ! is.matrix(x))
#     stop("'x' must be a dataframe or a matrix.")
#   
#   n <- nrow(x)
#   if (n == 0) {
#     integer(0)
#   } else {
#     1L:n
#   }
# }
# 
# 
# #' @export
# #' @rdname seq-along
# #' 
# seq_along_cols <- function(x)
# {
#   if (!is.data.frame(x) & ! is.matrix(x))
#     stop("'x' must be a dataframe or a matrix.")
#   
#   p <- ncol(x)
#   if (p == 0) {
#     integer(0)
#   } else {
#     1L:p
#   }
# }
# 
# 
# #' Add NA columns for missing variables
# #' 
# #' @param x A dataframe.
# #' @param expected Vector of expected column names.
# #' @keywords internal
# #'  
# missing_vars_na <- function(x, expected)
# {
#   vv <- setdiff(expected, names(x))
#   for (v in vv) {
#     x[[v]] <- NA
#   }
#   x
# }
# 
# 
# #' Replace values in matching rows
# #' 
# #' Replaces the values for a set of selected variables in rows where all values 
# #' match for a set of criterion variables. If several rows match, only the data
# #' in the first row is used for replacement.
# #' 
# #' @param x Dataframe where values are replaced.
# #' @param y Dataframe with replacement values.
# #' @param by Variables to match rows by
# #' @param replace Variables to replace values for. Can be an existing on
# #'   non-existing variable (in \code{x}).
# #' @return A dataframe.
# #' @export
# #' @keywords internal
# #' @example examples/example-match-replace.R
# #' 
# match_replace <- function(x, y, by=NULL, replace=NULL)
# { 
#   nx <- names(x)
#   ny <- names(y)
#   
#   # default: use all common vars except the ones passed in replace
#   if (is.null(by)) {
#     nxy <- intersect(nx, ny)      # all common ones
#     by <- setdiff(nxy, replace)   # except the ones given in replace
#   }
#     
#   # default: vars not used for matching(by) by but in y
#   if (is.null(replace))
#     replace <- setdiff(ny, by)
# 
#   # match rows in x on y for criterion vars
#   ii <- row_match(x[ , by, drop=FALSE], 
#                   y[ , by, drop=FALSE], nomatch = NA)   
#   notna <- !is.na(ii)           # remove non-matching entries
#   ss <- seq_along(ii)           # row sequence
#   ii <- ii[notna]               # keep row sequence and indexes that match
#   ss <- ss[notna]
#   
#   x <- missing_vars_na(x, replace)   # add missing 
#   x[ss, replace] <- y[ii, replace]   # replace values 
#   x
# }
# 
# 
# xy <- function(x)
# {
#   x[ , c("x", "y"), drop=FALSE]
# }
# 
# 
# xyz <- function(x)
# {
#   stop_if_not_has_xyz(x)
#   x[ , c("x", "y", "z"), drop=FALSE]
# }
# 
# 
# "xyz<-" <- function(x, value)
# {
#   stop_if_not_has_xyz(x)
#   x[ , c("x", "y", "z")] <- value
#   x
# }
# 
# 
# 
# #### +------------- List functions  --------------- ####
# 
# 
# #' Are all list entries equal?
# #' 
# #' Comparison is performed by \code{all.equal}.
# #' @param x A list.
# #' @return Logical.
# #' @export
# #' @keywords internal
# #' @examples
# #' x <- replicate(3, diag(2), simplify = FALSE)
# #' all_list_entries_equal(x)
# #' x[[1]][1,1] <- 2
# #' all_list_entries_equal(x)
# #' 
# all_list_entries_equal <- function(x)
# {
#   if (!is.list(x))
#     stop("'x' must be a list.")
#   d <- mapply(all.equal, x[-1], x[1])
#   isTRUE(all(as.logical(d)))
# }
# 
# 
# #' All matrices have identical rownames?
# #' 
# #' Are all rownames identical? If the matrices do not have rownames (i.e.
# #' \code{NULL}) they are also considered to be identical.
# #' 
# #' @param x A list of matrices.
# #' @keywords internal
# #' 
# all_rownames_equal <- function(x)
# {
#   stop_if_not_is_list_of_matrices(x)
#   l <- lapply(x, rownames)
#   all_list_entries_equal(l)
# }
# 
# 
# list_to_dataframe <- function(x)
# {
#   m <- do.call(cbind, x)            # dataframe to join
#   as.data.frame(m)   
# }
# 
# 
# #### +------------- Positioning   --------------- ####
# 
# 
# # position labels by quadrant
# pos_by_quadrant <- function(x, y=NULL, ori=c(0,0)) 
# {
#   xy = xy.coords(rbind(x), y)
#   x = xy$x
#   y = xy$y
#   pos = rep(NA, length(x))
#   pos[x <= ori[1]]  <- 2   # left quadrants
#   pos[x > ori[1]] <- 4     # right quadrants
#   pos
# }
# 
# 
# 
# #### +------------- Plotting --------------- ####
# 
# 
# ptext <- function(label, ...) 
# {
#   usr <- par()$usr
#   text(usr[ 1 ], usr[ 4 ], label,  adj = c( 0, 1 ), ...)
# }
# 
# 
# one_segment <- function(from=c(0,0), to=NULL, stretch=1.1, ...) 
# {
#   if (is.null(to))
#     to <- from * stretch
#   if (any(to != 0) & all(!is.na(to))) {
#     segments(from[1], from[2], to[1], to[2], ...)
#     segments(-from[1], -from[2], -to[1], -to[2], ...)
#   }
# }
# 
# 
# all_segments <- function(x)
# {
#   for (i in 1L:nrow(x))
#     one_segment(x[i, ])
# }
# 
# 
# one_segment2 <- function(from=c(0,0), to=NULL, stretch=1.1, ...) 
# {
#   if (is.null(to))
#     to <- from * stretch
#   if (any(to != 0) & all(!is.na(to))) {
#     segments(from[1], from[2], to[1], to[2], ...)
#   }
# }
# 
# 
# draw_labels <- function(x, pole1="", pole2="", stretch=1.1) 
# {
#   x <- x * stretch
#   if (any(x != 0) & all(!is.na(x))) {
#     pos <- pos_by_quadrant(x)
#     text(x[1], x[2], labels=pole1, pos=pos)
#     x <- -x
#     pos <- pos_by_quadrant(x)
#     text(x[1], x[2], labels=pole2, pos=pos)
#   }
# }
# 
# 
# draw_label2 <- function(x, label="", stretch=1.1, ...) 
# {
#   x <- x * stretch
#   if (any(x != 0) & all(!is.na(x))) {
#     pos <- pos_by_quadrant(x)
#     text(x[1], x[2], labels=label, pos=pos, ...)
#   }
# }
# 
# 
# draw_label_background <- function(x, label="", stretch=1.1, col="#FFFFFF30", ...) 
# {
#   x <- x * stretch
#   if (any(x != 0) & all(!is.na(x))) {
#     pos <- pos_by_quadrant(x)
#     textbox2(x[1], x[2], label=label, pos=pos, col=col, ...)
#   }
# }
# 
# 
# 
