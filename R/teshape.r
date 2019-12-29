#' Interconvert data structures
#'
#' Interconvert an array, a raw data frame, and frequency distribution
#' data.frame.
#'
#' Multivariate categorical data can be represented in several ways. Three comon
#' ways are : a contingency table, a data frame of raw observations (1 row = 1
#' subject), and a long data frame with a variable containing the counts in the
#' contingency table.
#'
#' @param data a data frame or array.
#' @param out the output format, see examples.
#' @param var the name of the frequency variable in the dataset, if not
#'   \code{"freq"}.
#' @param ... ...
#' @return a matrix containing the Markov basis as its columns (for easy
#'   addition to tables)
#' @seealso [stats::xtabs()]
#' @export
#' @examples
#'
#' # converting a talbe to a data frame
#' (tab <- structure(
#'   array(1:8, c(2,2,2)),
#'   .Dimnames = list(
#'     A = c("a1", "a2"),
#'     B = c("b1", "b2"),
#'     C = c("c1", "c2")
#'   )
#' ))
#'
#' teshape(tab, "freq")
#' teshape(tab, "raw") # nrow = sum(1:8)
#'
#'
#' # converting a summarized data frame into a table or raw data frame
#' (data <- teshape(tab, "freq"))
#' teshape(data, "tab")
#' stats::xtabs(freq ~ ., data = data)
#' teshape(data, "tab") == stats::xtabs(freq ~ ., data = data)
#' teshape(data, "raw")
#'
#'
#' # converting a raw data frame into a table or summarized data frame
#' (data <- teshape(tab, "raw"))
#' teshape(data, "tab")
#' teshape(data, "freq")
#'
#' 
teshape <- function(data, out = c("freq", "tab", "raw"), var){
  
  if(is.array(data)){
    input <- "tab"
  } else if(any(tolower(names(data)) %in% c("freq", "count", "frequency")) || !missing(var)){
    input   <- "freq"
    if(!missing(var) && !is.character(var)) var <- deparse(substitute(var))
  } else if(is.data.frame(data)){
    input <- "raw"
  } else {
    stop("input must be either \"tab\", \"raw\", or \"freq\"",
      call. = FALSE)	  	
  }

  
  out <- match.arg(out)
  
  
  
  if (input == out) return(data)
  
  
  
  
  ## if the input is a table...
  if(input == "tab"){
    df <- do.call(expand.grid, dimnames(data))
    p <- ncol(df)
    df$freq <- as.vector(data)
    if(out == "freq") return(df)
    
    # continuing if out == "raw"
    raw_df <- NULL
    for(k in 1:nrow(df)){
      if(df[k,]$freq > 0){
        raw_df <- rbind(raw_df, suppressMessages(
          df[rep(k, df[k,]$freq),1:p]
        ))
      }
    }
    raw_df <- raw_df[sample(nrow(raw_df)),]
    row.names(raw_df) <- 1:nrow(raw_df)
    return(raw_df)
  }
  
  
  
  ## if the input is a data frame of observations...  
  if(input == "raw"){
    if(out == "tab"){
      return(table(data))
    } else { # out == "freq"
      df <- melt(table(data))
      p <- ncol(df) - 1
      names(df) <- c(names(df)[-(p+1)], "freq")
      return(df)
    }
  }  
  
  
  
  ## if the input is a data frame of freqs...  
  if(input == "freq"){
    
    # make sure a frequency variable is found
  	if(!any(c("freq","count","frequency") %in% tolower(names(data))) && missing(var)){
  	  stop("frequency variable not found,  see ?teshape")
  	}
    
  	# move freq to end
  	if("freq" %in% tolower(names(data))){ 
  	  c <- which("freq" == tolower(names(data)))
  	} else if("count" %in% tolower(names(data))){ 
  	  c <- which("count" == tolower(names(data)))
  	} else if("frequency" %in% tolower(names(data))){ 
  	  c <- which("frequency" == tolower(names(data)))
  	} else if(!missing(var) && var %in% names(data)){
      c <- which(var == names(data))
  	}  	    		
  	df <- cbind(data[,(1:ncol(data))[-c]], freq = data[,c])
    p  <- ncol(data) - 1  	
  	
    if(out == "tab"){
      tab <- table(df[,1:p])
      tab[1:nrow(df)] <- df$freq
      return(tab)
    } else { # out == "raw"
      raw_df <- NULL
      for(k in 1:nrow(df)){
        if(df[k,]$freq > 0){
          raw_df <- rbind(raw_df, suppressMessages(
            df[rep(k, df[k,]$freq),1:p]
          ))
        }
      }
      raw_df <- raw_df[sample(nrow(raw_df)),]
      row.names(raw_df) <- 1:nrow(raw_df)
      return(raw_df)      
    }
    
  }   
  
}

