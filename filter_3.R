#' Compute filter_3 values
filter_3 <- function(data, AI_max = NULL, AI_min = NULL) {
  
  if(is.null(AI_max)){
    stop("You should insert highest AI value: AI_max")
  }
  
  if(is.null(AI_min)){
    stop("You should insert lowest AI value: AI_min")
  }
  
  attach(data)
  filter3 <- ifelse((C-O-N)==0, 1,
                  ifelse(((1+C-O-0.5*HIon)/(C-O-N))<AI_max,
                         ifelse(((1+C-O-0.5*HIon)/(C-O-N))>AI_min,
                                ifelse(N/C<=NtoC_max,1
                                       ,0),0),0))
  detach(data)
  
  return(filter3)
}
