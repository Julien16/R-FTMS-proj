#' Compute filter_1 values

filter_1 <- function(data, OtoC_min = NULL, OtoC_max = NULL, Intensity_Min = NULL, 
                     Intensity_Max = NULL, ExpMass_Min = NULL, 
                     ExpMass_Max = NULL, OC = NULL, 
                     Intensity = NULL) {
  
  if(is.null(OtoC_min)){
    stop("You should insert lowest O/C ratio threshold value: OtoC_min")
  }
  
  if(is.null(OtoC_max)){
    stop("You should insert highest O/C ratio threshold value: OtoC_max")
  }
  
  if(is.null(Intensity_Min)){
    stop("You should insert lowest intensity threshold value: Intensity_Min")
  }
  
  if(is.null(Intensity_Max)){
    stop("You should insert highest intensity threshold value: Intensity_Max")
  }
  
  if(is.null(ExpMass_Min)){
    stop("You should insert lowest mass threshold value: ExpMass_Min")
  }
  
  if(is.null(ExpMass_Max)){
    stop("You should insert highest mass threshold value: ExpMass_Max")
  }
  
  if(is.null(Intensity)){
    stop("You should insert Intensity value: Intensity")
  }
  
  if(is.null(OC)){
    stop("You should insert Oxygen/Carbon value: OC")
  }
  
  filter1 <- ifelse(OC<OtoC_max,
                    ifelse(OC>OtoC_min,
                           ifelse(formula_list$Intensity<Intensity_Max,
                                  ifelse(formula_list$Intensity>Intensity_Min,
                                         ifelse(formula_list$ExpMass<ExpMass_Max,
                                                ifelse(formula_list$ExpMass>ExpMass_Min,1
                                                       ,0),0),0),0),0),0)
  
  return (filter1)
  
}
