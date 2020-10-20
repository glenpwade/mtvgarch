## --- Multivariate TV-GARCH Model Class --- ####
##


## --- mtvgarch_class Definition --- ####
mtvgarch <- setClass(Class = "mtvgarch_class",
               slots = c(Tobs="integer",N="integer"),
               contains = c("namedList")
               )

## === Initialise  ===####
setMethod("initialize","mtvgarch_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            # Default initial values
            .Object@N <- as.integer(0)
            .Object@Tobs <- as.integer(0)
            # Return:
            .Object
          })

## === Constructor  ===####
setGeneric(name="mtvgarch",
           valueClass = "mtvgarch_class",
           signature = c("tvgarch_list","series.names"),
           def = function(tvgarch_list,series.names){
             this <- new("mtvgarch_class")

             # Do basic validation checks:
             if(length(tvgarch_list) < 2) {
               warning("multivariate mtvgarch objects require a minimum of 2 data series.")
               return(this)
             }
             if(length(tvgarch_list) != length(series.names) ) {
               warning("The number of series.names provided must match the number of data series in the list.")
               return(this)
             }
             # End validation

             this@N <- as.integer(length(tvgarch_list))
             this@Tobs <- tvgarch_list[[1]]$tvObj@Tobs

             for(n in 1:this@N){
               # Validate the class_type:
               objType <- class(tvgarch_list[[n]])
               if(objType[1] != "tvgarch_class") {
                 warning("multivariate mtvgarch objects require data to be a list of tvgarch_class objects")
                 return(this)
               }

               this[[series.names[n]]] <- tvgarch_list[[n]]
             }

             return(this)
           }
)

## ===   FilterData ===####
setGeneric(name="filterData",
           valueClass = "matrix",
           signature = c("e","mtvgarchObj"),
           def=function(e,mtvgarchObj){
             this <- mtvgarchObj

             # Do basic validation checks:
             objType <- class(e)
             if(objType[1] != "matrix"){
               warning("data must be a matrix")
               return(matrix())
             }
             objType <- class(mtvgarchObj)
             if(objType[1] != "mtvgarch_class"){
               warning("mtvgarchObj must be an instance of the mtvgarch_class")
               return(matrix())
             }
             if(this@N != NCOL(e)){
               warning("The number of data series in mtvgarch must match the number of series in 'e'")
               return(matrix())
             }
             # End validation

             filteredData <- matrix(nrow = this@Tobs, ncol = this@N)

             for(n in 1:this@N){

               filteredData[,n] <- e[,n]/sqrt(this[[n]]$Estimated$tv@g * this[[n]]$Estimated$garch@h)
             }
             colnames(filteredData) <- names(this)

             return(filteredData)

           }
)

