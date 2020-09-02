## --- Multivariate TV-GARCH Model Class --- ####
##


## --- mtvgarch_class Definition --- ####

CORRtype = list(CCC=1,CEC=2,STCC1=3,STEC1=4)
CORRshape = list(single=1,double=2,double1loc=3)
CORRspeedopt = list(gamma=1,gamma_std=2,eta=3)

mtvgarch <- setClass(Class = "mtvgarch_class",
               slots = c(tvgarch_list="list",corrType="integer"),
               contains = c("namedList")
               )

## Initialise with no params
setMethod("initialize","mtvgarch_class",
          function(.Object,...){
            .Object <- callNextMethod(.Object,...)
            # Return:
            .Object
          })

setGeneric(name="mtvgarch",
           valueClass = "mtvgarch_class",
           signature = c("tvgarch_list","corrObj"),
           def = function(tvgarch_list,corrObj){
             this <- new("mtvgarch_class")
             this@tvgarch_list <- tvgarch_list
             this@corrType <- corrObj$type
             this$corr <- corrObj
             # Do validation checks:


             # End validation

             if (this@corrType == CORRtype$CCC){

             }
             if (this@corrType == CORRtype$STCC1){

             }



             return(this)
           }
)


setGeneric(name="estimateMTVGARCH",
           valueClass = "mtvgarch_class",
           signature = c("z","mtvgarchObj","estimationControl"),
           def=function(z,mtvgarchObj,estimationControl){
             this <- mtvgarchObj

             if (this@corrType == CORRtype$CCC){

             }
             if (this@corrType == CORRtype$STCC1){
               this$corr <- estimateSTCC1(z,this$corr,estimationControl)

             }

             return(this)

           })

