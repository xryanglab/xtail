dispersionMatrix.DESeqDataSet <- function(object){
  if (!"dispersionMatrix" %in% names(assays(object))) return (NULL)
  disp <- assays(object)[["dispersionMatrix"]]
  colnames(disp) <- colnames(object)
  disp
}

#' @export
setMethod("dispersionMatrix", signature(object="DESeqDataSet"),
          dispersionMatrix.DESeqDataSet)

#' @name dispersions
#' @rdname dispersions
#' @exportMethod "dispersions<-"
setReplaceMethod("dispersionMatrix", signature(object="DESeqDataSet", value="matrix"),
                 function(object, value) {
                  assays(object)[["dispersionMatrix"]] <- value
                  validObject( object )
                  object
                 })
