#' @rdname getGenomicBins
#' @param ... Parameters passed to object-specific methods.
#' @export
setGeneric("getGenomicBins", function( object, ... ) standardGeneric("getGenomicBins"))

#' @rdname plotDiff
#' @param ... Parameters passed to object-specific methods.
#' @export
setGeneric("plotVolcano", function( object, ... ) standardGeneric("plotVolcano"))

## `referenceGenome` is already a generic in BSgenome, which chraw imports.
## Defining a second one here would mask it (or be masked by it) depending on
## the attach order, so the method below is registered on BSgenome's generic.

#' @rdname referenceGenome
#' @param ... Parameters passed to object-specific methods.
#' @export
setGeneric("pipeline", function( object, ... ) standardGeneric("pipeline"))
