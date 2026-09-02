## Every test works on the example objects shipped with the package, whose
## file paths point at "inst/files" and have to be rewritten to the installed
## location. Loading them once here keeps that boilerplate out of the tests.
##
## The objects are loaded into the testthat environment under the same names
## the tests already use; call `freshExample()` when a test mutates one and a
## later test needs the original back.

freshExample <- function( name ){
    e <- new.env(parent = emptyenv())
    utils::data(list = name, package = "chraw", envir = e)
    rewrite_paths(get(name, envir = e))
}

ce_examples <- freshExample("ce_examples")
ce_chipseq  <- freshExample("ce_chipseq")
ce_atac     <- freshExample("ce_atac")
ce_rat      <- freshExample("ce_rat")

## The bam indexes are not shipped with the package (they are build
## artefacts), but several functions read the bam files directly.
indexBam( ce_examples )
indexBam( ce_atac )
indexBam( ce_rat )
