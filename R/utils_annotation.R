## Loads an annotation package that is listed in 'Suggests' and returns the
## object it exports under its own name. `pkg` is NULL whenever the reference
## genome of the object has no annotation package associated with it.
loadAnnotationPackage <- function( pkg, what, genome ){
    if( is.null( pkg ) ){
        stop(sprintf(
            "No %s is available for the reference genome '%s'.", what, genome ))
    }
    if( !requireNamespace( pkg, quietly=TRUE ) ){
        stop(sprintf(
            paste("The package '%s' is needed to obtain the %s for the",
                  "reference genome '%s', but it is not installed.",
                  "Install it with:\n\tBiocManager::install(\"%s\")"),
            pkg, what, genome, pkg ))
    }
    getExportedValue( pkg, pkg )
}

#' Selects the corresponding BSgenome for a ChrawExperiment
#' @description This function inputs a ChrawExperiment and returns
#' the respective BSgenome object that corresponds to the experiment.
#'
#' @param x A ChrawExperiment object.
#'
#' @return A [BSgenome::BSgenome] object with the sequence of the reference
#' genome of `x`.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' bsgenome <- selectBSgenome( ce_examples )
#'
#' @importClassesFrom BSgenome BSgenome
#' @export
selectBSgenome <- function( x ){
    if( !is(x, "ChrawExperiment") )
        stop("Input needs to be a 'ChrawExperiment' object")
    genome <- referenceGenome( x )
    bsGenome <- switch( genome,
                       "hg38" = "BSgenome.Hsapiens.UCSC.hg38",
                       "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
                       "mm10" = "BSgenome.Mmusculus.UCSC.mm10",
                       "mm9" = "BSgenome.Mmusculus.UCSC.mm9",
                       "rn6" = "BSgenome.Rnorvegicus.UCSC.rn6",
                       "canFam3" = "BSgenome.Cfamiliaris.UCSC.canFam3"
                       )
    loadAnnotationPackage( bsGenome, "BSgenome package", genome )
}

#' Selects the corresponding org.**.**.db package for a ChrawExperiment object
#' @description This function inputs a ChrawExperiment and returns
#' the corresponding org.**.**.db object for the experiment.
#'
#' @param x A ChrawExperiment object.
#'
#' @return An `OrgDb` object, as defined by the \pkg{AnnotationDbi} package,
#' with the gene annotation of the organism of `x`.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' org <- selectOrgDb( ce_examples )
#'
#' @importClassesFrom AnnotationDbi OrgDb
#' @export
selectOrgDb <- function( x ){
    if( !is(x, "ChrawExperiment") )
        stop("Input needs to be a 'ChrawExperiment' object")
    genome <- referenceGenome( x )
    orgDb <- switch( genome,
                    "hg38" = "org.Hs.eg.db",
                    "hg19" = "org.Hs.eg.db",
                    "mm10" = "org.Mm.eg.db",
                    "mm9" = "org.Mm.eg.db",
                    "rn6" = "org.Rn.eg.db",
                    "canFam3" = "org.Cf.eg.db")
    loadAnnotationPackage( orgDb, "organism annotation package", genome )
}


#' Selects the corresponding TxDb package for a ChrawExperiment object
#' @description This function inputs a ChrawExperiment and returns
#' the corresponding TxDb object for the experiment.
#'
#' @param x A ChrawExperiment object.
#'
#' @return A [GenomicFeatures::TxDb] object with the transcript annotation of
#' the reference genome of `x`.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' txdb <- selectTxDb( ce_examples )
#'
#' @importClassesFrom GenomicFeatures TxDb
#' @importFrom txdbmaker makeTxDbFromGFF
#' @export
selectTxDb <- function( x ){
    if( !is(x, "ChrawExperiment") )
        stop("Input needs to be a 'ChrawExperiment' object")
    genome <- referenceGenome( x )
    txDb <- switch( genome,
                   "hg38" = "TxDb.Hsapiens.UCSC.hg38.refGene",
                   "hg19" = "TxDb.Hsapiens.UCSC.hg19.refGene",
                   "mm10" = "TxDb.Mmusculus.UCSC.mm10.knownGene",
                   "mm9" = "TxDb.Mmusculus.UCSC.mm9.knownGene",
                   "rn6" = "TxDb.Rnorvegicus.UCSC.rn6.refGene",
                   "canFam3" = "TxDb.Cfamiliaris.UCSC.canFam3.refGene")
    loadAnnotationPackage( txDb, "TxDb package", genome )
}

#' @importClassesFrom IRanges CompressedCharacterList
#' @importMethodsFrom S4Vectors %in%
#' @importFrom ChIPseeker annotatePeak
#' @importFrom ChIPpeakAnno annotatePeakInBatch
subsetTxDbByGenes <- function( txDb, activeGeneList ){
    if( !all(grepl("^[0-9]*$", activeGeneList)) &
        metadata(txDb)[metadata(txDb)$name == "Organism","value"] != "Macaca fascicularis" ){
        stop("At least one identifier from the 'activeGeneList' is not an Entrez ID")
    }
    gffFile <- asGFF(txDb)
    gffFile$phase <- NA
    gffFile$phase[gffFile$type == "CDS"] <- 0
    geneParents <- gffFile[gffFile$Name %in% activeGeneList]
    txParents <- gffFile[any(gffFile$Parent %in% geneParents$ID)]
    featureSons <- gffFile[any(gffFile$Parent %in% txParents$ID)]
    featureSons$Parent <- featureSons$Parent[featureSons$Parent %in% txParents$ID]
    gffFile <- c(geneParents, txParents, featureSons)
    tmpGFF <- tempfile()
    export( gffFile, con=tmpGFF, format="gff3" )
    txDbSub <- makeTxDbFromGFF( tmpGFF, format="gff3" )
    txDbSub
}

#' Downloads the default ChromHMM annotation for a ChrawExperiment object
#' @description This function inputs a ChrawExperiment and downloads the
#' default ChromHMM segmentation that the ENCODE project distributes for the
#' reference genome of the experiment. Downloads are cached with
#' [BiocFileCache::BiocFileCache], so the file is retrieved only once.
#'
#' @param x A ChrawExperiment object.
#'
#' @return A character string with the path to the cached ChromHMM file, or
#' `NULL` if no default annotation is available for the reference genome of
#' `x`.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' # Requires internet access the first time it is run.
#' \donttest{
#' chromhmm <- selectChromHMMData( ce_examples )
#' }
#'
#' @importFrom BiocFileCache BiocFileCache bfcrpath
#' @export
selectChromHMMData <- function( x ){
    if( !is(x, "ChrawExperiment") )
        stop("Input needs to be a 'ChrawExperiment' object")
    genome <- referenceGenome( x )
    message("As no ChromHMM file was specified, one will be selected from the ENCODE project")
    if( !genome %in% c("hg38", "mm10") ){
        warning(sprintf("No default ChromHMM data was found for the species %s.", genome))
        return(NULL)
    }
    SP <- switch(
        genome,
        "hg38" = "https://www.encodeproject.org/files/ENCFF409CGA/@@download/ENCFF409CGA.bed.gz",
        "mm10" = "https://www.encodeproject.org/files/ENCFF391RNO/@@download/ENCFF391RNO.bed.gz")
    bfc <- BiocFileCache( ask=FALSE )
    unname( bfcrpath( bfc, SP ) )
}

simplifyChromHMMLabs <- function( x ){
    x2 <- x
    for( i in seq_len(nrow(chromHMMDict)) ){
        x2 <- gsub(
            chromHMMDict$chromHmmLabs[i],
            chromHMMDict$simpleLabs[i],
            x2 )
    }
    x2[!x2 %in% chromHMMDict$simpleLabs] <- NA
    x2
}

createUCSCLinks <- function( object, experimentName ){
    ## Create UCSC links ##
    genome <- referenceGenome( object )
    SP <- switch(
        genome,
        "hg38" = "human",
        "hg19" = "human",
        "mm10" = "mouse",
        "mm9" = "mouse",
        "rn6" = "rat",
        "canFam3" = "dog" )
    regionRanges <- rowRanges(experiments(object)[[experimentName]])
    lnks <- paste0(
        "https://genome.ucsc.edu/cgi-bin/hgTracks?org=", SP,
        "&db=", genome,
        "&position=",
        as.character(seqnames(regionRanges)), ":",
        start(regionRanges), "-",
        end(regionRanges) )
    regionRanges$UCSC_link <- lnks
    elementMetadata(mcols(regionRanges))$type[colnames(mcols(regionRanges)) %in% "UCSC_link"] <- "UCSC link"
    rowRanges(experiments(object)[[experimentName]]) <- regionRanges
    object
}

