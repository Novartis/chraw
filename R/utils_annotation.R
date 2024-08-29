#' Selects the corresponding BSgenome for a ChrawExperiment
#' @description This function inputs a ChrawExperiment and returns
#' the respective BSgenome object that corresponds to the experiment.
#'
#' @param x A ChrawExperiment object.
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
    bsGenome <- switch( x@referenceGenome,
                       "hg38" = "BSgenome.Hsapiens.UCSC.hg38",
                       "hg19" = "BSgenome.Hsapiens.UCSC.hg19",
                       "mm10" = "BSgenome.Mmusculus.UCSC.mm10",
                       "rn6" = "BSgenome.Rnorvegicus.UCSC.rn6",
                       "canFam3" = "BSgenome.Cfamiliaris.UCSC.canFam3"
                       )
    require( bsGenome, character.only=TRUE )
    bsGenome <- get( bsGenome )
    bsGenome
}

#' Selects the corresponding org.**.**.db package for a ChrawExperiment object
#' @description This function inputs a ChrawExperiment and returns
#' the corresponding org.**.**.db object for the experiment.
#'
#' @param x A ChrawExperiment object.
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
    orgDb <- switch( x@referenceGenome,
                    "hg38" = "org.Hs.eg.db",
                    "hg19" = "org.Hs.eg.db",
                    "mm10" = "org.Mm.eg.db",
                    "rn6" = "org.Rn.eg.db",
                    "canFam3" = "org.Cf.eg.db")
    if( is.null( orgDb ) ){
        return(NULL)
    }
    require( orgDb, character.only=TRUE )
    orgDb <- get( orgDb )
    orgDb
}


#' Selects the corresponding TxDb package for a ChrawExperiment object
#' @description This function inputs a ChrawExperiment and returns
#' the corresponding TxDb object for the experiment.
#'
#' @param x A ChrawExperiment object.
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
    txDb <- switch( x@referenceGenome,
                   "hg38" = "TxDb.Hsapiens.UCSC.hg38.refGene",
                   "hg19" = "TxDb.Hsapiens.UCSC.hg19.refGene",
                   "mm10" = "TxDb.Mmusculus.UCSC.mm10.knownGene",
                   "rn6" = "TxDb.Rnorvegicus.UCSC.rn6.refGene",
                   "canFam3" = "TxDb.Cfamiliaris.UCSC.canFam3.refGene")
    require( txDb, character.only=TRUE )
    txDb <- get( txDb )
    txDb
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

#' Selects the corresponding TxDb package for a ChrawExperiment object
#' @description This function inputs a ChrawExperiment and returns
#' the corresponding TxDb object for the experiment.
#'
#' @param x A ChrawExperiment object.
#'
#' @examples
#' data(ce_examples)
#' ce_examples <- rewrite_paths(ce_examples)
#'
#' chromhmm <- selectChromHMMData( ce_examples )
#'
#' @export
selectChromHMMData <- function( x ){
    if( !is(x, "ChrawExperiment") )
        stop("Input needs to be a 'ChrawExperiment' object")
    referenceGenome <- x@referenceGenome
    message("As no ChromHMM file was specified, one will be selected from the ENCODE project")
    if( !referenceGenome %in% c("hg38", "mm10") ){
        warning(sprintf("No default ChromHMM data was found for the species %s.", referenceGenome))
        return(NULL)
    }
    SP <- switch(
        x@referenceGenome,
        "hg38" = "https://www.encodeproject.org/files/ENCFF409CGA/@@download/ENCFF409CGA.bed.gz",
        "mm10" = "https://www.encodeproject.org/files/ENCFF391RNO/@@download/ENCFF391RNO.bed.gz")
    tmpFile <- tempfile()
    download.file(SP, destfile=tmpFile)
    tmpFile
}

simplifyChromHMMLabs <- function( x ){
    data( chromHMMDict )
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
    SP <- switch(
        object@referenceGenome,
        "hg38" = "human",
        "hg19" = "human",
        "mm10" = "mouse",
        "rn6" = "rat",
        "canFam3" = "dog" )
    regionRanges <- rowRanges(experiments(object)[[experimentName]])
    lnks <- paste0(
        "http://ucsc.na.novartis.net/cgi-bin/hgTracks?org=", SP,
        "&db=", object@referenceGenome,
        "&position=",
        as.character(seqnames(regionRanges)), ":",
        start(regionRanges), "-",
        end(regionRanges) )
    regionRanges$UCSC_link <- lnks
    elementMetadata(mcols(regionRanges))$type[colnames(mcols(regionRanges)) %in% "UCSC_link"] <- "UCSC link"
    rowRanges(experiments(object)[[experimentName]]) <- regionRanges
    object
}

