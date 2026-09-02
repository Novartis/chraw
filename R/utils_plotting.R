## Shared look for every ggplot produced by the package. Keeping it in one
## place means the sample labels on the x axis are rotated the same way
## everywhere, instead of each function (and each vignette) picking its own
## angle.
##
## @param baseSize Base font size passed to `cowplot::theme_cowplot()`.
## @param rotateXLabels Logical, whether the x axis text should be rotated.
##   Set to `FALSE` for plots whose x axis is a continuous quantity.
#' @import cowplot
#' @import ggplot2
themeChraw <- function( baseSize = 12, rotateXLabels = TRUE ){
    th <- theme_cowplot( font_size = baseSize )
    if( rotateXLabels ){
        th <- th + theme(
            axis.text.x = element_text( angle=35, hjust=1, vjust=1 ) )
    }
    th
}

## rtracklayer reads bigwig files through the UCSC library, which parses the
## drive letter of a Windows path as a URL protocol and then fails with
## "Unrecognized protocol C in udcProtNew". A path that carries no colon works,
## so on Windows the absolute path is rewritten relative to the working
## directory whenever that is possible (same drive). Everywhere else, and
## whenever no relative path exists, the original path is returned unchanged.
bigwigPath <- function( path ){
    if( .Platform$OS.type != "windows" ) return( path )
    vapply( path, function( p ){
        if( is.na(p) || !file.exists(p) ) return( p )
        target <- normalizePath( p, winslash="/", mustWork=FALSE )
        base <- normalizePath( getwd(), winslash="/", mustWork=FALSE )
        ts <- strsplit( target, "/", fixed=TRUE )[[1L]]
        bs <- strsplit( base, "/", fixed=TRUE )[[1L]]
        ## Different drives: no relative path can be built.
        if( !identical( tolower(ts[1L]), tolower(bs[1L]) ) ) return( p )
        i <- 1L
        while( i <= min(length(ts), length(bs)) &&
               identical( tolower(ts[i]), tolower(bs[i]) ) ) i <- i + 1L
        rel <- paste( c( rep("..", length(bs) - i + 1L),
                         ts[seq(i, length(ts))] ), collapse="/" )
        if( file.exists( rel ) ) rel else p
    }, character(1), USE.NAMES=FALSE )
}

## Calls `fun` with bigwig paths the UCSC parser can cope with. When the files
## cannot be reached by a relative path (a different drive, say), the working
## directory is moved to the directory holding them and only the file names are
## passed, which is the only other way to get rid of the colon. The original
## working directory is always restored.
withBigwigPaths <- function( paths, fun ){
    if( .Platform$OS.type != "windows" ) return( fun( paths ) )
    rel <- bigwigPath( paths )
    if( !any( grepl(":", rel, fixed=TRUE) ) ) return( fun( rel ) )
    dirs <- unique( dirname( paths ) )
    if( length( dirs ) == 1L && dir.exists( dirs ) ){
        oldDir <- setwd( dirs )
        on.exit( setwd( oldDir ), add=TRUE )
        return( fun( basename( paths ) ) )
    }
    fun( paths )
}

## Axis and legend titles that are shared by several plots, so that the same
## quantity is never given two different names.
chrawLabels <- list(
    sample      = "Sample (replicate)",
    fragmentsM  = expression("Fragments (" * 10^6 * ")"),
    readsM      = expression("Reads (" * 10^6 * ")"),
    log2fc      = expression(log[2] ~ "fold change"),
    fragLength  = "Fragment length (bp)"
)
