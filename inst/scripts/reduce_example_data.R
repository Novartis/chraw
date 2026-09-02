## Reduces the example files shipped in inst/files to the genomic windows that
## the vignette and the unit tests actually use. The bam and bigwig files were
## already subset when they were created, but the narrowPeak files still held
## genome-wide peak calls, which alone accounted for ~41 MB of the package.
##
## Run from the root of the package:
##   Rscript inst/scripts/reduce_example_data.R
## and then re-create the data objects with inst/scripts/write_test_objects.R.

suppressPackageStartupMessages({
    library(GenomicRanges)
    library(Rsamtools)
})

## Windows exercised downstream, with a margin on either side.
margin <- 1e6
keepRegions <- list(
    atacseq = GRanges("chr8", IRanges(115000152 - margin, 115999730 + margin)),
    chipseq = GRanges("chr6", IRanges(67090000 - margin, 67170000 + margin)),
    ## The rat data is only used to check that the workflow runs on a third
    ## species; a single window of chr1 is enough.
    rat     = GRanges("chr1", IRanges(1, 5e6))
)

narrowPeakCols <- c("seqname", "start", "end", "name", "score", "strand",
                    "signalValue", "pValue", "qValue", "peak")

reduceOne <- function( file, region ){
    before <- file.info(file)$size
    peaks <- utils::read.table(file, sep = "\t")
    colnames(peaks) <- narrowPeakCols[seq_len(ncol(peaks))]
    gr <- GRanges(peaks$seqname, IRanges(peaks$start + 1L, peaks$end))
    keep <- overlapsAny(gr, region)
    con <- gzfile(file, "w")
    on.exit(close(con), add = TRUE)
    utils::write.table(peaks[keep, , drop = FALSE], con, sep = "\t",
                       quote = FALSE, row.names = FALSE, col.names = FALSE)
    close(con)
    on.exit()
    after <- file.info(file)$size
    cat(sprintf("%-58s %7.2f MB -> %6.1f kB (%d of %d peaks)\n",
                file, before/1024^2, after/1024, sum(keep), length(keep)))
}

## The bam files were subset to a whole chromosome when they were made; keeping
## only the window above removes most of the remaining reads.
reduceBam <- function( file, region ){
    before <- file.info(file)$size
    idxFile <- paste0(file, ".bai")
    createdIndex <- !file.exists(idxFile)
    if( createdIndex ) indexBam(file)
    tmp <- tempfile(fileext = ".bam")
    filterBam(file, tmp,
              param = ScanBamParam(which = region, what = scanBamWhat()),
              indexDestination = FALSE)
    file.copy(tmp, file, overwrite = TRUE)
    unlink(tmp)
    ## Index files are build artefacts; the package rebuilds them with
    ## indexBam(), so they are not shipped.
    unlink(idxFile)
    after <- file.info(file)$size
    cat(sprintf("%-58s %7.2f MB -> %6.1f kB\n",
                file, before/1024^2, after/1024))
}

## The rat dataset only exists to show that the workflow runs on a third
## species, so it can be thinned further. A fixed seed keeps this reproducible.
downsampleBam <- function( file, fraction, seed = 20240101 ){
    before <- file.info(file)$size
    idxFile <- paste0(file, ".bai")
    if( !file.exists(idxFile) ) indexBam(file)
    qnames <- unique(scanBam(file, param = ScanBamParam(what = "qname"))[[1]]$qname)
    set.seed(seed)
    keep <- sample(qnames, ceiling(length(qnames) * fraction))
    tmp <- tempfile(fileext = ".bam")
    filterBam(file, tmp,
              param = ScanBamParam(what = scanBamWhat()),
              filter = FilterRules(list(sub = function(x) x$qname %in% keep)),
              indexDestination = FALSE)
    file.copy(tmp, file, overwrite = TRUE)
    unlink(c(tmp, idxFile))
    cat(sprintf("%-58s %7.2f MB -> %6.1f kB (%d of %d read names)\n",
                file, before/1024^2, file.info(file)$size/1024,
                length(keep), length(qnames)))
}

for (dataset in names(keepRegions)) {
    files <- list.files(file.path("inst/files", dataset, "peaks"),
                        pattern = "[.]narrowPeak[.]gz$", full.names = TRUE)
    for (f in files) reduceOne(f, keepRegions[[dataset]])
    bams <- list.files(file.path("inst/files", dataset, "bam"),
                       pattern = "[.]bam$", full.names = TRUE)
    for (f in bams) reduceBam(f, keepRegions[[dataset]])
    if (dataset == "rat") for (f in bams) downsampleBam(f, 0.35)
}

cat(sprintf("\ninst/files is now %.1f MB\n",
            sum(file.info(list.files("inst/files", recursive = TRUE,
                                     full.names = TRUE))$size) / 1024^2))
