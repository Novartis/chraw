
test_that("Test constructors and validators of the ChrawExperiment object", {

    sampleFile <- system.file("files", "sampleSheet.txt", package="chraw", mustWork=TRUE)
    sampleDF <- as.data.frame(read.csv(sampleFile, sep = "\t", header=TRUE))

    inst_dir <- system.file("files/", package="chraw")
    sampleDF$bamFile <- gsub("inst/files", inst_dir, sampleDF$bamFile)
    sampleDF$bigwigFile <- gsub("inst/files", inst_dir, sampleDF$bigwigFile)
    sampleDF$peakFile <- gsub("inst/files", inst_dir, sampleDF$peakFile)
    sampleDF$qcFile <- gsub("inst/files", inst_dir, sampleDF$qcFile)

    ce <- ChrawExperimentFromDataFrame( sampleDF , genome = 'rn6')

    expect_s4_class( ce, "ChrawExperiment")

    colData(ce)$peakFile <- NULL
    colData(ce)$bamFile <- NULL
    expect_error(
        validObject(ce),
        "are missing" )

    expect_error(
      ChrawExperimentFromDataFrame(sampleDF),
      '"genome" is missing'
    )

    expect_error(ChrawExperimentFromDataFrame( sampleDF = '', genome = 'rn6' ),
                 'At least')

    expect_error(
      ChrawExperimentFromDataFrame( sampleDF[,"condition"], genome = 'rn6'),
      'At least one of the following required columns' )

    # expect_error(
    #     ChrawExperimentFromSampleFile( sampleFile ),
    #     '"genome" is missing' )
    #
    # expect_error(ChrawExperimentFromSampleFile( sampleFile = '', genome = 'rn6' ),
    #              'The sample sheet was not found.')
    #
    #
    # ce3 <- ChrawExperimentFromSampleFile( sampleFile, genome = 'mm10' )
    # expect_s4_class(ce3, "ChrawExperiment")

})
