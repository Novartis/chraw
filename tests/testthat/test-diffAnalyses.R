test_that("Test for differential functions", {
    data("ce_examples")
    ce_examples <- rewrite_paths(ce_examples)

    # expect_error(
    #     pullDiffResults( ce_examples, experimentName="FragLengthDist" ),
    #     "No differential results were found" )

    expect_error( testForDiffSignal(
        object="error", experimentName ="Peaks",
        design = ~condition,
        contrasts = list(
            contrast1=c("condition", "agonist_3h", "CTRL_3h"),
            contrast2=c("condition", "agonist_27h", "CTRL_27h") ) ),
        "parameter must be a" )

    expect_error( testForDiffSignal(
        object=ce_examples, experimentName ="notexistent",
        design = ~condition,
        contrasts = list(
          contrast1=c("condition", "agonist_3h", "CTRL_3h"),
          contrast2=c("condition", "agonist_27h", "CTRL_27h") ) ),
        "The count experiment 'notexistent' was not found" )

    expect_error( testForDiffSignal(
        object=ce_examples, experimentName ="Peaks",
        design = ~nonsense,
        contrasts = list(
          contrast1=c("condition", "agonist_3h", "CTRL_3h"),
          contrast2=c("condition", "agonist_27h", "CTRL_27h") ) ),
        "All variables in the formula must be column" )

    expect_error( testForDiffSignal(
        object=ce_examples, experimentName ="Peaks",
        design = ~condition,
        contrasts = list(
            contrast1=c("nonsense", "agonist_3h", "CTRL_3h"),
            contrast2=c("condition", "agonist_27h", "CTRL_27h") ) ),
        "The following variables" )

    expect_error( testForDiffSignal(
        object=ce_examples, experimentName="Peaks",
        design = ~condition,
        contrasts="nonsense"),
        "The parameter 'contrasts'")

    expect_error( testForDiffSignal(
        object=ce_examples, experimentName="Peaks",
        design = ~condition,
        contrasts=list(c("nonsense", "A", "B")),
        "The list passed to the'") )


    lvs <- levels(factor(colData(ce_examples)[colnames(experiments(ce_examples)[["Peaks"]]),"condition"]))
    stopifnot(length(lvs) == 4)

    expect_error(
        testForDiffSignal(
            object=ce_examples, experimentName ="Peaks",
            design = ~condition,
            contrasts = list(
                contrast01=c("condition", "nonsense", lvs[2]),
                contrast02=c("condition", lvs[1], lvs[2])) ),
        "were not found as values" )


    expect_error(
        testForDiffSignal(
            object=ce_examples, experimentName ="Peaks",
            design = ~condition,
            contrasts = list(
                contrast01=c("condition", lvs[1], lvs[2], "nonsense"),
                contrast02=c("condition", lvs[2], lvs[1])) ),
        "each contrast must be" )

    expect_message( tst <- testForDiffSignal(
        object=ce_examples, experimentName ="Peaks",
        design = ~condition,
        contrasts = list(
            contrast01=c("condition", lvs[1], lvs[3]),
            contrast02=c("condition", lvs[2], lvs[4])),
        dropLevels=FALSE ), "Fitting all data together" )

    expect_error(
        pullDiffResults( tst, experimentName="nonsense" ),
        "The experiment 'nonsense' was not found" )

    expect_error(
        pullDiffResults( tst, experimentName="Peaks", contrastName=c("contrast01", "non", "sense") ),
        "was/were not found in the ChrawExperiment" )

    diffRes1 <- pullDiffResults( tst, experimentName="Peaks", contrastName="contrast01" )
    expect_s4_class( diffRes1, "GRanges" )
    diffRes2 <- pullDiffResults( tst, experimentName="Peaks" )
    diffRes3 <- pullDiffResults( tst, experimentName="Peaks", outputFormat="df_wide" )
    diffRes4 <- pullDiffResults( tst, experimentName="Peaks", outputFormat="df_long" )
    expect_identical(mcols(diffRes1)$contrast01_pvalue, diffRes2$contrast01_pvalue)
    expect_identical(mcols(diffRes1)$contrast01_pvalue, diffRes3$contrast01_pvalue)
    expect_identical(mcols(diffRes1)$contrast01_pvalue, diffRes4$pvalue[diffRes4$contrastName == "contrast01"])

    expect_s4_class( rowData(experiments(tst)[["Peaks"]]), "DataFrame" )
    expect_true(all(colnames(elementMetadata( rowData(experiments(tst)[["Peaks"]]) )) %in% c("type", "description", "contrastName", "design")))

    expect_s3_class( getDiffResultSummary( tst ), "data.frame" )
    #expect_error( getDiffResultSummary( ce_examples ), "No results were found" )

    tst <- testForDiffSignal(
      object=ce_examples, experimentName ="Peaks",
      design = ~condition,
      contrasts = list(
        contrast01=c("condition", lvs[2], lvs[1]),
        contrast02=c("condition", lvs[3], lvs[1]),
        contrast03=c("condition", lvs[4], lvs[1])),
      dropLevels=FALSE )

    tst2 <- testForDiffSignal(
        object=ce_examples, experimentName ="Peaks",
        design = ~condition, dropLevels=FALSE )

    expect_identical(tst, tst2)

    ceSub <- ce_examples[,colData(ce_examples)$condition %in% c("agonist_3h", "agonist_27h")]

    #rs0 <- rowData(experiments(tst3)[["Peaks"]])
    rs1 <- rowData(experiments(tst2)[["Peaks"]])#[,elementMetadata(rowData(experiments(tst2)[["Peaks"]]))$contrastName == "contrast01"]

    expect_warning( tst <- testForDiffSignal(
        object=tst, experimentName ="Peaks",
        design = ~condition,
        contrasts = list(
          contrast01=c("condition", lvs[2], lvs[1]),
          contrast02=c("condition", lvs[3], lvs[1]),
          contrast03=c("condition", lvs[4], lvs[1])),
        dropLevels=FALSE ),
        "The ChrawExperiment object has already" )

    expect_identical(
        rowData(experiments(tst)[["Peaks"]])$contrast01_pvalue,
        rs1$contrast01_pvalue )

    suppressWarnings(tst <- testForDiffSignal(
        object=tst, experimentName ="Peaks",
        design = ~condition,
        contrasts = list( contrast01=c("condition", lvs[2], lvs[1]) ),
        replaceAll=TRUE )
    )

    colData(tst)$sizeFactors <- 1

    expect_message( tst2 <- testForDiffSignal(
        object=tst, experimentName ="Peaks",
        design = ~condition,
        contrasts = list( contrast01b=c("condition", lvs[2], lvs[1]) ),
        replaceAll=TRUE, useColSizeFactors=TRUE ),
        "Using size factors defined in the ChrawExperiment" )

    rs3 <- rowData(experiments(tst)[["Peaks"]])
    expect_true(all(elementMetadata(rs3)$contrastName[elementMetadata(rs3)$type == "differential_test"] == "contrast01"))

    xx <- plotVolcano( tst, experimentName="Peaks" )
    expect_s3_class(xx, "gg")

    xx <- plotMA( tst, experimentName="Peaks" )
    expect_s3_class( xx, "gg" )
} )
