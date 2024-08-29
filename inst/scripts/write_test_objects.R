#CREATE CE_CHIPSEQ
sampleName <- c("H3K27ac_NR1I3_CTRL_3h_rep1","H3K27ac_NR1I3_CTRL_3h_rep2",
                "H3K27ac_NR1I3_CTRL_27h_rep1","H3K27ac_NR1I3_CTRL_27h_rep2",
                "H3K27ac_NR1I3_agonist_3h_rep1","H3K27ac_NR1I3_agonist_3h_rep2",
                "H3K27ac_NR1I3_agonist_27h_rep1","H3K27ac_NR1I3_agonist_27h_rep2")

condition <- c("H3K27ac_NR1I3_CTRL_3h","H3K27ac_NR1I3_CTRL_3h",
          "H3K27ac_NR1I3_CTRL_27h","H3K27ac_NR1I3_CTRL_27h",
          "H3K27ac_NR1I3_agonist_3h","H3K27ac_NR1I3_agonist_3h",
          "H3K27ac_NR1I3_agonist_27h","H3K27ac_NR1I3_agonist_27h")

bamFile <- c("inst/files/chipseq/bam/ce_chipseq_ctrl_3h_rep1.bam",
             "inst/files/chipseq/bam/ce_chipseq_ctrl_3h_rep2.bam",
             "inst/files/chipseq/bam/ce_chipseq_ctrl_27h_rep1.bam",
             "inst/files/chipseq/bam/ce_chipseq_ctrl_27h_rep2.bam",
             "inst/files/chipseq/bam/ce_chipseq_agonist_3h_rep1.bam",
             "inst/files/chipseq/bam/ce_chipseq_agonist_3h_rep2.bam",
             "inst/files/chipseq/bam/ce_chipseq_agonist_27h_rep1.bam",
             "inst/files/chipseq/bam/ce_chipseq_agonist_27h_rep2.bam")

peakFile <- c("inst/files/chipseq/peaks/ce_chipseq_ctrl_3h_rep1.narrowPeak.gz",
              "inst/files/chipseq/peaks/ce_chipseq_ctrl_3h_rep2.narrowPeak.gz",
              "inst/files/chipseq/peaks/ce_chipseq_ctrl_27h_rep1.narrowPeak.gz",
              "inst/files/chipseq/peaks/ce_chipseq_ctrl_27h_rep2.narrowPeak.gz",
              "inst/files/chipseq/peaks/ce_chipseq_agonist_3h_rep1.narrowPeak.gz",
              "inst/files/chipseq/peaks/ce_chipseq_agonist_3h_rep2.narrowPeak.gz",
              "inst/files/chipseq/peaks/ce_chipseq_agonist_27h_rep1.narrowPeak.gz",
              "inst/files/chipseq/peaks/ce_chipseq_agonist_27h_rep2.narrowPeak.gz")

bigwigFile <- c("inst/files/chipseq/bigwig/ce_chipseq_ctrl_3h_rep1.bigwig",
                "inst/files/chipseq/bigwig/ce_chipseq_ctrl_3h_rep2.bigwig",
                "inst/files/chipseq/bigwig/ce_chipseq_ctrl_27h_rep1.bigwig",
                "inst/files/chipseq/bigwig/ce_chipseq_ctrl_27h_rep2.bigwig",
                "inst/files/chipseq/bigwig/ce_chipseq_agonist_3h_rep1.bigwig",
                "inst/files/chipseq/bigwig/ce_chipseq_agonist_3h_rep2.bigwig",
                "inst/files/chipseq/bigwig/ce_chipseq_agonist_27h_rep1.bigwig",
                "inst/files/chipseq/bigwig/ce_chipseq_agonist_27h_rep2.bigwig")

qcFile <- c("inst/files/chipseq/qc/ctrl_3h/ctrl_3h.qc.json",
            "inst/files/chipseq/qc/ctrl_3h/ctrl_3h.qc.json",
            "inst/files/chipseq/qc/ctrl_27h/ctrl_27h.qc.json",
            "inst/files/chipseq/qc/ctrl_27h/ctrl_27h.qc.json",
            "inst/files/chipseq/qc/agonist_3h/agonist_3h.qc.json",
            "inst/files/chipseq/qc/agonist_3h/agonist_3h.qc.json",
            "inst/files/chipseq/qc/agonist_27h/agonist_27h.qc.json",
            "inst/files/chipseq/qc/agonist_27h/agonist_27h.qc.json")

sampleDF_ <- data.frame(sampleName,bamFile,bigwigFile,peakFile,qcFile,condition,row.names = sampleName)
ce_chipseq <- ChrawExperimentFromDataFrame(sampleDF_,genome="mm10",minimalMetadata = FALSE)
save(ce_chipseq, file="data/ce_chipseq.RData",compress=TRUE)
############################################################################################################################################################################################
#CREATE CE_ATACSEQ
sampleName <- c("TRPS1_CTRL_rep1",
                "TRPS1_CTRL_rep2",
                "TRPS1_KO_rep1",
                "TRPS1_KO_rep2")

condition <- c("TRPS1_CTRL",
          "TRPS1_CTRL",
          "TRPS1_KO",
          "TRPS1_KO")

bamFile <- c("inst/files/atacseq/bam/ce_atacseq_ctrl_rep1.bam",
             "inst/files/atacseq/bam/ce_atacseq_ctrl_rep2.bam",
             "inst/files/atacseq/bam/ce_atacseq_ko_rep1.bam",
             "inst/files/atacseq/bam/ce_atacseq_ko_rep2.bam")

peakFile <- c("inst/files/atacseq/peaks/ce_atacseq_ctrl_rep1.narrowPeak.gz",
              "inst/files/atacseq/peaks/ce_atacseq_ctrl_rep2.narrowPeak.gz",
              "inst/files/atacseq/peaks/ce_atacseq_ko_rep1.narrowPeak.gz",
              "inst/files/atacseq/peaks/ce_atacseq_ko_rep2.narrowPeak.gz")

bigwigFile <- c("inst/files/atacseq/bigwig/ce_atacseq_ctrl_rep1.bigwig",
                "inst/files/atacseq/bigwig/ce_atacseq_ctrl_rep2.bigwig",
                "inst/files/atacseq/bigwig/ce_atacseq_ko_rep1.bigwig",
                "inst/files/atacseq/bigwig/ce_atacseq_ko_rep2.bigwig")

qcFile <- c("inst/files/atacseq/qc/ctrl/ctrl.qc.json",
            "inst/files/atacseq/qc/ctrl/ctrl.qc.json",
            "inst/files/atacseq/qc/ko/ko.qc.json",
            "inst/files/atacseq/qc/ko/ko.qc.json")

sampleDF_ <- data.frame(sampleName,bamFile,bigwigFile,peakFile,qcFile,condition,row.names = sampleName)
ce_atac <- ChrawExperimentFromDataFrame(sampleDF_,genome="hg38",minimalMetadata = FALSE)
save(ce_atac, file="data/ce_atac.RData",compress = TRUE)
############################################################################################################################################################################################
#CREATE CE_RAT
sampleName <- c("rat_ts_h3k4me1_1",
                "rat_ts_h3k27me3_1",
                "rat_ts_input_1")

condition <- c("h3k4me1",
          "h3k27me3",
          "input")

bamFile <- c("inst/files/rat/bam/rat_h3k4me1_1.bam",
             "inst/files/rat/bam/rat_h3k27me3_1.bam",
             "inst/files/rat/bam/rat_input_1.bam")

peakFile <- c("inst/files/rat/peaks/rat_h3k4me1_1.narrowPeak.gz",
              "inst/files/rat/peaks/rat_h3k27me3_1.narrowPeak.gz",
              "inst/files/rat/peaks/rat_input_1.narrowPeak.gz")

bigwigFile <- c("inst/files/rat/bigwig/rat_h3k4me1_1.bigwig",
                "inst/files/rat/bigwig/rat_h3k27me3_1.bigwig",
                "inst/files/rat/bigwig/rat_input_1.bigwig")

qcFile <- c("inst/files/rat/qc/h3k4me1/h3k4me1.qc.json",
            "inst/files/rat/qc/h3k27me3/h3k27me3.qc.json",
            "inst/files/rat/qc/input/input.qc.json")

sampleDF_ <- data.frame(sampleName,bamFile,bigwigFile,peakFile,qcFile,condition,row.names = sampleName)
ce_rat <- ChrawExperimentFromDataFrame(sampleDF_,genome="rn6",minimalMetadata = FALSE)


Peaks <- importNarrowPeaks( ce_rat,
                            includeSamples=rownames(colData(ce_rat)), merge=TRUE )
names(Peaks) <- sprintf("Peak%0.5d", seq_len(length(Peaks)))

ce_rat <- addCountExperiment(
  ce_rat, regions=Peaks,
  includeSamples=rownames(colData(ce_rat)),
  name="Peaks", BPPARAM=MulticoreParam(2))

save(ce_rat, file="data/ce_rat.RData",compress=TRUE)
##################################################################################################################################################################################################
#CREATE CE_EXAMPLES
ce_examples <- ce_chipseq

metaDat <- t(data.frame(strsplit(colData(ce_examples)$condition, "_|\\+")))
rownames(metaDat) <- NULL
metaDat <- as.data.frame(metaDat)
colnames(metaDat) <- c("donor", "target", "condition")
colData(ce_examples) <- cbind(colData(ce_examples), DataFrame( metaDat ))
colData(ce_examples)$sample_alias <- colData(ce_examples)$sampleName
colData(ce_examples)$condition <- c("CTRL_3h","CTRL_3h",
                                   "CTRL_27h","CTRL_27h",
                                   "agonist_3h","agonist_3h",
                                   "agonist_27h","agonist_27h")
colnames(colData(ce_examples))[10] <- "Time"
ce_examples$sample_alias <- gsub("H3K27ac_NR1I3_","",ce_examples$sample_alias)

ce_examples <- addFragmentLengthDist(
  ce_examples,
  param=csaw::readParam(pe = "both", restrict="chr6"),
  BPPARAM=MulticoreParam(2))

Samples <- rownames(colData(ce_examples))

Peaks <- importNarrowPeaks( ce_examples,
                            includeSamples=Samples, merge=TRUE )

#Peaks <- Peaks %>% plyranges::filter(seqnames == 'chr6', start >= 67090000, end <= 67170000)

names(Peaks) <- sprintf("Peak%0.5d", seq_len(length(Peaks)))

ce_examples <- addCountExperiment(
  ce_examples, regions=Peaks,
  includeSamples=Samples,
  name="Peaks", BPPARAM=MulticoreParam(2))

ce_examples <- testForDiffSignal(
  ce_examples, experimentName="Peaks",
  design=~condition,
  contrasts=list(agonist_3h=c("condition", "agonist_3h", "CTRL_3h"),
                 agonist_27h=c("condition","agonist_27h","CTRL_27h")))

ce_examples <- annotateExperimentRegions( ce_examples, "Peaks" )

save(ce_examples, file="data/ce_examples.RData",compress=TRUE)
##################################################################################################################################################################################################
tools::checkRdaFiles("data/")
tools::resaveRdaFiles("data/", compress = "xz")
tools::checkRdaFiles("data/")
