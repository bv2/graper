#' CLL multi-omic data
#'
#' A list containing differen omic measurements in CLL patient samples
#'
#' @format A list containing omic measurements for 92 patients
#' \itemize{
#'   \item{expr_lincRNA}{expression of lincRNA, p=1566}
#'   \item{expr_mRNA}{expression of mRNA, p=10910}
#'   \item{expr_mRNA}{expression of miRNA, p=1096}
#'   \item{met_3utr}{methylation in 3utr region, p=8661}
#'   \item{met_5utr}{methylation in 5utr region, p=8661}
#'   \item{met_cds_genebody}{methylation in coding region, p=9265}
#'   \item{met_intergenic}{methylation in intergenic region, p=11785}
#'   \item{met_noncds_genebody}{methylation in non-coding region, p=6884}
#'   \item{met_prom2k}{methylation in promotor region, p=12007}
#'   \item{surv1}{drug screen viability data, highest concentration, p=64}
#'   \item{surv2}{drug screen viability data, second-highest concentration, p=64}
#'   \item{surv3}{drug screen viability data, intermediate concentration, p=64}
#'   \item{surv4}{drug screen viability data, second-lowest concentration, p=64}
#'   \item{surv5}{drug screen viability data, lowest concentration, p=64}
#' }
#' @name CLLOmics
#' @usage data(CLLOmics)
NULL

# inputdir<-'~/Documents/LassoVariants/DataCompendium/CLL_views/R' allfiles<-list.files(inputdir)
# allfiles<-allfiles[grepl('.rds',allfiles)] CLLOmics<-lapply(allfiles, function(file) { print(file)
# readRDS(file.path(inputdir,file)) }) #mislabeled mRNA and linRNA adn overlapping miRNA and lncRNa
# load('/Users/bvelten/Documents/cll/Preprocessing_omics/proc170116/expr_miRNA.RData')
# load('/Users/bvelten/Documents/cll/Preprocessing_omics/proc170116/expr_mRNA.RData')
# load('/Users/bvelten/Documents/cll/Preprocessing_omics/proc170116/expr_lincRNA.RData') names(CLLOmics)<-sapply(allfiles,
# function(str) sub('.rds', '',str)) sapply(CLLOmics, dim) CLLOmics$expr_lincRNA<-lincRNA[rownames(CLLOmics$surv1),]
# CLLOmics$expr_mRNA<-mRNA[rownames(CLLOmics$surv1),] CLLOmics$expr_miRNA<-miRNA[rownames(CLLOmics$surv1),] sapply(CLLOmics,
# dim) save(CLLOmics, file='~/svn/huber/users/bvelten/LassoVariants/grpRR/data/CLLOmics.RData')


#' GSDC multi-omic data
#'
#' A Rdata object containing differen omic measurements in cancer cell lines (see Iorio et al)
#'
#' @format Including for n=962 samples
#' \itemize{
#'   \item{auc}{Aera under the curve for p=265 drugs}
#'   \item{BEM}{Binary Event Matrix, p=1250 (hypermethylation, mutations, CNA)}
#'   \item{CellLinesCommon}{Info on cell lines}
#'   \item{DrugMeta}{Info on drugs}
#'   \item{expr}{RMA normalized expression, p=17737}
#'   \item{exprgenesMeta}{Info on Genes in Expression data}
#'   \item{ic50}{IC50 values for p=265 drugs}
#'   \item{speedScoreBinary}{Speed Score from expression data, p=22}
#'   \item{X}{joint design matrix, p=19009}
#'   \item{annot}{annotation to 5 groups 'mut'   'meth'  'CNA'   'expr'  'speed' }
#' }
#' @name gdsc1000
#' @usage data(gdsc1000)
#'
NULL


# created in datagdsc1000.Rmd
