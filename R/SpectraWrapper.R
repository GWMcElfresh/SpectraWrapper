#' @title RunSpectra
#' @description Wrapper function for the Spectra package
#' @param seuratObj Input Seurat object
#' @param seuratToAdataDir Directory to house files used in generating Ann Data
#' @param assayName Assay used for generating Ann Data
#' @param inputJSONfile Input Seurat object
#' @param immportfile Input Seurat object

RunSpectra <- function(seuratObj, seuratToAdataDir, assayName, inputJSONfile, immportfile){

  SeuratToAData(seuratObj, seuratToAdataDir, assayName)
  gex_datafile = paste0(seuratToAdataDir, "/GEX.h5")
  metadata_file = paste0(seuratToAdataDir, "/metadata.tsv")
  genesets_json = paste0(seuratToAdataDir, "/genesets.json")
  variable_features_file = paste0(seuratToAdataDir, "varfeats.csv")

  ProduceGenesetJSONfile(inputJSONfile, immportfile, genesets_json)

  #copy run_Spectra.py in inst/scripts and supply custom arguments
  str <- readr::read_file(system.file("scripts/run_Spectra.py", package = "CellMembrane"))
  script <- tempfile()
  readr::write_file(str, script)

  # run_Spectra(features_file, metadata_file, gex_datafile, genesets_json)
  newstr <- paste0("run_Spectra(features_file = '", variable_features_file,
                   "', gex_datafile = '", gex_datafile,
                   "', metadata_file = '", metadata_file,
                   "', genesets_json = '", genesets_json,
                   "', seuratToAdataDir = '", seuratToAdataDir,
                   "')")

  #write the new arguments into the script and execute
  readr::write_file(newstr, script, append = TRUE)
  system2(reticulate::py_exe(), script)
  seuratObj <- SpectraToSeuratReduction(seuratObj, seuratToAdataDir, assayName, reduction.name = 'spectra', reduction.key = 'Spectra_')
  seuratObj <- SpectraToSeuratMetadata(seuratObj, paste0(seuratToAdataDir, "/SPECTRA_cell_scores.csv"))
  return(seuratObj)
}


#' @title SeuratToAData
#' @param seuratObj The Seurat object to be written into a 10x file format
#' @param seuratToAdataDir The local path to write the output files.
#' @param assayName The name of the gene expression data assay.
#' @description A wrapper function to prepare a Seurat object for Conga.
#' @export
SeuratToAData <- function(seuratObj,
                          seuratToAdataDir,
                          assayName = 'RNA') {
  if (!dir.exists(seuratToAdataDir)) {
    dir.create(seuratToAdataDir, recursive = T)
  }
  seuratObj <- Seurat::FindVariableFeatures(seuratObj)
  write.table(Seurat::VariableFeatures(seuratObj),
              R.utils::getAbsolutePath(paste0(seuratToAdataDir, "/varfeats.csv"),
                                       mustWork = FALSE), row.names = FALSE, col.names = FALSE)

  DropletUtils::write10xCounts(x = Seurat::GetAssayData(seuratObj, assay = assayName, slot = 'data'),
                               path = R.utils::getAbsolutePath(paste0(seuratToAdataDir, "/GEX.h5"), mustWork = FALSE),
                               overwrite = TRUE)
  data.table::fwrite(seuratObj@meta.data,
                     R.utils::getAbsolutePath(paste0(seuratToAdataDir, "/metadata.tsv"),
                                              mustWork = FALSE), row.names = TRUE, col.names = TRUE)

}


ProduceGenesetJSONfile <- function(inputpath, inputpath2, outputpath){
  # read in gene sets
  jsonin <- jsonlite::read_json(inputpath, simplifyVector = TRUE)
  immport <- readRDS(inputpath2)
  genesets <- RIRA::ListGeneSets()

  geneset_auto <- list()
  geneset_auto$global <- jsonin$global
  for (i in c(1:length(genesets))){
    # to use cell-type-specific genesets, replace global here with the cell-type label, e.g.:
    # geneset_auto$T_NK[[genesets[i]]] <- RIRA::GetGeneSet(genesets[i])
    geneset_auto$global[[genesets[i]]] <- RIRA::GetGeneSet(genesets[i])
  }
  for (i in 1:length(immport)){
    geneset_auto$global[[names(import[i])]] <- immport[[i]]
  }
  jsonlite::write_json(geneset_auto, outputpath)
}

SpectraToSeuratMetadata <- function(seuratObj, cellresults_file, plotComponents = TRUE){
  # Note: this adds NAs for missing cells. We could in theory change this to zeros if NAs are a problem.
  SpectraScores <- data.table::fread(cellresults_file) |> as.data.frame()
  rownames(SpectraScores) <- SpectraScores$V1
  SpectraScores <- SpectraScores |> select(-"V1")
  colnames(SpectraScores) <- paste0('Spectra_', 1:ncol(SpectraScores))

  for (cmp in  colnames(SpectraScores)){
    seuratObj <- Seurat::AddMetaData(seuratObj, SpectraScores[,cmp], col.name = cmp)
    if (plotComponents){
      print(Seurat::FeaturePlot(seuratObj, features = cmp, order = T) & ggplot2::scale_colour_gradientn(colours = c('navy', 'dodgerblue', 'white', 'gold', 'red')))
    }
  }
  return(seuratObj)
}

SpectraToSeuratReduction <- function(seuratObj, seuratToAdataDir, assayName = 'RNA', reduction.name = 'spectra', reduction.key = 'Spectra_') {
  SpectraScores <- data.table::fread(paste0(seuratToAdataDir, "/SPECTRA_cell_scores.csv"), header = T) |> as.data.frame()
  rownames(SpectraScores) <- SpectraScores$V1
  SpectraScores <- SpectraScores |> select(-"V1")
  colnames(SpectraScores) <- paste0('Spectra_', colnames(SpectraScores))

  geneweights <- data.table::fread(paste0(seuratToAdataDir, "/geneweights.csv"), header = T) |> as.data.frame()
  rownames(geneweights) <- geneweights$V1
  geneweights <- geneweights |> select(-V1)

  overlap <- data.table::fread(paste0(seuratToAdataDir, "/SPECTRA_overlap.csv"), header = T) |> as.data.frame()
  rownames(overlap) <- overlap$V1
  overlap <- overlap |> select(-V1)

  markers <- data.table::fread(paste0(seuratToAdataDir, "/SPECTRA_markers.csv"), header = T) |> as.data.frame()
  rownames(markers) <- markers$V1
  markers <- markers |> select(-V1)

  jsonin <- jsonlite::read_json(paste0(seuratToAdataDir,"/SPECTRA_L.json"), simplifyVector = TRUE)

  embeddings <- SpectraScores
  colnames(embeddings) <- paste0(reduction.key, 1:ncol(embeddings))

  # Note: since SDA could drop cells, add back in the missing cells with zeros
  extraCells <- setdiff(rownames(embeddings), colnames(seuratObj))
  if (length(extraCells) > 0) {
    stop(paste0('There were ', length(extraCells), ' with data in the Spectra results but not present in the seurat object.  Top barcodes: ', paste0(head(extraCells), collapse = ',')))
  }

  missingCells <- setdiff(colnames(seuratObj), rownames(embeddings))
  if (length(missingCells) > 0) {
    toAdd <- matrix(rep(0, ncol(embeddings)*length(missingCells)), ncol = ncol(embeddings))
    colnames(toAdd) <- colnames(embeddings)
    rownames(toAdd) <- missingCells
    embeddings <- rbind(embeddings, toAdd)
  }

  embeddings <- embeddings[colnames(seuratObj),] |> as.matrix()

  # rows = features. This needs to use the projected slot, since it probably has fewer features than the assay
  loadings <- t(geneweights)
  colnames(loadings) <- paste0(reduction.key, 1:ncol(loadings))

  # See: https://satijalab.org/seurat/archive/v3.0/dim_reduction_vignette.html
  spectra.reduction <- Seurat::CreateDimReducObject(
    embeddings = embeddings,
    loadings  = loadings,
    key = reduction.key,
    assay = assayName
  )

  SpectraResults <- list()
  SpectraResults$overlap <- overlap
  SpectraResults$markers <- markers
  SpectraResults$jsonin <- jsonin

  for (x in c('overlap', 'markers', 'jsonin')) {
    if (x %in% names(SpectraResults)) {
      spectra.reduction@misc[[x]] <- SpectraResults[[x]]
    }
  }

  seuratObj[[reduction.name]] <- spectra.reduction

  return(seuratObj)
}
