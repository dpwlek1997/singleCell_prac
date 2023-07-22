# This script had been generated for ArchRtoSignac with more varaible Matrixes not just peakMatrix and GeneScoreMatrix.
# Such as geneExpressionMatrix, tileMatrix and other matrixes. The codes were modifed from https://github.com/swaruplabUCI/ArchRtoSignac/blob/main/R/ArchRtoSignac.R. 

#############################
# Install and load packages #
#############################
devtools::install_github("swaruplabUCI/ArchRtoSignac")
packages <- c("ArchRtoSignac", "ArchR","Seurat", "Signac","stringr") # required packages
loadinglibrary(packages)

#################################################################
# Step 1: obtain ArchRproject peak matrix for object conversion #
#################################################################
pkm <- getPeakMatrix(proj)

#####################################################################
# Step 2: extract appropriate ensembl gene annotation -> UCSC style #
#####################################################################
library(EnsDb.Hsapiens.v86)
annotations <- getAnnotation(reference = EnsDb.Hsapiens.v86, refversion = "hg38")

#######################################################
# Step 3: convert ArchRproejct to Signac SeuratObject #
#######################################################
## peak added
fragments_dir <- "../multiome/fragments/"
seurat_object <- ArchR2Signac(  # modified ArchR2Signac function a bit, please refer line ___. 
  ArchRProject = proj3_integ_90,
  refversion = "hg38",
  fragments_dir = fragments_dir,
  pm = pkm, 
  fragments_fromcellranger = "Yes", 
  annotation = annotations 
)
seurat_object

##############################
# Step 4: add other matrixes #
############################## 
# Here, multiome data was used. So ArchRproject has GeneExprssionMatirx, GeneScoreMatrix, GeneIngtegrationMatrix.
## GeneExpressionMatirx: gene expression from snRNAseq of multiome
## GeneScoreMatrix: assumed gene expression based on chromatin accessibility from snATACseq of multiome
## GeneIntegrationMatrix: 

get_another_Matrix <- function(ArchRProject, SeuratObject, matrix_name){     # modified of getGeneScoreMatrix() in ArchRtoSignac
  matrix <- ArchR::getMatrixFromProject(ArchRProject, useMatrix=matrix_name)
  assay <- assays(matrix)[[1]]
  
  GeneFeatures <- getFeatures(
    ArchRProj = ArchRProject,
    useMatrix = matrix_name,
    select = NULL,
    ignoreCase = TRUE
  )
  
  # set the column/row names of assay so that it matches with the Seurat object 
  colnames(assay) <- gsub(".*#", "", colnames(assay)) # i had samplename#barcode-1, so removed the sample name
  ix <- match(colnames(SeuratObject), colnames(assay))
  assay <- assay[,ix]
  rownames(assay) <- GeneFeatures
  
  return(assay)
}

# check which matrix you got
getAvailableMatrices(ArchRProject)

# all expression matrix 
assay_GeneExpressionMatrix <- get_another_Matrix(ArchRProject, seurat_object, "GeneExpressionMatrix")
assay_GeneIntegrationMatrix <- get_another_Matrix(ArchRProject, seurat_object, "GeneIntegrationMatrix")
assay_GeneScoreMatrix <- get_another_Matrix(ArchRProject, seurat_object, "GeneScoreMatrix")

# add to seurat object 
seurat_object[['RNA']] <- CreateAssayObject(counts = assay_GeneExpressionMatrix) # 36438
seurat_object[["integrated"]] <- CreateAssayObject(counts = assay_GeneIntegrationMatrix) #  20397 
seurat_object[['genescore']] <- CreateAssayObject(counts = gsm)

###################################
# Step 5: add dimension reduction #
###################################




# modified ArchR2Signac #
ArchR2Signac <- function(
  ArchRProject,
  refversion, # write the EnsDb version
  samples = NULL, # Provide a list of unique sample
  fragments_dir = NULL, # directory of the cellranger output, the folder that contains all samples
  pm, # geting peak martix
  fragments_fromcellranger = NULL, # "NO" | "N" | "No" or "YES" | "Y" | "Yes"
  fragments_file_extension = NULL, #  '.tsv.gz' or '.fragments.tsv.gz'
  # output_dir = '/outs/', # removal due to the input format for snapATAC, added when fragments_fromcellranger == "NO" | "N" | "No"
  annotation # annotation from getAnnotation()
 ){
   if (is.null(samples)){
     samples <- unique(ArchRProject@cellColData$Sample)
   }

   if(fragments_fromcellranger == "YES" | fragments_fromcellranger == "Y" | fragments_fromcellranger == "Yes") {
     print("In Progress:")
     print("Prepare Seurat list for each sample")

     output_dir = '/outs/'
      
     seurat_list <- lapply(samples, function(cur_sample){
       print(cur_sample)
       #print out the sample name in progress
       
       print(is.list(fragments_dir))
       print(paste0(fragments_dir, cur_sample, output_dir, 'fragments.tsv.gz'))
       cur_fragments <- ifelse(is.list(fragments_dir),
                               paste0(fragments_dir[[which(samples == cur_sample)]], output_dir, 'atac_fragments.tsv.gz'), # yeji added atac for my case
                               paste0(fragments_dir, cur_sample, output_dir, 'atac_fragments.tsv.gz'))

       # seeking the pattern matched in colnames(pm); metadata of the corresponding sample
       cur_pm <- pm[,grepl(paste0(cur_sample, '#'), colnames(pm))]
       cur_meta <- ArchRProject@cellColData %>% as.data.frame %>% subset(Sample == cur_sample)

       # change colnames and rowname format:
       colnames(cur_pm) <- do.call(rbind, str_split(colnames(cur_pm), '#'))[,2]
       rownames(cur_meta) <- do.call(rbind, str_split(rownames(cur_meta), '#'))[,2]
       print(dim(cur_pm))
       # create chromatin assay
       cur_chromatin <- Signac::CreateChromatinAssay(
         counts=cur_pm, # should we add data instead counts
         sep = c('-', '-'),
         fragments=cur_fragments, # do we need this?
         ranges=ArchRProject@peakSet,
         genome=refversion,
         annotation = annotation
       )

       # create a new Seurat obj with only the archR peaks:
       cur_atac <- Seurat::CreateSeuratObject(
         cur_chromatin,
         assay='peaks',
         meta.data = cur_meta,
       )

     })
   }

   if(fragments_fromcellranger == "NO" | fragments_fromcellranger == "N" | fragments_fromcellranger == "No") {

     print("IF selecting NO, please make sure to provide fragments_file_extension")
     print("In Progress:")
     print("Prepare Seurat list for each sample")

     seurat_list <- lapply(samples, function(cur_sample){
       print(cur_sample)
       #print out the sample name in progress
       cur_fragments <- ifelse(is.list(fragments_dir),
                        paste0(fragments_dir[[which(samples == cur_sample)]], fragments_file_extension),
                        paste0(fragments_dir, cur_sample, fragments_file_extension))
       # seeking the pattern matched in colnames(pm); metadata of the corresponding sample
       cur_pm <- pm[,grepl(paste0(cur_sample, '#'), colnames(pm))]
       cur_meta <- ArchRProject@cellColData %>% as.data.frame %>% subset(Sample == cur_sample)

       # change colnames and rowname format:
       colnames(cur_pm) <- do.call(rbind, str_split(colnames(cur_pm), '#'))[,2]
       rownames(cur_meta) <- do.call(rbind, str_split(rownames(cur_meta), '#'))[,2]
       print(dim(cur_pm))
       # create chromatin assay
       cur_chromatin <- Signac::CreateChromatinAssay(
         counts=cur_pm, # should we add data instead counts
         sep = c('-', '-'),
         fragments=cur_fragments, # do we need this?
         ranges=ArchRProject@peakSet,
         genome=refversion,
         annotation = annotation
       )

       # create a new Seurat obj with only the archR peaks:
       cur_atac <- Seurat::CreateSeuratObject(
         cur_chromatin,
         assay='peaks',
         meta.data = cur_meta,
       )

     })

   }

   print("In Progress:")
   
   if (length(seurat_list) == 1){
     SeuratObject <- seurat_list[[1]]
   }else{
      SeuratObject <- merge(
        print("Merge Seurat list")
        x=seurat_list[[1]], y=seurat_list[2:length(seurat_list)],
        add.cell.ids = samples
      )

      print("Return SeuratObject")
   }
   
   SeuratObject
}






