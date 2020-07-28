## 3 CNN Training

#! Should CNN training have its own pickCells (classification-based selection) function or share the same one from the SVM section? 

#' 3.1 Generation and upload of single cell images for CNN training
#' @param x Image 
#' @param y Cell segmentation mask generated from segmentCyto
#' @param pheno0 String. Name of control (normal) phenotype
#' @param pheno1 String. Name of positive phenotype
#' @return 
#' @export
genImageSet <- function(x, y, base_dir, pheno0, pheno1){
  #Break the 20x magnification image into single cell images
  stack<-EBImage::stackObjects(y, x*intens,ext = c(satck_size, satck_size)) #! do we need the ext option? Where does satck_size come from?
  
  #Save individual files in phenotype dir
  
  
  #! I don't think we should include this preprocessing function. Move to generator creation if we still want that option
  # Fimg <- function(x) {
  #   x1<-readImage(x, type="tiff")
  #   dx <- (dim(x1)[1])
  #   dx1 <- dx%/%8
  #   ft <- matrix(c(-0.05,-0.02,0, -0.02, -0.05), nrow=5, ncol=5)
  #   ft[5,5] <- 2
  #   y1 <- filter2(x1, ft)
  #   y2 <- y1[(dx1-10):(dx-dx1),(dx1-10):(dx-dx1)]
  #   writeImage(y2, paste(x),type = "tiff", bits.per.sample = 8L, compression = "LZW")
  # }
  
}

#' 3.2 Create testing, training directories (split based off of user-inputed ratio) and move files to appropriate destinations
#' Create files into testing and training sets, create test directory. Remaining files in original base_dir=train_dir 
#' @param base_dir 
#' @param pheno0 String. Name of control (normal) phenotype
#' @param pheno1 String. Name of positive phenotype
#' @param test_ratio Number between 0-1. Proportion of images to relocate to testing directory
#' @export
splitDirs <- function(base_dir, pheno0, pheno1, test_ratio) {
  cls_list = c(pheno0, pheno1)
  for (cls in cls_list) {
    createdir(file.path(base_dir,'test',cls))
    createdir(file.path(base_dir,'train',cls))
    src = file.path(base_dir, cls) 
    allfiles = list.files(path=src, full.names=TRUE)
    allfiles = sample(allfiles)  #random shuffling 
    test_n= as.integer(length(allfiles)*test_ratio)
    test_files=allfiles[0:test_n]
    train_files=allfiles[test_n:as.integer(length(allfiles))]
    print(paste0("Testing-",cls,": ",length(test_files)))
    print(paste0("Training-",cls,": ",length(train_files)))
  }
  for (i in test_files) {
    path = file.path(base_dir,"test",cls) 
    file.copy(i, path)
  }
}

##Internal Functions##
#' #' Directory create
#' createdir <- function(path, subpath) {
#'   dir <- file.path(path, subpath)
#'   dir.create(dir)
#'   return(dir)
#' }
