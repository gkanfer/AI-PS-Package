## 3 CNN Training
#' 3.1 Generation and upload of single cell images for CNN training
#' @param y Cell segmentation mask generated from segmentCyto
#' @param x Image 
#' @param 
#' @return 
#' @export
genImageSet <- function(y, x, base_dir, train_dir, validation_dir, test_dir, train_norm_dir, train_swol_dir, validation_norm_dir, original_dataset_dir_norm, original_dataset_dir_swol, tr.mx, val.mx, tst.mx){
  #Break the 20x magnification image into single cell images
  stack<-stackObjects(y,x*intens,ext = c(satck_size,satck_size))
  #Folders creation 
  train_dir <- createdir(base_dir, "train")
  val_dir <- file.path(base_dir, "val")
  dir.create(train_dir)
  validation_dir <- file.path(base_dir, "validation")
  dir.create(validation_dir)
  test_dir <- file.path(base_dir, "test")
  dir.create(test_dir)
  train_norm_dir <- file.path(train_dir, "Norm")
  dir.create(train_norm_dir)
  train_swol_dir <- file.path(train_dir, "Swol")
  dir.create(train_swol_dir)
  validation_norm_dir <- file.path(validation_dir, "Norm")
  dir.create(validation_norm_dir)
  validation_swol_dir <- file.path(validation_dir, "Swol")
  dir.create(validation_swol_dir)
  test_swol_dir <- file.path(test_dir, "Swol")
  dir.create(test_swol_dir)
  test_norm_dir <- file.path(test_dir, "Norm")
  dir.create(test_norm_dir)
  # Transfer all images to the apropriete destenation
  setwd(original_dataset_dir_norm)
  fnames <- dir()[grep(".tif",dir())]
  fnames = sample(fnames)
  length(fnames)
  file.copy(file.path(original_dataset_dir_norm, fnames[1:(tr.mx+100)]),
            file.path(train_norm_dir))
  file.copy(file.path(original_dataset_dir_norm, fnames[(tr.mx+100+1):val.mx+200]),
            file.path(validation_norm_dir))
  file.copy(file.path(original_dataset_dir_norm, fnames[(val.mx+200+1):tst.mx+200]),
            file.path(test_norm_dir))
  setwd(original_dataset_dir_swol)
  fnames <- dir()[grep(".tif",dir())]
  fnames = sample(fnames)
  length(fnames)
  file.copy(file.path(original_dataset_dir_swol, fnames[1:tr.mx]),
            file.path(train_swol_dir))
  file.copy(file.path(original_dataset_dir_swol, fnames[(tr.mx+1):val.mx]),
            file.path(validation_swol_dir))
  file.copy(file.path(original_dataset_dir_swol, fnames[(val.mx+1):(tst.mx)]),
            file.path(test_swol_dir))
  
  Fimg <- function(x) {
    x1<-readImage(x, type="tiff")
    dx <- (dim(x1)[1])
    dx1 <- dx%/%8
    ft <- matrix(c(-0.05,-0.02,0, -0.02, -0.05), nrow=5, ncol=5)
    ft[5,5] <- 2
    y1 <- filter2(x1, ft)
    y2 <- y1[(dx1-10):(dx-dx1),(dx1-10):(dx-dx1)]
    writeImage(y2, paste(x),type = "tiff", bits.per.sample = 8L, compression = "LZW")
  }
  
  setwd(train_dir)
  setwd(paste(getwd(),"Swol", sep="/"))
  #dir()
  tifsS<-dir()[grep(".tif",dir())]
  lapply(tifsS, Fimg)
  
  setwd(train_dir)
  setwd(paste(getwd(),"Norm",sep="/"))
  #dir()
  tifsS<-dir()[grep(".tif",dir())]
  lapply(tifsS, Fimg)
  
  setwd(validation_dir)
  setwd(paste(getwd(),"Swol",sep="/"))
  #dir()
  tifsS<-dir()[grep(".tif$",dir())]
  lapply(tifsS, Fimg)
  
  setwd(validation_dir)
  setwd(paste(getwd(),"Norm",sep="/"))
  #dir()
  tifsS<-dir()[grep(".tif$",dir())]
  lapply(tifsS, Fimg)
}

##Internal Functions##
#' Directory create
createdir <- function(path, subpath) {
  dir <- file.path(path, subpath)
  dir.create(dir)
  return(dir)
}