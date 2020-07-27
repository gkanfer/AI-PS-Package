#' 2.1 Measures cell features
#' @param
#' @return
#' @export

extractFeatures <- function(x, mask_cyto, label_class="Positive") {
  table_shape = EBImage::computeFeatures.shape(mask_cyto,x)
  table_moment = EBImage::computeFeatures.moment(mask_cyto,x)
  table_basic = EBImage::computeFeatures.basic(mask_cyto,x)
  table_test <- as.data.frame(cbind(table_basic,table_moment,table_shape))
  #! what is purpose of rownames, table restructuring ?
  rownameTable<-row.names(table_test)
  table_test$predict<-label_class
  feature_table<-data.frame(cbind(rownameTable,table_test))
  #Ts.mix<-feature_table[,2:20]
  rowNameTable<-feature_table[,1]
  #Ts.mix$predict<-label_class
  #Features_Table<-list()
  #Features_Table[["Ts.mix"]]<-Ts.mix
  #Features_Table[["rowNameTable"]]<-rowNameTable
  return(list(table=feature_table, names=rowNameTable))
}

#' 2.2 Choose cells for classification
#' @param
#' @return
#' @export
pickCells<-function(mask_nuc, xy, Ts.mix,int, parameters_nuc, font_size=0.7, label_class="Postive", Mdisplay_select=TRUE) {
  seg_disp = paintObjects(mask_nuc, toRGB(x*intens),opac=c(1,0.8),col=c("Green",NA),thick=TRUE,closed=FALSE)
  display(seg_disp,"raster")
  celltext = text(x=xy[,1], y= parameters_nuc[,2], labels="", col="yellow", cex=font_size)
  c<-0
  readline(paste0("Select", label_class, "cells"))
  temp<-locator()
  c<-c(c,temp)
  xy<-computeFeatures.moment(mask_nuc)[,c('m.cx','m.cy')]
  df.c <- cbind(c$x,c$y)
  knn.out <- ann(as.matrix(xy), as.matrix(df.c), k=2)
  row_numb<-knn.out$knnIndexDist
  class(row_numb)
  row_numb<-as.data.frame(row_numb)
  Ts.mix$predict<-0
  if (label_class=="Postive"){
    Ts.mix[row_numb$V1,20]<-"P"
    Ts.mix[-row_numb$V1,20]<-"N"
  }else{
    Ts.mix[row_numb$V1,20]<-"N"  
    Ts.mix[-row_numb$V1,20]<-"P"
  }
  nr = which(Ts.mix$predict %in% "N")
  
  Table_class<-list()
  if (display_select) {
    seg.temp = rmObjects(mask_nuc, nr)
    seg_class_sel = paintObjects(seg.temp,toRGB(x*intens),opac=c(1,0.8),col=c("Green",NA),thick=TRUE,closed=FALSE)
  }
  return(list(table_train=Ts.mix, segsel=seg_class_sel))
} 

#'
