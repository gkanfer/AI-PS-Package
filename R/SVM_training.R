#! Do we need this function? #!

#' 2.1 Measures cell features
#' @param
#' @return
#' @export

extractFeatures <- function(x, mask_cyto, label_class="Positive") {
  table_shape = computeFeatures.shape(mask_cyto,x)
  table_moment = computeFeatures.moment(mask_cyto,x)
  table_basic = computeFeatures.basic(mask_cyto,x)
  table_test <- as.data.frame(cbind(table_basic,table_moment,table_shape))
  #! what is purpose of rownames, table restructuring ?
  rownameTable<-row.names(table_test)
  feature_table<-data.frame(cbind(rownameTable,table_test))
  Ts.mix<-feature_table[,2:20]
  rowNameTable<-feature_table[,1]
  Ts.mix$predict<-label_class
  Features_Table<-list()
  Features_Table[["Ts.mix"]]<-Ts.mix
  Features_Table[["rowNameTable"]]<-rowNameTable
  return(Features_Table)
}

#' 2.2 Choose cells for classification
#' @param
#' @return
#' @export
pickCells <- function(mask_nuc, x, Ts.mix, intens, param) {

}
