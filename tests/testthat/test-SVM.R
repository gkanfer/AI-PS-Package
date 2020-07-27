test_that("Extract features works", {
  img = LoadImage("~/Documents/test_tifs/HAB135_WT_Tom20-CAT-dapi_1.tif")
  nseg = segmentNucleus(img, index=1, minmaxnorm=TRUE,
                        int=2, filter_size=25, offset=0.01, opensize=3,
                        small_obj=30, use_watershed=FALSE, distmap_value=2,
                        rm_outliers=FALSE, out_p=0.95, displaymasks=TRUE)
  cseg = segmentCyto(img, nseg$seg, index=3, int=10, filter_size=100,
                     offset=0.01, size_smooth=19, opensize=7, largeobj=100000, minmaxnorm=TRUE)
  table = extractFeatures(img, cseg$seg)
  expect(nrow(table$table)>1)
})