# AI-PS Packages

# - -----------------------------------------------------------------------
#' 1. Segmentation
#' 1.1 Nucleus segmentation 
#' 1.2 Cytosol segmentation 
#' 2. SVM training
#' 2.1 Features measures 
#' 2.2 Cell picking
#' 2.3 SVM model device
#' 2.4 Features reduction using PCA analysis
#' 2.5 SVM Model tuning 
#' 3. CNN training
#' 4. Deployment example
#' 4.1 SVM
#' 4.2 CNN
#' 5.confusion matrix and ROC curve 

# -------------------------------------------------------------------------
# 1. Segmentation
# -------------------------------------------------------------------------

# Image load
load_image<-function(x,path,minmax_nrom=TRUE,...){
          x<-readImage()
          dim.y<-dim(x)[1]
          dim.X<-dim(x)[2]
          CH1<-x[1:dim.y,1:dim.X,1]
          CH2<-x[1:dim.y,1:dim.X,2]
          colorMode(CH1)<-"Grayscale"
          colorMode(CH2)<-"Grayscale"
          if (minmax_nrom){
              minCH1<-min(as.vector(CH1))
              maxCH1<-max(as.vector(CH1))
              minCH2<-min(as.vector(CH2))
              maxCH2<-max(as.vector(CH2))
              CH1<-normalize(CH1, ft=c(0,1),c(minCH1,maxCH1))
              CH2<-normalize(CH2, ft=c(0,1) ,c(minCH2,minCH2)) 
                          }
          channel_set<-list(CH1,CH2)
          return(channel_set)
          }


#' 1.1 Nucleus segmentation 
#'  
Nucleus_segmentation<-function(x,intensity=1,filter_size=16,offset=0.04,opensize=3,
                               rmObjects_small=30,Use_watershed=TRUE,distmap_value=2,rm_outliers=TRUE,Out_P=0.95,
                               Masks_display=TRUE){
          CH1_seg<- x*2
          CH1mask2 = thresh(CH1_seg, filter_size, filter_size,offset)
          CH1mask2 = opening(CH1mask2, makeBrush(opensize, shape= "diamond")) 
          CH1mask2 = fillHull(CH1mask2) 
          #extraction
          CH1mask2 = bwlabel(CH1mask2)  #binary to object
          if (Use_watershed){
            CH1mask2 = watershed( distmap(CH1mask2), distmap_value )
          } 
          chackpoint<-computeFeatures.shape(CH1mask2)
          nmask<-CH1mask2
          if (length(chackpoint) < 3 ){
            stop("no segmented object was detected, change setup parameters")
          }
          if (is.null(nrow(chackpoint))){
            stop("no segmented object was detected, change setup parameters")
          }
          if (nrow(chackpoint) > 2000 ){
            stop("poor segentation,change setup parameters")
          }
          nf = computeFeatures.shape(nmask)
          nr = which(nf[,2] < rmObjects_small) 
          nseg = rmObjects(nmask, nr) 
          nn = max(nseg) 
          colorMode(nseg)<-Grayscale
          chackpoint<-computeFeatures.shape(nseg)
          if (length(chackpoint) < 3 ){
            stop("no segmented object was detected, change setup parameters")
          }
          if (is.null(nrow(chackpoint))){
            stop("no segmented object was detected, change setup parameters")
          }
          if (nrow(chackpoint) > 2000 ){
            stop("poor segentation,change setup parameters")
          }
          #remove outliers 
          if (rm_outliers){
                  int.dapi<-computeFeatures.basic(nseg,x)
                  y<-which(scores(int.dapi[,1], type="z",prob = 0.95))
                  tp<-as.numeric(attr(y,"names"))
                  if (length(tp) < 1 ){
                    stop("Outlires detection failed")
                  }
                  nseg<-rmObjects(nseg,tp)
          }
          chackpoint<-computeFeatures.shape(nseg)
          df<-as.data.frame(chackpoint)
          xy<-computeFeatures.moment(nseg)[,c('m.cx','m.cy')]
          if (length(xy) < 3 ){
            stop("no segmented object was detected, change setup parameters")
          }
          if (is.null(nrow(xy))){
            stop("no segmented object was detected, change setup parameters")
          }
          if (nrow(xy) > 500 ){
            stop("poor segentation,change setup parameters")
          }
          gsegg=nseg
          mask<-list()
          mask["mask_nuc"]<-gsegg
          masks[["parameters_nuc"]]<-xy
          if (Masks_display){
                seg_CH1mask2 = paintObjects(CH1mask2,toRGB(x),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)  
                seg_mask<-paintObjects(gsegg,toRGB(x),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
                mask["mask_start"]<-seg_CH1mask2
                mask["mask_final"]<-seg_mask
          }
          return(masks)
          
}

# 1.2 GFP segmentation 
Cytosol_segmentation <- function (x,y,intensity=40,filter_size=10,offset=0.1,mat_size_smooth=19,
                                  opensize=7,rmObjects_large=30000,Use_watershed=TRUE,Masks_display=TRUE){
          Cyto<- x*40
          Cyto<-filter2(Cyto,makeBrush(mat_size_smooth,shape = "disc") , boundary = c("circular", "replicate"))
          thr<-thresh(Cyto, filter_size, filter_size, offset)
          colorMode(thr)<-Grayscale
          cmask = opening(thr, kern=makeBrush(opensize,shape="disc"))
          ar<-as.vector(Cyto)
          ar.sum<-as.numeric(summary(ar))
          open2<-opening(Cyto>(ar.sum[2]))
          colorMode(open2)<-Grayscale
          combine<-cmask
          combine[open2 > cmask]<-open2[open2 > cmask]
          combine[y > combine]<-y[y > combine]
          cseg = propagate(Cyto, y, lambda=1.0e-2, mask=cmask)
          cseg <- fillHull(cseg)
          colorMode(csegpink)<-Grayscale
          cseg = propagate(Cyto, gsegg, lambda=1.0e-2, mask=combine)
          cseg <- fillHull(cseg)
          colorMode(cseg)<-Grayscale
          xy<-computeFeatures.moment(csegpink)[,c('m.cx','m.cy')]
          if (is.null(nrow(xy))){
            stop("no segmented object was detected, change setup parameters")
          }
          if (nrow(xy) > 2000 ){
            stop("poor segentation,change setup parameters")
          }
          cf = computeFeatures.shape(cseg)
          cfErea<-data.frame(cf[,1])
          cfErea$num<-row.names(cfErea)
          ci = which(cf[,1] > rmObjects_large) 
          csegpink = rmObjects(cseg, ci,reenumerate = F) 
          xy.gsegg<-as.numeric(row.names(computeFeatures.moment(y)[,c('m.cx','m.cy')]))
          xy.cseg<- as.numeric(row.names(computeFeatures.moment(cseg)[,c('m.cx','m.cy')])) 
          ind.deff<-setdiff(xy.gsegg,xy.cseg)
          gsegg<-rmObjects(gsegg,ind.deff,reenumerate=F)
          gsegg<-reenumerate(gsegg)
          cseg<-reenumerate(cseg)
          masks<-list()
          masks[["mask_nuc"]]<-gsegg
          masks[["parameters_nuc"]]<-xy.gsegg
          masks[["mask_cyto"]]<-cseg
          masks[["parameters_cyto"]]<-xy.cseg
          if (Masks_display){
              seg_gsegg = paintObjects(CH1mask2,toRGB(Cyto),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)  
              seg_cseg<-paintObjects(gsegg,toRGB(Cyto),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
              mask["display_mask_nuc"]<-seg_gsegg
              mask["display_mask_cyto"]<-seg_mask
          }
          
          return(masks)
          }


# -------------------------------------------------------------------------
# 2. SVM training
# -------------------------------------------------------------------------
# 2.1 Features measures 
Feature_extract<- function(x,mask_cyto,label_class="Postive"){
          table_shape = computeFeatures.shape(mask_cyto,x)
          table_moment = computeFeatures.moment(mask_cyto,x)
          table_basic = computeFeatures.basic(mask_cyto,x)
          table_test<-data.frame(cbind(table_basic,table_moment,table_shape))
          rownameTable<-row.names(table_test_pink)
          table_test_pink<-data.frame(cbind(rownameTable,table_test_pink))
          Ts.mix<-table_test_pink[,2:20]
          rowNameTable<-table_test_pink[,1]
          Ts.mix$predict<-label_class
          Features_Table<-list()
          Features_Table[["Ts.mix"]]<-Ts.mix
          Features_Table[["rowNameTable"]]<-rowNameTable
          return(Features_Table)
}

# 2.2 Cell picking 
Pick_cells<-function(mask_nuc,x,Ts.mix,intens,parameters_nuc,font_size=0.7,label_class="Postive",
                     Masks_display){
          seg_disp = paintObjects(mask_nuc,toRGB(x*intens),opac=c(1,0.8),col=c("Green",NA),thick=TRUE,closed=FALSE)
          display(seg_pink_positive,"raster")
          celltext = text(x= parameters_nuc[,1], y= parameters_nuc[,2] , labels="", col="yellow", cex = font_size)
          c<-0
          readline("loop")
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
              Ts.mix[row_numb$V1,21]<-"P"
          }else{
              Ts.mix[row_numb$V1,21]<-"N"  
              }
          nr = which(Ts.training$predict %in% 0)
          
          Table_class<-list()
          if (Masks_display){
            gsegg.temp = rmObjects(mask_nuc, nr)
            seg_class_sel = paintObjects(gsegg.temp,toRGB(x*intens),opac=c(1,0.8),col=c("Green",NA),thick=TRUE,closed=FALSE)
            Table_class[["seg_class_sel"]]<-seg_class_sel
            
          }
          Table_class[["Table_class_train"]]<-Ts.mix
          return(Table_class)
        } 

# 2.3 SVM model device 

Model_svm<-function(Table_class_train,kernel_linear=TRUE,cost= 10, degree = 45){
          x<-(Table_class_train[,2:20])
          y<-Table_class_train[,21]
          
          if (kernel_linear){
              acc<-rep(0,100)
              for (i in 1:length(acc)){
                TestIndex<-sample(nrow(x),round(nrow(x)/2))
                model<-svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="linear",type = "C",cost= 10, degree = 45,probability = TRUE)
                y.pred<-predict(model,x[-TestIndex,])
                acc[i]=length(which(y.pred==y[-TestIndex]))/length(y.pred)
                }
              }else{
              acc<-rep(0,100)
              for (i in 1:length(acc)){
                TestIndex<-sample(nrow(x),round(nrow(x)/2))
                model<-svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="radial",type = "C",cost= 10, degree = 45,probability = TRUE)
                y.pred<-predict(model,x[-TestIndex,])
                acc[i]=length(which(y.pred==y[-TestIndex]))/length(y.pred)}
                }
          mean(acc)
          return(acc)
}
# 2.4 Features reduction using PCA analysis

Fetures_PCA<-function(Table_class_train,cartesian_lim_x=c(-12,10),cartesian_lim_y =c(-10,5),font_size=14){
          myPr <- prcomp(Table_class_train[, 2:19], scale = TRUE) 
          PCA_PLOT<-ggbiplot(myPr, pc.biplot = T,obs.scale = 2, var.scale = 1, groups = F, ellipse = F, circle = F,alpha = 0.02,varname.adjust = 5)+
                    coord_cartesian(xlim = cartesian_lim_x,ylim = cartesian_lim_y)+
                    theme(axis.title.x =element_text(size = font_size))+
                    theme(axis.title.y =element_text(size = font_size))+
                    theme(axis.text = element_text(size=font_size))+
                    theme(axis.line = element_line(colour = "black",size=1), panel.border = element_blank(),panel.background = element_blank())+
                    labs(title="PCA - feature selection") + 
                    theme(plot.title = element_text(hjust = 0.5))+
                    theme(panel.border = element_blank(),panel.background = element_blank())+
                    theme(panel.grid.minor = element_line(colour = "white"))+
                    theme(panel.border = element_blank(),panel.background = element_blank())+
                    theme(panel.grid.minor = element_line(colour = "white"))+
                    theme(legend.position="none")
          return(PCA_PLOT)
}

# 2.5 Features selection and SVM Model tuning
SVM_FeatureSel<-function(Table_class_train,kernel_linear=TRUE,cost= 10, degree = 45){
        vec<-NULL  
        for (i in 2:20){
            print(colnames(Table_class_train[i]))
            vec.temp<-readLines("Include Feture ? Y/N")
            vec.temp<-as.character(vec.temp)
            if (vec.temp=="Y"){
                vec<-c(vec,i)
            }
            rm(vec.temp)  
        }
        Table_class_sel<-Table_class_train[,c(vec,-21)]
        x<-Table_class_sel
        y<-Table_class_train[,21]
        if (kernel_linear){
          model<-svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="linear",type = "C",cost= 10, degree = 45,probability = TRUE)
          }else{
          model<-svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="radial",type = "C",cost= 10, degree = 45,probability = TRUE)
          }
        return(model)
}


# -------------------------------------------------------------------------
#' 3. CNN training
# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
#' 4. Deployment example
# -------------------------------------------------------------------------
#' 4.1 SVM
maskGen_deploySVM<-function(x,y,model,Table_class_train,label_class="Postive",erode_mask=T,opensize=9){
        y.pred<-predict(model,Table_class_train, decision.values = T)
        if (label_class=="Postive"){
                Target_class<-"P"
                noTarget_class<-"N"
        }else{
                Target_class<-"N"
                noTarget_class<-"P"    
        }
        d=attr(y.pred,"decision.values")
        new.y.pred=rep(noTarget_class,length(y.pred))
        NewCutoff=0
        new.y.pred[d<NewCutoff]=Target_class
        d<-round(d,1)
        Table_class_train$pred<-as.array(new.y.pred)
        ir = which(Table_class_train$pred %in% noTarget_class) 
        x = rmObjects(x, ir) 
        nr = which(Ts.mix$pred %in% noTarget_class) 
        y = rmObjects(y, nr) 
        if (length(nr) == length(Table_class_train$pred)){
               stop("no classifed cells were detected")
        }
        if (erode_mask){
                y.mask<-erode(y,makeBrush(opensize, shape= "diamond"))
        y.mask_bin<-thresh(y.mask)
        }else{
                y.mask_bin<- thresh(y) 
        }
        return(y.mask_bin)
}

#' 4.2 CNN
maskGen_deployCNN<-function(x,y,z,intens=20,satck_size=200,input_size=61,prediction_th=0.998,model,erode_mask=T,opensize=17){
        stack<-stackObjects(y,x*intens,ext = c(satck_size,satck_size))
        px_st = stack
        n = numberOfFrames(px_st)
        px = resize(px_st, input_size, input_size)
        pxlist = getFrames(px, i=1:n)
        px_A = abind(pxlist, along=0)
        tensor = array_reshape(px_A, c(n, 61, 61, 1))
        preds = predict_on_batch(model, tensor)
        df<-as.data.frame(xy)
        df$Prediction<-as.vector(preds)
        ind.na<-which(is.na(df[,3]))
        df[ind.na,3]<-0
        ind<-which(df[,3] < prediction_th)
        if (length(ind) < 1 ){
                stop("no classifed cells were detected")
        }
        y = rmObjects(y, ind,reenumerate = T)
        z = rmObjects(z, ind,reenumerate = T)
        if (length(table(y)) < 2){
                stop("no classifed cells were detected")
        }
        stack.pred<-stackObjects(y,x*intens,ext = c(satck_size,satck_size))
        indd<-which(stack.pred>0)
        stack.pred[indd]<-1
        xy<-computeFeatures.moment(y)[,c('m.cx','m.cy')]
        if (is.null(nrow(xy))){
                stop("no classifed cells were detected")
                
        }
        df.temp<-as.data.frame(xy)
        df.temp$index<-c(1:dim(stack.pred)[3])
        df.temp$singel_object
        for (m in 1:dim(stack.pred)[3]){
                st.temp<-stack.pred[,,m]
                bg<-length(which(st.temp==0))
                fg<-length(which(st.temp > 0))
                if ((bg/fg) > 10){
                        df.temp[m,4]<-"yes"
                } else {
                        df.temp[m,4]<-"no"
                }
        } 
        rm(ind)
        ind<-which(df.temp[,4]=="no")
        if (length(ind) < 1 ){
                stop("no classifed cells were detected")
                
        }
        y = rmObjects(y, ind,reenumerate = T)
        z = rmObjects(z, ind,reenumerate = T)
        z<-bwlabel(z)
        if (erode_mask){
                z.mask<-erode(z,makeBrush(opensize, shape= "diamond"))
                z.mask_bin<-thresh(z.mask)
        }else{
                z.mask_bin<- thresh(z) 
        }
        return(z.mask_bin)
}


# -------------------------------------------------------------------------
#' 5.confusion matrix and ROC curve 
# -------------------------------------------------------------------------









