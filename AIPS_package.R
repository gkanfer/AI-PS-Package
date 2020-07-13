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
load_image<-function(Image_File_name,path,minmax_nrom=TRUE,...){
          x<-readImage(paste0(path,Image_File_name))
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
          channel_set<-list()
          channel_set[["CH1"]]<-CH1
          channel_set[["CH2"]]<-CH2
          return(channel_set)
          }


#' 1.1 Nucleus segmentation 
Nucleus_segmentation<-function(x,intens=2,filter_size=16,offset=0.04,opensize=3,
                               rmObjects_small=30,Use_watershed=TRUE,distmap_value=2,rm_outliers=TRUE,Out_P=0.95,
                               Masks_display=TRUE){
          CH1_seg<- x*intens
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
          mask[["mask_nuc"]]<-gsegg
          mask[["parameters_nuc"]]<-xy
          if (Masks_display){
                seg_CH1mask2 = paintObjects(CH1mask2,toRGB(x),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)  
                seg_mask<-paintObjects(gsegg,toRGB(x),opac=c(1, 1),col=c("red",NA),thick=TRUE,closed=TRUE)
                mask[["mask_start"]]<-seg_CH1mask2
                mask[["mask_final"]]<-seg_mask
          }
          return(mask)
          
}

# 1.2 GFP segmentation 
Cytosol_segmentation <- function (x,y,intens=40,filter_size=10,offset=0.1,mat_size_smooth=19,
                                  opensize=7,rmObjects_large=30000,Use_watershed=TRUE){
          Cyto<- x*intens
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
          colorMode(cseg)<-Grayscale
          cseg = propagate(Cyto, y, lambda=1.0e-2, mask=combine)
          cseg <- fillHull(cseg)
          colorMode(cseg)<-Grayscale
          xy<-computeFeatures.moment(cseg)[,c('m.cx','m.cy')]
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
          cseg = rmObjects(cseg, ci,reenumerate = F) 
          xy.gsegg<-as.numeric(row.names(computeFeatures.moment(y)[,c('m.cx','m.cy')]))
          xy.cseg<- as.numeric(row.names(computeFeatures.moment(cseg)[,c('m.cx','m.cy')])) 
          ind.deff<-setdiff(xy.gsegg,xy.cseg)
          gsegg<-rmObjects(y,ind.deff,reenumerate=F)
          gsegg<-reenumerate(gsegg)
          cseg<-reenumerate(cseg)
          xy.gsegg_table<-computeFeatures.moment(gsegg)[,c('m.cx','m.cy')]
          xy.cseg_table<- computeFeatures.moment(cseg)[,c('m.cx','m.cy')]
          masks<-list()
          masks[["mask_nuc"]]<-gsegg
          masks[["parameters_nuc"]]<-xy.gsegg_table
          masks[["mask_cyto"]]<-cseg
          masks[["parameters_cyto"]]<-xy.cseg_table
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
        table_test <- as.data.frame(cbind(table_basic,table_moment,table_shape))
        rownameTable<-row.names(table_test)
        table_test_Feture<-data.frame(cbind(rownameTable,table_test))
        Ts.mix<-table_test_Feture[,2:20]
        rowNameTable<-table_test_Feture[,1]
        Ts.mix$predict<-label_class
        Features_Table<-list()
        Features_Table[["Ts.mix"]]<-Ts.mix
        Features_Table[["rowNameTable"]]<-rowNameTable
        return(Features_Table)
}

# 2.2 Cell picking 
Pick_cells<-function(mask_nuc,x,Ts.mix,intens,parameters_nuc,font_size=0.7,label_class="Postive",
                     Masks_display=TRUE){
          seg_disp = paintObjects(mask_nuc,toRGB(x*intens),opac=c(1,0.8),col=c("Green",NA),thick=TRUE,closed=FALSE)
          display(seg_disp,"raster")
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
              Ts.mix[row_numb$V1,20]<-"P"
              Ts.mix[-row_numb$V1,20]<-"N"
          }else{
              Ts.mix[row_numb$V1,20]<-"N"  
              Ts.mix[-row_numb$V1,20]<-"P"
              }
          nr = which(Ts.mix$predict %in% "N")
          
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
        ind0<-which(is.na(Table_class_train$predict))
        if (length(ind0)>1) Table_class_train<-Table_class_train[-ind0,]
        ind<-which(is.na(Table_class_train$predict))
        if (length(ind)>1) Table_class_train<-Table_class_train[-ind,]
        ind1<-grep("\\d",Table_class_train$predict)
        if (length(ind)>1) Table_class_train<-Table_class_train[-ind1,]  
        x<-(Table_class_train[,2:19])
        y<-Table_class_train[,20]
          if (kernel_linear){
              acc<-rep(0,100)
               pb <- progress_bar$new(total = 100)
              for (i in 1:length(acc)){
                pb$tick()
                TestIndex<-sample(nrow(x),round(nrow(x)/2))
                model<-svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="linear",type = "C",cost= 10, degree = 45,probability = TRUE)
                y.pred<-predict(model,x[-TestIndex,])
                acc[i]=length(which(y.pred==y[-TestIndex]))/length(y.pred)
                }
              }else{
              acc<-rep(0,100)
              pb <- progress_bar$new(total = 100)
              for (i in 1:length(acc)){
                pb$tick()
                TestIndex<-sample(nrow(x),round(nrow(x)/2))
                model<-svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="radial",type = "C",cost= 10, degree = 45,probability = TRUE)
                y.pred<-predict(model,x[-TestIndex,])
                acc[i]=length(which(y.pred==y[-TestIndex]))/length(y.pred)}
                }
          mean(acc)
          return(acc)
}
# 2.4 Features reduction using PCA analysis

Fetures_PCA<-function(Table_class_train,cartesian_lim_x=c(-12,10),cartesian_lim_y=c(-10,5),font_size=14){
        ind0<-which(is.na(Table_class_train$predict))
        if (length(ind0)>1) Table_class_train<-Table_class_train[-ind0,]
        ind<-which(is.na(Table_class_train$predict))
        if (length(ind)>1) Table_class_train<-Table_class_train[-ind,]
        ind1<-grep("\\d",Table_class_train$predict)
        if (length(ind)>1) Table_class_train<-Table_class_train[-ind1,]   
        Table_class_train<-Table_class_train %>% filter_at(vars(1:20), all_vars(!is.infinite(.)))
        x<-1
        z<-sample(1:nrow(Table_class_train),nrow(Table_class_train))
        while (x<10){
                myPr <- try(prcomp(Table_class_train[z, 2:19], scale = TRUE))
                if (inherits(myPr, "try-error")) {
                        x<-x+1
                        z<-sample(1:nrow(Table_class_train),round(nrow(Table_class_train)/x,2))
                }else{
                        break
                }
        }
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
        ind0<-which(is.na(Table_class_train$predict))
        if (length(ind0)>1) Table_class_train<-Table_class_train[-ind0,]
        ind<-which(is.na(Table_class_train$predict))
        if (length(ind)>1) Table_class_train<-Table_class_train[-ind,]
        ind1<-grep("\\d",Table_class_train$predict)
        if (length(ind)>1) Table_class_train<-Table_class_train[-ind1,]   
        Table_class_train<-Table_class_train %>% filter_at(vars(1:20), all_vars(!is.infinite(.)))
        vec<-NULL  
        for (i in 2:19){
                print(colnames(Table_class_train[i]))
                vec.temp<-readline("Include Feture ? Y/N")
                vec.temp<-as.character(vec.temp)
                if (vec.temp=="Y"){
                vec<-c(vec,i)
                }
                rm(vec.temp)  
        }
        Table_class_sel<-Table_class_train[,c(vec,which(colnames(Table_class_train)=="predict"))]
        x<-Table_class_sel[,-(ncol(Table_class_sel))]
        y<-Table_class_train$predict
        pb <- progress_bar$new(total = 100)
        if (kernel_linear){
                acc<-rep(0,100)
                for (i in 1:length(acc)){
                        pb$tick()
                        TestIndex<-sample(nrow(x),round(nrow(x)/2))
                        model<-svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="linear",type = "C",cost= 10, degree = 45,probability = TRUE)
                        y.pred<-predict(model,x[-TestIndex,])
                        acc[i]=length(which(y.pred==y[-TestIndex]))/length(y.pred)
                }
        }else{
                acc<-rep(0,100)
                for (i in 1:length(acc)){
                        pb$tick()
                        TestIndex<-sample(nrow(x),round(nrow(x)/2))
                        model<-svm(x=x[TestIndex,], y=as.factor(y[TestIndex]),kernel="radial",type = "C",cost= 10, degree = 45,probability = TRUE)
                        y.pred<-predict(model,x[-TestIndex,])
                        acc[i]=length(which(y.pred==y[-TestIndex]))/length(y.pred)}
        }
        model_svm<-list()
        model_svm[["ACC"]]<-acc
        model_svm[["Selected_Features"]]<-colnames(x)
        model_svm[["model"]]<-model
        return(model_svm)
}


# -------------------------------------------------------------------------
# 3.CNN training
# -------------------------------------------------------------------------
# 3.1 single cell images generation and upload
Image_set_gen<-function(y,x,base_dir,train_dir,validation_dir,test_dir,train_norm_dir,train_swol_dir,validation_norm_dir,
                        original_dataset_dir_norm,original_dataset_dir_swol,tr.mx,val.mx,tst.mx){
        #Break the 20x magnification image into single cell images
        stack<-stackObjects(y,x*intens,ext = c(satck_size,satck_size))
        #Folders creation 
        train_dir <- file.path(base_dir, "train")
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


# 3.2 CNN build model
CNN_build_model<-function(batch.size,img,lr,train_dir,validation_dir,epoch, img_size, checkpoint_dir){
        model <- keras_model_sequential() %>% 
                layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu", 
                              input_shape = c(img, img, 1), padding="same") %>%
                layer_conv_2d(filters = 64, kernel_size = c(3, 3), activation = "relu") %>%
                layer_batch_normalization(momentum = 0.9) %>%
                layer_max_pooling_2d(pool_size = c(2, 2)) %>% 
                layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>% 
                layer_conv_2d(filters = 128, kernel_size = c(3, 3), activation = "relu") %>% 
                layer_batch_normalization(momentum = 0.9) %>%
                layer_max_pooling_2d(pool_size = c(2, 2)) %>% 
                layer_conv_2d(filters = 256, kernel_size = c(3, 3), activation = "relu") %>% 
                layer_conv_2d(filters = 256, kernel_size = c(3, 3), activation = "relu") %>% 
                layer_conv_2d(filters = 256, kernel_size = c(3, 3), activation = "relu") %>%
                layer_batch_normalization(momentum = 0.9) %>%
                layer_max_pooling_2d(pool_size = c(2, 2)) %>% 
                layer_conv_2d(filters = 512, kernel_size = c(3, 3), activation = "relu") %>% 
                layer_conv_2d(filters = 512, kernel_size = c(3, 3), activation = "relu") %>% 
                layer_conv_2d(filters = 512, kernel_size = c(3, 3), activation = "relu") %>% 
                layer_conv_2d(filters = 512, kernel_size = c(3, 3), activation = "relu") %>% 
                layer_batch_normalization(momentum = 0.9) %>%
                layer_max_pooling_2d(pool_size = c(2, 2)) %>% 
                layer_flatten() %>% 
                layer_dense(units = 2048, activation = "relu") %>%
                layer_dropout(0.2) %>%
                layer_dense(units = 1, activation = "sigmoid")
        
        model %>% compile(
                loss = "binary_crossentropy",
                optimizer = optimizer_adam(lr),
                metrics = c("acc")
        )
        
        train_datagen <- image_data_generator(
                rescale = 1/255,
                rotation_range = 90, 
                fill_mode = "constant",
                horizontal_flip = TRUE, 
                zoom_range = 0.25)
        
        validation_datagen <- image_data_generator(rescale = 1/255)
        
        train_generator <- flow_images_from_directory(
                train_dir,
                train_datagen,
                color_mode = "grayscale",
                target_size = c(img_size, img_size),
                batch_size = batch.size,
                class_mode = "binary"
        )
        validation_generator <- flow_images_from_directory(
                validation_dir,
                validation_datagen,
                color_mode = "grayscale",
                target_size = c(img_size, img_size),
                batch_size = batch.size,
                class_mode = "binary"
        )
        
        train.n = train_generator$n
        train_step = (train.n%/%batch.size)
        val.n = validation_generator$n
        val_step = (val.n%/%batch.size)
        
        dir.create(checkpoint_dir, showWarnings = FALSE)
        filepath <- file.path(checkpoint_dir, "weights.{epoch:02d}-{val_loss:.2f}.hdf5")
        
        # Create checkpoint callback
        cp_callback <- callback_model_checkpoint(
                filepath = filepath,
                save_weights_only = FALSE,
                period = 3,
                verbose = 1
        )
        
        history <- model %>% fit_generator(
                train_generator,
                steps_per_epoch =train_step,
                epochs = epoch,
                validation_data = validation_generator,
                validation_steps = val_step,
                callbacks = list(cp_callback), 
                verbose = 1, 
                workers = 10
        )
        
        ##Add save history plot (as png) or just let them save from Rstudio plot output during training?
        
        model_h5<-save_model_hdf5(model, "mito_model_4.h5")    
        ##Do we need this if they are saving the model at every X checkpoint? If we are keepint the checkpoint funcitonality I also think we need to include a test for choosing the best model--see section 3.3 
        return(model_h5) ##If we add the test for best model I think this function would return nothing? Since the checkpoint files are being saved in the function.. there's nothing to return ?
     
}

#' 3.3 Choose best model from checkpoints 
best_model <- function(test_dir, img_size, model_dir) {
        test_datagen <- image_data_generator(rescale = 1/255)
        test_generator <- flow_images_from_directory(
                test_dir,
                test_datagen,
                target_size = c(img_size, img_size), 
                color_mode = "grayscale",
                batch_size = 1, 
                class_mode = NULL, 
                shuffle=FALSE
        )
        setwd(model_dir)
        fnames <- dir()
        dfL <- list()
        modelfiles<- function(x) {
                model <- load_model_hdf5(x)
                pred = predict_generator(model, test_generator, steps=(n))
                cm <- CM(pred)
                acc <- cm$overall[[1]]
                nam <- paste(fnames[i])
                df <- data.frame(nam, acc)
                return(df)
        }
        for (i in 1:length(fnames)) {
                model <- load_model_hdf5(fnames[i], custom_objects = c("sensitivity"=custom, "specificity"= custom2))
                pred = predict_generator(model, test_generator, steps=(n))
                cm <- CM(pred)
                acc <- cm$overall[[1]]
                nam <- paste(fnames[i])
                df <- data.frame(nam, acc)
                dfL[[i]] <- df
                print(i)
        }
        
}


# -------------------------------------------------------------------------
#' 4. Deployment example
# -------------------------------------------------------------------------
#' 4.1 SVM
maskGen_deploySVM<-function(x,y,model,Table_class_train,label_class="Postive",erode_mask=T,opensize=9,Sel=TRUE,Selected_Features,TH){
        ind<- which(colnames(Table_class_train) %in% "predict")
        if (length(ind) > 0) {
        Table_class_train<-Table_class_train[,-ind]
        }
        if (Sel){
        Table_class_train<-Table_class_train[,which(colnames(Table_class_train) %in% Selected_Features)]
        }
        y.pred<-predict(model,Table_class_train, decision.values = T)
        if (label_class=="Postive"){
                Target_class<-"P"
                noTarget_class<-"N"
        }else{
                Target_class<-"N"
                noTarget_class<-"P"    
        }
        d=attr(y.pred,"decision.values")[1:nrow(Table_class_train),1]
        new.y.pred=rep(noTarget_class,length(y.pred))
        NewCutoff=TH
        new.y.pred[d>NewCutoff]=Target_class
        d<-round(d,1)
        table_pred<-Table_class_train
        table_pred$pred<-as.array(new.y.pred)
        ir = which(table_pred$pred %in% noTarget_class) 
        x = rmObjects(x, ir) 
        nr = which(table_pred$pred %in% noTarget_class) 
        y = rmObjects(y, nr) 
        if (length(nr) == length(table_pred$pred)){
               stop("no classifed cells were detected")
        }
        if (erode_mask){
                y.mask<-erode(y,makeBrush(opensize, shape= "diamond"))
        y.mask_bin<-thresh(y.mask)
        }else{
                y.mask_bin<- thresh(y) 
        }
        call_mask<-list()
        call_mask[["decision_values"]]<-d
        call_mask[["table"]]<-table_pred
        call_mask[["mask"]]<-y.mask_bin
        return(call_mask)
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

#' 5.1 Test data prep 
test_data <- function(test_dir, img_size) {
        test_datagen <- image_data_generator(rescale = 1/255)
        test_generator <- flow_images_from_directory(
                test_dir,
                test_datagen,
                target_size = c(img_size, img_size), 
                color_mode = "grayscale",
                batch_size = 1, 
                class_mode = NULL, 
                shuffle = FALSE
        )
        return(test_generator)
}

#' 5.2 Load model, Make predictions
test_predict <- function(model_file) {
        
}

#' 5.2 CM

CM <- function(x, ) {
        
        
        
        df <- as.tibble(cbind(pred, test_generator$filenames)) %>%
                rename(
                        predict_proba = V1,
                        filename = V2
                ) %>%
                mutate(predict_proba = as.numeric(predict_proba)) %>%
                mutate(predicted_label = ifelse(predict_proba > 0.5, 1, 0)) %>%
                mutate(predicted_label = as.integer(predicted_label)) %>%
                mutate(predicted_label_name = ifelse(predicted_label == 0, "Norm", "Swol")) %>%
                separate(filename, into=c("true_label","fname"), sep = "[//]" )
        cm <- confusionMatrix(as.factor(df$predicted_label_name), as.factor(df$true_label), positive = "Swol")
        return(cm)
}



library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(caret)


CM <- function(x){
        df <- as.tibble(cbind(pred, test_generator$filenames)) %>%
                rename(
                        predict_proba = V1,
                        filename = V2
                ) %>%
                mutate(predict_proba = as.numeric(predict_proba)) %>%
                mutate(predicted_label = ifelse(predict_proba > 0.5, 1, 0)) %>%
                mutate(predicted_label = as.integer(predicted_label)) %>%
                mutate(predicted_label_name = ifelse(predicted_label == 0, "Norm", "Swol")) %>%
                separate(filename, into=c("true_label","fname"), sep = "[//]" )
        cm <- confusionMatrix(as.factor(df$predicted_label_name), as.factor(df$true_label), positive = "Swol")
        return(cm)
}

fnames <- dir()
dfL <- list()
for (i in 37:length(fnames)) {
        model <- load_model_hdf5(fnames[i], custom_objects = c("sensitivity"=custom, "specificity"= custom2))
        pred = predict_generator(model, test_generator, steps=(n))
        cm <- CM(pred)
        acc <- cm$overall[[1]]
        nam <- paste(fnames[i])
        df <- data.frame(nam, acc)
        dfL[[i]] <- df
        print(i)
}


#' 5.4 ROC 
library("e1071")
library("dplyr")
library("data.table")
library("gtools")
library("yaImpute")
library("tidyr")
library(stringr)
library(caret)
library(ROCR)
pred = prediction(df$predict_proba, df$true_label)
perf = performance(pred, measure="tpr", x.measure="fpr")
plot(perf, colorize=TRUE, add=FALSE)
perf2 = performance(pred, measure="acc", x.measure="cutoff")
plot(perf2, add=FALSE)
perf3 = performance(pred, measure="prec", x.measure="rec")
plot(perf3, add=FALSE, colorize=TRUE)
perf4 <- performance(pred, "acc")
plot(perf4, avg= "vertical", spread.estimate="boxplot", show.spread.at= seq(0.1, 0.9, by=0.1))
perf5 <- performance(pred, "cost")
plot(perf5)
perf5 <- performance(pred, "ecost")
plot(perf5)
threshold1 <- function(predict, response) {
        perf <- performance(pred, "sens", "spec")
        df <- data.frame(cut = perf@alpha.values[[1]], sens = perf@x.values[[1]], spec = perf@y.values[[1]])
        df[which.max(df$sens + df$spec), "cut"]
}
thr = threshold1(pred)
print(thr)

plot(0,0,type="n", xlim= c(0,1), ylim= c(0,500), xlab="Prediction", ylab="Density")
for (i in 1:length(pred@predictions)) { 
        lines(density(subset(pred@predictions[[i]], pred@labels[[i]] == "0")), col= "red")
        lines(density(subset(pred@predictions[[i]], pred@labels[[i]] == "1")), col="green") }






