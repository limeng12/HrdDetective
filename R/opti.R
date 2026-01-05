#' Optimized the model, get best cellularity and best ploidy, then segments each chromosome
#'
#' @param data.file the input file, can be got from sequenza-utils
#' @param number_of_cores the number of cores used to optimized the model.
#' @param CNt.max the upper limit of ploidy.
#' @param min.tumor.reads the mininal number of reads in tumor required to support the SNP.
#' @param cellularity the range of the input cellularity.
#' @param ploidy the range of the input ploidy.
#' @return list of below elements.
#' @return \code{seg.tab} list of all segments.
#' @return \code{seqz_list} the input data.
#' @return \code{celluPloidyMa} cellularity and ploidy matrix.
#' @return \code{best_ploidy} optimized ploidy.
#' @return \code{best_cellularity} optimized cellularity.
#' @examples
#' data.file <-  system.file("extdata", "example.seqz", package = "HrdDetective");
#'
#' seg_obj<-opti_and_segs(data.file=data.file);
#'
#' output_plot(seg_obj$seg.tab, seg_obj$seqz_list, seg_obj$celluPloidyMat,
#'             seg_obj$best_ploidy,
#'             seg_obj$best_cellularity,
#'             output_dir_name="./");
#'
#' @export



opti_and_segs<-function(data.file,
                        number_of_cores=6,
                        CNt.max=8,
                        min.tumor.reads=20,
                        cellularity = seq(0.1, 1, by = 0.01*grid_density),
                        ploidy = seq(1, 6, by = 0.1*grid_density) ){
  
  window_size=1e5;
  grid_density=1;
  extract_data<-sequenza.extract_custom(file=data.file,
                                        window=window_size,
                                        min.reads.baf=min.tumor.reads);
  
  seqz_list<-extract_data[[1]];
  seqz_hetero_list<-extract_data[[2]];
  avg_depth_ratio<-extract_data[[3]];
  seqz_list$X<-NULL
  seqz_hetero_list$X<-NULL
  seqz_hetero_list_processed<-list();
  seqz_list_processed<-list();
  for(i in names(seqz_hetero_list) ){
    seqz_hetero_list_processed[[i]]$df_data_copynumber<-winso_poten(seqz_hetero_list[[i]])$df_data_copynumber;
    seqz_hetero_list_processed[[i]]$potential_pos<-(winso_poten(seqz_hetero_list[[i]])$potential_pos);
    seqz_hetero_list_processed[[i]]$sd<-(est_sd(seqz_hetero_list_processed[[i]]$df_data_copynumber,
                                                seqz_hetero_list_processed[[i]]$potential_pos) );
    seqz_list_processed[[i]]$df_data_copynumber<-winso_poten(seqz_list[[i]])$df_data_copynumber;
    seqz_list_processed[[i]]$potential_pos<-winso_poten(seqz_list[[i]])$potential_pos;
    seqz_list_processed[[i]]$sd<-(est_sd(seqz_list_processed[[i]]$df_data_copynumber,
                                         seqz_list_processed[[i]]$potential_pos) );
  }
  print("finished processed")

  # return(seqz_hetero_list_processed);
  celluPloidyMat <- expand.grid(ploidy = ploidy, cellularity = cellularity,
                        KEEP.OUT.ATTRS = FALSE);
  best_cellularity<-0;
  best_ploidy<-0;
  fit_cp_best<-Inf;
  celluPloidyMat$fit_v<-rep(NA, nrow(celluPloidyMat));
  q <- task_q$new();
  q$initialize(number_of_cores);
  for(i in 1:nrow(celluPloidyMat) ){
    print(paste0("cellularity and ploidy: ", celluPloidyMat[i,]$cellularity," ", celluPloidyMat[i,]$ploidy));
    q$push(function(t_process_data_list, tt_cellularity, tt_ploidy, avg_depth_ratio, t_CNt.max){
      library(HrdDetective)
      init_log_table();
      fit_cp_one<-get_best_seg_multi_chr_simple(t_process_data_list,
                                                tt_cellularity=tt_cellularity,
                                                tt_ploidy=tt_ploidy,
                                                avg_depth_ratio=avg_depth_ratio,
                                                t_CNt.max=t_CNt.max);
      
      fit_cp_one
    }, list(t_process_data_list=seqz_hetero_list_processed,
            tt_cellularity=celluPloidyMat[i,]$cellularity,
            tt_ploidy=celluPloidyMat[i,]$ploidy,
            avg_depth_ratio=avg_depth_ratio,
            t_CNt.max=CNt.max) );
  }
  
  while (!q$is_idle()) {
    cat("c");
    task_result <- q$pop(Inf);
    cur_ploidy<-task_result$result$tt_ploidy
    cur_cellurity<-task_result$result$tt_cellularity
    cur_fit_value<-task_result$result$fit_value_total;
    celluPloidyMat$fit_v[intersect(which(celluPloidyMat[,1]==cur_ploidy),
                           which( celluPloidyMat[,2]==cur_cellurity) ) ]<-cur_fit_value
    print(paste0("cur_ploidy:", cur_ploidy," ",
                 "cur_cellurity:",cur_cellurity," ",
                 "cur_fit:",cur_fit_value ) );
  }
  
  best_fit<-min(celluPloidyMat[,"fit_v"]);
  best_ploidy<-celluPloidyMat$ploidy[which.min(celluPloidyMat[,"fit_v"]) ];
  best_cellularity<-celluPloidyMat$cellularity[which.min(celluPloidyMat[,"fit_v"]) ];
  segments_list<-list();
  fit_cp_best_list<-get_best_seg_multi_chr(seqz_list_processed,
                                           tt_cellularity=best_cellularity,
                                           tt_ploidy=best_ploidy,
                                           t_CNt.max=CNt.max);
  for( chr in names(fit_cp_best_list$result_segments)){
    fit_cp_best_list$result_segments$chr$chrom=fit_cp_best_list$result_segments$chr$chromosome
    filtered_one_seg<-segment.breaks (seqz_list[[chr]],
                                      fit_cp_best_list$result_segments[[chr]], min.reads.baf = 1);
    segments_list[[chr]]<-inner_join( filtered_one_seg,
      fit_cp_best_list$result_segments[[chr]][,c("chromosome","start.pos","end.pos","CNt","A","B","LPP")],
        c("chromosome"="chromosome","start.pos"="start.pos","end.pos"="end.pos") );
    segments_list[[chr]]$end.pos<-segments_list[[chr]]$end.pos-1;
  }
  seg.tab <- do.call(rbind, segments_list);
  seg.tab<-seg.tab[seg.tab$N.BAF>=5,];
  return(list(seg.tab=seg.tab,
              seqz_list=seqz_list,
              celluPloidyMat=celluPloidyMat,
              best_ploidy=best_ploidy,
              best_cellularity=best_cellularity,
              seqz_hetero_list_processed=seqz_hetero_list_processed,
              seqz_list_processed=seqz_list_processed
              ) );
}
