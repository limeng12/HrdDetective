get_best_seg_multi_chr_simple<-function(process_data_list,tt_cellularity=0.8, tt_ploidy=6,
                                        avg_depth_ratio=avg_depth_ratio,
                                        t_CNt.max=t_CNt.max){
  fit_value_total<-0;
  result_list<-list();
  for(chr in names(process_data_list)){
    print(paste0("chromosome: ", chr) );
    ratio_sd_v<-process_data_list[[chr]]$sd[1]
    bf_sd_v<-process_data_list[[chr]]$sd[2]
    t_bf_sd_logit<-process_data_list[[chr]]$sd[3]
    one_fit<-get_best_seg(process_data_list[[chr]]$df_data_copynumber,
                          process_data_list[[chr]]$potential_pos,
                          tt_cellularity=tt_cellularity, tt_ploidy=tt_ploidy,
                          avg_depth_ratio=avg_depth_ratio,
                          t_ratio_sd=ratio_sd_v,
                          t_bf_sd=bf_sd_v,
                          t_bf_sd_logit=t_bf_sd_logit,
                          t_CNt.max=t_CNt.max
                          );
    fit_value_total<-fit_value_total+one_fit$one_fit_v;
  }
  return(list(fit_value_total=fit_value_total,
  tt_cellularity=tt_cellularity,
  tt_ploidy=tt_ploidy) );
}
get_best_seg_multi_chr<-function(process_data_list,tt_cellularity=0.8, tt_ploidy=6 ,t_CNt.max=t_CNt.max){
  ratio_sd_v<-process_data_list[[1]]$sd[1]
  bf_sd_v<-process_data_list[[1]]$sd[2]
  t_bf_sd_logit<-process_data_list[[1]]$sd[3]
  fit_value_total<-0;
  result_list<-list();
  for(chr in names(process_data_list) ){
    print(paste0("chromosome: ", chr) );
    one_fit<-get_best_seg(process_data_list[[chr]]$df_data_copynumber,
                          process_data_list[[chr]]$potential_pos,
                          tt_cellularity=tt_cellularity, tt_ploidy=tt_ploidy,
                          t_ratio_sd=ratio_sd_v,t_bf_sd=bf_sd_v ,
                          t_bf_sd_logit=t_bf_sd_logit,
                          t_CNt.max=t_CNt.max
    );
    tt_df_data<-one_fit$df_data
    
    model.pts<-one_fit$model.pts;
    fit_value_total<-fit_value_total+one_fit$one_fit_v;
    one_fit$cp_final[length(one_fit$cp_final)]<-one_fit$cp_final[length(one_fit$cp_final)];
    one_fit$cp_final[1]<-0;
    start_pos<-c();end_pos<-c();seg_pos<-c();seg_id<-c();Cnt<-c();A<-c();B<-c();LPP<-c();
    if(length(one_fit$cp_final)>=2){
      for(i in 2:(length(one_fit$cp_final )) ){
        cnt<-baf_bayes_one_fit_cpp(tt_df_data$adjusted.ratio[(one_fit$cp_final[i-1]+1):(one_fit$cp_final[i])],
                                   tt_df_data$Bf[(one_fit$cp_final[i-1]+1):(one_fit$cp_final[i])],
                                   tt_df_data$zygosity.normal[(one_fit$cp_final[i-1]+1):(one_fit$cp_final[i])],
                                   model.pts$CNt, model.pts$B,
                                   model.pts$depth.ratio,model.pts$BAF,
                                   model.pts$trunca_ratio_t,model.pts$trunca_baf_t,
                                   ratio_sd=ratio_sd_v,
                                   bf_sd=bf_sd_v,
                                   bf_sd_logit=t_bf_sd_logit
                                   );
        if(length( (one_fit$cp_final[i-1]+1):(one_fit$cp_final[i]) )<3){
          next;
        }
        start_pos<-c(start_pos, tt_df_data$position[one_fit$cp_final[i-1]+1 ] );
        end_pos<-c(end_pos, tt_df_data$position[(one_fit$cp_final[i]) ] );
        Cnt<-c(Cnt, cnt[1]);
        A<-c(A, cnt[2]);
        B<-c(B, cnt[3]);
        LPP<-c(LPP, cnt[4]);
      }
    }
    result_list[[chr]]<-data.frame(
                                   start.pos=start_pos,
                                   end.pos=end_pos,
                                   chromosome=rep(chr,length(start_pos) ),
                                   CNt=Cnt,
                                   A=A,
                                   B=B,
                                   LPP=LPP,
                                   stringsAsFactors = FALSE);
  }
  return(list(fit_value_total=fit_value_total,
              result_segments=result_list) );
}
