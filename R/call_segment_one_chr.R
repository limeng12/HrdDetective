get_best_seg<-function(t_df_data,
                       potential_pos,
                       tt_cellularity=0.8, tt_ploidy=3,
                       avg_depth_ratio=1,
                       t_ratio_sd=t_ratio_sd,
                       t_bf_sd=t_bf_sd,
                       t_bf_sd_logit=t_bf_sd_logit,
                       t_CNt.max=t_CNt.max){
  n<-nrow(t_df_data);
  beta<-2*log(n)*( (tt_ploidy*( 2.2/(1+exp(-1*10*(tt_cellularity-0.95) ) )+1) ) );
  t_ratio<-t_df_data$adjusted.ratio;
  t_bf<-t_df_data$Bf;
  t_seg_pos<-t_df_data$position;
  zygosity<-t_df_data$zygosity.normal;
  model.pts <- baf.model.points_h(cellularity = tt_cellularity,
                                  ploidy = tt_ploidy,
                                  avg.depth.ratio = avg_depth_ratio,
                                  sd.Bf=t_bf_sd,
                                  sd.ratio=t_ratio_sd,
                                  t_CNt.max=t_CNt.max);
  best_seg<-baf_bayes_fit_cpp(t_ratio, t_bf, t_seg_pos, zygosity,
                              model.pts$CNt, model.pts$B,
                              model.pts$BAF,
                              model.pts$depth.ratio,
                              model.pts$trunca_ratio_t,
                              model.pts$trunca_baf_t,
                              potential_pos,
                              beta,
                              t_ratio_sd,
                              t_bf_sd,
                              t_bf_sd_logit);
  cp<-best_seg$cp;
  fv<-best_seg$fv;
  cp_final<-c(length(cp));
  index<-length(cp);
  while(index!=1){
    cp_final<-c(cp_final,cp[index]);
    index<-cp[index];
  }
  the_dis_threshold<-0.03*max(t_seg_pos);
  if(the_dis_threshold<5*1e6){
    the_dis_threshold<-5*1e6
  }
  cp_final<-c(cp_final, which(diff(t_seg_pos)>the_dis_threshold), which(diff(t_seg_pos)>the_dis_threshold)+1);
  one_fit_v<-fv[length(fv)]
  list(one_fit_v=one_fit_v,
       cp_final=sort(unique(cp_final), decreasing = FALSE),
       t_seg_pos=t_seg_pos[cp_final],
       df_data=t_df_data,
       model.pts=model.pts,
       potential_pos=potential_pos
       ) ;
}
