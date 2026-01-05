est_sd<-function(df_data_copynumber, potential_pos){
  ratio_list<-c();
  bf_list<-c();
  adjust_ratios<-df_data_copynumber$adjusted.ratio;
  bfs<-df_data_copynumber$Bf;
  sapply( slide( adjust_ratios,
                function(y){sd(y, na.rm=TRUE) }, .after= ceiling(length(adjust_ratios)/20) ), "[", 1)->ratio_list;
  sapply( slide(bfs,
                function(y){sd(y, na.rm=TRUE) }, .after= ceiling(length(bfs)/20) ), "[", 1)->bf_list;
  ratio_list<-sort(ratio_list);
  bf_list<-sort(bf_list);
  t_ratio_sd<-( quantile(ratio_list, 0.25, na.rm=TRUE) );
  t_bf_sd<-( quantile(bf_list, 0.25, na.rm=TRUE) );
  c(t_ratio_sd, t_bf_sd, 0);
}
winso_poten<-function(t_df_data){
  df_data_copynumber<-as.data.frame(t_df_data[,c("chromosome", "position",
                                                 "adjusted.ratio", "Bf", "zygosity.normal")],
                                    stringsAsFactors =FALSE);
  df_data_copynumber[,1:3]<-copynumber::winsorize(df_data_copynumber[,1:3]);
  df_data_copynumber[,3]<-pmin(df_data_copynumber[,3], 2.5);
  df_data_copynumber<-df_data_copynumber[order(df_data_copynumber[,2]),];
  potential_pos_pcr1<-copynumber::pcf(df_data_copynumber[,c(1,2,3)], gamma=15)
  potential_pos<-c(which(df_data_copynumber$position %in% potential_pos_pcr1$start.pos),
                   which(df_data_copynumber$position %in% potential_pos_pcr1$end.pos) );
  potential_pos_pcr2<-copynumber::pcf(df_data_copynumber[,c(1,2,4)], gamma=25)
  potential_pos<-c(potential_pos,which(df_data_copynumber$position %in% potential_pos_pcr2$start.pos),
                   which(df_data_copynumber$position %in% potential_pos_pcr2$end.pos) );
  potential_pos<-c(potential_pos,nrow(df_data_copynumber) );
  potential_pos<-c(1,potential_pos);
  potential_pos<- setdiff( potential_pos, potential_pos+1 );
  potential_pos<-sort(unique(potential_pos) );
  list(df_data_copynumber=df_data_copynumber, potential_pos=potential_pos)
}
