########CKD-pQTL孟德尔随机化分析######
library(data.table)
library(TwoSampleMR)
exposure_dat<- read.csv('CKD_exposure_dat.csv',check.names = F)
table(exposure_dat$id.exposure)
exposure_dat=exposure_dat[exposure_dat$id.exposure=='CKD_overall_EA_JW_20180223_nstud23',]

vector<- c('11543_84_LIMA1_LIMA1',
           '4563_61_PLCG1_PLCG1',
           '2681_23_HGF_HGF',
           '6580_29_PZP_Pregnancy_zone_protein',
           '19154_41_SERPINE2_Protease_nexin_I',
           '6366_38_TXNDC15_TXD15'
)

pp = c('LIMA1','PLCG1','HGF',
       'PZP','SERPINE2','TXNDC15')

for (i in 1:6) {
data = fread(paste0('7.CKD-pqtl/',vector[i],'.txt.gz'),data.table = F)  
outcome_dat=data[data$rsids%in%exposure_dat$SNP,]  
colnames(outcome_dat) <- c("chr.outcome","pos.outcome",'name','SNP',
                           "effect_allele.outcome",
                           "other_allele.outcome", 
                          "beta.outcome",         
                           "pval.outcome",'minus_log10_pval',
                          "se.outcome",           
                           "samplesize.outcome", "eaf.outcome")

outcome_dat$id.outcome = vector[i]
outcome_dat$outcome = pp[i]
harm_rt<-harmonise_data(exposure_dat=exposure_dat,
                        outcome_dat=outcome_dat)
write.csv(harm_rt,file = paste0('7.CKD-pqtl/CKD_', pp[i],'_harm_rt.csv'),row.names = F)
}

df_MR_pleio <- list()
df_MR_hetero <- list()
resall<- list()
harm_rt_a <- list()
for (i in 1:6) {
  
  harm_rt1 = read.csv(paste0('7.CKD-pqtl/CKD_',pp[i],'_harm_rt.csv'),
                  check.names = F)
  harm_rt_a=rbind(harm_rt_a,harm_rt1)
  res_hete <- mr_heterogeneity(harm_rt1)
  if (res_hete$Q_pval[2]<0.05) {
  res=mr(harm_rt1, method_list = c("mr_egger_regression",
                                  "mr_weighted_median", 
                                  "mr_ivw_mre", 
                                  "mr_simple_mode", 
                                  "mr_weighted_mode"))
  
  } else{
  res=mr(harm_rt1, method_list = c("mr_egger_regression",
                                  "mr_weighted_median", 
                                  "mr_ivw_fe", 
                                  "mr_simple_mode", 
                                  "mr_weighted_mode"))
  }
  resall = rbind(resall,res)
  id = resall
  id$pval1=round(id$pval,digits=3)
  id$OR=exp(id$b)
  id$or_uci95 = exp(id$b+1.96*id$se)
  id$or_lci95 = exp(id$b-1.96*id$se)
} 
  write.csv(id,'7.CKD-pqtl/CKD-pQTL_mr.csv',row.names = F)
  
  
#######剔除异常SNP########
  library(TwoSampleMR)

  dat2 <- list()
  for (i in c('PZP','SERPINE2')) {

    dat=read.csv(paste0('7.CKD-pqtl/CKD_',i,'_harm_rt.csv'),
                 check.names = F)
    
    res_single <- mr_singlesnp(dat,all_method =c('mr_egger_regression',
                                                 'mr_ivw_mre'))
    res_single=res_single[!duplicated(res_single$SNP),]
    res_single=res_single[res_single$p>0.05,]
    dat = dat[dat$SNP%in%res_single$SNP,]
    dat2 = rbind(dat2,dat)
  }
  
  table(harm_rt_a$outcome)
  
  dat1 = harm_rt_a[harm_rt_a$outcome%in%c('LIMA1','PLCG1','TXNDC15','HGF'),]
  data1 = rbind(dat1,dat2)
  write.csv(data1,'7.CKD-pqtl/harm_rt2.csv',row.names = F)
  ####
###############MR#################
data1 <- read.csv('7.CKD-pqtl/harm_rt2.csv',check.names = F)
  df_MR_pleio <- list()
  df_MR_hetero <- list()
  resall<- list()
  for (i in 1:6) {
    
    harm_rt1 = data1[data1$outcome==pp[i],]
    res_hete <- mr_heterogeneity(harm_rt1)
    if (res_hete$Q_pval[2]<0.05) {
      res=mr(harm_rt1, method_list = c("mr_egger_regression",
                                       "mr_weighted_median", 
                                       "mr_ivw_mre", 
                                       "mr_simple_mode", 
                                       "mr_weighted_mode"))
      
    } else{
      res=mr(harm_rt1, method_list = c("mr_egger_regression",
                                       "mr_weighted_median", 
                                       "mr_ivw_fe", 
                                       "mr_simple_mode", 
                                       "mr_weighted_mode"))
    }
    resall = rbind(resall,res)
    id = resall
    id$pval1=round(id$pval,digits=3)
    id$OR=exp(id$b)
    id$or_uci95 = exp(id$b+1.96*id$se)
    id$or_lci95 = exp(id$b-1.96*id$se)
  } 
  
  write.csv(id,'7.CKD-pqtl/CKD-pQTL_mr2.csv',row.names = F)
###
############7.森林图#######

library(rms)
library(pROC)
library(ormPlot)
library(regplot)
library(dplyr)
library(forestploter)
library(grid)
source("d:/public data/代码2024/survival_model.R")
library(stringr)
library(readxl)


id <- read.csv('7.CKD-pqtl/CKD-pQTL_mr2_0.05.csv')###手动筛选CKD-pQTL_mr2.csv文件所获得

id$pval1=round(id$pval,digits=3)
id$OR=exp(id$b)
id$or_uci95 = exp(id$b+1.96*id$se)
id$or_lci95 = exp(id$b-1.96*id$se)
# id$pval=id$pval1
id$OR=round(id$OR,3)
id$or_uci95 = round(id$or_uci95,3)
id$or_lci95 = round(id$or_lci95,3)

df_tb = id

df_tb$`OR(95% CI)` <- paste0(df_tb$OR,
                             "(", df_tb$or_lci95,"-",
                             df_tb$or_uci95, ")")
# 森林图
dt <- df_tb
dt$` `<- paste(rep(" ", 20), collapse = " ")

forest_dt <- dt[, c("exposure", 'outcome'," ", 
                    'OR(95% CI)','pval',"nsnp")]

p <- forest(forest_dt,
            est = dt$OR,       #效应值
            lower = dt$or_lci95,     #可信区间下限
            # upper = dt$or_uci95,      #可信区间上限
            upper=ifelse(dt$or_uci95<13,dt$or_uci95,20),
            sizes = 0.5,     #黑框的大小
            ci_column = 3,   #在那一列画森林图，要选空的那一列
            ref_line = 1,
            # xlim = c(0, 2),
            # ticks_at = c(0.95, 1, 1.05),
            theme = tm)
plot(p)
g <- edit_plot(p,
               row = which(dt$OR < 1),
               col = 3,
               which = "ci",
               gp = gpar(fill = "#ABDAEC"))
plot(g)

dev.copy2pdf(file = "7.CKD-pqtl/forestplot.pdf", 
             width = 9, height = 4)
dev.off()
###


############7.敏感性分析画图#######

library(TwoSampleMR)
data <- read.csv('7.CKD-pqtl/harm_rt2.csv',check.names = F)
id <- read.csv('7.CKD-pqtl/CKD-pQTL_mr2_0.05.csv')
cpg = id$id.outcome
name = id$outcome
df_MR_pleio <- list()
df_MR_hetero <- list()

for (i in 1:6) {
  # i=1
  dat=data[data$id.outcome==cpg[i],]
  
  res_hete <- mr_heterogeneity(dat)
  # res_pleio <- mr_pleiotropy_test(dat)
  #   df_MR_pleio = rbind(df_MR_pleio,res_pleio)
  #   df_MR_hetero = rbind(df_MR_hetero,res_hete)
  # write.csv(df_MR_pleio,'4.pQTL-CAD/MR_pleio.csv',row.names = F)
  # write.csv(df_MR_hetero,'4.pQTL-CAD/MR_hetero.csv',row.names = F)
  # }
  
  if (res_hete$Q_pval[2]<0.05) {
    
    res=mr(dat, method_list = c("mr_egger_regression",
                                "mr_weighted_median",
                                "mr_ivw_mre",
                                "mr_simple_mode",
                                "mr_weighted_mode"))
  } else{
    
    res <- mr(dat,method_list = c("mr_egger_regression",
                                  "mr_weighted_median", 
                                  "mr_ivw_fe",
                                  "mr_simple_mode",
                                  "mr_weighted_mode"))
    }
    
    
    library(patchwork)
    library(ggplot2)
    library(ggpubr)
    
    single <- mr_leaveoneout(dat)
    single=single[!duplicated(single$SNP),]
    ##rs6441169
    pa1 <- mr_leaveoneout_plot(single)[[1]]
    
    pa1 <- pa1+scale_color_brewer(palette = 'Set1')+
      theme_classic(base_size = 11)+
      theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.position = 'none')+
      update_geom_defaults("line", list(size = 5))+
      labs(x=name[i])
    pa1
    
    # 7. Visualization results
    # 7.1 Scatter Plot
    ###########    IVW这一条线的截距指的是混淆因素的影响，越小越好。
    
    pa2 <- mr_scatter_plot(dat,mr_results = res)[[1]]+
      scale_color_brewer(palette = 'Set1')+
      # theme_classic(base_size = 18)+
      theme(panel.grid=element_blank())+
      bgcolor('white')+
      theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
            axis.ticks = element_line(size = 1))+
      update_geom_defaults("line", list(size = 5))+
      labs(x=name[i])
    pa2$labels$x <- 'CKD'
    pa2$labels$y <- name[i]
    pa2
    # 7.2 Forest Plot
    if (res_hete$Q_pval[2]<0.05) {
    res_single <- mr_singlesnp(dat,all_method =c('mr_egger_regression',
                                                 'mr_ivw_mre'))
    }else{
      res_single <- mr_singlesnp(dat,all_method =c('mr_egger_regression',
                                                   'mr_ivw_fe'))
    }
    res_single=res_single[!duplicated(res_single$SNP),]
    x =  res_single[nrow(res_single),]
    x$SNP='All - IVW'
    res_single = rbind(res_single[1:(nrow(res_single)-1),],x)
    # res_single = res_single[res_single$p>0.05,]
    
    pa3 <- mr_forest_plot(res_single)[[1]]+
      scale_color_brewer(palette = 'Set1')+
      theme_classic(base_size = 11)+
      
      theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
            axis.ticks = element_line(size = 1),
            legend.position = 'none')+
      update_geom_defaults("line", list(size = 5))+
      labs(x=name[i])
    pa3
    # 7.3 Funnel Plot
    
    pa4 <- mr_funnel_plot(res_single)[[1]]+
      scale_color_brewer(palette = 'Set1')+
      # theme_classic(base_size = 18)+
      theme(panel.grid=element_blank())+
      bgcolor('white')+
      theme(panel.border = element_rect(size = 1.7,fill = 'transparent'),
            axis.ticks = element_line(size = 1))+   ##legend.position = 'none'
      update_geom_defaults("line", list(size = 5))+
      labs(x=name[i])
    
    # pa3+pa1+pa2+pa4
    # pa2+pa4
    # ggsave(paste0('7.CKD-pqtl/敏感性检验/',name[i],'散点图-漏斗图.pdf'),
    #        width = 8.5,height=5)
    # print(pa2)
    # ggsave(paste0('7.CKD-pqtl/敏感性检验/',name[i],'散点图.pdf'),
    #        width = 4.5,height=4.5)
    # print(pa4)
    # ggsave(paste0('7.CKD-pqtl/敏感性检验/',name[i],'漏斗图.pdf'),
    #        width = 4,height=4)
    pa3+pa1
    ggsave(paste0('7.CKD-pqtl/敏感性检验/',name[i],'森林图.pdf'),
           width = 8.5,height=length(dat)/10*4)
  }
  
  ##
  
  ######中介效应计算#######
  
  mr1 <- read.csv('1.CKD-CAD/CKD-CAD_mr2.csv',check.names = F)
  mr1 = mr1[mr1$method=='Inverse variance weighted (fixed effects)',]
  mr1$b = round(mr1$b,3)
  mr1$uci95 = round(mr1$b+1.96*mr1$se,3)
  mr1$lci95 = round(mr1$b-1.96*mr1$se,3)
  mr1$`b(95% CI)` <- paste0(mr1$b,
                            "(", mr1$lci95,"-",
                            mr1$uci95, ")")
  
  mr2<- read.csv('7.CKD-pqtl/CKD-pQTL_mr2_0.05.csv',check.names = F)
  mr2$b = round(mr2$b,3)
  mr2$uci95 = round(mr2$b+1.96*mr2$se,3)
  mr2$lci95 = round(mr2$b-1.96*mr2$se,3)
  mr2$`b1(95% CI)` <- paste0(mr2$b,
                             "(", mr2$lci95,"-",
                             mr2$uci95, ")")
  
  mr12 = merge(mr1[,c(1,7,14,15,16)],mr2[,c(1,3,7,14,15,16)],by='id.exposure')
  mr12$Mediator=mr12$outcome
  
  mr3 <- read.csv('4.pQTL-CAD/pQTL_MR_0.05_2.csv',check.names = F)
  mr3$b = round(mr3$b,3)
  mr3$uci95 = round(mr3$b+1.96*mr3$se,3)
  mr3$lci95 = round(mr3$b-1.96*mr3$se,3)
  mr3$`b2(95% CI)` <- paste0(mr3$b,
                             "(", mr3$lci95,"-",
                             mr3$uci95, ")")
  mr3$Mediator=mr3$exposure
  mr123 = merge(mr12[,c(2,3,4,5,7,8,9,10,11)],mr3[,c(7,17,14,15,16)],by='Mediator')
  mr123$b12 = mr123$b.y*mr123$b
  mr123$uci95_12 =mr123$uci95.y*mr123$uci95
  mr123$lci95_12 =mr123$lci95.y*mr123$lci95
  mr123$`b_12(95% CI)` <- paste0(mr123$b12,
                                 "(", mr123$lci95_12,"-",
                                 mr123$uci95_12, ")")
  
  mr123$proportion =round(mr123$b12/mr123$b.x*100,3)
  mr123$proportion_u =round(mr123$uci95_12/mr123$uci95.x*100,3)
  mr123$proportion_l =round(mr123$lci95_12/mr123$lci95.x*100,3)
  mr123$`proportion(95% CI)` <- paste0(mr123$proportion,
                                       "(", mr123$proportion_l,"-",
                                       mr123$proportion_u, ")")
  
  mr123$b_direct = round(mr123$b.x-mr123$b12,3)
  mr123$uci95_direct =round(mr123$uci95.x-mr123$uci95_12,3)
  mr123$lci95_direct =round(mr123$lci95.x-mr123$lci95_12,3)
  mr123$`b_direct(95% CI)` <- paste0(mr123$b_direct,
                                     "(", mr123$lci95_direct,"-",
                                     mr123$uci95_direct, ")")
  
  write.csv(mr123,'7.CKD-pqtl/table5.csv',row.names = F)
  ####
  
  
  