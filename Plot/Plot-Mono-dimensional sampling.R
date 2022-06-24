### Script to produce the figures associated to Mono-dimensional Sampling : Figure 2 and Figures S5 to S13 in SI I.


# Set working directory to the script location.  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


install.packages(c("ggridges", "ggplot2", "reshape2", "plyr", "dplyr", "tidyverse","plotly","GGaly","rio","ggthemes","foreach","parallel","doParallel","foreach"))

install.packages(c("GGally","pastecs"))

install.packages(c("rstatix", "ggpubr"))
install.packages("psych")

library(readxl)
library(ggplot2)
library(reshape2) 
library(plyr)
library(dplyr)
library(ggthemes) # Charger

library(plotly)
library(GGally)
library(data.table)

library(remotes)
library(tidyverse)
library(ggpubr)
library(rstatix)


remotes::install_github('rpkgs/gg.layers')
library(gg.layers)
library(ggplot2)


library(psych) 
library(pastecs)


c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)





### Load data

list_skarka=c('ES','SE','IT','PT','UK','FR','EL','CY','IE','DE')



# To keep the columns' names
to_get_columns <- fread("../Outputs-Mono/results_table_df_pointUK_1.csv")


list_result = list()




for (country in list_skarka){
  print(country)
  for(count in 0:8 ){
    
    name_file=paste("../Outputs-Mono/results_table_df_point",country,'_',as.character(count),".csv",sep='')
    
    if(file.exists(name_file)){
      #results_table =read.csv(name_file,sep=";")
      results_table <- fread(name_file)
      
      list_result=append(list_result,list(results_table))
      
      #Total_results_table=rbind(Total_results_table,results_table)
      
    }else{print(c(name_file,"does not exist"))}
    
  }}

# combine list of dataframes in a big dataframe
Total_results_table=data.frame(rbindlist(list_result))


colnames(Total_results_table)=colnames(to_get_columns)

# Delete to save memory

list_result=0
to_get_columns= 0

# Set some columns as factors

Total_results_table$night_monitoring<-as.factor(Total_results_table$night_monitoring)
Total_results_table$CNTR<-as.factor(Total_results_table$CNTR)
Total_results_table$Bio_class<-as.factor(Total_results_table$Bio_class)
Total_results_table$market_for_substitution<-as.factor(Total_results_table$market_for_substitution)
Total_results_table$Nsource<-as.factor(Total_results_table$Nsource)



names(Total_results_table)[names(Total_results_table) == "WDP"] <- "WD, m3 water"
names(Total_results_table)[names(Total_results_table) == "GWP100"] <- "GW100, kg CO2-eq"
names(Total_results_table)[names(Total_results_table) == "FEP"] <- "FE, kg P-eq"
names(Total_results_table)[names(Total_results_table) == "TETPinf"] <- "TETinf, kg 1.4-DC"






# Create quantile for the content in target molecule

quantile_numbers = 3
subtotal <- within(Total_results_table, Quantile_conc <- as.integer(cut(bioact_molec_dbio, quantile(bioact_molec_dbio, probs=0:quantile_numbers/quantile_numbers), include.lowest=TRUE)))

#The column "Quantile_con" is a factor
subtotal$Quantile_conc<-as.factor(subtotal$Quantile_conc)




###
#Chose a few random points from the sample
###


set.seed(1)
number_points = 100000

row_numbers=dim(subtotal)[1]
randomindexes=sample.int(row_numbers,number_points , replace = FALSE)


micro_sample=subtotal[randomindexes,]









#Melt 



# NOT ENOUGH MEMORY TO PROCESS A WHOLE TABLE 
Total_melted=melt(subtotal,measure.vars = c("WD, m3 water","GW100, kg CO2-eq","FE, kg P-eq","TETinf, kg 1.4-DC"),variable.name = 'Impact_Category',value.name='Score')

Total_melted$Impact_Category=as.factor(Total_melted$Impact_Category)

# To label IC with exponents and subscripts, we create a new column 
# identical to IC but with an expression

Total_melted$facets = factor(Total_melted$Impact_Category, labels = c(
  "WD ~(m^{3}~water)", 
  "GW100~(kg~CO[2]-eq)", 
  "FE~(kg~P-eq)",
  "TETinf~(kg~1.4-DC)"))



# Add mean column for plot vercital lines

Total_melted <- Total_melted %>% group_by(facets,night_monitoring,CNTR,Bio_class,Quantile_conc)%>% 
  mutate(mean = mean(Score),mean=mean(Score)) %>% # calculate mean for plotting as well
  ungroup()






# Run regressions
model_TET <- lm(`TETinf, kg 1.4-DC` ~ bioact_molec_dbio, data = micro_sample)
summary(model_TET)

coefs_TET<- coef(lm(`TETinf, kg 1.4-DC` ~ as.numeric(Quantile_conc), data = micro_sample))


model_FE <- lm(`FE, kg P-eq` ~ bioact_molec_dbio, data = micro_sample)
summary(model_FE)

coefs_FE<- coef(lm(`FE, kg P-eq` ~ as.numeric(Quantile_conc), data = micro_sample))


model_WD <- lm(`WD, m3 water` ~ bioact_molec_dbio, data = micro_sample)
summary(model_WD)

coefs_WD<- coef(lm(`WD, m3 water` ~ as.numeric(Quantile_conc), data = micro_sample))




model_GW <- lm(`GW100, kg CO2-eq` ~ bioact_molec_dbio, data = micro_sample)
summary(model_GW)



coefs_GW<- coef(lm(`GW100, kg CO2-eq` ~ as.numeric(Quantile_conc), data = micro_sample))


# Correlations


cor(micro_sample$bioact_molec_dbio, micro_sample$`GW100, kg CO2-eq`,method =  "spearman")
cor(micro_sample$bioact_molec_dbio, micro_sample$`TETinf, kg 1.4-DC`,method =  "spearman")
cor(micro_sample$bioact_molec_dbio, micro_sample$`FE, kg P-eq`,method =  "spearman")
cor(micro_sample$bioact_molec_dbio, micro_sample$`WD, m3 water`,method =  "spearman")







#with regression

coeffs=c(coefs_GW[2], coefs_WD[2], coefs_TET[2],coefs_FE[2])

col=c()
for(e in coeffs){if(e < 0){
  b="Negative"}
  else{b="Positive"}
  col=c(col,b)
  print(b)}

regr <- data.frame(Impact_Category = c("GW100, kg CO2-eq", "WD, m3 water", "TETinf, kg 1.4-DC","FE, kg P-eq"),
                   inter = c(coefs_GW[1], coefs_WD[1], coefs_TET[1],coefs_FE[1]),
                   slo=coeffs,
                   Coefficient=col)


# Boundaries of the quantiles

Quantile0=unname(quantile(Total_melted$bioact_molec_dbio, probs = seq(0, 1, 1/3))[1])
Quantile1=unname(quantile(Total_melted$bioact_molec_dbio, probs = seq(0, 1, 1/3))[2])
Quantile2=unname(quantile(Total_melted$bioact_molec_dbio, probs = seq(0, 1, 1/3))[3])
Quantile3=unname(quantile(Total_melted$bioact_molec_dbio, probs = seq(0, 1, 1/3))[4])

Quantile0
Quantile1
Quantile2
Quantile3




#FIGURE 2 


regr$facets = factor(regr$Impact_Category, 
                     labels = c( "FE~(kg~P-eq)",
                                 "GW100~(kg~CO[2]-eq)",
                                 "TETinf~(kg~1.4-DC)",
                                 "WD ~(m^{3}~water)"
                     ))



Total_melted %>%
  mutate(class = fct_reorder(CNTR, Score, .fun='median')) %>%
  ggplot( aes(x=Quantile_conc, y=Score,fill=reorder(CNTR,Score))) +facet_grid(rows = vars(`facets`),scales="free",labeller = label_parsed)+scale_fill_brewer(palette="Spectral")+
  geom_boxplot2(width = 0.8, width.errorbar = 0.05)+stat_summary(fun=mean, geom="point", shape=21, size=3,alpha=0.7)+
  xlab("Quantile for the content of bioprospected compound in the biomass") + 
  ylab("Score")+guides(fill=guide_legend(title="Country"))+theme(axis.text=element_text(size=12),
                                                                 legend.title = element_text(size=12,face="bold"))+
  
  geom_abline(aes(intercept =inter, slope = slo,color=Coefficient),regr,size=1.5,alpha=0.5)+theme_calc()


ggsave('figure2.jpeg', width = 6, height = 6,dpi=600)













################################
#Supplementary Information
#############################













ggplot(micro_sample[micro_sample$Quantile_conc=="1",], aes(x=bioact_molec_dbio, y=`TETinf, kg 1.4-DC`,color=market_for_substitution),alpha=0.2)+scale_colour_brewer(palette="Spectral")+
  geom_point(size=2,alpha=0.6)+labs(x = bquote("bioact_molec_dbio"~(g.(g[dbio])^-1)), y = bquote(`TETinf, kg 1.4-DC`))
ggsave('FigureS13.jpeg', width = 10, height = 6,dpi=600)





ggplot(micro_sample, aes(x=bioact_molec_dbio, y=log(`GW100, kg CO2-eq`),color=CNTR))+scale_colour_brewer(palette="Spectral")+
  geom_point(size=2,alpha=0.6)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+labs(x = bquote("bioact_molec_dbio"~(g.(g[dbio])^-1)), y = bquote(log("GW100 ,"~kg~CO^2-eq)))

ggsave('FigureS10.jpeg', width = 10, height = 7,dpi=600)


model_GW_log <- lm(log(`GW100, kg CO2-eq`) ~ bioact_molec_dbio, data = micro_sample)
summary(model_GW_log)


# 
# ggplot(micro_sample, aes(x=bioact_molec_dbio, y=log(`WD, m3 water`),color=CNTR))+scale_colour_brewer(palette="Spectral")+
#   geom_point(size=2,alpha=0.6)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+labs(x = bquote("bioact_molec_dbio"~(g.(g[dbio])^-1)), y = bquote(log("WD,"~m^3~water)))
# 
# ggsave('Figure_pairplot_WD_bioact-molec_dbio_ok_log_correc.jpeg', width = 10, height = 7,dpi=600)

# 
# ggplot(micro_sample, aes(x=bioact_molec_dbio, y=log(`FE, kg P-eq`),color=CNTR))+scale_colour_brewer(palette="Spectral")+
#   geom_point(size=2,alpha=0.6)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+labs(x = bquote("bioact_molec_dbio"~(g.(g[dbio])^-1)), y = bquote(log(`FE, kg P-eq`)))
# 
# ggsave('Figure_pairplot_FE_bioact-molec_dbio_ok_log_correc.jpeg', width = 10, height = 7,dpi=600)
# 








#####
#Other pairlpots


# Tplateau # NO
# ggplot(micro_sample, aes(x=T_plateau, y=log(`GW100, kg CO2-eq`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)
# 
# cor(micro_sample$T_plateau, micro_sample$`GW100, kg CO2-eq`,method =  "spearman")
# 
# 
# # Areal productivity ok
# ggplot(micro_sample, aes(x=`Areal productivity kg.m-2.d-1`, y=log(`GW100, kg CO2-eq`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Areal productivity"~(kg.~m^-2~.~d^-1)), y = bquote("log(GW100"~(kg~CO^2-eq)~")"))
# ggsave('Figure_pairplot_GW_arealprod_log.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# 
# ggplot(micro_sample, aes(x=`Areal productivity kg.m-2.d-1`, y=log(`FE, kg P-eq`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Areal productivity"~(kg.~m^-2~.~d^-1)), y = bquote("log(FE, kg P-eq)"))
# ggsave('Figure_pairplot_FE_arealprod_log.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# 
# ggplot(micro_sample, aes(x=`Areal productivity kg.m-2.d-1`, y=`TETinf, kg 1.4-DC`,color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Areal productivity"~(kg.~m^-2~.~d^-1)), y = bquote(`TETinf, kg 1.4-DC`))
# ggsave('Figure_pairplot_TET_arealprod_nolog.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# 
# ggplot(micro_sample, aes(x=`Areal productivity kg.m-2.d-1`, y=log(`TETinf, kg 1.4-DC`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Areal productivity"~(kg.~m^-2~.~d^-1)), y = bquote(log(`TETinf, kg 1.4-DC`)))
# ggsave('Figure_pairplot_TET_arealprod_log.jpeg', width = 10, height = 7,dpi=600)
# 
# ggplot(micro_sample, aes(x=`Areal productivity kg.m-2.d-1`, y=log(`WD, m3 water`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Areal productivity"~(kg.~m^-2~.~d^-1)), y = bquote(log(`WD, m3 water`)))
# ggsave('Figure_pairplot_WD_arealprod_log.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# 
# # Lipid ok
# micro_sample$lipid_af_dwggplot(DF, aes(A, B)) + geom_point() + 
#   ylab(bquote(Vc[max](mu~mol ~CO[2]~ m^-2~s^-1)))
# 
# 
# 
# ggplot(micro_sample, aes(x=lipid_af_dw, y=log(`GW100, kg CO2-eq`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lipid_af_dw"~(g~.~(g["ash-free dw"])^-1)), y =  bquote("log(GW100"~(kg~CO^2-eq)~")"))
# ggsave('Figure_pairplot_GW_lipidafdw_log.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# ggplot(micro_sample, aes(x=lipid_af_dw, y=log(`FE, kg P-eq`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lipid_af_dw"~(g~.~(g["ash-free dw"])^-1)), y = bquote(log(`FE, kg P-eq`)))
# ggsave('Figure_pairplot_FE_lipidafdw_log.jpeg', width = 10, height = 7,dpi=600)
# 
# ggplot(micro_sample, aes(x=lipid_af_dw, y=log(`WD, m3 water`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lipid_af_dw"~(g~.~(g["ash-free dw"])^-1)), y = bquote(log(`WD, m3 water`)))
# ggsave('Figure_pairplot_WD_lipidafdw_log.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# ggplot(micro_sample, aes(x=lipid_af_dw, y=log(`TETinf, kg 1.4-DC`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lipid_af_dw"~(g~.~(g["ash-free dw"])^-1)), y = bquote(log(`TETinf, kg 1.4-DC`)))
# ggsave('Figure_pairplot_logTET_lipidafdw.jpeg', width = 10, height = 7,dpi=600)
# 
# ggplot(micro_sample, aes(x=lipid_af_dw, y=`TETinf, kg 1.4-DC`,color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lipid_af_dw"~(g~.~(g["ash-free dw"])^-1)), y = bquote(`TETinf, kg 1.4-DC`))
# ggsave('Figure_pairplot_TET_lipidafdw_nolog.jpeg', width = 10, height = 7,dpi=600)
# 
# 





# lat ok


# log


# ggplot(micro_sample, aes(x=lat, y=log(`GW100, kg CO2-eq`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lat (°)"), y =  bquote("log(GW100"~(kg~CO^2-eq)~")"))
# ggsave('Figure_pairplot_logGW_latitude_correc.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# ggplot(micro_sample, aes(x=lat, y=log(`FE, kg P-eq`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lat (°)"), y =  bquote(log(`FE, kg P-eq`)))
# ggsave('Figure_pairplot_logFE_latitude.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# ggplot(micro_sample, aes(x=lat, y=log(`WD, m3 water`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lat (°)"), y =  bquote(log(`WD, m3 water`)))
# ggsave('Figure_pairplot_logWD_latitude.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# ggplot(micro_sample, aes(x=lat, y=log(`TETinf, kg 1.4-DC`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lat (°)"), y =  bquote(log(`TETinf, kg 1.4-DC`)))
# ggsave('Figure_pairplot_logTET_latitude.jpeg', width = 10, height = 7,dpi=600)



# No log

ggplot(micro_sample, aes(x=lat, y=`GW100, kg CO2-eq`,color=CNTR)) +
  geom_point(size=2,alpha=0.5)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lat (°)"), y =  bquote("GW100,"~kg~CO^2-eq))
ggsave('FigureS5.jpeg', width = 10, height = 7,dpi=600)

# 
# ggplot(micro_sample, aes(x=lat, y=`FE, kg P-eq`,color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lat (°)"), y =  bquote(`FE, kg P-eq`))
# ggsave('Figure_pairplot_nologFE_latitude_ok.jpeg', width = 10, height = 7,dpi=600)
# 
# ggplot(micro_sample, aes(x=lat, y=`WD, m3 water`,color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lat (°)"), y =  bquote("WD,"~m^3~water))
# ggsave('Figure_pairplot_nologWD_latitude_ok.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# ggplot(micro_sample, aes(x=lat, y=`TETinf, kg 1.4-DC`,color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("lat (°)"), y =  bquote(`TETinf, kg 1.4-DC`))
# ggsave('Figure_pairplot_nologTET_latitude_ok.jpeg', width = 10, height = 7,dpi=600)


# Topt  ok

# log

ggplot(micro_sample, aes(x=Topt, y=log(`GW100, kg CO2-eq`),color=CNTR)) +
  geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Topt (°C)"), y =  bquote("log(GW100"~(kg~CO^2-eq)~")"))
ggsave('FigureS6.jpeg', width = 10, height = 7,dpi=600)


# Regression lines

for(country in list_skarka){
  print(country)
  
  model_GW_topt <- lm(log(`GW100, kg CO2-eq`) ~ Topt,
                      data = micro_sample[micro_sample$CNTR==country,])
  print(summary(model_GW_topt))
  
}


# ggplot(micro_sample, aes(x=Topt, y=log(`TETinf, kg 1.4-DC`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Topt (°C)"), y =  bquote(log(`TETinf, kg 1.4-DC`)))
# ggsave('Figure_pairplot_TET_topt_log.jpeg', width = 10, height = 7,dpi=600)
# 
# ggplot(micro_sample, aes(x=Topt, y=log(`FE, kg P-eq`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Topt (°C)"), y =  bquote(log(`FE, kg P-eq`)))
# ggsave('Figure_pairplot_FE_topt_log.jpeg', width = 10, height = 7,dpi=600)
# 
# ggplot(micro_sample, aes(x=Topt, y=log(`WD, m3 water`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Topt (°C)"), y =  bquote(log(`WD, m3 water`)))
# ggsave('Figure_pairplot_WD_topt_log.jpeg', width = 10, height = 7,dpi=600)
# 




#Volumetric.productivity.kg.L.2.d.1

# ggplot(micro_sample, aes(x=`Volumetric productivity kg.L-2.d-1`, y=log(`GW100, kg CO2-eq`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Volumetric productivity"~(kg.~L^-2~.~d^-1)), y =  bquote("log(GW100"~(kg~CO^2-eq)~")"))
# ggsave('Figure_pairplot_GW_Volumetric.productivity.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# 
# ggplot(micro_sample, aes(x=`Volumetric productivity kg.L-2.d-1`, y=log(`TETinf, kg 1.4-DC`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Volumetric productivity"~(kg.~L^-2~.~d^-1)), y =  bquote(log(`TETinf, kg 1.4-DC`)))
# ggsave('Figure_pairplot_TET_Volumetric.productivity_log.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# ggplot(micro_sample, aes(x=`Volumetric productivity kg.L-2.d-1`, y=log(`FE, kg P-eq`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Volumetric productivity"~(kg.~L^-2~.~d^-1)), y =  bquote(log(`FE, kg P-eq`)))
# ggsave('Figure_pairplot_FE_Volumetric.productivity_log.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# ggplot(micro_sample, aes(x=`Volumetric productivity kg.L-2.d-1`, y=log(`WD, m3 water`),color=CNTR)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_brewer(palette="Spectral")+labs(x = bquote("Volumetric productivity"~(kg.~L^-2~.~d^-1)), y =  bquote(log(`WD, m3 water`)))
# ggsave('Figure_pairplot_WD_Volumetric.productivity_log.jpeg', width = 10, height = 7,dpi=600)
# 
# 




# Correct name IC
micro_sample$GW100_1kWh <- micro_sample$GWP100_1kWh



ggplot(micro_sample, aes(x=`Volumetric productivity kg.L-2.d-1`, y=log(`GW100, kg CO2-eq`),color=GW100_1kWh)) +
  geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_gradient(low = "yellow", high = "red", na.value = NA)+labs(x = bquote("Volumetric productivity"~(kg.~L^-2~.~d^-1)), y = bquote("log(GW100"~(kg~CO^2-eq)~")"))

ggsave('FigureS8.jpeg', width = 10, height = 7,dpi=600)



# 
# ggplot(micro_sample, aes(x=Topt, y=log(`GW100, kg CO2-eq`),color=GW100_1kWh)) +
#   geom_point(size=2,alpha=0.5)+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+scale_colour_gradient(low = "yellow", high = "red", na.value = NA)+labs(x = bquote("Topt (°C)"), y = bquote("log(GW100"~(kg~CO^2-eq)~")"))
# ggsave('Figure_pairplot_GW_topt_color=GW100_kwh_continuous_log.jpeg', width = 10, height = 7,dpi=600)
# 
# 



### Impact= f(Impact)
jpeg(file="FigureS7.jpg",height=2200, width=2600) 

pairs(micro_sample[,c("GW100, kg CO2-eq","WD, m3 water","FE, kg P-eq","TETinf, kg 1.4-DC")], pch = 15,  cex = 3,cex.labels = 6,cex.axis=4,
      col = c25[micro_sample$CNTR],
      lower.panel=NULL,oma=c(6,6,6,6))
par(xpd = TRUE)
legend("bottomleft",fill = c25[unique(micro_sample$CNTR)], legend = c( unique(micro_sample$CNTR)), title="Country",cex = 5, lwd = 5)

dev.off() 







# Create summary table


stat_complete = stat.desc(Total_results_table)

write.csv(stat_complete,"Stats_tot_mono_complete.csv", row.names = TRUE)


# Find maximum and minimum

# max_GW_line = Total_results_table[Total_results_table$`GW100, kg CO2-eq`==max(Total_results_table$`GW100, kg CO2-eq`),]
# min_GW_line = Total_results_table[Total_results_table$`GW100, kg CO2-eq`==min(Total_results_table$`GW100, kg CO2-eq`),]
# 
# max_FE_line = Total_results_table[Total_results_table$`FE, kg P-eq`==max(Total_results_table$`FE, kg P-eq`),]
# min_FE_line = Total_results_table[Total_results_table$`FE, kg P-eq`==min(Total_results_table$`FE, kg P-eq`),]
# 
# max_WD_line = Total_results_table[Total_results_table$`WD, m3 water`==max(Total_results_table$`WD, m3 water`),]
# min_WD_line = Total_results_table[Total_results_table$`WD, m3 water`==min(Total_results_table$`WD, m3 water`),]
# 
# max_TET_line = Total_results_table[Total_results_table$`TETinf, kg 1.4-DC`==max(Total_results_table$`TETinf, kg 1.4-DC`),]
# min_TET_line = Total_results_table[Total_results_table$`TETinf, kg 1.4-DC`==min(Total_results_table$`TETinf, kg 1.4-DC`),]
# 
# min_max_table=rbind(max_TET_line ,
#                     min_TET_line,
#                     max_WD_line,
#                     min_WD_line,
#                     max_FE_line,
#                     min_FE_line,
#                     max_GW_line,
#                     min_GW_line) 
# 
# rownames(min_max_table)=c("Max_TETinf" ,
#                           "Min_TETinf",
#                           "Max_WD_",
#                           "Min_WD",
#                           "Max_FE",
#                           "Min_FE",
#                           "Max_GW",
#                           "Min_GW")
# 
# 
# 
# 
# write.csv(min_max_table,"min_max_individuals.csv", row.names = TRUE)




### Influence categorical variables









# night monitoring

# 
# 
# ggplot(Total_melted, aes(x=night_monitoring, y=Score)) +facet_grid(rows = vars(`facets`),scales="free",labeller = label_parsed)+
#   geom_boxplot2(width = 0.8, width.errorbar = 0.05)+stat_summary(fun=mean, geom="point", shape=21, size=3,alpha=0.7)
# ggsave('boxplot_category_night monitoring_correc.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# 
# 
# #Stats table
# 
# 
# 
# Stats_night_mon=Total_melted %>%
#   group_by(market_for_substitution,`Impact_Category`) %>%
#   get_summary_stats(Score, show = c("mean","median","sd"))
# 
# write.csv(Stats_night_mon,"Stats_night_mon_sampling1.csv", row.names = FALSE)
# 
# 
# #Bio class
# 
# ggplot(Total_melted, aes(x=Bio_class, y=Score)) +facet_grid(rows = vars(`facets`),scales="free",labeller = label_parsed)+
#   geom_boxplot2(width = 0.8, width.errorbar = 0.05)+stat_summary(fun=mean, geom="point", shape=21, size=3,alpha=0.7)
# ggsave('boxplot_category_bioclass.jpeg', width = 10, height = 7,dpi=600)
# 
# 
# #Stats table
# 
# 
# Stats_bio=Total_melted %>%
#   group_by(Bio_class,`Impact_Category`) %>%
#   get_summary_stats(Score, show = c("mean","median","sd"))
# 
# write.csv(Stats_bio,"Stats_bio_sampling1.csv", row.names = FALSE)
# 
# 
# 
# 
# #GW
# 
# kruskal.test(Score ~ Bio_class, data = Total_melted[Total_melted$Impact_Category=='GW100, kg CO2-eq',])
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=='GW100, kg CO2-eq',]$Score, 
#                      Total_melted[Total_melted$Impact_Category=='GW100, kg CO2-eq',]$Bio_class,
#                      p.adjust.method = "BH")
# 
# 
# #WD
# kruskal.test(Score ~ Bio_class, data = Total_melted[Total_melted$Impact_Category=="WD, m3 water",])
# 
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=="WD, m3 water",]$Score, 
#                      Total_melted[Total_melted$Impact_Category=="WD, m3 water",]$Bio_class,
#                      p.adjust.method = "BH")
# 
# #FE
# kruskal.test(Score ~ Bio_class, data = Total_melted[Total_melted$Impact_Category=="FE, kg P-eq",])
# 
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=="FE, kg P-eq",]$Score, 
#                      Total_melted[Total_melted$Impact_Category=="FE, kg P-eq",]$Bio_class,
#                      p.adjust.method = "BH")
# 
# #TETinf
# kruskal.test(Score ~ Bio_class, data = Total_melted[Total_melted$Impact_Category=="TETinf, kg 1.4-DC",])
# 
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=="TETinf, kg 1.4-DC",]$Score, 
#                      Total_melted[Total_melted$Impact_Category=="TETinf, kg 1.4-DC",]$Bio_class,
#                      p.adjust.method = "BH")
# 
# 



# Substition

ggplot(Total_melted, aes(x=market_for_substitution, y=Score)) +facet_grid(rows = vars(`facets`),scales="free",labeller = label_parsed)+
  geom_boxplot2(width = 0.8, width.errorbar = 0.05)+stat_summary(fun=mean, geom="point", shape=21, size=3,alpha=0.7)

ggsave('FigureS9.jpeg', width = 10, height = 7,dpi=600)

#Stats table

# 
# 
# Stats_subst=Total_melted %>%
#   group_by(market_for_substitution,`Impact_Category`) %>%
#   get_summary_stats(Score, show = c("mean","median","sd"))
# 
# write.csv(Stats_subst,"Stats_substmark_sampling1.csv", row.names = FALSE)
# 
# 
# 
# #GW
# 
# kruskal.test(Score ~ market_for_substitution, data = Total_melted[Total_melted$Impact_Category=='GW100, kg CO2-eq',])
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=='GW100, kg CO2-eq',]$Score, 
#                      Total_melted[Total_melted$Impact_Category=='GW100, kg CO2-eq',]$market_for_substitution,
#                      p.adjust.method = "BH")
# 
# 
# #WD
# kruskal.test(Score ~ market_for_substitution, data = Total_melted[Total_melted$Impact_Category=="WD, m3 water",])
# 
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=="WD, m3 water",]$Score, 
#                      Total_melted[Total_melted$Impact_Category=="WD, m3 water",]$market_for_substitution,
#                      p.adjust.method = "BH")
# 
# #FE
# kruskal.test(Score ~ market_for_substitution, data = Total_melted[Total_melted$Impact_Category=="FE, kg P-eq",])
# 
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=="FE, kg P-eq",]$Score, 
#                      Total_melted[Total_melted$Impact_Category=="FE, kg P-eq",]$market_for_substitution,
#                      p.adjust.method = "BH")
# 
# #TETinf
# kruskal.test(Score ~ market_for_substitution, data = Total_melted[Total_melted$Impact_Category=="TETinf, kg 1.4-DC",])
# 
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=="TETinf, kg 1.4-DC",]$Score, 
#                      Total_melted[Total_melted$Impact_Category=="TETinf, kg 1.4-DC",]$market_for_substitution,
#                      p.adjust.method = "BH")
# 
# 
# 
# 
# Total_melted$Nsource
# 
# 
# #Nutrients
# 
# ggplot(Total_melted, aes(x=Nsource, y=Score)) +facet_grid(rows = vars(`facets`),scales="free",labeller = label_parsed)+
#   geom_boxplot2(width = 0.8, width.errorbar = 0.05)+stat_summary(fun=mean, geom="point", shape=21, size=3,alpha=0.7)
# 
# ggsave('boxplot_Nsource.jpeg', width = 10, height = 7,dpi=600)
# 
# #Stats table
# 
# 
# 
# Stats_subst=Total_melted %>%
#   group_by(Nsource,`Impact_Category`) %>%
#   get_summary_stats(Score, show = c("mean","median","sd"))
# 
# write.csv(Stats_subst,"Stats_nsource_sampling1.csv", row.names = FALSE)
# 
# 
# 
# #GW
# 
# kruskal.test(Score ~ market_for_substitution, data = Total_melted[Total_melted$Impact_Category=='GW100, kg CO2-eq',])
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=='GW100, kg CO2-eq',]$Score, 
#                      Total_melted[Total_melted$Impact_Category=='GW100, kg CO2-eq',]$market_for_substitution,
#                      p.adjust.method = "BH")
# 
# 
# #WD
# kruskal.test(Score ~ market_for_substitution, data = Total_melted[Total_melted$Impact_Category=="WD, m3 water",])
# 
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=="WD, m3 water",]$Score, 
#                      Total_melted[Total_melted$Impact_Category=="WD, m3 water",]$market_for_substitution,
#                      p.adjust.method = "BH")
# 
# #FE
# kruskal.test(Score ~ market_for_substitution, data = Total_melted[Total_melted$Impact_Category=="FE, kg P-eq",])
# 
# 
# pairwise.wilcox.test(Total_melted[Total_melted$Impact_Category=="FE, kg P-eq",]$Score, 
#                      Total_melted[Total_melted$Impact_Category=="FE, kg P-eq",]$market_for_substitution,
#                      p.adjust.method = "BH")
# 
# #TETinf
# kruskal.test(Score ~ market_for_substitution, data = Total_melted[Total_melted$Impact_Category=="TETinf, kg 1.4-DC",])
# 
#                      