### Script to process results and produce Figure 3 and Figures S14 to S21.



# Set working directory to the script location.  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


install.packages(c("ggridges", "ggplot2", "reshape2", "plyr", "dplyr", "tidyverse","plotly","GGaly","rio","ggthemes","foreach","parallel","doParallel","foreach"))
library(readxl)
library(ggplot2)
library(reshape2) 
library(plyr)
library(dplyr)
library(ggthemes) # Charger
library(tidyverse)
library(bigmemory)
library(plotly)
library(GGally)
library(rio)
library(data.table)
library(ggridges)
library(doParallel)




# Load data
list_skarka=c('ES','SE','IT','PT','UK','FR','EL','CY','IE','DE')



# To keep the columns' names

to_get_columns <- fread(file="../Outputs-Multi/results_table_df_pointUK_2.csv")



# Load data



list_result = list()




for (country in list_skarka){
  print(country)
  for(count in 0:10 ){
    
    name_file=paste("../Outputs-Multi/results_table_df_point",country,'_',as.character(count),".csv",sep='')
    
    if(file.exists(name_file)){
      #results_table =read.csv(name_file,sep=";")
      results_table =fread(name_file)
      
      list_result=append(list_result,list(results_table))
      
      #Total_results_table=rbind(Total_results_table,results_table)
      
    }else{print(c(name_file,"does not exist"))}
    
  }}

#combine list of dataframes in a big dataframe

Total_results_table=data.frame(rbindlist(list_result))


colnames(Total_results_table)=colnames(to_get_columns)



# delete list_result to save memory and to_get_columns to save memory

list_result=none
to_get_columns= none


# Set some columns as factors

Total_results_table$night_monitoring<-as.factor(Total_results_table$night_monitoring)
Total_results_table$CNTR<-as.factor(Total_results_table$CNTR)
Total_results_table$Bio_class<-as.factor(Total_results_table$Bio_class)
Total_results_table$market_for_substitution<-as.factor(Total_results_table$market_for_substitution)
Total_results_table$Nsource<-as.factor(Total_results_table$Nsource)
Total_results_table$alga<-as.factor(Total_results_table$alga)


names(Total_results_table)[names(Total_results_table) == "WDP"] <- "WD, m3 water"
names(Total_results_table)[names(Total_results_table) == "GWP100"] <- "GW100, kg CO2-eq"
names(Total_results_table)[names(Total_results_table) == "FEP"] <- "FE, kg P-eq"
names(Total_results_table)[names(Total_results_table) == "TETPinf"] <- "TETinf, kg 1.4-DC"




Total_melted=melt(Total_results_table,measure.vars = c("WD, m3 water","GW100, kg CO2-eq","FE, kg P-eq","TETinf, kg 1.4-DC"),variable.name = 'Impact_Category',value.name='Score')

Total_melted$Impact_Category=as.factor(Total_melted$Impact_Category)





# Calculate quantiles for mean impact per strain-compound

Total_results_table=Total_results_table %>% group_by(alga)%>% 
  mutate(Mean_alga_GW = mean(`GW100, kg CO2-eq`)) %>% # calculate mean for plotting as well
  ungroup()

Total_results_table=Total_results_table %>% group_by(alga)%>% 
  mutate(Mean_alga_FE = mean(`FE, kg P-eq`)) %>% # calculate mean for plotting as well
  ungroup()

Total_results_table=Total_results_table %>% group_by(alga)%>% 
  mutate(Mean_alga_TET = mean(`TETinf, kg 1.4-DC`)) %>% # calculate mean for plotting as well
  ungroup()

Total_results_table=Total_results_table %>% group_by(alga)%>% 
  mutate(Mean_alga_WD = mean(`WD, m3 water`)) %>% # calculate mean for plotting as well
  ungroup()

quantile_number=200

colnames(Total_results_table)
Total_results_table <- within(Total_results_table, Quantile_GW_alga <- as.integer(cut(Mean_alga_GW, quantile(Mean_alga_GW, probs=0:quantile_number/quantile_number), include.lowest=TRUE)))
Total_results_table <- within(Total_results_table, Quantile_FE_alga <- as.integer(cut(Mean_alga_FE, quantile(Mean_alga_FE, probs=0:quantile_number/quantile_number), include.lowest=TRUE)))
Total_results_table <- within(Total_results_table, Quantile_TET_alga <- as.integer(cut(Mean_alga_TET, quantile(Mean_alga_TET, probs=0:quantile_number/quantile_number), include.lowest=TRUE)))
Total_results_table <- within(Total_results_table, Quantile_WD_alga <- as.integer(cut(Mean_alga_WD, quantile(Mean_alga_WD, probs=0:quantile_number/quantile_number), include.lowest=TRUE)))


####GW 


# Pick 200 random strain-compound pairs with 1 per quantile

number_points_per_quantile = 1
list_random_algae=list()


# To obtain same plot as Figure 2
set.seed(1)

#  Random sample within quantiles
for(i in 1:quantile_number){
  
  A=unique(Total_results_table[which(Total_results_table[,"Quantile_GW_alga"]==i),]$alga)
  
  
  random_algae_in_this_quant=sample(A,number_points_per_quantile , replace = FALSE)
  
  list_random_algae=append(list_random_algae,random_algae_in_this_quant)
  
}




# subset with only the algae to plot ( chosen in different quantiles)
sample_algae=Total_results_table[Total_results_table$alga %in% list_random_algae,]



# Figure3


ggplot(mutate(sample_algae, alga = reorder(alga, Quantile_GW_alga))) +
  geom_density_ridges(stat = "density",trim = TRUE,aes(x = `GW100, kg CO2-eq`, y = alga, group = interaction(CNTR,alga),fill = CNTR,height = ..density..), alpha = 0.7, scale = 150)+coord_cartesian(xlim = c(-20, 1500))+scale_fill_brewer(palette="Spectral")+
  labs(x = bquote("GW100"~(kg~CO^2-eq)), y =  bquote("Strain-compound tandem"))+theme(axis.text.y=element_text(size=0),
                                                                                       axis.text.x=element_text(size=14),
                                                                                       axis.title.x=element_text(size=16),
                                                                                       axis.title.y=element_text(size=16),
                                                                                       legend.title = element_text(size=13),
                                                                                       legend.text = element_text(size=13),
                                                                                       legend.key.size = unit(2, 'cm'),
                                                                                       legend.key.height = unit(0.2, 'cm'),
                                                                                       legend.key.width = unit(0.4, 'cm'))


ggsave('Figure_3_GW.pdf', width = 10, height = 7,dpi=600)






# With bioac_molec_dbio as color


ggplot(mutate(sample_algae, alga = reorder(alga, Quantile_GW_alga))) +
  geom_density_ridges(stat = "density",trim = TRUE,aes(x = `GW100, kg CO2-eq`, y = alga, group = interaction(CNTR,alga),fill = bioact_molec_dbio,height = ..density..), alpha = 0.7, scale = 150)+coord_cartesian(xlim = c(-20, 1500))+
  scale_fill_gradient2(low="yellow", high="red")+
  labs(x = bquote("GW100"~(kg~CO^2-eq)), y =  bquote("Strain-compound tandem"))+theme(axis.text.y=element_text(size=0),
                                                                                       axis.text.x=element_text(size=14),
                                                                                       axis.title.x=element_text(size=16),
                                                                                       axis.title.y=element_text(size=16),
                                                                                       legend.title = element_text(size=13),
                                                                                       legend.text = element_text(size=13),
                                                                                       legend.key.size = unit(2, 'cm'),
                                                                                       legend.key.height = unit(1, 'cm'),
                                                                                       legend.key.width = unit(0.5, 'cm')) + labs(fill =expression(paste(bioact_molec_dbio~(g.(g[dbio])^-1)))) 
ggsave('Figure_S14.pdf', width = 13, height = 7,dpi=600)



# With Topt as color


breaks=c(20,25,30)

ggplot(mutate(sample_algae, alga = reorder(alga, Quantile_GW_alga))) +
  geom_density_ridges(stat = "density",trim = TRUE,aes(x = `GW100, kg CO2-eq`, y = alga, group = interaction(CNTR,alga),fill = Topt,height = ..density..), alpha = 0.7, scale = 150)+coord_cartesian(xlim = c(-20, 1500))+
  scale_fill_gradientn(colours = c("blue","white","red"),
                       breaks= breaks, labels = format(breaks))+
  labs(x = bquote("GW100"~(kg~CO^2-eq)), y =  bquote("Strain-compound tandem"))+theme(axis.text.y=element_text(size=0),
                                                                                       axis.text.x=element_text(size=14),
                                                                                       axis.title.x=element_text(size=16),
                                                                                       axis.title.y=element_text(size=16),
                                                                                       legend.title = element_text(size=13),
                                                                                       legend.text = element_text(size=13),
                                                                                       legend.key.size = unit(2, 'cm'),
                                                                                       legend.key.height = unit(1, 'cm'),
                                                                                       legend.key.width = unit(0.5, 'cm')) + labs(fill ="Topt (°)") 


ggsave('Figure_S15.pdf', width = 13, height = 7,dpi=600)





# 
# #T_plateau 
# 
# ggplot(mutate(sample_algae, alga = reorder(alga, Quantile_GW_alga))) +
#   geom_density_ridges(stat = "density",trim = TRUE,aes(x = `GW100, kg CO2-eq`, y = alga, group = interaction(CNTR,alga),fill = T_plateau,height = ..density..), alpha = 0.7, scale = 150)+coord_cartesian(xlim = c(-20, 1500))+
#   labs(x = bquote("GW100"~(kg~CO^2-eq)), y =  bquote("Strain-compound tandem"))+theme(axis.text.y=element_text(size=0),
#                                                                                        axis.text.x=element_text(size=14),
#                                                                                        axis.title.x=element_text(size=16),
#                                                                                        axis.title.y=element_text(size=16),
#                                                                                        legend.title = element_text(size=13),
#                                                                                        legend.text = element_text(size=13),
#                                                                                        legend.key.size = unit(2, 'cm'),
#                                                                                        legend.key.height = unit(0.2, 'cm'),
#                                                                                        legend.key.width = unit(0.4, 'cm'))+ guides(fill=guide_legend(title="Tplateau (°)"))
# 
# ggsave('Figure_3_GW_appendix_Tpalteau.pdf', width = 10, height = 7,dpi=600)
# 
# 


# WD


ggplot(mutate(sample_algae, alga = reorder(alga, Quantile_WD_alga))) +
  geom_density_ridges(stat = "density",trim = TRUE,aes(x = `WD, m3 water`, y = alga, group = interaction(CNTR,alga),fill = CNTR,height = ..density..), alpha = 0.7, scale = 150)+coord_cartesian(xlim = c(-5, 100))+scale_fill_brewer(palette="Spectral")+
  labs(x = bquote("WD"~(m^3~"water")), y =  bquote("Strain-compound tandem"))+theme(axis.text.y=element_text(size=0),
                                                                                     axis.text.x=element_text(size=14),
                                                                                     axis.title.x=element_text(size=16),
                                                                                     axis.title.y=element_text(size=16),
                                                                                     legend.title = element_text(size=13),
                                                                                     legend.text = element_text(size=13),
                                                                                     legend.key.size = unit(2, 'cm'),
                                                                                     legend.key.height = unit(0.2, 'cm'),
                                                                                     legend.key.width = unit(0.4, 'cm'))
ggsave('Figure_3_WD.pdf', width = 10, height = 7,dpi=600)



# Bioact dbio 

ggplot(mutate(sample_algae, alga = reorder(alga, Quantile_WD_alga))) +
  geom_density_ridges(stat = "density",trim = TRUE,aes(x = `WD, m3 water`, y = alga, group = interaction(CNTR,alga),fill = bioact_molec_dbio,height = ..density..), alpha = 0.7, scale = 150)+coord_cartesian(xlim = c(-5, 100))+
  labs(x = bquote("WD"~(m^3~"water")), y =  bquote("Strain-compound tandem"))+theme(axis.text.y=element_text(size=0),
                                                                                     axis.text.x=element_text(size=14),
                                                                                     axis.title.x=element_text(size=16),
                                                                                     axis.title.y=element_text(size=16),
                                                                                     legend.title = element_text(size=13),
                                                                                     legend.text = element_text(size=13),
                                                                                     legend.key.size = unit(2, 'cm'),
                                                                                     legend.key.height = unit(0.2, 'cm'),
                                                                                     legend.key.width = unit(0.4, 'cm'))+ guides(fill=guide_legend(title=expression(paste(bioact_molec_dbio~(g.(g[dbio])^-1)))))
ggsave('Figure_3_WD_bioact_dbio.pdf', width = 10, height = 7,dpi=600)

# 
# # Topt 
# 
# colnames(sample_algae)
# ggplot(mutate(sample_algae, alga = reorder(alga, Quantile_WD_alga))) +
#   geom_density_ridges(stat = "density",trim = TRUE,aes(x = `WD, m3 water`, y = alga, group = interaction(CNTR,alga),fill = Topt,height = ..density..), alpha = 0.7, scale = 150)+coord_cartesian(xlim = c(-5, 100))+
#   labs(x = bquote("WD"~(m^3~"water")), y =  bquote("Strain-compound tandem"))+theme(axis.text.y=element_text(size=0),
#                                                                                      axis.text.x=element_text(size=14),
#                                                                                      axis.title.x=element_text(size=16),
#                                                                                      axis.title.y=element_text(size=16),
#                                                                                      legend.title = element_text(size=13),
#                                                                                      legend.text = element_text(size=13),
#                                                                                      legend.key.size = unit(2, 'cm'),
#                                                                                      legend.key.height = unit(0.2, 'cm'),
#                                                                                      legend.key.width = unit(0.4, 'cm'))+ guides(fill=guide_legend(title="Topt (°)"))
# ggsave('Figure_3_WD_topt.pdf', width = 10, height = 7,dpi=600)
# 


#TET

ggplot(mutate(sample_algae, alga = reorder(alga, Quantile_TET_alga))) +
  geom_density_ridges(stat = "density",trim = TRUE,aes(x = `TETinf, kg 1.4-DC`, y = alga, group = interaction(CNTR,alga),fill = CNTR,height = ..density..), alpha = 0.7, scale = 50)+coord_cartesian(xlim = c(-2, 2))+scale_fill_brewer(palette="Spectral")+
  labs(x = bquote("TETinf (kg 1.4-DC)"), y =  bquote("Strain-compound tandem"))+theme(axis.text.y=element_text(size=0),
                                                                                       axis.text.x=element_text(size=14),
                                                                                       axis.title.x=element_text(size=16),
                                                                                       axis.title.y=element_text(size=16),
                                                                                       legend.title = element_text(size=13),
                                                                                       legend.text = element_text(size=13),
                                                                                       legend.key.size = unit(2, 'cm'),
                                                                                       legend.key.height = unit(0.2, 'cm'),
                                                                                       legend.key.width = unit(0.4, 'cm'))

ggsave('Figure_3_TET.pdf', width = 10, height = 7,dpi=600)







#FE

ggplot(mutate(sample_algae, alga = reorder(alga, Quantile_FE_alga))) +
  geom_density_ridges(stat = "density",trim = TRUE,aes(x = `FE, kg P-eq`, y = alga, group = interaction(CNTR,alga),fill = CNTR,height = ..density..), alpha = 0.7, scale = 150)+coord_cartesian(xlim = c(-0.5, 1))+scale_fill_brewer(palette="Spectral")+
  labs(x = bquote("FE (kg P-eq)"), y =  bquote("Strain-compound tandem"))+theme(axis.text.y=element_text(size=0),
                                                                                 axis.text.x=element_text(size=14),
                                                                                 axis.title.x=element_text(size=16),
                                                                                 axis.title.y=element_text(size=16),
                                                                                 legend.title = element_text(size=13),
                                                                                 legend.text = element_text(size=13),
                                                                                 legend.key.size = unit(2, 'cm'),
                                                                                 legend.key.height = unit(0.2, 'cm'),
                                                                                 legend.key.width = unit(0.4, 'cm'))

ggsave('Figure_3_FE.pdf', width = 10, height = 7,dpi=600)

# 
# # Bioact dbio
# 
# ggplot(mutate(sample_algae, alga = reorder(alga, Quantile_FE_alga))) +
#   geom_density_ridges(stat = "density",trim = TRUE,aes(x = `FE, kg P-eq`, y = alga, group = interaction(CNTR,alga),fill = bioact_molec_dbio,height = ..density..), alpha = 0.7, scale = 150)+coord_cartesian(xlim = c(-0.5, 1))+
#   labs(x = bquote("FE (kg P-eq)"), y =  bquote("Strain-compound tandem"))+theme(axis.text.y=element_text(size=0),
#                                                                                  axis.text.x=element_text(size=14),
#                                                                                  axis.title.x=element_text(size=16),
#                                                                                  axis.title.y=element_text(size=16),
#                                                                                  legend.title = element_text(size=13),
#                                                                                  legend.text = element_text(size=13),
#                                                                                  legend.key.size = unit(2, 'cm'),
#                                                                                  legend.key.height = unit(0.2, 'cm'),
#                                                                                  legend.key.width = unit(0.4, 'cm'))+ guides(fill=guide_legend(title=expression(paste(bioact_molec_dbio~(g.(g[dbio])^-1)))))
# ggsave('Figure_3_FE_bioact_dbio.pdf', width = 10, height = 7,dpi=600)
# 
# 
# 






# Heteroskedasticity


# Calculate standard deviations for each strain-compound

Total_results_table=Total_results_table %>% group_by(alga)%>% 
  mutate(sd_alga_WD = sd(`WD, m3 water`)) %>% # calculate mean for plotting as well
  ungroup()

Total_results_table=Total_results_table %>% group_by(alga)%>% 
  mutate(sd_alga_GW = sd(`GW100, kg CO2-eq`)) %>% # calculate mean for plotting as well
  ungroup()

Total_results_table=Total_results_table %>% group_by(alga)%>% 
  mutate(sd_alga_TET = sd(`TETinf, kg 1.4-DC`)) %>% # calculate mean for plotting as well
  ungroup()

Total_results_table=Total_results_table %>% group_by(alga)%>% 
  mutate(sd_alga_FE = sd(`FE, kg P-eq`)) %>% # calculate mean for plotting as well
  ungroup()




Total_results_table_sd_plot = distinct_at(as_tibble(Total_results_table),vars(alga,sd_alga_GW,sd_alga_WD,sd_alga_TET,sd_alga_FE,Mean_alga_WD,Mean_alga_TET,Mean_alga_FE,Mean_alga_GW,bioact_molec_dbio,Topt))


# Calculate Coefficient of variation for each strain-compound

Total_results_table_sd_plot$CV_GW100 = Total_results_table_sd_plot$sd_alga_GW/Total_results_table_sd_plot$Mean_alga_GW
Total_results_table_sd_plot$CV_TET = Total_results_table_sd_plot$sd_alga_TET/Total_results_table_sd_plot$Mean_alga_TET
Total_results_table_sd_plot$CV_FE = Total_results_table_sd_plot$sd_alga_FE/Total_results_table_sd_plot$Mean_alga_FE
Total_results_table_sd_plot$CV_WD = Total_results_table_sd_plot$sd_alga_WD/Total_results_table_sd_plot$Mean_alga_WD



# PLOT


breaks <- c(20,25,30)



ggplot(Total_results_table_sd_plot, aes(x=Mean_alga_GW, y=sd_alga_GW,col=Topt)) + labs(x = bquote("Mean(GW100) per strain-compound tandem"~(kg~CO^2-eq)), y =  bquote("Sd(GW100) per strain-compound tandem"~(kg~CO^2-eq)))+
  geom_point(size=3)+scale_colour_gradientn(colours = c("blue","white","red"),
                                            breaks= breaks, labels = format(breaks))+ guides(fill=guide_legend(title="Topt (°)"))+theme(axis.text.y=element_text(size=14),
                                                                                                                                        axis.text.x=element_text(size=14),
                                                                                                                                        axis.title.x=element_text(size=17),
                                                                                                                                        axis.title.y=element_text(size=17),
                                                                                                                                        legend.title = element_text(size=13),
                                                                                                                                        legend.text = element_text(size=13),
                                                                                                                                        legend.key.size = unit(2, 'cm'),
                                                                                                                                        legend.key.height = unit(1.5, 'cm'),
                                                                                                                                        legend.key.width = unit(1, 'cm'))+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+geom_abline(intercept = 0, slope = 1, col="red")+
  labs(col="Topt (°)")



ggsave('FigureS16.jpg', width = 14, height = 10,dpi=600)






# X =meanGW  Y = coefficient of variation


ggplot(Total_results_table_sd_plot, aes(x=Mean_alga_GW, y=CV_GW100,colour=Topt)) + labs(x = bquote("Mean(GW100) per strain-compound tandem"~(kg~CO^2-eq)), y =  bquote("CV(GW100) per strain-compound tandem"))+
  geom_point()+theme(axis.text.y=element_text(size=14),
                     axis.text.x=element_text(size=14),
                     axis.title.x=element_text(size=17),
                     axis.title.y=element_text(size=17),
                     legend.title = element_text(size=13),
                     legend.text = element_text(size=13),
                     legend.key.size = unit(2, 'cm'),
                     legend.key.height = unit(1.5, 'cm'),
                     legend.key.width = unit(1, 'cm'))+labs(col="Topt (°)")
ggsave('FigureS19.jpg', width = 14, height = 10,dpi=600)





# X =TOpt Y = sd



ggplot(Total_results_table_sd_plot, aes(x=Topt, y=sd_alga_GW)) + labs(x = bquote("Topt (°)"), y =  bquote("Sd(GW100) per strain-compound tandem"~(kg~CO^2-eq)))+
  geom_point()+theme(axis.text.y=element_text(size=14),
                     axis.text.x=element_text(size=14),
                     axis.title.x=element_text(size=17),
                     axis.title.y=element_text(size=17))
ggsave('FigureS18.jpg', width = 14, height = 10,dpi=600)


# X =TOpt Y = CV



ggplot(Total_results_table_sd_plot, aes(x=Topt, y=CV_GW100)) + labs(x = bquote("Topt (°)"), y= bquote("CV(GW100) per strain-compound tandem"))+
  geom_point()+theme(axis.text.y=element_text(size=14),
                     axis.text.x=element_text(size=14),
                     axis.title.x=element_text(size=17),
                     axis.title.y=element_text(size=17))
ggsave('FigureS20.jpg', width = 14, height = 10,dpi=600)




# X = bioact-dbio   Y = Sd




ggplot(Total_results_table_sd_plot, aes(x=bioact_molec_dbio, y=sd_alga_GW)) + labs(x = bquote("bioact_molec_dbio"~(g.(g[dbio])^-1)), y =  bquote("Sd(GW100) per strain-compound tandem"~(kg~CO^2-eq)))+
  geom_point()+theme(axis.text.y=element_text(size=14),
                     axis.text.x=element_text(size=14),
                     axis.title.x=element_text(size=17),
                     axis.title.y=element_text(size=17))
ggsave('FigureS17.jpg', width = 14, height = 10,dpi=600)


# X =bioact-dbio Y = CV



ggplot(Total_results_table_sd_plot, aes(x=bioact_molec_dbio, y=CV_GW100,colour=Topt)) + labs(x = bquote("bioact_molec_dbio"~(g.(g[dbio])^-1)), y= bquote("CV(GW100) per strain-compound tandem"))+
  geom_point()+theme(axis.text.y=element_text(size=14),
                     axis.text.x=element_text(size=14),
                     axis.title.x=element_text(size=17),
                     axis.title.y=element_text(size=17),
                     legend.title = element_text(size=13),
                     legend.text = element_text(size=13),
                     legend.key.size = unit(2, 'cm'),
                     legend.key.height = unit(1.5, 'cm'),
                     legend.key.width = unit(1, 'cm'))+labs(col="Topt (°)")
ggsave('FigureS21.jpg', width = 14, height = 10,dpi=600)









# 
# 
# #FE
# 
# # X =meanFE Y = sd
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=Mean_alga_FE, y=sd_alga_WD)) + labs(x = bquote("Mean FE per strain-compound tandem FE (kg P-eq)"), y =  bquote("Standard deviation FE per strain-compound tandem (kg P-eq)"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=16),
#                      axis.title.y=element_text(size=16))
# ggsave('hetero_sdFE_meanFE_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# 
# 
# 
# # X =meanFE  Y = coefficient of variation
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=Mean_alga_FE, y=CV_FE,colour=Topt)) + labs(x = bquote("Mean FE per strain-compound tandem (kg P-eq)"), y =  bquote("Coefficient of Variariation for FE per strain-compound tandem"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=16),
#                      axis.title.y=element_text(size=16))
# ggsave('hetero_meanFE_cvFE_toptcolor_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# 
# 
# 
# # X =TOpt Y = sd
# 
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=Topt, y=sd_alga_FE)) + labs(x = bquote("Topt (°)"), y =  bquote("Standard deviation FE per strain-compound tandem (kg P-eq)"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=16),
#                      axis.title.y=element_text(size=16))
# ggsave('hetero_topt_sdFE.jpg', width = 10, height = 10,dpi=600)
# 
# 
# # X =TOpt Y = CV
# 
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=Topt, y=CV_FE)) + labs(x = bquote("Topt (°)"), y= bquote("Coefficient of Variariation for FE per strain-compound tandem"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=16),
#                      axis.title.y=element_text(size=16))
# ggsave('hetero_2_Topt_CVFE_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# 
# 
# # X = bioact-dbio   Y = Sd
# 
# 
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=bioact_molec_dbio, y=sd_alga_FE)) + labs(x = bquote("bioact_molec_dbio"~(g.(g[dbio])^-1)), y =  bquote("Standard deviation FE per strain-compound tandem (kg P-eq)"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=15),
#                      axis.title.y=element_text(size=15))
# ggsave('hetero_bioact_sdFE_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# # X =bioact-dbio Y = CV
# 
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=bioact_molec_dbio, y=CV_FE,colour=Topt)) + labs(x = bquote("bioact_molec_dbio"~(g.(g[dbio])^-1)), y= bquote("Coefficient of Variariation for FE per strain-compound tandem"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=15),
#                      axis.title.y=element_text(size=15))
# ggsave('hetero_2_bioact_dbio_CVFE_toptcol_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #WD
# 
# withoutoutlier_WD=Total_results_table_sd_plot[which(Total_results_table_sd_plot$alga!=764),]
# 
# 
# # X =meanFE Y = sd
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=Mean_alga_WD, y=sd_alga_WD)) + labs(x = bquote("Mean WD per strain-compound tandem" ~(m^3~"water")), y =  bquote("Standard deviation WD per strain-compound tandem"~(m^3~"water")))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=16),
#                      axis.title.y=element_text(size=16))
# ggsave('hetero_sdWD_meanWD_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# 
# 
# 
# # X =meanWD  Y = coefficient of variation
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=Mean_alga_WD, y=CV_WD,colour=Topt)) + labs(x = bquote("Mean WD per strain-compound tandem" ~(m^3~"water")), y =  bquote("Coefficient of Variariation for WD per strain-compound tandem"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=14),
#                      axis.title.y=element_text(size=14))
# ggsave('hetero_meanWD_cvWD_toptcolor_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# 
# 
# ggplot(withoutoutlier_WD, aes(x=Mean_alga_WD, y=CV_WD,colour=Topt)) + labs(x = bquote("Mean WD per strain-compound tandem" ~(m^3~"water")), y =  bquote("Coefficient of Variariation for WD per strain-compound tandem"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=12),
#                      axis.title.y=element_text(size=12))
# ggsave('hetero_meanWD_cvWD_toptcolor_nooutlier_goodsize.jpg', width = 10, height = 10,dpi=600)
# 

# 
# 
# 
# # X =TOpt Y = sd
# 
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=Topt, y=sd_alga_WD)) + labs(x = bquote("Topt (°)"), y =  bquote("Standard deviation WD per strain-compound tandem (kg P-eq)"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=16),
#                      axis.title.y=element_text(size=16))
# ggsave('hetero_topt_sdWD_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# # X =TOpt Y = CV
# 
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=Topt, y=CV_WD)) + labs(x = bquote("Topt (°)"), y= bquote("Coefficient of Variariation for WD per strain-compound tandem"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=16),
#                      axis.title.y=element_text(size=16))
# ggsave('hetero_2_Topt_CVWD_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# 
# 
# # X = bioact-dbio   Y = Sd
# 
# 
# 
# 
# ggplot(withoutoutlier_WD, aes(x=bioact_molec_dbio, y=sd_alga_WD)) + labs(x = bquote("bioact_molec_dbio"~(g.(g[dbio])^-1)), y =  bquote("Standard deviation WD per strain-compound tandem (kg P-eq)"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=15),
#                      axis.title.y=element_text(size=15))
# ggsave('hetero_bioact_sdWD_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# # X =bioact-dbio Y = CV
# 
# 
# 
# ggplot(Total_results_table_sd_plot, aes(x=bioact_molec_dbio, y=CV_WD,colour=Topt)) + labs(x = bquote("bioact_molec_dbio"~(g.(g[dbio])^-1)), y= bquote("Coefficient of Variariation for WD per strain-compound tandem"))+
#   geom_point()+theme(axis.text.y=element_text(size=14),
#                      axis.text.x=element_text(size=14),
#                      axis.title.x=element_text(size=14),
#                      axis.title.y=element_text(size=14))
# ggsave('hetero_2_bioact_dbio_CVWD_toptcol_goodsize.jpg', width = 10, height = 10,dpi=600)
# 
# 
# 
# 












