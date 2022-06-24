### Script to process the sensitivity indices and plot Figure 5 and Figure S22 to S23.

### You need to replace the ID of the sensitivy indices files

# Set working directory to the script location.  

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
library(ggpattern)

#Execute the whole script to print figure 2, 3 and 4 from the article.



getwd()

#Load data
list_skarka=c('ES','SE','IT','PT','UK','FR','EL','CY','IE','DE')





# Load Multidimensional indices

## Change the ID of the xlsx file with non-weighted indice to the match the file in Outputs-Multi

Total_sens_table_2 <- read_excel("../Outputs-Multi/sensi_multi_sampling_2_non_weighted_5_19_058642.xlsx")



# Load Mono-dimensional indices


to_get_columns_1 <- read_excel("../Outputs-Mono/sensi_multi_melt_pointEL_1.xlsx")



list_result = list()


for (country in list_skarka){
  print(country)
  for(count in 0:8 ){
    
    name_file=paste("../From_ubuntu_FinalSampling1_Unc_AD_1105/Outputs/sensi_multi_melt_point",country,'_',as.character(count),".xlsx",sep='')
    
    if(file.exists(name_file)){
      sens_table =read_excel(name_file)
      
      list_result=append(list_result,list(sens_table))
      
    }else{print(c(name_file,"does not exist"))}
    
  }}



Total_sens_table_1=data.frame(rbindlist(list_result))

colnames(Total_sens_table_1)=colnames(to_get_columns_1)








# Refine names and reorganize results for Multidimensional sampling


Total_sens_table_2$Impact_Category = as.factor(Total_sens_table_2$Impact_Category)
Total_sens_table_2$Impact_Category = as.factor(Total_sens_table_2$Impact_Category)



Total_sens_table_2$Parameter = as.factor(Total_sens_table_2$Parameter)

names(Total_sens_table_2)[names(Total_sens_table_2) == "Impact_Category"] <- "Impact Category"




ST_2 = subset(Total_sens_table_2,`Type indice`== "ST")
S1_2 = subset(Total_sens_table_2,`Type indice`== "S1")
conf_S1_2 =  subset(Total_sens_table_2,`Type indice`== "S1_conf")
conf_ST_2 =  subset(Total_sens_table_2,`Type indice`== "ST_conf")





# Recombine indices and their confidence intervals

final_S1_2=cbind(S1_2,"conf"=conf_S1_2$value)
final_ST_2=cbind(ST_2,"conf"=conf_ST_2$value)






# Add column with sum of indices for each category



final_ST_unique_2 <- final_ST_2%>% group_by(`Impact Category`)%>% 
  mutate(sumST = sum(value)) %>% # calculate mean for plotting as well
  ungroup()



final_S1_unique_2 <- final_S1_2 %>% group_by(`Impact Category`)%>% 
  mutate(sumS1 = sum(value)) %>% # calculate mean for plotting as well
  ungroup()


# Create column with % of sum of indices



final_S1_unique_2$Percent_of_total = final_S1_unique_2$value/final_S1_unique_2$sumS1
final_ST_unique_2$Percent_of_total = final_ST_unique_2$value/final_ST_unique_2$sumST














# Refine names and reorganize results for Monodimensional sampling







Total_sens_table_1$CNTR = as.factor(Total_sens_table_1$CNTR)
Total_sens_table_1$Impact_Category = as.factor(Total_sens_table_1$Impact_Category)
Total_sens_table_1$Impact_Category = as.factor(Total_sens_table_1$Impact_Category)


Total_sens_table_1$Parameter=as.factor(Total_sens_table_1$Parameter)

names(Total_sens_table_1)[names(Total_sens_table_1) == "Impact_Category"] <- "Impact Category"





ST_1 = subset(Total_sens_table_1,`Type indice`== "ST")
S1_1 = subset(Total_sens_table_1,`Type indice`== "S1")
conf_S1_1 =  subset(Total_sens_table_1,`Type indice`== "S1_conf")
conf_ST_1 =  subset(Total_sens_table_1,`Type indice`== "ST_conf")




# Add mean per country

S1_1 <- S1_1 %>% group_by(CNTR,Parameter,`Impact Category`)%>% 
  mutate(Mean_cntr = mean(value)) %>% # calculate mean for plotting as well
  ungroup()

conf_S1_1 <- conf_S1_1 %>% group_by(CNTR,Parameter,`Impact Category`)%>% 
  mutate(Mean_cntr = mean(value)) %>% # calculate mean for plotting as well
  ungroup()


ST_1 <- ST_1 %>% group_by(CNTR,Parameter,`Impact Category`)%>% 
  mutate(Mean_cntr = mean(value)) %>% # calculate mean for plotting as well
  ungroup()

conf_ST_1 <- conf_ST_1 %>% group_by(CNTR,Parameter,`Impact Category`)%>% 
  mutate(Mean_cntr = mean(value)) %>% # calculate mean for plotting as well
  ungroup()




S1_only_national_averages_1 = unique(S1_1[,c('Mean_cntr','Parameter','CNTR','Impact Category')])
conf_S1_only_national_averages_1 = unique(conf_S1_1[,c('Mean_cntr','Parameter','CNTR','Impact Category')])
ST_only_national_averages_1 = unique(ST_1[,c('Mean_cntr','Parameter','CNTR','Impact Category')])
conf_ST_only_national_averages_1 = unique(conf_ST_1[,c('Mean_cntr','Parameter','CNTR','Impact Category')])




S1_total_average_1 <- S1_only_national_averages_1 %>% group_by(Parameter,`Impact Category`)%>% 
  mutate(Mean_tot = mean(Mean_cntr)) %>% # calculate mean for plotting as well
  ungroup()

conf_S1_total_average_1 <- conf_S1_only_national_averages_1 %>% group_by(Parameter,`Impact Category`)%>% 
  mutate(Mean_tot = mean(Mean_cntr)) %>% # calculate mean for plotting as well
  ungroup()

ST_total_average_1 <- ST_only_national_averages_1 %>% group_by(Parameter,`Impact Category`)%>% 
  mutate(Mean_tot = mean(Mean_cntr)) %>% # calculate mean for plotting as well
  ungroup()

conf_ST_total_average_1 <- conf_ST_only_national_averages_1 %>% group_by(Parameter,`Impact Category`)%>% 
  mutate(Mean_tot = mean(Mean_cntr)) %>% # calculate mean for plotting as well
  ungroup()


# Recombine indices and their confidence intervals

final_S1_1=cbind(S1_total_average_1,"conf" = conf_S1_total_average_1$Mean_tot)
final_ST_1=cbind(ST_total_average_1,"conf" = conf_ST_total_average_1$Mean_tot)



final_S1_unique_1 = unique(final_S1_1[,c('Mean_tot','Parameter','Impact Category','conf')])


final_ST_unique_1 = unique(final_ST_1[,c('Mean_tot','Parameter','Impact Category','conf')])





final_S1_unique_IC_1 = subset(final_S1_unique_1,`Impact Category`==	'GWP100')
final_ST_unique_IC_1 = subset(final_ST_unique_1,`Impact Category`==	'GWP100')





# Add column with sum of indices for each category



final_ST_unique_1 <- final_ST_unique_1 %>% group_by(`Impact Category`)%>% 
  mutate(sumST = sum(Mean_tot)) %>% # calculate mean for plotting as well
  ungroup()



final_S1_unique_1 <- final_S1_unique_1 %>% group_by(`Impact Category`)%>% 
  mutate(sumS1 = sum(Mean_tot)) %>% # calculate mean for plotting as well
  ungroup()


# Create column with % of sum of indices



final_S1_unique_1$Percent_of_total_1 = final_S1_unique_1$Mean_tot/final_S1_unique_1$sumS1
final_ST_unique_1$Percent_of_total_1 = final_ST_unique_1$Mean_tot/final_ST_unique_1$sumST







# Combine Monodimensional and multidimensional sampling



# Sampling 1
subset_finalST_1 = final_ST_unique_1[,c("Mean_tot","Parameter","Impact Category","Percent_of_total_1")]
colnames(subset_finalST_1) = c("value","Parameter","Impact Category","Percent_of_total")
subset_finalST_1$Strategy = "Mono-dim"
subset_finalST_1$Strategy = as.factor(subset_finalST_1$Strategy)

# Sampling 2
subset_finalST_2 = final_ST_unique_2[,c("value","Parameter","Impact Category","Percent_of_total")]
subset_finalST_2$Strategy = "Multi-dim"
subset_finalST_2$Strategy = as.factor(subset_finalST_2$Strategy)


final_ST_combined = rbind(subset_finalST_1,subset_finalST_2)


c29 <- c(
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
  "darkorange4", "brown","grey95","mediumpurple2","cyan3","aliceblue"
)

# To %

final_ST_combined$Percent_of_total=final_ST_combined$Percent_of_total/100

# FIGURE 5



(ggplot(final_ST_combined, aes(fill=reorder(Parameter, -Percent_of_total), y=`Impact Category`, x=Percent_of_total))
  + geom_bar(position="fill",stat="identity")
  +facet_grid(vars(Strategy))
  +theme_bw()
  +scale_fill_manual(values=c29)
  +theme(legend.title = element_text(size=16))
  +theme(legend.text = element_text(size=11))
  +theme(legend.position="bottom")
  +theme(axis.title=element_text(size=14))
  +theme(axis.text=element_text(size=14))
  +theme(strip.text.y = element_text(size = 13))
  +labs(fill="Parameter",face="bold")
  +xlab("% Total-order Sensitivy")  
  +ylab("Impact Category")
  +scale_x_continuous(labels=c("0"="0","0.25"= "25", "0.50"="50","0.75"="75","1.00"="100"))
  +guides(fill=guide_legend(nrow=7,byrow=TRUE))
)
ggsave('Figur5.jpeg', width = 13, height = 6,dpi=700)









# Supplementary Information




# With confidence intervals

ggplot(final_ST_unique_2, aes(x=reorder(Parameter, -value), y=value))+
  geom_bar(stat='identity', fill="skyblue", alpha=0.7)+
  geom_linerange( aes(x=Parameter, ymin=value-conf, ymax=value+conf), colour="orange", alpha=0.9, size=1.3)+
  facet_wrap(~`Impact Category`)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c25)+
  theme(legend.title = element_text(size=14))+
  theme(legend.text = element_text(size=11))+
  theme(legend.position="bottom")+
  theme(axis.title=element_text(size=14))+
  theme(axis.text=element_text(size=14))+
  xlab("Parameter") +
  ylab("Total-order sensitivy indice")

ggsave('FigureS23.jpeg', width = 8, height = 10,dpi=600)



# With confidence intervals


ST_with_conf_1 = cbind(ST_1,"conf" = conf_ST_1$value)


ggplot(data = ST_with_conf_1, aes(x=Parameter, y= value)) +
  geom_point(size = 1, position=position_dodge(width=0.5),color="red") +
  geom_errorbar(
    aes(ymin = value-conf, ymax = value+conf),
    width = 0.1,
    linetype = "dotted",
    position=position_dodge(width=0.5)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  geom_boxplot(outlier.shape = NA)+facet_grid(vars(`Impact Category`))+coord_cartesian(ylim = c(0, 1.5))+
  theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=14))
ggsave('FigureS22.jpeg', width = 9, height = 15,dpi=700)




