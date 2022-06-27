### Script to generate the production mixes and Figure 4.


# Set working directory to the script location.  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))



install.packages(c("ggridges", "ggplot2", "reshape2", "plyr", "dplyr", "tidyverse","plotly","GGaly","rio","ggthemes","foreach","parallel","doParallel","foreach","EnvStats","rstatix"))
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

remotes::install_github('rpkgs/gg.layers')
library(gg.layers)
library(ggplot2)

devtools::install_github("associatedpress/aptheme")

library(aptheme)


# Load data


list_skarka=c('ES','SE','IT','PT','UK','FR','EL','CY','IE','DE')


# To keep the columns' names
to_get_columns <- fread(file="../Outputs-Multi/results_table_df_pointUK_2.csv")



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

Total_results_table = data.frame(rbindlist(list_result))


colnames(Total_results_table) = colnames(to_get_columns)



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

Total_melted$Impact_Category = as.factor(Total_melted$Impact_Category)


# Production mix uncertainty


tibble_totalmelted = as_tibble(Total_melted)
# To weigh the probalibilites for 1 location ot be part of a mix. 

unique_locations_tibble = distinct_at(tibble_totalmelted, vars(lat,long,alga,`Areal productivity kg.m-2.d-1`))

# We can weigh with the ranking for 1 strain only  as the ratio between locations regarding areal productivities are the same for all strains by construction. 
unique_locations_1alga = unique_locations_tibble[unique_locations_tibble$alga==0,]



set.seed(1) # to sample the same production mixes 

# Create random production mixes

list_of_random_markets_5 = list()
list_of_random_markets_15 = list()
list_of_random_markets_25 = list()



number_market=400
for(i in 1:number_market){
  
  samp_locations_5 <- sample(unique_locations_1alga$long, 5, prob=unique_locations_1alga$`Areal productivity kg.m-2.d-1`)
  samp_locations_15 <- sample(unique_locations_1alga$long, 10, prob=unique_locations_1alga$`Areal productivity kg.m-2.d-1`)
  samp_locations_25 <- sample(unique_locations_1alga$long, 25, prob=unique_locations_1alga$`Areal productivity kg.m-2.d-1`)
  
  
  list_of_random_markets_5=append(list_of_random_markets_5,list(samp_locations_5))
  list_of_random_markets_15=append(list_of_random_markets_15,list(samp_locations_15))
  list_of_random_markets_25=append(list_of_random_markets_25,list(samp_locations_25))
}


# # Average score 
# Total_melted_average=Total_melted %>% group_by(alga,Impact_Category,long)%>%
#   mutate(mean_alga_location = mean(Score)) %>% # calculate mean for plotting as well
#   ungroup()
# # 
# # Total_melted_average$median_alga_location=Total_melted$median_alga_location
# # total_comp = melt(Total_melted_average,measure.vars = c("mean_alga_location","median_alga_location"),variable.name = 'stat',value.name='value stat')
# # unique(Total_melted_average$Impact_Category)
# # mean(Total_melted_average[Total_melted_average$Impact_Category=="GW100, kg CO2-eq",]$median_alga_location)
# # mean(Total_melted_average[Total_melted_average$Impact_Category=="GW100, kg CO2-eq",]$mean_alga_location)
# # 
# # 
# # ggplot(total_comp[total_comp$Impact_Category=="GW100, kg CO2-eq",], aes(x=`value stat`, fill=`stat`))+geom_density(alpha=0.4)+coord_cartesian(xlim = c(0, 1000))
# #                                                                               




# Aggregate techno-operational uncertainty per strain-compound per location
# one long per location, no need to take lat, long together

#median score
Total_melted=Total_melted %>% group_by(alga,Impact_Category,long)%>% 
  mutate(median_alga_location = median(Score)) %>% 
  ungroup()


Total_melted_per_location_tibble=distinct_at(as_tibble(Total_melted),vars(long,alga,Impact_Category,median_alga_location))






registerDoParallel(cores = 50)


rank_markets=1:400


number_alga=0:5375


rank_mark =rep(rank_markets,each=length(number_alga))
alga = rep(number_alga,length(rank_markets))
length(rank_mark)
length(alga)


# Convert tibble back to data.table
Total_melted_per_location_table <- data.table(Total_melted_per_location_tibble)


list_datatables=list()



# Calculate impacts for all strains in all mixes.
for(index_mix in 1:length(list_of_random_markets_5)){
  
  
  prod_mix_5=list_of_random_markets_5[[index_mix]]
  prod_mix_15=list_of_random_markets_15[[index_mix]]
  prod_mix_25=list_of_random_markets_25[[index_mix]]
  # print(prod_mix)
  
  
  
  mean_for_market_mix = Total_melted_per_location_table[long %in% prod_mix_5,  .(Impact_mix_5=mean(median_alga_location)), by=.(alga,Impact_Category)]
  
  
  mean_for_market_mix[,"Impact_mix_15"] = Total_melted_per_location_table[long %in% prod_mix_15,  .(Impact_mix_15=mean(median_alga_location)), by=.(alga,Impact_Category)][,"Impact_mix_15"]
  
  mean_for_market_mix[,"Impact_mix_25"] = Total_melted_per_location_table[long %in% prod_mix_25,  .(Impact_mix_25=mean(median_alga_location)), by=.(alga,Impact_Category)][,"Impact_mix_25"]
  
  mean_for_market_mix[,"mix"]=index_mix
  #print(mean_for_market_mix)
  #append(list_datatables[index_mix] = mean_for_market_mix
  list_datatables=append(list_datatables,list(mean_for_market_mix ))        
  
}






# Create new table per mix

per_market=rbindlist(list_datatables)


class(per_market) <- class(as.data.frame(per_market))

per_market_melt = melt(per_market,measure.vars = c("Impact_mix_5","Impact_mix_15","Impact_mix_25"),variable.name = 'Mix size',value.name='Impact_per_mix')





# FIGURE  4 Boxplot




# without outliers





levels(per_market_melt$`Mix size`)=c("5 Locations","15 Locations","25 Locations")

per_market_melt=per_market_melt %>% group_by(Impact_Category,`Mix size`)%>%
  mutate(`Mean of mixes` = mean(Impact_per_mix),`Median of mixes` = median(Impact_per_mix)) %>% # calculate mean for plotting as well
  ungroup()

per_market_melt=melt(per_market_melt,measure.vars = c("Median of mixes","Mean of mixes"),variable.name = 'Stat',value.name='Score Stat')

per_market_melt$mix = as.factor(per_market_melt$mix)




# Refine axis names 



per_market_melt$facets_mixsize <- recode_factor(per_market_melt$`Mix size`,
                                                `5 Locations` = "5~Locations",
                                                `15 Locations`="15~Locations",
                                                `25 Locations`="25~Locations")



per_market_melt$facets_Impact_Category <- recode_factor(per_market_melt$Impact_Category, `WD, m3 water` = "WD ~(m^{3}~water)",
                                                        `GW100, kg CO2-eq`="GW100~(kg~CO[2]-eq)",
                                                        `FE, kg P-eq`="FE~(kg~P-eq)",
                                                        `TETinf, kg 1.4-DC`="TETinf~(kg~1.4-DC)")









ggplot(per_market_melt, aes(x=`mix`, y=`Impact_per_mix`))+
  geom_boxplot2(color="red",alpha=0.3)+stat_summary(fun.y=mean, geom="point", shape=20, size=1, color="black", fill="black")+
  stat_summary(fun.y=median, geom="point", shape=20, size=1, color="blue", fill="blue")+
  facet_grid(cols =vars(facets_mixsize) ,rows=vars(facets_Impact_Category),scales="free",labeller = label_parsed)+geom_hline(aes(yintercept = `Score Stat`,colour=Stat))+scale_color_manual(values=c("blue", "black"))+
  xlab("Production mixes") +ylab("Score for production mix")+theme(legend.position = "none")

ggsave('Figure_4.jpeg', width = 10, height = 7,dpi=600)







# Export Statistics


tibble_per_market_melt=as_tibble(per_market_melt)

unique_stats_tibble=distinct_at(tibble_per_market_melt, vars(Stat,Impact_Category,`Mix size`,`Score Stat`))

write.csv(unique_stats_tibble,"stats Production mix.csv",row.names = TRUE)


per_market_melt_stats=per_market_melt %>%
  group_by(`Mix size`,`Impact_Category`) %>%
  get_summary_stats(Impact_per_mix, show = c("mean","median","sd"))




per_market_melt_stats=per_market_melt %>% group_by(`Mix size`,`Impact_Category`) %>%
  summarise("Geometric Mean" = psych::geometric.mean(Impact_per_mix, na.rm = T),
            "Geometric Standard deviation" = EnvStats::geoSD(Impact_per_mix, na.rm = FALSE, sqrt.unbiased = TRUE),
            "Standard deviation"= sd(Impact_per_mix, na.rm = T),
            "Mean" = mean(Impact_per_mix, na.rm = T),
            "Median"= median(Impact_per_mix, na.rm = T))


write.csv(per_market_melt_stats,"stats Production mix_extended.csv",row.names = TRUE)
