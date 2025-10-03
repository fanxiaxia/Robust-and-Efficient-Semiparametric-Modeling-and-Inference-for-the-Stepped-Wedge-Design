data <- read.csv("dataforpreisser.csv")
data2 <- data[,c("timex","treatx","ct","hdist")]

data2 %>% filter(hdist == 3) %>% group_by(timex) %>% summarise(meanx = mean(treatx))

library(dplyr)
# data.naomit %>% group_by(hdist,hdist) %>% summarise(n = n()) -> uniqueclusters
# length(unique(data.naomit$hdist))
# length(unique(data.naomit$hdist))

data.naomit <- data2 %>% mutate(timewave = case_when(
  timex == 0 ~ 0,
  timex > 0 & timex <=3 ~ 1,
  timex > 3 & timex <=6 ~ 2,
  timex > 6 & timex <=8 ~ 3,
  timex > 8 & timex <=14 ~ 4,
  TRUE ~ NA
))

data.use <- data.naomit[,c("hdist","timewave","treatx","ct")]
colnames(data.use) <- c("cluster","time","x","y")
#data.use %>% group_by(cluster) %>% arrange(cluster,time) %>% 
#  mutate(change_point = cumsum(lag(x, default = 0) == 0 & x == 1)) %>%
#  mutate(l = ifelse(change_point > 0, time - time[which.max(change_point)] +1, 0)) %>%
#  ungroup() -> data.use

sapply(data.use, function(x) sum(is.na(x)))

data.use <- na.omit(data.use)

data.use <- data.use %>% mutate(time = time + 1) #time period ranges from 1 to 5 instead of 0 to 4
  
data.use %>% arrange(cluster) %>% mutate(ind = 1:n()) -> data.use

data.use %>% group_by(cluster,time) %>% mutate(size = n()) -> data.use

data.use %>% group_by(cluster) %>% mutate(basesize = size[which.min(time)]) -> data.use


sampled_data <- data.use %>%
  group_by(cluster, time) %>%
  slice_sample(prop = 1, replace = FALSE) %>%
  ungroup()

dataset <- sampled_data

#calculate crossover time

df_crosstime <- dataset %>%
  group_by(cluster) %>%
  arrange(time) %>%
  mutate(x_lag = lag(x, default = 0)) %>%  # Get previous x value
  filter(x == 1 & x_lag == 0) %>%  # Keep only the first transition from 0 → 1
  summarise(cross_time = min(time), 
            basesize = mean(basesize),
            .groups = "drop")  # Get first occurrence

cor(df_crosstime$cross_time,df_crosstime$basesize)
#-0.3113959

summary_data <- dataset %>%
  group_by(cluster, time) %>%
  summarise(sample_size = n(), .groups = 'drop')

library(tidyr)
sample_matrix <- summary_data %>%
  pivot_wider(names_from = time, values_from = sample_size, values_fill = 0)

diffsize <- as.data.frame(sample_matrix[,-1])
