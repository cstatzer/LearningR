# Install additional libraries to be able to generate all plots (run install functions once and then you can delete them):
# install.packages("viridis")
# install.packages("ggforce")
# install.packages("ggrepel")
# install.packages("survminer")
# install.packages("survival")
# install.packages("magrittr")
# install.packages("scatterpie")
library(ggExtra)


#################################
########### R basics ############ 
#################################
library(tidyverse)
Kelvin <- 273.15
tempsC <- c(5,20,30)
feeling <- c("cold","medium","hot")

mean(tempsC)
max(tempsC)
length(feeling)

# option A
temp <- mean(tempsC)
round(temp)
# option B
tempsC %>% mean() %>% round()





#################################
######## Data workflow ########## 
#################################
world <- read_csv("./data/gapminder.csv")
world
View(world)


world %>%
  filter(continent == "Americas")


world %>%
  filter(continent == "Asia") %>%
  select(country,year,lifeExp) %>%
  arrange(desc(lifeExp))

world %>%
  group_by(continent, year) %>%
  summarize(mean_life_exp = mean(lifeExp)) %>%
  arrange(year,continent)





##################################################
######## Data visualization with ggplot ########## 
##################################################
ggplot(data = world, mapping = aes(x = lifeExp)) +
  geom_histogram()

# The LEGO principle
plot <- ggplot(data = world,mapping = aes(x = continent, y = lifeExp))
plot
plot + geom_point()
plot + geom_boxplot()
plot + geom_violin() + geom_boxplot(width = 0.1)


# Data exploration
plot <- ggplot(data = world, mapping = aes(x = gdpPercap, y = lifeExp)) +
geom_point(alpha = 0.5) +
geom_smooth() +
ggtitle("Do rich countries have a higher life expectancy?") 
plot
ggMarginal(plot,type = "boxplot")

ggplot(data = world, mapping = aes(x = gdpPercap, y = lifeExp, color = continent)) +
  geom_point() +
  geom_smooth() +
  ggtitle("Do rich countries have a higher life expectancy?") 


# ggplot(data = world, mapping = aes(x = year, y = country, fill = gdpPercap)) +
#   geom_tile()


# The Economist
library(viridis)
library(ggforce)
library(ggrepel)

world_2002 <- world %>%
  filter(year == 2002) %>%
  mutate(pop=pop/1000000)

africa_leading <- world_2002 %>% 
  filter(continent == "Africa", 
         gdpPercap > 3000)

plot <- ggplot(data = world_2002, aes(x=gdpPercap, y=lifeExp, size = pop, color = continent, label = country)) +
  geom_point(alpha=0.8) +
  scale_size(range = c(1, 20), name="Population (mio)") +
  scale_color_viridis(discrete=TRUE, guide=FALSE) +
  facet_zoom(x = continent == "Africa") +
  geom_label_repel(data = africa_leading,size = 3) +
  theme_bw() +
  labs(title = "World demographics in 2002",x = "GDP per capita [US dollars per inhabitant]", y = "Life expectancy  [years]")

plot






#################################
########### Survival ############ 
#################################

library(tidyverse)
library(survminer)
library(survival)

data <- read_csv(file = "./data/lifespan.csv")
fit<- survfit(Surv(Time, Censored) ~ Strain, data = data)
# option A
ggsurvplot(fit = fit, data = data)
# option B
ggsurvplot(fit = fit, data = data,conf.int = TRUE,risk.table = TRUE)
# option C
ggsurvplot(fit = fit, data = data,conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv")

fit_repeat<- survfit(Surv(Time, Censored) ~ Strain + Repeat, data = data)
ggsurvplot_facet(fit = fit_repeat, data = data,facet.by = "Repeat",conf.int = TRUE,risk.table = TRUE,surv.median.line = "hv")




#################################
########## Microscopy ########### 
#################################
library(tidyverse)
library(waffle)

colors <- c("grey50","bisque1", "darkorange2", "red","darkred")

data <- read_csv(file = "./data/microscopy.csv")
absolute <- data$Frequency
names(absolute) <- data$category

waffle(absolute, rows = 10, size = 1, colors = colors, legend_pos = "bottom",title = "Absolute observations")

relative <- round(absolute/sum(absolute) * 25)
waffle(relative, rows = 1, size = 1, colors = colors, legend_pos = "bottom",title = "Relative fractions")





#################################
####### Merging datasets ######## 
#################################
chip <- read_csv("./data/chip.csv")
rnaseq <- read_csv("./data/rna_seq.csv")
full_join(chip,rnaseq, by = "gene")





#################################
######## Parsing strings ######## 
#################################
library(tidyverse)
input <- read_csv("./data/aa_sequence.csv") 
input
str <- input %>% unlist() %>% paste(., collapse=" ")
str_proc <- str_replace_all(string = str,pattern = "[[:digit:],[:space:]]",replacement = "")
str_proc





#################################
#### Data cleaning RNA seq ###### 
#################################
library(tidyverse)

data <- read_csv(file = "./data/RNAseq.csv")
data <- data %>% mutate(logFC = as.numeric(logFC)) %>% filter(!is.na(logFC))
readable <- data %>% separate(col = Gene_id,into = c("Junk","ID"),sep = "_")
readable <- readable %>% filter(!is.na(ID))
ggplot(data = data, aes(x = logFC,y = -log10(padj))) +
  geom_point(alpha = 0.5)





#################################
##### 96-well plate scoring ##### 
#################################
library(tidyverse)
library(ggplot2)
library(magrittr)
library(scatterpie)

lets<- c("H","G","F","E","D","C","B","A")
df <- read_csv("./data/96well.csv")
df <- df %>% mutate(region = paste0(Row,Column), Row = match(Row, lets))
df[,4:9] <- df[,4:9] / 80

p <- ggplot() + 
  geom_point(data=expand.grid(seq(1, 12), seq(1,8)), aes(x=Var1, y=Var2),
             color="grey80", fill="white", shape=21, size=6) +
  geom_scatterpie(aes(x= Column, y= Row, r = Total, group = region), data=df, 
                  cols=c("Extremely high","High", "Intermediate", "Low", "None" ))

p + scale_x_discrete(name ="Columns", limits=c(1:12)) +
  scale_y_discrete(name ="Rows", limits=lets) +
  scale_fill_manual(values = c("red","darkgreen","darkolivegreen3","darkseagreen1","grey")) +
  coord_fixed(ratio=12/12) + facet_grid(Day ~ .) +
  ggtitle("Promotor activity (day 1 and 8")

