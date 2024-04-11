### Script for the analyses of the manuscript entitled
### "Manipulating the strength of organism-environment feedback increases nonlinearity and apparent hysteresis of ecosystem response to environmental change"
### Authors: Aurelie Garnier, Florence D. Hulot and Owen L. Petchey

#####################
##### Libraries #####
#####################

library(extrafont) 
library(plyr)
library(dplyr)
library(tidyr)

library(rootSolve)
library(numDeriv)
library(deSolve)
library(vegan)
library(caret)
library(dtwclust)

library(ggplot2)
library(ggfortify)
library(grid)
library(cowplot)
library(gridExtra)
library(ggdendro)
library(gtable)


#### ACTIONS: 
#### 1) change accordingly where the data and script are stored
#### 2) set this directory as the working directory 
data_location <- "~/Documents/GitHub/OEF-strength/data/"
setwd("~/Documents/GitHub/OEF-strength/data/")


####################################
#####                          #####
#####  Minimal Ecosystem Model #####
#####                          #####
####################################

# A minimal model of an ecosystem showing hysteresis

# $\frac{dx}{dt} = a - bx + rf(x)$

# PARAMETERS
# $x$ is an "unwanted" ecosystem property.
# $a$ is an environmental factor that "promotes" $x$.
# $b$ is the rate of decay of $x$.
# $r$ is the rate at which $x$ recovers, as a function of $x$.

# "For a lake, one can think of $x$ as nutrients suspended in phytoplankton causing turbidity, 
# of $a$ as nutrient loading, and $b$ as nutrient removal rate, 
# and of $r$ as internal nutrient recycling. For desertification, one could interpret $x$ as barren soil, 
# $a$ as vegetation destruction, $b$ as recolonization of barren soil by plants and $r$ as erosion by wind and runoff."

# set some global parameter values
xs <- seq(0, 2, 0.01)
a = 1
as <- seq(0, 1, 0.001)
b = 1
r = 1


# Hysteresis and alternative stable states can occur if $f(x)$ is a function with a threshold, 
# e.g. the hill function:
#  $f(x) = \frac{x^p}{x^p + h^p}$

# $h$ is the hill coefficient set to 1
h = 1
# $p$ is the self feedback strength
p = 1
ps <- c(seq(1, 5, by=1), 10, 20, 50)


# the functions that implement this model
dx_dt <- function(x, a, b, r, p, h) {
  a - b*x + r*f_of_x(x, p, h)
}

f_of_x <- function(x, p, h) {
  x^p / (x^p + h^p)
}

dx_dt2 <- function(x, p=1) {
  a - b*x + r * x^p / (x^p + h^p)
}


###Explore f(x), the function that determines the internal feedback

# plot f_of_x vs x:
simul1 <- expand.grid(ps=ps, xs=xs)

results1 <- mutate(simul1, f_of_x=f_of_x(x=xs, p=ps, h=h))

qplot(xs, f_of_x, data=results1, col=as.factor(ps), geom="path") +
  xlab("System state (X)") +
  ylab("Value of self feedback (f(X))") +
  scale_colour_discrete(guide = guide_legend(title = "Nonlinearity of self feedback (p)"))


# set the range of values of a, for the x-axis of plot

simul2 <- expand.grid(as=as, ps=ps)
results2 <- simul2 %>%
  group_by(as, ps) %>%
  do(roots = as.data.frame(uniroot.all(dx_dt, interval=c(0,10),
                                       a=.$as, b=b, r=r, p=.$ps, h=h))) %>%
  tidyr::unnest()
names(results2) <- c("as", "ps", "roots")
results2 <- arrange(results2, ps, roots)
results2$stability <- ifelse(grad(dx_dt, x=results2$roots, a=results2$as, b=b, r=r, p=results2$ps, h=h)<0, "Stable", "Unstable")

col_p <- c("grey50", "#7AD151FF", "#22A884FF", "#2A788EFF", "#414487FF", "#440154FF")
ggplot(data=subset(results2, stability=="Stable" & roots<1 & ps<15), aes(x=as, y=roots, color=as.factor(ps)))  + geom_line(size=1) + 
  geom_line(data=subset(results2, stability=="Stable" & roots>1 & ps<15), aes(x=as, y=roots, color=as.factor(ps)),size=1) +
  geom_line(data=subset(results2, stability=="Unstable" & ps<15 ), aes(x=as, y=roots, color=as.factor(ps)), linetype="dashed",size=1) +
  labs(x=expression(paste("Environmental condition (",italic("a"),")")), y=expression(paste("Equilibrium state (",italic("Y"),")")) , title="(a)") +
  scale_colour_manual(values=col_p, guide = guide_legend(title = expression(paste("Stregnth of \n self feedback (",italic("p"),")")))) +
  ylim(0,2)

###############
#### MEASURING AMOUNT OF HYSTERESIS

## get non-linearity and hysteresis for a value of p
get_nl_hyst_by_p <- function(p, inv_rate=20000) {
  
  this.p <- p
  parameters <- c(a = a, b = b, r = r, p = this.p, h = h)
  state <- c(x = 0.1)
  
  model <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dx <- dx_dt(x, a_forcing2(t), b, r, p, h)
      list(c(dx), a=a_forcing2(t))
    })
  }
  
  times <- seq(0, inv_rate, length = 2001)
  a_forcing1 <- matrix(ncol=2, byrow=T, data=c(0,0, mean(times), 1, max(times), 0))
  a_forcing2 <- approxfun(x = a_forcing1[,1], y = a_forcing1[,2], method = "linear", rule = 2)
  out <- as.data.frame(ode(y = state, times = times, func = model, parms = parameters))
  
  ## nonlinearity
  get_nl <- function(x, y)
  {
    lin_pred <- predict(lm(y ~ x))
    gam_pred <- predict(mgcv::gam(y ~ s(x)))
    L <- sqrt(sum((lin_pred-gam_pred)^2)/length(x))
    L
  }
  nl_up <- get_nl(out$a[1:1000], out$x[1:1000])
  nl_down <- get_nl(out$a[1001:2000], out$x[1001:2000])
  
  ## hysteresis
  out2 <- data.frame(up_x=out$x[1:1000],
                     down_x=out$x[2000:1001])
  hyst = mean(abs(out2$up_x - out2$down_x))
  
  ## hysteresis between 0.25 and 0.75
  out3 <- data.frame(up_x=out$x[251:750],
                     down_x=out$x[1750:1251])
  hyst2 = mean(abs(out3$up_x - out3$down_x))
  
  return(list(nl_up=nl_up, nl_down=nl_down, hyst=hyst, hyst2=hyst2))
}





ps <- 10^seq(0, log10(10), length=20)
inv_rates <- c(200, 20000)

simul3 <- expand.grid(ps=ps,
                     inv_rates=inv_rates)

temp_res3 <- simul3 %>%
  group_by(ps, inv_rates) %>%
  do(rez = as.data.frame(get_nl_hyst_by_p(p=.$ps, inv_rate = .$inv_rates)))

results3 <- temp_res3 %>%
  tidyr::unnest() %>%
  gather(key=Variable, value=Value, 3:6) %>%
  mutate(Variable=case_when(Variable=="nl_up" ~ "Non-linearity (up phase)",
                            Variable=="nl_down" ~ "Non-linearity (down phase)",
                            Variable=="hyst" ~ "Hysteresis",
                            Variable=="hyst2" ~ "Hysteresis2"),
         `Simulation duration`=as.character(inv_rates))

results3 <- as.data.frame(results3)



### Comparison between the method mentioned by the reviewer (selection of environmental conditions to measure hysteresis) and our method (no selection of specific environmental conditions)

results4 <- temp_res3 %>%
  tidyr::unnest() %>%
  gather(key=Variable, value=Value, 3:6)
results5 <- subset(results4, Variable == "hyst" | Variable =="hyst2")
results5$inv_rates <- ifelse(results5$inv_rates==200, paste0(200," (fast rate)"), paste0(20000," (slow rate)"))

#png("Compare_hysteresis.png", res=200, width=1200, height=800)
ggplot(results5, aes(x=ps, y=Value, fill=Variable, colour=Variable, shape=inv_rates))+
  geom_line(colour="black") + geom_point(size=3) +
  scale_fill_manual(values=c("black","blue"), 
                    labels=c("Figure 2 \n- overall environmental conditions [0-1]",
                             "Comment [R1.5] from first reviews\n- selected environmental conditions [0.25-0.75]"), 
                    name="Method")+
  scale_colour_manual(values=c("black","blue"),
                      labels=c("Figure 2 \n- overall environmental conditions [0-1]",
                               "Comment [R1.5] from first reviews\n- selected environmental conditions [0.25-0.75]"), 
                      name="Method")+
  scale_shape_manual(values=c(21,22), 
                     name="Simulation duration")+
  labs(x=expression(paste("Strength of self feedback (",italic("p"),")")), y="", title="Hysteresis (absolute mean difference)")+ 
  guides(colour=guide_legend(order=1),
         shape=guide_legend(order=2),
         fill="none")+
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(), 
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key=element_blank())
#dev.off()


#######################################################################
#####                                                             #####
##### Analyses of the Experimental Feedback Strength Manipulation #####
#####                                                             #####
#######################################################################

# experimental design
expt_design <- data.frame(microcosm=seq(101,130,by=1),
                          gas=c(rep("Open",6), rep("Closed",6), rep("+ 50 mL",6), rep("+ 100 mL",6), rep("+ 200 mL",6)),
                          replicate=rep(seq(1,6,by=1),5))

# import the data
full_data <- read.csv(file=paste0(data_location,"Data_FeedbackStrengthManipulation.csv"), header=T)

# correct and order the levels of the gas treatment
full_data$gas <- gsub("200mL","+ 200 mL", gsub("100mL","+ 100 mL", gsub("50mL","+ 50 mL", full_data$gas)))
gas_level_ordered <- c("Open","+ 200 mL","+ 100 mL", "+ 50 mL", "Closed")
full_data$gas <- factor(full_data$gas, levels=gas_level_ordered, ordered = T)

# add a variable $name_microcosm for the figures
full_data$name_microcosm <- paste0(full_data$gas," #",full_data$replicate)

# log10 transformation
full_data$PRED <- log10(full_data$Density_spathidium + 1)
full_data$PREY <- log10(full_data$Density_prey + 1)
full_data$BACT <- log10(full_data$Density_bacteria + 1)
full_data$BIOM <- log10(full_data$Total_biomass + 1)

#####################
## PCA on time series

# scale and center some variables

# variables to transform
var <- c("Density_prey", "Density_spathidium", "Density_bacteria", "Total_biomass", "PREY", "PRED", "BACT","BIOM", "Oxygen_L")
t_var <- paste0("t_",c("Density_prey", "Density_spathidium", "Density_bacteria", "Total_biomass", "PREY", "PRED", "BACT","BIOM", "Oxygen_L"))

sub.data <- full_data[,c("microcosm","gas","replicate","day","Temperature_L",var)]


trans <- caret::preProcess(as.data.frame(full_data[,var]), method=c("center", "scale"))
trans_data1 <- predict(trans, as.data.frame(full_data[,var]))
names(trans_data1) <- paste0("t_",names(trans_data1))
trans_data <- cbind(sub.data, trans_data1)

# PCA on transformed variables
select_t.var <- c("t_PRED", "t_PREY", "t_BACT", "t_Total_biomass", "t_Oxygen_L")

trans_dd2 <- na.omit(trans_data)  
PCA.vegan1 <- rda(trans_dd2[,select_t.var])
PCA.vegan.eigen <- PCA.vegan1$CA$eig
#summary(PCA.vegan1)

# variance explained by the first two axes of the PCA
var.expl.PC1 <- paste0("PC1 (", (round(summary(PCA.vegan1)$cont$importance[2,1], digits=3))*100, "% variance explained)")
var.expl.PC2 <- paste0("PC2 (", (round(summary(PCA.vegan1)$cont$importance[2,2], digits=3))*100, "% variance explained)")

### PCA scores for the variables
variable.scores_a <- as.data.frame(summary(PCA.vegan1)$species)
variable.scores_a$variable <- c("Predator","Prey","Bacteria","Total biomass","DO")

### PCA scores over time
time.scores_a<- summary(PCA.vegan1)$sites
time.scores <- cbind(trans_dd2[,c("microcosm","gas","day")],time.scores_a)
time.scores$up_or_down <- ifelse(time.scores$day<34, "increasing", ifelse(time.scores$day>34, "decreasing", NA))
time.scores$gas <- factor(time.scores$gas, levels=gas_level_ordered, ordered = T)
time.scores$up_or_down <- factor(time.scores$up_or_down, levels=c("increasing","decreasing"))


data_day1 <- data.frame(microcosm=0, day=0, PC1_t=0, PC2_t=0)

time.scores2 <- time.scores
for (i in unique(time.scores2$microcosm)) {
  #i=101
  df <- subset(time.scores2, microcosm==i)
  df <- df[order(df$day),]
  
  PC1 <- c(df[2:nrow(df),"PC1"],NA)
  PC2 <- c(df[2:nrow(df),"PC2"],NA)
  
  add_df <- data.frame(microcosm=i, day=df$day, PC1_t=PC1, PC2_t=PC2)
  data_day1 <- rbind(data_day1, add_df)
}

time.scores2 <- merge(time.scores, data_day1[-1,],by=c("microcosm","day"))
time.scores3 <- merge(expt_design, time.scores2, by=c("microcosm"))



###################################################
## Hierarchical clustering for the dissolved oxygen

# transform the data into list
list_oxygen <- NULL
mm <- unique(trans_dd2$microcosm)
for(i_m in 1:length(mm)){
  # i_m = 1 
  df <- subset(trans_dd2, microcosm==mm[i_m])
  df <- df[order(df$day),]
  
  list_oxygen[[i_m]] <- df[,c("t_Oxygen_L")]
}

hc_oxygen <- tsclust(list_oxygen, type = "hierarchical", k = 2L,
                     distance = "dtw", 
                     control = hierarchical_control(method = "ward.D"),
                     preproc = NULL, # no pre-processing because of the direct use of transformed (z-score: centered and scaled) data (i.e. trans_dd2)
                     seed = 899)

# By default, the dendrogram is plotted in hierarchical clustering
plot(hc_oxygen)

# preparation for the figure 
dd.row <- as.dendrogram(hc_oxygen)
ddata_oxygen <- dendro_data(dd.row)
labs <- label(ddata_oxygen)
labs$label <- as.numeric(as.character(labs$label)) + 100
labs <- merge(labs, unique(full_data[,c("microcosm","name_microcosm","gas")]), by.x="label" , by.y="microcosm")
labs$gas <- factor(labs$gas, levels = gas_level_ordered, ordered = T)
labs_oxygen <- labs

# get the id of microcosm belonging to cluster
clusterCut_oxygen <- data.frame(cluster=cutree(hc_oxygen, 2))
clusterCut_oxygen$mm <- as.numeric(row.names(clusterCut_oxygen))+100

clusterCut_oxygen2 <- merge(expt_design, clusterCut_oxygen, by.x="microcosm", by.y="mm")
clusterCut_oxygen2$cluster <- as.factor(clusterCut_oxygen2$cluster - 1)
clusterCut_oxygen2$gas <- factor(clusterCut_oxygen2$gas, levels=gas_level_ordered, ordered = T)

### Test if there is a match between cluster and gas treatment
mod_cluster_gas_oxygen <- glm(cluster ~ gas, clusterCut_oxygen2, family=binomial)
anova(mod_cluster_gas_oxygen, test="LRT")

################################################
## Hierarchical clustering for the total biomass

# transform the data into list for TSclust
list_biomass <- NULL
mm <- unique(trans_dd2$microcosm)
for(i_m in 1:length(mm)){
  # i_m = 1 
  df <- subset(trans_dd2, microcosm==mm[i_m])
  df <- df[order(df$day),]
  
  list_biomass[[i_m]] <- df[,c("t_Total_biomass")]
}

hc_biomass <- tsclust(list_biomass, 
                      type = "hierarchical", k = 2L,
                      distance = "dtw", 
                      control = hierarchical_control(method = "ward.D"),
                      preproc = NULL, # no pre-processing because of the direct use of transformed (z-score: centered and scaled) data (i.e. trans_dd2)
                      seed = 899)
# By default, the dendrogram is plotted in hierarchical clustering
plot(hc_biomass)

# preparation for the figure
dd.row <- as.dendrogram(hc_biomass)
ddata_biomass <- dendro_data(dd.row)
labs <- label(ddata_biomass)
labs$label <- as.numeric(as.character(labs$label)) + 100
labs <- merge(labs, unique(full_data[,c("microcosm","name_microcosm","gas")]), by.x="label" , by.y="microcosm")
labs$gas <- factor(labs$gas, levels = gas_level_ordered, ordered = T)
labs_biomass <- labs

# get the id of microcosm belonging to cluster
clusterCut_biomass <- data.frame(cluster=cutree(hc_biomass, 2))
clusterCut_biomass$mm <- as.numeric(row.names(clusterCut_biomass))+100

clusterCut_biomass2 <- merge(expt_design, clusterCut_biomass, by.x="microcosm", by.y="mm")
clusterCut_biomass2$cluster <- as.factor(clusterCut_biomass2$cluster - 1)
clusterCut_biomass2$gas <- factor(clusterCut_biomass2$gas, levels=gas_level_ordered, ordered = T)

#### Test if there is a match between cluster and gas treatment
mod_cluster_gas_biomass <- glm(cluster ~ gas, clusterCut_biomass2, family=binomial)
anova(mod_cluster_gas_biomass, test="LRT")

# get the three clusters (A,B,C) combining the two ts clustering
clusterCut_biomass3 <- clusterCut_biomass
names(clusterCut_biomass3) <- c("cluster_biomass", "mm")

clusterCut_oxygen3 <- clusterCut_oxygen
names(clusterCut_oxygen3) <- c("cluster_oxygen", "mm")

clusterCut <- merge(clusterCut_biomass3, clusterCut_oxygen3, by="mm")
clusterCut$all <- ifelse(clusterCut$cluster_biomass==1 & clusterCut$cluster_oxygen==1, "A",
                         ifelse(clusterCut$cluster_biomass==1 & clusterCut$cluster_oxygen==2,"B","C"))

time.scores4 <- merge(time.scores3, clusterCut, by.x="microcosm",by.y="mm")
names(time.scores4)[which(names(time.scores4)=="all")] <- "overall_cluster"

# add the information about the cluster (A, B or C) for the hierarchical clustering
labs_oxygen <- merge(labs_oxygen, clusterCut[,c("mm","all")], by.x="label",by.y="mm")
names(labs_oxygen)[which(names(labs_oxygen)=="all")] <- "overall_cluster"

labs_biomass <- merge(labs_biomass, clusterCut[,c("mm","all")], by.x="label",by.y="mm")
names(labs_biomass)[which(names(labs_biomass)=="all")] <- "overall_cluster"

## Test in response to a Reviewer's comment [R2.XX] "the addition of air - via bubbling - might affect the nutrient release" which would affect primarily on bacteria
mean_bacteria <- ddply(full_data, .(microcosm, gas),
                       summarise,
                       mean = mean(Density_bacteria, na.rm =T))
mean_bacteria$gas <- as.factor(mean_bacteria$gas)

mod_bacteria <- lm(mean ~ gas, data=mean_bacteria)
anova(mod_bacteria)




#################################
## Nonlinearity and hysteresis 


# Transform the data
## different day range than for hysteresis analysis
hh <- dplyr::mutate(full_data, up_or_down=ifelse(day<34, "increasing", ifelse(day>34, "decreasing", NA)),
                    Temperature=round(Temperature_L, 0)) %>%
  dplyr::select(microcosm, gas, Oxygen_L, Temperature, Temperature_L, Total_biomass, up_or_down, replicate) %>%
  na.omit()


####################################
## Nonlinearity for Dissolved Oxygen

# function
get_nl <- function(x, y)
{
  lin_pred <- predict(lm(y ~ x))
  gam_pred <- predict(mgcv::gam(y ~ s(x)))
  L <- sqrt(sum((lin_pred-gam_pred)^2)/length(x))
  L
}

nl_oxygen <- dplyr::group_by(hh, microcosm, up_or_down, gas) %>%
  dplyr::summarise(get_nl=get_nl(Temperature_L, Oxygen_L))
nl_oxygen$gas <- factor(nl_oxygen$gas, levels = gas_level_ordered, ordered = T)
nl_oxygen$up_or_down <- factor(nl_oxygen$up_or_down, levels = sort(unique(nl_oxygen$up_or_down))[c(2,1)])

# merge with the cluster of coefficients
nl_oxygen <- merge(nl_oxygen, clusterCut_oxygen, by.x="microcosm", by.y="mm")
nl_oxygen$cluster <- as.factor(nl_oxygen$cluster -1 )

### Analyse the **up** phase with ANOVA

mod_nl_oxygen_up <- lm(get_nl ~ gas, filter(nl_oxygen, up_or_down=="increasing"))
autoplot(mod_nl_oxygen_up)
anova(mod_nl_oxygen_up)

#### + get_nl ~ gas * cluster
# with interaction
mod_nl_oxygen_up_int <- lm(get_nl ~ gas * cluster, filter(nl_oxygen, up_or_down=="increasing"))
# without interaction
mod_nl_oxygen_up_NOint <- lm(get_nl ~ gas + cluster, filter(nl_oxygen, up_or_down=="increasing"))

print(anova(mod_nl_oxygen_up_int,mod_nl_oxygen_up_NOint)) ## there is no difference between the 2 models - let's choose the simplest without interaction
anova(mod_nl_oxygen_up_NOint)


### Analyse the **down** phase with ANOVA

mod_nl_oxygen_down <- lm(get_nl ~ gas, filter(nl_oxygen, up_or_down=="decreasing"))
autoplot(mod_nl_oxygen_down)
anova(mod_nl_oxygen_down)

#### + get_nl ~ gas * cluster
# with interaction
mod_nl_oxygen_down_int <- lm(get_nl ~ gas * cluster, filter(nl_oxygen, up_or_down=="decreasing"))
# without interaction
mod_nl_oxygen_down_NOint <- lm(get_nl ~ gas + cluster, filter(nl_oxygen, up_or_down=="decreasing"))

print(anova(mod_nl_oxygen_down_int,mod_nl_oxygen_down_NOint)) ## there is no difference between the 2 models - let's choose the simplest without interaction
anova(mod_nl_oxygen_down_NOint)


##################################
## Hysteresis for Dissolved oxygen

# Match oxygen levels when increasing and decreasing temperature
hh_oxygen_up <- dplyr::filter(hh, up_or_down=="increasing") %>%
  dplyr::rename(oxy_L_up=Oxygen_L) %>%
  dplyr::mutate(Temperature=round(Temperature_L, 0)) 
hh_oxygen_up <- hh_oxygen_up[,-which(names(hh_oxygen_up) %in% c("up_or_down","Temperature_L","replicate","gas","Total_biomass"))]

hh_oxygen_down <- dplyr::filter(hh, up_or_down=="decreasing") %>%
  dplyr::rename(oxy_L_down=Oxygen_L) %>%
  dplyr::mutate(Temperature=round(Temperature_L, 0)) 
hh_oxygen_down <- hh_oxygen_down[,-which(names(hh_oxygen_down) %in% c("up_or_down","Temperature_L","replicate","gas","Total_biomass"))]


hh_oxygen_all <- full_join(hh_oxygen_up, hh_oxygen_down)


## Measure hysteresis by difference

hh_diff2 <- hh_oxygen_all %>% mutate(diff=abs(oxy_L_up - oxy_L_down)) 
hh_oxygen_diff <- ddply(hh_diff2, .(microcosm), summarise, mean_diff=mean(diff, na.rm=T))
hh_oxygen_diff <- merge(expt_design, hh_oxygen_diff, by="microcosm")


hh_oxygen_diff <- merge(hh_oxygen_diff, clusterCut_oxygen, by.x="microcosm", by.y="mm")
hh_oxygen_diff$gas <- factor(hh_oxygen_diff$gas, levels=gas_level_ordered, ordered = T)


### Analyse with ANOVA:
mod_hh_oxygen <- lm(mean_diff ~ gas, hh_oxygen_diff)
autoplot(mod_hh_oxygen)
anova(mod_hh_oxygen)

#### + get_nl ~ gas * cluster
hh_oxygen_diff$cluster <- as.factor(hh_oxygen_diff$cluster - 1)
# with interaction
mod_hh_oxygen_int <- lm(mean_diff ~ gas * cluster, hh_oxygen_diff)
# without interaction
mod_hh_oxygen_NOint <- lm(mean_diff ~ gas + cluster, hh_oxygen_diff)

print(anova(mod_hh_oxygen_int,mod_hh_oxygen_NOint)) ## there is no difference between the 2 models - let's choose the simplest without interaction
anova(mod_hh_oxygen_NOint)



#################################
## Nonlinearity for Total biomass

nl_biomass <- dplyr::group_by(hh, microcosm, up_or_down, gas) %>%
  dplyr::summarise(get_nl=get_nl(Temperature_L, Total_biomass))
nl_biomass$gas <- factor(nl_biomass$gas, levels = gas_level_ordered, ordered = T)
nl_biomass$up_or_down <- factor(nl_biomass$up_or_down, levels = sort(unique(nl_biomass$up_or_down))[c(2,1)])

# merge with the cluster of coefficients
nl_biomass <- merge(nl_biomass, clusterCut_biomass, by.x="microcosm", by.y="mm")
nl_biomass$cluster <- as.factor(nl_biomass$cluster -1 )

### Analyse the up phase with ANOVA
mod_nl_biomass_up <- lm(get_nl ~ gas, filter(nl_biomass, up_or_down=="increasing"))
autoplot(mod_nl_biomass_up)
anova(mod_nl_biomass_up)

#### + get_nl ~ gas * cluster
# with interaction
mod_nl_biomass_up_int <- lm(get_nl ~ gas * cluster, filter(nl_biomass, up_or_down=="increasing"))
# without interaction
mod_nl_biomass_up_NOint <- lm(get_nl ~  gas + cluster, filter(nl_biomass, up_or_down=="increasing"))

print(anova(mod_nl_biomass_up_int,mod_nl_biomass_up_NOint)) ## there is no difference between the 2 models - let's choose the simplest without interaction
anova(mod_nl_biomass_up_NOint)

### Analyse the down phase with ANOVA

mod_nl_biomass_down <- lm(get_nl ~ gas, filter(nl_biomass, up_or_down=="decreasing"))
autoplot(mod_nl_biomass_down)
anova(mod_nl_biomass_down)

#### + get_nl ~ gas * cluster
# with interaction
mod_nl_biomass_down_int <- lm(get_nl ~ gas * cluster, filter(nl_biomass, up_or_down=="decreasing"))
# without interaction
mod_nl_biomass_down_NOint <- lm(get_nl ~ gas + cluster, filter(nl_biomass, up_or_down=="decreasing"))

print(anova(mod_nl_biomass_down_int, mod_nl_biomass_down_NOint)) ## there is a difference between the 2 models - let's choose the model with interaction
anova(mod_nl_biomass_down_int)


###############################
## Hysteresis for Total biomass

# Match oxygen levels when increasing and decreasing temperature
hh_biomass_up <- dplyr::filter(hh, up_or_down=="increasing") %>%
  dplyr::rename(biom_up=Total_biomass) %>%
  dplyr::mutate(Temperature=round(Temperature_L, 0)) 
hh_biomass_up <- hh_biomass_up[,-which(names(hh_biomass_up) %in% c("up_or_down","Temperature_L","replicate","gas","Oxygen_L"))]


hh_biomass_down <- dplyr::filter(hh, up_or_down=="decreasing") %>%
  dplyr::rename(biom_down=Total_biomass) %>%
  dplyr::mutate(Temperature=round(Temperature_L, 0)) 
hh_biomass_down <- hh_biomass_down[,-which(names(hh_biomass_down) %in% c("up_or_down","Temperature_L","replicate","gas","Oxygen_L"))]


hh_biomass_all <- full_join(hh_biomass_up, hh_biomass_down)

## Measure hysteresis by difference
hh_diff2 <- hh_biomass_all %>% mutate(diff=abs(biom_up - biom_down)) 

hh_biomass_diff <- ddply(hh_diff2, .(microcosm), summarise, mean_diff=mean(diff, na.rm=T))
hh_biomass_diff <- merge(expt_design, hh_biomass_diff, by="microcosm")

hh_biomass_diff <- merge(hh_biomass_diff, clusterCut_biomass, by.x="microcosm", by.y="mm")
hh_biomass_diff$gas <- factor(hh_biomass_diff$gas, levels=gas_level_ordered, ordered = T)

### Analyse with ANOVA:
mod_hh_biomass <- lm(mean_diff ~ gas, hh_biomass_diff)
autoplot(mod_hh_biomass)
anova(mod_hh_biomass)

#### + get_nl ~ gas * cluster
hh_biomass_diff$cluster <- as.factor(hh_biomass_diff$cluster - 1)
# with interaction
mod_hh_biomass_int <- lm(mean_diff ~ gas * cluster, hh_biomass_diff)
# without interaction
mod_hh_biomass_NOint <- lm(mean_diff ~ gas + cluster, hh_biomass_diff)

print(anova(mod_hh_biomass_int,mod_hh_biomass_NOint)) ## there is no difference between the 2 models - let's choose the simplest without interaction
# summary(mod_hh_biomass_NOint)
anova(mod_hh_biomass_NOint)




############################################### REPRODUCTION OF  THE FIGURES
color_gas <- c("gold","lightgreen","lightblue","royalblue3","darkorchid4")
treatment <- gas_level_ordered


##################################
##### Preliminary experiment #####
##################################

preli_expt <- read.csv(file=paste0(data_location,"Data_Preliminary_experiment.csv"), header=T, sep = ",")
preli_expt$design <- factor(preli_expt$design, levels = c("Open","Closed","+ 25mL","+ 50mL","+ 100mL","+ 200mL","+ 400mL"))

preli_expt_plot <- ggplot(preli_expt, aes(x=deltaDay,y=Value, group=SensorName, fill=Mode, shape=Mode)) +
  geom_line(aes(linetype=Mode, color=Mode)) + geom_point() + 
  facet_grid(~ design) +
  scale_fill_manual(name="Phase",breaks=c("dry","humid"), labels=c("head","liquid"),values = c("white","black"))+
  scale_shape_manual(name="Phase",breaks=c("dry","humid"), labels=c("head","liquid"),values = c(22,21))+
  scale_color_manual(name="Phase",breaks=c("dry","humid"), labels=c("head","liquid"),values = c("grey50","black"))+
  scale_linetype_manual(name="Phase",breaks=c("dry","humid"), labels=c("head","liquid"),values = c(2,1))+
  labs(x="Time (days)", y= "Oxygen concentration (%O2)")+
  ylim(0,20)+
  theme_bw() +
  theme(text = element_text(family="Times"),
        panel.grid = element_blank(),
        strip.background = element_rect(fill="white", color="white"))

tiff("Figure2_preli_experiment.tiff", res=600, height = 1500, width = 3000)
preli_expt_plot
dev.off()


#############################################
### Reproduction of the figure 3 - Simulation

col_p <- c("grey50", "#7AD151FF", "#22A884FF", "#2A788EFF", "#414487FF", "#440154FF")
scheffer <- ggplot(data=subset(results2, stability=="Stable" & roots<1 & ps<15), aes(x=as, y=roots, color=as.factor(ps)))  + geom_line(size=1) + 
  geom_line(data=subset(results2, stability=="Stable" & roots>1 & ps<15), aes(x=as, y=roots, color=as.factor(ps)),size=1) +
  geom_line(data=subset(results2, stability=="Unstable" & ps<15 ), aes(x=as, y=roots, color=as.factor(ps)), linetype="dashed",size=1) +
  labs(x=expression(paste("Environmental condition (",italic("a"),")")), y=expression(paste("Equilibrium state (",italic("Y"),")")) , title="(a)") +
  scale_colour_manual(values=col_p, guide = guide_legend(title = expression(paste("Strength of \nself-feedback (",italic("p"),")")))) +
  ylim(0,2)+ 
  theme(axis.line = element_line(colour = "black"),
        legend.background = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key=element_blank(),
        text = element_text(size=15, family="Times"))

results3$`Simulation duration` <- ifelse(results3$inv_rates==200, paste0(200," (fast rate)"), paste0(20000," (slow rate)"))

NL_up <- ggplot(subset(results3, Variable=="Non-linearity (up phase)"), aes(x=ps, y=Value, shape=`Simulation duration`, fill=`Simulation duration`)) +
  geom_line() + geom_point(size=3) +
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=c("white","black"))+
  labs(x="", y="Value", title="(b) Nonlinearity", subtitle = "(up phase)")+ 
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position="none",
        text = element_text(size=13, family="Times"))


NL_down <- ggplot(subset(results3, Variable=="Non-linearity (down phase)"), aes(x=ps, y=Value, shape=`Simulation duration`, fill=`Simulation duration`)) +
  geom_line() + geom_point(size=3) +
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=c("white","black"))+
  labs(x=expression(paste("Strength of self-feedback (",italic("p"),")")), y="", title="(c) Nonlinearity", subtitle = "(down phase)")+ 
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position="none",
        axis.title.x = element_text(hjust=-0.1),
        text = element_text(size=13, family="Times"))

#expression(paste("Strength of self-feedback (",italic("p"),")"))

Hyst <- ggplot(subset(results3, Variable=="Hysteresis"), aes(x=ps, y=Value, shape=`Simulation duration`, fill=`Simulation duration`)) +
  geom_line() + geom_point(size=3) +
  scale_shape_manual(values=c(21,22))+
  scale_fill_manual(values=c("white","black"))+
  labs(x="", y="", title="(d) Hysteresis", subtitle = "(absolute mean difference)")+ 
  theme(axis.line = element_line(colour = "black"),
        panel.background = element_blank(), 
        legend.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size=13, family="Times"))

measure_simul <- ggpubr::ggarrange(NL_up, NL_down, Hyst, nrow=1, common.legend = TRUE, legend = "right")

tiff("Figure3_simulation.tiff", res=300, height = 3000, width = 3000)
gridExtra::grid.arrange(scheffer, measure_simul, 
                        ncol=1, heights=c(2,1))
dev.off()


###########################################
### Reproduction of the figure 4 - Dynamics


## function to change the color background of strip according to a conditional value (i.e. overall cluster)
## solution found at https://stackoverflow.com/questions/19440069/ggplot2-facet-wrap-strip-color-based-on-variable-in-data-set
gtable_select <- function (x, ...) 
{
  matches <- c(...)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  x
}

### add the information of the cluster in a new data set
full_data_cl <- merge(full_data, labs_oxygen[,c("label","overall_cluster")], by.x="microcosm", by.y="label")
color_cluster <- c("A"="#3973A4", "B"="#5E3AAD", "C"="#923638")

fig3_list <- list()

for(ttt in 1:length(treatment)){
  # ttt=1
  
   plot1 <- ggplot(data=subset(full_data_cl, gas==treatment[ttt]), aes(x=day, y=Oxygen_H) ) +
    geom_point(colour="black", alpha=.5) + 
    geom_line(colour="black")+ 
    ylim(0,20)+
    geom_point(aes(x=day, y=Oxygen_L), colour="royalblue", alpha=.5) + 
    stat_smooth(inherit.aes = F, data=subset(full_data_cl, gas==treatment[ttt]), aes(x=day, y=Oxygen_L), method = "loess", size = 1.5, se=F, colour="royalblue")+ # method = "lm", formula = y ~ x + I(x^2)+I(x^3)
    facet_wrap(~name_microcosm, ncol=6)+
    theme_bw() + labs(x="", y="% O2", title=treatment[ttt]) +
    theme(axis.text.x = element_blank(), panel.grid = element_blank(), strip.background = element_blank(),
          plot.margin = unit(c(0,0.1,0,0.38), "cm"), 
          text = element_text(size=30, family="Times"), 
          plot.title = element_text(size=30, face="bold")) 
   
   dummy <- ggplot(data=subset(full_data_cl, gas==treatment[ttt]), aes(x=day, y=Oxygen_H)) + 
     facet_wrap(~name_microcosm, ncol=6)+
     scale_fill_manual(values=color_cluster)+
     geom_rect(aes(fill=overall_cluster), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+ theme_minimal()+
     theme(strip.text.x = element_text(size=30, color = "white", family="Times"))
   
   g1 <- ggplotGrob(plot1)
   g2 <- ggplotGrob(dummy)
   panels <- grepl(pattern="panel", g2$layout$name)
   strips <- grepl(pattern="strip-t", g2$layout$name)
   g2$layout$t[panels] <- g2$layout$t[panels] - 1
   g2$layout$b[panels] <- g2$layout$b[panels] - 1
   new_strips <- gtable_select(g2, panels | strips)
   # grid.newpage()
   # grid.draw(new_strips)
   gtable_stack <- function(g1, g2){
     g1$grobs <- c(g1$grobs, g2$grobs)
     g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
     g1$layout <- rbind(g1$layout, g2$layout)
     g1
   }
   
  
  fig3_list[[ttt]] <- gtable_stack(g1, new_strips)

  fig3_list[[length(treatment)+ttt]] <- ggplot(data=subset(full_data, gas==treatment[ttt]), aes(x=day, y=PREY) ) +
    geom_point(alpha=.5, colour="darkgreen") + 
    stat_smooth(inherit.aes = F, data=subset(full_data, gas==treatment[ttt]), aes(x=day, y=PREY), method = "loess", size = 1.5, se=F, colour="darkgreen")+ #"lm", formula = y ~ x + I(x^2)+I(x^3)+I(x^4)
    geom_point(data=subset(full_data, gas==treatment[ttt]), aes(x=day, y=PRED), colour="red" , alpha=.5) + 
    stat_smooth(inherit.aes = F, data=subset(full_data, gas==treatment[ttt]), aes(x=day, y=PRED), method = "loess", size = 1.5, se=F, colour="red")+ #, formula = y ~ x + I(x^2)
    geom_point(data=subset(full_data, gas==treatment[ttt]), aes(x=day, y=BACT), colour="brown" , alpha=.5)+
    stat_smooth(inherit.aes = F, data=subset(full_data, gas==treatment[ttt]), aes(x=day, y=BACT), method = "loess", size = 1.5, se=F, colour="brown")+ #, formula = y ~ x + I(x^2)+I(x^3)+I(x^4)+I(x^5)
    facet_wrap(~name_microcosm, ncol=6)+
    ylim(0,8) +
    theme_bw() + labs(x="Time (days)", y="Organisms (log10)", title="") +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(), panel.grid = element_blank(),
          plot.margin = unit(c(-0.9,0.1,0,.8), "cm"),
          text = element_text(size=30, family="Times"))
}

blankPlot <- ggplot() + geom_blank(aes(1,1)) + cowplot::theme_nothing()

tiff("Figure4_dynamics.tiff",res=300, width=8000, height=11000)
grid.arrange(blankPlot,
             fig3_list[[1]],fig3_list[[length(treatment)+1]], blankPlot,
             fig3_list[[2]],fig3_list[[length(treatment)+2]], blankPlot,
             fig3_list[[3]],fig3_list[[length(treatment)+3]], blankPlot, 
             fig3_list[[4]],fig3_list[[length(treatment)+4]], blankPlot, 
             fig3_list[[5]],fig3_list[[length(treatment)+5]], blankPlot,
             ncol=1, nrow=16,widths=c(5), heights=c(.2, 3, 4,.2, 3,4,.2, 3,4,.2, 3,4,.2, 3, 4, .2))
dev.off()

#######################################################
### Reproduction of the figure 5 - time series analysis

# some parameter for the figure
min.d.incr <- 1
max.d.decr <- 67

# 4a
plot_variable <- ggplot(variable.scores_a, aes(x=PC1, y=PC2))+
  geom_segment(data = variable.scores_a, aes(xend = variable.scores_a[ ,"PC1"], yend=variable.scores_a[ ,"PC2"]),
               x=0, y=0, colour="black",
               arrow=arrow(angle=25, length=unit(0.25, "cm"),type = "closed")) +
  labs(x=var.expl.PC1, y=var.expl.PC2, title="(a) Variable scores") +
  geom_text(data=variable.scores_a, aes(x= PC1, y= PC2,label=variable), size=5, colour="black", hjust=0, vjust=1.3, family="Times") +
  geom_vline(xintercept = 0, linetype=2, color="grey") + geom_hline(yintercept = 0, linetype=2, color="grey")+
  xlim(-4,9) + ylim(-3.3,2)+
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        plot.margin=unit(c(0.1,0.1,.7,0.1),"cm"),
        panel.grid = element_blank(),
        text = element_text(size=15, family="Times"))

# 4b
plot_time1 <- ggplot(data=time.scores4, aes(x=PC1, y=PC2, group=microcosm, color=overall_cluster)) +
  geom_vline(xintercept = 0, linetype=2, color="grey") + geom_hline(yintercept = 0, linetype=2, color="grey")+
  geom_segment(aes(x=PC1, y=PC2, xend=PC1_t, yend=PC2_t), alpha=.5)+
  scale_color_manual(values=color_cluster)+
  scale_fill_manual(values=color_cluster)+
  stat_ellipse(geom="polygon", aes(group=overall_cluster, fill=overall_cluster), alpha=.3, level=.95)+
  scale_x_discrete(limits=c(0,1))+ scale_y_discrete(limits=c(-1,0,1))+
  labs(x="PC1", y="PC2", title="(b) Temporal scores", fill="Cluster", colour="Cluster") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid = element_blank(),
        plot.margin=unit(c(0.1,0.1,.7,0.1),"cm"),
        strip.background = element_blank(),
        text = element_text(size=15, family="Times"),
        legend.position = c(.85,.25))
# 4c
plot_time2 <- ggplot(data=time.scores4, aes(x=PC1, y=PC2, group=microcosm, color=overall_cluster)) +
  geom_vline(xintercept = 0, linetype=2, color="grey") + geom_hline(yintercept = 0, linetype=2, color="grey")+
  geom_segment(aes(x=PC1, y=PC2, xend=PC1_t, yend=PC2_t))+
  scale_color_manual(values=color_cluster)+
  facet_grid(up_or_down~gas.y) +
  geom_point(data=subset(time.scores4, day==min.d.incr), shape=16, color="black") +
  geom_point(data=subset(time.scores4, day==max.d.decr), shape=15, color="black") +
  scale_x_discrete(limits=c(0,1))+ scale_y_discrete(limits=c(-1,0,1))+
  labs(x="PC1", y="PC2", title="(c) Decomposed temporal scores") +
  theme(panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid = element_blank(),
        plot.margin=unit(c(0.1,0.1,.7,0.1),"cm"),
        strip.background = element_blank(),
        text = element_text(size=15, family="Times"),
        legend.position = "none")


# 4d 
plot_clust_oxygen <- ggplot(segment(ddata_oxygen)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))+
  geom_text(data=labs_oxygen,aes(label=name_microcosm, x=x, y=-.5, color=overall_cluster), angle=70, hjust=1, size=5, family="Times") +
  scale_color_manual(values=color_cluster)+
  xlim(-1,33)+
  scale_y_continuous(limits=c(-210,360), breaks = seq(0,350,by=50))+
  geom_hline(aes(yintercept=250), colour="black", linetype=2)+ annotate("text", x=33, y=270, label="2 clusters", hjust=1, size=5, family="Times")+
  labs(colour="Gas", x="", y="hierarchical DTW distance", title="(d) Dissolved oxygen")+
  geom_segment(aes(x=-Inf, xend=-Inf, y=0, yend=Inf))+
  annotate("rect", xmin=-1,xmax=20,ymin=-170,ymax=-140, fill="white", color="black")+
  annotate("rect", xmin=20.2,xmax=30,ymin=-170,ymax=-140, fill="white", color="black")+
  annotate("text", x=10, y=-154, label="cluster 1", family="Times")+
  annotate("text", x=25, y=-154, label="cluster 2", family="Times")+
  guides(colour = guide_legend(override.aes = list(angle=0)))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        plot.margin=unit(c(0.1,0.1,0,0.1),"cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(hjust = 0.8),
        text = element_text(size=15, family="Times"))

# 4e
plot_clust_biomass <- ggplot(segment(ddata_biomass)) + 
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))+
  geom_text(data=labs_biomass,aes(label=name_microcosm, x=x, y=-.5, color=overall_cluster), angle=70, hjust=1, size=5, family="Times") +
  scale_color_manual(values=color_cluster)+
  xlim(-1,33)+
  scale_y_continuous(limits=c(-105,180), breaks = seq(0,175,by=25))+
  geom_hline(aes(yintercept=125), colour="black", linetype=2)+ annotate("text", x=33, y=135, label="2 clusters", hjust=1, size=5, family="Times")+
  labs(colour="Gas", x="", y="hierarchical DTW distance", title="(e) Total biomass")+
  geom_segment(aes(x=-Inf, xend=-Inf, y=0, yend=Inf))+
  annotate("rect", xmin=-1,xmax=25,ymin=-85,ymax=-70, fill="white", color="black")+
  annotate("rect", xmin=25.2,xmax=30,ymin=-85,ymax=-70, fill="white", color="black")+
  annotate("text", x=13, y=-77, label="cluster 1", family="Times")+
  annotate("text", x=27.7, y=-77, label="cluster 2", family="Times")+
  guides(colour = guide_legend(override.aes = list(angle=0)))+
  theme(legend.position = "none",
        panel.background = element_blank(),
        plot.margin=unit(c(0.1,0.1,0,0.1),"cm"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),       
        axis.title.y = element_text(hjust = 0.8),
        text = element_text(size=15, family="Times"))


blankPlot <- ggplot() + geom_blank(aes(.1,.1)) + cowplot::theme_nothing()

lay <- rbind(c(1,2,3,3),
             c(1,2,3,3),
             c(4,4,5,5),
             c(4,4,5,5),
             c(4,4,5,5))

tiff("Figure5_TimeSeriesAnalysis.tiff",res=300, width=4000, height=3000)
grid.arrange(plot_variable, plot_time1, plot_time2, 
             plot_clust_oxygen, plot_clust_biomass,
             ncol=4, nrow=5,  #,widths=c(10,10),heights=c(10,.2,10),
             layout_matrix = lay)
dev.off()


####################################################################
### Reproduction of the figure 6 - DO and T.biomass over temperature

## Prepare the data (link two consecutive points to create a segment)
data_day1 <- data.frame(microcosm=0, day=0, 
                        oxy_t=0, biom_t=0, temp_t=0)

full_data2 <- as.data.frame(full_data)
for (i in unique(full_data2$microcosm)) {
  #i=101
  df <- subset(full_data2, microcosm==i)
  df <- df[order(df$day),]
  
  oxy_t2 <- c(df[2:nrow(df),"Oxygen_L"],NA)
  biom_t2 <- c(df[2:nrow(df),"Total_biomass"],NA)
  temp_t2 <- c(df[2:nrow(df),"Temperature_L"],NA)
  
  add_df <- data.frame(microcosm=i, day=df$day, 
                       oxy_t=oxy_t2, biom_t=biom_t2, temp_t=temp_t2)
  data_day1 <- rbind(data_day1, add_df)
}

full_data3 <- merge(full_data2, data_day1[-1,],by=c("microcosm","day"))
full_data3_cl <- merge(full_data3, labs_oxygen[,c("label","overall_cluster")], by.x="microcosm", by.y="label")


fig5_list <- list()

for(ttt in 1:length(treatment)){
  # ttt=1
  
  # Oxygen
  plot1 <- ggplot(data=subset(full_data3_cl, gas==treatment[ttt]), aes(x=Temperature_L, y=Oxygen_L, colour=day) ) +
    scale_color_gradient(low = "grey", high = "darkblue")+
    geom_point(size=1) +
    geom_segment(aes(x=Temperature_L, y=Oxygen_L, xend=temp_t, yend=oxy_t), size=2)+
    geom_point(data=subset(full_data3_cl, day==min.d.incr & gas==treatment[ttt]), shape=16, color="black", size=3) +
    geom_point(data=subset(full_data3_cl, day==max.d.decr & gas==treatment[ttt]), shape=15, color="black", size=3) +
    facet_wrap(~name_microcosm, ncol=6)+
    labs(x="",y="% O2", title=treatment[ttt])+  
    ylim(0,20)+
    theme_bw()+ 
    theme(legend.position = "none", 
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_line(size=.2),
          axis.ticks = element_line(size=.2),
          strip.background = element_blank(),
          panel.spacing = unit(2, "pt"), 
          axis.text.x = element_blank(),
          plot.margin = unit(c(0,0.1,0,1.47), "cm"),
          text = element_text(size=30, family="Times"),
          plot.title = element_text(size=30, face="bold"))
  
  
  dummy <- ggplot(data=subset(full_data_cl, gas==treatment[ttt]), aes(x=day, y=Oxygen_H)) + 
    facet_wrap(~name_microcosm, ncol=6)+
    scale_fill_manual(values=color_cluster)+
    geom_rect(aes(fill=overall_cluster), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)+ theme_minimal()+
    theme(strip.text.x = element_text(size=30, color = "white", family="Times"))
  
  g1 <- ggplotGrob(plot1)
  g2 <- ggplotGrob(dummy)
  panels <- grepl(pattern="panel", g2$layout$name)
  strips <- grepl(pattern="strip-t", g2$layout$name)
  g2$layout$t[panels] <- g2$layout$t[panels] - 1
  g2$layout$b[panels] <- g2$layout$b[panels] - 1
  new_strips <- gtable_select(g2, panels | strips)
  # grid.newpage()
  # grid.draw(new_strips)
  gtable_stack <- function(g1, g2){
    g1$grobs <- c(g1$grobs, g2$grobs)
    g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
    g1$layout <- rbind(g1$layout, g2$layout)
    g1
  }
  
  
  fig5_list[[ttt]] <- gtable_stack(g1, new_strips)
  
  
  
  
  fig5_list[[length(treatment)+ttt]] <- ggplot(data=subset(full_data3_cl, gas==treatment[ttt]), aes(x=Temperature_L, y=Total_biomass, colour=day) ) +
    scale_color_gradient(low = "grey", high = "maroon4")+
    geom_point(size=1) +
    geom_segment(aes(x=Temperature_L, y=Total_biomass, xend=temp_t, yend=biom_t), size=2)+
    geom_point(data=subset(full_data3_cl, day==min.d.incr & gas==treatment[ttt]), shape=16, color="black", size=3) +
    geom_point(data=subset(full_data3_cl, day==max.d.decr & gas==treatment[ttt]), shape=15, color="black", size=3) +
    scale_y_continuous(labels = function(x) sprintf("%.2f", x))+
    facet_wrap(~name_microcosm, ncol=6)+#, scales="free"
    labs(x="Temperature (°C)", y= expression(paste("Total biomass (mm"^{3}, "/mL)"))) +  
    theme_bw()+ 
    theme(legend.position = "none", 
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_line(size=.2),
          axis.ticks = element_line(size=.2),
          strip.background = element_blank(),
          strip.text = element_blank(),
          panel.spacing = unit(2, "pt"), 
          plot.margin = unit(c(-0.9,0.1,0,0.57), "cm"),
          text = element_text(size=27, family="Times"))
  
}

blankPlot <- ggplot() + geom_blank(aes(1,1)) + cowplot::theme_nothing()

tiff("Figure6_versusTemperature.tiff",res=300, width=8000, height=11000)
grid.arrange(blankPlot,
             fig5_list[[1]],fig5_list[[length(treatment)+1]], blankPlot,
             fig5_list[[2]],fig5_list[[length(treatment)+2]], blankPlot,
             fig5_list[[3]],fig5_list[[length(treatment)+3]], blankPlot, 
             fig5_list[[4]],fig5_list[[length(treatment)+4]], blankPlot, 
             fig5_list[[5]],fig5_list[[length(treatment)+5]], blankPlot,
             ncol=1, nrow=16,widths=c(5), heights=c(.2, 4, 4,.2, 4, 4,.2, 4, 4,.2, 4, 4,.2, 4, 4, .2))
dev.off()





#######################################################################################
### Reproduction of the figure 7 - Nonlinearity and hysteresis for DO and Total biomass

## OXYGEN

plot_nl_oxygen_incr <- ggplot(subset(nl_oxygen,up_or_down=="increasing"), aes(x=gas, y=get_nl, shape=as.factor(cluster), fill=as.factor(cluster))) + 
  scale_fill_manual(values=c("black","white"))+
  scale_shape_manual(values=c(21,24))+
  geom_jitter(size=2,width = 0.15, alpha=.8) + 
  ylim(-0.1,3)+
  labs(y="Nonlinearity", x="", title="Dissolved oxygen", subtitle = "(a) Increasing temperature")+
  theme_bw()+
  theme(text = element_text(family="Times"),
        axis.text.x = element_text(angle=45, vjust=.7),
        strip.background = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

plot_nl_oxygen_decr <- ggplot(subset(nl_oxygen,up_or_down=="decreasing"), aes(x=gas, y=get_nl, fill=as.factor(cluster), shape=as.factor(cluster))) + 
  scale_fill_manual(values=c("black","white"))+
  scale_shape_manual(values=c(21,24))+
  geom_jitter(size=2,width = 0.15, alpha=.8) + 
  ylim(-0.1,3)+
  labs(y="Nonlinearity", x="", title="", subtitle = "(b) Decreasing temperature")+
  theme_bw()+
  theme(text = element_text(family="Times"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=45, vjust=.7),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

Hyst_oxygen <- ggplot(hh_oxygen_diff, aes(x=gas, y=mean_diff, fill=as.factor(cluster), shape=as.factor(cluster))) +
  geom_jitter(size=2,width = 0.15, alpha=.8) + 
  ylim(0,9)+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values=c("black","white"))+
  labs(x="", y="Mean difference", title="", subtitle="(c) Hysteresis")+
  theme_bw()+
  theme(text = element_text(family="Times"),
        axis.text.x = element_text(angle=45, vjust=.7),
        strip.background = element_blank(),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())


## BIOMASS

plot_nl_biomass_incr <- ggplot(subset(nl_biomass,up_or_down=="increasing"), aes(x=gas, y=get_nl, shape=as.factor(cluster), fill=as.factor(cluster))) + 
  scale_fill_manual(values=c("black","white"))+
  scale_shape_manual(values=c(21,24))+
  geom_jitter(size=2,width = 0.15, alpha=.8) + 
  ylim(-0.005,0.065)+
  labs(y="Nonlinearity", x="", title="Total biomass", subtitle = "(d) Increasing temperature")+
  theme_bw()+
  theme(text = element_text(family="Times"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=45, vjust=.7),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

plot_nl_biomass_decr <- ggplot(subset(nl_biomass,up_or_down=="decreasing"), aes(x=gas, y=get_nl, fill=as.factor(cluster), shape=as.factor(cluster))) + 
  scale_fill_manual(values=c("black","white"))+
  scale_shape_manual(values=c(21,24))+
  geom_jitter(size=2,width = 0.15, alpha=.8) + 
  ylim(-0.005,0.065)+
  labs(y="Nonlinearity", x="Gas exchange treatment", title="", subtitle = "(e) Decreasing temperature")+
  theme_bw()+
  theme(text = element_text(family="Times"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=45, vjust=.7),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())

Hyst_biomass <- ggplot(hh_biomass_diff, aes(x=gas, y=mean_diff, fill=as.factor(cluster), shape=as.factor(cluster))) +
  geom_jitter(size=2, alpha=.8, width=0.15) + 
  ylim(0,.15)+
  scale_shape_manual(values = c(21,24))+
  scale_fill_manual(values=c("black","white"))+
  labs(x="", y="Mean difference", title="", subtitle="(f) Hysteresis")+
  theme_bw()+
  theme(text = element_text(family="Times"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=45, vjust=.7),
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())



tiff("Figure7_Results.tiff", res=600, width=4400, height=3200)
gridExtra::grid.arrange(plot_nl_oxygen_incr, plot_nl_oxygen_decr, Hyst_oxygen,
                        plot_nl_biomass_incr, plot_nl_biomass_decr, Hyst_biomass,
                        ncol=3, nrow=2)
dev.off()








