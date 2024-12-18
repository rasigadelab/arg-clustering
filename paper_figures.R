
library(ggplot2)
library(ggpubr)

#### Figures specs ####
# A4 210 * 297mm
fig_specs <- list(
  path = "./Figures",
  width = 9,
  height = 7,
  dpi = 200
)
fig_specs_4plots <- list(
  path = "./Figures",
  width = 20,
  height = 6.5,
  dpi = 300
)
fig_specs_4heatmaps <- list(
  path = "./Figures",
  width = 22,
  height = 6,
  dpi = 300
)

alpha_heatmap <- 0.8
t=simLength

## Colors ##
#color_palette <- c("steelblue1","gold","firebrick1") # heatmaps (continuous)
#color_palette <- c("blue2","yellow1","red1") # heatmaps (continuous)
#color_palette_reverse <- c("red1","yellow1","blue2")
color_lines <- c("blue2","green3","yellow2","orange2","red","blue") # lines (discrete)
#color_lines_3 <- c("blue2","green3","orange2","red","blue") # lines (discrete)
color_palette <- c("navy","purple2","gold","darkorange","red3") # heatmaps (continuous)
color_palette_reverse <- c("red3","darkorange","gold","purple2","navy")
# color_lines <- c("blue2","dodgerblue","green2","yellow2","red","blue")

# P <- P + scale_colour_viridis_d()
# P <- P + scale_fill_gradient(low="skyblue", high="red1", limits=c(0,0.002))
# P <- P + scale_fill_gradient2(low="blue2", mid="yellow1", midpoint=0.001, high="red1", limits=c(0,0.002)) # ok+
# P <- P + scale_fill_gradientn(colours = c("blue2","yellow1","red1"), limits=c(0,0.002)) # ok+
# P <- P + scale_fill_gradientn(colours = c("blue2","green1","yellow1","orange1","red1"), limits=c(0,0.002)) # ok
# P <- P + scale_fill_viridis_c(option = "plasma", limits=c(0,0.002)) # ok
# P <- P + scale_fill_distiller(palette = "RdPu") # bof
# P <- P + scale_fill_distiller(palette = "YlOrBr") # bof
# P <- P + scale_fill_distiller(palette = "RdYlBu") # ok+

# color test #
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x00_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x00_env1_sp1/mean_N_env1_sp1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("test.tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)


#### Data preparation - 1env 1species ####

### Data preparation ###
nbPop <- 1 # 1 species
out_dt <- setDT(as.data.frame(out_parallel))
setnames(out_dt, "rUpt","rUpt_1")
## Add N and Pab ##
for(i in 1:nbPop){
  out_dt[, paste0("m_",i) := get(paste0("NEf1_",i)) / get(paste0("x11_",i))] # mean arc distance
  out_dt[, paste0("S_",i) := get(paste0("NEf2_",i)) / get(paste0("x11_",i)) - get(paste0("m_",i))^2] # variance arc distance
  out_dt[, paste0("N_",i) := get(paste0("x00_",i)) + get(paste0("x01_",i)) + get(paste0("x10_",i)) + get(paste0("x11_",i))] # size population species 1
  out_dt[, paste0("ms_",i) := paste0(get(paste0("m_",i)),";",get(paste0("S_",i)))]
  out_dt[, paste0("alpha_",i) := sapply(get(paste0("ms_",i)), function(x){
    m=as.numeric(strsplit(x,";")[[1]][1]);
    s=as.numeric(strsplit(x,";")[[1]][2]);
    return(absolve(m, s)[1])})]
  out_dt[, paste0("beta_",i) := sapply(get(paste0("ms_",i)), function(x){
    m=as.numeric(strsplit(x,";")[[1]][1]);
    s=as.numeric(strsplit(x,";")[[1]][2]);
    return(absolve(m, s)[2])})]
  out_dt[, paste0("ab_",i) := paste0(get(paste0("alpha_",i)),";",get(paste0("beta_",i)))]
  out_dt[, paste0("pab_",i) := sapply(get(paste0("ab_",i)), function(x){
    a=as.numeric(strsplit(x,";")[[1]][1]);
    b=as.numeric(strsplit(x,";")[[1]][2]);
    return(comobg(50,a,b))})]
}

## Mean last time steps ##
out_dt_last_steps <- out_dt[time >= (simLength-4000)] # subset last time steps
out_dt_last_steps[, mean_N_1 := mean(N_1), by=param_set]
out_dt_last_steps[, mean_x00_1 := mean(x00_1), by=param_set]
out_dt_last_steps[, mean_x01_1 := mean(x01_1), by=param_set]
out_dt_last_steps[, mean_x11_1 := mean(x11_1), by=param_set]
out_dt_last_steps[, mean_Dna11_1 := mean(Dna11_1), by=param_set]
out_dt_last_steps[, mean_Dna01_1 := mean(Dna01_1), by=param_set]
out_dt_last_steps[, mean_pab_1 := mean(pab_1), by=param_set]

## Parameters to log (if needed) ##
out_dt_last_steps[, param_1 := log10(param_1)]
out_dt_last_steps[, param_2 := log10(param_2)]

## Antibiotic dynamics ##
# 0=no pressure, 1=pressure
dp <- pars["deathPeriod_1"]
# Treatment 1 : relaxed pressure
out_dt[treatment==1, atb1 := 0]
out_dt[treatment==1, atb2 := 0]
# Treatment 2 : constant pressure
out_dt[treatment==2, atb1 := 1]
out_dt[treatment==2, atb2 := 1]
# Treatment 3 : periodic pressure
out_dt[treatment==3, atb1 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, dp, trt=3)[1]})]
out_dt[treatment==3, atb2 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, dp, trt=3)[1]})]
# Treatment 4 : alternating pressure
out_dt[treatment==4, atb1 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, dp, trt=4)[2]})]
out_dt[treatment==4, atb2 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, dp, trt=4)[3]})]
out_dt[atb1==0, atb1 := NA]
out_dt[atb2==0, atb2 := NA]


#### Data preparation - 1env 2species ####
## Init ##
nbPop <- 2
out_dt <- setDT(as.data.frame(out_parallel))
## Add N & Pab ##
for(i in 1:nbPop){ # add columns
  out_dt[, paste0("m_",i) := get(paste0("NEf1_",i)) / get(paste0("x11_",i))] # mean arc distance
  out_dt[, paste0("S_",i) := get(paste0("NEf2_",i)) / get(paste0("x11_",i)) - get(paste0("m_",i))^2] # variance arc distance
  out_dt[, paste0("N_",i) := get(paste0("x00_",i)) + get(paste0("x01_",i)) + get(paste0("x10_",i)) + get(paste0("x11_",i))] # size population species 1
  out_dt[, paste0("ms_",i) := paste0(get(paste0("m_",i)),";",get(paste0("S_",i)))]
  out_dt[, paste0("alpha_",i) := sapply(get(paste0("ms_",i)), function(x){
    m=as.numeric(strsplit(x,";")[[1]][1]);
    s=as.numeric(strsplit(x,";")[[1]][2]);
    return(absolve(m, s)[1])})]
  out_dt[, paste0("beta_",i) := sapply(get(paste0("ms_",i)), function(x){
    m=as.numeric(strsplit(x,";")[[1]][1]);
    s=as.numeric(strsplit(x,";")[[1]][2]);
    return(absolve(m, s)[2])})]
  out_dt[, paste0("ab_",i) := paste0(get(paste0("alpha_",i)),";",get(paste0("beta_",i)))]
  out_dt[, paste0("pab_",i) := sapply(get(paste0("ab_",i)), function(x){
    a=as.numeric(strsplit(x,";")[[1]][1]);
    b=as.numeric(strsplit(x,";")[[1]][2]);
    return(comobg(50,a,b))})]
}

## Mean last time steps species1 ##
out_dt_last_steps <- out_dt[time >= (simLength-4000)] # subset last time steps
out_dt_last_steps[, mean_N_1 := mean(N_1), by=param_set]
out_dt_last_steps[, mean_x00_1 := mean(x00_1), by=param_set]
out_dt_last_steps[, mean_x01_1 := mean(x01_1), by=param_set]
out_dt_last_steps[, mean_x10_1 := mean(x10_1), by=param_set]
out_dt_last_steps[, mean_x11_1 := mean(x11_1), by=param_set]
out_dt_last_steps[, mean_Dna01_1 := mean(Dna01_1), by=param_set]
out_dt_last_steps[, mean_Dna10_1 := mean(Dna10_1), by=param_set]
out_dt_last_steps[, mean_Dna11_1 := mean(Dna11_1), by=param_set]
## Mean last time steps species2 ##
out_dt_last_steps[, mean_N_2 := mean(N_2), by=param_set]
out_dt_last_steps[, mean_x00_2 := mean(x00_2), by=param_set]
out_dt_last_steps[, mean_x01_2 := mean(x01_2), by=param_set]
out_dt_last_steps[, mean_x10_2 := mean(x10_2), by=param_set]
out_dt_last_steps[, mean_x11_2 := mean(x11_2), by=param_set]
out_dt_last_steps[, mean_Dna01_2 := mean(Dna01_2), by=param_set]
out_dt_last_steps[, mean_Dna10_2 := mean(Dna10_2), by=param_set]
out_dt_last_steps[, mean_Dna11_2 := mean(Dna11_2), by=param_set]
#View(out_dt_last_steps[time==simLength])

## Params to log (if needed) ##
out_dt_last_steps[, param_1 := log10(param_1)]
out_dt_last_steps[, param_2 := log10(param_2)]

## Antibiotic dynamics ##
# 0=no pressure, 1=pressure
dp <- pars["deathPeriod_1"]
# Treatment 1 : relaxed pressure
out_dt[treatment==1, atb1 := 0]
out_dt[treatment==1, atb2 := 0]
# Treatment 2 : constant pressure
out_dt[treatment==2, atb1 := 1]
out_dt[treatment==2, atb2 := 1]
# Treatment 3 : periodic pressure
out_dt[treatment==3, atb1 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, dp, trt=3)[1]})]
out_dt[treatment==3, atb2 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, dp, trt=3)[1]})]
# Treatment 4 : alternating pressure
out_dt[treatment==4, atb1 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, dp, trt=4)[2]})]
out_dt[treatment==4, atb2 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, dp, trt=4)[3]})]
out_dt[atb1==0, atb1 := NA]
out_dt[atb2==0, atb2 := NA]

# Treatment 5 : periodic uneven recovery period (life period) !ONGOING!
out_dt[treatment==5, atb1 := deathRateFunc_tris(time, 0, 1, dp=deathPeriod, trt=5, lp=lifePeriod)[2]]




#### Figures 2.1 - dynamics=f(mob)/f(dup) - 1env 1species ####
trt <- 3
param_subset <- c(0,1E-3,1E-2,1E-1)
### Style 1 ###
# x11=f(t) #
for(trt in 1:4){
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% param_subset])
  P <- P + ggtitle(paste0("x11=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=log10(x11_1), color=factor(param_1)))
  P <- P + geom_line(aes(x=time, y=atb1+0.1, color="atb1"))
  P <- P + geom_line(aes(x=time, y=atb2+0.3, color="atb2"))
  P <- P + scale_color_manual(values=color_lines)
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="x11", limits=c(1,10), breaks=seq(0,10,2))
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("x11=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# Conditional Pab=f(t) log #
Es <- as.double(muvar(1, 50)[1])
for(trt in 1:4){
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% param_subset])
  P <- P + ggtitle(paste0("pab_cond_log=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=log10(pab_1 / Es), color=factor(param_1)))
  P <- P + geom_line(aes(x=time, y=atb1-3.15, color="atb1"))
  P <- P + geom_line(aes(x=time, y=atb2-3.1, color="atb2"))
  P <- P + scale_color_manual(values=color_lines)
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="log pab cond", limits=c(-2.2,0), breaks=seq(-2,0,0.5))
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("pab_cond_log=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}

### Style 2 ###
# x11=f(t) #
for(trt in 1:4){
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% param_subset][x11_1>1E2])
  P <- P + ggtitle(paste0("x11=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=log10(x11_1), color=factor(param_1)))
  P <- P + geom_line(aes(x=time, y=atb1+0.4, color="Drug A"))
  P <- P + geom_line(aes(x=time, y=atb2+0.05, color="Drug B"))
  P <- P + scale_color_manual(values=color_lines)
  P <- P + scale_x_continuous(name="generation", limits=c(0,16E3), expand=c(0,0), breaks=seq(0,simLength,4000), labels=c("0", paste0(seq(4,simLength/1E3,4),"E3")))
  P <- P + scale_y_continuous(name="x11", limits=c(1,10), breaks=seq(0,10,2))
  P <- P + labs(color="param_1")
  P <- P + theme_minimal()
  P <- P + theme(panel.grid.minor.y=element_line(color="white")) # remove minor y grid
  P <- P + geom_hline(yintercept = 2, linetype = "solid", color = "grey30", size=0.3)
  P <- P + geom_point(data=data.frame(x=seq(0,16000,4000),y=2), aes(x=x, y=y), shape=3, color="grey30")
  P <- P + geom_segment(x=0,y=2,xend=0,yend=10,color='grey30')
  P <- P + geom_point(data=data.frame(x=0,y=seq(2,10,2)), aes(x=x, y=y), shape=3, color="grey30")
  plot(P)
  assign(paste0("P_x11_trt",trt), P)
  ggsave(paste0("x11=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# Conditional Pab=f(t) log #
Es <- as.double(muvar(1, 50)[1])
for(trt in 1:4){
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% param_subset])
  P <- P + ggtitle(paste0("pab_cond_log=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=log10(pab_1 / Es), color=factor(param_1)))
  P <- P + geom_line(aes(x=time, y=atb1-3.3, color="Drug A"))
  P <- P + geom_line(aes(x=time, y=atb2-3.4, color="Drug B"))
  P <- P + scale_color_manual(values=color_lines)
  P <- P + scale_x_continuous(name="generation", limits=c(0,16E3), breaks=seq(0,simLength,4000), labels=c("0", paste0(seq(4,simLength/1E3,4),"E3")), expand=c(0,0)) # expand removes margin
  P <- P + scale_y_continuous(name="log pab cond", limits=c(-2.4,0), breaks=seq(-2,0,0.5))
  P <- P + labs(color="param_1")
  #P <- P + theme_light()
  P <- P + theme_minimal()
  P <- P + theme(panel.grid.minor.y=element_line(color="white")) # remove minor y grid
  #P <- P + theme(axis.text.x = element_text(vjust = 25)) # move axis text
  P <- P + geom_hline(yintercept = -2, linetype = "solid", color = "grey30", size=0.3)
  P <- P + geom_point(data=data.frame(x=seq(0,16000,4000),y=-2), aes(x=x, y=y), shape=3, color="grey30")
  P <- P + geom_segment(x=0,y=-2,xend=0,yend=0,color='grey30')
  P <- P + geom_point(data=data.frame(x=0,y=seq(-2,0,0.5)), aes(x=x, y=y), shape=3, color="grey30")
  plot(P)
  assign(paste0("P_pab_trt",trt), P)
  ggsave(paste0("pab_cond_log=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}


### Fuse plots ###
## x11 ##
for(i in 1:4){
  P <- get(paste0("P_x11_trt",i))
  # remove y axis titles
  P <- P + theme(axis.title.y=element_blank())
  # x axis text size
  P <- P + theme(axis.title.x=element_text(size=9))
  assign(paste0("P_x11_trt",i), P)
}
P <- ggarrange(P_x11_trt1, P_x11_trt2, P_x11_trt3, P_x11_trt4,
               labels=c("A","B","C","D"), ncol=4, nrow=1,
               common.legend=TRUE, legend="right")
plot(P)
ggsave(paste0("x11=f(t)_param1_all_trt.tiff"), width=fig_specs_4plots$width, height=fig_specs_4plots$height, units="cm", path=fig_specs_4plots$path, dpi=fig_specs_4plots$dpi)

## Pab cond ##
for(i in 1:4){
  P <- get(paste0("P_pab_trt",i))
  # remove y axis titles
  P <- P + theme(axis.title.y=element_blank())
  # x axis text size
  P <- P + theme(axis.title.x=element_text(size=9))
  assign(paste0("P_pab_trt",i), P)
}
P <- ggarrange(P_pab_trt1, P_pab_trt2, P_pab_trt3, P_pab_trt4,
               labels=c("A","B","C","D"), ncol=4, nrow=1,
               common.legend=TRUE, legend="right")
plot(P)
ggsave(paste0("pab=f(t)_param1_all_trt.tiff"), width=fig_specs_4plots$width, height=fig_specs_4plots$height, units="cm", path=fig_specs_4plots$path, dpi=fig_specs_4plots$dpi)




### Other possible plots (not included in paper figures) ###
# m=f(t) #
for(trt in 1:4){
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% param_subset])
  P <- P + ggtitle(paste0("m=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=m_1, color=factor(param_1)))
  P <- P + scale_color_manual(values=color_lines)
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="m", limits=c(0,0.25001))
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("m=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# Pab=f(t) #
for(trt in 1:4){
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% param_subset])
  P <- P + ggtitle(paste0("pab=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=pab_1, color=factor(param_1)))
  P <- P + scale_color_manual(values=color_lines)
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="pab", limits=c(0,0.01))
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("pab=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# Conditional Pab=f(t) #
Es <- as.double(muvar(1, 50)[1])
for(trt in 1:4){
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% param_subset])
  P <- P + ggtitle(paste0("pab_cond=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=pab_1 / Es, color=factor(param_1)))
  P <- P + geom_line(aes(x=time, y=atb1-1.08, color="atb1"))
  P <- P + geom_line(aes(x=time, y=atb2-1.06, color="atb2"))
  P <- P + scale_color_manual(values=color_lines)
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="pab cond", limits=c(-0.1,1), breaks=seq(0,1,0.2))
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("pab_cond=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# x00=f(t) #
for(trt in 1:4){
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% param_subset])
  P <- P + ggtitle(paste0("x00=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=log10(x00_1), color=factor(param_1)))
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="x00", limits=c(0,10))
  P <- P + scale_color_manual(values=color_lines)
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("x00=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# x01=f(t) #
for(trt in 1:4){
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% param_subset])
  P <- P + ggtitle(paste0("x01=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=log10(x01_1), color=factor(param_1)))
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="x01", limits=c(0,10))
  P <- P + scale_color_manual(values=color_lines)
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("x01=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}




#### Figures 2.1 - competition dynamics=f(mob)/f(dup) - 1env 2species ####

param_subset <- c(0,1E-3,1E-2,1E-1)
trt <- 3
### Style 2 ###
## N_1/N_tot=f(t) ##
# Supplementary figure specs
time_limit <- simLength
time_break <- 1E4
time_breaks <- seq(0,simLength,time_break)
# Figure 2.1 specs
time_limit <- 16000
time_break <- 4000
time_breaks <- seq(0,simLength,time_break)
time_labels <- c("0", paste0(seq(4,simLength/1E3,4),"E3"))
for(trt in 1:4){
  P <- ggplot(out_dt[treatment==trt][param_1 %in% param_subset])#[x11_1>1E-6]
  P <- P + ggtitle(paste0("N_1/N_tot=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=N_1/(N_1+N_2), color=factor(param_1)))
  P <- P + geom_line(aes(x=time, y=atb1-1.14, color="Drug A"))
  P <- P + geom_line(aes(x=time, y=atb2-1.09, color="Drug B"))
  P <- P + scale_color_manual(values=color_lines)
  P <- P + scale_x_continuous(name="generation", limits=c(0,time_limit), expand=c(0,0), breaks=time_breaks, labels=time_labels)
  #P <- P + scale_x_continuous(name="generation", limits=c(0,time_limit), expand=c(0,0), breaks=time_breaks, labels = function(x) format(x, scientific = TRUE))
  P <- P + scale_y_continuous(name="N1/Ntot",limits=c(-0.15,1), breaks=seq(0,1,0.2))
  P <- P + labs(color="param_1")
  P <- P + theme_minimal()
  P <- P + theme(panel.grid.minor.y=element_line(color="white")) # remove minor y grid
  P <- P + geom_hline(yintercept=0, linetype = "solid", color = "grey30", size=0.3)
  P <- P + geom_point(data=data.frame(x=seq(0,time_limit,time_break),y=0), aes(x=x, y=y), shape=3, color="grey30")
  P <- P + geom_segment(x=0,y=0,xend=0,yend=1,color='grey30')
  P <- P + geom_point(data=data.frame(x=0,y=seq(0,1,0.2)), aes(x=x, y=y), shape=3, color="grey30")
  assign(paste0("P_N1_trt",trt), P)
  plot(P)
  ggsave(paste0("N1=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}

### Fuse plots ###
## N1 ##
for(i in 1:4){
  P <- get(paste0("P_N1_trt",i))
  # remove y axis titles
  P <- P + theme(axis.title.y=element_blank())
  # x axis text size
  P <- P + theme(axis.title.x=element_text(size=9))
  assign(paste0("P_N1_trt",i), P)
}
P <- ggarrange(P_N1_trt1, P_N1_trt2, P_N1_trt3, P_N1_trt4,
               labels=c("A","B","C","D"), ncol=4, nrow=1,
               common.legend=TRUE, legend="right")
plot(P)
ggsave(paste0("N1=f(t)_param1_all_trt.tiff"), width=fig_specs_4plots$width, height=fig_specs_4plots$height, units="cm", path=fig_specs_4plots$path, dpi=fig_specs_4plots$dpi)




### Other plots ###
# N_1/N_tot=f(t) #
for(trt in 1:4){
  P <- ggplot(out_dt[treatment==trt][param_1 %in% param_subset])#[x11_1>1E-6]
  P <- P + ggtitle(paste0("N_1/N_tot=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=N_1/(N_1+N_2), color=factor(param_1)))
  P <- P + geom_line(aes(x=time, y=atb1-1.09, color="atb1"))
  P <- P + geom_line(aes(x=time, y=atb2-1.08, color="atb2"))
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="N1/Ntot",limits=c(-0.1,1), breaks=seq(0,1,0.2))
  P <- P + scale_color_manual(values=color_lines, breaks=c(as.character(param_subset),"atb1","atb2"))
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("N1=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# N_2/N_tot=f(t) #
for(trt in 1:4){
  P <- ggplot(out_dt[treatment==trt][param_1 %in% param_subset])#[x11_1>1E-6]
  P <- P + ggtitle(paste0("N_2/N_tot=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=N_2/(N_1+N_2), color=factor(param_1)))
  P <- P + geom_line(aes(x=time, y=atb1-1.09, color="atb1"))
  P <- P + geom_line(aes(x=time, y=atb2-1.08, color="atb2"))
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="N2/Ntot",limits=c(-0.1,1), breaks=seq(0,1,0.2))
  P <- P + scale_color_manual(values=color_lines, breaks=c(as.character(param_subset),"atb1","atb2"))
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("N2=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# Cells 1sim #
color_cells <- c("green1","blue1","red1","magenta","blue","red")
for(trt in 1:4){
  P <- ggplot(out_dt[treatment==trt][param_1==1E-3])#[x11_1>1E-6]
  P <- P + ggtitle(paste0("Cell dynamics trt",trt))
  P <- P + geom_line(aes(x=time, y=log10(x00_1), color="x00", linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=log10(x10_1), color="x10", linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=log10(x01_1), color="x01", linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=log10(x11_1), color="x11", linetype="sp1"))
  
  P <- P + geom_line(aes(x=time, y=log10(x00_2), color="x00", linetype="sp2"))
  P <- P + geom_line(aes(x=time, y=log10(x10_2), color="x10", linetype="sp2"))
  P <- P + geom_line(aes(x=time, y=log10(x01_2), color="x01", linetype="sp2"))
  P <- P + geom_line(aes(x=time, y=log10(x11_2), color="x11", linetype="sp2"))
  
  P <- P + geom_line(aes(x=time, y=atb1-1.15, color="atb1"))
  P <- P + geom_line(aes(x=time, y=atb2-1.05, color="atb2"))
  
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="cells",limits=c(-0.2,10), breaks=seq(0,10,2))
  
  P <- P + scale_color_manual(values=color_cells, breaks=c("x00","x10","x01","x11","atb1","atb2"))
  
  P <- P + labs(color="", linetype="species")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("cells=f(t)_param1_trt",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# m species 1 #
for(trt in 1:4){
  P <- ggplot(out_dt[treatment==trt][param_1 %in% param_subset][x11_1>1E-6])#
  P <- P + ggtitle(paste0("m trt",trt))
  P <- P + geom_line(aes(x=time, y=m_1, color=factor(param_1)))
  P <- P + geom_line(aes(x=time, y=atb1-1.09, color="atb1"))
  P <- P + geom_line(aes(x=time, y=atb2-1.08, color="atb2"))
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="N1/Ntot",limits=c(0,0.251), breaks=seq(0,0.251,0.05))
  P <- P + scale_color_manual(values=color_lines, breaks=c(as.character(param_subset),"atb1","atb2"))
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("m1=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}





#### Figure 2.2 - heatmaps mob/dup - 1env 1species ####
trt <- 3
# x11 #
for(trt in 1:4){
  P <- ggplot(out_dt_last_steps[treatment==trt])
  P <- P + ggtitle(paste0("mean_x11_1/mean_N_1 trt=",trt))
  P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_1/mean_N_1), alpha=alpha_heatmap)
  P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1), breaks=seq(0,1,0.2))
  P <- P + scale_x_continuous(name="deletion rate")
  P <- P + scale_y_continuous(name="transfer rate")
  P <- P + labs(fill="")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("Heatmap_x11_t=",t,"_trt",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
  assign(paste0("P_x11_trt",trt), P)
}
# x00 #
for(trt in 1:4){
  P <- ggplot(out_dt_last_steps[treatment==trt])
  P <- P + ggtitle(paste0("mean_x00_1/mean_N_1 trt=",trt))
  P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x00_1/mean_N_1), alpha=alpha_heatmap)
  P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1), breaks=seq(0,1,0.2))
  P <- P + scale_x_continuous(name="deletion rate")
  P <- P + scale_y_continuous(name="transfer rate")
  P <- P + labs(fill="")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("Heatmap_x00_t=",t,"_trt",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
  assign(paste0("P_x00_trt",trt), P)
}
# Conditional Pab #
Es <- as.double(muvar(1,50)[1])
out_dt_last_steps[time==simLength & x11_1<=1E-6, pab_1:=NA]
for(trt in 1:4){
  #P <- ggplot(out_dt_last_steps[time==simLength][x11_1>1E-6][treatment==trt])
  P <- ggplot(out_dt_last_steps[time==simLength][treatment==trt])
  P <- P + ggtitle(paste0("Pab_cond trt=",trt))
  P <- P + geom_raster(aes(x=param_1, y=param_2, fill=pab_1/Es), alpha=alpha_heatmap)
  if(pars[["fcost_1"]]==0.03){
    P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,0.1501), breaks=seq(0,0.2,0.05)) # for cost 3%
  }
  if(pars[["fcost_1"]]==0.05){
    P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,0.25), breaks=seq(0,0.25,0.05)) # for cost 5%
  }
  P <- P + scale_x_continuous(name="deletion rate")#, limits=c(min(out_dt_last_steps$param_1),max(out_dt_last_steps$param_1)))
  P <- P + scale_y_continuous(name="transfer rate")#, limits=c(min(out_dt_last_steps$param_2),max(out_dt_last_steps$param_2)))
  P <- P + labs(fill="")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("Heatmap_pab_cond_t=",t,"_trt",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
  assign(paste0("P_pab_trt",trt), P)
}


### Fuse plots ###
## x11 ##
for(i in 1:4){
  P <- get(paste0("P_x11_trt",i))
  # remove y axis titles
  P <- P + theme(axis.title.y=element_blank())
  # x axis text size
  P <- P + theme(axis.title.x=element_text(size=9))
  assign(paste0("P_x11_trt",i), P)
}
P <- ggarrange(P_x11_trt1, P_x11_trt2, P_x11_trt3, P_x11_trt4,
               labels=c("A","B","C","D"), ncol=4, nrow=1,
               common.legend=TRUE, legend="right")
plot(P)
ggsave(paste0("x11=f(t)_param1_all_trt.tiff"), width=fig_specs_4heatmaps$width, height=fig_specs_4heatmaps$height, units="cm", path=fig_specs_4heatmaps$path, dpi=fig_specs_4heatmaps$dpi)

## x00 ##
for(i in 1:4){
  P <- get(paste0("P_x00_trt",i))
  # remove y axis titles
  P <- P + theme(axis.title.y=element_blank())
  # x axis text size
  P <- P + theme(axis.title.x=element_text(size=9))
  assign(paste0("P_x00_trt",i), P)
}
P <- ggarrange(P_x00_trt1, P_x00_trt2, P_x00_trt3, P_x00_trt4,
               labels=c("A","B","C","D"), ncol=4, nrow=1,
               common.legend=TRUE, legend="right")
plot(P)
ggsave(paste0("x00=f(t)_param1_all_trt.tiff"), width=fig_specs_4heatmaps$width, height=fig_specs_4heatmaps$height, units="cm", path=fig_specs_4heatmaps$path, dpi=fig_specs_4heatmaps$dpi)

## Pab cond ##
for(i in 1:4){
  P <- get(paste0("P_pab_trt",i))
  # remove y axis titles
  P <- P + theme(axis.title.y=element_blank())
  # x axis text size
  P <- P + theme(axis.title.x=element_text(size=9))
  assign(paste0("P_pab_trt",i), P)
}
P <- ggarrange(P_pab_trt1, P_pab_trt2, P_pab_trt3, P_pab_trt4,
               labels=c("A","B","C","D"), ncol=4, nrow=1,
               common.legend=TRUE, legend="right")
plot(P)
ggsave(paste0("pab=f(t)_param1_all_trt.tiff"), width=fig_specs_4heatmaps$width, height=fig_specs_4heatmaps$height, units="cm", path=fig_specs_4heatmaps$path, dpi=fig_specs_4heatmaps$dpi)





### Other plots not included in the paper figure ###
# x00 #
#P <- ggplot(out_dt_last_steps)
P <- ggplot(out_dt_last_steps[time==simLength][treatment==trt])
P <- P + ggtitle(paste0("mean_x00_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x00_1/mean_N_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x00_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
# x01 #
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x01_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x01_1/mean_N_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x01_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
# Arc distance m #
P <- ggplot(out_dt_last_steps[time==simLength][x11_1>1E-6])
P <- P + ggtitle(paste0("m t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette_reverse)
P <- P + scale_fill_gradientn(colours=color_palette_reverse, limits=c(0.2,0.25001))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_m_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
# Pab #
P <- ggplot(out_dt_last_steps[time==simLength][x11_1>1E-6])
P <- P + ggtitle(paste0("Pab t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=pab_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,0.0015))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_pab_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
# Dna11 #
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna11_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna11_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna11_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
# Dna01 #
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna01_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna01_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna01_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)





#### Figure 2.3 - 1env 2species ####

## Heatmap colors species competition ##
#color_palette <- c("navy","purple2","gold","darkorange","red3") # heatmaps (continuous)
color_palette_compete <- c("navy","purple2","ivory","darkorange","red3") # heatmaps (continuous)
trt=3
# N_1 / (N_1 + N_2)
for(trt in 1:4){
  P <- ggplot(out_dt_last_steps[time==simLength][treatment==trt])
  P <- P + ggtitle(paste0("N_1 / (N_1 + N_2)"))
  P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_N_1 / (mean_N_1 + mean_N_2)), alpha=alpha_heatmap)
  P <- P + scale_fill_gradientn(colours=color_palette_compete, limits=c(0,1), breaks=seq(0,1,0.2))
  P <- P + scale_x_continuous(name="param1")
  P <- P + scale_y_continuous(name="param2")
  P <- P + labs(fill="")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("Heatmap_N_1_t=",t,"_trt",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
  assign(paste0("P_N1_trt",trt), P)
}

# Cell dynamics - genotypes=f(t) #
color_pop <- c(susceptible="green2",single_ARG_A="blue2",single_ARG_B='red2',double_ARG="magenta2",DrugA="blue",DrugB="red")
# subset oen simulation from heatmap data
out_dt_bis <- out_dt[param_1==1E-3][param_2==1E-1]
min_cell <- 1E0
out_dt_bis[x10_1<min_cell, x10_1:=NA]
out_dt_bis[x10_2<min_cell, x10_2:=NA]
out_dt_bis[x01_1<min_cell, x01_1:=NA]
out_dt_bis[x01_2<min_cell, x01_2:=NA]
out_dt_bis[x11_1<min_cell, x11_1:=NA]
out_dt_bis[x11_2<min_cell, x11_2:=NA]
# make figure
for(trt in 1:4){
  P <- ggplot(out_dt_bis[treatment==trt])
  P <- P + ggtitle(paste0("Cells trt",trt))
  #P <- P + geom_line(aes(x=time, y=N_1, color="Total", linetype="sp1"))
  #P <- P + geom_line(aes(x=time, y=N_2, color="Total", linetype="sp2"))
  P <- P + geom_line(aes(x=time, y=log10(x00_1), color="susceptible", linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=log10(x00_2), color="susceptible", linetype="sp2"))
  P <- P + geom_line(aes(x=time, y=log10(x10_1), color="single_ARG_A", linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=log10(x10_2), color="single_ARG_A", linetype="sp2"))
  P <- P + geom_line(aes(x=time, y=log10(x01_1), color="single_ARG_B", linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=log10(x01_2), color="single_ARG_B", linetype="sp2"))
  P <- P + geom_line(aes(x=time, y=log10(x11_1), color="double_ARG", linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=log10(x11_2), color="double_ARG", linetype="sp2"))
  P <- P + geom_line(aes(x=time, y=atb1+log10(min_cell)-1.5, color="DrugA"))
  P <- P + geom_line(aes(x=time, y=atb2+log10(min_cell)-1.9, color="DrugB"))
  P <- P + scale_y_continuous(name="cells", limits=c(log10(min_cell)-1,10), breaks=seq(log10(min_cell),10,2))
  #P <- P + scale_x_continuous(name="generation", limits=c(0,16E3), breaks=seq(0,16E3,4E3), expand=c(0,0))
  P <- P + scale_x_continuous(name="generation", limits=c(0,16E3), expand=c(0,0), breaks=seq(0,simLength,4000), labels=c("0", paste0(seq(4,simLength/1E3,4),"E3")))
  P <- P + scale_color_manual(values=color_pop, name="")
  P <- P + scale_linetype_discrete(name="", labels=c("sp1", "sp2"))
  P <- P + theme_minimal()
  P <- P + theme(panel.grid.minor.y=element_line(color="white")) # remove minor y grid
  P <- P + geom_hline(yintercept=log10(min_cell), linetype = "solid", color = "grey30", size=0.3)
  P <- P + geom_point(data=data.frame(x=seq(0,16000,4000),y=log10(min_cell)), aes(x=x, y=y), shape=3, color="grey30")
  P <- P + geom_segment(x=0,y=log10(min_cell),xend=0,yend=10,color='grey30')
  P <- P + geom_point(data=data.frame(x=0,y=seq(log10(min_cell),10,2)), aes(x=x, y=y), shape=3, color="grey30")
  plot(P)
  assign(paste0("P_Cells_trt",trt), P)
  ggsave(paste0("Cell_dynamics_trt",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}

### Fuse plots ###
## N1 ##
for(i in 1:4){
  P <- get(paste0("P_N1_trt",i))
  # remove y axis titles
  P <- P + theme(axis.title.y=element_blank())
  # x axis text size
  P <- P + theme(axis.title.x=element_text(size=9))
  assign(paste0("P_N1_trt",i), P)
}
P <- ggarrange(P_N1_trt1, P_N1_trt2, P_N1_trt3, P_N1_trt4,
               labels=c("A","B","C","D"), ncol=4, nrow=1,
               common.legend=TRUE, legend="right")
plot(P)
ggsave(paste0("N1_heatmaps_param1_all_trt.tiff"), width=fig_specs_4heatmaps$width, height=fig_specs_4heatmaps$height, units="cm", path=fig_specs_4heatmaps$path, dpi=fig_specs_4heatmaps$dpi)

## Cell dynamics ##
for(i in 1:4){
  P <- get(paste0("P_Cells_trt",i))
  # remove y axis titles
  P <- P + theme(axis.title.y=element_blank())
  # x axis text size
  P <- P + theme(axis.title.x=element_text(size=9))
  assign(paste0("P_Cells_trt",i), P)
}
P <- ggarrange(P_Cells_trt1, P_Cells_trt2, P_Cells_trt3, P_Cells_trt4,
               labels=c("A","B","C","D"), ncol=4, nrow=1,
               common.legend=TRUE, legend="right")
plot(P)
ggsave(paste0("Cells=f(t)_param1_all_trt.tiff"), width=fig_specs_4heatmaps$width, height=fig_specs_4heatmaps$height, units="cm", path=fig_specs_4heatmaps$path, dpi=fig_specs_4heatmaps$dpi)




### Other heatmaps not in paper ###
# N_1 / (N_1 + N_2)
P <- ggplot(out_dt_last_steps[time==simLength])
P <- P + ggtitle(paste0("N_1 / (N_1 + N_2)"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_N_1 / (mean_N_1 + mean_N_2)), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette_compete, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_N_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# N_2 / (N_1 + N_2)
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("N_2 / (N_1 + N_2)"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_N_2 / (mean_N_1 + mean_N_2)), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette_compete, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_N_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x00_1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x00_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x00_1/mean_N_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x00_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x00_2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x00_2/mean_N_2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x00_2/mean_N_2), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x00_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x11_1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_1/mean_N_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x11_2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_2/mean_N_2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_2/mean_N_2), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x01_1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x01_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x01_1/mean_N_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x01_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x01_2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x01_2/mean_N_2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x01_2/mean_N_2), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x01_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Arc distance m_1
P <- ggplot(out_dt_last_steps[time==simLength & x11_1>1E-6])
P <- P + ggtitle(paste0("m_1 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette_reverse)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_m_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Pab m_1
Es <- as.double(muvar(1, 50)[1])
P <- ggplot(out_dt_last_steps[time==simLength & x11_1>1E-6])
P <- P + ggtitle(paste0("pab_1 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=pab_1/Es), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_pab_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Arc distance m_2
P <- ggplot(out_dt_last_steps[time==simLength & x11_1>1E-6])
P <- P + ggtitle(paste0("m_2 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_2), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_m_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Pab m_1
Es <- as.double(muvar(1, 50)[1])
P <- ggplot(out_dt_last_steps[time==simLength & x11_1>1E-6])
P <- P + ggtitle(paste0("pab_2 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=pab_2/Es), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_pab_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna11_1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna11_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna11_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna11_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna11_2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna11_2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna11_2), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna11_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna01_1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna01_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna01_1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna01_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna01_2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna01_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna01_2), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna01_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)





#### Supplementary Figure 1 - example dynamics ####

# load heatmap data DNAdelation/HorizontalTransfer 1patch/1species

# Cell dynamics - genotypes=f(t) #
# subset oen simulation from heatmap data
out_dt_bis <- out_dt[param_1==1E-4][param_2==1E-1]
min_cell <- 1E0
out_dt_bis[x10_1<min_cell, x10_1:=NA]
out_dt_bis[x01_1<min_cell, x01_1:=NA]
out_dt_bis[x10_1<min_cell, x10_1:=NA]
out_dt_bis[x11_1<min_cell, x11_1:=NA]
# make figure
color_pop <- c(susceptible="green2",single_ARG_A="blue2",single_ARG_B='red2',double_ARG="magenta2",DrugA="blue",DrugB="red")
for(trt in 1:4){
  P <- ggplot(out_dt_bis[treatment==trt])
  P <- P + ggtitle(paste0("Cells trt",trt))
  #P <- P + geom_line(aes(x=time, y=N_1, color="Total", linetype="sp1"))
  #P <- P + geom_line(aes(x=time, y=N_2, color="Total", linetype="sp2"))
  P <- P + geom_line(aes(x=time, y=log10(x00_1), color="susceptible"))#, linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=log10(x10_1), color="single_ARG_A"))#, linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=log10(x01_1), color="single_ARG_B"))#, linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=log10(x11_1), color="double_ARG"))#, linetype="sp1"))
  P <- P + geom_line(aes(x=time, y=atb1+log10(min_cell)-1.8, color="DrugA"))
  P <- P + geom_line(aes(x=time, y=atb2+log10(min_cell)-2.3, color="DrugB"))
  P <- P + scale_y_continuous(name="cells", limits=c(log10(min_cell)-1.5,10), breaks=seq(log10(min_cell),10,2))
  P <- P + scale_x_continuous(name="generation", limits=c(0,16E3), expand=c(0,0), breaks=seq(0,simLength,4000), labels=c("0", paste0(seq(4,simLength/1E3,4),"E3")))
  P <- P + scale_color_manual(values=color_pop, name="")
  P <- P + theme_minimal()
  P <- P + theme(panel.grid.minor.y=element_line(color="white")) # remove minor y grid
  P <- P + geom_hline(yintercept=log10(min_cell), linetype = "solid", color = "grey30", size=0.3)
  P <- P + geom_point(data=data.frame(x=seq(0,16000,4000),y=log10(min_cell)), aes(x=x, y=y), shape=3, color="grey30")
  P <- P + geom_segment(x=0,y=log10(min_cell),xend=0,yend=10,color='grey30')
  P <- P + geom_point(data=data.frame(x=0,y=seq(log10(min_cell),10,2)), aes(x=x, y=y), shape=3, color="grey30")
  plot(P)
  assign(paste0("P_Cells_trt",trt), P)
  ggsave(paste0("Cell_dynamics_trt",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}

# Conditional Pab=f(t) log #
Es <- as.double(muvar(1, 50)[1])
color_pop <- c(pab="orange2",DrugA="blue",DrugB="red")
for(trt in 1:4){
  P <- ggplot(out_dt_bis[x11_1>1E-6][treatment==trt])
  P <- P + ggtitle(paste0("pab_cond_log=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=log10(pab_1 / Es), color="pab"))
  P <- P + geom_line(aes(x=time, y=atb1-3.3, color="DrugA"))
  P <- P + geom_line(aes(x=time, y=atb2-3.4, color="DrugB"))
  P <- P + scale_color_manual(values=color_pop, name="")
  P <- P + scale_x_continuous(name="generation", limits=c(0,16E3), breaks=seq(0,simLength,4000), labels=c("0", paste0(seq(4,simLength/1E3,4),"E3")), expand=c(0,0)) # expand removes margin
  P <- P + scale_y_continuous(name="log pab cond", limits=c(-2.4,0), breaks=seq(-2,0,0.5))
  P <- P + labs(color="param_1")
  P <- P + theme_minimal()
  P <- P + theme(panel.grid.minor.y=element_line(color="white")) # remove minor y grid
  P <- P + geom_hline(yintercept = -2, linetype = "solid", color = "grey30", size=0.3)
  P <- P + geom_point(data=data.frame(x=seq(0,16000,4000),y=-2), aes(x=x, y=y), shape=3, color="grey30")
  P <- P + geom_segment(x=0,y=-2,xend=0,yend=0,color='grey30')
  P <- P + geom_point(data=data.frame(x=0,y=seq(-2,0,0.5)), aes(x=x, y=y), shape=3, color="grey30")
  plot(P)
  assign(paste0("P_pab_trt",trt), P)
  ggsave(paste0("Pab_dynamics_trt",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}

# # DNA=f(t) - ratios
# P <- ggplot(out_dt)
# P <- P + geom_line(aes(x=time, y=log10(Dna01_1), color="Dna01"))
# P <- P + geom_line(aes(x=time, y=log10(Dna10_1), color="Dna10"), linetype="dashed")
# P <- P + geom_line(aes(x=time, y=log10(Dna11_1), color="Dna11"))
# #P <- P + scale_y_continuous(name="DNA (log)", limits=c(-4,0), breaks=seq(-4,0,1))
# P <- P + scale_color_manual(values=c("blue2",'red2',"magenta2"))
# P <- P + labs(color="")
# P <- P + theme_light()
# plot(P)
# ggsave("Dna=f(t)_1sp.tiff", width=fig_specs$width, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)



### Fuse plots ###
## Cell dynamics ##
for(i in 1:4){
  P <- get(paste0("P_Cells_trt",i))
  # remove y axis titles
  P <- P + theme(axis.title.y=element_blank())
  # x axis text size
  P <- P + theme(axis.title.x=element_text(size=9))
  assign(paste0("P_Cells_trt",i), P)
}
P <- ggarrange(P_Cells_trt1, P_Cells_trt2, P_Cells_trt3, P_Cells_trt4,
               labels=c("A","B","C","D"), ncol=4, nrow=1,
               common.legend=TRUE, legend="right")
plot(P)
ggsave(paste0("Cells=f(t)_all_trt.tiff"), width=fig_specs_4heatmaps$width, height=fig_specs_4heatmaps$height, units="cm", path=fig_specs_4heatmaps$path, dpi=fig_specs_4heatmaps$dpi)

## Pab cond ##
for(i in 1:4){
  P <- get(paste0("P_pab_trt",i))
  # remove y axis titles
  P <- P + theme(axis.title.y=element_blank())
  # x axis text size
  P <- P + theme(axis.title.x=element_text(size=9))
  assign(paste0("P_pab_trt",i), P)
}
P <- ggarrange(P_pab_trt1, P_pab_trt2, P_pab_trt3, P_pab_trt4,
               labels=c("A","B","C","D"), ncol=4, nrow=1,
               common.legend=TRUE, legend="right")
plot(P)
ggsave(paste0("pab=f(t)_all_trt.tiff"), width=fig_specs_4plots$width, height=fig_specs_4plots$height, units="cm", path=fig_specs_4plots$path, dpi=fig_specs_4plots$dpi)



#### 2env 1species ####
nbPop <- 1
nbEnvironment <- 2
out_dt <- setDT(as.data.frame(out_parallel))
for(i in 1:nbPop){ # add columns
  print(paste0("pop ",i))
  for(j in 1:nbEnvironment){
    print(paste0("env ",j))
    out_dt[, paste0("m_env",j,"_sp",i) := get(paste0("NEf1_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i))] # mean arc distance
    out_dt[, paste0("S_env",j,"_sp",i) := get(paste0("NEf2_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i)) - get(paste0("m_env",j,"_sp",i))^2] # variance arc distance
    out_dt[, paste0("N_env",j,"_sp",i) := get(paste0("x00_env",j,"_sp",i)) + get(paste0("x01_env",j,"_sp",i)) + get(paste0("x10_env",j,"_sp",i)) + get(paste0("x11_env",j,"_sp",i))] # size population species i
    
    out_dt[, paste0("ms_env",j,"_sp",i) := paste0(get(paste0("m_env",j,"_sp",i)),";",get(paste0("S_env",j,"_sp",i)))]
    out_dt[, paste0("alpha_env",j,"_sp",i) := sapply(get(paste0("ms_env",j,"_sp",i)), function(x){
      m=as.numeric(strsplit(x,";")[[1]][1]);
      s=as.numeric(strsplit(x,";")[[1]][2]);
      return(absolve(m, s)[1])})]
    out_dt[, paste0("beta_env",j,"_sp",i) := sapply(get(paste0("ms_env",j,"_sp",i)), function(x){
      m=as.numeric(strsplit(x,";")[[1]][1]);
      s=as.numeric(strsplit(x,";")[[1]][2]);
      return(absolve(m, s)[2])})]
    out_dt[, paste0("ab_env",j,"_sp",i) := paste0(get(paste0("alpha_env",j,"_sp",i)),";",get(paste0("beta_env",j,"_sp",i)))]
    out_dt[, paste0("pab_env",j,"_sp",i) := sapply(get(paste0("ab_env",j,"_sp",i)), function(x){
      a=as.numeric(strsplit(x,";")[[1]][1]);
      b=as.numeric(strsplit(x,";")[[1]][2]);
      return(comobg(50,a,b))})]
  }
}
## check ##
out_dt[1:2]
## Mean last time steps ##
colnames(out_dt)
out_dt_last_steps <- out_dt[time >= (simLength-4000)] # subset last 2000 time steps
out_dt_last_steps[, mean_N_env1_sp1 := mean(N_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x00_env1_sp1 := mean(x00_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x01_env1_sp1 := mean(x01_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x10_env1_sp1 := mean(x10_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x11_env1_sp1 := mean(x11_env1_sp1), by=param_set]
out_dt_last_steps[, mean_Dna01_env1_sp1 := mean(Dna01_env1_sp1), by=param_set]
out_dt_last_steps[, mean_Dna10_env1_sp1 := mean(Dna10_env1_sp1), by=param_set]
out_dt_last_steps[, mean_Dna11_env1_sp1 := mean(Dna11_env1_sp1), by=param_set]
## Params to log (if needed) ##
out_dt_last_steps[, param_1 := log10(param_1)]
out_dt_last_steps[, param_2 := log10(param_2)]


### Heatmaps ###
# x00
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x00_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x00_env1_sp1/mean_N_env1_sp1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x00_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x11
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_env1_sp1/mean_N_env1_sp1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0.9,1), breaks=seq(0,1,0.05)) # c5 t2E4
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0.95,1), breaks=seq(0,1,0.01))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x01
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x01_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x01_env1_sp1/mean_N_env1_sp1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x01_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Arc distance m
P <- ggplot(out_dt_last_steps[time==simLength & x11_env1_sp1>1E-6])
P <- P + ggtitle(paste0("m_env1_sp1 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_env1_sp1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,0.25), breaks=seq(0,0.25,0.05))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_m_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# conditional Pab m
Es <- as.double(muvar(1, 50)[1])
P <- ggplot(out_dt_last_steps[time==simLength & x11_env1_sp1>1E-6])
P <- P + ggtitle(paste0("Pab_cond_env1_sp1 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=pab_env1_sp1/Es), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,0.22), breaks=seq(0,0.25,0.05))
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,0.25), breaks=seq(0,0.25,0.05))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_pab_cond_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna11
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna11_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna11_env1_sp1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna11_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna01
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna01_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna01_env1_sp1), alpha=alpha_heatmap)
P <- P + scale_fill_gradientn(colours=color_palette)#, limits=c(0,1), breaks=seq(0,1,0.2))
P <- P + scale_x_continuous(name="param1")
P <- P + scale_y_continuous(name="param2")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna01_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)



#### LHS 1env 1sp ####

# Get result file
out_dt <- setDT(as.data.frame(out_parallel))
# Add columns
nbPop <- 1
setnames(out_dt, "rUpt","rUpt_1")
for(i in 1:nbPop){
  out_dt[, paste0("m_",i) := get(paste0("NEf1_",i)) / get(paste0("x11_",i))] # mean arc distance
  out_dt[, paste0("S_",i) := get(paste0("NEf2_",i)) / get(paste0("x11_",i)) - get(paste0("m_",i))^2] # variance arc distance
  out_dt[, paste0("N_",i) := get(paste0("x00_",i)) + get(paste0("x01_",i)) + get(paste0("x10_",i)) + get(paste0("x11_",i))] # size population species 1
  out_dt[, paste0("ms_",i) := paste0(get(paste0("m_",i)),";",get(paste0("S_",i)))]
  out_dt[, paste0("alpha_",i) := sapply(get(paste0("ms_",i)), function(x){
    m=as.numeric(strsplit(x,";")[[1]][1]);
    s=as.numeric(strsplit(x,";")[[1]][2]);
    return(absolve(m, s)[1])})]
  out_dt[, paste0("beta_",i) := sapply(get(paste0("ms_",i)), function(x){
    m=as.numeric(strsplit(x,";")[[1]][1]);
    s=as.numeric(strsplit(x,";")[[1]][2]);
    return(absolve(m, s)[2])})]
  out_dt[, paste0("ab_",i) := paste0(get(paste0("alpha_",i)),";",get(paste0("beta_",i)))]
  out_dt[, paste0("pab_",i) := sapply(get(paste0("ab_",i)), function(x){
    a=as.numeric(strsplit(x,";")[[1]][1]);
    b=as.numeric(strsplit(x,";")[[1]][2]);
    return(comobg(50,a,b))})]
  out_dt[, Es := sapply(gamma, function(x){as.double(muvar(1, x)[1])})]
  out_dt[, paste0("pab_cond_",i) := get(paste0("pab_",i)) / Es]
}
setDT(lhs)
#lhs[, Es := as.double(muvar(1, 50)[1])] # gamma=50
#lhs[, stressPower := deathPeriod * deathIncrease / lifePeriod]
out_dt[1:2]

# Mean last time steps #
out_dt_last_steps <- out_dt[time >= (simLength-4000)] # subset last time steps
out_dt_last_steps[, mean_N_1 := mean(N_1), by=param_set]
out_dt_last_steps[, mean_x00_1 := mean(x00_1), by=param_set]
out_dt_last_steps[, mean_x01_1 := mean(x01_1), by=param_set]
out_dt_last_steps[, mean_x11_1 := mean(x11_1), by=param_set]
out_dt_last_steps[, mean_Dna11_1 := mean(Dna11_1), by=param_set]
out_dt_last_steps[, mean_Dna01_1 := mean(Dna01_1), by=param_set]
out_dt_last_steps[, mean_pab_1 := mean(pab_1), by=param_set]
out_dt_last_steps[, mean_pab_cond_1 := mean(pab_cond_1), by=param_set]

# Data subsets #
cluster_set <- out_dt[time==simLength][m_1<0.24][x11_1>1E-6][,param_set]
scattered_set <- out_dt[time==simLength][m_1>0.24][x11_1>1E-6][,param_set]

strong_cluster_set <- out_dt[time==simLength][m_1<0.15][x11_1>1E-6][,param_set]
weak_cluster_set <- out_dt[time==simLength][m_1>0.15 & m_1<0.20][x11_1>1E-6][,param_set]
mild_cluster_set <- out_dt[time==simLength][m_1>0.20 & m_1<0.24][x11_1>1E-6][,param_set]

extinction_set <- out_dt[x11_1<1E-6 & time==simLength][,param_set]
failed_simu <- seq(1,nrow(lhs),1)[-out_dt[time==simLength][,param_set]]


# Parameter distribution - classic #
P <- ggplot(as.data.frame(lhs[scattered_set,]))
P <- P + ggtitle(paste0("scattered sets, n=", length(scattered_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[cluster_set,]))
P <- P + ggtitle(paste0("cluster sets, n=", length(cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[weak_cluster_set,]))
P <- P + ggtitle(paste0("weak cluster sets, n=", length(weak_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[mild_cluster_set,]))
P <- P + ggtitle(paste0("mild cluster sets, n=", length(mild_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[strong_cluster_set,]))
P <- P + ggtitle(paste0("strong cluster sets, n=", length(strong_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[failed_simu,]))
P <- P + ggtitle(paste0("failed sets, n=", length(failed_simu),"/",nrow(lhs)))

P <- P + geom_boxplot(aes(x=1, y=(log10(rMob)-log10(min(rMob)))/(log10(max(rMob))-log10(min(rMob)))))
P <- P + geom_boxplot(aes(x=2, y=(log10(rDup)-log10(min(rDup)))/(log10(max(rDup))-log10(min(rDup)))))
P <- P + geom_boxplot(aes(x=3, y=(log10(rUpt)-log10(min(rUpt)))/(log10(max(rUpt))-log10(min(rUpt)))))
P <- P + geom_boxplot(aes(x=4, y=(rClear-min(rClear))/(max(rClear)-min(rClear))))
P <- P + geom_boxplot(aes(x=5, y=(log10(fcost)-log10(min(fcost)))/(log10(max(fcost))-log10(min(fcost)))))
P <- P + geom_boxplot(aes(x=6, y=(deathBasal-min(deathBasal))/(max(deathBasal)-min(deathBasal))))
P <- P + geom_boxplot(aes(x=7, y=(deathIncrease-min(deathIncrease))/(max(deathIncrease)-min(deathIncrease))))
P <- P + geom_boxplot(aes(x=8, y=(deathPeriod-min(deathPeriod))/(max(deathPeriod)-min(deathPeriod))))
P <- P + geom_boxplot(aes(x=9, y=(lifePeriod-min(lifePeriod))/(max(lifePeriod)-min(lifePeriod))))
P <- P + scale_y_continuous(name="range ratio")
P <- P + scale_x_continuous(name="",breaks=1:9,labels=c("rMob","rDup","rUpt","rClear","fcost","deathBasal","deathIncrease","deathPeriod","lifePeriod"))
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
P <- P + labs(color="Parameter")
plot(P)

ggsave("lhs_1env_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)







#### LHS 2env 1sp ####
nbPop <- 1
nbEnvironment <- 2
out_dt <- setDT(as.data.frame(out_parallel))
for(i in 1:nbPop){ # add columns
  print(paste0("pop ",i))
  for(j in 1:nbEnvironment){
    print(paste0("env ",j))
    out_dt[, paste0("m_env",j,"_sp",i) := get(paste0("NEf1_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i))] # mean arc distance
    out_dt[, paste0("S_env",j,"_sp",i) := get(paste0("NEf2_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i)) - get(paste0("m_env",j,"_sp",i))^2] # variance arc distance
    out_dt[, paste0("N_env",j,"_sp",i) := get(paste0("x00_env",j,"_sp",i)) + get(paste0("x01_env",j,"_sp",i)) + get(paste0("x10_env",j,"_sp",i)) + get(paste0("x11_env",j,"_sp",i))] # size population species i
    
    out_dt[, paste0("ms_env",j,"_sp",i) := paste0(get(paste0("m_env",j,"_sp",i)),";",get(paste0("S_env",j,"_sp",i)))]
    out_dt[, paste0("alpha_env",j,"_sp",i) := sapply(get(paste0("ms_env",j,"_sp",i)), function(x){
      m=as.numeric(strsplit(x,";")[[1]][1]);
      s=as.numeric(strsplit(x,";")[[1]][2]);
      return(absolve(m, s)[1])})]
    out_dt[, paste0("beta_env",j,"_sp",i) := sapply(get(paste0("ms_env",j,"_sp",i)), function(x){
      m=as.numeric(strsplit(x,";")[[1]][1]);
      s=as.numeric(strsplit(x,";")[[1]][2]);
      return(absolve(m, s)[2])})]
    out_dt[, paste0("ab_env",j,"_sp",i) := paste0(get(paste0("alpha_env",j,"_sp",i)),";",get(paste0("beta_env",j,"_sp",i)))]
    out_dt[, paste0("pab_env",j,"_sp",i) := sapply(get(paste0("ab_env",j,"_sp",i)), function(x){
      a=as.numeric(strsplit(x,";")[[1]][1]);
      b=as.numeric(strsplit(x,";")[[1]][2]);
      return(comobg(50,a,b))})]
    
    Es <- as.double(muvar(1, 50)[1])
    out_dt[, paste0("pab_cond_env",j,"_sp",i) := get(paste0("pab_env",j,"_sp",i)) / Es]
    
  }
}

## check ##
out_dt[1:2]
## Mean last time steps ##
colnames(out_dt)
out_dt_last_steps <- out_dt[time >= (simLength-4000)] # subset last 2000 time steps
out_dt_last_steps[, mean_N_env1_sp1 := mean(N_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x00_env1_sp1 := mean(x00_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x01_env1_sp1 := mean(x01_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x10_env1_sp1 := mean(x10_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x11_env1_sp1 := mean(x11_env1_sp1), by=param_set]
out_dt_last_steps[, mean_Dna01_env1_sp1 := mean(Dna01_env1_sp1), by=param_set]
out_dt_last_steps[, mean_Dna10_env1_sp1 := mean(Dna10_env1_sp1), by=param_set]
out_dt_last_steps[, mean_Dna11_env1_sp1 := mean(Dna11_env1_sp1), by=param_set]










