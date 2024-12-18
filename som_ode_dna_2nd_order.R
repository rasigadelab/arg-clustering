library(deSolve)
library(data.table)
library(ggplot2)
library(lhs)
library(foreach)
library(doParallel)

#source("210407_somfuncs.R")
source("210929_somfuncs.R")
source("220209_gsom2.R")

#### Make cluster for parallel computing ####
# Make cluster
detectCores() # number of cores of the computer
cl <- makeCluster(7) # number of cores allowed to parallel computation
registerDoParallel(cl)
# Get info
getDoParWorkers()
getDoParName()
getDoParVersion()
# Stop cluster
stopImplicitCluster()
stopCluster(cl)

#### Simulation ####

## PARAMETERS ##
parsList <- list()
## Population parameters ##
b <- log(2) # Base birth rate, if =log(2) -> doubling time = time step
fcost <- 0.03 # Fitness cost per ARG
parsList["b00"] <- b # Growth rate
parsList["b01"] <- b - fcost
parsList["b10"] <- b - fcost
parsList["b11"] <- b - 2*fcost
parsList["K"] <- 1E10 # Carrying capacity
parsList["deathBasal"] <- 0.1 # basal cell death rate
parsList["x00_niche"] <- 0 # constant x00 population (interface niche of wild type cells)
parsList["WTin"] <- 0 # Wild type input (replace part of the population by WT cells)
## 2 patch model ##
parsList["K1"] <- 1E10 # Carrying capacity environment 1 (for 2patch model)
parsList["K2"] <- 1E10 # Carrying capacity environment 2 (for 2patch model)
parsList["G"] <- 0.05 # migration rate (for 2patch model)
## HGT/DNA Parameters ##
parsList["gamma"] <- 50 # manage s distribution (size transferred DNA)
parsList["Es"]  <- as.double(muvar(1, parsList$gamma)[1]) # size transferred DNA
parsList["rMob"] <- 0.001 # gene mobilization (excision)
parsList["rDup"] <- 0.01 # gene duplication
parsList["rUpt"] <- parsList$rDup / 50
parsList["Dna_cst"] <- 1 # limits DNA uptake when eDNA low
#parsList["rUpt"] <- (parsList[["rMob"]] + parsList[["rDup"]])/10 # eDNA uptake
parsList["rClear"] <- .1 # eDNA clearance rate
## Antibiotic treatment ##
parsList["deathPeriod"] <- 2000 # duration of the treatment periods
parsList["deathIncrease"] <- 0.1 # death rate increase during treatments (added to basal death rate for susceptible cells)
parsList["treatmentType"] <- 2 # 1: no atb, 2: cst atb, 3:synchro atb, 4:altern atb, 5: synchro atb / uneven periods
parsList["lifePeriod"] <- 200 # if treatmentType=3: recovery time between treatments
rm(b,fcost)

## VARIABLES ##
yiniList <- list()
occup <- .5 # Initial cell occupancy (K ratio)
yiniList["x00"] <- parsList$K * occup * 0.25
yiniList["x01"] <- parsList$K * occup * 0.25
yiniList["x10"] <- parsList$K * occup * 0.25
yiniList["x11"] <- parsList$K * occup * 0.25
alpha <- 1
beta <- 1 # initial arc distance (1: uniform arc distance distribution, >1: gene clustering)
Ef1 <- as.double(muvar(alpha, beta))[1] # 1st moment -> mean
Ef2 <- as.double(muvar(alpha, beta))[2] + Ef1^2 # 2nd moment about zero -> allow to retrieve variance (sum squarred error from zero and not from the mean)
yiniList["NEf1"] <- as.double(yiniList$x11 * Ef1)
yiniList["NEf2"] <- as.double(yiniList$x11 * Ef2)
yiniList["Dna01"] <- 0 # extracellular DNA
yiniList["Dna10"] <- 0 # extracellular DNA
yiniList["Dna11"] <- 0 # extracellular DNA
rm(occup,Ef1,Ef2,alpha,beta)

## Reshaping parsList and yiniList ##
# 1 species #
names(parsList) <- paste0(names(parsList), "_1") # add species suffix 
pars <- c(nSpecies = 1, unlist(parsList)) # unlist for ode function
names(yiniList) <- paste0(names(yiniList), "_1") # add species suffix
yini <- c(unlist(yiniList))# ode() takes vector argument not list

# 2 species #
parsSp1 <- parsList
names(parsSp1) <- paste0(names(parsSp1), "_1") # add species suffix
parsSp2 <- parsList
names(parsSp2) <- paste0(names(parsSp2), "_2") # add species suffix
yiniSp1 <- yiniList
names(yiniSp1) <- paste0(names(yiniSp1), "_1") # add species suffix
yiniSp2 <- yiniList
names(yiniSp2) <- paste0(names(yiniSp2), "_2") # add species suffix
# specify gene cluster species2
alpha <- 1
beta <- 1E3 # 1:scattered genes, 1E4:clustered genes
Ef1 <- as.double(muvar(alpha, beta))[1] # 1st moment -> mean
Ef2 <- as.double(muvar(alpha, beta))[2] + Ef1^2 # 2nd moment about zero -> allow to retrieve variance (sum squarred error from zero and not from the mean)
yiniSp2["NEf1_2"] <- as.double(yiniSp2$x11 * Ef1)
yiniSp2["NEf2_2"] <- as.double(yiniSp2$x11 * Ef2)
rm(Ef1,Ef2,alpha,beta)
# specify non-transformable species1
parsSp1["rDup_1"] <- 0
parsSp1["rMob_1"] <- 0
parsSp1["rUpt_1"] <- 0
# initialize pars & yini
pars <- c(nSpecies = 2, unlist(parsSp1), unlist(parsSp2)) # unlist for ode function
yini <- c(unlist(yiniSp1), unlist(yiniSp2))

# 2 patches #
names(parsList) <- paste0(names(parsList), "_sp1") # add species suffix
pars <- c(nSpecies = 1, unlist(parsList)) # unlist for ode function
yiniList_env1 <- yiniList
yiniList_env2 <- yiniList
names(yiniList_env1) <- paste0(names(yiniList_env1), "_env1") # add environment suffix
names(yiniList_env2) <- paste0(names(yiniList_env2), "_env2") # add environment suffix
yiniList <- c(yiniList_env1, yiniList_env2)
names(yiniList) <- paste0(names(yiniList), "_sp1") # add species suffix
yini <- c(unlist(yiniList))# ode() takes vector argument not list
rm(yiniList_env1, yiniList_env2)


## SIMULATION ##
simLength <- 40000 # simulation duration
simLength <- 16000
simLength <- 5000
simStep   <- 1 # 1 step = doubling time if b=log(2)
times <- seq(0, simLength, by = simStep) # simulation time points
max_dt <- 1
max_dt <- 1E-1 # max time step used by the ode solver
max_dt <- 1E-2
max_dt <- 1E-3

out <- lsoda(yini, times, LVcompet, pars,
             rtol = 1E-6, # relative error tolerance, either a scalar or an array as long as y (default 1E-6)
             atol = 1E-6, # absolute error tolerance, either a scalar or an array as long as y (default 1E-6)
             verbose = FALSE,
             hmin = 0, # minimum value of the integration stepsize
             hmax = max_dt, # maximum value of the integration stepsize (default NULL)
             hini = 0, # initial step size to be attempted; if 0, the initial step size is determined by the solver
             maxsteps = 5000)

out <- lsoda(yini, times, LVcompet_wtinput, pars, DEBUG = TRUE, # x00 input (input of WT cells)
             verbose = FALSE,
             hmin = 0, 
             hmax = max_dt, 
             hini = 0, 
             maxsteps = 5000)

out <- lsoda(yini, times, LVcompet_2patch, pars, DEBUG = TRUE, # 2 environments
             verbose = FALSE,
             hmin = 0, 
             hmax = max_dt, 
             hini = 0, 
             maxsteps = 5000)

out <- lsoda(yini, times, LVcompet_2patchReg, pars, DEBUG = TRUE, # 2 environments - DNA uptake regulation
             verbose = FALSE,
             hmin = 0, 
             hmax = max_dt, 
             hini = 0, 
             maxsteps = 5000)

out <- lsoda(yini, times, LVcompet_uptReg, pars, DEBUG = TRUE, # DNA uptake = rUpt * DNAtot / (DNAtot + Dna_cst)
             verbose = FALSE,
             hmin = 0, 
             hmax = max_dt, 
             hini = 0, 
             maxsteps = 5000)

out <- lsoda(yini, times, LVcompet_sharedDNA, pars, DEBUG=TRUE, # shared eDNA pool
             verbose = FALSE,
             hmin = 0, 
             hmax = max_dt, 
             hini = 0, 
             maxsteps = 5000)

out <- lsoda(yini, times, LVcompet_loss, pars, DEBUG=TRUE, # gene loss instead of excision
             verbose = FALSE,
             hmin = 0, 
             hmax = max_dt, 
             hini = 0, 
             maxsteps = 5000)

out <- lsoda(yini, times, LVcompet_2patchReg_loss, pars, DEBUG=TRUE, # gene loss instead of excision
             verbose = FALSE,
             hmin = 0, 
             hmax = max_dt, 
             hini = 0, 
             maxsteps = 5000)

out <- lsoda(yini, times, LVcompet_2patchReg_loss_spaceBis, pars, DEBUG=TRUE, # gene loss instead of excision
             verbose = FALSE,
             hmin = 0, 
             hmax = max_dt, 
             hini = 0, 
             maxsteps = 5000)

deSolve::diagnostics.deSolve(out)
save(out, file="~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures/out_parallel.Rdata")


#### Plots - 1 env ####
figure_path <- "./Figures"
fig_specs <- list(
  "width"=8,
  "heigth"=6,
  "dpi"=150
)
nbPop <- 1
nbPop <- 2
out_dt <- setDT(as.data.frame(out))
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

## 1 species ##
# Genotypes=f(t) - ratios
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=log10(x00_1), color="x00"))
P <- P + geom_line(aes(x=time, y=log10(x01_1), color="x01"))
P <- P + geom_line(aes(x=time, y=log10(x10_1), color="x10"), linetype="dashed")
P <- P + geom_line(aes(x=time, y=log10(x11_1), color="x11"))
#P <- P + scale_y_continuous(name="genotypes (log)", limits=c(6, log10(parsList$K_1)), breaks=seq(4,10,2))
P <- P + scale_y_continuous(name="genotypes (log)", limits=c(6,10),breaks=seq(0,10,2))
P <- P + scale_color_manual(values=c("green2","blue2",'red2',"magenta2"))
P <- P + labs(color="")
P <- P + theme_light()
plot(P)
ggsave("Genotypes=f(t)_1sp.tiff", width=fig_specs$width, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)

# DNA=f(t) - ratios
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=log10(Dna01_1), color="Dna01"))
P <- P + geom_line(aes(x=time, y=log10(Dna10_1), color="Dna10"), linetype="dashed")
P <- P + geom_line(aes(x=time, y=log10(Dna11_1), color="Dna11"))
#P <- P + scale_y_continuous(name="DNA (log)", limits=c(-4,0), breaks=seq(-4,0,1))
P <- P + scale_color_manual(values=c("blue2",'red2',"magenta2"))
P <- P + labs(color="")
P <- P + theme_light()
plot(P)
ggsave("Dna=f(t)_1sp.tiff", width=fig_specs$width, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)

# m=f(t)
P <- ggplot(out_dt[x11_1>1E-6])
P <- P + geom_line(aes(x=time, y=m_1))
P <- P + scale_x_continuous(name="time", limits=c(0,simLength), breaks=seq(0,simLength,5000))
P <- P + scale_y_continuous(name="mean arc distance", limits=c(0, 0.25), breaks=seq(0,0.25,0.05))
P <- P + labs(color="")
P <- P + theme_light()
plot(P)
ggsave("m=f(t)_1sp.tiff", width=fig_specs$width-2, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)

# Probability comobilization
Es <- as.double(muvar(1, 50)[1])
P <- ggplot(out_dt[x11_1>1E-6])
P <- P + ggtitle(paste0("pab_cond_log=f(t) trt",trt))
P <- P + geom_line(aes(x=time, y=log10(pab_1 / Es)))
#P <- P + geom_line(aes(x=time, y=atb1-3.15, color="atb1"))
#P <- P + geom_line(aes(x=time, y=atb2-3.1, color="atb2"))
#P <- P + scale_color_manual(values=color_lines)
P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
P <- P + scale_y_continuous(name="log pab cond")#, limits=c(-2.2,0), breaks=seq(-2,0,0.5))
P <- P + labs(color="param_1")
P <- P + theme_light()
plot(P)
ggsave("pab_cond_log=f(t).tiff", width=fig_specs$width-2, height=fig_specs$heigth, units="cm", path=figure_path, dpi=fig_specs$dpi)




# Check intermediate values #
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=pab, color="pab"))
P <- P + geom_line(aes(x=time, y=eta, color="eta"))
P <- P + geom_line(aes(x=time, y=Ef1, color="Ef1"))
P <- P + geom_line(aes(x=time, y=Ed1, color="Ed1"))
P <- P + geom_line(aes(x=time, y=Es1, color="Es1"))
P <- P + geom_line(aes(x=time, y=Eu1, color="Eu1"))

P <- P + geom_line(aes(x=time, y=Ef2, color="Ef2"))
P <- P + geom_line(aes(x=time, y=Ed2, color="Ed2"))
P <- P + geom_line(aes(x=time, y=Es2, color="Es2"))
P <- P + geom_line(aes(x=time, y=Eu2, color="Eu2"))

P <- P + geom_line(aes(x=time, y=dNf_, color="dNf_"))
P <- P + geom_line(aes(x=time, y=dNd_, color="dNd_"))
P <- P + geom_line(aes(x=time, y=dNs_, color="dNs_"))
P <- P + geom_line(aes(x=time, y=dNu_, color="dNu_"))

P <- P + theme_light()
plot(P)


## 2 species ##
# genotypes=f(t)
P <- ggplot(out_dt[param_1==1E-5 & param_2==1])
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=log10(N_1), color="Total", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(N_2), color="Total", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x00_1), color="WT", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x00_2), color="WT", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x01_1+x10_1), color="single-gene", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x01_2+x10_2), color="single-gene", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x11_1), color="double-gene", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x11_2), color="double-gene", linetype="sp2"))
P <- P + scale_y_continuous(name="cells")
P <- P + scale_color_discrete(name="")
P <- P + scale_linetype_discrete(name="", labels=c("sp1", "sp2"))
P <- P + theme_light()
plot(P)
ggsave("N&x=f(t)_2sp.tiff", width=15, height=8, units="cm", path=figure_path, dpi=150)
ggsave("N&x=f(t)_2sp_bis.tiff", width=15, height=8, units="cm", path=figure_path, dpi=150)

# genotypes=f(t) - debugging #
for(i in 1:7){
  p1 <- 1E-4
  p2 <- 1
  trt <- 4
  P <- ggplot(out_dt[param_1==p1 & param_2==p2][treatment==trt])
  P <- P + ggtitle(paste0("p1=",p1," p2=",p2," trt=",trt))
  if(i==1){
    P <- P + geom_line(aes(x=time, y=log10(N_1), color="Total1", linetype="sp1"))
    P <- P + geom_line(aes(x=time, y=log10(N_2), color="Total2", linetype="sp2"))
  }
  if(i==2){
    P <- P + geom_line(aes(x=time, y=log10(x00_1), color="WT1", linetype="sp1"))
    P <- P + geom_line(aes(x=time, y=log10(x00_2), color="WT2", linetype="sp2"))
  }
  if(i==3){
    P <- P + geom_line(aes(x=time, y=log10(x10_1), color="single-gene-A-1", linetype="sp1"))
    P <- P + geom_line(aes(x=time, y=log10(x10_2), color="single-gene-A-2", linetype="sp2"))
  }
  if(i==4){
    P <- P + geom_line(aes(x=time, y=log10(x01_1), color="single-gene-B-1", linetype="sp1"))
    P <- P + geom_line(aes(x=time, y=log10(x01_2), color="single-gene-B-2", linetype="sp2"))
  }
  if(i==5){
    P <- P + geom_line(aes(x=time, y=log10(x11_1), color="double-gene-1", linetype="sp1"))
    P <- P + geom_line(aes(x=time, y=log10(x11_2), color="double-gene-2", linetype="sp2"))
  }
  if(i==6){
    P <- P + geom_line(aes(x=time, y=log10(Dna10_1), color="DNA-A-1", linetype="sp1"))
    P <- P + geom_line(aes(x=time, y=log10(Dna10_2), color="DNA-A-2", linetype="sp2"))
  }
  if(i==7){
    P <- P + geom_line(aes(x=time, y=log10(Dna01_1), color="DNA-B-1", linetype="sp1"))
    P <- P + geom_line(aes(x=time, y=log10(Dna01_2), color="DNA-B-2", linetype="sp2"))
  }
  P <- P + scale_x_continuous(limit=c(0,8E3))
  P <- P + scale_y_continuous(name="cells", limits=c(0,10))
  P <- P + scale_color_discrete(name="")
  P <- P + scale_linetype_discrete(name="", labels=c("sp1", "sp2"))
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("Debug_fig_p1=",p1,"_p2=",p2,"_",i,".tiff"), width=15, height=8, units="cm", path=figure_path, dpi=150)
}




# arc distance = f(t)
P <- ggplot(out_dt[x11_1>1E-6]) # m value not true under x11<1E-6
P <- P + geom_line(aes(x=time, y=m_1, color="m", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=m_2, color="m", linetype="sp2"))
P <- P + scale_y_continuous(name="cells")
P <- P + scale_color_discrete(name="")
P <- P + scale_linetype_discrete(name="", labels=c("sp1", "sp2"))
P <- P + theme_light()
plot(P)
ggsave("m=f(t)_2sp.tiff", width=15, height=8, units="cm", path=figure_path, dpi=150)
ggsave("m=f(t)_2sp_bis.tiff", width=15, height=8, units="cm", path=figure_path, dpi=150)

# Dna=f(t)
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=Dna01_1, color="Dna01", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=Dna10_1, color="Dna10", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=Dna11_1, color="Dna11", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=Dna01_2, color="Dna01", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=Dna10_2, color="Dna10", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=Dna11_2, color="Dna11", linetype="sp2"))
P <- P + scale_y_continuous(name="DNA")
P <- P + scale_color_discrete(name="")
P <- P + scale_linetype_discrete(name="", labels=c("sp1", "sp2"))
P <- P + theme_light()
#P <- P + geom_vline(xintercept=2333.17, linetype="dotted", color="red")
plot(P)
ggsave("Dna=f(t)_2sp.tiff", width=15, height=8, units="cm", path=figure_path, dpi=150)

# Dna=f(t) Dna_total
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=Dna01_1+Dna01_2, color="Dna01 total"))
P <- P + geom_line(aes(x=time, y=Dna10_1+Dna10_2, color="Dna10 total"))
P <- P + geom_line(aes(x=time, y=Dna11_1+Dna11_2, color="Dna11 total"))
P <- P + scale_y_continuous(name="DNA")
P <- P + scale_color_discrete(name="")
P <- P + theme_light()
#P <- P + geom_vline(xintercept=2333.17, linetype="dotted", color="red")
plot(P)
ggsave("Dna_total=f(t)_2sp.tiff", width=15, height=8, units="cm", path=figure_path, dpi=150)


#### Plots - 1 env - param variations ####
nbPop <- 1
out_dt <- setDT(as.data.frame(out_parallel))
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
  out_dt[, paste0("pab_cond_",i) := get(paste0("pab_",i)) * Es]
}

# Checks #
m = 0.2166306
s = 0.02503318
muvar(absolve(m, s)[1], absolve(m, s)[2]) # ok
muvar(as.double(absolve(m, s-m^2))[1], as.double(absolve(m, s-m^2))[2]) # nope
a = 3.516429e-06
b = 1
comobg(50,a,b)
out_dt[time==simLength][1:3]
out_dt[time==simLength][param_2==0.1]
View(out_dt[time==simLength])

# Mean last time steps #
out_dt_last_steps <- out_dt[time >= (simLength-4000)] # subset last time steps
out_dt_last_steps[, mean_N_1 := mean(N_1), by=param_set]
out_dt_last_steps[, mean_x00_1 := mean(x00_1), by=param_set]
out_dt_last_steps[, mean_x01_1 := mean(x01_1), by=param_set]
out_dt_last_steps[, mean_x11_1 := mean(x11_1), by=param_set]
out_dt_last_steps[, mean_Dna11_1 := mean(Dna11_1), by=param_set]
out_dt_last_steps[, mean_Dna01_1 := mean(Dna01_1), by=param_set]
out_dt_last_steps[, mean_pab_1 := mean(pab_1), by=param_set]

out_dt_last_steps[, param_1 := log10(param_1)]
out_dt_last_steps[, param_2 := log10(param_2)]

#View(out_dt_last_steps[time==simLength][x11_1>1E-6])

fig_specs <- list(
  path = "./Figures",
  width = 9,
  height = 7,
  dpi = 200
)

color_palette <- c("blue2","yellow1","red1")
color_palette_reverse <- c("red1","yellow1","blue2")
# P <- P + scale_fill_gradient(low="skyblue", high="red1", limits=c(0,0.002))
# P <- P + scale_fill_gradient2(low="blue2", mid="yellow1", midpoint=0.001, high="red1", limits=c(0,0.002)) # ok+
# P <- P + scale_fill_gradientn(colours = c("blue2","yellow1","red1"), limits=c(0,0.002)) # ok+
# P <- P + scale_fill_gradientn(colours = c("blue2","green1","yellow1","orange1","red1"), limits=c(0,0.002)) # ok
# P <- P + scale_fill_viridis_c(option = "plasma", limits=c(0,0.002)) # ok
# P <- P + scale_fill_distiller(palette = "RdPu") # bof
# P <- P + scale_fill_distiller(palette = "YlOrBr") # bof
# P <- P + scale_fill_distiller(palette = "RdYlBu") # ok+

## param_1 & param_2 variations ##
t=simLength

### Heatmaps ###
# x00
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x00_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x00_1/mean_N_1))
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x00_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x11
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_1/mean_N_1))
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x01
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x01_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x01_1/mean_N_1))
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x01_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Arc distance m
P <- ggplot(out_dt_last_steps[time==simLength][x11_1>1E-6])
P <- P + ggtitle(paste0("m t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_1))
P <- P + scale_fill_gradientn(colours=color_palette_reverse, limits=c(0.2,0.25001))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_m_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Pab
P <- ggplot(out_dt_last_steps[time==simLength][x11_1>1E-6])
P <- P + ggtitle(paste0("Pab t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=pab_1))
P <- P + scale_fill_gradientn(colours=color_palette, limits=c(0,0.0015))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_pab_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna11
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna11_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna11_1))
P <- P + scale_fill_gradientn(colours=color_palette)
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna11_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna01
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna01_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna01_1))
P <- P + scale_fill_gradientn(colours=color_palette)
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna01_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)


### geom lines ###
## param_2 vs param_1 ##
# upt_x00_dna11_1
P <- ggplot(out_dt_last_steps[time==simLength])
P <- P + ggtitle(paste0("upt_x00_dna11_1"))
P <- P + geom_line(aes(x=param_2, y=mean_upt_x00_dna11_1, color=factor(param_1)))
P <- P + labs(color="param_1")
P <- P + theme_light()
plot(P)
ggsave(paste0("upt_x00_dna11_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x11
P <- ggplot(out_dt_last_steps[time==simLength])
P <- P + ggtitle(paste0("mean_x11_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_line(aes(x=param_2, y=mean_x11_1/mean_N_1, color=factor(param_1)))
P <- P + labs(color="param_1")
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_mean_x11_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# m
P <- ggplot(out_dt_last_steps[time==simLength & x11_1>1E-6])
P <- P + ggtitle(paste0("m_1"))
P <- P + geom_line(aes(x=param_2, y=m_1, color=factor(param_1)))
P <- P + labs(color="param_1")
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_m_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

## param_1 vs param_2 ##
# x11
P <- ggplot(out_dt_last_steps[time==simLength])
P <- P + ggtitle(paste0("mean_x11_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_line(aes(x=param_1, y=mean_x11_1/mean_N_1, color=factor(param_2)))
#P <- P + geom_line(aes(x=param_1, y=mean_x11_1/mean_N_1, color=factor(param_2, levels=c("1_200","2_200","3_200","3_500","3_1000","4_200","4_500","4_1000"))))
P <- P + labs(color="param_2")
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_mean_x11_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# m
P <- ggplot(out_dt_last_steps[time==simLength & x11_1>1E-6])
P <- P + ggtitle(paste0("m_1"))
P <- P + geom_line(aes(x=param_1, y=m_1, color=factor(param_2)))
#P <- P + geom_line(aes(x=param_1, y=m_1, color=factor(param_2, levels=c("1_200","2_200","3_200","3_500","3_1000","4_200","4_500","4_1000"))))
P <- P + labs(color="param_2")
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_m_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

unique(out_parallel$param_2)


### Dynamics ###
out_dt[treatment==1, atb1 := 0]
out_dt[treatment==1, atb2 := 0]
out_dt[treatment==2, atb1 := 1]
out_dt[treatment==2, atb2 := 1]
out_dt[treatment==3, atb1 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, 2000, trt=3)[1]})]
out_dt[treatment==3, atb2 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, 2000, trt=3)[1]})]
out_dt[treatment==4, atb1 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, 2000, trt=4)[2]})]
out_dt[treatment==4, atb2 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, 2000, trt=4)[3]})]
out_dt[atb1==0, atb1 := NA]
out_dt[atb2==0, atb2 := NA]


## Param_1 & dynamics ##
trt <- 4
#color_lines <- c("blue2","dodgerblue","green2","yellow2","red","blue")
#P <- P + scale_colour_viridis_d()
color_lines <- c("blue2","green3","yellow2","orange2","red","blue")
param_subset <- c(0,1E-3,1E-2,1E-1)
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


## param_1 / treatments ##
## Lines ##
## x11_1 / N_1 ##
P <- ggplot(out_dt_last_steps[time==simLength])
P <- P + ggtitle(paste0("x11_1/N_1"))
P <- P + geom_line(aes(x=param_1, y=mean_x11_1/mean_N_1, color=factor(treatment)))
P <- P + labs(color="environment")
P <- P + scale_color_discrete(labels=c("constant","periodic","alternating"))
P <- P + scale_y_continuous(limits=c(0,1))
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_x11N11_1_trt.tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

## x11_1 ##
P <- ggplot(out_dt_last_steps[time==simLength])
P <- P + ggtitle(paste0("x11_1"))
P <- P + geom_line(aes(x=param_1, y=mean_x11_1, color=factor(treatment)))
P <- P + labs(color="environment")
P <- P + scale_color_discrete(labels=c("constant","periodic","alternating"))
#P <- P + scale_y_continuous(limits=c(0,1))
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_x11_1_trt.tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

## x00_1 ##
P <- ggplot(out_dt_last_steps[time==simLength])
P <- P + ggtitle(paste0("x00_1"))
P <- P + geom_line(aes(x=param_1, y=mean_x00_1, color=factor(treatment)))
P <- P + labs(color="environment")
P <- P + scale_color_discrete(labels=c("constant","periodic","alternating"))
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_x00_1_trt.tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

## x11_1 ##
P <- ggplot(out_dt_last_steps[time==simLength])
P <- P + ggtitle(paste0("x01_1"))
P <- P + geom_line(aes(x=param_1, y=mean_x01_1, color=factor(treatment)))
P <- P + labs(color="environment")
P <- P + scale_color_discrete(labels=c("constant","periodic","alternating"))
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_x01_1_trt.tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

## m_1 ##
P <- ggplot(out_dt_last_steps[time==simLength][mean_x11_1>1E-6])
P <- P + ggtitle(paste0("m_1 "))
P <- P + geom_line(aes(x=param_1, y=m_1, color=factor(treatment)))
P <- P + labs(color="environment")
P <- P + scale_color_discrete(labels=c("constant","periodic","alternating"))
P <- P + scale_y_continuous(limits=c(0,0.25))
P <- P + scale_x_continuous(limits=c(-5,0))
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_m_1_trt.tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)




#### Plots - 1 env & 2 species - param variations ####
fig_specs <- list(
  path = "./Figures",
  width = 9,
  height = 7,
  dpi = 200
)
nbPop <- 2
out_dt <- setDT(as.data.frame(out_parallel))
for(i in 1:nbPop){ # add columns
  out_dt[, paste0("m_",i) := get(paste0("NEf1_",i)) / get(paste0("x11_",i))] # mean arc distance
  out_dt[, paste0("S_",i) := get(paste0("NEf2_",i)) / get(paste0("x11_",i)) - get(paste0("m_",i))^2] # variance arc distance
  out_dt[, paste0("N_",i) := get(paste0("x00_",i)) + get(paste0("x01_",i)) + get(paste0("x10_",i)) + get(paste0("x11_",i))] # size population species 1
}

# Mean last time steps
colnames(out_dt)
out_dt_last_steps <- out_dt[time >= (simLength-4000)] # subset last time steps
out_dt_last_steps[, mean_N_1 := mean(N_1), by=param_set]
out_dt_last_steps[, mean_x00_1 := mean(x00_1), by=param_set]
out_dt_last_steps[, mean_x01_1 := mean(x01_1), by=param_set]
out_dt_last_steps[, mean_x10_1 := mean(x10_1), by=param_set]
out_dt_last_steps[, mean_x11_1 := mean(x11_1), by=param_set]
out_dt_last_steps[, mean_Dna01_1 := mean(Dna01_1), by=param_set]
out_dt_last_steps[, mean_Dna10_1 := mean(Dna10_1), by=param_set]
out_dt_last_steps[, mean_Dna11_1 := mean(Dna11_1), by=param_set]

out_dt_last_steps[, mean_N_2 := mean(N_2), by=param_set]
out_dt_last_steps[, mean_x00_2 := mean(x00_2), by=param_set]
out_dt_last_steps[, mean_x01_2 := mean(x01_2), by=param_set]
out_dt_last_steps[, mean_x10_2 := mean(x10_2), by=param_set]
out_dt_last_steps[, mean_x11_2 := mean(x11_2), by=param_set]
out_dt_last_steps[, mean_Dna01_2 := mean(Dna01_2), by=param_set]
out_dt_last_steps[, mean_Dna10_2 := mean(Dna10_2), by=param_set]
out_dt_last_steps[, mean_Dna11_2 := mean(Dna11_2), by=param_set]
#View(out_dt[time==simLength])

out_dt[param_1==0, param_1:=10^-4.5] # mob -5 -2
out_dt[param_2==0, param_2:=10^-4.5] # dup -4 -1

out_dt[, param_1 := log10(param_1)]
out_dt[, param_2 := log10(param_2)]

## param_1 & param_2 variations ##
t=simLength

# N_1 / (N_1 + N_2)
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("N_1 / (N_1 + N_2)"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_N_1 / (mean_N_1 + mean_N_2)))
P <- P + scale_fill_viridis_c(limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_N_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# N_2 / (N_1 + N_2)
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("N_2 / (N_1 + N_2)"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_N_2 / (mean_N_1 + mean_N_2)))
P <- P + scale_fill_viridis_c(limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_N_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x00_1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x00_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x00_1/mean_N_1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x00_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x00_2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x00_2/mean_N_2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x00_2/mean_N_2))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x00_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x11_1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_1/mean_N_1))
P <- P + scale_fill_viridis_c(limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x11_2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_2/mean_N_2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_2/mean_N_2))
P <- P + scale_fill_viridis_c(limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x01_1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x01_1/mean_N_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x01_1/mean_N_1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x01_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x01_2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x01_2/mean_N_2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x01_2/mean_N_2))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x01_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Arc distance m_1
P <- ggplot(out_dt_last_steps[time==simLength & x11_1>1E-6])
P <- P + ggtitle(paste0("m_1 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_m_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Arc distance m_2
P <- ggplot(out_dt_last_steps[time==simLength & x11_1>1E-6])
P <- P + ggtitle(paste0("m_2 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_2))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_m_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna11_1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna11_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna11_1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna11_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna11_2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna11_2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna11_2))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna11_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna01_1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna01_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna01_1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna01_1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna01_2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna01_1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna01_2))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna01_2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)


### variable=f(param_1) ###
## N_1 / N_tot ##
P <- ggplot(out_dt_last_steps[time==simLength])#[param_2==-4.5]
P <- P + ggtitle(paste0("N_1/(N_1+N_2)"))
P <- P + geom_line(aes(x=param_1, y=mean_N_1/(mean_N_1+mean_N_2), color=factor(treatment)))
P <- P + labs(color="environment")
P <- P + scale_color_discrete(labels=c("constant","periodic","alternating"))
P <- P + scale_y_continuous(limits=c(0,1))
P <- P + geom_hline(yintercept=0.5, linetype="dashed", alpha=0.5)
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_N_1_trt.tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

## x11_1 / N_1 ##
P <- ggplot(out_dt[time==simLength])
P <- P + ggtitle(paste0("x11_1/N_1"))
P <- P + geom_line(aes(x=param_1, y=mean_x11_1/mean_N_1, color=factor(treatment)))
P <- P + labs(color="environment")
P <- P + scale_color_discrete(labels=c("constant","periodic","alternating"))
P <- P + scale_y_continuous(limits=c(0,1))
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_x11_1_trt.tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

## m_1 ##
P <- ggplot(out_dt[time==simLength][mean_x11_1>1E-6])
P <- P + ggtitle(paste0("m_1 "))
P <- P + geom_line(aes(x=param_1, y=m_1, color=factor(treatment)))
P <- P + labs(color="environment")
P <- P + scale_color_discrete(labels=c("constant","periodic","alternating"))
P <- P + scale_y_continuous(limits=c(0,0.25))
P <- P + scale_x_continuous(limits=c(-5,0))
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_m_1_trt.tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)


### Dynamics ###
out_dt[treatment==1, atb1 := 0]
out_dt[treatment==1, atb2 := 0]
out_dt[treatment==2, atb1 := 1]
out_dt[treatment==2, atb2 := 1]
out_dt[treatment==3, atb1 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, 2000, trt=3)[1]})]
out_dt[treatment==3, atb2 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, 2000, trt=3)[1]})]
out_dt[treatment==4, atb1 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, 2000, trt=4)[2]})]
out_dt[treatment==4, atb2 := sapply(time, function(x){deathRateFunc_tris(x, 0, 1, 2000, trt=4)[3]})]
out_dt[atb1==0, atb1 := NA]
out_dt[atb2==0, atb2 := NA]

color_lines <- c("blue2","green3","yellow2","orange2","red","blue")
param_subset <- c(0,1E-3,1E-2,1E-1)
trt <- 2

# m=f(t) #
for(trt in 1:4){
  #P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% c(0,1E-5,1E-4,1E-3,1E-2,1E-1)])
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt])
  P <- P + ggtitle(paste0("m_1=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=m_1, color=factor(param_1)))
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(limits=c(0,0.25))
  P <- P + scale_colour_viridis_d()
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("m=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# N_1/N_tot=f(t) #
for(trt in 1:4){
  P <- ggplot(out_dt[treatment==trt][param_1 %in% param_subset])#[x11_1>1E-6]
  P <- P + ggtitle(paste0("N_1/N_tot=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=N_1/(N_1+N_2), color=factor(param_1)))
  P <- P + geom_line(aes(x=time, y=atb1-1.09, color="atb1"))
  P <- P + geom_line(aes(x=time, y=atb2-1.08, color="atb2"))
  P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
  P <- P + scale_y_continuous(name="N1/Ntot",limits=c(-0.1,1), breaks=seq(0,1,0.2))
  P <- P + scale_color_manual(values=color_lines)
  #P <- P + scale_colour_viridis_d()
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("N1=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# x11_1/N_tot=f(t) #
for(trt in 1:4){
  #P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% c(0,1E-5,1E-4,1E-3,1E-2,1E-1)])
  P <- ggplot(out_dt[x11_1>1E-6][treatment==trt])
  P <- P + ggtitle(paste0("x11_1/N_1=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=x11_1/(N_1), color=factor(param_1)))
  P <- P + scale_y_continuous(limits=c(0,1))
  P <- P + scale_colour_viridis_d()
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("x11=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# x00_1/N_tot=f(t) #
for(trt in 1:4){
  #P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% c(0,1E-5,1E-4,1E-3,1E-2,1E-1)])
  P <- ggplot(out_dt[treatment==trt])
  P <- P + ggtitle(paste0("x00_1/N_1=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=x00_1/(N_1), color=factor(param_1)))
  P <- P + scale_y_continuous(limits=c(0,1))
  P <- P + scale_colour_viridis_d()
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("x00=f(t)_param1_",trt,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# N_1/N_tot=f(t) color bis #
for(trt in 1:4){
  #P <- ggplot(out_dt[x11_1>1E-6][treatment==trt][param_1 %in% c(0,1E-5,1E-4,1E-3,1E-2,1E-1)])
  P <- ggplot(out_dt[treatment==trt])
  P <- P + ggtitle(paste0("N_1/N_tot=f(t) trt",trt))
  P <- P + geom_line(aes(x=time, y=N_1/(N_1+N_2), color=factor(param_1)))
  P <- P + scale_y_continuous(limits=c(0,1))
  P <- P + scale_colour_manual(values=c("black","royalblue2","chartreuse3","gold2","darkorange2","red3"))
  P <- P + labs(color="param_1")
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("N1=f(t)_param1_",trt,"_color2.tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}


### Dynamics 1 param ###

P <- ggplot(out_dt[treatment==trt][param_2 == 0.03])
P <- ggplot(out_dt[treatment==trt][param_1 == 2000])
P <- ggplot(out_dt[treatment==trt][param_2 == 0.03][param_1 == 1000])
P <- P + ggtitle(paste0("N_1/N_tot=f(t) trt",trt))
P <- P + geom_line(aes(x=time, y=N_1/(N_1+N_2), color=factor(param_1)))
P <- P + geom_line(aes(x=time, y=N_1/(N_1+N_2), color=factor(param_2)))
P <- P + scale_x_continuous(name="generation", limits=c(0,simLength))
P <- P + scale_y_continuous(name="N1/Ntot",limits=c(-0.1,1), breaks=seq(0,1,0.2))
P <- P + scale_color_manual(values=color_lines)
P <- P + labs(color="param_1")
P <- P + theme_light()
plot(P)



#### Plots - 1 env - LHS results ####

# Formate results #
figure_path <- "~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures"
# Get result file
out_dt <- setDT(as.data.frame(out_parallel))
# Add columns
nbPop <- 1
setnames(out_dt, "rUpt","rUpt_1")
for(i in 1:nbPop){
  out_dt[, paste0("m_",i) := get(paste0("NEf1_",i)) / get(paste0("x11_",i))] # mean arc distance
  out_dt[, paste0("S_",i) := get(paste0("NEf2_",i)) / get(paste0("x11_",i)) - get(paste0("m_",i))^2] # variance arc distance
  out_dt[, paste0("N_",i) := get(paste0("x00_",i)) + get(paste0("x01_",i)) + get(paste0("x10_",i)) + get(paste0("x11_",i))] # size population species 1
  out_dt[, paste0("alpha_",i) := absolve(get(paste0("m_",i)), get(paste0("S_",i)))[1]] # as.double(absolve(Ef1, Ef2 - Ef1^2))
  out_dt[, paste0("beta_",i) := absolve(get(paste0("m_",i)), get(paste0("S_",i)))[2]]
  out_dt[, paste0("pab_",i) := comobg(gamma, get(paste0("alpha_",i)), get(paste0("beta_",i)))]
  out_dt[, paste0("x11_extinct_",i) := sapply(get(paste0("x11_",i)), function(x) if(x>1E-6){1}else{0})]
  out_dt[, paste0("pab_x11_",i) := get(paste0("x11_",i)) * get(paste0("pab_",i)) * get(paste0("x11_extinct_",i))] # Dna11 release rate
  out_dt[, paste0("pab_x11_x00_",i) := get(paste0("pab_x11_",i)) * get(paste0("x00_",i))]
  out_dt[, paste0("upt_x00_dna11_",i) := get(paste0("Dna11_",i)) * get(paste0("x00_",i)) * get(paste0("rUpt_",i))] # x00 -> x11 rate
}
setDT(lhs)
lhs[, Es := as.double(muvar(1, gamma)[1])]
#lhs[, rMob := rUpt * 10 * rMob_ratio]
#lhs[, rDup := rUpt * 10 * (1 - rMob_ratio)]
#lhs[, rRelease := (rMob + rDup) * Es]
#lhs[, rLoss := rMob * Es]
#lhs[, stressPower := deathPeriod * deathIncrease / lifePeriod]

# Mean last time steps #
out_dt_last_steps <- out_dt[time >= (simLength-2000)] # subset last 2000 time steps
out_dt_last_steps[, mean_N_1 := mean(N_1), by=param_set]
out_dt_last_steps[, mean_x00_1 := mean(x00_1), by=param_set]
out_dt_last_steps[, mean_x01_1 := mean(x01_1), by=param_set]
out_dt_last_steps[, mean_x11_1 := mean(x11_1), by=param_set]
out_dt_last_steps[, mean_pab_1 := mean(pab_1), by=param_set]
out_dt_last_steps[, mean_pab_x11_1 := mean(pab_x11_1), by=param_set]
out_dt_last_steps[, mean_pab_x11_x00_1 := mean(pab_x11_x00_1), by=param_set]
out_dt_last_steps[, mean_upt_x00_dna11_1 := mean(upt_x00_dna11_1), by=param_set]

# Data subsets #
cluster_set <- out_dt[m_1<0.24 & time==simLength][x11_1>1E-6][,param_set]
scattered_set <- out_dt[m_1>0.247 & time==simLength][x11_1>1E-6][,param_set]
weak_cluster_set <- out_dt[m_1<0.247 & m_1>0.24 & time==simLength][x11_1>1E-6][,param_set]
mild_cluster_set <- out_dt[m_1<0.24 & m_1>0.15 & time==simLength][x11_1>1E-6][,param_set]
strong_cluster_set <- out_dt[m_1<0.15 & time==simLength][x11_1>1E-6][,param_set]
extinction_set <- out_dt[x11_1<1E-6 & time==simLength][,param_set]
failed_simu <- seq(1,nrow(lhs),1)[-out_dt[time==simLength][,param_set]]

out_dt[m_1<0.24][time==simLength][x11_1>1E-6]

# Heat map m #
P <- ggplot(out_dt[x11_1>1E-6])
P <- P + ggtitle("")
P <- P + geom_raster(aes(x=time, y=param_set, fill=m_1))
P <- P + scale_fill_viridis_c()#limits=c(0,0.3)
P <- P + labs(color="")
P <- P + theme_light()
plot(P)
ggsave("lhs_1env.tiff", width=10, height=8, units="cm", path=figure_path, dpi=150)

# m dynamics #
P <- ggplot(out_dt[x11_1>1E-6])
P <- ggplot(out_dt[param_set %in% weak_cluster_set])
P <- ggplot(out_dt[param_set %in% mild_cluster_set])
P <- ggplot(out_dt[param_set %in% strong_cluster_set])
P <- P + geom_line(aes(x=time, y=m_1, color=as.character(param_set)))
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_cluster_m.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# x00 dynamics #
P <- ggplot(out_dt)
P <- ggplot(out_dt[param_set %in% weak_cluster_set])
P <- ggplot(out_dt[param_set %in% mild_cluster_set])
P <- ggplot(out_dt[param_set %in% strong_cluster_set])
P <- P + geom_line(aes(x=time, y=x00_1, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_cluster_x00.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# x11 dynamics #
P <- ggplot(out_dt)
P <- ggplot(out_dt[param_set %in% weak_cluster_set])
P <- ggplot(out_dt[param_set %in% mild_cluster_set])
P <- ggplot(out_dt[param_set %in% strong_cluster_set])
P <- P + geom_line(aes(x=time, y=x11_1, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_cluster_x11.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# Dna11 dynamics #
P <- ggplot(out_dt)
P <- ggplot(out_dt[param_set %in% weak_cluster_set])
P <- ggplot(out_dt[param_set %in% mild_cluster_set])
P <- ggplot(out_dt[param_set %in% strong_cluster_set])
P <- P + geom_line(aes(x=time, y=Dna11_1, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_cluster_Dna11.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)


# Parameter distribution - rMob_ratio #
P <- ggplot(as.data.frame(lhs[cluster_set,]))
P <- P + ggtitle(paste0("cluster sets, n=", length(cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[mild_cluster_set,]))
P <- P + ggtitle(paste0("mild cluster sets, n=", length(mild_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[strong_cluster_set,]))
P <- P + ggtitle(paste0("strong cluster sets, n=", length(strong_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[extinction_set,]))
P <- P + ggtitle(paste0("extinction sets, n=", length(failed_simu),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[failed_simu,]))
P <- P + ggtitle(paste0("failed sets, n=", length(failed_simu),"/",nrow(lhs)))

#P <- P + geom_boxplot(aes(x=1, y=(gamma-min(gamma))/(max(gamma)-min(gamma))))
P <- P + geom_boxplot(aes(x=1, y=(log10(rUpt)-log10(min(rUpt)))/(log10(max(rUpt))-log10(min(rUpt)))))
P <- P + geom_boxplot(aes(x=2, y=(rMob_ratio-min(rMob_ratio))/(max(rMob_ratio)-min(rMob_ratio))))
#P <- P + geom_boxplot(aes(x=2, y=(log10(rMob)-log10(min(rMob)))/(log10(max(rMob))-log10(min(rMob)))))
#P <- P + geom_boxplot(aes(x=3, y=(log10(rDup)-log10(min(rDup)))/(log10(max(rDup))-log10(min(rDup)))))
#P <- P + geom_boxplot(aes(x=3, y=(log10(rRelease)-log10(min(rRelease)))/(log10(max(rRelease))-log10(min(rRelease)))))
#P <- P + geom_boxplot(aes(x=4, y=(log10(rLoss)-log10(min(rLoss)))/(log10(max(rLoss))-log10(min(rLoss)))))
P <- P + geom_boxplot(aes(x=3, y=(rClear-min(rClear))/(max(rClear)-min(rClear))))
P <- P + geom_boxplot(aes(x=4, y=(log10(fcost)-log10(min(fcost)))/(log10(max(fcost))-log10(min(fcost)))))
P <- P + geom_boxplot(aes(x=5, y=(deathBasal-min(deathBasal))/(max(deathBasal)-min(deathBasal))))
P <- P + geom_boxplot(aes(x=6, y=(deathIncrease-min(deathIncrease))/(max(deathIncrease)-min(deathIncrease))))
P <- P + geom_boxplot(aes(x=7, y=(deathPeriod-min(deathPeriod))/(max(deathPeriod)-min(deathPeriod))))
P <- P + geom_boxplot(aes(x=8, y=(lifePeriod-min(lifePeriod))/(max(lifePeriod)-min(lifePeriod))))
#P <- P + geom_boxplot(aes(x=13, y=(log10(WTin)-log10(min(WTin)))/(log10(max(WTin))-log10(min(WTin)))))
#P <- P + geom_boxplot(aes(x=14, y=(log10(stressPower)-log10(min(stressPower)))/(log10(max(stressPower))-log10(min(stressPower)))))
#P <- P + geom_boxplot(aes(x=14, y=(stressPower-min(stressPower))/(max(stressPower)-min(stressPower))))
P <- P + scale_y_continuous(name="range ratio")
#P <- P + scale_x_continuous(name="parameters",breaks=1:13,labels=c("gamma","rMob","rDup","rUpt","rRelease","rLoss","rClear","fcost","deathBasal","deathIncrease","deathPeriod","lifePeriod","WTin"))
P <- P + scale_x_continuous(name="",breaks=1:8,labels=c("rUpt","rMob_ratio","rClear","fcost","deathBasal","deathIncrease","deathPeriod","lifePeriod"))
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
P <- P + labs(color="Parameter")
plot(P)
ggsave("lhs_1env_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_strong_cluster_param_bis.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_mild_cluster_param_bis.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_fail_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)


# Parameter distribution - rMob_ratio - 2 species #
P <- ggplot(as.data.frame(lhs[cluster_set,]))
P <- P + ggtitle(paste0("cluster sets, n=", length(cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[mild_cluster_set,]))
P <- P + ggtitle(paste0("mild cluster sets, n=", length(mild_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[strong_cluster_set,]))
P <- P + ggtitle(paste0("strong cluster sets, n=", length(strong_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[extinction_set,]))
P <- P + ggtitle(paste0("extinction sets, n=", length(extinction_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[scattered_set,]))
P <- P + ggtitle(paste0("scattered sets, n=", length(scattered_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[failed_simu,]))
P <- P + ggtitle(paste0("failed sets, n=", length(failed_simu),"/",nrow(lhs)))

P <- P + geom_boxplot(aes(x=1, y=(gamma-min(gamma))/(max(gamma)-min(gamma))))
#P <- P + geom_boxplot(aes(x=2, y=(log10(rMob_ratio)-log10(min(rMob_ratio)))/(log10(max(rMob_ratio))-log10(min(rMob_ratio)))))
P <- P + geom_boxplot(aes(x=2, y=(rMob_ratio-min(rMob_ratio))/(max(rMob_ratio)-min(rMob_ratio))))
P <- P + geom_boxplot(aes(x=3, y=(log10(rMob)-log10(min(rMob)))/(log10(max(rMob))-log10(min(rMob)))))
P <- P + geom_boxplot(aes(x=4, y=(log10(rDup)-log10(min(rDup)))/(log10(max(rDup))-log10(min(rDup)))))
P <- P + geom_boxplot(aes(x=5, y=(log10(rUpt)-log10(min(rUpt)))/(log10(max(rUpt))-log10(min(rUpt)))))
P <- P + geom_boxplot(aes(x=6, y=(rClear-min(rClear))/(max(rClear)-min(rClear))))
P <- P + geom_boxplot(aes(x=7, y=(log10(fcost)-log10(min(fcost)))/(log10(max(fcost))-log10(min(fcost)))))
P <- P + geom_boxplot(aes(x=8, y=(deathBasal-min(deathBasal))/(max(deathBasal)-min(deathBasal))))
P <- P + geom_boxplot(aes(x=9, y=(deathIncrease-min(deathIncrease))/(max(deathIncrease)-min(deathIncrease))))
P <- P + geom_boxplot(aes(x=10, y=(deathPeriod-min(deathPeriod))/(max(deathPeriod)-min(deathPeriod))))
P <- P + geom_boxplot(aes(x=11, y=(lifePeriod-min(lifePeriod))/(max(lifePeriod)-min(lifePeriod))))
#P <- P + geom_boxplot(aes(x=14, y=(log10(stressPower)-log10(min(stressPower)))/(log10(max(stressPower))-log10(min(stressPower)))))
#P <- P + geom_boxplot(aes(x=14, y=(stressPower-min(stressPower))/(max(stressPower)-min(stressPower))))
P <- P + scale_y_continuous(name="range ratio")
P <- P + scale_x_continuous(name="",breaks=1:11,labels=c("gamma","rMob_ratio","rMob","rDup","rUpt","rClear","fcost","deathBasal","deathIncrease","deathPeriod","lifePeriod"))
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
P <- P + labs(color="Parameter")
plot(P)
ggsave("lhs_1env_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_mild_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_strong_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_extinction_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_scattered_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_fail_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_species2_win.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_species2_win_last_steps.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_species1_win_last_steps.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)


# Parameter distribution - classic #
cluster_set <- out_dt[m_1<0.20 & time==simLength][x11_1>1E-6][,param_set]

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
ggsave("lhs_1env_weak_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_mild_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_strong_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_fail_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)

# Variable distribution - classic #
cluster_set <- out_dt[m_1<0.23 & time==simLength][x11_1>1E-6][,param_set]
scattered_set <- out_dt[m_1>0.245 & time==simLength][x11_1>1E-6][,param_set]
intermediary_set <- out_dt[m_1>=0.23 & m_1<=0.245 & time==simLength][x11_1>1E-6][,param_set]
setDT(out_dt_last_steps)
out_dt_last_steps[param_set %in% scattered_set, category := "scattered"]
out_dt_last_steps[param_set %in% intermediary_set, category := "intermediary"]
out_dt_last_steps[param_set %in% cluster_set, category := "clustered"]
out_dt_last_steps[param_set %in% extinction_set, category := "extinction"]
setDF(out_dt_last_steps)
out_dt_last_steps$category <- factor(out_dt_last_steps$category, levels = c("scattered", "intermediary","clustered","extinction"))

P <- ggplot(out_dt_last_steps[out_dt_last_steps$time==simLength,])
P <- P + ggtitle("")
#P <- P + geom_boxplot(aes(x=1, y=mean_N_1, color=category))
P <- P + geom_boxplot(aes(x=1, y=mean_x00_1/mean_N_1, color=category))
P <- P + geom_boxplot(aes(x=2, y=mean_x01_1/mean_N_1, color=category))
P <- P + geom_boxplot(aes(x=3, y=mean_x11_1/mean_N_1, color=category))
#P <- P + geom_boxplot(aes(x=4, y=mean_x00_x11, color=category))
P <- P + scale_y_continuous(name="cell ratio")
P <- P + scale_x_continuous(name="genotypes",breaks=1:3,labels=c("x00","x01","x11"))
P <- P + theme_light()
#P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
P <- P + labs(color="")
plot(P)
ggsave("lhs_1env_var.tiff", width=12, height=8, units="cm", path=figure_path, dpi=150)


# Check hmax #
nrow(out_dt[time==simLength][hmax==1E-1]) # 
nrow(out_dt[time==simLength][hmax==1E-2]) # 
nrow(out_dt[time==simLength][hmax==1E-3]) # 
P <- ggplot(out_dt[time==simLength])
P <- P + geom_bar(aes(x = log10(hmax)))
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_hmax_tf.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)

# DNA distribution #
P <- ggplot(out_dt[time==simLength])
P <- P + geom_histogram(aes(x = log10(Dna10_1), fill="DNA10"), alpha=0.5)
P <- P + geom_histogram(aes(x = log10(Dna01_1), fill="DNA01"), alpha=0.5)
P <- P + geom_histogram(aes(x = log10(Dna11_1), fill="DNA11"), alpha=0.5)
P <- P + scale_x_continuous(name="[eDNA]")
P <- P + labs(fill="eDNA type")
P <- P + theme_light()
plot(P)
ggsave("lhs_DNA_distrib.tiff", width=15, height=10, units="cm", path=figure_path, dpi=150)

# m distribution #
P <- ggplot(out_dt[time==simLength][x11_1>1E-6])
P <- P + geom_histogram(aes(x = m_1, fill="m"), alpha=0.8)
P <- P + scale_x_continuous(name="m")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave("lhs_m_distrib.tiff", width=15, height=10, units="cm", path=figure_path, dpi=150)

# N_1 distribution
P <- ggplot(out_dt_last_steps[time==simLength])
P <- P + geom_histogram(aes(x = mean_N_1, fill="N"), alpha=0.8)
P <- P + scale_x_continuous(name="N")
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave("lhs_N_distrib.tiff", width=15, height=10, units="cm", path=figure_path, dpi=150)



#### Plots - 1 env & 2 species - LHS results ####
# Formate results #
figure_path <- "~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures"
# Get result file
out_dt <- setDT(as.data.frame(out_parallel))
# Add columns
nbPop <- 2
colnames(out_dt)
setnames(out_dt, "rUpt","rUpt_1")
for(i in 1:nbPop){
  out_dt[, paste0("m_",i) := get(paste0("NEf1_",i)) / get(paste0("x11_",i))] # mean arc distance
  out_dt[, paste0("S_",i) := get(paste0("NEf2_",i)) / get(paste0("x11_",i)) - get(paste0("m_",i))^2] # variance arc distance
  out_dt[, paste0("N_",i) := get(paste0("x00_",i)) + get(paste0("x01_",i)) + get(paste0("x10_",i)) + get(paste0("x11_",i))] # size population species 1
  out_dt[, paste0("alpha_",i) := absolve(get(paste0("m_",i)), get(paste0("S_",i)))[1]] # as.double(absolve(Ef1, Ef2 - Ef1^2))
  out_dt[, paste0("beta_",i) := absolve(get(paste0("m_",i)), get(paste0("S_",i)))[2]]
  out_dt[, paste0("pab_",i) := comobg(gamma, get(paste0("alpha_",i)), get(paste0("beta_",i)))]
  out_dt[, paste0("x11_extinct_",i) := sapply(get(paste0("x11_",i)), function(x) if(x>1E-6){1}else{0})]
  out_dt[, paste0("pab_x11_",i) := get(paste0("x11_",i)) * get(paste0("pab_",i)) * get(paste0("x11_extinct_",i))] # Dna11 release rate
  out_dt[, paste0("pab_x11_x00_",i) := get(paste0("pab_x11_",i)) * get(paste0("x00_",i))]
  out_dt[, paste0("upt_x00_dna11_",i) := get(paste0("Dna11_",i)) * get(paste0("x00_",i)) * get(paste0("rUpt_",i))] # x00 -> x11 rate
}

#out_dt[, "rMob_ratio_2" := rMob_2 / (rMob_2 + rDup_2)] # rMob_ratio = rMob / (rMob + rDup)
setDT(lhs)
lhs[, Es := as.double(muvar(1, gamma)[1])]
#lhs[, rMob := rUpt * 10 * rMob_ratio]
#lhs[, rDup := rUpt * 10 * (1-rMob_ratio)]
# lhs[, stressPower := deathPeriod * deathIncrease / lifePeriod] # try to find better indicator

# Mean last 2000 steps #
# with periodic stresses, taking the mean of the last time steps is more representative
out_dt_last_steps <- out_dt[time >= (simLength-5000)] # subset last 2000 time steps
out_dt_last_steps[, mean_N_1 := mean(N_1), by=param_set]
out_dt_last_steps[, mean_x00_1 := mean(x00_1), by=param_set]
out_dt_last_steps[, mean_x01_1 := mean(x01_1), by=param_set]
out_dt_last_steps[, mean_x11_1 := mean(x11_1), by=param_set]
out_dt_last_steps[, mean_Dna11_1 := mean(Dna11_1), by=param_set]
out_dt_last_steps[, mean_Dna01_1 := mean(Dna01_1), by=param_set]
out_dt_last_steps[, mean_pab_1 := mean(pab_1), by=param_set]
out_dt_last_steps[, mean_pab_x11_1 := mean(pab_x11_1), by=param_set]
out_dt_last_steps[, mean_pab_x11_x00_1 := mean(pab_x11_x00_1), by=param_set]
out_dt_last_steps[, mean_upt_x00_dna11_1 := mean(upt_x00_dna11_1), by=param_set]

out_dt_last_steps[, mean_N_2 := mean(N_2), by=param_set]
out_dt_last_steps[, mean_x00_2 := mean(x00_2), by=param_set]
out_dt_last_steps[, mean_x01_2 := mean(x01_2), by=param_set]
out_dt_last_steps[, mean_x11_2 := mean(x11_2), by=param_set]
out_dt_last_steps[, mean_Dna11_2 := mean(Dna11_2), by=param_set]
out_dt_last_steps[, mean_Dna01_2 := mean(Dna01_2), by=param_set]
out_dt_last_steps[, mean_pab_2 := mean(pab_2), by=param_set]
out_dt_last_steps[, mean_pab_x11_2 := mean(pab_x11_2), by=param_set]
out_dt_last_steps[, mean_pab_x11_x00_2 := mean(pab_x11_x00_2), by=param_set]
out_dt_last_steps[, mean_upt_x00_dna11_2 := mean(upt_x00_dna11_2), by=param_set]

# Result groups (1 species) #
mild_cluster_set <- out_dt[m_1<0.24 & m_1>0.1 & time==simLength][x11_1>1][,param_set]
strong_cluster_set <- out_dt[m_1<0.1 & time==simLength][x11_1>1][,param_set]
extinction_set <- out_dt[x11_1<1 & time==simLength][,param_set]
failed_simu <- seq(1,nrow(lhs),1)[-out_dt[time==simLength][,param_set]]
scattered_set <- out_dt[m_1>0.24 & time==simLength][,param_set]
colnames(out_dt)
# Result groups (2 species)
species1_win_50 <- out_dt_last_steps[time==simLength & mean_N_1 / (mean_N_1 + mean_N_2) > 0.5][,param_set]
species1_win_60 <- out_dt_last_steps[time==simLength & mean_N_1 / (mean_N_1 + mean_N_2) > 0.6][,param_set]
species2_win_50 <- out_dt_last_steps[time==simLength & mean_N_2 / (mean_N_1 + mean_N_2) > 0.5][,param_set]
species2_win_60 <- out_dt_last_steps[time==simLength & mean_N_2 / (mean_N_1 + mean_N_2) > 0.6][,param_set]
species2_win_80 <- out_dt_last_steps[time==simLength & mean_N_2 / (mean_N_1 + mean_N_2) > 0.8][,param_set]

# Heat map % N_2 #
P <- ggplot(out_dt)
P <- P + ggtitle("")
P <- P + geom_raster(aes(x=time, y=param_set, fill=N_2/(N_1+N_2)))
P <- P + scale_fill_viridis_c()#limits=c(0,0.3)
P <- P + labs(color="")
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_N2.tiff", width=10, height=8, units="cm", path=figure_path, dpi=150)

# Heat map % N_1 #
P <- ggplot(out_dt)
P <- P + ggtitle("")
P <- P + geom_raster(aes(x=time, y=param_set, fill=N_1/(N_1+N_2)))
P <- P + scale_fill_viridis_c()#limits=c(0,0.3)
P <- P + labs(color="")
P <- P + theme_light()
plot(P)

# Heat map m_1 #
P <- ggplot(out_dt[x11_2>1E-6])
P <- P + ggtitle("")
P <- P + geom_raster(aes(x=time, y=param_set, fill=m_1))
P <- P + scale_fill_viridis_c()#limits=c(0,0.3)
P <- P + labs(color="")
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_m1.tiff", width=10, height=8, units="cm", path=figure_path, dpi=150)

# Heat map m_2 #
P <- ggplot(out_dt[x11_2>1E-6])
P <- P + ggtitle("")
P <- P + geom_raster(aes(x=time, y=param_set, fill=m_2))
P <- P + scale_fill_viridis_c()#limits=c(0,0.3)
P <- P + labs(color="")
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_m2.tiff", width=10, height=8, units="cm", path=figure_path, dpi=150)

# m dynamics #
P <- ggplot(out_dt[x11_2>1E-6])
P <- ggplot(out_dt[param_set %in% cluster_set])
P <- ggplot(out_dt[param_set %in% species2_win_strong][x11_2>1])
P <- P + geom_line(aes(x=time, y=m_2, color=as.character(param_set)))
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_cluster_m_2.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# N_2 dynamics #
P <- ggplot(out_dt)
P <- ggplot(out_dt[param_set %in% species2_win_80])
P <- P + geom_line(aes(x=time, y=N_2, color=as.character(param_set)))
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_2sp_N2.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# x11_2 dynamics #
P <- ggplot(out_dt)
P <- ggplot(out_dt[param_set %in% species2_win_60])
P <- ggplot(out_dt[param_set %in% species2_win_80])
P <- P + geom_line(aes(x=time, y=x11_2, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_2sp_x11_sp2.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# x00_2 dynamics #
P <- ggplot(out_dt)
P <- ggplot(out_dt[param_set %in% species2_win_60])
P <- ggplot(out_dt[param_set %in% species2_win_80])
P <- P + geom_line(aes(x=time, y=x00_2, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_2sp_sp2.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# x11_1 dynamics #
P <- ggplot(out_dt)
P <- ggplot(out_dt[param_set %in% species2_win_60])
P <- ggplot(out_dt[param_set %in% species2_win_80])
P <- P + geom_line(aes(x=time, y=x11_1, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_2sp_x11_sp1.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# x00_1 dynamics #
P <- ggplot(out_dt)
P <- ggplot(out_dt[param_set %in% species2_win_60])
P <- ggplot(out_dt[param_set %in% species2_win_80])
P <- P + geom_line(aes(x=time, y=x00_1, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_2sp_sp1.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# Dna11 dynamics #
P <- ggplot(out_dt)
P <- ggplot(out_dt[param_set %in% cluster_set])
P <- P + geom_line(aes(x=time, y=Dna11_2, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_cluster_Dna11_sp2.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# Parameters distribution 1 #
P <- ggplot(as.data.frame(lhs[cluster_set,]))
P <- P + ggtitle(paste0("cluster sets, n=", length(cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[mild_cluster_set,]))
P <- P + ggtitle(paste0("mild cluster sets, n=", length(mild_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[strong_cluster_set,]))
P <- P + ggtitle(paste0("strong cluster sets, n=", length(strong_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[extinction_set,]))
P <- P + ggtitle(paste0("extinction sets, n=", length(extinction_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[scattered_set,]))
P <- P + ggtitle(paste0("scattered sets, n=", length(scattered_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[failed_simu,]))
P <- P + ggtitle(paste0("failed sets, n=", length(failed_simu),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[species2_win_last_steps,]))
P <- P + ggtitle(paste0("species2 win (>50%, last 2000 steps), n=", length(species2_win_last_steps),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[species2_win_strong,]))
P <- P + ggtitle(paste0("species2 win (>80%, last 2000 steps), n=", length(species2_win_last_steps),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[species1_win_last_steps,]))
P <- P + ggtitle(paste0("species1 win (last 2000 steps), n=", length(species2_win_strong),"/",nrow(lhs)))

P <- P + geom_boxplot(aes(x=1, y=(gamma-min(gamma))/(max(gamma)-min(gamma))))
P <- P + geom_boxplot(aes(x=2, y=(rMob_ratio-min(rMob_ratio))/(max(rMob_ratio)-min(rMob_ratio))))
P <- P + geom_boxplot(aes(x=3, y=(log10(rMob)-log10(min(rMob)))/(log10(max(rMob))-log10(min(rMob)))))
P <- P + geom_boxplot(aes(x=4, y=(log10(rDup)-log10(min(rDup)))/(log10(max(rDup))-log10(min(rDup)))))
P <- P + geom_boxplot(aes(x=5, y=(log10(rUpt)-log10(min(rUpt)))/(log10(max(rUpt))-log10(min(rUpt)))))
P <- P + geom_boxplot(aes(x=6, y=(rClear-min(rClear))/(max(rClear)-min(rClear))))
P <- P + geom_boxplot(aes(x=7, y=(log10(fcost)-log10(min(fcost)))/(log10(max(fcost))-log10(min(fcost)))))
P <- P + geom_boxplot(aes(x=8, y=(deathBasal-min(deathBasal))/(max(deathBasal)-min(deathBasal))))
P <- P + geom_boxplot(aes(x=9, y=(deathIncrease-min(deathIncrease))/(max(deathIncrease)-min(deathIncrease))))
P <- P + geom_boxplot(aes(x=10, y=(deathPeriod-min(deathPeriod))/(max(deathPeriod)-min(deathPeriod))))
P <- P + geom_boxplot(aes(x=11, y=(lifePeriod-min(lifePeriod))/(max(lifePeriod)-min(lifePeriod))))
P <- P + scale_y_continuous(name="range ratio")
P <- P + scale_x_continuous(name="",breaks=1:11,labels=c("gamma","rMob_ratio","rMob","rDup","rUpt","rClear","fcost","deathBasal","deathIncrease","deathPeriod","lifePeriod"))
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
P <- P + labs(color="Parameter")
plot(P)
ggsave("lhs_1env_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_mild_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_strong_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_extinction_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_scattered_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_fail_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_species2_win_last_2000steps.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_species2_win_last_2000steps_strong.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)


# Parameters distribution 2 #
P <- ggplot(as.data.frame(lhs[cluster_set,]))
P <- P + ggtitle(paste0("cluster sets, n=", length(cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[mild_cluster_set,]))
P <- P + ggtitle(paste0("mild cluster sets, n=", length(mild_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[strong_cluster_set,]))
P <- P + ggtitle(paste0("strong cluster sets, n=", length(strong_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[extinction_set,]))
P <- P + ggtitle(paste0("extinction sets, n=", length(extinction_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[species2_win_last_steps,]))
P <- P + ggtitle(paste0("species2 win (>50%, last 2000 steps), n=", length(species2_win_last_steps),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[species2_win_strong,]))
P <- P + ggtitle(paste0("species2 win (>80%, last 2000 steps), n=", length(species2_win_strong),"/",nrow(lhs)))

#P <- P + geom_boxplot(aes(x=1, y=(gamma-min(gamma))/(max(gamma)-min(gamma))))
P <- P + geom_boxplot(aes(x=1, y=(log10(rUpt)-log10(min(rUpt)))/(log10(max(rUpt))-log10(min(rUpt)))))
P <- P + geom_boxplot(aes(x=2, y=(rMob_ratio-min(rMob_ratio))/(max(rMob_ratio)-min(rMob_ratio))))
#P <- P + geom_boxplot(aes(x=3, y=(log10(rMob)-log10(min(rMob)))/(log10(max(rMob))-log10(min(rMob)))))
#P <- P + geom_boxplot(aes(x=4, y=(log10(rDup)-log10(min(rDup)))/(log10(max(rDup))-log10(min(rDup)))))
P <- P + geom_boxplot(aes(x=3, y=(rClear-min(rClear))/(max(rClear)-min(rClear))))
P <- P + geom_boxplot(aes(x=4, y=(log10(fcost)-log10(min(fcost)))/(log10(max(fcost))-log10(min(fcost)))))
P <- P + geom_boxplot(aes(x=5, y=(deathBasal-min(deathBasal))/(max(deathBasal)-min(deathBasal))))
P <- P + geom_boxplot(aes(x=6, y=(deathIncrease-min(deathIncrease))/(max(deathIncrease)-min(deathIncrease))))
P <- P + geom_boxplot(aes(x=7, y=(deathPeriod-min(deathPeriod))/(max(deathPeriod)-min(deathPeriod))))
P <- P + geom_boxplot(aes(x=8, y=(lifePeriod-min(lifePeriod))/(max(lifePeriod)-min(lifePeriod))))
P <- P + scale_y_continuous(name="range ratio")
P <- P + scale_x_continuous(name="",breaks=1:8,labels=c("rUpt","rMob_ratio","rClear","fcost","deathBasal","deathIncrease","deathPeriod","lifePeriod"))
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
P <- P + labs(color="Parameter")
plot(P)
ggsave("lhs_1env_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_mild_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_strong_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_extinction_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_species2_win_last_2000steps.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_species2_win_last_2000steps_strong.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)



# No rMob_ratio #
# Parameters distribution 1 #
P <- ggplot(as.data.frame(lhs[species1_win_60,]))
P <- P + ggtitle(paste0("species1 win (>60%, last 2000 steps), n=", length(species1_win_60),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[species2_win_50,]))
P <- P + ggtitle(paste0("species2 win (>50%, last 2000 steps), n=", length(species2_win_50),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[species2_win_60,]))
P <- P + ggtitle(paste0("species2 win (>60%, last 2000 steps), n=", length(species2_win_60),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[species2_win_80,]))
P <- P + ggtitle(paste0("species2 win (>80%, last 2000 steps), n=", length(species2_win_80),"/",nrow(lhs)))

P <- P + geom_boxplot(aes(x=1, y=(gamma-min(gamma))/(max(gamma)-min(gamma))))
#P <- P + geom_boxplot(aes(x=2, y=(rMob_ratio-min(rMob_ratio))/(max(rMob_ratio)-min(rMob_ratio))))
P <- P + geom_boxplot(aes(x=2, y=(log10(rMob)-log10(min(rMob)))/(log10(max(rMob))-log10(min(rMob)))))
P <- P + geom_boxplot(aes(x=3, y=(log10(rDup)-log10(min(rDup)))/(log10(max(rDup))-log10(min(rDup)))))
P <- P + geom_boxplot(aes(x=4, y=(log10(rUpt)-log10(min(rUpt)))/(log10(max(rUpt))-log10(min(rUpt)))))
P <- P + geom_boxplot(aes(x=5, y=(rClear-min(rClear))/(max(rClear)-min(rClear))))
P <- P + geom_boxplot(aes(x=6, y=(log10(fcost)-log10(min(fcost)))/(log10(max(fcost))-log10(min(fcost)))))
P <- P + geom_boxplot(aes(x=7, y=(deathBasal-min(deathBasal))/(max(deathBasal)-min(deathBasal))))
P <- P + geom_boxplot(aes(x=8, y=(deathIncrease-min(deathIncrease))/(max(deathIncrease)-min(deathIncrease))))
P <- P + geom_boxplot(aes(x=9, y=(deathPeriod-min(deathPeriod))/(max(deathPeriod)-min(deathPeriod))))
P <- P + geom_boxplot(aes(x=10, y=(lifePeriod-min(lifePeriod))/(max(lifePeriod)-min(lifePeriod))))
P <- P + scale_y_continuous(name="range ratio")
P <- P + scale_x_continuous(name="",breaks=1:10,labels=c("gamma","rMob","rDup","rUpt","rClear","fcost","deathBasal","deathIncrease","deathPeriod","lifePeriod"))
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
P <- P + labs(color="Parameter")
plot(P)
ggsave("lhs_1env_species2_win_50.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_species1_win_60.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_species1_win_80.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)



#### Plots - 2 env ####
figure_path <- "~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures"
nbPop <- 1
nEnvironment <- 2
out_dt <- setDT(as.data.frame(out))
for(i in 1:nbPop){ # add columns
  for(j in 1:nEnvironment){
    out_dt[, paste0("m_env",j,"_sp",i) := get(paste0("NEf1_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i))] # mean arc distance
    out_dt[, paste0("S_env",j,"_sp",i) := get(paste0("NEf2_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i)) - get(paste0("m_env",j,"_sp",i))^2] # variance arc distance
    out_dt[, paste0("N_env",j,"_sp",i) := get(paste0("x00_env",j,"_sp",i)) + get(paste0("x01_env",j,"_sp",i)) + get(paste0("x10_env",j,"_sp",i)) + get(paste0("x11_env",j,"_sp",i))] # size population species i
  }
}
fig_specs <- list(
  path = "~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures",
  width = 15,
  height = 12,
  dpi = 150
)

# Cells
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=log10(x00_env1_sp1), color="x00_sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=log10(x01_env1_sp1), color="x01_sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=log10(x10_env1_sp1), color="x10_sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=log10(x11_env1_sp1), color="x11_sp1", linetype="env1"))
#P <- P + geom_line(aes(x=time, y=N_env1_sp1, color="Ntotal", linetype="env1"))
P <- P + geom_line(aes(x=time, y=log10(x00_env2_sp1), color="x00_sp1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=log10(x01_env2_sp1), color="x01_sp1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=log10(x10_env2_sp1), color="x10_sp1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=log10(x11_env2_sp1), color="x11_sp1", linetype="env2"))
#P <- P + geom_line(aes(x=time, y=N_env2_sp1, color="Ntotal", linetype="env2"))
P <- P + scale_y_continuous(name="cells")
#P <- P + scale_y_continuous(name="cells", limits=c(-K/2, K*2))
P <- P + labs(colour="genotypes", linetype="environment")
P <- P + theme_light()
plot(P)
ggsave("Cells_2env.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Ratio genotypes
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=x00_env1_sp1/N_env1_sp1, color="x00_sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=x01_env1_sp1/N_env1_sp1, color="x01_sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=x10_env1_sp1/N_env1_sp1, color="x10_sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=x11_env1_sp1/N_env1_sp1, color="x11_sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=x00_env2_sp1/N_env2_sp1, color="x00_sp1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=x01_env2_sp1/N_env2_sp1, color="x01_sp1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=x10_env2_sp1/N_env2_sp1, color="x10_sp1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=x11_env2_sp1/N_env2_sp1, color="x11_sp1", linetype="env2"))
P <- P + scale_y_continuous(name="ratio")
P <- P + labs(colour="genotypes", linetype="environment")
P <- P + theme_light()
plot(P)
ggsave("Genotypes_ratio_2env.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# DNA
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=Dna01_env1_sp1, color="Dna01_sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=Dna10_env1_sp1, color="Dna10_sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=Dna11_env1_sp1, color="Dna11_sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=Dna01_env2_sp1, color="Dna01_sp1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=Dna10_env2_sp1, color="Dna10_sp1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=Dna11_env2_sp1, color="Dna11_sp1", linetype="env2"))
P <- P + scale_y_continuous(name="DNA")
#P <- P + scale_y_continuous(name="DNA", limits=c(-1, 2))
P <- P + labs(colour="DNA", linetype="environment")
P <- P + theme_light()
plot(P)
ggsave("Dna_2env.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Arc distance
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=m_env1_sp1, color="m", linetype="env1"))
P <- P + geom_line(aes(x=time, y=m_env2_sp1, color="m", linetype="env2"))
#P <- P + geom_line(aes(x=time, y=S_env1_sp1, color="S", linetype="env1"))
#P <- P + geom_line(aes(x=time, y=S_env2_sp1, color="S", linetype="env2"))
P <- P + scale_y_continuous(name="mean arc distance")
#P <- P + scale_y_continuous(name="mean arc distance", limits=c(0,0.5))
P <- P + labs(colour="arc distance", linetype="environment")
P <- P + theme_light()
plot(P)
ggsave("m_2env.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Moments
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=NEf1_env1_sp1, color="NEf1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=NEf1_env2_sp1, color="NEf1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=NEf2_env1_sp1, color="NEf2", linetype="env1"))
P <- P + geom_line(aes(x=time, y=NEf2_env2_sp1, color="NEf2", linetype="env2"))
P <- P + labs(colour="moments", linetype="environment")
P <- P + theme_light()
plot(P)
ggsave("Moments.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

## param_1 variation ##
P <- ggplot(out_dt[time %in% seq(0,5000,100)])
P <- ggplot(out_dt[time %in% seq(0,10000,100)])
P <- P + ggtitle(paste0("env1 m"))
P <- P + geom_raster(aes(x=time, y=param_1, fill=m_env1_sp1))
P <- P + geom_raster(aes(x=time, y=log10(param_1), fill=m_env1_sp1))
P <- P + scale_x_continuous(name="time")
P <- P + scale_y_continuous(name="param 1")
P <- P + scale_fill_viridis_c() #
P <- P + scale_fill_viridis_c(limits=c(0,0.260)) #
P <- P + theme_light()
plot(P)
ggsave("Env1_m.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# m env2
P <- ggplot(out_dt[time %in% seq(0,5000,100)])
P <- ggplot(out_dt[time %in% seq(0,10000,100)])
P <- P + ggtitle(paste0("env2 m"))
P <- P + geom_raster(aes(x=time, y=param_1, fill=m_env2_sp1))
P <- P + geom_raster(aes(x=time, y=log10(param_1), fill=m_env2_sp1))
P <- P + scale_x_continuous(name="time")
P <- P + scale_y_continuous(name="param 1")
P <- P + scale_fill_viridis_c() #
P <- P + scale_fill_viridis_c(limits=c(0,0.260)) #
P <- P + theme_light()
plot(P)
ggsave("Env2_m.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x00 env1
P <- ggplot(out_dt[time %in% seq(0,5000,100)])
P <- ggplot(out_dt[time %in% seq(0,10000,100)])
P <- P + ggtitle(paste0("env1 x00"))
P <- P + geom_raster(aes(x=time, y=param_1, fill=x11_env2_sp1/K1))
P <- P + geom_raster(aes(x=time, y=log10(param_1), fill=x00_env1_sp1/K1))
P <- P + scale_x_continuous(name="time")
P <- P + scale_y_continuous(name="param 1")
P <- P + scale_fill_viridis_c() #
P <- P + scale_fill_viridis_c(limits=c(0,1)) #
P <- P + theme_light()
plot(P)
ggsave("Env1_x00.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x11 env2
P <- ggplot(out_dt[time %in% seq(0,5000,100)])
P <- ggplot(out_dt[time %in% seq(0,10000,100)])
P <- P + ggtitle(paste0("env2 x11"))
P <- P + geom_raster(aes(x=time, y=param_1, fill=x11_env2_sp1/K2))
P <- P + geom_raster(aes(x=time, y=log10(param_1), fill=x11_env2_sp1/K2))
P <- P + scale_x_continuous(name="time")
P <- P + scale_y_continuous(name="param 1")
P <- P + scale_fill_viridis_c() #
P <- P + scale_fill_viridis_c(limits=c(0,1)) #
P <- P + theme_light()
plot(P)
ggsave("Env2_x11.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

## param_1 param_2 variation ##
# Arc distance env1
for(t in seq(0,simLength,1000)){
  P <- ggplot(out_dt[time==t])
  P <- P + ggtitle(paste0("m t=",t))
  P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_env1_sp1))
  #P <- P + scale_x_continuous(name="Treatment period")
  #P <- P + scale_y_continuous(name="Treatment intensity")
  P <- P + scale_fill_viridis_c(limits=c(0,0.260)) #
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("Heatmap_env1_m_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# Arc distance env2
for(t in seq(0,simLength,1000)){
  P <- ggplot(out_dt[time==t])
  P <- P + ggtitle(paste0("m t=",t))
  P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_env2_sp1))
  #P <- P + scale_x_continuous(name="Treatment period")
  #P <- P + scale_y_continuous(name="Treatment intensity")
  P <- P + scale_fill_viridis_c(limits=c(0,0.260)) #
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("Heatmap_env2_m_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}
# x11 env2
for(t in seq(0,simLength,1000)){
  P <- ggplot(out_dt[time==t])
  P <- P + ggtitle(paste0("x11 t=",t))
  P <- P + geom_raster(aes(x=param_1, y=param_2, fill=x11_env2_sp1/(1E7*param_2)))
  #P <- P + scale_x_continuous(name="Treatment period")
  #P <- P + scale_y_continuous(name="Treatment intensity")
  P <- P + scale_fill_viridis_c(limits=c(0,1.1))
  P <- P + theme_light()
  plot(P)
  ggsave(paste0("Heatmap_env2_x11_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)
}


#### Plots - 2 env - param variations ####
fig_specs <- list(
  path = "./Figures",
  width = 10,
  height = 8,
  dpi = 150
)
nbPop <- 1
nbEnvironment <- 2
out_dt <- setDT(as.data.frame(out_parallel))
for(i in 1:nbPop){ # add columns
  for(j in 1:nbEnvironment){
    out_dt[, paste0("m_env",j,"_sp",i) := get(paste0("NEf1_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i))] # mean arc distance
    out_dt[, paste0("S_env",j,"_sp",i) := get(paste0("NEf2_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i)) - get(paste0("m_env",j,"_sp",i))^2] # variance arc distance
    out_dt[, paste0("N_env",j,"_sp",i) := get(paste0("x00_env",j,"_sp",i)) + get(paste0("x01_env",j,"_sp",i)) + get(paste0("x10_env",j,"_sp",i)) + get(paste0("x11_env",j,"_sp",i))] # size population species i
  }
}

# Mean last time steps
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
out_dt_last_steps[, mean_upt_x00_dna11_env1_sp1 := mean_Dna11_env1_sp1 * mean_x00_env1_sp1 * rUpt] # x00 -> x11 rate

out_dt_last_steps[param_1==0, param_1:=10^-4.5]
out_dt_last_steps[param_2==0, param_2:=10^-4.5]
unique(out_dt_last_steps$param_1)

out_dt_last_steps[, param_1 := log10(param_1)]
out_dt_last_steps[, param_2 := log10(param_2)]

## param_1 & param_2 variations ##
t=simLength
# x00
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x00_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x00_env1_sp1/mean_N_env1_sp1))
P <- P + scale_fill_viridis_c() # limits=c(0,1)
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x00_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x11
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_env1_sp1/mean_N_env1_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x01
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x01_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x01_env1_sp1/mean_N_env1_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x01_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Arc distance m
P <- ggplot(out_dt_last_steps[time==simLength & x11_env1_sp1>1E-6])
P <- P + ggtitle(paste0("m_env1_sp1 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_env1_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_m_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna11
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna11_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna11_env1_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna11_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna01
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna01_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna01_env1_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna01_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)


## geom lines ##
# mean_upt_x00_dna11_env1_sp1
P <- ggplot(out_dt[time==simLength])
P <- P + ggtitle(paste0("mean_upt_x00_dna11_env1_sp1"))
P <- P + geom_line(aes(x=param_2, y=mean_upt_x00_dna11_env1_sp1, color=factor(param_1)))
P <- P + labs(color="param_1")
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_mean_upt_x00_dna11_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# x11
P <- ggplot(out_dt[time==simLength])
P <- P + ggtitle(paste0("mean_x11_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_line(aes(x=param_2, y=mean_x11_env1_sp1/mean_N_env1_sp1, color=factor(param_1)))
P <- P + labs(color="param_1")
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_mean_x11_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# m
P <- ggplot(out_dt[time==simLength & x11_env1_sp1>1E-6])
P <- P + ggtitle(paste0("m_env1_sp1"))
P <- P + geom_line(aes(x=param_2, y=m_env1_sp1, color=factor(param_1)))
P <- P + labs(color="param_1")
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_m_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna11
P <- ggplot(out_dt[time==simLength])
P <- P + ggtitle(paste0("mean_Dna11_env1_sp1"))
P <- P + geom_line(aes(x=param_2, y=mean_Dna11_env1_sp1, color=factor(param_1)))
P <- P + labs(color="param_1")
P <- P + theme_light()
plot(P)
ggsave(paste0("Line_Dna11_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)


#### Plots - 2 env - LHS results ####
figure_path <- "~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures"
nbPop <- 1
nbEnvironment <- 2
out_dt <- setDT(as.data.frame(out_parallel))
for(i in 1:nbPop){ # add columns
  for(j in 1:nbEnvironment){
    out_dt[, paste0("m_env",j,"_sp",i) := get(paste0("NEf1_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i))] # mean arc distance
    out_dt[, paste0("S_env",j,"_sp",i) := get(paste0("NEf2_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i)) - get(paste0("m_env",j,"_sp",i))^2] # variance arc distance
    out_dt[, paste0("N_env",j,"_sp",i) := get(paste0("x00_env",j,"_sp",i)) + get(paste0("x01_env",j,"_sp",i)) + get(paste0("x10_env",j,"_sp",i)) + get(paste0("x11_env",j,"_sp",i))] # size population species i
  }
}
setDT(lhs)
lhs[, Es := as.double(muvar(1, gamma)[1])]
lhs[, rRelease := (rMob + rDup) * Es]
lhs[, rLoss := rMob * Es]
lhs[, stressPower := deathPeriod * deathIncrease / lifePeriod]

fig_specs <- list(
  path = "~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures",
  width = 13,
  height = 10,
  dpi = 150
)
mild_cluster_set <- out_dt[m_env1_sp1<0.24 & m_env1_sp1>0.1 & time==simLength][x11_env1_sp1>1][,param_set]
strong_cluster_set <- out_dt[m_env1_sp1<0.1 & time==simLength][x11_env1_sp1>1][,param_set]
extinction_set <- out_dt[x11_env1_sp1<1 & time==simLength][,param_set]
failed_simu <- seq(1,nrow(lhs),1)[-out_dt[time==simLength][,param_set]]
scattered_set <- out_dt[m_env1_sp1>0.24 & time==simLength][,param_set]

out_dt[m_env1_sp1<0.24][time==simLength][x11_env1_sp1>1E-6]

# Heat map m env1 #
P <- ggplot(out_dt)
P <- ggplot(out_dt[x11_env1_sp1>1E-6])
P <- P + ggtitle(paste0("env1 m"))
P <- P + geom_raster(aes(x=time, y=param_set, fill=m_env1_sp1))
#P <- P + geom_raster(aes(x=time, y=param_set, fill=m_env2_sp1))
P <- P + scale_x_continuous(name="time")
P <- P + scale_y_continuous(name="param set")
P <- P + scale_fill_viridis_c() #limits=c(0,0.260)
P <- P + theme_light()
plot(P)
ggsave("lhs_env1_m.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Heat map m env2 #
P <- ggplot(out_dt)
P <- ggplot(out_dt[x11_env1_sp1>0.1])
P <- P + ggtitle(paste0("env2 m"))
P <- P + geom_raster(aes(x=time, y=param_set, fill=m_env2_sp1))
P <- P + scale_x_continuous(name="time")
P <- P + scale_y_continuous(name="param set")
P <- P + scale_fill_viridis_c() #
P <- P + theme_light()
plot(P)
ggsave("lhs_env2_m.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# m dynamics env1 #
P <- ggplot(out_dt[param_set %in% strong_cluster_set])
P <- P + geom_line(aes(x=time, y=m_env1_sp1, color=as.character(param_set)))
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_env1_cluster_m.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# x11 dynamics #
P <- ggplot(out_dt[param_set %in% strong_cluster_set])
P <- P + geom_line(aes(x=time, y=x11_env1_sp1, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_env1_cluster_x11.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# Dna11 dynamics #
P <- ggplot(out_dt[param_set %in% strong_cluster_set])
P <- P + geom_line(aes(x=time, y=Dna11_env1_sp1, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_env1_cluster_Dna11.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# Parameters distribution #
P <- ggplot(as.data.frame(lhs))
P <- P + ggtitle(paste0("all sets, n=", nrow(lhs),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[cluster_set,]))
P <- P + ggtitle(paste0("cluster sets, n=", length(cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[mild_cluster_set,]))
P <- P + ggtitle(paste0("mild cluster sets, n=", length(mild_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[strong_cluster_set,]))
P <- P + ggtitle(paste0("strong cluster sets, n=", length(strong_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[extinction_set,]))
P <- P + ggtitle(paste0("extinction sets, n=", length(extinction_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[scattered_set,]))
P <- P + ggtitle(paste0("scattered sets, n=", length(scattered_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[failed_simu,]))
P <- P + ggtitle(paste0("failed sets, n=", length(failed_simu),"/",nrow(lhs)))

P <- P + geom_boxplot(aes(x=1, y=(gamma-min(gamma))/(max(gamma)-min(gamma))))
P <- P + geom_boxplot(aes(x=2, y=(log10(rMob)-log10(min(rMob)))/(log10(max(rMob))-log10(min(rMob)))))
P <- P + geom_boxplot(aes(x=3, y=(log10(rDup)-log10(min(rDup)))/(log10(max(rDup))-log10(min(rDup)))))
P <- P + geom_boxplot(aes(x=4, y=(log10(rUpt)-log10(min(rUpt)))/(log10(max(rUpt))-log10(min(rUpt)))))
P <- P + geom_boxplot(aes(x=5, y=(log10(rRelease)-log10(min(rRelease)))/(log10(max(rRelease))-log10(min(rRelease)))))
P <- P + geom_boxplot(aes(x=6, y=(log10(rLoss)-log10(min(rLoss)))/(log10(max(rLoss))-log10(min(rLoss)))))
P <- P + geom_boxplot(aes(x=7, y=(rClear-min(rClear))/(max(rClear)-min(rClear))))
P <- P + geom_boxplot(aes(x=8, y=(log10(fcost)-log10(min(fcost)))/(log10(max(fcost))-log10(min(fcost)))))
P <- P + geom_boxplot(aes(x=9, y=(deathBasal-min(deathBasal))/(max(deathBasal)-min(deathBasal))))
P <- P + geom_boxplot(aes(x=10, y=(deathIncrease-min(deathIncrease))/(max(deathIncrease)-min(deathIncrease))))
P <- P + geom_boxplot(aes(x=11, y=(log10(G)-log10(min(G)))/(log10(max(G))-log10(min(G)))))
P <- P + scale_y_continuous(name="range ratio")
P <- P + scale_x_continuous(name="",breaks=1:11,labels=c("gamma","rMob","rDup","rUpt","rRelease","rLoss","rClear","fcost","deathBasal","deathIncrease","G"))
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
P <- P + labs(color="Parameter")
plot(P)
ggsave("lhs_1env_all_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_mild_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_strong_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_extinction_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_scattered_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_fail_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)


# Check hmax #
nrow(out_dt[time==simLength][hmax==1E-1]) # 
nrow(out_dt[time==simLength][hmax==1E-2]) # 
nrow(out_dt[time==simLength][hmax==1E-3]) # 
P <- ggplot(out_dt[time==simLength])
P <- P + geom_bar(aes(x = log10(hmax)))
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_hmax_tf.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)





#### Plots - 2 env - 2species ####
nbPop <- 2
nEnvironment <- 2
out_dt <- setDT(as.data.frame(out))
for(i in 1:nbPop){ # add columns
  for(j in 1:nEnvironment){
    out_dt[, paste0("m_env",j,"_sp",i) := get(paste0("NEf1_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i))] # mean arc distance
    out_dt[, paste0("S_env",j,"_sp",i) := get(paste0("NEf2_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i)) - get(paste0("m_env",j,"_sp",i))^2] # variance arc distance
    out_dt[, paste0("N_env",j,"_sp",i) := get(paste0("x00_env",j,"_sp",i)) + get(paste0("x01_env",j,"_sp",i)) + get(paste0("x10_env",j,"_sp",i)) + get(paste0("x11_env",j,"_sp",i))] # size population species i
  }
}
fig_specs <- list(
  path = "./Figures",
  width = 15,
  height = 12,
  dpi = 150
)

# Cells env1 #
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=log10(x00_env1_sp1), color="x00", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x01_env1_sp1), color="x01", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x10_env1_sp1), color="x10", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x11_env1_sp1), color="x11", linetype="sp1"))
#P <- P + geom_line(aes(x=time, y=log10(x00_env1_sp1+x10_env1_sp1+x01_env1_sp1+x11_env1_sp1), color="total", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x00_env1_sp2), color="x00", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x01_env1_sp2), color="x01", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x10_env1_sp2), color="x10", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x11_env1_sp2), color="x11", linetype="sp2"))
#P <- P + geom_line(aes(x=time, y=log10(x00_env1_sp2+x10_env1_sp2+x01_env1_sp2+x11_env1_sp2), color="total", linetype="sp2"))
P <- P + scale_y_continuous(name="cells")
P <- P + labs(colour="genotypes", linetype="species")
P <- P + theme_light()
plot(P)
ggsave("Cells_env1.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Cells env2 #
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=log10(x00_env2_sp1), color="x00", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x01_env2_sp1), color="x01", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x10_env2_sp1), color="x10", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x11_env2_sp1), color="x11", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x00_env2_sp2), color="x00", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x01_env2_sp2), color="x01", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x10_env2_sp2), color="x10", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x11_env2_sp2), color="x11", linetype="sp2"))
P <- P + scale_y_continuous(name="cells")
P <- P + labs(colour="genotypes", linetype="species")
P <- P + theme_light()
plot(P)
ggsave("Cells_env2.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Arc distance
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=m_env1_sp1, color="sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=m_env2_sp1, color="sp1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=m_env1_sp2, color="sp2", linetype="env1"))
P <- P + geom_line(aes(x=time, y=m_env2_sp2, color="sp2", linetype="env2"))
P <- P + scale_y_continuous(name="mean arc distance")
P <- P + labs(colour="arc distance", linetype="environment")
P <- P + theme_light()
plot(P)
ggsave("m_2env.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)


#### Plots - 2 env - 2species - param variations ####
fig_specs <- list(
  path = "./Figures",
  width = 10,
  height = 8,
  dpi = 150
)
nbPop <- 2
nbEnvironment <- 2
out_dt <- setDT(as.data.frame(out_parallel))
for(i in 1:nbPop){ # add columns
  for(j in 1:nbEnvironment){
    out_dt[, paste0("m_env",j,"_sp",i) := get(paste0("NEf1_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i))] # mean arc distance
    out_dt[, paste0("S_env",j,"_sp",i) := get(paste0("NEf2_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i)) - get(paste0("m_env",j,"_sp",i))^2] # variance arc distance
    out_dt[, paste0("N_env",j,"_sp",i) := get(paste0("x00_env",j,"_sp",i)) + get(paste0("x01_env",j,"_sp",i)) + get(paste0("x10_env",j,"_sp",i)) + get(paste0("x11_env",j,"_sp",i))] # size population species i
  }
}

# Mean last time steps
colnames(out_dt)
out_dt_last_steps <- out_dt[time >= (simLength-4000)] # subset last time steps
# env1 sp1
out_dt_last_steps[, mean_N_env1_sp1 := mean(N_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x00_env1_sp1 := mean(x00_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x01_env1_sp1 := mean(x01_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x10_env1_sp1 := mean(x10_env1_sp1), by=param_set]
out_dt_last_steps[, mean_x11_env1_sp1 := mean(x11_env1_sp1), by=param_set]
out_dt_last_steps[, mean_Dna01_env1_sp1 := mean(Dna01_env1_sp1), by=param_set]
out_dt_last_steps[, mean_Dna10_env1_sp1 := mean(Dna10_env1_sp1), by=param_set]
out_dt_last_steps[, mean_Dna11_env1_sp1 := mean(Dna11_env1_sp1), by=param_set]
# env1 sp2
out_dt_last_steps[, mean_N_env1_sp2 := mean(N_env1_sp2), by=param_set]
out_dt_last_steps[, mean_x00_env1_sp2 := mean(x00_env1_sp2), by=param_set]
out_dt_last_steps[, mean_x01_env1_sp2 := mean(x01_env1_sp2), by=param_set]
out_dt_last_steps[, mean_x10_env1_sp2 := mean(x10_env1_sp2), by=param_set]
out_dt_last_steps[, mean_x11_env1_sp2 := mean(x11_env1_sp2), by=param_set]
out_dt_last_steps[, mean_Dna01_env1_sp2 := mean(Dna01_env1_sp2), by=param_set]
out_dt_last_steps[, mean_Dna10_env1_sp2 := mean(Dna10_env1_sp2), by=param_set]
out_dt_last_steps[, mean_Dna11_env1_sp2 := mean(Dna11_env1_sp2), by=param_set]
# env2 sp1
out_dt_last_steps[, mean_N_env2_sp1 := mean(N_env2_sp1), by=param_set]
out_dt_last_steps[, mean_x00_env2_sp1 := mean(x00_env2_sp1), by=param_set]
out_dt_last_steps[, mean_x01_env2_sp1 := mean(x01_env2_sp1), by=param_set]
out_dt_last_steps[, mean_x10_env2_sp1 := mean(x10_env2_sp1), by=param_set]
out_dt_last_steps[, mean_x11_env2_sp1 := mean(x11_env2_sp1), by=param_set]
out_dt_last_steps[, mean_Dna01_env2_sp1 := mean(Dna01_env2_sp1), by=param_set]
out_dt_last_steps[, mean_Dna10_env2_sp1 := mean(Dna10_env2_sp1), by=param_set]
out_dt_last_steps[, mean_Dna11_env2_sp1 := mean(Dna11_env2_sp1), by=param_set]
# env2 sp2
out_dt_last_steps[, mean_N_env2_sp2 := mean(N_env2_sp2), by=param_set]
out_dt_last_steps[, mean_x00_env2_sp2 := mean(x00_env2_sp2), by=param_set]
out_dt_last_steps[, mean_x01_env2_sp2 := mean(x01_env2_sp2), by=param_set]
out_dt_last_steps[, mean_x10_env2_sp2 := mean(x10_env2_sp2), by=param_set]
out_dt_last_steps[, mean_x11_env2_sp2 := mean(x11_env2_sp2), by=param_set]
out_dt_last_steps[, mean_Dna01_env2_sp2 := mean(Dna01_env2_sp2), by=param_set]
out_dt_last_steps[, mean_Dna10_env2_sp2 := mean(Dna10_env2_sp2), by=param_set]
out_dt_last_steps[, mean_Dna11_env2_sp2 := mean(Dna11_env2_sp2), by=param_set]

out_dt_last_steps[param_1==0, param_1:=10^-4.5]
out_dt_last_steps[param_2==0, param_2:=10^-4.5]

out_dt_last_steps[, param_1 := log10(param_1)]
out_dt_last_steps[, param_2 := log10(param_2)]

## param_1 & param_2 variations ##
t=simLength
# Env1 species1 ratio
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("N_env1_sp1/Ntot_env1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_N_env1_sp1/(mean_N_env1_sp1+mean_N_env1_sp2)))
P <- P + scale_fill_viridis_c(limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_N_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Env2 species1 ratio
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("N_env2_sp1/Ntot_env2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_N_env2_sp1/(mean_N_env2_sp1+mean_N_env2_sp2)))
P <- P + scale_fill_viridis_c(limits=c(0,1))
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_N_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Env1 x11_sp1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_env1_sp1/mean_N_env1_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Env1 x11_sp2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_env1_sp2/mean_N_env1_sp2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_env1_sp2/mean_N_env1_sp2))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_env1_sp2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Env2 x11_sp1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_env2_sp1/mean_N_env2_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_env2_sp1/mean_N_env2_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_env2_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Env2 x11_sp2
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x11_env2_sp2/mean_N_env2_sp2"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x11_env2_sp2/mean_N_env2_sp2))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x11_env2_sp2_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)



# Env1 x01_sp1
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_x01_env1_sp1/mean_N_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_x01_env1_sp1/mean_N_env1_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_x01_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Env1 m_sp1
P <- ggplot(out_dt_last_steps[time==simLength & x11_env1_sp1>1E-6])
P <- P + ggtitle(paste0("m_env1_sp1 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_env1_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_m_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Env1 m_sp2
P <- ggplot(out_dt_last_steps[time==simLength & x11_env1_sp2>1E-6])
P <- P + ggtitle(paste0("m_env1_sp2 t=",t))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=m_env1_sp2))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_m_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)


# Dna11
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna11_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna11_env1_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna11_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Dna01
P <- ggplot(out_dt_last_steps)
P <- P + ggtitle(paste0("mean_Dna01_env1_sp1"))
P <- P + geom_raster(aes(x=param_1, y=param_2, fill=mean_Dna01_env1_sp1))
P <- P + scale_fill_viridis_c()
P <- P + labs(fill="")
P <- P + theme_light()
plot(P)
ggsave(paste0("Heatmap_Dna01_env1_sp1_t=",t,".tiff"), width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)



#### Plots - 2 env - LHS results ####
figure_path <- "~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures"
nbPop <- 1
nbEnvironment <- 2
out_dt <- setDT(as.data.frame(out_parallel))
for(i in 1:nbPop){ # add columns
  for(j in 1:nbEnvironment){
    out_dt[, paste0("m_env",j,"_sp",i) := get(paste0("NEf1_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i))] # mean arc distance
    out_dt[, paste0("S_env",j,"_sp",i) := get(paste0("NEf2_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i)) - get(paste0("m_env",j,"_sp",i))^2] # variance arc distance
    out_dt[, paste0("N_env",j,"_sp",i) := get(paste0("x00_env",j,"_sp",i)) + get(paste0("x01_env",j,"_sp",i)) + get(paste0("x10_env",j,"_sp",i)) + get(paste0("x11_env",j,"_sp",i))] # size population species i
  }
}
setDT(lhs)
lhs[, Es := as.double(muvar(1, gamma)[1])]
lhs[, rRelease := (rMob + rDup) * Es]
lhs[, rLoss := rMob * Es]
lhs[, stressPower := deathPeriod * deathIncrease / lifePeriod]

fig_specs <- list(
  path = "~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures",
  width = 13,
  height = 10,
  dpi = 150
)
mild_cluster_set <- out_dt[m_env1_sp1<0.24 & m_env1_sp1>0.1 & time==simLength][x11_env1_sp1>1][,param_set]
strong_cluster_set <- out_dt[m_env1_sp1<0.1 & time==simLength][x11_env1_sp1>1][,param_set]
extinction_set <- out_dt[x11_env1_sp1<1 & time==simLength][,param_set]
failed_simu <- seq(1,nrow(lhs),1)[-out_dt[time==simLength][,param_set]]
scattered_set <- out_dt[m_env1_sp1>0.24 & time==simLength][,param_set]

out_dt[m_env1_sp1<0.24][time==simLength][x11_env1_sp1>1E-6]

# Heat map m env1 #
P <- ggplot(out_dt)
P <- ggplot(out_dt[x11_env1_sp1>1E-6])
P <- P + ggtitle(paste0("env1 m"))
P <- P + geom_raster(aes(x=time, y=param_set, fill=m_env1_sp1))
#P <- P + geom_raster(aes(x=time, y=param_set, fill=m_env2_sp1))
P <- P + scale_x_continuous(name="time")
P <- P + scale_y_continuous(name="param set")
P <- P + scale_fill_viridis_c() #limits=c(0,0.260)
P <- P + theme_light()
plot(P)
ggsave("lhs_env1_m.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Heat map m env2 #
P <- ggplot(out_dt)
P <- ggplot(out_dt[x11_env1_sp1>0.1])
P <- P + ggtitle(paste0("env2 m"))
P <- P + geom_raster(aes(x=time, y=param_set, fill=m_env2_sp1))
P <- P + scale_x_continuous(name="time")
P <- P + scale_y_continuous(name="param set")
P <- P + scale_fill_viridis_c() #
P <- P + theme_light()
plot(P)
ggsave("lhs_env2_m.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# m dynamics env1 #
P <- ggplot(out_dt[param_set %in% strong_cluster_set])
P <- P + geom_line(aes(x=time, y=m_env1_sp1, color=as.character(param_set)))
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_env1_cluster_m.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# x11 dynamics #
P <- ggplot(out_dt[param_set %in% strong_cluster_set])
P <- P + geom_line(aes(x=time, y=x11_env1_sp1, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_env1_cluster_x11.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# Dna11 dynamics #
P <- ggplot(out_dt[param_set %in% strong_cluster_set])
P <- P + geom_line(aes(x=time, y=Dna11_env1_sp1, color=as.character(param_set)))
P <- P + scale_y_continuous()#limits=c(0,100)
P <- P + labs(color="Param set")
P <- P + theme_light()
plot(P)
ggsave("lhs_env1_cluster_Dna11.tiff", width=15, height=12, units="cm", path=figure_path, dpi=150)

# Parameters distribution #
P <- ggplot(as.data.frame(lhs))
P <- P + ggtitle(paste0("all sets, n=", nrow(lhs),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[cluster_set,]))
P <- P + ggtitle(paste0("cluster sets, n=", length(cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[mild_cluster_set,]))
P <- P + ggtitle(paste0("mild cluster sets, n=", length(mild_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[strong_cluster_set,]))
P <- P + ggtitle(paste0("strong cluster sets, n=", length(strong_cluster_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[extinction_set,]))
P <- P + ggtitle(paste0("extinction sets, n=", length(extinction_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[scattered_set,]))
P <- P + ggtitle(paste0("scattered sets, n=", length(scattered_set),"/",nrow(lhs)))

P <- ggplot(as.data.frame(lhs[failed_simu,]))
P <- P + ggtitle(paste0("failed sets, n=", length(failed_simu),"/",nrow(lhs)))

P <- P + geom_boxplot(aes(x=1, y=(gamma-min(gamma))/(max(gamma)-min(gamma))))
P <- P + geom_boxplot(aes(x=2, y=(log10(rMob)-log10(min(rMob)))/(log10(max(rMob))-log10(min(rMob)))))
P <- P + geom_boxplot(aes(x=3, y=(log10(rDup)-log10(min(rDup)))/(log10(max(rDup))-log10(min(rDup)))))
P <- P + geom_boxplot(aes(x=4, y=(log10(rUpt)-log10(min(rUpt)))/(log10(max(rUpt))-log10(min(rUpt)))))
P <- P + geom_boxplot(aes(x=5, y=(log10(rRelease)-log10(min(rRelease)))/(log10(max(rRelease))-log10(min(rRelease)))))
P <- P + geom_boxplot(aes(x=6, y=(log10(rLoss)-log10(min(rLoss)))/(log10(max(rLoss))-log10(min(rLoss)))))
P <- P + geom_boxplot(aes(x=7, y=(rClear-min(rClear))/(max(rClear)-min(rClear))))
P <- P + geom_boxplot(aes(x=8, y=(log10(fcost)-log10(min(fcost)))/(log10(max(fcost))-log10(min(fcost)))))
P <- P + geom_boxplot(aes(x=9, y=(deathBasal-min(deathBasal))/(max(deathBasal)-min(deathBasal))))
P <- P + geom_boxplot(aes(x=10, y=(deathIncrease-min(deathIncrease))/(max(deathIncrease)-min(deathIncrease))))
P <- P + geom_boxplot(aes(x=11, y=(log10(G)-log10(min(G)))/(log10(max(G))-log10(min(G)))))
P <- P + scale_y_continuous(name="range ratio")
P <- P + scale_x_continuous(name="",breaks=1:11,labels=c("gamma","rMob","rDup","rUpt","rRelease","rLoss","rClear","fcost","deathBasal","deathIncrease","G"))
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
P <- P + labs(color="Parameter")
plot(P)
ggsave("lhs_1env_all_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_mild_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_strong_cluster_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_extinction_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_scattered_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)
ggsave("lhs_1env_fail_param.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)


# Check hmax #
nrow(out_dt[time==simLength][hmax==1E-1]) # 
nrow(out_dt[time==simLength][hmax==1E-2]) # 
nrow(out_dt[time==simLength][hmax==1E-3]) # 
P <- ggplot(out_dt[time==simLength])
P <- P + geom_bar(aes(x = log10(hmax)))
P <- P + theme_light()
plot(P)
ggsave("lhs_1env_hmax_tf.tiff", width=12, height=10, units="cm", path=figure_path, dpi=150)





#### Plots - 2 env - 2species ####
nbPop <- 2
nEnvironment <- 2
out_dt <- setDT(as.data.frame(out))
for(i in 1:nbPop){ # add columns
  for(j in 1:nEnvironment){
    out_dt[, paste0("m_env",j,"_sp",i) := get(paste0("NEf1_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i))] # mean arc distance
    out_dt[, paste0("S_env",j,"_sp",i) := get(paste0("NEf2_env",j,"_sp",i)) / get(paste0("x11_env",j,"_sp",i)) - get(paste0("m_env",j,"_sp",i))^2] # variance arc distance
    out_dt[, paste0("N_env",j,"_sp",i) := get(paste0("x00_env",j,"_sp",i)) + get(paste0("x01_env",j,"_sp",i)) + get(paste0("x10_env",j,"_sp",i)) + get(paste0("x11_env",j,"_sp",i))] # size population species i
  }
}
fig_specs <- list(
  path = "./Figures",
  width = 15,
  height = 12,
  dpi = 150
)

# Cells env1 #
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=log10(x00_env1_sp1), color="x00", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x01_env1_sp1), color="x01", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x10_env1_sp1), color="x10", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x11_env1_sp1), color="x11", linetype="sp1"))
#P <- P + geom_line(aes(x=time, y=log10(x00_env1_sp1+x10_env1_sp1+x01_env1_sp1+x11_env1_sp1), color="total", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x00_env1_sp2), color="x00", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x01_env1_sp2), color="x01", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x10_env1_sp2), color="x10", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x11_env1_sp2), color="x11", linetype="sp2"))
#P <- P + geom_line(aes(x=time, y=log10(x00_env1_sp2+x10_env1_sp2+x01_env1_sp2+x11_env1_sp2), color="total", linetype="sp2"))
P <- P + scale_y_continuous(name="cells")
P <- P + labs(colour="genotypes", linetype="species")
P <- P + theme_light()
plot(P)
ggsave("Cells_env1.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Cells env2 #
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=log10(x00_env2_sp1), color="x00", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x01_env2_sp1), color="x01", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x10_env2_sp1), color="x10", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x11_env2_sp1), color="x11", linetype="sp1"))
P <- P + geom_line(aes(x=time, y=log10(x00_env2_sp2), color="x00", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x01_env2_sp2), color="x01", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x10_env2_sp2), color="x10", linetype="sp2"))
P <- P + geom_line(aes(x=time, y=log10(x11_env2_sp2), color="x11", linetype="sp2"))
P <- P + scale_y_continuous(name="cells")
P <- P + labs(colour="genotypes", linetype="species")
P <- P + theme_light()
plot(P)
ggsave("Cells_env2.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)

# Arc distance
P <- ggplot(out_dt)
P <- P + geom_line(aes(x=time, y=m_env1_sp1, color="sp1", linetype="env1"))
P <- P + geom_line(aes(x=time, y=m_env2_sp1, color="sp1", linetype="env2"))
P <- P + geom_line(aes(x=time, y=m_env1_sp2, color="sp2", linetype="env1"))
P <- P + geom_line(aes(x=time, y=m_env2_sp2, color="sp2", linetype="env2"))
P <- P + scale_y_continuous(name="mean arc distance")
P <- P + labs(colour="arc distance", linetype="environment")
P <- P + theme_light()
plot(P)
ggsave("m_2env.tiff", width=fig_specs$width, height=fig_specs$height, units="cm", path=fig_specs$path, dpi=fig_specs$dpi)



#### Divers ####
# absolve: m & S to alpha & beta
# muvar: alpha & beta to m & S
absolve(0.2546625, 0.020080822) # 0, 65
muvar(0,65) # ! should be mean and variance?
plot(darc(seq(0,0.5,0.01),0,65)) # uniform
absolve(0.2933128, 0.013572511) # 0, 13
plot(darc(seq(0,0.5,0.01),0,13)) # uniform
absolve(0.2480782, 0.02113912) # 7.8E-3, 1.16E2
plot(darc(seq(0,0.5,0.01),0.00782,116)) # cluster

plot(darc(seq(0,0.5,0.01),1,0.5)) # afar genes
plot(darc(seq(0,0.5,0.01),1,1)) # uniform
plot(darc(seq(0,0.5,0.01),1,100)) # close genes

#### Plot LHS data_all ####

### Plots ###
figure_path <- "~/Documents/Article-fabric/Article - Multi-resistance/R_model_JP/Figures"
figure_path <- "C:\\Users\\Gabriel Carvalho\\Documents\\Multi-resistance\\Figures"
fig_specs <- list(
  "width" = 11,
  "heigh" = 9,
  "dpi" = 200
)

## 2 species ##
# N2 ratio #
P <- ggplot(data_all[time==simLength][x11_2>1E-6])
P <- ggplot(data_all[time==simLength][deathPeriod>500][rMob_2<10^(-0.5)])
P <- ggplot(data_all[time==simLength][deathPeriod>500])
P <- ggplot(data_all[time==simLength])
P <- P + ggtitle("Species 2")
P <- P + geom_boxplot(aes(x=as.character(scenario), y=mean_N_2/(mean_N_1+mean_N_2), fill=factor(treatmentType)), lwd=0.4)
P <- P + scale_x_discrete(name="", breaks=1:3, labels=c("scat vs clust","NT vs clust","NT vs scat"))
P <- P + scale_y_continuous(name="N2/(N1+N2)")
P <- P + scale_fill_hue(labels=c("relaxed","constant","periodic","alterning"))
P <- P + geom_hline(yintercept=0.5, linetype="dashed", color="grey", size=0.5)
P <- P + labs(fill="treatment")
P <- P + theme_light()
plot(P)
ggsave("all_1env_N2.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# N2 ratio - violin #
P <- ggplot(data_all[time==simLength])
P <- P + ggtitle("Species 2")
P <- P + geom_violin(aes(x=as.character(scenario), y=mean_N_2/(mean_N_1+mean_N_2), color=factor(treatmentType)))
P <- P + scale_x_discrete(name="", breaks=1:3, labels=c("scat vs clust","NT vs clust","NT vs scat"))
P <- P + scale_y_continuous(name="N2/(N1+N2)")
P <- P + scale_color_hue(labels=c("constant","synchro","altern"))
P <- P + labs(color="treatment")
P <- P + theme_light()
plot(P)
ggsave("all_1env_N2_violin.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# N2 ratio - density #
i <- 3 # treatment type
P <- ggplot(data_all[time==simLength][treatmentType==i])
P <- P + ggtitle(paste0("Species 2 - treatment type = ",i))
P <- P + geom_density(aes(x=mean_N_2/(mean_N_1+mean_N_2), fill=factor(scenario)), alpha=0.5)
P <- P + scale_x_continuous(name="N2/(N1+N2)")
P <- P + scale_fill_hue(labels=c("scat vs clust","NT vs clust","NT vs scat"))
P <- P + labs(fill="scenario")
P <- P + theme_light()
plot(P)
ggsave(paste0("all_1env_N2_density_trt_",i,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# x11 cells sp1 #
P <- ggplot(data_all[time==simLength])
P <- P + ggtitle("Species 1")
P <- P + geom_boxplot(aes(x=as.character(scenario), y=mean_x11_1, color=factor(treatmentType)), lwd=0.4)
P <- P + scale_x_discrete(name="", breaks=1:3, labels=c("scat vs clust","NT vs clust","NT vs scat"))
#P <- P + scale_y_continuous(name="")
P <- P + scale_color_hue(labels=c("constant","synchro","altern"))
P <- P + labs(color="treatment")
P <- P + theme_light()
plot(P)
ggsave("all_1env_x11_sp1.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)

# x11 cells sp2 #
P <- ggplot(data_all[time==simLength])
P <- P + ggtitle("Species 2")
P <- P + geom_boxplot(aes(x=as.character(scenario), y=mean_x11_2, color=factor(treatmentType)), lwd=0.4)
P <- P + scale_x_discrete(name="", breaks=1:3, labels=c("scat vs clust","NT vs clust","NT vs scat"))
#P <- P + scale_y_continuous(name="")
P <- P + scale_color_hue(labels=c("constant","synchro","altern"))
P <- P + labs(color="treatment")
P <- P + theme_light()
plot(P)
ggsave("all_1env_x11_sp2.tiff", width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)


# No rMob_ratio #
# Parameters distribution 1 #
species2_win_60 <- data_all[time==simLength & mean_N_2 / (mean_N_1 + mean_N_2) > 0.6][,param_set]
species2_win_80 <- data_all[time==simLength & mean_N_2 / (mean_N_1 + mean_N_2) > 0.8][,param_set]

scen <- 2

P <- ggplot(data_all[param_set %in% species2_win_60][scenario==scen])
P <- P + ggtitle(paste0("species2 win (>60%, last 5000 steps), scen", scen))

P <- ggplot(data_all[param_set %in% species2_win_80][scenario==scen])
P <- P + ggtitle(paste0("species2 win (>80%, last 5000 steps), scen", scen))

P <- P + geom_boxplot(aes(x=1, y=(log10(deathBasal)-log10(min(deathBasal)))/(log10(max(deathBasal))-log10(min(deathBasal))), fill=factor(treatmentType)))
P <- P + geom_boxplot(aes(x=2, y=(log10(deathIncrease)-log10(min(deathIncrease)))/(log10(max(deathIncrease))-log10(min(deathIncrease))), fill=factor(treatmentType)))
P <- P + geom_boxplot(aes(x=3, y=(log10(deathPeriod)-log10(min(deathPeriod)))/(log10(max(deathPeriod))-log10(min(deathPeriod))), fill=factor(treatmentType)))
P <- P + geom_boxplot(aes(x=4, y=(log10(fcost)-log10(min(fcost)))/(log10(max(fcost))-log10(min(fcost))), fill=factor(treatmentType)))
P <- P + geom_boxplot(aes(x=5, y=(log10(rClear)-log10(min(rClear)))/(log10(max(rClear))-log10(min(rClear))), fill=factor(treatmentType)))
P <- P + geom_boxplot(aes(x=6, y=(log10(rDup_2)-log10(min(rDup_2)))/(log10(max(rDup_2))-log10(min(rDup_2))), fill=factor(treatmentType)))
P <- P + geom_boxplot(aes(x=7, y=(log10(rMob_2)-log10(min(rMob_2)))/(log10(max(rMob_2))-log10(min(rMob_2))), fill=factor(treatmentType)))
P <- P + geom_boxplot(aes(x=8, y=(log10(rUpt_2)-log10(min(rUpt_2)))/(log10(max(rUpt_2))-log10(min(rUpt_2))), fill=factor(treatmentType)))
P <- P + scale_y_continuous(name="range ratio")
P <- P + scale_x_continuous(name="",breaks=1:8,labels=c("deathBasal","deathIncrease","deathPeriod","fcost","rClear","rDup","rMob","rUpt"))
P <- P + scale_fill_hue(labels=c("relaxed","constant","periodic","alterning"))
P <- P + geom_hline(yintercept=0.5, linetype="dashed", color="grey", size=0.5)
P <- P + theme_light()
P <- P + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
P <- P + labs(fill="pressure")
plot(P)
ggsave(paste0("lhs_1env_species2_win_60_scen",scen,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)
ggsave(paste0("lhs_1env_species2_win_80_scen",scen,".tiff"), width=fig_specs$width, height=fig_specs$heigh, units="cm", path=figure_path, dpi=fig_specs$dpi)







## 1 species ##
m_class <- sapply(data_all$m_1, function(x){
  if(x > 0.24){return("m_25_24")}
  if(x < 0.24 & x > 0.20){return("m_24_20")}
  if(x < 0.20 & x > 0.15){return("m_20_15")}
  if(x < 0.15){return("m_15_0")}
})
data_all[, "m_class" := m_class]
rm(m_class)

# m boxplot class #
P <- ggplot(data_all[time==simLength][x11_1>1E-6])
P <- P + ggtitle("")
P <- P + geom_bar(aes(x=m_class, fill=as.character(treatmentType)), position="dodge")
P <- P + scale_x_discrete(name="m class", labels=c("strong","intermediary","low","scattered"))
P <- P + scale_y_continuous(name="arc distance")
P <- P + scale_fill_hue(labels=c("constant","altern","synchro"))
P <- P + labs(fill="treatment")
P <- P + theme_light()
plot(P)
ggsave("all_1sp_m.tiff", width=11, height=8, units="cm", path=figure_path, dpi=150)

# % x00 ##
P <- ggplot(data_all[time==simLength])
P <- P + ggtitle("")
P <- P + geom_boxplot(aes(x=m_class, y=mean_x00_1, color=as.character(treatmentType)), lwd=0.4)
P <- P + scale_x_discrete(name="m class", labels=c("strong","intermediary","low","scattered"))
P <- P + scale_y_continuous(name="WT cells", limits=c(0,1E9))
P <- P + scale_color_hue(labels=c("constant","altern","synchro"))
P <- P + labs(color="treatment")
P <- P + theme_light()
plot(P)
ggsave("all_1sp_x00.tiff", width=11, height=8, units="cm", path=figure_path, dpi=150)

# % x11 ##
P <- ggplot(data_all[time==simLength])
P <- P + ggtitle("")
P <- P + geom_boxplot(aes(x=m_class, y=mean_x11_1, color=as.character(treatmentType)), lwd=0.4)
P <- P + scale_x_discrete(name="m class", labels=c("strong","intermediary","low","scattered"))
P <- P + scale_y_continuous(name="double-gene cells", limits=c(0,1E9))
P <- P + scale_color_hue(labels=c("constant","altern","synchro"))
P <- P + labs(color="treatment")
P <- P + theme_light()
plot(P)
ggsave("all_1sp_x11.tiff", width=11, height=8, units="cm", path=figure_path, dpi=150)

# % x01 ##
P <- ggplot(data_all[time==simLength])
P <- P + ggtitle("")
P <- P + geom_boxplot(aes(x=m_class, y=mean_x01_1, color=as.character(treatmentType)), lwd=0.4)
P <- P + scale_x_discrete(name="m class", labels=c("strong","intermediary","low","scattered"))
P <- P + scale_y_continuous(name="single-gene cells", limits=c(0,1E9))
P <- P + scale_color_hue(labels=c("constant","altern","synchro"))
P <- P + labs(color="treatment")
P <- P + theme_light()
plot(P)
ggsave("all_1sp_x01.tiff", width=11, height=8, units="cm", path=figure_path, dpi=150)













#### Draft ####


### Population step test ###
## Initialization ##
alpha <- 1
beta <- 1
gamma <- 50

rMob = 1
rDup = 1
rUpt = 1
rClear = 1
K = 1E9
b00 = log(2)
b01 = log(2)
b10 = log(2)
b11 = log(2)

x00 = 1E6
x01 = 1E6
x10 = 1E6
x11 = 1E6
Dna01 = 0.5
Dna10 = 0.5
Dna11 = 0.5

Es <- as.double(muvar(1, gamma)[1])

Ef1 <- as.double(muvar(alpha, beta))[1] # 1st moment -> mean
Ef2 <- as.double(muvar(alpha, beta))[2] + Ef1^2 # 2nd moment -> variance

NEf1 <- as.double(x11 * Ef1)
NEf2 <- as.double(x11 * Ef2)



# Current moments
Ef1 <- NEf1 / x11
Ef2 <- NEf2 / x11
# Probability of co-mobilization
pab <- comobg(gamma, alpha, beta)
# Conditional moments: double-gene event
Md <- as.double(muvarg(gamma, alpha, beta))
Ed1 <- Md[1]
Ed2 <- Md[2] + Ed1^2
# Conditional moments: single-gene insertion (uniform arc distance)
Eu1 <- 1/4
Eu2 <- 1/12
# Eta P(A|B)
eta <- pab / Es
eta1 <- 1 - eta
# Conditional moments: single-gene deletion, CORRECTION 21-07-30, P(!A,B)
Es1 <- (Es * Ef1 - pab * Ed1) / (Es * eta1)
Es2 <- (Es * Ef2 - pab * Ed2) / (Es * eta1)
# Introduce population dynamics
N <- x00 + x01 + x10 + x11
rRelease <- Es * rDup # rMob remove but not release DNA
rLoss <- rUpt * N + rClear
# eDNA
dDna01 <- rRelease * (x01 + eta1 * x11 /2) - rLoss * Dna01
dDna10 <- rRelease * (x10 + eta1 * x11 /2) - rLoss * Dna10
dDna11 <- rRelease * (eta * x11)        - rLoss * Dna11
# Gains
tau00_01 <- rUpt * x00 * Dna01
tau00_10 <- rUpt * x00 * Dna10
tau00_11 <- rUpt * x00 * Dna11
tau01_11 <- rUpt * x01 * Dna10
tau10_11 <- rUpt * x10 * Dna01
# Losses
l <- Es * rMob
tau01_00 <- l * x01
tau10_00 <- l * x10
tau11_00 <- l * eta  * x11
tau11_01 <- l * eta1 * x11
tau11_10 <- l * eta1 * x11

# Variable death rates
dR <- deathRateFunc_tris(Time,db=0.2,di=0.,dp=200,trt=2,lp=200) #

damp <- 1 - N / K # logistic growth rate
r00 <- (b00 * damp - dR[1])
r01 <- (b01 * damp - dR[2])
r10 <- (b10 * damp - dR[3])
r11 <- (b11 * damp - dR[4])

# Populations
dx00 <- r00 * x00 + (tau01_00 + tau10_00 + tau11_00) - (tau00_01 + tau00_10 + tau00_11)
dx01 <- r01 * x01 + (tau00_01 + tau11_01) - (tau01_11 + tau01_00)
dx10 <- r10 * x10 + (tau00_10 + tau11_10) - (tau10_11 + tau10_00)
dx11 <- r11 * x11 + (tau00_11 + tau01_11 + tau10_11) - (tau11_00 + tau11_01 + tau11_10)

# Component derivatives (intermediate forms have underscore suffix)
dNf_ <- r11 * x11
dNd_ <- tau00_11 - tau11_00
dNs_ <- -(tau11_01 + tau11_10)
dNu_ <- +(tau01_11 + tau10_11)

dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1
dNEf2 <- dNf_ * Ef2 + dNd_ * Ed2 + dNs_ * Es2 + dNu_ * Eu2


# Tests #
dNs_ <- tau11_01 + tau11_10
dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1

dNEf1 <- dNf_ * Ef1 + tau00_11 * Ed1 + dNu_ * Eu1


dNf_ <- r11 * x11 - (tau11_00 + tau11_01 + tau11_10)
dNd_ <- tau00_11
dNs_ <- 0
dNu_ <- tau01_11 + tau10_11
dNEf1 <- dNf_ * Ef1 + dNd_ * Ed1 + dNs_ * Es1 + dNu_ * Eu1



dNf_ * Ef1
dNd_ * Ed1
dNs_ * Es1
dNu_ * Eu1

(NEf1 + dNEf1) / (x11 + dx11) # Ef1
(NEf2 + dNEf2) / (x11 + dx11) # Ef2
