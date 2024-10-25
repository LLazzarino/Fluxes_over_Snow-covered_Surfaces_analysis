#
library(wavelets)
library(plyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(RColorBrewer)
library(raster)
library(rasterVis)
library(reshape2)
library("pracma")
library(gtools)
library(chron)
library(grid)
library(zoo)
library(data.table) 
##########################################################################
# funzioni
##########################################################################
source("/disk1/luca/FLOSS/Completi_SUB_FLOSS.R")
source("/disk1/luca/FLOSS/Spettri_FLOSS_SUB_60Hz.R")
source("/disk1/luca/FLOSS/Sub_FLOSS.R")
#source("/Users/luca/Documents/DatiAmazzonia/Linhares/Sub_MRA_Linhares.R")
#source("/Users/luca/Documents/DatiAmazzonia/Linhares/Sub_read_Linhares.R")
##########################################################################
# parametri
##########################################################################
options(digits.secs=1)
#
ndati = 216000
numero = 1:ndati
deltat = 1/60
lagtouse = ndati/10
lagtofit = lagtouse/2
Nit<-30                    # parametri
controllo=1/(2^20)         # per i fit 
nums <- as.numeric(1:10 %o% 10 ^ (-11:2))
# costanti per log_block_average (valide per spettri di 72000 dati)
dln =(log(10)-log(1))/20
edln = exp(dln)
fstart = .01
nint = 60
# costanti per il fit degli spettri
twothird=2/3
fourthird=4/3
fivethird=5/3
alfa1 = 0.5
numerofit=20
#
rip = 2048 # lunghezza di ogni spettro a bassa frequenza
ripuniti = 2090
intervallotot   = 1 : ripuniti
intervallostart = 1 : 41 
intervallomedia = 42 : 52
intervalloalta  = 1 : rip
intervallo      = 1 : rip
nn = 4096*52 # lunghezza massima dei dati per lo spettro...
numtmp = 0
#
plotta_auto_spettri = TRUE
##########################################################################
# setto la directory di lavoro
##########################################################################
setwd("/disk1/luca/FLOSS/")
#ls()
##########################################################################
# indici di lettura
##########################################################################
zs = c(1,2,5,10,15,20,30) # livelli FLOSS
uidx = seq(1,28,4) # indici u
vidx = uidx + 1    # indici v
widx = uidx + 2    # indici w
Tidx = uidx + 3    # indici T
##########################################################################
# lettura FLOSS
#
# ciclo_file_giornalieri (ci sono solo le ore notturne)
#########################################################################
file_orari = list.files("Dati")
numero_orari = length(file_orari)
mesi = c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
ore = c(20,21,22,23,0,1,2,3,4,5)
#
for(j in 1:numero_orari){ # ciclo sui giorni
  print(paste("inizio file",j,"di",numero_orari))
  #j = 66
  ggiorno = as.numeric(substr(file_orari[j],4,5))
  gmese   = substr(file_orari[j],1,3)
  gmese   = which(gmese == mesi)
  file_read = readLines(paste("Dati/",file_orari[j],sep=""))
  temp = tempfile()
#  writeLines(gsub(pattern = "([*])\\1+", replace = " NA",x = file_read, perl = TRUE), 
#             temp, sep="\n") ## Write it to a temporary file
  writeLines(gsub("\\*{7}", replace = " NA",x = file_read, perl = TRUE), 
             temp, sep="\n") ## Write it to a temporary file
  dati_giorno = fread(temp,header=FALSE,nrow = 2160000)
  dati_giorno = as.data.frame(dati_giorno)
  colnames(dati_giorno)[uidx] = paste("u",zs,sep="_")
  colnames(dati_giorno)[vidx] = paste("v",zs,sep="_")
  colnames(dati_giorno)[widx] = paste("w",zs,sep="_")
  colnames(dati_giorno)[Tidx] = paste("T",zs,sep="_")
  dati_giorno[,Tidx] = dati_giorno[,Tidx] + 273.16 # porto la temperatura in Kelvin
  rm(temp)
#  
  for (i in 1:10)  { #ciclo sulle ore
  # i = 9
  gora    = i-1
  gminuti = 30
  if (mesi[gmese] %in% c("nov","dec")) ganno = "2002"
  if (mesi[gmese] %in% c("jan","feb","mar","apr")) ganno = "2003"
#  data = paste(ganno,"-",gmese,"-0",ggiorno," ",gora,":",gminuti,sep="")
  data = paste(ganno,"-",gmese,"-",ggiorno," ",gora,":",gminuti,sep="")
  data = strptime(data,"%Y-%m-%d %H:%M",tz="UTC")-(3600*4)
  indici_ora = ((i-1)*ndati+1):(ndati*i)
  ####################
  # ciclo_livelli
  ####################
  data_char = as.character(format(data,"%Y-%m-%d %H:%M:%S",tz="UTC"))
  # 
  for(kk in 1:7){
    #kk = 4
    dati_ora = dati_giorno[indici_ora,c(uidx[kk],vidx[kk],widx[kk],Tidx[kk])]
#    str(dati_ora)
    colnames(dati_ora) = c("u","v","w","T")
    #
    if(((all(is.na(dati_ora$u)) |  all(is.na(dati_ora$v)) |
         all(is.na(dati_ora$w)) |  all(is.na(dati_ora$T))))) {
      tutto_lev = data.frame(GooN=FALSE)
      scrivi_lev = t(c(data_char,rep(NA,51)))
#      scriviacf_lev = data.frame(rep(data_char,lagtouse+1),matrix(NA,nrow = lagtouse+1,ncol = 5))
#      scrivispettri_lev = data.frame(rep(data_char,rip),matrix(NA,nrow = rip,ncol = 10))
    } else {
      tutto_lev = calcoli_FLOSS(dati_ora,zs[kk],FALSE)
      print(kk)
      print(summary(dati_ora))
      #      print("qui")
      tutto_lev_1min = calcoli_5min(dati_ora,zs[kk],1,FALSE)
      scrivi_lev = data.frame(data_char,tutto_lev$NA_percent,tutto_lev$angolo,tutto_lev$mean_u,tutto_lev$mean_u_H,
                              tutto_lev$mean_T,
                              tutto_lev$sigma_u,tutto_lev$sigma_v,tutto_lev$sigma_w,tutto_lev$sigma_T, 
                              tutto_lev$uw,tutto_lev$vw,tutto_lev$uT,tutto_lev$vT,tutto_lev$wT,
                              tutto_lev$ustar_u,tutto_lev$UnoSuL,tutto_lev$ZsuL,tutto_lev$thetastar,
                              tutto_lev$p_u,tutto_lev$q_u,tutto_lev$p_v,tutto_lev$q_v, 
                              tutto_lev$p_w,tutto_lev$q_w,tutto_lev$p_T,tutto_lev$q_T, 
                              tutto_lev$Rminu,tutto_lev$Rminv,tutto_lev$Rminw,tutto_lev$RminT,
                              tutto_lev$epsilon_u,tutto_lev$epsilon_v,tutto_lev$epsilon_w,
                              tutto_lev$Kol_u,tutto_lev$Kol_v,tutto_lev$Kol_w,tutto_lev_1min,
                              tutto_lev$thetamed,tutto_lev$sigma_theta)
      #
      # scriviacf_lev = data.frame(data_char,tutto_lev$acfu$lag,tutto_lev$acfu$acf,tutto_lev$acfv$acf,
      #                            tutto_lev$acfw$acf,tutto_lev$acfT$acf)
      # #
      # scrivispettri_lev = data.frame(data_char,tutto_lev$spettro_u$low_freq,tutto_lev$spettro_u$low_spec,
      #                                tutto_lev$spettro_v$low_spec,
      #                                tutto_lev$spettro_w$low_spec,tutto_lev$spettro_T$low_spec,
      #                                tutto_lev$spettro_u$high_freq,tutto_lev$spettro_u$high_spec,
      #                                tutto_lev$spettro_v$high_spec,tutto_lev$spettro_w$high_spec,
      #                                tutto_lev$spettro_T$high_spec)
      if(plotta_auto_spettri){
        ######## autocorrelazioni
        plotta_acf(data_char,tutto_lev$acfu,tutto_lev$acfv,tutto_lev$acfw,tutto_lev$acfT,
                   tutto_lev$mean_u,tutto_lev$mean_T,tutto_lev$sigma_u,
                   tutto_lev$sigma_v,tutto_lev$sigma_w,tutto_lev$sigma_T,
                   tutto_lev$q_u,tutto_lev$p_u,tutto_lev$q_v,tutto_lev$p_v,tutto_lev$q_w,
                   tutto_lev$p_w,tutto_lev$q_T,tutto_lev$p_T,
                   zs[kk],paste(zs[kk],"m"))
        ######## spettri
        if(tutto_lev$GooN) {
          plotta_spettri(data_char,tutto_lev$spettro_u,tutto_lev$spettro_v,tutto_lev$spettro_w,tutto_lev$spettro_T,
                         tutto_lev$epsilon_u,tutto_lev$epsilon_v,tutto_lev$epsilon_w,
                         tutto_lev$p_u,tutto_lev$q_u,tutto_lev$p_v,tutto_lev$q_v,tutto_lev$p_w,tutto_lev$q_w,tutto_lev$p_T,tutto_lev$q_T,
                         tutto_lev$ZsuL,tutto_lev$mean_u,
                         zs[kk],paste(zs[kk],"m"))
        }
      } #fine plotta
    }  
    ########### scrivo file
     filewrite = paste("Output/Riassunto_FLOSS_",zs[kk],"m.txt",sep="")
     if(j == 1 & i == 1){
       write.table(file=filewrite,scrivi_lev,
                   row.names=FALSE,col.names=c("data","NA_percent","angolo","mean_u","mean_u_H","mean_T",
                                               "sigma_u","sigma_v","sigma_w","sigma_T",
                                               "uw","vw","uT","vT","wT",
                                               "ustar_u","UnoSuL","ZsuL","thetastar",
                                               "p_u","q_u","p_v","q_v",
                                               "p_w","q_w","p_T","q_T",
                                               "Rminu","Rminv","Rminw","RminT",
                                               "epsilonu","epsilonv","epsilonw",
                                               "Kol_u","Kol_v","Kol_w",
                                               "sigma_u_1min","sigma_v_1min","sigma_w_1min","sigma_T_1min",
                                               "uw_1min","vw_1min","uT_1min","vT_1min",
                                               "wT_1min","ustar_1min","H0_1min","ZsuL_1min",
                                               "thetastar_1min","thetamed","sigma_theta"),
                   quote=FALSE,sep="\t",append=FALSE)
     }else{
       write.table(file=filewrite,scrivi_lev,row.names=FALSE,
                   col.names=FALSE,quote=FALSE,sep="\t",append=TRUE)
     } 
  } #end do kk sui 7 livelli
  ##
  print(paste("fine ora",i))
  } # fine ciclo ore
  print(data_char)
  print(paste("fine file",j,"di",numero_orari))
} #fine ciclo j






