setwd("/Users/haeun/OneDrive - University of North Carolina at Chapel Hill/Rotation3_Palmer/Work/data/NSCLC/")
fileprefix = "NSNSCLC_atezo_chemo_Socinski2018"

TSD <- read.csv(paste0(fileprefix, ".csv"))
digizeit <- DIGI.CLEANUP(TSD)

# if survival is in percentage, convert it to 0-1
if (max(digizeit$S) > 2){
  digizeit$S = digizeit$S / 100
}

#DF <- read.csv(file.path("at_risk_table.csv"))
#AR <-filter(DF, Slide == 'C', Arm == 'TC3/IC3')
AR <- read.csv(paste0(fileprefix, '_atrisk.csv'))
pub.risk <- K.COORDINATES(AR, digizeit)

IPD <- GENERATEINDIVIDUALDATA(tot.events = unique(AR$TE), arm.id = unique(AR$Arm), digizeit = digizeit, pub.risk = pub.risk)
IPD <- data.frame(Time=IPD$Time, Event=IPD$Event, Arm=IPD$Arm, Subpop=rep(unique(AR$Subpop), nrow(IPD)), Slide = rep(unique(AR$Slide), nrow(IPD)))
nm <- paste(unique(AR$Slide),  paste(unique(AR$Arm), unique(AR$Subpop), sep='_'),sep='_')
write.csv(IPD, paste0(fileprefix, "_indiv.csv"), row.names = FALSE)

#IPD$Time <- as.numeric(as.character(IPD$Time))
#IPD$Event <- as.numeric(as.character(IPD$Event))
#survfit(Surv(Time, Event) ~ 1, data=IPD)