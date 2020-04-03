TSD <- read.csv("D:/Dropbox/ExtractIPD/exampledata/input/F_OS_CPS_1_pembrolizumab.csv")
digizeit <- DIGI.CLEANUP(TSD)

DF <- read.csv("D:/Dropbox/ExtractIPD/exampledata/input/AR.csv")
AR <-filter(DF, Slide == 'A', Arm == 'pembrolizumab')
pub.risk <- K.COORDINATES(AR, digizeit)



IPD <- GENERATEINDIVIDUALDATA(tot.events = unique(AR$TE), arm.id = unique(AR$Arm), digizeit = digizeit, pub.risk = pub.risk)
IPD <- data.frame(Time=IPD$Time, Event=IPD$Event, Arm=IPD$Arm, Subpop=rep(unique(AR$Subpop), nrow(IPD)), Slide = rep(unique(AR$Slide), nrow(IPD)))
nm <- paste(unique(AR$Slide),  paste(unique(AR$Arm), unique(AR$Subpop), sep='_'),sep='_')
write.csv(IPD, paste(paste('D:/Dropbox/ExtractIPD/exampledata/output', nm, sep='/'), 'csv', sep='.'), row.names = FALSE)

IPD$Time <- as.numeric(as.character(IPD$Time))
IPD$Event <- as.numeric(as.character(IPD$Event))
#survfit(Surv(Time, Event) ~ 1, data=IPD)