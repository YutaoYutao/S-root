# For: titling assay GWAS screening
# Date: 13-10-2017
# CalcDirection: Calculates direction based on change between next and previous//based on lastNode direction

#RITEtoData: Transform Dataset
path = "/Users/zyt/Desktop/test-pst2/" #path where csvfiles are
setupFile <- read.csv("pst2plate.csv", sep = ",") #make a file explain the 8 plates.
timeStep <- 20 #minutes
ExpName <- "Direction"
SorW <- "W" #Summertime=S, wintertime=W
View(setupFile)

RITEtoData <- function(path, setupFile, timeStep, ExpName, SorW){
  storewd <- getwd()
  setwd(path)
  out.file<- matrix(, nrow = 0, ncol = 0)
  file.names <- dir(path, pattern =".csv")
  for(i in 1:length(file.names)){
    file <- read.csv(file.names[i],header=TRUE, sep=",", stringsAsFactors=FALSE)
    filename <- t(as.data.frame(strsplit(as.character(file.names[i]),"\\_|\\."))) #split imagename
    colnames(filename) <- c("sr", "Plate", "image") #Add column names
    file$Plate <- filename[,2]


    DP <- t(as.data.frame(strsplit(as.character(file$image),SorW))) #split imagename
    colnames(DP) <- c("Date", "Image") #Add column names
    file <- cbind(file, DP[,1:2]) #bind date to the current dataset
    file$Date <- strptime(file$Date, "%Y%m%d_%H%M%S")
    file$Date <- as.POSIXct(file$Date)

    out.file <- rbind(out.file, file)
  }

  Raw_main <- matrix(ncol=ncol(out.file), nrow = 0)
  colnames(Raw_main) <- colnames(out.file) ##make matrix for main roots

  LatRoot <- FALSE
  if(nlevels(as.factor(out.file$root_order)) > 1){
  Raw_lat <- matrix(ncol=ncol(out.file), nrow = 0)
  colnames(Raw_lat) <- colnames(out.file) #Make matrix for lateral roots
  LatRoot <- TRUE
  }

  #loop over lines SRoutputFile to seperate main and lateral roots
  if(LatRoot == TRUE){
  for(i in 1:nrow(out.file)){
    if(out.file$root_order[i] == 0){
      Raw_main <- rbind(Raw_main, out.file[i,])
    }
    if(out.file$root_order[i] > 0){
      Raw_lat <- rbind(Raw_lat, out.file[i,])
    }
  }
  }
  if(LatRoot == FALSE) Raw_main <- out.file

  NR <- t(as.data.frame(strsplit(as.character(Raw_main$root_name), "\\_")))
  colnames(NR) <-  c("ontology", "Root")
  rownames(NR) <- NULL
  Raw_main <- cbind(Raw_main, NR[,1:2])

  Raw_main$Treatment = "X"
  Raw_main$Genotype = "X"
  Raw_main$Root <- as.numeric(levels(Raw_main$Root)[Raw_main$Root])
  Raw_main$Plate <- as.integer(Raw_main$Plate)

  #Extract info
  for(a in 1:nrow(Raw_main)){
    for (b in 1:nrow(setupFile)){
      if (Raw_main$Plate[a] == setupFile$Plate[b]){
      if (Raw_main$Root[a] == setupFile$Root[b]){
        Raw_main$Treatment[a] <- as.character(setupFile$Treatment[b])
        Raw_main$Genotype[a] <- as.character(setupFile$genotypes[b])
      }
    }
    }
  }

  Raw_main$root <- as.factor(Raw_main$root)
  Raw_main$Time <- NA

  #Sort on time/root
  require(dplyr)
  require(plyr)
  Raw_main <- arrange(Raw_main, Plate, Root, root, Date)

  time =0
  StartDate = min(Raw_main$Date)
  Day = 24*60
    for (i in 1:nrow(Raw_main)){
      if(Raw_main$Date[i] < (StartDate+60*60)) Raw_main$Time[i] = as.numeric(Raw_main$Date[i] - StartDate)
      else if(Raw_main$Date[i] < (StartDate+Day*60)) Raw_main$Time[i] = as.numeric(Raw_main$Date[i] - StartDate)*60
      else Raw_main$Time[i] = as.numeric(Raw_main$Date[i] - StartDate)*Day
      }

  setwd(storewd)
  write.csv(Raw_main, file= paste(ExpName, "_Raw_main.csv", sep=""))


  if(LatRoot == TRUE){
  NR2 <- t(as.data.frame(strsplit(as.character(Raw_lat$parent_name), "\\_")))
  colnames(NR2) <-  c("ontology", "Parent")
  rownames(NR2) <- NULL
  Raw_lat <- cbind(Raw_lat, NR2[,1:2])

  Raw_lat$Treatment = "X"
  Raw_lat$Genotype = "X"
  Raw_lat$Plate <- as.integer(Raw_lat$Plate)

  for(a in 1:nrow(Raw_lat)){
    for (b in 1:nrow(setupFile)){
      if (Raw_lat$Plate[a] == setupFile$Plate[b]){
        if (Raw_lat$Parent[a] == setupFile$Root[b]){
          Raw_lat$Treatment[a] <- as.character(setupFile$Treatment[b])
          Raw_lat$Genotype[a] <- as.character(setupFile$genotypes[b])
        }
      }
    }
  }


  Raw_lat <- arrange(Raw_lat, Plate, Parent, root_name, root, Date)
  Raw_lat$root <- as.factor(Raw_lat$root)
  Raw_lat$Time <- NA
  time =0
  minDay = min(Raw_main$Date)
  Day = 60*24
  Raw_lat$EmergeDay = NA
  for (i in 1:nrow(Raw_lat)){
    Raw_lat$Time[i] = time

    if (Raw_lat$Time[i] == 0){
      if(Raw_lat$Date[i] < (minDay+Day*60)) Raw_lat$EmergeDay[i] = 1
      else Raw_lat$EmergeDay[i] = ceiling(as.numeric(Raw_lat$Date[i]-minDay))
    }
    else Raw_lat$EmergeDay[i] = Raw_lat$EmergeDay[i-1]

    if(i != nrow(Raw_lat)){
      if(Raw_lat$root[i] == Raw_lat$root[i+1]) time = time +timeStep
      else time =0
    }
  }

  setwd(storewd)
  write.csv(Raw_lat, file=paste(ExpName, "_Raw_lat.csv", sep=""))
  }
}

RITEtoData(path, setupFile, timeStep, ExpName, SorW)

#DataToElongation: use dataset from previous to calculate elongation
Data_lat <- NA
Data_main <- read.csv("Direction_Raw_main.csv")
timeResolution <- 20
timeStep <- 20
expName <- "Direction"

DataToElongation <- function(Data_lat, Data_main, timeResolution, expName){
  Data_main <- Data_main[,c("Root", "root", "Plate", "Treatment", "Genotype", "Date", "Time", "length", "direction")]

  #Add timeCategories
  require(dplyr)
  require(plyr)
  Data_main$timeCategory <- 0
  Data_main<- mutate(Data_main, timeCategory = round(Time/timeResolution))
  Data_main$GrowthPerHour <- with(Data_main, ave(length, Plate, Root, root,
                                                 FUN=function(x) c(NA,NA,NA,NA,NA,NA,diff(x,12)/((timeStep/60)*12))))

  Data_main$DirDiff <- with(Data_main, ave(direction, Plate, Root, root,
                                                 FUN=function(x) c(NA,diff(x,2))))

  Main_summary <- ddply(Data_main, c("Treatment", "Genotype", "timeCategory"), summarise,
                        N = sum(!is.na(length)),
                        Length = mean(length, na.rm=TRUE),
                        Dir = mean(direction, na.rm = TRUE),
                        Dirdiff = mean(DirDiff, na.rm=TRUE),
                        GR = mean(GrowthPerHour, na.rm = TRUE),
                        SE_Length =  sd(length, na.rm = TRUE)/ sqrt(N),
                        SE_GR =  sd(GrowthPerHour, na.rm = TRUE)/ sqrt(N),
                        SE_Dir = sd(direction, na.rm=TRUE)/sqrt(N),
                        SE_Dirdiff = sd(DirDiff,na.rm=TRUE)/sqrt(N),
                        HAS = mean(timeCategory) * (timeResolution/60))
  Main_summary$timeCategory <- as.numeric(Main_summary$timeCategory)


  if(!is.na(Data_lat)){
  Data_lat <- Data_lat[,c("root_name", "root", "parent", "Plate", "Treatment","Genotype","Date", "EmergeDay", "Time", "length")]
  Data_lat$timeCategory <-0
  Data_lat <- mutate(Data_lat, timeCategory = ceiling(Time/timeResolution))
  Data_lat$GrowthPerHour <- with(Data_lat, ave(length, root, Plate, parent,
                                               FUN=function(x) c(NA,NA,NA,NA,NA,NA,diff(x,12)/((timeStep/60)*12))))

  Lat_summary <- ddply(Data_lat, c("Treatment", "Genotype", "timeCategory", "EmergeDay"), summarise,
                        N = sum(!is.na(length)),
                        Length = mean(length, na.rm=TRUE),
                        GR = mean(GrowthPerHour, na.rm = TRUE),
                        SE_Length =  sd(length, na.rm = TRUE)/ sqrt(N),
                        SE_GR =  sd(GrowthPerHour, na.rm = TRUE)/ sqrt(N),
                       HAS = mean(timeCategory) * (timeResolution/60))
  Lat_summary$timeCategory <- as.numeric(Lat_summary$timeCategory)
  write.csv(Data_lat, file=paste(expName, "_Elongation_lat.csv", sep=""))
  write.csv(Lat_summary, file=paste(expName, "_ElongationSummary_lat.csv", sep=""))
  }

  write.csv(Data_main, file=paste(expName, "_Elongation_main.csv", sep=""))
  write.csv(Main_summary, file=paste(expName, "_ElongationSummary_main.csv", sep=""))

}

DataToElongation(Data_lat, Data_main, timeResolution, expName)

#DataToDirection: use dataset from previous to calculate direction

Data_main <- read.csv("Direction_Raw_main.csv")
timeResolution <- 20
timeStep <- 20
expName <- "Direction"

DataToDirection <- function(Data_lat, Data_main, timeResolution, expName){
  Data_main <- Data_main[,c("Root", "root", "Plate", "Treatment", "Genotype", "Date", "Time", "length", "direction", "vectorAngle", "nodeDirection5", "nodeDirection10", "nodeDirectionAv")]

  #Add timeCategories
  require(dplyr)
  require(plyr)
  Data_main$timeCategory <- 0
  Data_main<- mutate(Data_main, timeCategory = round(Time/timeResolution))
  Data_main$GrowthPerHour <- with(Data_main, ave(length, Plate, Root, root,
                                                 FUN=function(x) c(NA,NA,NA,diff(x,6)/((timeStep/60)*6))))

  Data_main$DirDiff <- with(Data_main, ave(direction, Plate, Root, root,
                                           FUN=function(x) c(NA,diff(x,2))))
  Data_main$AngDiff <- with(Data_main, ave(vectorAngle, Plate, Root, root,
                                           FUN=function(x) c(NA,diff(x,2))))
  Data_main$DirAvDiff <- with(Data_main, ave(nodeDirectionAv, Plate, Root, root,
                                           FUN=function(x) c(NA, NA, NA ,diff(x,6))))

  Main_summary <- ddply(Data_main, c("Treatment", "Genotype", "timeCategory"), summarise,
                        N = sum(!is.na(length)),
                        Length = mean(length, na.rm=TRUE),
                        Dir = mean(direction, na.rm = TRUE),
                        Ang = mean(vectorAngle, na.rm = TRUE),
                        Dir5 = mean(nodeDirection5, na.rm = TRUE),
                        Dir10 = mean(nodeDirection10, na.rm = TRUE),
                        DirAv = mean(nodeDirectionAv, na.rm = TRUE),
                        Dirdiff = mean(DirDiff, na.rm=TRUE),
                        Angdiff = mean(AngDiff, na.rm=TRUE),
                        DirAvdiff =mean(DirAvDiff, na.rm = TRUE),
                        GR = mean(GrowthPerHour, na.rm = TRUE),
                        SE_Length =  sd(length, na.rm = TRUE)/ sqrt(N),
                        SE_GR =  sd(GrowthPerHour, na.rm = TRUE)/ sqrt(N),
                        SE_Dir = sd(direction, na.rm=TRUE)/sqrt(N),
                        SE_Ang = sd(vectorAngle, na.rm = TRUE)/sqrt(N),
                        SE_Dir5 = sd(nodeDirection5, na.rm = TRUE)/sqrt(N),
                        SE_Dir10 = sd(nodeDirection10, na.rm = TRUE)/sqrt(N),
                        SE_DirAv = sd(nodeDirectionAv, na.rm = TRUE)/sqrt(N),
                        SE_Dirdiff = sd(DirDiff,na.rm=TRUE)/sqrt(N),
                        SE_Angdiff = sd(AngDiff,na.rm=TRUE)/sqrt(N),
                        SE_DirAvdiff = sd(DirAvDiff,na.rm=TRUE)/sqrt(N),
                        HAS = mean(timeCategory) * (timeResolution/60))
  Main_summary$timeCategory <- as.numeric(Main_summary$timeCategory)


  if(!is.na(Data_lat)){
    Data_lat <- Data_lat[,c("root_name", "root", "parent", "Plate", "Treatment","Genotype","Date", "EmergeDay", "Time", "length")]
    Data_lat$timeCategory <-0
    Data_lat <- mutate(Data_lat, timeCategory = ceiling(Time/timeResolution))
    Data_lat$GrowthPerHour <- with(Data_lat, ave(length, root, Plate, parent,
                                                 FUN=function(x) c(NA,NA,NA,NA,NA,NA,diff(x,12)/((timeStep/60)*12))))

    Lat_summary <- ddply(Data_lat, c("Treatment", "Genotype", "timeCategory", "EmergeDay"), summarise,
                         N = sum(!is.na(length)),
                         Length = mean(length, na.rm=TRUE),
                         GR = mean(GrowthPerHour, na.rm = TRUE),
                         SE_Length =  sd(length, na.rm = TRUE)/ sqrt(N),
                         SE_GR =  sd(GrowthPerHour, na.rm = TRUE)/ sqrt(N),
                         HAS = mean(timeCategory) * (timeResolution/60))
    Lat_summary$timeCategory <- as.numeric(Lat_summary$timeCategory)
    write.csv(Data_lat, file=paste(expName, "_Elongation_lat.csv", sep=""))
    write.csv(Lat_summary, file=paste(expName, "_ElongationSummary_lat.csv", sep=""))
  }

  write.csv(Data_main, file=paste(expName, "_Elongation_main.csv", sep=""))
  write.csv(Main_summary, file=paste(expName, "_ElongationSummary_main.csv", sep=""))

}

DataToDirection(Data_lat, Data_main, timeResolution, expName)

#ElongationToGraph:

#Data_lat <- read.csv("_Elongation_lat.csv")
Data_main<- read.csv("Direction_Elongation_main.csv")
#Summary_lat<- read.csv("7d_ElongationSummary_lat.csv")
Summary_main<- read.csv("Direction_ElongationSummary_main.csv")
#SavePath <- "~/Iko/PhD/Data/Camera project/20170303_B2B3/Graphs/" #path where csvfiles are
#timeResolution <- 120

#ElongationToGraph <- function(Data_main, Summary_main, SavePath, timeResolution){

  require(dplyr)
  require(plyr)
  require(dplyr)
  require(plyr)
  Main_smooth <- ddply(Data_main, c("Treatment", "Genotype", "timeCategory", "root"), summarise,
                       N = sum(!is.na(length)),
                       Length = mean(length, na.rm=TRUE),
                       Dir = mean(direction, na.rm = TRUE),
                       Ang = mean(vectorAngle, na.rm = TRUE),
                       Dir5 = mean(nodeDirection5, na.rm = TRUE),
                       Dir10 = mean(nodeDirection10, na.rm = TRUE),
                       DirAv = mean(nodeDirectionAv, na.rm = TRUE),
                       Dirdiff = mean(DirDiff, na.rm=TRUE),
                       Angdiff = mean(AngDiff, na.rm=TRUE),
                       GR = mean(GrowthPerHour, na.rm = TRUE),
                       SE_Length =  sd(length, na.rm = TRUE)/ sqrt(N),
                       SE_GR =  sd(GrowthPerHour, na.rm = TRUE)/ sqrt(N),
                       SE_Dir = sd(direction, na.rm=TRUE)/sqrt(N),
                       SE_Ang = sd(vectorAngle, na.rm = TRUE)/sqrt(N),
                       SE_Dir5 = sd(nodeDirection5, na.rm = TRUE)/sqrt(N),
                       SE_Dir10 = sd(nodeDirection10, na.rm = TRUE)/sqrt(N),
                       SE_DirAv = sd(nodeDirectionAv, na.rm = TRUE)/sqrt(N),
                       SE_Dirdiff = sd(DirDiff,na.rm=TRUE)/sqrt(N),
                       SE_Angdiff = sd(AngDiff,na.rm=TRUE)/sqrt(N),
                       HAS = mean(timeCategory) * (timeResolution/60))
  require(ggplot2)
  setwd(SavePath)
  for(i in 1:nlevels(Main_smooth$root)){
    Data <- subset(Main_smooth, Main_smooth$root == levels(Main_smooth$root)[i])
    Data <- droplevels(Data)
    if(nrow(Data) > 24){
    TreatmentName <- levels(Data$Treatment)
    timeCategoryName <- levels(Data$timeCategory)
    GenotypeName <- levels(Data$Genotype)
    PLOT <- paste(SavePath, levels(Main_smooth$root)[i], TreatmentName, GenotypeName, timeCategoryName,  sep="_")
    p <- ggplot(Data, aes(x=Genotype, y=Dirdiff))
    p + geom_line() + geom_point(size = 0.5) + xlab("accessions")
    ggsave(PLOT, device = "jpeg")
    }
  }

#ElongationToGraph(Data_main, Summary_main, SavePath, timeResolution)

#Overall graphs
require(ggplot2)

#nodeDirectionAv on timepionts
p <- ggplot(Summary_main[grep("69",Summary_main$timeCategory),], aes(x=Genotype, y=(DirAv*-1), colour = Treatment))
p  + geom_point(size = 3.5) + xlab("accessions") + ylab("nodeDirectionAv 24hrs") +
  geom_errorbar(aes(ymin= (DirAv*-1)- SE_DirAv, ymax=(DirAv*-1)+SE_DirAv, width=0.4))



#nDAD on time points
p <- ggplot(Summary_main[grep("66",Summary_main$timeCategory),], aes(x=Genotype, y=(DirAvdiff), colour = Treatment))
p  + geom_point(size = 3.5) + xlab("accessions") + ylab("nDirAvdiff 23hrs") +
  geom_errorbar(aes(ymin= (DirAvdiff)- SE_DirAvdiff, ymax=(DirAvdiff)+SE_DirAvdiff, width=0.4))



#Growth rate on time points
p <- ggplot(Summary_main[grep("49",Summary_main$timeCategory),], aes(x=Genotype, y=GR, colour = Treatment))
p  + geom_point(size = 3.5) + xlab("accessions") + ylab("GR 15hrs") +
  geom_errorbar(aes(ymin= GR- SE_GR, ymax=GR+SE_GR, width=0.4))

#nodeDirectionAv
  p <- ggplot(Summary_main[grep("Col",Summary_main$Genotype),], aes(x=timeCategory*20/60, y=(DirAv*-1), colour = Treatment))
  p  + geom_point(size = 0.1) + xlab("hours after start treatment") + ylab("Col0 nodeDirectionAv") + xlim(0,24) +geom_line(aes(group=Treatment))+
    geom_errorbar(aes(ymin= (DirAv*-1)- SE_DirAv, ymax=(DirAv*-1)+SE_DirAv, width=.1))

  #nodeDirection10
  p <- ggplot(Summary_main[grep("Col",Summary_main$Genotype),], aes(x=timeCategory*20/60, y=(Dir10*-1), colour = Treatment))
  p  + geom_point(size = 0.1) + xlab("hours after start treatment") + ylab("Col0 nodeDirection10") + xlim(0,24) +geom_line(aes(group=Treatment))+
    geom_errorbar(aes(ymin= (Dir10*-1)- SE_Dir10, ymax=(Dir10*-1)+SE_Dir10, width=.1))

  #nodeDirection5
  p <- ggplot(Summary_main[grep("Col",Summary_main$Genotype),], aes(x=timeCategory*20/60, y=(Dir5*-1), colour = Treatment))
  p  + geom_point(size = 0.1) + xlab("hours after start treatment") + ylab("Col0 nodeDirection5") + xlim(0,24) +geom_line(aes(group=Treatment))+
    geom_errorbar(aes(ymin= (Dir5*-1)- SE_Dir5, ymax=(Dir5*-1)+SE_Dir5, width=.1))

  #nodeDirAvDiff
  p <- ggplot(Summary_main[grep("Col",Summary_main$Genotype),], aes(x=timeCategory*20/60, y=(DirAvdiff), colour = Treatment))
  p  + geom_point(size = 0.1) + xlab("hours after start treatment") + ylab("Col0 nodeDirAv Change") + xlim(0,24) +geom_line(aes(group=Treatment))+
    geom_errorbar(aes(ymin= (DirAvdiff)- SE_DirAvdiff, ymax=(DirAvdiff)+SE_DirAvdiff, width=.1))

  #Elongation
  p <- ggplot(Summary_main[grep("Col",Summary_main$Genotype),], aes(x=timeCategory*20/60, y=GR, colour = Treatment))
  p  + geom_point(size = 0.1) + xlab("hours after start treatment") + ylab("Col0 Growth rate") + xlim(0,24) +geom_line(aes(group=Treatment))+
    geom_errorbar(aes(ymin= GR- SE_GR, ymax=GR+ SE_GR, width=.1))
