setwd("")
SetupFile <- read.csv("Tilting_setup.csv")
SRFile <- read.csv("testTilting.csv")
ExpName <- "exad_titling"


DataInSegments <- function(SetupFile, SRFile, ExpName){

  #Prepare markerfile
  DP<- (data.frame (date= substr(SRFile$image, start=1, stop = 8),
                    Plate= as.numeric(substr(SRFile$image, start=10, stop = 12)) ))#Separate date & plate no
  SRFile <- cbind(SRFile, DP[,1:2]) #bind date & plate to the current dataset

  #determine root number (NR)
  NR <- t(as.data.frame(strsplit(as.character(SRFile$root_name), "\\_")))
  colnames(NR) <-  c("ontology", "Root") #Add column names
  rownames(NR) <- NULL
  SRFile <- cbind(SRFile, NR) #bind root number to the current dataset

  #determine time (H= hours)
  require(car)
  H<- recode(SRFile$date, "c(2)='24'; c(3)='48'; c(4)='72'")
  SRFile <- cbind(SRFile, H) #bind time to the current dataset

  #match marker
  for(i in 1:nrow(SRFile)){
    for (b in 1:nrow(SetupFile)){
      if (as.numeric(SRFile$Plate[i]) == SetupFile$Plate[b]){
        if(SRFile$Root[i] == SetupFile$Root[b]) {
          SRFile$Genotype[i] <- as.character(SetupFile$genotypes[b])
          SRFile$Line[i] <- as.character(SetupFile$lines[b])
          SRFile$Condition[i] <- as.character(SetupFile$conditions[b])
          SRFile$Treatment[i] <- as.character(SetupFile$Treatment[b])
        }
      }
    }
  }

  #summarize node data
  write.csv(SRFile, file=paste(ExpName, "_FullDataSet.csv", sep=""))

  require(plyr)
  Summary <- ddply(SRFile, c("H","Treatment", "Genotype"), summarise,
                   N=sum(!is.na( growth)),
                   length = mean( growth, na.rm = TRUE),
                   vect_angle = mean( vector_angle, na.rm = TRUE),
                   node_angle = mean( av_node_angle, na.rm = TRUE),
                   se_length = sd( growth, na.rm = TRUE)/sqrt(N),
                   se_vectAngle = sd( vector_angle, na.rm = TRUE)/ sqrt(N),
                   se_nodeAngle = sd( av_node_angle, na.rm = TRUE)/ sqrt(N))

  write.csv(Summary, file=paste(ExpName, "_SummaryTimepoints.csv", sep=""))


}
 DataInSegments(SetupFile, SRFile, ExpName)



###########################################################################################################

#Graphs and Statistics
setwd("~/Desktop")
Data <- read.csv("exad_titling_FullDataSet.csv")
Data$Genotype <- factor(Data $ Genotype, level=c ("Col-0", "exad-1", "exad-2"))
Data$H <- factor(Data$H, level= c( "24", "48","72"))
Data$Treatment <- factor(Data $ Treatment, level=c("Control", "NaCl100mM"))
Data$interaction <- interaction(Data $ Genotype, Data$H, Data$Treatment, Data$Plate)

require(ggplot2)

p <- ggplot(Data, aes(x=Treatment, y= vector_angle*-1, Fill = Genotype)) + geom_boxplot(outlier.colour = NA) + scale_fill_discrete() + ylab("Root VectorAngle") +  xlab("Hours post-stress")+ scale_x_discrete(limits= c("Control", "NaCl100mM"))
p  + theme(axis.text.x = element_text(angle=0, hjust = 0.5, vjust = 0, size = 10)) + theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(angle=90, size = 12)) + facet_grid(H~.)
#+ geom_boxplot(notch = F)

# new
require(ggplot2)
p <- ggplot(Data, aes(x = Treatment, y= vector_angle* -1,fill = Genotype)) + geom_boxplot(outlier.colour = NA)   + scale_fill_manual(values = c("Col-0" = "palegreen3", "exad-1" = "steelblue1", "exad-2" = "steelblue3", "exad-CP"= "steelblue4"))  + ylab("Root Vector Angle (degree)") + scale_x_discrete(limits= c("Control", "NaCl100mM"))
p + scale_x_discrete()  + theme(axis.text.x = element_text(angle=0, hjust = .5, vjust = .5, size = 10)) +theme(axis.title.x = element_blank()) + facet_grid(H ~ . )+expand_limits(y= c(-60, 60)) #+ geom_boxplot(notch = FALSE)

require(ggplot2)
p <- ggplot(Data, aes(x = Treatment, y= growth, fill = Genotype)) + geom_boxplot(outlier.colour = NA)   + scale_fill_manual(values = c("Col-0" = "palegreen3", "exad-1" = "steelblue1", "exad-2" = "steelblue3", "exad-CP"= "steelblue4"))  + ylab("Root Growth Rate (cm/h)") + scale_x_discrete(limits= c("Control", "NaCl100mM"))
p + scale_x_discrete()  + theme(axis.text.x = element_text(angle=0, hjust = .5, vjust = .5, size = 10)) +theme(axis.title.x = element_blank()) +expand_limits(y= c(0, 0.75)) + facet_grid(H ~ . ) #+ geom_boxplot(notch = FALSE)

p <- ggplot(Data, aes(x=Treatment, y= growth*-1, fill = Genotype )) + geom_boxplot(outlier.colour = NA) + scale_fill_discrete() + ylab("Root GrowthRate") +  xlab("Hours post-stress")+ scale_x_discrete(limits= c("Control", "NaCl100mM"))
p  + theme(axis.text.x = element_text(angle=0, hjust = 0.5, vjust = 0, size = 10)) + theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_text(angle=90, size = 12)) + facet_grid(H~.)



#Significance
Data <- read.csv("exad_titling_FullDataSet.csv")
Data$Condition <- as.factor(Data$Condition)
Data$Line <- as.factor(Data$Line)
Data$Plate <- as.factor(Data$Plate)
Data$H <- as.factor(Data$H)

qqnorm(Data$vector_angle)
qqline(Data$vector_angle)

require(nlme)
Model <- lme(fixed=growth~Line*Condition,random=~1|Plate, data = Data)
anova(Model)
summary(Model)

int <- interaction(Data$Line, Data$Condition, Data$H)
Data <- cbind(Data, int)
Data$int <- droplevels(Data$int)
head(Data)

Model <- lm(growth~int, Data)
anova(Model)

library(multcomp)
require(multcomp)
contrasts<- c("Col0.0.24 - exad1.0.24 = 0",
              "Col0.0.24 - exad2.0.24 = 0",
              "Col0.100.24 - exad1.100.24 = 0",
              "Col0.100.24 - exad2.100.24 = 0",
              "Col0.0.48 - exad1.0.48 = 0",
              "Col0.0.48 - exad2.0.48 = 0",
              "Col0.100.48 - exad1.100.48 = 0",
              "Col0.100.48 - exad2.100.48 = 0",
              "Col0.0.72 - exad1.0.72 = 0",
              "Col0.0.72 - exad2.0.72 = 0",
              "Col0.100.72 - exad1.100.72 = 0",
              "Col0.100.72 - exad2.100.72 = 0")

summary(glht(Model, linfct = mcp(int = contrasts)))


