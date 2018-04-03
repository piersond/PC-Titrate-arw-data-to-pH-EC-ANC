###########################
# File: COSA_compiler_2016.R
# Description: Compile COSA data (pH, ANC, EC) from raw export files
# Date: January 7, 2016
# Revised:  
# Author: dpierson02 -> Derek Pierson
# Notes: 
# 
###########################

#Load packages
library("ggplot2")

#Find and list ANC raw data files
setwd("C:/Users/dpierson02/Google Drive/Biogeochemistry Lab/Run Files/Raw ANC/")

#######################################
#Number of previous runs to analyze (default = 1)
#Increase number if multiple sets of run data have been added since last running script
Runs <- 1
#######################################  



#Find all folders in the working directory and remove all non-.txt files
  ANCFile.list <- list.files()
  Files.ANC <- ANCFile.list[grepl(".TXT",ANCFile.list) == TRUE];
  Files.ANC <- as.data.frame(Files.ANC)
  Files.ANC$Date <- strtrim(Files.ANC[,1], 6) 
  Files.ANC$Modified <- file.info(as.character(Files.ANC[,1]))$mtime  
  Files.ANC$Modified <- as.POSIXct(Files.ANC$Modified, format="%Y/%m/%d %I:%M:%S")
  #Files.ANC$Modified <- Files.ANC$Modified - (60*60)  #<-- Fixes time offest

  
  #######################################
  #Build Samples List?
  Answer <- "Y" # Y or N
  #######################################
  
  if(Answer == "N") {
    
    samples <- NULL
    for (i in 4500:nrow(Files.ANC))
    {
      tbl <- read.table(as.character(Files.ANC[i,1]), sep=" ", header = FALSE, fill = TRUE,) 
      x <- NULL
      x <- as.numeric(is.na(tbl[1,2]))
      y <- as.character(tbl[1,1])
      tbl[,1] <- ifelse(x>0,y,paste0(y,tbl[,2]))
      ifelse(x>0,"ok",tbl[,2] <- NA)
      nm <- NULL
      nm <- as.data.frame(as.character(tbl[1,1]))
      nm[1,2] <- i
      samples <- rbind(samples,nm)
      
    }
    colnames(samples) <- c("Sample ID", "i")
  }  
  
  
  
   
#Set ANC analysis volume
Vol <- 35  #35 ml is default method (25 ml used as alternate method for high EC samples)



#Find and list raw EC-pH export files
setwd("C:/Users/dpierson02/Google Drive/Biogeochemistry Lab/Run Files/Raw EC/")

  #Find all folders in the working directory and remove all non-.txt files
  FileCount.EC <- list.files()
  Files.EC <- FileCount.EC[grepl(".TXT",FileCount.EC) == TRUE];
  #Compile.Record <- read.csv("COSA_Compiler_Record.csv")

#Set variables for summary files
RunSummary <- NULL
Summary <- NULL

#Works back through specified number of "runs" 
for (c in (length(Files.EC)):(length(Files.EC)-(runs-1)) 
  {

    #Load EC-pH summary file
      setwd("C:/Users/dpierson02/Google Drive/Biogeochemistry Lab/Run Files/Raw EC/")
      file <- Files.EC[c]
      #ifelse(grepl(file,Compile.Record[,2]), next, paste0("Running Analysis on file number: ",c))
      
    #Import EC and pH data from COSA summary export file
    tbl.EC <- read.fwf(
      file=file,
      skip=0,
      widths=c(42, 20, 11, 26, 10, 20))
    
    ifelse(nrow(tbl.EC)<5,next,"ok")
    
    tbl.EC <- tbl.EC[,2:6]
    colnames(tbl.EC) <- c("SampleID", "Date", "Time", "EC", "pH")
    
    tbl.EC$Date<- strtrim(tbl.EC$Date, 8) 
    tbl.EC$Date <- gsub(" ","",tbl.EC$Date)
    #tbl.EC$Date <- ifelse(strtrim(tbl.EC$Date, 2) == "12" | strtrim(tbl.EC$Date, 2) == "11" | strtrim(tbl.EC$Date, 2) == "10", tbl.EC$Date, paste0("0",tbl.EC$Date))                 
    
    tbl.EC$Time <- gsub(" ","",tbl.EC$Time) 
    tbl.EC$Time <- gsub("PM"," PM",tbl.EC$Time) 
    tbl.EC$Time <- gsub("AM"," AM",tbl.EC$Time) 
    #tbl.EC$Time <- paste0("0",tbl.EC$Time)
    tbl.EC$Time <- gsub("00:","12:",tbl.EC$Time) 
  
    tbl.EC$modified <- as.POSIXct(as.character(paste(tbl.EC$Date, tbl.EC$Time)), format="%m/%d/%y %I:%M:%S %p")
  

    #Exceptions for missing samples in EC table
      tbl.EC <- subset(tbl.EC, tbl.EC[,6] != as.POSIXct("2015-11-23 11:33:10", format="%Y-%m-%d %I:%M:%S"))
      if(tbl.EC[nrow(tbl.EC),2] == "9/23/15") {tbl.EC <- tbl.EC[1:(nrow(tbl.EC)-1),]}
      if(tbl.EC[1,2] == "1/5/16") {tbl.EC <- tbl.EC[c(1:20,22:26,28:nrow(tbl.EC)),]}
      
      
    #Time change adjustments  
      #if(c <= 136) {tbl.EC$modified <- tbl.EC$modified - 60*60}   
       
    #Set titration volume from first sample ID
      Vol <- 35  #35 ml is default method (25 ml used as alternate method for high EC samples)
      if(as.character(gsub(" ", "",tbl.EC[1,1])) == "V25") {Vol <- 25}

    
    #RUN ANC ANALYSIS 
      setwd("C:/Users/dpierson02/Google Drive/Biogeochemistry Lab/Run Files/Raw ANC/")
      ANCsummary <- NULL
   
    for (i in 1:nrow(tbl.EC))
      {
        

      
        file.time1 <- tbl.EC[i,6] - 300  # - (60*60*5)
        file.time2 <- tbl.EC[i,6] + 300  # - (60*60*5)
        
        file.df <- Files.ANC[which(Files.ANC$Modified > file.time1 & Files.ANC$Modified < file.time2),]
        if(nrow(file.df) < 1 & i != 1 & i != nrow(tbl.EC)) {     
          file.time1 <- tbl.EC[i-1,6]
          file.time2 <- tbl.EC[i+1,6]
          file.df <- Files.ANC[which(Files.ANC$Modified > file.time1 & Files.ANC$Modified < file.time2),]
        } 

        
        tbl <- read.table(as.character(file.df[1,1]), sep=" ", header = FALSE, fill = TRUE) 
        tbl[,1] <- paste(tbl[,1],tbl[,2],tbl[,3],tbl[,4])
        tbl[,1] <- gsub(" NA","", tbl[,1])
        tbl <- tbl[,c(1,7:ncol(tbl))]
        x <- NULL
        x <- as.numeric(is.na(tbl[1,2]))
        y <- as.character(tbl[1,1])
        tbl[,1] <- ifelse(x>0,y,paste0(y,tbl[,2]))
        ifelse(x>0,"ok",tbl[,2] <- NA)
        
              if(nrow(tbl) < 5){
                data <- NULL
                data <- as.data.frame(tbl[1,1])
                colnames(data)[1] <- "SampleID"
                data$pH <- "0"
                data$ANC <- "0"
                data$CaCO3 <- "0"
                data$simpleCaCO3 <- "0"
                data$Rsqd <- "0"
                data$SampleVol <- "0"
                data$AcidVol <- "0"
                data$Intercept <- "0"
                data$Slope <- "0"
                data$Points <- "0"
                data$i <- i
                ANCsummary <- rbind(ANCsummary,data)
                next
              }
        
        
        
        #Table Corrections
        for(j in 1:nrow(tbl)){tbl[j,41] <- ifelse(is.na(tbl[j,41]) == TRUE, tbl[j,40], tbl[j,41])}
        for(j in 1:nrow(tbl)){tbl[j,25] <- ifelse(is.na(tbl[j,25]) == TRUE, tbl[j,24], tbl[j,25])}
        for(j in 1:nrow(tbl)){tbl[j,29] <- ifelse(is.na(tbl[j,29]) == TRUE, tbl[j,28], tbl[j,29])}
        for(j in 1:nrow(tbl)){tbl[j,42] <- ifelse(is.na(tbl[j,42]) == TRUE, tbl[j,41], tbl[j,42])}
        for(j in 1:nrow(tbl)){tbl[j,26] <- ifelse(is.na(tbl[j,26]) == TRUE, tbl[j,25], tbl[j,26])}
        
        if(ncol(tbl) > 44) {for(j in 1:nrow(tbl)){tbl[j,45] <- ifelse(is.na(tbl[j,45]) == TRUE, tbl[j,44], tbl[j,45])}}
        if(ncol(tbl) > 49) {for(j in 1:nrow(tbl)){tbl[j,50] <- ifelse(is.na(tbl[j,50]) == TRUE, tbl[j,49], tbl[j,50])}}

        tbl <- tbl[,1:(ncol(tbl)-2)]
        
        tbl <- tbl[colSums(is.na(tbl))<5]
        tbl <- tbl[complete.cases(tbl),]
        
        if(ncol(tbl) < 7){
          ifelse(tbl[2,2] == tbl[2,3] , tbl <- tbl[,c(1,2,4)], "OK")
        
          ifelse(tbl[1,2] == "HALF", tbl <- tbl[,c(1,3,4)], "OK")
        
          ifelse(ncol(tbl) == 4, colnames(tbl)[2:3] <- c("Va","pH"), colnames(tbl)[2:3] <- c("Va","pH"))
        
          ifelse(ncol(tbl) == 4,df <- tbl[,2:3],df <- tbl[,2:3])
        }
        
        if(ncol(tbl) > 6)
          {
            tbl <- tbl[,c(1,2,6)]
            colnames(tbl)[2:3] <- c("Va","pH")
            df <- tbl[,2:3]
            }
        
        #Grans Analysis
        row.names(df) <- NULL
        df$TVol <- df[nrow(df),1] + Vol
        df$F1 <- df$TVol*(10^-df[,2]) 
        
        #Titration Plot
        ggplot(data=df, aes(x=Va,y=pH,group=1)) + geom_point() + geom_line()
        
        #Grans Plot 
        ggplot(data=df, aes(x=Va,y=F1,group=1)) + geom_point() + geom_line() 
        
        dfz <- NULL
        
        #Best at 5.55 & 4.6 as of June 10th, 2014
        dfz <- df[which(df$pH < 5.7 & df$pH > 4.5),]
        
              if(nrow(dfz) < 3){
                data <- NULL
                data <- as.data.frame(tbl[1,1])
                colnames(data)[1] <- "SampleID"
                data$pH <- "0"
                data$ANC <- "0"
                data$CaCO3 <- "0"
                data$simpleCaCO3 <- "0"
                data$Rsqd <- "0"
                data$SampleVol <- "0"
                data$AcidVol <- "0"
                data$Intercept <- "0"
                data$Slope <- "0"
                data$Points <- "0"
                data$i <- i
                ANCsummary <- rbind(ANCsummary,data)
                next
              }
        
        #Point included for Grans analysis
        ggplot(data=dfz, aes(x=Va,y=F1,group=1)) + geom_point() + geom_line() 
        
        rsq <- cor(x=dfz$Va, y=dfz$F1, use="pairwise.complete.obs") 
        
        lmz = lm(F1 ~ Va, data=dfz) 
        
        coeffs <- coef(lmz)
        
        #ANC calculation
        ANC <- (((coeffs[1]/coeffs[2])*1000*0.01)/Vol)*-1000  #<-- ANC in ueq/L
        ANC
        
        CaCO3 <- ((coeffs[1]/coeffs[2])*-50044*0.01)/Vol      #<-- Alk as mg/L CaCO3
        CaCO3
        simpleCaCO3 <- ANC/20
        
        data <- NULL
        data <- as.data.frame(tbl[1,1])
        colnames(data)[1] <- "SampleID"
        data$pH <- (tbl$pH[1])
        data$ANC <- round(ANC,3)
        data$CaCO3 <- round(CaCO3,3)
        data$simpleCaCO3 <- simpleCaCO3
        data$Rsqd <- rsq
        data$SampleVol <- Vol
        data$AcidVol <- df[nrow(df),1]
        data$Intercept <- coeffs[1]
        data$Slope <- coeffs[2]
        data$Points <- nrow(dfz)
        data$i <- i
        ANCsummary <- rbind(ANCsummary,data)
      }
    RunSummary <- tbl.EC
    RunSummary$ANC <- as.numeric(ANCsummary$ANC) #* 0.95  #<--NOTICE: % Modifier for ANC values based on standard shift
    RunSummary$RunFile <- file 
    RunSummary <- RunSummary[,c(1:3,5,7,4,8)]
    RunSummary$ANCvol <- Vol
    RunSummary$i <- ANCsummary$i
    RunSummary$ECfile <- c 
    
    #Export each COSA run file
    write.csv(RunSummary,file=paste0("C:/Users/dpierson02/Google Drive/Biogeochemistry Lab/Run Files/COSA Run Files/","COSA_",strtrim(file,6),".csv"))    
    
  Summary <- rbind(Summary,RunSummary)
        
}


#Export summary file
#write.csv(Summary,file=paste0("C:/Users/dpierson02/Google Drive/Biogeochem Lab/Water Chemistry/Run Files/COSA/Runs/","COSA_Summary_",Sys.Date(),".csv"))
