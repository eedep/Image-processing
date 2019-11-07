# prior to analysis, tif images were converted using imageJ - command: run("Enhance Contrast", "saturated=0.35");, and saved as JPG files

install.packages("tidyverse")
install.packages("imager")
install.packages("ggsci")
install.packages("ggpubr")

library(tidyverse)
library(ggsci)
library(ggpubr)
library(imager)


path = "F:/Charite-Bernhard (Backup 20190930)/Bernhard/20190918_Pascales Revision/Converted Images"
setwd(path)


animal.file <- paste(path,"/","APP23p40KO Experimental list_histo.csv",sep="")
animal.info <- read.csv2(animal.file)
animal.info[,1] <- as.character(animal.info[,1])
animal.info[,2] <- as.character(animal.info[,2])
animal.info[,3] <- as.character(animal.info[,3])
animal.info[,4] <- as.character(animal.info[,4])
animal.info[,5] <- as.character(animal.info[,5])
animal.info[,6] <- as.character(animal.info[,6])
animal.info[,7] <- as.numeric(animal.info[,7])


Numextract <- function(string){
  unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}

hist.bins<-255
hist.breaks <- c(seq(0,1,1/hist.bins))


file.names <- dir(path, pattern=".jpg")

threshold <- NULL
image.info.df <- data.frame()
h.data <- data.frame()
row.names <- c()


      for (i in 1:length(file.names)) {
        
        actual.file <- load.image(file.names[i])
        actual.file.gr <- grayscale(actual.file)
        
        file.name <- strsplit(file.names[i], ".jpeg")

        thrshld <- round((1-mean(actual.file.gr))*hist.bins,digits = 0)
        threshold <- c(threshold,  thrshld)
        
        histogram <- hist(actual.file.gr, hist.breaks, plot=FALSE, right=TRUE, include.lowest=TRUE)
        image.info <- as.numeric(Numextract(file.name))
        histogram$density <- as.numeric(histogram$density/2.5)

        if(length(image.info)!=4){
          
          if(length(image.info) == 3){
            image.info <- c(image.info, "n.a.")
          } else
            if(length(image.info) == 2){
              image.info <- c(image.info, "n.a.", "n.a.")
            } else
              if(length(image.info) == 1) {
                image.info <- c(image.info, "n.a.", "n.a.", "n.a.")
                } else
                  image.info <- c("n.a.","n.a.", "n.a.", "n.a.")
        }

        h.data <- rbind.data.frame(h.data, rev(histogram$density))
        image.info.df <- rbind.data.frame(image.info.df, image.info)
        
        row.names <- c(row.names, file.names[i])
   

}


h.data <- cbind.data.frame(file.names,image.info.df,threshold,h.data)
rownames(h.data) <- row.names
Xnames <- paste("X",str_pad(seq(1,hist.bins),3,pad="0"),sep="")
colnames(h.data) <- c("File Name","Animal ID", "Slice No.", "Image No.", "Channel","Inverted Image Mean", Xnames)

 for (k in 1:length(h.data$`Channel`)){
         if(h.data[k,5] == "5"){
           h.data[k,5] <- "4G8 - Cy5 (647 nm)"
         } else
          if(h.data[k,5] == "3"){
            h.data[k,5] <- "BACE - Cy3 (532 nm)"
          } 
 
 }

h.data <- h.data %>%
  arrange(`Animal ID`,`Slice No.`,`Image No.`,`Channel`)

h.data <- cbind(h.data[1:5],`Experimental Group`="",`APP23 Genotype`="",`p40 Genotype`="", `Age in Days`=0, `Gender`="",h.data[6:length(h.data[1,])])

h.data$`Experimental Group` <- as.character(h.data$`Experimental Group`)
h.data$`APP23 Genotype` <- as.character(h.data$`APP23 Genotype`)
h.data$`p40 Genotype` <- as.character(h.data$`p40 Genotype`)
h.data$Gender <- as.character(h.data$Gender)

for(l in 1:length(h.data[,1])){
  
  for(ll in 1:length(animal.info$Animal.ID)){
    
    if(as.character(h.data$`Animal ID`[l]) == as.character(animal.info$Animal.ID[ll])){
      
      h.data$`Experimental Group`[l] <- paste(animal.info$Genotype.APP23[ll],animal.info$Genotype.p40[ll],sep=" ")
      h.data$`APP23 Genotype`[l] <- animal.info$Genotype.APP23[ll]
      h.data$`p40 Genotype`[l] <- animal.info$Genotype.p40[ll]
      h.data$`Gender`[l] <- animal.info$Gender[ll]
      h.data$`Age in Days`[l] <- animal.info$Age.in.Days[ll]
      if(str_detect(h.data$`File Name`[l],"cond")==TRUE){
        h.data$`Experimental Group`[l] <- "Secondary AB Control"
      }
    } 
  }
}

write.csv2(h.data,file= paste("h.data.csv"))
                      
                      #hist.bins<-255
                      #Xnames <- paste("X",str_pad(seq(1,hist.bins),3,pad="0"),sep="")
                      #h.data <- read.csv2("h.data.csv",sep=";")
                      #h.data <- h.data[,-1]
                      #colnames(h.data) <- c("File Name","Animal ID", "Slice No.", "Image No.", "Channel","Experimental Group","APP23 Genotype","p40 Genotype", "Age in Days", "Gender","Inverted Image Mean",Xnames) # optional for reading in the csv file

h.data.BACE <- subset(h.data, `Channel`=="BACE - Cy3 (532 nm)")
h.data.4G8 <- subset(h.data,`Channel`=="4G8 - Cy5 (647 nm)")

h.data.BACE <- h.data.BACE %>%
  gather(key, density, matches("^X\\d")) %>%
  group_by(`Animal ID`, `Slice No.`,`Image No.`) %>%
  arrange(`Animal ID`, `Slice No.`,`Image No.`)%>%
  ungroup()

integral.data.BACE <- h.data.BACE %>%
  gather(key, density, matches("^X\\d")) %>%
  group_by(`Animal ID`, `Slice No.`,`Image No.`,`Experimental Group`) %>%
  mutate(integral_hist=cumsum(density)) %>%
  select(-c(`density`)) %>%
  ungroup() 

h.data.4G8 <- h.data.4G8 %>%
  gather(key, density, matches("^X\\d")) %>%
  group_by(`Animal ID`, `Slice No.`,`Image No.`) %>%
  arrange(`Animal ID`, `Slice No.`,`Image No.`)%>%
  ungroup()

integral.data.4G8 <- h.data.4G8 %>%
  gather(key, density, matches("^X\\d")) %>%
  group_by(`Animal ID`, `Slice No.`,`Image No.`,`Experimental Group`) %>%
  mutate(integral_hist=cumsum(density)) %>%
  select(-c(`density`)) %>%
  ungroup() 

analysis.threshold.4G8 <- round(median(integral.data.4G8$`Inverted Image Mean`)-2*sd(integral.data.4G8$`Inverted Image Mean`),digits=0)
analysis.threshold.BACE <- round(median(integral.data.BACE$`Inverted Image Mean`)-2*sd(integral.data.BACE$`Inverted Image Mean`),digits=0)


quality.check1.4G8 <- subset(h.data[,1:11],`Channel`=="4G8 - Cy5 (647 nm)")
quality.check1.4G8$`Animal ID` <- factor(quality.check1.4G8$`Animal ID`)
quality.check1.4G8$`Slice No.` <- factor(quality.check1.4G8$`Slice No.`)
quality.check1.4G8$Channel <- factor(quality.check1.4G8$Channel)

quality.check1.4G8 <- quality.check1.4G8 %>%
  group_by(`Animal ID`,`Channel`) %>%
  mutate(`median_Animal`= median(`Inverted Image Mean`)) %>%
  mutate(`SD_Animal`= sd(`Inverted Image Mean`)) %>%
  group_by(`Animal ID`,`Slice No.`, `Channel`) %>%
  mutate(`median_Slice`= median(`Inverted Image Mean`)) %>%
  mutate(`SD_Slice`= sd(`Inverted Image Mean`)) %>%
  ungroup()
  
quality.check1.4G8 <- cbind(quality.check1.4G8, 'Quality_Test_A_passed'="",'Quality_Test_B_passed'="",'Overall Sample Quality'=0)
quality.check1.4G8$Quality_Test_A_passed <- as.character(quality.check1.4G8$Quality_Test_A_passed)
quality.check1.4G8$Quality_Test_B_passed <- as.character(quality.check1.4G8$Quality_Test_B_passed)
quality.check1.4G8$`Overall Sample Quality` <- as.numeric(quality.check1.4G8$`Overall Sample Quality`)


for(qc in 1:NROW(quality.check1.4G8)){

  if(is.na(quality.check1.4G8$SD_Animal[qc])){quality.check1.4G8$SD_Animal[qc] <- 0}
  if(is.na(quality.check1.4G8$SD_Slice[qc])){quality.check1.4G8$SD_Slice[qc] <- 0}
  
  }

for(qc in 1:NROW(quality.check1.4G8)){
  if(quality.check1.4G8$`Inverted Image Mean`[qc] >= (quality.check1.4G8$median_Slice[qc]-(2*quality.check1.4G8$SD_Slice[qc])) & quality.check1.4G8$`Inverted Image Mean`[qc] <= (quality.check1.4G8$median_Slice[qc]+(2*quality.check1.4G8$SD_Slice[qc]))){
    quality.check1.4G8$Quality_Test_A_passed[qc] <- "Yes"  
    quality.check1.4G8$`Overall Sample Quality`[qc] <- (quality.check1.4G8$`Overall Sample Quality`[qc]+50)}else
      quality.check1.4G8$Quality_Test_A_passed[qc] <- "No" 

 if(quality.check1.4G8$`Inverted Image Mean`[qc] >= (quality.check1.4G8$median_Animal[qc]-(2*quality.check1.4G8$SD_Animal[qc])) & quality.check1.4G8$`Inverted Image Mean`[qc] <= (quality.check1.4G8$median_Animal[qc]+(2*quality.check1.4G8$SD_Animal[qc]))){
    quality.check1.4G8$Quality_Test_B_passed[qc] <- "Yes"  
    quality.check1.4G8$`Overall Sample Quality`[qc] <- (quality.check1.4G8$`Overall Sample Quality`[qc]+50)}else
      quality.check1.4G8$Quality_Test_B_passed[qc] <- "No"
}

quality.check2.4G8 <- quality.check1.4G8 %>%
  group_by(`Animal ID`,`Channel`) %>%
  summarize(`Overall Quality`=mean(`Overall Sample Quality`)) %>%
  ungroup

write.csv2(quality.check1.4G8,file= paste("quality.check1.4G8.csv"))
write.csv2(quality.check2.4G8,file= paste("quality.check2.4G8.csv"))
Mean.Sample.Quality.4G8 <- mean(quality.check1.4G8$`Overall Sample Quality`)

plot.quality1.4G8 <- ggplot(quality.check1.4G8,aes(x=`Animal ID`,y=`Inverted Image Mean`, color=`Slice No.`, fill=`Channel`)) +
  geom_boxplot(notch=FALSE, width=0.2,size=1,shape="x",color="black",outlier.shape = NA)  +
  geom_jitter(width=0.2,alpha=0.35) + 
  geom_hline(yintercept = median(quality.check1.4G8$`Inverted Image Mean`)-1.96*sd(quality.check1.4G8$`Inverted Image Mean`), linetype="F1", color = "black")+
  geom_hline(yintercept = median(quality.check1.4G8$`Inverted Image Mean`)+1.96*sd(quality.check1.4G8$`Inverted Image Mean`), linetype="F1", color = "black")+
  scale_color_npg() +
  theme_bw()+ 
  ylab("mean image intensity (4G8)") +
  xlab("Animal ID")+
  facet_grid(rows=vars(`Channel`))+
  ylim(0,255) +
  theme(axis.text.x = element_text(angle = 90))+
  scale_size(range=c(1,3.5))

plot.quality2.4G8 <- ggplot(quality.check1.4G8,aes(x=`Slice No.`,y=`Inverted Image Mean`, color=`Slice No.`, fill=`Channel`)) +
  geom_jitter(width=0.2,alpha=0.5) + 
  geom_boxplot(notch=FALSE, width=0.2,size=1,shape="x",color="black",outlier.shape = NA)  +
  geom_hline(yintercept = median(quality.check1.4G8$`Inverted Image Mean`)-1.96*sd(quality.check1.4G8$`Inverted Image Mean`), linetype="F1", color = "black")+
  geom_hline(yintercept = median(quality.check1.4G8$`Inverted Image Mean`)+1.96*sd(quality.check1.4G8$`Inverted Image Mean`), linetype="F1", color = "black")+
  scale_color_npg() +
  theme_light()+ 
  ylab("mean image intensity (4G8)") +
  xlab("Slice No.")+
  facet_grid(cols=vars(`Channel`))+
  ylim(0,255)+
  scale_size(range=c(1,3.5))


png(paste(Sys.Date(),"plot.quality1.4G8.png",sep="_"),width = 10, height = 7, units = 'in', res = 300)
print(plot.quality1.4G8)     
dev.off() 

png(paste(Sys.Date(),"plot.quality2.4G8.png",sep="_"),width = 10, height = 7, units = 'in', res = 300)
print(plot.quality2.4G8)     
dev.off() 


quality.check1.BACE <- subset(h.data[,1:11],`Channel`=="BACE - Cy3 (532 nm)")
quality.check1.BACE$`Animal ID` <- factor(quality.check1.BACE$`Animal ID`)
quality.check1.BACE$`Slice No.` <- factor(quality.check1.BACE$`Slice No.`)
quality.check1.BACE$Channel <- factor(quality.check1.BACE$Channel)

quality.check1.BACE <- quality.check1.BACE %>%
  group_by(`Animal ID`,`Channel`) %>%
  mutate(`median_Animal`= median(`Inverted Image Mean`)) %>%
  mutate(`SD_Animal`= sd(`Inverted Image Mean`)) %>%
  group_by(`Animal ID`,`Slice No.`, `Channel`) %>%
  mutate(`median_Slice`= median(`Inverted Image Mean`)) %>%
  mutate(`SD_Slice`= sd(`Inverted Image Mean`)) %>%
  ungroup()

quality.check1.BACE <- cbind(quality.check1.BACE, 'Quality_Test_A_passed'="",'Quality_Test_B_passed'="",'Overall Sample Quality'=0)
quality.check1.BACE$Quality_Test_A_passed <- as.character(quality.check1.BACE$Quality_Test_A_passed)
quality.check1.BACE$Quality_Test_B_passed <- as.character(quality.check1.BACE$Quality_Test_B_passed)
quality.check1.BACE$`Overall Sample Quality` <- as.numeric(quality.check1.BACE$`Overall Sample Quality`)


for(qc in 1:NROW(quality.check1.BACE)){
  
  if(is.na(quality.check1.BACE$SD_Animal[qc])){quality.check1.BACE$SD_Animal[qc] <- 0}
  if(is.na(quality.check1.BACE$SD_Slice[qc])){quality.check1.BACE$SD_Slice[qc] <- 0}
  
}

for(qc in 1:NROW(quality.check1.BACE)){
  if(quality.check1.BACE$`Inverted Image Mean`[qc] >= (quality.check1.BACE$median_Slice[qc]-(2*quality.check1.BACE$SD_Slice[qc])) & quality.check1.BACE$`Inverted Image Mean`[qc] <= (quality.check1.BACE$median_Slice[qc]+(2*quality.check1.BACE$SD_Slice[qc]))){
    quality.check1.BACE$Quality_Test_A_passed[qc] <- "Yes"  
    quality.check1.BACE$`Overall Sample Quality`[qc] <- (quality.check1.BACE$`Overall Sample Quality`[qc]+50)}else
      quality.check1.BACE$Quality_Test_A_passed[qc] <- "No" 
    
    if(quality.check1.BACE$`Inverted Image Mean`[qc] >= (quality.check1.BACE$median_Animal[qc]-(2*quality.check1.BACE$SD_Animal[qc])) & quality.check1.BACE$`Inverted Image Mean`[qc] <= (quality.check1.BACE$median_Animal[qc]+(2*quality.check1.BACE$SD_Animal[qc]))){
      quality.check1.BACE$Quality_Test_B_passed[qc] <- "Yes"  
      quality.check1.BACE$`Overall Sample Quality`[qc] <- (quality.check1.BACE$`Overall Sample Quality`[qc]+50)}else
        quality.check1.BACE$Quality_Test_B_passed[qc] <- "No"
}

quality.check2.BACE <- quality.check1.BACE %>%
  group_by(`Animal ID`,`Channel`) %>%
  summarize(`Overall Quality`=mean(`Overall Sample Quality`)) %>%
  ungroup

write.csv2(quality.check1.BACE,file= paste("quality.check1.BACE.csv"))
write.csv2(quality.check2.BACE,file= paste("quality.check2.BACE.csv"))
Mean.Sample.Quality.BACE <- mean(quality.check1.BACE$`Overall Sample Quality`)


plot.quality1.BACE <- ggplot(quality.check1.BACE,aes(x=`Animal ID`,y=`Inverted Image Mean`, color=`Slice No.`, fill=`Channel`)) +
  geom_boxplot(notch=FALSE, width=0.2,size=1,shape="x",color="black",outlier.shape = NA)  +
  geom_jitter(width=0.2,alpha=0.35) + 
  geom_hline(yintercept = median(quality.check1.BACE$`Inverted Image Mean`)-1.96*sd(quality.check1.BACE$`Inverted Image Mean`), linetype="F1", color = "black")+
  geom_hline(yintercept = median(quality.check1.BACE$`Inverted Image Mean`)+1.96*sd(quality.check1.BACE$`Inverted Image Mean`), linetype="F1", color = "black")+
  scale_color_npg() +
  theme_bw()+ 
  ylab("mean image intensity (BACE)") +
  xlab("Animal ID")+
  facet_grid(rows=vars(`Channel`))+
  ylim(0,255) +
  theme(axis.text.x = element_text(angle = 90))+
  scale_size(range=c(1,3.5))

plot.quality2.BACE <- ggplot(quality.check1.BACE,aes(x=`Slice No.`,y=`Inverted Image Mean`, color=`Slice No.`, fill=`Channel`)) +
  geom_jitter(width=0.2,alpha=0.5) + 
  geom_boxplot(notch=FALSE, width=0.2,size=1,shape="x",color="black",outlier.shape = NA)  +
  geom_hline(yintercept = median(quality.check1.BACE$`Inverted Image Mean`)-1.96*sd(quality.check1.BACE$`Inverted Image Mean`), linetype="F1", color = "black")+
  geom_hline(yintercept = median(quality.check1.BACE$`Inverted Image Mean`)+1.96*sd(quality.check1.BACE$`Inverted Image Mean`), linetype="F1", color = "black")+
  scale_color_npg() +
  theme_light()+ 
  ylab("mean image intensity (BACE)") +
  xlab("Slice No.")+
  facet_grid(cols=vars(`Channel`))+
  ylim(0,255)+
  scale_size(range=c(1,3.5))

png(paste(Sys.Date(),"plot.quality1.BACE.png",sep="_"),width = 10, height = 7, units = 'in', res = 300)
print(plot.quality1.BACE)     
dev.off() 

png(paste(Sys.Date(),"plot.quality2.BACE.png",sep="_"),width = 10, height = 7, units = 'in', res = 300)
print(plot.quality2.BACE)     
dev.off() 


select.data.all <- subset(integral.data.4G8,key==paste("X",as.character(analysis.threshold.4G8),sep=""))

select.data.BACE <- subset(integral.data.BACE,key==paste("X",as.character(analysis.threshold.BACE),sep=""))


select.data.all <- select.data.all %>%
  arrange(`Animal ID`,`Channel`)
select.data.BACE <- select.data.BACE %>%
  arrange(`Animal ID`,`Channel`)

select.data.all <- cbind(select.data.all[,1],`File Name (II)`="",select.data.all[,-1],`BACE_Av_integral_hist`=NA,`Ratio 4G8/BACE`=NA,`Ratio BACE/4G8`=NA,`Quality Overall`="O.K.")
select.data.all$`File Name (II)` <- as.character(select.data.all$`File Name (II)`)
select.data.all$`Quality Overall` <- as.character(select.data.all$`Quality Overall`)


for(n in 1:NROW(select.data.all)){
  for(nn in 1:NROW(select.data.BACE)){
    if(select.data.all$`Animal ID`[n]==select.data.BACE$`Animal ID`[nn] & select.data.all$`Slice No.`[n]==select.data.BACE$`Slice No.`[nn]&select.data.all$`Image No.`[n]==select.data.BACE$`Image No.`[nn]& select.data.all$`Experimental Group`[n]==select.data.BACE$`Experimental Group`[nn] ){
      select.data.all$`Ratio 4G8/BACE`[n] <- select.data.all$integral_hist[n]/select.data.BACE$integral_hist[nn]
      select.data.all$`BACE_Av_integral_hist`[n] <- select.data.BACE$integral_hist[nn]
      select.data.all$`Ratio BACE/4G8`[n] <- select.data.BACE$integral_hist[nn]/select.data.all$integral_hist[n]
      select.data.all$`File Name (II)`[n] <- as.character(select.data.BACE$`File Name`[nn])
        }
  }
}

for(n in 1:NROW(select.data.all)){
  for(nn in 1:NROW(select.data.BACE)){
    if(select.data.all$`Animal ID`[n]==select.data.BACE$`Animal ID`[nn] & select.data.all$`Slice No.`[n]==select.data.BACE$`Slice No.`[nn]&select.data.all$`Image No.`[n]==select.data.BACE$`Image No.`[nn]& select.data.all$`Experimental Group`[n]==select.data.BACE$`Experimental Group`[nn] ){
      select.data.all$`Ratio 4G8/BACE`[n] <- select.data.all$integral_hist[n]/select.data.BACE$integral_hist[nn]
      select.data.all$`BACE_Av_integral_hist`[n] <- select.data.BACE$integral_hist[nn]
      select.data.all$`Ratio BACE/4G8`[n] <- select.data.BACE$integral_hist[nn]/select.data.all$integral_hist[n]
      select.data.all$`File Name (II)`[n] <- as.character(select.data.BACE$`File Name`[nn])
    }
  }
}

colnames(select.data.all)[1] <- "File Name (I)"
select.data.all <- select.data.all[,-6]
colnames(select.data.all)[12] <- "Analysis Bin/Threshold"
colnames(select.data.all)[13] <- "4G8_Av_integral_hist"


for(s in 1:NROW(quality.check1.4G8)){
  if(quality.check1.4G8$Quality_Test_B_passed[s]=="No"){
    for(ss in 1:NROW(select.data.all)){
      if(select.data.all$`File Name (I)`[ss] == as.character(quality.check1.4G8$`File Name`[s])){
        select.data.all$`Quality Overall`[ss] <- "Quality Flagged"
      }
    }
  }

  if(quality.check1.4G8$Quality_Test_A_passed[s]=="No"){
    for(ss in 1:NROW(select.data.all)){
      if(select.data.all$`File Name (I)`[ss] == as.character(quality.check1.4G8$`File Name`[s])){
        select.data.all$`Quality Overall`[ss] <- "Quality Flagged"
      }
    }
  }
}

for(s in 1:NROW(quality.check1.BACE)){
  if(quality.check1.BACE$Quality_Test_B_passed[s]=="No"){
    for(ss in 1:NROW(select.data.all)){
      if(select.data.all$`File Name (II)`[ss] == as.character(quality.check1.BACE$`File Name`[s])){
        select.data.all$`Quality Overall`[ss] <- "Quality Flagged"
      }
    }
  }
  if(quality.check1.BACE$Quality_Test_A_passed[s]=="No"){
    for(ss in 1:NROW(select.data.BACE)){
      if(select.data.all$`File Name (II)`[ss] == as.character(quality.check1.BACE$`File Name`[s])){
        select.data.all$`Quality Overall`[ss] <- "Quality Flagged"
      }
    }
  }
}

write.csv2(select.data.all,"select.data.all.csv")
#write.csv2(select.data.BACE,"select.data.BACE.csv")


select.data.all <- subset(select.data.all,!is.na(`BACE_Av_integral_hist`))
select.data.all <- subset(select.data.all,!is.na(`4G8_Av_integral_hist`))

#select.data.all <- subset(select.data.all,`Quality Overall`!="Quality Flagged")


select.data.all <- select.data.all %>%
  arrange(`Ratio BACE/4G8`)

num.obs.4G8 <- count(select.data.all, vars=`Animal ID`,`Experimental Group`)
colnames(num.obs.4G8) <- c("Animal ID","Experimental Group", "No.Obs.Units")


select.data.all <- aggregate(select.data.all[,13:16], by=list(`Animal ID`=select.data.all$`Animal ID`,`Experimental Group`=select.data.all$`Experimental Group`,`Genotype APP23`=select.data.all$`APP23 Genotype`,`Genotype p40`=select.data.all$`p40 Genotype`,`Age in Days`=select.data.all$`Age in Days`,`Gender`=select.data.all$`Gender`),FUN=mean)

select.data.all <- cbind(select.data.all,"No. Obs. Units"=0)

for(o in 1:length(select.data.all$`Animal ID`)){
  for(oo in 1: length(num.obs.4G8$`Animal ID`)){
    if(select.data.all$`Animal ID`[o]==num.obs.4G8$`Animal ID`[oo] & select.data.all$`Experimental Group`[o]==num.obs.4G8$`Experimental Group`[oo]){
      select.data.all$`No. Obs. Units`[o] <- num.obs.4G8$No.Obs.Units[oo]
    }
  }
}

select.data.all <- subset(select.data.all,`Experimental Group` != "Secondary AB Control")

ymax.ratio <- max(select.data.all$`Ratio BACE/4G8`)*1.3

plot.4G8.ratio <- ggplot(select.data.all,aes(size=`No. Obs. Units` ,x=`Experimental Group`,y=`Ratio BACE/4G8`, color=`Experimental Group`, fill=`Experimental Group`, shape=`Gender`)) +
  geom_boxplot(color="black",width=0.2,size=1,shape="x", outlier.shape = NA,alpha=0.5) +
  geom_jitter(width=0.2,alpha=0.5,color="black") + 
  scale_color_npg() +
  theme_classic()+ 
  ylab("Ratio BACE/4G8") +
  xlab("put in p value")+
  scale_size(range=c(1,3.5))+ 
  facet_grid(cols=vars(`Gender`))+
  ylim(0,ymax.ratio)+
  theme(axis.text.x = element_blank())+  
  scale_shape_manual(values=c(21,22,23))
  

plot.bubble.ratio <- ggplot(select.data.all,aes(size=`No. Obs. Units`,x=`4G8_Av_integral_hist`,y=`BACE_Av_integral_hist`, group=interaction(`Experimental Group`,`Animal ID`),color=`Experimental Group`, fill=`Experimental Group`, shape=`Gender`)) +
  geom_point()+
  scale_color_npg() +
  theme_classic()+ 
  ylab("BACE Staining Signal") +
  xlab("4G8 Plaque Load")+
  scale_size(range=c(1,3)) + 
  scale_shape_manual(values=c(21,22,23))+
  geom_text(aes(label=`Animal ID`),hjust=1.25, vjust=0,color="black",size=3)


write.csv2(select.data.all, "select.data.all.means.csv")

png(paste(Sys.Date(),"plot.4G8.ratio.png",sep="_"),width = 10, height = 7, units = 'in', res = 300)
print(plot.4G8.ratio)     
dev.off() 
png(paste(Sys.Date(),"plot.bubble.ratio.png",sep="_"),width = 10, height = 7, units = 'in', res = 300)
print(plot.bubble.ratio)     
dev.off() 

integral.data.4G8 <- aggregate(integral.data.4G8[,13],by=list(`Animal ID`=integral.data.4G8$`Animal ID`,`Experimental Group`=integral.data.4G8$`Experimental Group`,`Genotype APP23`=integral.data.4G8$`APP23 Genotype`,`Genotype p40`=integral.data.4G8$`p40 Genotype`,`Age in Days`=integral.data.4G8$`Age in Days`,`Gender`=integral.data.4G8$`Gender`,`key`=integral.data.4G8$key),FUN=median)

integral.data.BACE <- aggregate(integral.data.BACE[,13],by=list(`Animal ID`=integral.data.BACE$`Animal ID`,`Experimental Group`=integral.data.BACE$`Experimental Group`,`Genotype APP23`=integral.data.BACE$`APP23 Genotype`,`Genotype p40`=integral.data.BACE$`p40 Genotype`,`Age in Days`=integral.data.BACE$`Age in Days`,`Gender`=integral.data.BACE$`Gender`,`key`=integral.data.BACE$key),FUN=mean)

plot.integral.4G8 <- ggplot(data=integral.data.4G8,aes(x=key, y=`integral_hist`, group=`Experimental Group`)) +
  #geom_line()+ 
  geom_point(aes(color=`Experimental Group`, shape=`Gender`)) + 
  geom_vline(xintercept = analysis.threshold.4G8, linetype="dotted", color = "black")+
  scale_color_npg() +
  theme_bw()+
  ylab("Cumulative Signal [%]")+ 
  xlab("Thresholding: Signal vs. Background") +
  theme(axis.text.x = element_blank()) #+
#  facet_grid(rows = vars(`Experimental Group`))
  
plot.integral.BACE <- ggplot(data=integral.data.BACE,aes(x=key, y=`integral_hist`, group=interaction(`Experimental Group`,`Gender`))) +
  #geom_line()+ 
  geom_point(aes(color=`Experimental Group`, shape=`Gender`)) + 
  geom_vline(xintercept = analysis.threshold.4G8, linetype="dotted", color = "black")+
  scale_color_npg() +
  theme_bw()+
  ylab("Cumulative Signal [%]")+ 
  xlab("Thresholding: Signal vs. Background") +
  theme(axis.text.x = element_blank())
  #facet_grid(rows = vars(`Experimental Group`))


# if you need example images being thresholded, here´s code for a random selection from all your images. 
#You can plot the thresholded images with plot(threshold(ex.i1,thr=a(nalysis.threshold/hist.bins))):

#ex.i <- NULL
#for(ei in 1:round(length(file.names)/10,digits=0)){
#  ex.i <- c(ex.i,round((ei*length(file.names)/10),digits=0))
#}
#
#ex.i1 <- grayscale(load.image(file.names[ex.i[1]]))
#ex.i2 <- grayscale(load.image(file.names[ex.i[2]]))
#ex.i3 <- grayscale(load.image(file.names[ex.i[3]]))
#ex.i4 <- grayscale(load.image(file.names[ex.i[4]]))
#ex.i5 <- grayscale(load.image(file.names[ex.i[5]]))
#ex.i6 <- grayscale(load.image(file.names[ex.i[6]]))
#ex.i7 <- grayscale(load.image(file.names[ex.i[7]]))
#ex.i8 <- grayscale(load.image(file.names[ex.i[8]]))
#ex.i9 <- grayscale(load.image(file.names[ex.i[9]]))
#ex.i10 <- grayscale(load.image(file.names[ex.i[10]]))




