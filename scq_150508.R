#2015-11-16 Single cell qPCR y1 + Actb. 

#Load data:
#load("~/Desktop/PhD/Results/2015/Single cell PCR/2015-05-08/2015-05-08 Single cell data.RData")

View(`2015-05-08 Single cell qPCR`)

#Process data:
dataA <- `2015-05-08 Single cell qPCR`

#Replace those values equal to 999 or 50 (non-amplified) for 40:
dataA[,2:ncol(dataA)] <- apply(dataA[,2:ncol(dataA)],2,function(x) ifelse(x==50,40,x))

#Subset by timepoint:
t_data_day3.0 <- t(dataA[,grep("Day7|no.cell|bulk|gene",colnames(dataA),invert = T)])
t_data_day7.0 <- t(dataA[,grep("Day3|no.cell|bulk|gene",colnames(dataA),invert = T)])
colnames(t_data_day3.0)<-dataA$gene
colnames(t_data_day7.0)<-dataA$gene


#Boxplots: 
new.mfrow <- par(mfrow=c(2,4))
       for (i in 1:22){
             zz1 <- list(t_data_day3.0[,i],t_data_day7.0[,i])
           if(i==20) new.lim<-c(12,30)
           else new.lim<-c(20,40) 
 
 boxplot(zz1,    
         outline=F,border="grey 61", col="grey 91", boxwex=0.5,
         ylab="Cq value",xlab="Subset",
         cex.lab=1.2,
         ylim=new.lim,
         names=c("3.0d","7.0d"))
       
 stripchart(zz1[1],
            vertical=T,method="jitter",add=T,
            pch=21,cex=1.2,col="black",bg="black")
 
 stripchart(zz1[2],at=2,
            vertical=T,method="jitter",add=T,
            pch=21,cex=1.2,col="black",bg="blue")
 
 
 plotT<-dataA$gene[i]
 title(main=plotT,cex.main=1.7)
}

par(new.mfrow)
par(pin=c(5,5))

#Save as PDF: size 8 x 15.5 inches


#Subset germline transcripts + cMyc + Apex1 + Actb:
#Targets: c(1,3,5,8,11,14,17,20)
dataB <- dataA[c(1,3,5,8,11,14,17,15,9,20),]
dataC <- dataB[c(1,5,6,2,3,7,4,8,9,10),]

#Subset data from day 3.0:
 data_3.0dA <- dataC[c(2:59)]
 data_3.0d  <- as.matrix(data_3.0dA[-c(7,15,19,21,24,26,37,38,40,41,48,49,50,55)])

 rownames(data_3.0d) <- dataC[,1]
 colnames(data_3.0d) <- c(paste("c0",1:9,sep=""),paste("c",10:42,sep=""),"NTC","Bulk")

tdata_3.0d<-t(data_3.0d)

#Subset data from day 7.0:
data_7.0d<-as.matrix(dataC[60:103])
rownames(data_7.0d) <- dataC[,1]
colnames(data_7.0d) <- c(paste("c0",1:9,sep=""),paste("c",10:42,sep=""),"NTC","Bulk")

tdata_7.0d<-t(data_7.0d)


x <- as.matrix(cbind(tdata_3.0d,tdata_7.0d))
y <- as.matrix(rbind(data_3.0d,data_7.0d))



#Activate gplots package:
library(gplots)


#Set up heatmap-colors:
   colorRB1 <- colorRampPalette(c("gray51","yellow","darkviolet","black"))

target_lab <- c(expression(gamma[1]),expression(gamma[3]),expression(gamma["2b"]),
                expression(gamma["2c"]),expression(epsilon),expression(mu),
                expression(mu - gamma[1]),"c-Myc","Apex1","Actb")

#Heatmap:
heatmap.2(y,Rowv=F,Colv=F,dendrogram ="none",scale="none",
          col=colorRB1,
          trace="none",colsep=(1:(ncol(y))),rowsep=(1:(nrow(y))),sepcolor="black",
          srtCol=60,cexCol=1,cexRow=1,sepwidth=c(0.07,0.07),
          adjRow = c(0,0.5),
          #labRow=c(target_lab,target_lab),
          breaks=seq(0,40,0.5),
          lmat=rbind(c(0,3,3,2,0),
                     c(0,1,1,0,0),
                     c(0,1,1,0,0),
                     c(4,4,4,0,0)),
          lwid=c(0.5,1,1,1,0.05), lhei=c(0.25,0.7,0.7,0.5),
          density.info="none",keysize=0.5)













######################
target_lab <- c(expression(gamma[1]),expression(gamma[3]),expression(gamma["2b"]),
                 expression(gamma["2c"]),expression(epsilon),expression(mu),
                 expression(mu - gamma[1]),"Actb")

#Heatmap:
heatmap.2(x,Rowv=F,Colv=F,dendrogram ="none",scale="none",
          col=colorRB1,
          trace="none",colsep=(1:(ncol(x))),rowsep=(1:(nrow(x))),sepcolor="black",
          srtCol=60,cexCol=0.8,cexRow=0.8,sepwidth=c(0.07,0.07),
          adjRow = c(0,0.5),
          #labRow=c(rep("1 cell",50),data[51:77,2]),
          #labCol=c(target_lab,target_lab),
          breaks=seq(0,40,0.5),
          lmat=rbind(c(0,3,3,2,0),
                     c(0,1,1,0,0),
                     c(0,1,1,0,0),
                     c(4,4,4,0,0)),
          lwid=c(0.5,0.8,0.8,3,0.05), lhei=c(0.25,0.7,0.7,0.5),
          density.info="none",keysize=0.5)

title("Single cell qPCR", cex.main = 2)

#Save 5000x35000
