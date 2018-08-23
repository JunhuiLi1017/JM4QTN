calGenoProb <- function(GeneticMap,GenoData,method,croType=NULL,steps=0,Gn=2){
    if(nrow(GeneticMap) != ncol(GenoData)){
      stop("Marker No. in geneticMap and GenoData must be equal\n")
    }
    if(method=="AM"){
      GenoData[is.na(GenoData)] <- 9
      tGeneticMapt <- t(GeneticMap[,-1])
      colnames(tGeneticMapt) <- colnames(GenoData)
      RetValue <- rbind(tGeneticMapt,GenoData)
    }else if(method=="LM"){
      if(steps < 0){
        stop("Parameter steps must be postive\n")
      }
      if(croType=="RIL" | croType=="DH") {
        if(!all(is.na(GenoData[!GenoData==2 & !GenoData==0] ))){
          warning("GenoData numeric should be 0 or 2 for RIL and DH population, other numeric has been transformed to NA\n")
          GenoData[!GenoData==2 & !GenoData==0] <- NA
        }
      }
      if(croType=="F2") {
        if(!all(is.na(GenoData[!GenoData==2 & !GenoData==1 & !GenoData==0] ))){
          warning("GenoData numeric should be 0 1 2 for F2 population, other numeric has been transformed to NA\n")
          GenoData[!GenoData==2 & !GenoData==1 & !GenoData==0] <- NA
        }
      }
      if(croType=="BCP1" | croType=="BCP2") {
        if(!all(is.na(GenoData[!GenoData==2 & !GenoData==1 & !GenoData==0] ))){
          warning("GenoData numeric should be 1 2 and 1 0 for BCP1 and BCP2 population respectively, other numeric has been transformed to NA\n")
          GenoData[!GenoData==2 & !GenoData==1 & !GenoData==0] <- NA
        }
      }
      NoChr <- nlevels(as.factor(GeneticMap$chr))
      ConChr <- levels(as.factor(GeneticMap$chr))
      rownames(GeneticMap) <- GeneticMap[,1]
      #RMapList
      RMapList <- NULL
      for (i in 1:NoChr){
        tempV <- subset(GeneticMap,GeneticMap$chr==ConChr[i])
        RMapList[[i]] <- tempV$pos
        names(RMapList[[i]]) <- tempV$marker
      }
      if(steps==0){
        VMmatrix <- GenoData
        nObs <- nrow(GenoData)
        for(i in 1:NoChr){
          Chr_R <- RMapList[[i]]
          for(nS in 1:nObs){
            for(j in 1:length(Chr_R)){
              if(GenoData[nS,names(Chr_R)[j]] %in% NA){
                jL <- j
                jR <- j
                while(jR < length(Chr_R) & GenoData[nS,names(Chr_R)[jR]] %in% NA){
                  jR <- jR+1
                }
                while(jL > 1 & GenoData[nS,names(Chr_R)[jL]] %in% NA){
                  jL <- jL-1
                }
                if(jL==1 & GenoData[nS,names(Chr_R)[jL]] %in% NA){
                  xb <- as.numeric(Chr_R[jR]-Chr_R[j])
                  rb <- HaldaneMap(xb/100)
                  if(croType != "BCP2" & GenoData[nS,names(Chr_R)[jR]] %in% 2){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("N2",croType,Gn,rb),digits=3)+1
                  }
                  if(croType != "DH" & croType != "RIL"& GenoData[nS,names(Chr_R)[jR]] %in% 1){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("N1",croType,Gn,rb),digits=3)+1
                  }
                  if(croType != "BCP1"& GenoData[nS,names(Chr_R)[jR]] %in% 0){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("N0",croType,Gn,rb),digits=3)+1
                  }
                }
                if(jR==length(Chr_R) & GenoData[nS,names(Chr_R)[jR]] %in% NA){
                  xa <- as.numeric(Chr_R[j]-Chr_R[jL])
                  ra <- HaldaneMap(xa/100)
                  if (croType != "BCP2" & GenoData[nS,names(Chr_R)[jL]] %in% 2){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("2N",croType,Gn,ra),digits=3)+1
                  }
                  if (croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 1){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("1N",croType,Gn,ra),digits=3)+1
                  }
                  if (croType != "BCP1" & GenoData[nS,names(Chr_R)[jL]] %in% 0){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("0N",croType,Gn,ra),digits=3)+1
                  }
                }
                #No 1st and last missing marker
                if(!GenoData[nS,names(Chr_R)[jL]] %in% NA & !GenoData[nS,names(Chr_R)[jR]] %in% NA){
                  #else{
                  xa <- as.numeric(Chr_R[j]-Chr_R[jL])
                  xb <- as.numeric(Chr_R[jR]-Chr_R[j])
                  ra <- HaldaneMap(xa/100)
                  rb <- HaldaneMap(xb/100)
                  if (croType != "BCP2" & GenoData[nS,names(Chr_R)[jL]] %in% 2 & GenoData[nS,names(Chr_R)[jR]] %in% 2){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("22",croType,Gn,ra,rb),digits=3)+1
                  }
                  if (croType != "BCP2" & croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 2 & GenoData[nS,names(Chr_R)[jR]] %in% 1){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("21",croType,Gn,ra,rb),digits=3)+1
                  }
                  if (croType != "BCP2" & croType != "BCP1" & GenoData[nS,names(Chr_R)[jL]] %in% 2 & GenoData[nS,names(Chr_R)[jR]] %in% 0){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("20",croType,Gn,ra,rb),digits=3)+1
                  }
                  if (croType != "BCP2" & croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 1 & GenoData[nS,names(Chr_R)[jR]] %in% 2){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("12",croType,Gn,ra,rb),digits=3)+1
                  }
                  if (croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 1 & GenoData[nS,names(Chr_R)[jR]] %in% 1){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("11",croType,Gn,ra,rb),digits=3)+1
                  }
                  if (croType != "BCP1" & croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 1 & GenoData[nS,names(Chr_R)[jR]] %in% 0){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("10",croType,Gn,ra,rb),digits=3)+1
                  }
                  if (croType != "BCP2" & croType != "BCP1" & GenoData[nS,names(Chr_R)[jL]] %in% 0 & GenoData[nS,names(Chr_R)[jR]] %in% 2){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("02",croType,Gn,ra,rb),digits=3)+1
                  }
                  if (croType != "BCP1" & croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 0 & GenoData[nS,names(Chr_R)[jR]] %in% 1){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("01",croType,Gn,ra,rb),digits=3)+1
                  }
                  if (croType != "BCP1" & GenoData[nS,names(Chr_R)[jL]] %in% 0 & GenoData[nS,names(Chr_R)[jR]] %in% 0){
                    VMmatrix[nS,names(Chr_R)[j]] <- round(Eagd("00",croType,Gn,ra,rb),digits=3)+1
                  }
                } 
              }
            }#j
          }#nObs
        }#NoChr
        VMmatrix <- round(VMmatrix,digits = 3)
        RetValue <- rbind(t(GeneticMap[,-1]),VMmatrix)
      }else{
        #VMapList
        VMapList <- NULL
        geneticMap <- NULL
        for(i in 1:NoChr){
          #i=1
          chi <- RMapList[[i]]
          newCh <- chi[1]
          n <- 0
          for(j in 2:length(chi)){
            #j=8
            gap <- chi[j]-chi[j-1]
            cushion <- chi[j-1]
            if(gap > steps){
              nVM <- trunc(gap/steps)
              perdist <- gap/(nVM+1)
              vecVM <- round(cushion + perdist*c(1:nVM),digits=3)
              names(vecVM) <- paste("V",i,"_",c((n+j):(n+j+nVM-1)),sep="") 
              newCh <- append(newCh,vecVM)
              n <- n+nVM
            }#gap
            newCh <- append(newCh,chi[j])
          }#j
          subMap <- matrix(NA,2,length(newCh))
          subMap <- as.data.frame(subMap)
          colnames(subMap) <- names(newCh)
          #subMap <- as.data.frame(subMap)
          rownames(subMap) <- c("chr","pos")
          subMap[1,] <- ConChr[i]
          subMap[2,] <- newCh
          if(i==1){
            geneticMap <- subMap
          }else{
            geneticMap <- cbind(geneticMap,subMap)
          }
          VMapList[i] <- list(newCh)
        }#i
        #VMapLine
        NoMM <- 0      #No. of Marker on the Map
        LineDis <- 0   #Line distance
        VMapLine <- NULL     #mapline including Virtual markers
        FMCh <- NULL   #vector of First marker of each chromosome 
        EMCh <- NULL   #vector of End marker of each chromosome 
        for (i in 1:NoChr) {
          chi <- VMapList[[i]]
          StartP <- chi[1]
          for (j in 1:length(chi)) {
            NoMM <- NoMM + 1
            VMapLine <- append(VMapLine, (chi[j]+LineDis))
            names(VMapLine)[NoMM] <- names(chi)[j]
            if (j == 1) {
              FMCh <- append(FMCh,names(chi)[j])
            }else if (j == length(chi)) {
              EMCh <- append(EMCh,names(chi)[j])
              LineDis <- VMapLine[NoMM]
            }
          }
          #LineDis <- LineDis + 20
        }
        nM_V <- NoMM          #nM_V: No. of Markers
        nObs <- nrow(GenoData)
        VMmatrix <- matrix(NA,nObs,nM_V)
        colnames(VMmatrix) <- names(VMapLine)
        rownames(VMmatrix) <- rownames(GenoData)
        for(j in 1:ncol(GenoData)){
          VMmatrix[,colnames(GenoData)[j]] <- GenoData[,j]
        }
        for(i in 1:NoChr){
          #i=1
          Chr_R <- RMapList[[i]]
          Chr_V <- VMapList[[i]]
          for(jr in 2:length(Chr_R)){
            #jr=8
            nL <- which(names(Chr_V) %in% names(Chr_R)[jr-1])
            nR <- which(names(Chr_V) %in% names(Chr_R)[jr])
            for(nS in 1:nObs){
              #nS=2
              jR <- jr  #j=1
              jL <- jr-1
              while(jR < length(Chr_R) & GenoData[nS,names(Chr_R)[jR]] %in% NA){
                jR <- jR+1
              }
              while(jL > 1 & GenoData[nS,names(Chr_R)[jL]] %in% NA){
                jL <- jL-1
              }
              for(jv in nL:nR){
                #jv = 66
                if(VMmatrix[nS,names(Chr_V)[jv]] %in% NA){
                  if(jL==1 & jR ==length(Chr_R)){
                    stop("marker data without valid marker")
                  }
                  #the 1st missing marker 
                  if(jL==1 & GenoData[nS,names(Chr_R)[jL]] %in% NA){
                    xb <- as.numeric(Chr_R[jR]-Chr_V[jv])
                    rb <- HaldaneMap(xb/100)
                    if(croType != "BCP2" & GenoData[nS,names(Chr_R)[jR]] %in% 2){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("N2",croType,Gn,rb),digits=3)+1
                    }
                    if(croType != "DH" & croType != "RIL"& GenoData[nS,names(Chr_R)[jR]] %in% 1){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("N1",croType,Gn,rb),digits=3)+1
                    }
                    if(croType != "BCP1"& GenoData[nS,names(Chr_R)[jR]] %in% 0){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("N0",croType,Gn,rb),digits=3)+1
                    }
                  }
                  #the last missing marker
                  if(jR==length(Chr_R) & GenoData[nS,names(Chr_R)[jR]] %in% NA){
                    xa <- as.numeric(Chr_V[jv]-Chr_R[jL])
                    ra <- HaldaneMap(xa/100)
                    if (croType != "BCP2" & GenoData[nS,names(Chr_R)[jL]] %in% 2){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("2N",croType,Gn,ra),digits=3)+1
                    }
                    if (croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 1){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("1N",croType,Gn,ra),digits=3)+1
                    }
                    if (croType != "BCP1" & GenoData[nS,names(Chr_R)[jL]] %in% 0){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("0N",croType,Gn,ra),digits=3)+1
                    }
                  }
                  #No 1st and last missing marker
                  if(!GenoData[nS,names(Chr_R)[jL]] %in% NA & !GenoData[nS,names(Chr_R)[jR]] %in% NA){
                    #else{
                    xa <- as.numeric(Chr_V[jv]-Chr_R[jL])
                    xb <- as.numeric(Chr_R[jR]-Chr_V[jv])
                    ra <- HaldaneMap(xa/100)
                    rb <- HaldaneMap(xb/100)
                    if (croType != "BCP2" & GenoData[nS,names(Chr_R)[jL]] %in% 2 & GenoData[nS,names(Chr_R)[jR]] %in% 2){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("22",croType,Gn,ra,rb),digits=3)+1
                    }
                    if (croType != "BCP2" & croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 2 & GenoData[nS,names(Chr_R)[jR]] %in% 1){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("21",croType,Gn,ra,rb),digits=3)+1
                    }
                    if (croType != "BCP2" & croType != "BCP1" & GenoData[nS,names(Chr_R)[jL]] %in% 2 & GenoData[nS,names(Chr_R)[jR]] %in% 0){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("20",croType,Gn,ra,rb),digits=3)+1
                    }
                    if (croType != "BCP2" & croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 1 & GenoData[nS,names(Chr_R)[jR]] %in% 2){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("12",croType,Gn,ra,rb),digits=3)+1
                    }
                    if (croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 1 & GenoData[nS,names(Chr_R)[jR]] %in% 1){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("11",croType,Gn,ra,rb),digits=3)+1
                    }
                    if (croType != "BCP1" & croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 1 & GenoData[nS,names(Chr_R)[jR]] %in% 0){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("10",croType,Gn,ra,rb),digits=3)+1
                    }
                    if (croType != "BCP2" & croType != "BCP1" & GenoData[nS,names(Chr_R)[jL]] %in% 0 & GenoData[nS,names(Chr_R)[jR]] %in% 2){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("02",croType,Gn,ra,rb),digits=3)+1
                    }
                    if (croType != "BCP1" & croType != "DH" & croType != "RIL" & GenoData[nS,names(Chr_R)[jL]] %in% 0 & GenoData[nS,names(Chr_R)[jR]] %in% 1){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("01",croType,Gn,ra,rb),digits=3)+1
                    }
                    if (croType != "BCP1" & GenoData[nS,names(Chr_R)[jL]] %in% 0 & GenoData[nS,names(Chr_R)[jR]] %in% 0){
                      VMmatrix[nS,names(Chr_V)[jv]] <- round(Eagd("00",croType,Gn,ra,rb),digits=3)+1
                    }
                  }
                }
              }#jv
            }#nObs
          }#jr
        }#NoChr
        VMmatrix <- round(VMmatrix,digits = 3)
        RetValue <- rbind(geneticMap,VMmatrix)
      }#steps
    }else{
      stop("method must be AM or LM")
    }
    
    #* --------------------------
    # create and set working dir
    #*---------------------------
    mainDir=getwd()
    subDir="OUTPUT_GenoProb"
    ifelse(!file.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
    if(method=="AM"){
      write.table(RetValue,file=paste(file.path(mainDir, subDir),"/GenoData_",method,".xls",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")
    }else if(method=="LM"){
      write.table(RetValue,file=paste(file.path(mainDir, subDir),"/GenoData_step",steps,"_",method,".xls",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")
    }  
    return(as.data.frame(RetValue))
  }
