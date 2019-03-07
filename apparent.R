apparent = function (InputFile, MaxIdent=0.10, alpha=0.01, nloci=300, self=TRUE, plot=TRUE, files=TRUE, Dyad=FALSE) {
  
  ################################################################################################
  # Parse the tab-delimited input file and convert genotypic states to numeric genotypic classes, 
  # based on primary and secondary alleles across the population.
  ################################################################################################

  # Some previous definitions
  pos = 1
  setwd(getwd())
  envir = as.environment(pos)
  options(warn=-1)
  ptm <- proc.time()
  
  print("Importing data...")
  GK <- cbind(as.data.frame(InputFile[,1]),as.data.frame(InputFile[,2]))
  colnames(GK) <- c("genos","key")
  Data <- t(as.data.frame(InputFile[,3:ncol(InputFile)]))
  ConvertedData <- data.frame (matrix (ncol = ncol(Data), nrow = 0))
  
  for (i in 1:nrow(Data)) {    
    alleles <- setdiff(strsplit(paste(Data[i,],collapse=""),"")[[1]],c("/","-"))
    pri <- paste(alleles[1],"/",alleles[1],sep="")
    het1 <- paste(alleles[1],"/",alleles[2],sep="")
    het2 <- paste(alleles[2],"/",alleles[1],sep="")
    sec <- paste(alleles[2],"/",alleles[2],sep="")
    
    Data[i,][Data[i,] == pri] <- 0
    Data[i,][Data[i,] == het1 | Data[i,] == het2] <- 0.5
    Data[i,][Data[i,] == sec] <- 1
    Data[i,][Data[i,] == "-/-"] <- NA
    
    line <- as.numeric(as.vector(Data[i,]))
    ConvertedData <- rbind(ConvertedData,line)
  }
  colnames(ConvertedData) <- as.vector(t(GK$genos))
  
  #################################################################################################
  # Create the working matrix for the population, based on the individual key assignments (Mothers, 
  # Fathers, Parents, Offspring, or All) provided in the input file.
  #################################################################################################
  
  Mothers <- vector(mode="numeric",length=0)
  MothersNames <- list()
  Fathers <- vector(mode="numeric",length=0)
  FathersNames <- list()
  Offs <- vector(mode="numeric",length=0)
  OffsNames <- list()
  
  for (i in 1:ncol(ConvertedData)){
    if (GK$key[i] == "Mo") {
      Mothers[i] <- as.data.frame(ConvertedData[,i])
      MothersNames <- append(MothersNames,as.name(as.matrix(GK[i,1])))
      next    
    } else if (GK$key[i] == "Fa") {
      Fathers[i] <- as.data.frame(ConvertedData[,i])
      FathersNames <- append(FathersNames,as.name(as.matrix(GK[i,1])))
      next  
    } else if (GK$key[i] == "Off") {
      Offs[i] <- as.data.frame(ConvertedData[,i])
      OffsNames <- append(OffsNames,as.name(as.matrix(GK[i,1])))
      next
    } else if (GK$key[i] == "Pa") {
      Mothers[i] <- as.data.frame(ConvertedData[,i])
      MothersNames <- append(MothersNames,as.name(as.matrix(GK[i,1])))
      Fathers[i] <- as.data.frame(ConvertedData[,i])
      FathersNames <- append(FathersNames,as.name(as.matrix(GK[i,1])))
      next
    } else if (GK$key[i] == "All") {
      Mothers[i] <- as.data.frame(ConvertedData[,i])
      MothersNames <- append(MothersNames,as.name(as.matrix(GK[i,1])))
      Fathers[i] <- as.data.frame(ConvertedData[,i])
      FathersNames <- append(FathersNames,as.name(as.matrix(GK[i,1])))
      Offs[i] <- as.data.frame(ConvertedData[,i])
      OffsNames <- append(OffsNames,as.name(as.matrix(GK[i,1])))
      next
    } else {
      stop("Please, check the format of the key column (Column 2) in your input file. Acceptable keys for each genotype are Mo, Fa, Off, Pa and All.")
    }
  }
  
  Mothers <- as.data.frame(Mothers[!sapply(Mothers,is.null)]) 
  colnames(Mothers) <- MothersNames 
  Fathers <- as.data.frame(Fathers[!sapply(Fathers,is.null)]) 
  colnames(Fathers) <- FathersNames 
  Offs <- as.data.frame(Offs[!sapply(Offs,is.null)]) 
  colnames(Offs) <- OffsNames
  
  WorkMatrix <- as.data.frame(cbind(Mothers,Fathers,Offs))
  
  MoN <- ncol(Mothers)
  MoS <- 1
  MoE <- ncol(Mothers)
  FaN <- ncol(Fathers)
  FaS <- MoE + 1
  FaE <- FaS + ncol(Fathers) - 1
  OfN <- ncol(Offs)  
  OfS <- FaE + 1
  OfE <- OfS + ncol(Offs) - 1    
  
  ##############################################################################################
  # Creating the genotypes of the Expected Progeny (EPij) for each pair of potential parents
  # (i and j), based only on parental homozygous loci. Then calculating the Gower Dissimilarity 
  # (GD) between each EPij and each potential Offspring (j) in the population (Offj).
  ###############################################################################################
  
  print("Creating the EP(ij)'s and estimating the GDij|POk...")
  
  # Load intermediate files, construct output vectors and initialize progress bar
  Parent1List <- vector(mode="numeric",length=0)
  Parent2List <- vector(mode="numeric",length=0)
  ObsProgList <- vector(mode="numeric",length=0)
  TypeList <- vector(mode="numeric",length=0)
  SNPsNumber <- vector(mode="numeric",length=0)
  GD <- vector(mode="numeric",length=0)
  CheckNames1 <- vector(mode="numeric",length=0)
  pb <- txtProgressBar(min=0,max=MoE,style=3)  
  
  # Loop through each pair of parents
  for (i in MoS:MoE) {
    for (j in FaS:FaE) {
      Sys.sleep(0.01)
      setTxtProgressBar(pb,i)
      
      # Skip parental pairs that were previously considered. 
      # Skip triads for which the same genotype is both parent and offspring (whenever "All" appears in the input file).  
      D1 <- paste(colnames(WorkMatrix[i]),colnames(WorkMatrix[j]),sep=".")
      D2 <- paste(colnames(WorkMatrix[j]),colnames(WorkMatrix[i]),sep=".")
      
      if (D1 %in% CheckNames1 && D2 %in% CheckNames1) {
        next        
      } else {
        CheckNames1 <- append(CheckNames1,c(D1,D2))
      }
      
      ParentsPair <- WorkMatrix[,c(i,j)]
      Offsprings <- WorkMatrix[,c(OfS:OfE)]
      
      OffsFinal <- data.frame(matrix(ncol=0,nrow=nrow(WorkMatrix)))
      FinalNames <- vector(mode="numeric",length=0)
      n <- 1
      for (k in 1:ncol(Offsprings)) {
        if ( identical(colnames(Offsprings[k]),colnames(ParentsPair[1])) | identical(colnames(Offsprings[k]),colnames(ParentsPair[2])) ) {
          next
        } else {
          OffsFinal[n] <- Offsprings[,k]
          FinalNames <- append(FinalNames,colnames(Offsprings[k]))
          n <- n + 1
        }
      }
      colnames(OffsFinal) <- FinalNames
      ParentsPairOffs <- cbind(ParentsPair,OffsFinal)
      OffsFinalN <- ncol(OffsFinal)
      
      # Defining the type of cross for each triad (self versus outcross)
      if ( colnames(WorkMatrix[i]) == colnames(WorkMatrix[j]) ) {
        Type <- "self"
        TypeList <- append(TypeList,rep(Type,OffsFinalN))
      } else {
        Type <- "outcross"
        TypeList <- append(TypeList,rep(Type,OffsFinalN))
      }
      
      # Creating the EP(ij)'s and calculating the GD(ij|k)'s...
      ExpProgName <- paste(c(colnames(ParentsPairOffs[1]),colnames(ParentsPairOffs[2])),sep=" vs. ",collapse=" vs. ")
      Parent1 <- colnames(ParentsPairOffs[1])
      Parent2 <- colnames(ParentsPairOffs[2])
      Parent1List <- append(Parent1List,rep(Parent1,OffsFinalN))
      Parent2List <- append(Parent2List,rep(Parent2,OffsFinalN))
      
      P1 <- ParentsPairOffs[,1]
      P2 <- ParentsPairOffs[,2]
      P1[P1 == 0.5] <- NA
      P2[P2 == 0.5] <- NA
      ExpProg <- as.data.frame((P1 + P2) / 2)
      colnames(ExpProg) <- ExpProgName
      
      for (l in 3:(ncol(ParentsPairOffs))) {
        Diff <- abs(ExpProg - ParentsPairOffs[,l])
        Diff <- Diff[!is.na(Diff)]
        SNPsNumber <- append (SNPsNumber,length(Diff))
        ObsProgList <- append(ObsProgList,colnames(ParentsPairOffs[l]))
        GD <- append(GD, (sum(Diff)) / length(Diff))
      }
    }  
  }
  close (pb)
  
  # Output 1 - All triads
  Out1 <- data.frame(Parent1List,Parent2List,ObsProgList,TypeList,SNPsNumber,GD)
  colnames(Out1) <- c("Parent1","Parent2","Offspring","Cross.Type","SNPs","GD")
  Triad_All <- Out1
  assign("Triad_All", Out1, envir = envir)

  if (files==TRUE){
    write.csv(Triad_All,"apparent-Triad-All.csv",row.names=F)
  }
  
  #########################################################################
  # Parsing and reporting the significatives of the Triad analysis.
  # Finding the triad GAP and testing its significance and calculates the 
  # p-values for the significant triads.
  #########################################################################
  print("Finding the Triad GAP, testing its significance and assign p-values for the significant Triad hits, if found....")

  # The log report files
  OMeanGD <- mean(Out1$GD)
  OSdGD <- sd (Out1$GD)
  OMeanSNPs <- mean(Out1$SNPs)
  OSdSNPs <- sd(Out1$SNPs)
  LogOut1 <- data.frame(OMeanGD,OSdGD,OMeanSNPs,OSdSNPs)
  colnames(LogOut1) <- c("Overall mean GDij|POk","Standard deviation GDij|POk","Overall mean usable loci","Standard deviation usable loci")
  assign("Triad_Log1", LogOut1, envir = envir)
  
  Geno <- vector(mode="numeric",length=0)
  MeanGDGeno <- vector(mode="numeric",length=0)
  MinGDGeno <- vector(mode="numeric",length=0)
  MaxGDGeno <- vector(mode="numeric",length=0)
  MeanSNPsGeno <- vector(mode="numeric",length=0)
  for (i in 1:nrow(GK)) {
    genoOut <- subset (Out1, Out1$Parent1 == GK$genos[i])
    Geno <- append(Geno,as.character(GK$genos[i]))
    MeanGDGeno <- append(MeanGDGeno,mean(genoOut$GD))
    MinGDGeno <- append(MinGDGeno,min(genoOut$GD))
    MaxGDGeno <- append(MaxGDGeno,max(genoOut$GD))
    MeanSNPsGeno <- append(MeanSNPsGeno,mean(genoOut$SNPs))
  }
  LogOut2 <- data.frame(Geno,MeanGDGeno,MinGDGeno,MaxGDGeno,MeanSNPsGeno)
  colnames(LogOut2) <- c("Genotype","Mean GDij|POk","Min GDij|POk","Max GDij|POk","Mean usable loci")
  assign("Triad_Log2", LogOut2, envir = envir)
  if (files==TRUE) {
    blank <- data.frame(" "," "," "," "," "," ")
    write.table(LogOut1,"apparent-Triad-Summary.csv",sep=",",col.names=T,row.names=F)
    write.table(blank,"apparent-Triad-Summary.csv",sep=",",col.names=F,row.names=F,append=T)
    write.table(blank,"apparent-Triad-Summary.csv",sep=",",col.names=F,row.names=F,append=T)
    write.table(LogOut2,"apparent-Triad-Summary.csv",sep=",",col.names=T,row.names=F,append=T)
  } 

  # Finding the GAP, testing its significance and calculating p-values.
  Out1a <- na.omit(Out1[order(Out1$GD),])
  Out1b <- subset (Out1a, Out1a$GD <= MaxIdent)
  Out1b <- subset (Out1b, Out1b$SNPs > nloci)
  
  if (self == FALSE) {
    Out1b <- subset (Out1b, Out1b$Cross.Type != 'self')
  }
  
  if (nrow(Out1b) > 0) {
    Tdiff <- vector(mode="numeric",length=0)  
    for (i in 1:nrow(Out1b)) {
      Tdiff[i] <- Out1b$GD[i+1] - Out1b$GD[i]
    }
  } else {
    print ("No triads (pair of parents + offspring) were found.")
  }
  Tpv <- vector(mode="numeric",length=0)
  TMax <- max(na.omit(Tdiff))
  TIndex <- match(TMax,Tdiff)
  Tvect1 <- c(sample(na.omit(Tdiff[-TIndex]),29,replace=T),TMax)
  TDtGap <- dixon.test(Tvect1)
  
  if (TDtGap$p.value < alpha ) {
    TCutoff <- Out1b$GD[TIndex]
    L <- subset (Out1b, Out1b$GD <= TCutoff)
    H <- subset (Out1b, Out1b$GD > TCutoff)
    
    if (nrow(L) > 0) {    
      if (nrow(H) > 29) {
        S <- H$GD[1:29]
      } else if (nrow(H) < 3) {
        S <- sample(H$GD,29,replace=T)
      } else if ( (nrow(H) > 3) & (nrow(H) < 29) ) {
        S <- H$GD[1:nrow(H)]
      } 
      for (i in 1:nrow(L)) { 
        Tvect2 <- c(S, L$GD[i])
        TDt <- dixon.test(Tvect2)
        Tpv <- append (Tpv,TDt$p.value)
      }
      
      # Significant triads
      Out2 <- data.frame(L, Tpv)
      colnames(Out2) <- c("Parent1","Parent2","Offspring","Cross.Type","SNPs","GD","p.value")
      Out2a <- subset(Out2,Out2$p.value < alpha)

      # Check hits with duplicated offspring and exclude them
      dupl.offs <- data.frame(unique(Out2a$Offspring[duplicated(Out2a$Offspring)]))
      if ((nrow(dupl.offs)) > 1) {
        for (i in 1:nrow(dupl.offs)){
          Out1c <- subset (Out1b, Out1b$Offspring != dupl.offs[i,1])
        }
        Tdiff <- vector(mode="numeric",length=0)  
        for (i in 1:nrow(Out1c)) {
          Tdiff[i] <- Out1c$GD[i+1] - Out1c$GD[i]
        }
        Tpv <- vector(mode="numeric",length=0)
        TMax <- max(na.omit(Tdiff))
        TIndex <- match(TMax,Tdiff)
        Tvect1 <- c(sample(na.omit(Tdiff[-TIndex]),29,replace=T),TMax)
        TDtGap <- dixon.test(Tvect1)
        
        if (TDtGap$p.value < alpha ) {
          TCutoff <- Out1c$GD[TIndex]
          L <- subset (Out1c, Out1c$GD <= TCutoff)
          H <- subset (Out1c, Out1c$GD > TCutoff)
          
          if (nrow(L) > 0) {    
            if (nrow(H) > 29) {
              S <- H$GD[1:29]
            } else if (nrow(H) < 3) {
              S <- sample(H$GD,29,replace=T)
            } else if ( (nrow(H) > 3) & (nrow(H) < 29) ) {
              S <- H$GD[1:nrow(H)]
            } 
            for (i in 1:nrow(L)) { 
              Tvect2 <- c(S, L$GD[i])
              TDt <- dixon.test(Tvect2)
              Tpv <- append (Tpv,TDt$p.value)
            }
          }
        }
      }
      # Output 2 - Unique and significant Triads
      Out2 <- data.frame(L, Tpv)
      colnames(Out2) <- c("Parent1","Parent2","Offspring","Cross.Type","SNPs","GD","p.value")
      Out2a <- subset(Out2,Out2$p.value < alpha)
      Out2a <- Out2a[order(Out2a$p.value),]
      Triad_Sig <- Out2a
      assign("Triad_Sig", Out2a, envir = envir)
      if (files==TRUE){
        write.csv(Triad_Sig,"apparent-Triad-Sig.csv",row.names=F)
      }
    }
  } else {
    print("The Triad analysis GAP was not significant at the declared alpha level. No triads (pair of parents + offspring) were found.")
  }
  
  # Print the Triad analysis plots
  SortGD <- as.data.frame(na.omit(sort(GD)))
  colnames(SortGD) <- "GD"
  ThresholdT <- mean(c(SortGD[TIndex + 1,1], SortGD[TIndex,1]))
  SortGD$Colour[SortGD$GD <= ThresholdT] = "red"
  SortGD$Colour[SortGD$GD > ThresholdT] = "black"
  
  if (plot==TRUE & files==FALSE) {
    plot (SortGD$GD,xlab="Test triads, ordered by GDij|POk",ylab="Gower Genetic Dissimilarity (GDij|k)",
          main="Triad analysis plot",pch=1,cex=.5,col=SortGD$Colour,xaxt="n")
    axis(1,at=c(1,nrow(SortGD)),labels=c("1",nrow(SortGD)),cex.axis=.7)
    abline(h = ThresholdT, lty = 2, col =  "tomato")
  
  } else if (plot==TRUE & files==TRUE) {
    plot (SortGD$GD,xlab="Test triads, ordered by GDij|POk",ylab="Gower Genetic Dissimilarity (GDij|k)",
          main="Triad analysis plot",pch=1,cex=.5,col=SortGD$Colour,xaxt="n")
    axis(1,at=c(1,nrow(SortGD)),labels=c("1",nrow(SortGD)),cex.axis=.7)
    abline(h = ThresholdT, lty = 2, col =  "tomato")
    
    png("apparent-Triad-Plot.png", width=6,height=6,units='in',res=400)
    par(mar = c(5,5,4,2), xpd = F,cex.axis=.6,cex.lab=.6,cex.main=.8)
    plot (SortGD$GD,xlab="Test triads, ordered by GDij|POk",ylab="Gower Genetic Dissimilarity (GDij|k)",
          main="Triad analysis plot",pch=1,cex=.5,col=SortGD$Colour,xaxt="n")
    axis(1,at=c(1,nrow(SortGD)),labels=c("1",nrow(SortGD)),cex.axis=.7)
    abline(h = ThresholdT, lty = 2, col =  "tomato")
    dev.off()
  }
  
  ##########################################################
  # Dyad analysis - for the Triad analysis non-assignments
  ##########################################################
  
  if (Dyad==TRUE) {
    print ("Starting the Dyad analysis...")
    AllOffs <- as.character(colnames(Offsprings))
    OffsAssigned <- as.character(unique(Out2a[,3]))
    OffsD <- as.data.frame(setdiff(AllOffs,OffsAssigned))
    
    if (nrow(OffsD) < 1) {
      stop ("No Dyad analysis is needed because the Triad analysis identified the parents for all tested offspring. Please refer to the output file apparent-Triad-Sig.txt.")
    } else {
      ParD <- as.data.frame(unique(Out1[,1]))
      colnames(OffsD) <- "X"
      colnames(ParD) <- "X"
      GDM <- data.frame(matrix(ncol=nrow(OffsD),nrow=nrow(ParD)))
      GDCV <- data.frame(matrix(ncol=nrow(OffsD),nrow=nrow(ParD)))
      
      for (i in 1:nrow(ParD)) {  
        Pa1 <- Out1[grep(ParD[i,1],Out1$Parent1),c(1,3,6)]
        Pa2 <- Out1[grep(ParD[i,1],Out1$Parent2),c(2,3,6)]
        colnames (Pa1) <- c("P","O","GD")
        colnames (Pa2) <- c("P","O","GD")
        Pa1Pa2 <- rbind(Pa1,Pa2)  
        for (j in 1:nrow(OffsD)) {
          Samp <- head(Pa1Pa2[grep(OffsD$X[j],Pa1Pa2$O),3],-1)
          avg <- mean(Samp)
          dev <- sd(Samp)
          
          # Calculating GDM and GDCV values
          dat <- ConvertedData[c(as.character(ParD[i,]),as.character(OffsD[j,]))]
          dat <- na.omit(dat)
          colnames(dat) <- c("geno1","geno2")
          dat$diff <- abs(dat$geno1 - dat$geno2)
          dist <- sum(dat$diff) / length(dat$diff)
          ratio <- dev/dist      
          GDM[i,j] <- avg
          GDCV[i,j] <- ratio
          colnames(GDM) <- as.character(OffsD$X)
          rownames(GDM) <- as.character(ParD$X)
          colnames(GDCV) <- as.character(OffsD$X)
          rownames(GDCV) <- as.character(ParD$X)
        } 
      }
      
      # Assign p-values to each significant parent-offspring pair   
      DPa <- vector(mode="numeric",length=0)
      DOf <- vector(mode="numeric",length=0)
      DMPv <- vector(mode="numeric",length=0)
      DRPv <- vector(mode="numeric",length=0)
      DCumPv <- vector(mode="numeric",length=0)
      
      for (i in 1:ncol(GDM)) {
        M <- cbind(as.matrix(rownames(GDM)),as.data.frame(GDM[,i]))
        colnames(M) <- c("P","Mean")
        M <- na.omit(M[order(M$Mean),])      
        R <- cbind(as.matrix(rownames(GDCV)),as.data.frame(GDCV[,i]))
        colnames(R) <- c("P","Ratio")
        R <- na.omit(R[order(R$Ratio),])
        
        # Testing GDMs
        Mp <- as.data.frame(pnorm ((M$Mean - mean(M$Mean)) / sd (M$Mean)))
        MpVect <- cbind(M$P,Mp)
        colnames(MpVect) <- c("P","Mpnorm")
        MpSig <- subset(MpVect,MpVect$Mpnorm < sqrt(alpha))
        if (nrow(MpSig) < 1) {
          next
        }
        
        # Testing GDCVs
        Rp <- as.data.frame(pnorm ((R$Ratio - mean(R$Ratio)) / sd (R$Ratio)))
        RpVect <- cbind(R$P,Rp)
        colnames(RpVect) <- c("P","Rpnorm")
        RpSig <- subset(RpVect,RpVect$Rpnorm > (1 - sqrt(alpha)))
        
        # Genotype(s) that passed on both GDM and GDCV tests
        MRSig <- merge(MpSig,RpSig, by = "P")
        
        if (nrow(MRSig) < 1) {
          next
        } else {
          for (j in 1: nrow(MRSig)) {
            DPa <- append(DPa,as.character(MRSig$P[j]))
            DOf <- append(DOf,colnames(GDM[i]))
            DMPv <- append(DMPv,MRSig$Mpnorm[j])
            DRPv <- append(DRPv,(1 - MRSig$Rpnorm[j]))
            DCumPv <- append(DCumPv,(MRSig$Mpnorm[j] * (1 - MRSig$Rpnorm[j])))
          }     
        }
      }
      
      # Dyad analysis output 
      Out3 <- data.frame(DPa,DOf,DMPv,DRPv,DCumPv)
      Out3 <- subset(Out3, ((Out3$DMPv < alpha) & (Out3$DCumPv < alpha)) | ((Out3$DRPv < alpha) & (Out3$DCumPv < alpha)) )
      Out3 <- Out3[order(Out3$DCumPv),]
      
      if (nrow(Out3) < 1) {
        stop ("At the declared alpha level, the Dyad analysis was unable to find any signficant parent-offspring associations.")
      
      } else if (nrow(Out3) == 1) {
        colnames(Out3) <- c("Parent","Offspring","GDM p-value","GDCV p-value","Cumulative p-value")
        assign("Dyad_Sig", Out3b, envir = envir)
        Out3b <- Out3
        
      } else if (nrow(Out3) > 1) {
        # Skip all parent-offspring pairs already considered, as well as all offspring successfuly assigned to parental pairs in the Triad analysis.     
        Out3a <- data.frame(matrix(ncol=7,nrow=0))
        Out3b <- data.frame(matrix(ncol=5,nrow=0))
        Out3$D1 <- paste(Out3$DPa,Out3$DOf,sep=".")
        Out3$D2 <- paste(Out3$DOf,Out3$DPa,sep=".")
        
        # If a duplicated analysis, keep the likely one (strongest evidence > 200)
        for (i in 1:nrow(Out3)) {
          if (Out3$D2[i] %in% Out3$D1) {
            j <- grep(Out3$D2[i], Out3$D1)
            if ( (Out3[j,5] / Out3[i,5] ) > 200 ) {
              Out3a <- rbind (Out3a,Out3[i,])
            } else {
              Out3a <- rbind (Out3a,Out3[j,])
            }
            Out3a <- rbind (Out3a,Out3[i,])
          }
        }
        Out3a <- Out3a[order(Out3a$DOf,Out3a$DCumPv),]
        
        for (i in 1:nrow(Out3a)) {
          if (i == nrow(Out3a)) {
            break
          } else {
            Occur <- length(grep(Out3a$DOf[i],Out3a$DOf))     
            if (Occur == 1) {
              Out3b <- rbind (Out3b,Out3a[i,c(1:5)])
            } else if (Occur > 1) {
              if ( (Out3a[i+1,2] != Out3a[i,2]) && ((Out3a[i+1,5] / Out3a[i,5]) >= 200) ) {
                Out3b <- rbind (Out3b,Out3a[i,c(1:5)])
              }
            }
          }
        }
        Out3b <- Out3b[,1:5]
        Out3b <- Out3b[order(Out3b$DCumPv),]
        colnames(Out3b) <- c("Parent","Offspring","GDM p-value","GDCV p-value","Cumulative p-value")
        Dyad_Sig <- Out3b
        assign("Dyad_Sig", Out3b, envir = envir)
        if (files==TRUE){
          write.csv(Dyad_Sig,"apparent-Dyad-Sig.csv",row.names=F)
        }
      }
      
      # Plotting the Dyad results
      if ( plot==TRUE & files==TRUE & nrow(Out3b > 0) ) {
        pdf("apparent-Dyad-Plot.pdf",width=7,height=4)
        for (i in 1:nrow(Out3b)) {
          par(mfrow=c(2,1))
          par(mar=c(2,1,1,1))
          par(mgp=c(.5,.5,0))
          par(oma=c(0,0,0,0))
          PlotM <- as.data.frame(GDM[,as.character(Out3b$Offspring[i])])
          PlotR <- as.data.frame(GDCV[,as.character(Out3b$Offspring[i])])
          PlotM <- cbind(rownames(GDM),PlotM)
          PlotR <- cbind(rownames(GDCV),PlotR)
          colnames(PlotM) <- c("P","Mean")
          colnames(PlotR) <- c("P","Ratio")
          PlotM <- na.omit(PlotM[order(PlotM$Mean),])
          PlotR <- na.omit(PlotR[order(PlotR$Ratio),])
          # Normal probabilities of Means (GDM)
          PlotMpnorm <- as.data.frame(pnorm((PlotM$Mean - mean(PlotM$Mean)) / sd (PlotM$Mean)))
          PlotM <- cbind(PlotM,PlotMpnorm)
          colnames(PlotM) <- c("P","Mean","Mpnorm")
          PlotM$Colour[PlotM$Mpnorm <= sqrt(alpha)] = "red"
          PlotM$Colour[PlotM$Mpnorm > sqrt(alpha)] = "black"
          M_axis <- PlotM[c(1,nrow(PlotM)),c(1,2)]
          cntM <- grep("red",PlotM$Colour)
          M_lower_bound <- -qnorm(1 - sqrt(alpha)) * sd(PlotM$Mean) + mean(PlotM$Mean)
          # # Normal probabilities of Ratio (GDCV)
          PlotRpnorm <- as.data.frame(pnorm((PlotR$Ratio - mean(PlotR$Ratio)) / sd (PlotR$Ratio)))
          PlotR <- cbind(PlotR,PlotRpnorm)
          colnames(PlotR) <- c("P","Ratio","Rpnorm")
          PlotR$Colour[PlotR$Rpnorm <= 1-(sqrt(alpha))] = "black"
          PlotR$Colour[PlotR$Rpnorm > 1-(sqrt(alpha))] = "red"
          R_axis = PlotR[c(1,nrow(PlotR)),c(1,2)]
          cntR <- grep("red",PlotR$Colour)
          R_upper_bound <- qnorm(1 - sqrt(alpha)) * sd(PlotR$Ratio) + mean(PlotR$Ratio)
          # GDM plot
          plot (PlotM$Mean,rep(1,nrow(PlotM)),xlab="GDM",col=PlotM$Colour,ylab="", 
                main=Out3b$Offspring[i],pch=1,cex=.5,cex.lab=.7,axes=F,line=-1)
          axis(side=1,at=M_axis$Mean,labels=round(M_axis$Mean,2),cex.axis=.5,line=-3)
          segments(x0=M_lower_bound, y0=.85, x1=M_lower_bound, y1=1,col='red',lty=3)
          #text(M_text,.7,labels="GDM",pos=3,cex=.4,offset=.1)
          text(PlotM$Mean[cntM],1,labels=PlotM$P[cntM],pos=4,cex=.4,srt=45,offset=.3)
          text(M_lower_bound,.7,labels="Lower bound\ncutoff",pos=3,cex=.4,offset=.3)
          # GDCV plot
          plot (PlotR$Ratio,rep(1,nrow(PlotR)),xlab="GDCV",col=PlotR$Colour,ylab="",
                pch=1,cex=.5,cex.lab=.7,axes=F,line=-1)
          axis(side=1,at=R_axis$Ratio,labels=round(R_axis$Ratio,2),cex.axis=.5,line=-3)
          segments(x0=R_upper_bound, y0=.85, x1=R_upper_bound, y1=1,col='red',lty=3)
          text(PlotR$Ratio[cntR],1,labels=PlotR$P[cntR],pos=2,cex=.4,srt=315,offset=.3)
          text(R_upper_bound,.7,labels="Upper bound\ncutoff",pos=3,cex=.4,offset=.3)
        }
        dev.off()
      }
    }
  }
  print("Done.")
  proc.time() - ptm
} 
