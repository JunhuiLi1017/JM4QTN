StatPheno <-
function(phenoData,defineForm=NULL,effNotation="G"){
  #Environment(E) effect: multiple year and location;  Block(B) effect: block effect nested in every Environment;  Repetition(R):
  if(ncol(phenoData)<5){
    stop("The number of column should be more than 4!")
  }else{
    colnames(phenoData)[1:4] <- c("E","B","R","G")
    phenoData[,"E"] <- as.factor(phenoData[,"E"])
    phenoData[,"B"] <- as.factor(phenoData[,"B"])
    phenoData[,"G"] <- as.factor(phenoData[,"G"])
  }
  Phelist <- list(list())
  #1 normality abalysis for phenotype
  for(n in 1:(ncol(phenoData)-4)){
    normtest <- shapiro.test(phenoData[,4+n])
    if(normtest$"p.value" < 0.05){
      warning(paste("P value of shapiro.test for ",colnames(phenoData)[4+n]," is ",normtest$"p.value",sep=""))
    }
    Phelist[[colnames(phenoData)[4+n]]]["normality_test"] <- list(normtest)
  }
  #2-3 anova and lsmeans
  if(is.null(defineForm)){
    Form <- NULL
    if(nlevels(phenoData[,"E"])>1){
      if(nlevels(phenoData[,"B"])>1){
        for(i in 1:(ncol(phenoData)-4)){
          Form[i] <- as.character(paste(colnames(phenoData)[4+i],"~ E + G + E:G + B%in%E",sep=""))
        }
      }else{
        Form[i] <- as.character(paste(colnames(phenoData)[4+i],"~ E*G",sep=""))
      }
      #}else if(nlevels(phenoData[,"E"]) == 1 && nlevels(phenoData[,"B"])>1){
      #Form[i] <- as.character(paste(colnames(phenoData)[4+i],"~ G + B",sep=""))
    }else{
      stop("Phenotype data is not multiple year-location data!")
    }
    for(j in 1:(ncol(phenoData)-4)){
      lm0 <- lm(as.formula(Form[j]),data=phenoData)
      ano0 <- anova(lm0)
      lms <- as.data.frame(predict(lsmeans(lm0,"G")))
      colnames(lms) <- colnames(phenoData)[4+j]
      Phelist[[colnames(lms)]]["formula"] <- list(defineForm[j])
      Phelist[[colnames(lms)]]["ANOVA"] <- list(ano0)
      Phelist[[colnames(lms)]]["lsmeans"] <- list(lms)
    }
  }else{
    for(j in 1:length(defineForm)){
      lm0 <- lm(as.formula(defineForm[j]),data=phenoData)
      ano0 <- anova(lm0)
      lms <- as.data.frame(predict(lsmeans(lm0,effNotation)))
      #colnames(lms) <- gsub('\\)|\\(', '', defineForm[j][2])
      colnames(lms) <- unlist(strsplit(defineForm[j], "~"))[1]
      Phelist[[colnames(lms)]]["formula"] <- list(defineForm[j])
      Phelist[[colnames(lms)]]["ANOVA"] <- list(ano0)
      Phelist[[colnames(lms)]]["lsmeans"] <- list(lms)
    }
  }
  #4 varcomponent
  return(Phelist)
}
