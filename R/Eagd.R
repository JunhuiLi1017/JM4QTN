Eagd <-
function(marType,croType,Gn=2,x,y=0){   
  if(Gn < 1){
    stop("Gn should > 0")
  }
  z <- x + y -2*x*y
  if(marType == "22"){
    if (croType=="BCP1"){
      rval <- FreGeno(croType,Gn+1,8,x,y)/(FreGeno(croType,Gn+1,8,x,y)+FreGeno(croType,Gn+1,6,x,y))
    }
    if (croType=="DH"){
      rval <- (1-x-y)/(1-z)
    }
    if (croType=="F2"|croType=="RIL"){
      rval <- ((1-x)*(1-x)*(1-y)*(1-y) - x*x*y*y)/((1-z)*(1-z))
    }
    if (croType=="Fn"){
      rval <- (FreGeno(croType,Gn,1,x,y) - FreGeno(croType,Gn,3,x,y))/(FreGeno(croType,Gn,1,x,y) + FreGeno(croType,Gn,6,x,y) + FreGeno(croType,Gn,3,x,y))
    }
  }else if(marType == "21"){
    if (croType=="BCP1"){
      rval <- FreGeno(croType,Gn+1,7,x,y)/(FreGeno(croType,Gn+1,7,x,y)+FreGeno(croType,Gn+1,5,x,y))
    }
    if (croType=="F2"){
      rval <- ((1-x)*(1-x)*(1-y)*y - x*x*y*(1-y))/((1-z)*z)
    }
    if (croType=="Fn"){
      rval <- (FreGeno(croType,Gn,5,x,y) - FreGeno(croType,Gn,15,x,y))/(FreGeno(croType,Gn,5,x,y) + FreGeno(croType,Gn,7,x,y) + FreGeno(croType,Gn,11,x,y) + FreGeno(croType,Gn,15,x,y))
    }
  }else if(marType == "20"){
    if (croType=="DH"){
      rval <- (y-x)/z
    }
    if (croType=="F2"|croType=="RIL"){
      rval <- ((1-x)*(1-x)*y*y - x*x*(1-y)*(1-y))/(z*z)
    }
    if (croType=="Fn"){
      rval <- (FreGeno(croType,Gn,2,x,y) - FreGeno(croType,Gn,4,x,y))/(FreGeno(croType,Gn,2,x,y) + FreGeno(croType,Gn,12,x,y) + FreGeno(croType,Gn,4,x,y))
    }
  }else if(marType == "12"){
    if (croType=="BCP1"){
      rval <- FreGeno(croType,Gn+1,4,x,y)/(FreGeno(croType,Gn+1,4,x,y)+FreGeno(croType,Gn+1,2,x,y))
    }
    if (croType=="F2"){
      rval <- ((1-y)*(1-y)*x*(1-x) - x*(1-x)*y*y)/(z*(1-z))
    }
    if (croType=="Fn"){
      rval <- (FreGeno(croType,Gn,8,x,y) - FreGeno(croType,Gn,14,x,y))/(FreGeno(croType,Gn,8,x,y) + FreGeno(croType,Gn,10,x,y) + FreGeno(croType,Gn,16,x,y) + FreGeno(croType,Gn,14,x,y))
    }
  }else if(marType == "11"){
    if (croType=="BCP1"){
      rval <- FreGeno(croType,Gn+1,3,x,y)/(FreGeno(croType,Gn+1,3,x,y)+FreGeno(croType,Gn+1,1,x,y))
    }
    if (croType=="BCP2"){
      rval <- -(FreGeno(croType,Gn+1,6,x,y)/(FreGeno(croType,Gn+1,8,x,y)+FreGeno(croType,Gn+1,6,x,y)))
    }
    if (croType=="F2"|croType=="Fn"){
      rval <- 0
    }
  }else if(marType == "10"){
    if (croType=="BCP2"){
      rval <- -(FreGeno(croType,Gn+1,5,x,y)/(FreGeno(croType,Gn+1,7,x,y)+FreGeno(croType,Gn+1,5,x,y)))
    }
    if (croType=="F2"){
      rval <- (x*(1-x)*y*y -x*(1-x)*(1-y)*(1-y))/(z*(1-z))
    }
    if (croType=="Fn"){
      rval <- (FreGeno(croType,Gn,14,x,y) - FreGeno(croType,Gn,8,x,y))/(FreGeno(croType,Gn,14,x,y) + FreGeno(croType,Gn,10,x,y) + FreGeno(croType,Gn,16,x,y) + FreGeno(croType,Gn,8,x,y))
    }
  }else if(marType == "02"){
    if (croType=="DH"){
      rval <- (x-y)/z
    }
    if (croType=="F2"|croType=="RIL"){
      rval <- ((1-y)*(1-y)*x*x - (1-x)*(1-x)*y*y)/(z*z)
    }
    if (croType=="Fn"){
      rval <- (FreGeno(croType,Gn,4,x,y) - FreGeno(croType,Gn,2,x,y))/(FreGeno(croType,Gn,4,x,y) + FreGeno(croType,Gn,12,x,y) + FreGeno(croType,Gn,2,x,y))
    }
  }else if(marType == "01"){
    if (croType=="BCP2"){
      rval <- -(FreGeno(croType,Gn+1,2,x,y)/(FreGeno(croType,Gn+1,4,x,y)+FreGeno(croType,Gn+1,2,x,y)))
    }
    if (croType=="F2"){
      rval <- (x*x*(1-y)*y - (1-x)*(1-x)*(1-y)*y)/((1-z)*z)
    }
    if (croType=="Fn"){
      rval <- (FreGeno(croType,Gn,15,x,y) - FreGeno(croType,Gn,5,x,y))/(FreGeno(croType,Gn,5,x,y) + FreGeno(croType,Gn,7,x,y) + FreGeno(croType,Gn,11,x,y) + FreGeno(croType,Gn,15,x,y))
    }
  }else if(marType == "00"){
    if (croType=="BCP2"){
      rval <- -(FreGeno(croType,Gn+1,1,x,y)/(FreGeno(croType,Gn+1,3,x,y)+FreGeno(croType,Gn+1,1,x,y)))
    }
    if (croType=="DH"){
      rval <- (x+y-1)/(1-z)
    }
    if (croType=="F2"|croType=="RIL"){
      rval <- (x*x*y*y - (1-x)*(1-x)*(1-y)*(1-y))/((1-z)*(1-z))
    }
    if (croType=="Fn"){
      rval <- (FreGeno(croType,Gn,3,x,y) - FreGeno(croType,Gn,1,x,y))/(FreGeno(croType,Gn,1,x,y) + FreGeno(croType,Gn,6,x,y) + FreGeno(croType,Gn,3,x,y))
    }
  }else if(marType == "2N"){
    if (croType=="BCP1"){
      rval <- FreGeno(croType,Gn,4,x)/(FreGeno(croType,Gn,3,x)+FreGeno(croType,Gn,4,x))
    }
    if (croType=="DH"){
      rval <- 1-2*x
    }
    if (croType=="Fn"){
      rval <- (FreGeno(croType,Gn,9,x)-FreGeno(croType,Gn,7,x))/(FreGeno(croType,Gn,9,x)+FreGeno(croType,Gn,8,x)+FreGeno(croType,Gn,7,x))
    }
    if (croType=="F2"){
      rval <- (1-x)*(1-x) - x*x
    }
    if (croType=="RIL"){
      rval <- ((1-x)*(1-x) - x*x)/((1-x)*(1-x) + x*x)
    }
  }else if(marType == "1N"){
    if (croType=="BCP1"){
      rval <- FreGeno(croType,Gn,2,x)/(FreGeno(croType,Gn,1,x)+FreGeno(croType,Gn,2,x))
    }
    if (croType=="BCP2"){
      rval <- -FreGeno(croType,Gn,3,x)/(FreGeno(croType,Gn,4,x)+FreGeno(croType,Gn,3,x))
    }
    if (croType=="F2"|croType=="Fn"){
      rval <- 0
    }
  }else if(marType == "0N"){
    if (croType=="BCP2"){
      rval <- -FreGeno(croType,Gn,1,x)/(FreGeno(croType,Gn,1,x)+FreGeno(croType,Gn,2,x))
    }
    if (croType=="DH"){
      rval <- 2*x-1
    }
    if(croType=="Fn"){
      rval <- (FreGeno(croType,Gn,3,x)-FreGeno(croType,Gn,1,x))/(FreGeno(croType,Gn,3,x)+FreGeno(croType,Gn,2,x)+FreGeno(croType,Gn,1,x))
    }
    if (croType=="F2"){
      rval <- x*x - (1-x)*(1-x)
    }
    if(croType=="RIL"){
      rval <- (x*x - (1-x)*(1-x))/((1-x)*(1-x) + x*x)
    }
  }else if(marType == "N2"){
    if (croType=="BCP1"){
      rval <- FreGeno(croType,Gn,4,x)/(FreGeno(croType,Gn,4,x)+FreGeno(croType,Gn,2,x))
    }
    if (croType=="DH"){
      rval <- 1-2*x
    }
    if (croType=="Fn"){
      rval <- (FreGeno(croType,Gn,9,x)-FreGeno(croType,Gn,3,x))/(FreGeno(croType,Gn,9,x)+FreGeno(croType,Gn,6,x)+FreGeno(croType,Gn,3,x))
    }
    if (croType=="F2"){
      rval <- (1-x)*(1-x) - x*x
    }
    if (croType=="RIL"){
      rval <- ((1-x)*(1-x) - x*x)/((1-x)*(1-x) + x*x)
    }
  }else if(marType == "N1"){
    if (croType=="BCP1"){
      rval <- FreGeno(croType,Gn,3,x)/(FreGeno(croType,Gn,3,x)+FreGeno(croType,Gn,1,x))
    }
    if (croType=="BCP2"){
      rval <- -FreGeno(croType,Gn,2,x)/(FreGeno(croType,Gn,2,x)+FreGeno(croType,Gn,4,x))
    }
    if (croType=="F2"|croType=="Fn"){
      rval <- 0
    }
  }else if(marType == "N0"){
    if (croType=="BCP2"){
      rval <- -FreGeno(croType,Gn,1,x)/(FreGeno(croType,Gn,1,x)+FreGeno(croType,Gn,3,x))
    }
    if (croType=="DH"){
      rval <- 2*x-1
    }
    if(croType=="Fn"){
      rval <- (FreGeno(croType,Gn,7,x)-FreGeno(croType,Gn,1,x))/(FreGeno(croType,Gn,1,x)+FreGeno(croType,Gn,4,x)+FreGeno(croType,Gn,7,x))
      #error inrval <- (x*x - (1-x)*(1-x))/((1-x)*(1-x) + x*x + (2/(2^(Gn-1)-1))*x*(1-x))
    }
    if (croType=="F2"){
      rval <- x*x - (1-x)*(1-x)
    }
    if(croType=="RIL"){
      rval <- (x*x - (1-x)*(1-x))/((1-x)*(1-x) + x*x)
    }
  }else{
    stop("Wrong mark type for argument marType")
  }
  return(rval)
}
