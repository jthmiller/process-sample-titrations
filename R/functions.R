## functions

run_titrator <- function(data_dir, sal_mass_file, date_cal, cal_file){
  
  Mass<-read.csv(sal_mass_file, header=T, sep=",", na.string="NA", as.is=T) 
  # only keep sample ID and weight and remove the extra stuff at the bottom
  Mass <- Mass[1:nrows,c('Sample.ID1','Weight..g.','Sample.Index', 'Salinity')]
  ## parse through all the data in the one file ###
  sample_names <- Mass$Sample.ID1
  
  #### pH Calibration #####
  pHCal<-read.csv(cal_file) # read in the pH Calibration file
  
  ## read in sample files
  files_to_parse <- list.files(data_dir, pattern = "*.csv")
  
  #select the calibration for the correct date
  pHData<-pHCal[pHCal$Date==date,]
  
  # calculate pH 3 and 3.5 based on the slope and intercept from pH 4, 7, and 10 calibration
  mod.pH<-lm(c(pHData$pH4, pHData$pH7, pHData$pH10)~c(4,7,10)) # linear model
  
  
  png(paste0(gsub('\\/','_',date),'_pHmvplot.png'), height = 400, width = 400)
  plot(c(4,7,10), c(pHData$pH4, pHData$pH7, pHData$pH10), xlab = 'pH', ylab = 'mv')
  lines(c(4,7,10), predict(mod.pH))
  R2<-summary(mod.pH)$r.squared
  legend('topright', legend = bquote(R^2 == .(format(R2, digits = 3))), bty='n')
  dev.off()
  
  # Select the mV for pH=3 and pH=3.5 based on your probe calibration
  pH35<-mod.pH$coefficients[1]+mod.pH$coefficients[2]*3.5
  pH3<-mod.pH$coefficients[1]+mod.pH$coefficients[2]*3
  
  tables_out <- lapply(files_to_parse, return_tables)
  ph <- lapply(tables_out,'[[',1)
  ph <- do.call(rbind,ph)
  
  set <- lapply(tables_out,'[[',2)
  set <- do.call(rbind,set)
  AllData <- set
  
  AllData <- AllData %>% 
    separate(sampname, c("sampname", "date_time"), "_", extra="merge")
  # Identifies rows starting with "Scope" in column 1
  colnames(AllData) <- c("ID1","ID2","Time", "mV", "Volume", "dV/dt", "Temperature")
  AllData$mV <- as.numeric(AllData$mV) ## supress the warnings since NA will be produced
  AllData$Temperature <- as.numeric(AllData$Temperature)
  AllData$Volume <- as.numeric(AllData$Volume)
  
  sample_names_list <- split( AllData , f = AllData$ID1)
  
  ##### titration##########
  #create an empty matrix to put the TA values in
  nrows<-nrow(Mass) # -1 because there is an extra line in the mass file
  TA <- data.frame(matrix(nrow = nrows, ncol = 4))
  rownames(TA)<-Mass$Sample.ID1[1:nrows]
  colnames(TA)<-c("SampleID",'TA','Mass','Sample.Index')
  
  out <- lapply(names(sample_names_list), function(name, data = sample_names_list, acid = 0.100025) {
    Data <- data[[name]]
    mV <- which(Data$mV<pH3 & Data$mV>pH35) 
    d_in <- (-0.00000400*mean(Data$Temperature[mV], na.rm=T)^2-0.0001108*mean(Data$Temperature[mV], na.rm=T)+1.02878) #A23
    s <- Mass[Mass$Sample.ID1==name,4]
    mass <- Mass[Mass$Sample.ID1==name,2]
    sample.index <- Mass[Mass$Sample.ID1==name,3]# this is the order that the sample was run
    tot_alk <- 1000000*at(S=s,T=mean(Data$Temperature[mV], na.rm=T), C = acid, d = d_in, pHTris = NULL, ETris = NULL, weight = mass, E=Data$mV[mV], volume=Data$Volume[mV])
    return(data.frame(name, tot_alk, mass, s, sample.index))
  })
  
  out <- data.frame(do.call(rbind, out) )
  rownames(out) = out$name
  
  Mass$salt <- seq(from = 31, length.out = length(sample_names_list), by = 0.02)
  rownames(Mass) <- Mass$Sample.ID1
  
  out[rownames(Mass),'evap'] <- out[rownames(Mass) ,'tot_alk']*31/Mass[ rownames(Mass),'salt']
  
  #exports your data as a CSV file
  write.table(out,"07212023/data_output/practiceoutput 20230721.csv",sep=",", row.names=FALSE)
  
}

##############################################################################################
# at function
#at function is based on code in saecarb package by Steeve Comeau, Heloise Lavigne and Jean-Pierre Gattuso
##############################################################################################

at<-function (S = 35, T = 25, C = 0.1, d = 1, pHTris = NULL, ETris = NULL, 
              weight, E, volume) 
{
  R <- 8.31447215
  F <- 96485.339924
  test <- length(T) != length(E)
  if (test) {
    T <- rep(T[1], length(E))
  }
  Tk <- T + 273.15
  p <- data.frame(E = E, volume = volume, Tk = Tk)
  z <- p
  if (!is.null(pHTris) & !is.null(ETris)) {
    pH <- pHTris + (ETris/1000 - E/1000)/(R * Tk * log(10)/F)
    p <- data.frame(p, pH = pH)
    iii <- which((3 <= p$pH) & (p$pH <= 3.5))
    z <- p[iii, ]
  }
  options(digits = 9)
  m <- z$volume * d
  m0 <- weight
  
  F1 <- (m0 + m) * exp((z$E/1000)/(R * (z$Tk)/F))
  print(F1)
  print(F)
  f <- lm(m ~ F1)
  TA <- f$coefficients[1] * C/m0[1]
  E0 <- z$E/1000 - (R * z$Tk/F) * log((-m0 * TA + m * C)/(m0 + 
                                                            m))
  Hprime <- exp((z$E/1000 - E0)/(R * z$Tk/F))
  Cl <- S/1.80655
  St <- 0.14 * Cl/96.062
  Ks <- Ks(S, T[1], 0)
  Z <- 1 + St/Ks
  Ft <- 6.7e-05 * Cl/18.9984
  Kf <- exp(874/z$Tk - 9.68 + 0.111 * S^(1/2))
  y <- (m/m0)
  regr <- nls(y ~ ((At + (St/(1 + Ks * Z/(f * Hprime))) + (Ft/(1 + 
                                                                 Kf/(f * Hprime))) + (f * Hprime/Z))/(C - f * Hprime/Z)), 
              start = list(At = TA, f = 1))
  ALK <- summary(regr)$parameters[1]
  attr(ALK, "unit") <- "mol/kg-soln"
  attr(ALK, "name") <- "Total Alkalinity"
  return(ALK)
}
##############################################################################################
##############################################################################################
##############################################################################################

return_tables <- function(file){
  mess <- readLines(paste(data_dir,"/", file, sep=""), encoding = 'iso-8859-1')
  mess <- iconv(mess,'iso-8859-1', 'utf-8')
  mess <- mess[!mess == ""]
  table_data_names <- grep('^SET|^pH',mess, value = T)
  
  mess <- mess[grep('^SET',mess, invert = T)]
  mess <- mess[grep('^pH',mess, invert = T)]
  
  headers <- grep('^Time', mess)
  sheet <- length(mess)
  
  data.ends <- lapply(headers,function(head){
    if(any(headers > head)){
      min(headers[which(headers > head)]) - 1
    } else {
      sheet
    }
  })
  data.ends <- unlist(data.ends)
  
  tables_of_data <- mapply(function(head, tail){
    data <- strsplit(mess[head:tail ],",")
    data <- data.frame(do.call(rbind,data))
    colnames(data) <- data[1,]
    data <- data[-1,]
    return(data)
  }, headers, data.ends)
  names(tables_of_data) <- unlist(lapply(strsplit(table_data_names, ","),'[[',1))
  
  tables_of_data <- lapply(tables_of_data, function(data){
    return(data.frame(sampname = gsub('.csv','',file), data))
  })
  
  return(tables_of_data)
}

##############################################################################################
##############################################################################################
##############################################################################################