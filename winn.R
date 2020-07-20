# WiNN version 0.1

library(hwwntest)
library(gam)

##################################################### 
# Utility function for test.wn
# Find greatest power of 2 that is closest to "n"
# ("closest.pw2", could be > or < "n") and the greatest 
# power of 2 that is closest to "n" and < "n"
# ("greatest.pw2")
#####################################################
greatest.closest.power.2 <- function(n){
  
  i <- 0;
  repeat{
    if( n >= 2^i & n < 2^(i+1)){break}
    else{i <- i+1}
  }
  
  i.closest <- i
  
  # Is this really the power od 2 closest to "n"?
  # (n=255 would return a "power of 2" i value = 7, while it is closer to 2^8)
  if(abs(n - 2^(i+1)) < (n - 2^i)){
    i.closest <- (i+1)  
  }
  
  return(list(closest.pw2=i.closest, greatest.pw2=i))
}


######################################
# White Noise test
# The test combines hwwntest and Box 
# tests
######################################
test.wn <- function(x,box.test.lag=20){
  # The p-value for WN to return
  p <- NA
  
  # Remove NAs
  x.no.na <- x[!is.na(x)]
  
  # The WN test can be applied for series lenghts >= 16  
  if(length(x.no.na) >= 16){
    
    # Determine greatest power of 2 which is < length of data
    z <- greatest.closest.power.2(length(x.no.na))[["greatest.pw2"]]
    
    # Compute delta
    delta <- length(x.no.na) - 2^z
    
    # If delta < 20 , delete "delta" values randomly from the series
    p <- c()
    if(delta < 20){
      
      if(delta > 0){
        x.no.na <- x.no.na[-sample(1:length(x.no.na), delta)]
      }
      
      # Perform the test and return p.value for w.noise (Ho: WN)
      p <- genwwn.test(x.no.na)[["p.value"]]
    }else{
      #print("here 1")
      
      # Create 2 series of length 2^z:
      #   - from position 1 to 2^z
      #   - from position (n - 2^z) to n
      x.no.na.1 <- x.no.na[1:2^z]
      n <- length(x.no.na)
      x.no.na.2 <- x.no.na[(n - 2^z + 1):n]
      
      p1 <- genwwn.test(x.no.na.1)[["p.value"]]
      p2 <- genwwn.test(x.no.na.2)[["p.value"]]
      
      p <- p1; if(p2 <= p1){p <- p2}
    }
  }else{
    print("Error: WN test cannot applied, length of series < 16")
  }
  
  # Add the Box test
  p.box <- Box.test(x, box.test.lag)$p.value

  # Return the smallest of the 2 p-values
  p.ret <- p
  if(!is.na(p.box) & !is.na(p) & p.box <= p){p.ret <- p.box} #print(" ===> Selected p.box")}
  
  if(is.na(p.ret)){
    stop("WN test cannot applied, length of series < 16")}
  
  return(p.ret)
}

###################
# smoothing spline
###################
apply.gam.hastie <- function(y, fun="spline",k){
  
  x <- 1:length(y)
  
  # Use spline or loess in gam
  fun.type <- "s"
  if(fun=="loess"){fun.type -> "lo"}
  #form <- paste("y ~ ", fun.type,"(x)", sep="") s(x, k = -1, bs = "cs"
  form <- paste("y ~ ", fun.type,"(x, ", k , ")", sep="")
  m <- gam(as.formula(form), na.action = na.exclude)
  
  return(predict(m, newdata=data.frame(x=x)))
}

######## 
# Plot
########
plot.uncorr.corr <- function(orig.m, 
                             unresid.metab, 
                             signal, 
                             corrected.signal,
                             seq.order,
                             plate.coord,
                             n.plates,
                             file.name,is.wn, p, k)
{
  correction        <- signal - corrected.signal
  range  <- c(0, as.numeric(plate.coord[length(plate.coord)]))

  pdf(file.name, width = 15, height = 10)
  layout(matrix(c(1,2,3), nrow = 3, ncol = 1, byrow = TRUE))
    
  plot(seq.order,unresid.metab, type="l", xlim=range,col="black",
        xlab="Reading Sequence", ylab="Measurement",
        main=paste(orig.m,": Untransformed")) 
  # Add vertical lines demarcating the plates
  draw.plates <- function(n, plate.coord, y.pos){
    start = 0
    for(j in 1:n){
      abline(v=plate.coord[j], lty=2)
      text(x=start + (plate.coord[j] - start)/2 , y=y.pos, j, col = "blue", cex = 1.2)
      # paste("Plate:", j), col = "blue", cex = 1.2)
      start <- plate.coord[j]
    }
  }
  draw.plates(n=n.plates, plate.coord=as.numeric(plate.coord), y.pos=max(unresid.metab)*0.85)
  
  plot(seq.order, unlist(signal), type="l", xlim=range,col="darkgrey",
        xlab="Reading Sequence", ylab="Residuals",
        main=paste(orig.m,": Normalized by plate sd + residualized by plate"))
  draw.plates(n=n.plates, plate.coord=as.numeric(plate.coord), y.pos=max(signal)*0.85)
  lines(seq.order, unlist(correction), col="red")
  legend("bottomleft", lty = c(1,1), col = c("darkgrey", "red"), 
          legend = c("Residualized", "Correction"))
  
  title <- paste("Corrected (WN Test p-val:", round(p,4), " - ", is.wn," - Best smooth. param.:", k,")")
  plot(seq.order, corrected.signal, type="l", xlim=range,
        xlab="Reading Sequence", ylab="Corrected Residuals",
        main=paste(orig.m,":", title))
  draw.plates(n=n.plates, plate.coord=as.numeric(plate.coord), y.pos=max(corrected.signal)*0.85)
  
  dev.off()
}

########
# WINN
########
winn <- function(met.dat, 
                 stdy.id, 
                 n.start.smooth = 1, 
                 n.stop.smooth = 10,
                 debug=T, 
                 runall=F, 
                 save.res.norm.data =T, 
                 save.corrected.data = T,
                 save.correction.summary = T,
                 save.pdfs=T){
  
  print("Running WINN with following parameters:")
  print(paste("met.dat =", deparse(substitute(met.dat))))
  print(paste("stdy.id =", stdy.id))
  print(paste("n.start.smooth =", n.start.smooth))
  print(paste("n.stop.smooth =", n.stop.smooth))
  print(paste("debug =", debug))
  print(paste("runall =", runall))
  print(paste("save.res.norm.data =", save.res.norm.data))
  print(paste("save.corrected.data =", save.corrected.data))
  print(paste("save.correction.summary =", save.correction.summary))
  print(paste("save.pdfs =", save.pdfs))
  
  # Create study directory if requested
  if(save.res.norm.data | save.corrected.data | save.correction.summary | save.pdfs){
  
    here <- paste(getwd(),"/", stdy.id, sep="")
    dir.create(here)
    
    # If save.pdf == TRUE, create also the plot directory in the study directory
    if(save.pdfs){
    pdf.dir <- paste(here,"/untransf_resid_correct_pdfs",sep="")
    dir.create(pdf.dir)
    }
  }
  
  ################################
  # Step 1: Normalize by plate SD
  ################################
  sd.plate <- function(plate.id, log.dat){
    grep.plates <- grep(paste("plate_",plate.id,"_",sep=""), row.names(log.dat))
    sd.plate <-  lapply(log.dat[grep.plates,], function(x){return(sd(x,na.rm=T))})
    
    sd.plate.df <- data.frame(sd.plate)
    sd.names    <- paste("sd.", names(sd.plate.df), sep="")
    names(sd.plate.df) <- sd.names
    
    sd.plate.df <- sd.plate.df[rep(1, each = length(grep.plates)), ]
    
    return(sd.plate.df)
  }
  
  # Get the plates
  plates <- unique(unlist(lapply(strsplit(rownames(met.dat), "_"), function(x){return(x[2])})))
  
  sd.plate.str <- c()
  
  for(i in 1:length(plates)){
    #print(i);print(plates[i])
    assign(paste("sd.plate.",i,sep=""), sd.plate(plates[i], met.dat))
    sd.plate.str <- c(sd.plate.str, paste("sd.plate.",i,sep=""))
  }
  
  sd.plate.str
  # "sd.plate.1"  "sd.plate.2"  "sd.plate.3" .... "sd.plate.25"
  
  # The rbind command:  "rbind(sd.plate.1 , sd.plate.2 , ... sd.plate.25)
  r.cmd <- paste("rbind(", paste(sd.plate.str, collapse=" , "), ")",sep="")
  
  # Evaluate and assign
  assign("sd.by.plate", eval(parse(text=r.cmd)))
  #str(sd.by.plate)
  
  # Finally normalize here
  norm.met.dt <- met.dat/sd.by.plate
  names(norm.met.dt) <- paste("norm.",names(norm.met.dt), sep="")
  
  #############################################
  # Step 2:  Residualize by plate number
  #############################################
  for(i in 1:(length(plates) - 1)){
    on.plate <- paste("on.plate.",i,sep="")
    norm.met.dt[[on.plate]] <- 0
    norm.met.dt[grep(paste("plate_",plates[i],"_",sep=""), rownames(norm.met.dt)), ][[on.plate]] <- 1
  }
  
  # Regress and get the residuals
  metbs <- names(norm.met.dt)[grep("norm.", names(norm.met.dt))]
  
  for(metb in metbs){
    form <- paste(metb," ~ ", paste(paste("as.factor(on.plate.", 1:(length(plates) - 1), ")",sep=""), collapse=" + "))
    tryCatch(
      {
        reg <- lm(as.formula(form), data=norm.met.dt, na.action=na.exclude)
        
        # Calling "resid" on the regression object created with the na.action=na.exclude
        # option takes care of the missing values. "resid(reg)" return NA if the
        # observation is also NA (remember that reg$resid will only return the residuals for the 
        # non-missing values) 
        norm.met.dt[[paste("res.", metb, sep="")]] <- resid(reg)
      },
      error = function(e){
        print(paste("Error processing", metb))
        print(e)
      })
  }
  
  res.norm.met.dt <- norm.met.dt[grep("res.norm.", names(norm.met.dt), value=T)]
  str(res.norm.met.dt, list.len=5)
  
  if(save.res.norm.data){ 
    save(list=c("res.norm.met.dt"), file=paste(here,"/",stdy.id,"_res.norm.met.RData",sep=""))
  }
  
  ####################################
  # Step 3: WN test (first)
  ####################################
  p.wn.tst <- sapply(names(res.norm.met.dt), function(x){
    ret <- NA; p <- NA;
    tryCatch({
      p <- test.wn(res.norm.met.dt[,x])
    },warning=function(w){
      print(paste("WARNING @ test.wn > ", x,":",w))
    },error = function(e) {
      p <- NA;
      print(paste("ERROR @ test.wn > ", x,":",e)) 
    })
    
    if(!is.na(p)){
      ret <- p
    }
    
    return(ret)})
  
  names(p.wn.tst)  <- names(res.norm.met.dt)
  
  # Determine number of independent dimensions by PCA
  is.wn.before.corr <- c(); 
  do.pca <- function(res.norm.met.dt){
    
    # PCA
    n.90.pc <- tryCatch({
      pc <- prcomp(res.norm.met.dt, scale=TRUE)
      cum.prop    <- summary(pc)$importance["Cumulative Proportion", ]
      print(paste("PCA - number of dimensions:", min(which(cum.prop >= 0.9))))
      min(which(cum.prop >= 0.9)) 
      },warning=function(w){
        print(paste("WARNING @ PCA:", w))
        return(NA)
      },error = function(e) {
        print(paste("ERROR @ PCA", e))
        return(NA)
      })
  
    ret.n.90.pc <- c()
    if(is.na(n.90.pc)){ret.n.90.pc <- dim(res.norm.met.dt)[2]}
    else{ret.n.90.pc <- n.90.pc}

    return(ret.n.90.pc)
  }
  
  # The number of independent dimensions
  print("Running first PCA ...")
  n.90.pc.before.corr <- do.pca(res.norm.met.dt=res.norm.met.dt)
  if(debug){print(paste("DBG 1 - number independent components:", n.90.pc.before.corr))}
  
  # Set the p-value threshold to Bonferroni:
  p.thr.before.corr <- 0.05/n.90.pc.before.corr
  if(debug){print(paste("DBG 2 - p-value threshold:",  p.thr.before.corr))}
    
  # Determine which metab is WN
  is.wn.before.corr <- ifelse(p.wn.tst <= p.thr.before.corr, FALSE, TRUE)
  names(is.wn.before.corr) <- names(res.norm.met.dt)
  
  # Select the metabolites that cannot be tested for WN because
  # errors in the test for white noise function (eg length < 16)
  mets.cannot.test <- names(is.wn.before.corr[is.na(is.wn.before.corr)])
  if(length(mets.cannot.test) > 0){
    print("Metabolites that could not be tested for WN:", log.file=cannot.correct)
    for(i in mets.cannot.test){print(i, log.file=cannot.correct)}
  }
  
  # Select the metabolites that do *not* pass the WN test and must be corrected
  mets.to.correct <- names(is.wn.before.corr[!is.na(is.wn.before.corr) & is.wn.before.corr == FALSE])
  
  print("Metabolites not passing WN to be corrected:")
  for(i in mets.to.correct){print(i)}
  
  # Perform a second PCA restricted to the metabs to be corrected
  print("Running second PCA ...")
  n.90.pc.after.corr <- do.pca(res.norm.met.dt=res.norm.met.dt[mets.to.correct])
  if(debug){print(paste("DBG 3 - number independent components among mets to correct:", n.90.pc.after.corr))}
  
  # n.90.pc.var <- 185
  p.thr.after.corr <- 0.05/n.90.pc.after.corr
  if(debug){print(paste("DBG 4 - p-value threshold for mets to correct:", p.thr.after.corr))}
  
  ####################################
  # Step 4: apply correction
  ####################################
  plates <- as.numeric(plates)

  # If this is true, just run them all regardless of WN test
  if(runall){
    mets.to.correct <- names(res.norm.met.dt)
    p.thr.after.corr <- 0.05/length(res.norm.met.dt)
    n.90.pc.after.corr <- length(res.norm.met.dt)
  }

  # Post-correction WN tets
  is.wn.post.corr <- c()
  count <- 0
  best.pvals <- c()
  
  # If we have metabolites to correct
  if(length(mets.to.correct) > 0){
    
    # Get plate coordinates. It will be used later for plotting purposes
    plate.coord <- c()
    for(i in 1:length(plates)){
      tmp.df <- met.dat[grep(paste("plate_",i,"_order_",sep=""),rownames(met.dat), value=T),]
      plate.coord <- c(plate.coord, tail(sapply(rownames(tmp.df), function(x){strsplit(x,"_order_")[[1]][2]}), n=1))
    }
    
    # Get also the sequence order for later
    seq.order <- as.numeric(sapply(rownames(met.dat), function(x){ strsplit(x,"_order_")[[1]][2]}))

    # Loop through the metabolites
    for(met in mets.to.correct){
     
      # Extract the metabolite
      met.dt<- res.norm.met.dt[met]
      
      print("##########################")
      print(paste("Correcting ", met))
      
      # Update count
      count <- count + 1; if(count %% 20 == 0){print(count)}
      
      pvals           <- c()
      z.best          <- c()
      p.best          <- 0
      k.best          <- c()
      pred.trend.best <- c()
    
      # Loop through the "smoothing parameter" values and select the 
      # one which gives the greatest p-value for WN test (i.e. the one that 
      # best detrends the signal to WN)
      for(k in n.start.smooth:n.stop.smooth){
        
        # The vector that stores the trends for this smoothing parameter
        pred.trend <- c()
      
        # Loop through the plates
        for(i in plates){
          # Subset to this plate 
          plate.dt.met <- met.dt[grep(paste("plate_",i,"_",sep=""), rownames(met.dt)),]
          
          # Get the metabolite measurements in this plate
          y <- plate.dt.met
        
          # Get the trend for this metabolite for this plate and append it to pred.trend
          # ret.trend <- apply.gam.hastie(y=y, k=k)
         
          ret.trend <- tryCatch({
            apply.gam.hastie(y=y, k=k)
          },warning=function(w, z=length(plate.dt.met), k=k, i=i){
            print(paste("WARNING @ apply.gam.hastie > ", met," - smooth=", k," - plate=", i,":", w))
            return(rep(0,z))
          },error = function(e, z=length(plate.dt.met,k=k,i=i)) {
            print(paste("ERROR @ apply.gam.hastie > ", met," - smooth=", k," - plate=", i,":", e))
            return(rep(0,z))
          })
          
          # append to pred.trend
          pred.trend <- c(pred.trend, ret.trend)
        } # closes loop through plates
        
        # Correct the metabolite by subtracting  the trend relative to this 
        # smoothing parameter
        z <- met.dt[[met]] - pred.trend
    
        # Get the white noise test pvalue relative to this smoothing parameter
        p.wn <- test.wn(z)
      
        # Add the WN p-value to the vector of WN test pvalues
        pvals <- c(pvals, p.wn)
        
        # Is this the "best" p-val so far? If so, set z.best, p.best, pred.trend.best
        if(p.wn >= p.best){
          z.best <- z; 
          p.best <- p.wn; 
          k.best <- k; 
          pred.trend.best <- pred.trend
        }
      }  # closes for k in 1:n.smooth
      
      # Define a new "corrected metabolite" variable set to the "z.best" value obtained 
      # after looping through the smothing parameters
      corr.met <- paste("corr.", met, sep="")
      met.dt[[corr.met]] <- z.best
      
      # Add the corrected metabolite to the original dset  
      res.norm.met.dt[[corr.met]] <- met.dt[[corr.met]]
      
      # Now compare the "p.best" value for this metabolite with the 
      # threshold p.value determined in step 3
      result <- "Keep Ho: WN" 
      if(p.best <= p.thr.after.corr){result <- "Reject Ho: not WN";} 
      print(paste(corr.met," :  p.wn=", p.best, " - ", result, sep=""))
      
      # Add the test result to the vector of WN test results 
      is.wn.post.corr <- c(is.wn.post.corr, p.best <= p.thr.after.corr)
      
      # Add the p.best value to the vector of "best" pvalues
      best.pvals <- c(best.pvals,p.best)
      
      ############# 
      # Plot
      #############
      if(save.pdfs){
        file.name <- paste(pdf.dir,"/",met,".pdf",sep="")
        orig.m <- strsplit(met,"res.norm.")[[1]][2]
    
        plot.uncorr.corr(orig.m = orig.m, 
                       unresid.metab=met.dat[[orig.m]], 
                       signal=res.norm.met.dt[met], 
                       corrected.signal=res.norm.met.dt[[corr.met]],
                       seq.order=seq.order,
                       plate.coord = plate.coord,
                       n.plates=length(plates), 
                       file.name=file.name,
                       is.wn=result,
                       p=p.best, 
                       k=k.best)
      } 
    } # Closes for(met in ...)
    
    tmp.dt <- data.frame(corrected.metabs=mets.to.correct, 
                         p.val.after.correction=best.pvals,
                         p.thr.corr=p.thr.after.corr,
                         n.90.pc.corr=n.90.pc.after.corr,
                         is.wn.after.correction=!(is.wn.post.corr),
                         stringsAsFactors = FALSE)
    
    p.wn.tst.before.correction <- data.frame(metab=names(p.wn.tst), 
                                             p.val.before.correction=p.wn.tst, 
                                             p.thr.before.correction=p.thr.before.corr,
                                             n.90.pc.before.corr=n.90.pc.before.corr, 
                                             stringsAsFactors = FALSE)
    
    is.wn.after.correction <- merge(tmp.dt, p.wn.tst.before.correction, by.x = "corrected.metabs", by.y="metab", sort=F)
    
    is.wn.after.correction <- is.wn.after.correction[c("corrected.metabs", "p.val.before.correction", 
                                                       "p.val.after.correction", "is.wn.after.correction",
                                                       "p.thr.before.correction", "n.90.pc.before.corr",
                                                       "p.thr.corr", "n.90.pc.corr")]
    
    is.wn.after.correction <- is.wn.after.correction[order(is.wn.after.correction$p.val.after.correction),]
    
    # Save
    res.norm.met.with.correction.dt <- res.norm.met.dt
    
    if(save.corrected.data){
      save(list=c("res.norm.met.with.correction.dt"), file=paste(here,"/res.norm.met.with.correction.dt.RData",sep=""))
    }
    
    if(save.correction.summary){
      save(list=c("is.wn.after.correction"), file=paste(here,"/is.wn.after.correction.RData",sep=""))
    }
  }else{
    print("No metabolites to correct")
    is.wn.after.correction <- NA
  } # closes in length(met.to.correct > 0
  
  # Return "res.norm.met.dt". If there were metabolites to correct, the dset will include"corr.res.norm."
  # corrected metabolite else only the normalized and residualized.
  return(list(corrected.dset=res.norm.met.dt, correction.summary=is.wn.after.correction))
}  
