library(tidyverse)
library(ribor)
library(psych)
library(wavethresh)
library(rstatix)
library(sgr)
library(parallel)

pkm.ribo <- Ribo("./pkm.ribo", rename = rename_default )
rc_CDS <- get_region_counts(ribo.object    = pkm.ribo,
                            tidy       = TRUE,
                            alias      = TRUE,
                            transcript = FALSE,
                            region     = "CDS",
                            compact    = FALSE)
highexp=rc_CDS%>%filter(count>3000, experiment=="20201209-shRNA-PKM-KD-1-ribo")
transcripts=highexp$transcript
experiments=get_experiments(pkm.ribo)
size=50
ssize=3
threshold=7
#Load all of these functions
  #applies the Bayesian wavelet thresholding function to the
  #coverage data passed to it using the optimal hyperparameters (Daubechie filter 4 (0, 0))
wavetransform=function(cov){
  cov2=cov%>%mutate(position=as.numeric(position))%>%arrange(desc(position))
  max=cov2[[1,2]]
  power=floor(log(max, 2))
  bcov=cov%>%arrange(experiment)
  hecc=c()
  for(i in 1:6){
    great=cov%>%filter(experiment %in% experiments[i])
    joy=great$counts
    great=joy[1:(2^power)]
    great2=joy[(max-(2^power-1)):max]
    blocks.thr1 <- wavethresh::BAYES.THR(great, plotfn=FALSE, filter.number=4, 
                                         family = "DaubExPhase", alpha=0, beta=0)
    blocks.thr2 <- wavethresh::BAYES.THR(great2, plotfn=FALSE, filter.number=4, 
                                         family = "DaubExPhase", alpha=0, beta=0)
    test=c(blocks.thr1, blocks.thr2[(2^(power+1)+1-max):(2^power)])
    hecc=c(hecc, test)
  }
  bcov[,4]=hecc
  return(bcov)
}
  #overlapping Z-score algorithm - does not remove prior peaks and only goes from 5' - 3'
zscoretransform=function(cov, size, threshold, ssize){
  acov=cov%>%mutate(position=as.numeric(position))%>%group_by(experiment)%>%
    mutate(avgcount=mean(count))%>%
    arrange(experiment)%>%mutate(countr=counts+0.0001, z=counts)
  print(acov)
  acov2=acov%>%arrange(desc(position))
  max=acov2[[1,2]]
  #calculates z - leaves out the beginning and end of each transcript
  window=c(1:size)
  subwindow=c(1:ssize)
  stuff=array(dim=c(4, max, 6))
  matrix=array(dim=c(3, 3))
  matrix[,1]=c(1, 0, 0)
  matrix[,2]=c(0, 1, -1)
  matrix[,3]=c(0, 1, 0)
  
  for(i in 0:5){ #experiment number
    print(i)
    for(j in 0:2){ #moving window position
      f=matrix[1, j+1]
      g=matrix[2, j+1]
      h=matrix[3, j+1]
      for(k in (j*(size/2)+f):(max-size+j*(size/2)+f)){#moving window range of nucleotide positions
        m=1
        save=c(1:2*ssize+1)
        counting=1
        for(l in (-1*j*size/2+g):((size-1-j*(size/2))+g+h)){#mini-loop to calculate mean and SD
          window[m]=as.numeric(acov[[i*max+k+l, 6]])
          m=m+1
        }
        if(j==0){
          window=window[-c(1:ssize)]
        }
        else if(j==1){
          window=window[-c((size/2-ssize+2):(size/2+ssize))]
        }
        else{
          window=window[-c((size-ssize+1):(size))]
        }
        mean=mean(window)
        sd=sd(window)
        zscore=(acov[[i*max+k, 7]]-mean)/(sd+0.0001)
        stuff[j+1,k,i+1]=zscore
      }
    }
  }
  acov2=acov%>%arrange(desc(position))
  max=acov2[[i,2]]
  values=c(1:3)
  
  #calculating 4th row of 3D array
  for(j in 1:6){
    for(i in 1+size:max-size){
      for(k in 1:3){
        values[k]=as.numeric(stuff[k, i, j])
      }
      mean=mean(values)
      if(is.na(mean)){
        mean=0
      }
    }
  }
  #setting last column of acov
  for(i in 0:5){
    for(j in (size):(max-size)){
      stuff[4, j, i+1]=(stuff[1, j, i+1]+stuff[2, j, i+1]+stuff[3, j, i+1])/3
      acov[i*max+j, 7]=stuff[4, j, i+1]
    }
  }
  #acov=acov%>%mutate(score=z*counts)
  acov2=acov%>%filter(z>6 & counts>7)
  newlist = list(acov, acov2)
  return (newlist)
}
  #isolates pause sites from the transformed coverage data and puts it in the proper form
  #for heatmaps and ICC/kappa calculation
makeset=function(acov, acov2){
  nrows=nrow(acov2)
  
  acov2=acov2%>%arrange(desc(position))
  max=acov2[[1,2]]
  
  vector3 <- c(1:max)
  for(i in 1:nrows){
    a=acov2[[i, 2]]
    vector3[a]=0
  }
  count=0
  for(i in 1:max){
    if(vector3[i]==0){
      count=count+1
    }
  }
  set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
                 "control1"=1:count, "control2"=1:count, "control3"=1:count)
  #vector4 records all nucleotide sites in the GAPDH transcript with score > 30 in at
  #least one of the experiment
  vector3 <- c(1:max)
  for(i in 1:nrow(acov2)){
    a=acov2[[i, 2]]
    vector3[a]=0
  }
  vector4=c(1:count)
  count=1
  for(i in 1:max){
    if(vector3[i]==0){
      vector4[count]=i
      count=count+1
    }
  }
  #constructing the skeleton for set
  count=count-1
  for(i in 1:count){
    set[i, 1]=vector4[i]
  }
  
  for(i in 1:count){
    for(j in 2:7){
      set[i,j]=0
    }
  }
  
  acov3=acov%>%filter(position%in%vector4)
  count=1
  acov3=acov3%>%arrange(position)
  
  rows=nrow(set)
  #fills in the empty score cells
  for(i in 1:rows){
    for(j in 2:7){
      set[i, j]=acov3[(i-1)*6+j-1, ncol(acov3)]
    }
  }
  #compute averages across all experiments and the ratio of these averages
  set=set%>%mutate(avgKD=(KD1+KD2+KD3)/3, avgcontrol=(control1+control2+control3)/3)%>%
    mutate(ratio=log(avgKD/avgcontrol, 2))%>%filter(!(is.na(position)))
  set2 <- set[,-1]
  rownames(set2) <- set[,1]
  return(set)
}
  #same as makeset, but replaces the scores with 1 for pause site and 0 for not pause site
makekappaset=function(acov, acov2, threshold){
  nrows=nrow(acov2)
  
  acov2=acov2%>%arrange(desc(position))
  max=acov2[[1,2]]
  
  vector3 <- c(1:max)
  count=0
  for(i in 1:max){
    if(vector3[i]==0){
      count=count+1
    }
  }
  vector4=c(1:count)
  for(i in 1:nrows){
    a=acov2[[i, 2]]
    vector3[a]=0
  }
  count=1
  for(i in 1:max){
    if(vector3[i]==0){
      vector4[count]=i
      count=count+1
    }
  }
  count=count-1
  set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
                 "control1"=1:count, "control2"=1:count, "control3"=1:count)
  #vector4 records all nucleotide sites in the GAPDH transcript with score > 30 in at
  #least one of the experiments
  #constructing the skeleton for set
  for(i in 1:count){
    set[i, 1]=vector4[i]
  }
  
  for(i in 1:count){
    for(j in 2:7){
      set[i,j]=0
    }
  }
  
  acov3=acov%>%filter(position%in%vector4)
  count=1
  acov3=acov3%>%arrange(position)
  
  rows=nrow(set)
  #fills in the empty score cells
  for(i in 1:rows){
    for(j in 2:7){
      if(acov3[[(i-1)*6+j-1, ncol(acov3)]]>threshold){
        set[i, j]=1
      }
    }
  }
  set2 <- set[,-1]
  rownames(set2) <- set[,1]
  return(set)
}
get_icc=function(set){
  subset=set[,c(2:4)]
  subset2=set[,c(5:7)]
  le=ICC(subset)
  le2=ICC(subset2)
  le=le[[1]]
  le2=le2[[1]]
  meanicc=(le[2,2]+le2[2,2])/2
  return(meanicc)
}
get_kappa=function(set){
  subset=set[,c(2:4)]
  subset2=set[,c(5:7)]
  le=cohen.kappa(subset)
  le2=cohen.kappa(subset2)
  le=le[[5]]
  le2=le2[[5]]
  meanicc=(le+le2)/2
  return(meanicc)
}
  #Stack Overflow algorithm with influence = 0. Removes prior peaks. 
laggedtransform=function(cov, window, threshold, subwindow){
  acov=cov%>%mutate(position=as.numeric(position))%>%group_by(experiment)%>%
    arrange(experiment)%>%mutate(scan1=0, scan2=0, scan3=0)
  max=acov%>%arrange(desc(position))
  max=max[[1,2]]
  ohboy=acov%>%arrange(desc(position))
  for(i in 1:6){
    print(i)
    values=acov$counts
    max=ohboy[[1,2]]
    max2=ohboy[[1,2]]
    j3=window
    for(j in window:max){
      position=(i-1)*max2+j3
      position2=(i-1)*max2+j
      mean=mean(values[(position-window+1):(position-4)])
      sd=sd(values[(position-window+1):(position-4)])
      acov[[position2, 5]]=(values[position]-mean)/sd
      if(is.na(acov[[position2, 5]] | is.infinite(acov[[position2, 5]]))){
        acov[[position2, 5]]=0
      }
      if(acov[[position2, 5]]>threshold & values[position]>7){
        print(j)
        values=values[-position]
        max=max-1
        j3=j3-1
      }
      j3=j3+1
    }
  }
  
  for(i in 1:6){
    values=acov$counts
    max=ohboy[[1,2]]
    max2=ohboy[[1,2]]
    print(i)
    for(j in (max-window+1):1){
      position=(i-1)*max2+j
      mean=mean(values[(position+4):(position+3+window)])
      sd=sd(values[(position+4):(position+3+window)])
      acov[[position, 6]]=(values[position]-mean)/sd
      if(is.na(acov[[position, 6]] | is.infinite(acov[[position, 6]]))){
        acov[[position, 6]]=0
      }
      if(acov[[position, 6]]>threshold & values[position]>threshold){
        print(j)
        values=values[-position]
      }
    }
  }
  for(i in 1:6){
    values=acov$counts
    max=ohboy[[1,2]]
    print(i)
    for(j in (max-window/2+1):(window/2)){
      position=(i-1)*max+j
      mean=mean(values[c((position-window/2+1):(position-3), (position+3):(position+window/2-1))])
      sd=sd(values[c((position-window/2+1):(position-3), (position+3):(position+window/2-1))])
      acov[[position, 7]]=(values[position]-mean)/sd
      if(is.na(acov[[position, 7]]) | is.infinite(acov[[position, 7]])){
        acov[[position, 7]]=0
      }
    }
  }
  acov=acov%>%mutate(z=(scan1+scan2+scan3)/3)
  acov2=acov%>%filter(z>threshold & counts>7)
  acov=acov[,c(1:6, 8, 7)]
  return(list(acov, acov2))
}
  #if two peaks are within a subwindow of each other, this function will remove the peak
  #that is smaller in at least 5 of the 6 conditions
removepeaks=function(acov2, subwindow){
  acov2=acov2%>%mutate(remove=0)
  ncols=ncol(acov2)
  for(i in 1:nrow(acov2)){
    if(!(i==1) & acov2[[i, ncols]]==0){
      if((acov2[[i,2]]-acov2[[i-1, 2]])<=subwindow & (acov2[[i,2]]-acov2[[i-1, 2]])>0){
        if(acov2[[i-1,3]]>acov2[[i,3]]){
          acov2[[i, ncols]]=1
        }
        else{
          acov2[[i, ncols]]=-4
        }
      }
    }
    if(!(i==nrow(acov2)) & acov2[[i, ncols]]==0){
      if((acov2[[i,2]]-acov2[[i+1, 2]])>=-subwindow & (acov2[[i,2]]-acov2[[i+1, 2]])<0){
        if(acov2[[i+1,3]]>acov2[[i,3]]){
          acov2[[i, ncols]]=1
        }
        else{
          acov2[[i, ncols]]=-6
        }
      }
    }
  }
  acov2=acov2%>%group_by(position)%>%mutate(sum=sum(remove))
  acov2=acov2%>%filter(sum<=0)
  return(acov2)
}
getzscorestats_o=function(transcript, size, threshold, ssize){
  cov <- get_coverage(ribo.object = pkm.ribo,
                      name        = transcript,
                      range.lower = 15,
                      range.upper = 40,
                      length      = TRUE,
                      alias       = TRUE,
                      compact = FALSE,
                      tidy=TRUE,
                      experiment  = get_experiments(pkm.ribo))
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=zscoretransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2)
  kappaset=makekappaset(acov, acov2, threshold)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  newlist=(c(get_icc(set), get_kappa(kappaset)))
}
getzscorestats_l=function(transcript, size, threshold, ssize){
  cov <- get_coverage(ribo.object = pkm.ribo,
                      name        = transcript,
                      range.lower = 15,
                      range.upper = 40,
                      length      = TRUE,
                      alias       = TRUE,
                      compact = FALSE,
                      tidy=TRUE,
                      experiment  = get_experiments(pkm.ribo))
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=laggedtransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2)
  kappaset=makekappaset(acov, acov2, threshold)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  newlist=(c(get_icc(set), get_kappa(kappaset)))
}
  #for wavelet thresholded data, scores at each nucleotide position are calculated by dividing
  #the count at that region by the global count average
getwavestats=function(transcript, threshold){
  cov <- get_coverage(ribo.object = pkm.ribo,
                      name        = transcript,
                      range.lower = 15,
                      range.upper = 40,
                      length      = TRUE,
                      alias       = TRUE,
                      compact = FALSE,
                      tidy=TRUE,
                      experiment  = get_experiments(pkm.ribo))
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  cov=wavetransform(cov)
  acov=cov%>%mutate(position=as.numeric(position))%>%group_by(experiment)%>%
    mutate(mean=counts/mean(counts))
  acov2=acov%>%filter(mean>7)
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2)
  kappaset=makekappaset(acov, acov2, threshold)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  get_icc(set)
  newlist=(c(get_icc(set), get_kappa(kappaset)))
  return (newlist)
}

#Code to run -- each function will return a 2 by 186 list with the ICC and kappa values for each
#algorithm on the 186 highest expressed transcripts in pkm.ribo. If the code throws any errors,
#please let me know and i can take a look
wavestats=mcmapply(function(x,y){
  return(getwavestats(x,y))
}, x = transcripts, y = threshold, mc.cores=48)

zscorestats_o=mcmapply(function(x,y,z,w){
  return(getzscorestats_o(x,y,z,w))
}, x = transcripts, y = size, z=threshold, w=ssize, 
mc.cores=48)

zscorestats_l=mcmapply(function(x,y,z,w){
  return(getzscorestats_l(x,y,z,w))
}, x = transcripts, y = size, z=threshold, w=ssize,
mc.cores=48)


