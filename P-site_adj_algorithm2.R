library(tidyverse)
library(ribor)
library(psych)
library(wavethresh)
library(rstatix)
library(sgr)
library(parallel)

## CC added library reshape; edger
library(reshape2)
library(edgeR)

pkm.ribo <- Ribo("./pkm.ribo", rename = rename_default )
# CC -> Added range lower/upper for consistency
rc_CDS <- get_region_counts(ribo.object    = pkm.ribo,
                            range.lower = 29,
                            range.upper = 36,
                            tidy       = TRUE,
                            alias      = TRUE,
                            transcript = FALSE,
                            normalize=TRUE,
                            region     = "CDS",
                            compact    = FALSE)
#retrieves the 350 most expressed transcripts
# CC removed deprecated function; updated selection of highexpression genes to include all experiments
region_lengths <- get_internal_region_lengths(ribo.object = pkm.ribo, alias = TRUE)
cds=region_lengths%>%select(transcript, CDS)
rc_CDS=rc_CDS%>%left_join(cds)%>%mutate(count=count/CDS)

rc_CDS_w = dcast(rc_CDS[,-5], transcript ~ experiment)
# A value of 500 gives 354 transcripts; 1500 gives 106. 
high_exp = rowSums( cpm(rc_CDS_w[,-1]) > 500 ) > 2
sum(high_exp)

#highexp=rc_CDS%>%filter(count>0.4125, experiment=="20201209-shRNA-PKM-KD-1-ribo")
transcripts=rc_CDS_w$transcript[high_exp]
experiments=get_experiments(pkm.ribo)

#Load all of these functions
#applies the Bayesian wavelet thresholding function to the
#coverage data passed to it using the optimal hyperparameters (Daubechie filter 4 (0, 0))
wavetransform=function(cov){
  cov2=cov%>%mutate(position=as.numeric(position))%>%arrange(desc(position))
  max=cov2[[1,2]]
  power=floor(log(max, 2))
  bcov=cov%>%arrange(experiment)
  hecc=c()
  noise=c(rep(0.00001,nrow(cov)/6))
  noise=jitter(noise, factor = 1, amount=NULL)
  cov=cov%>%mutate(noises=noise)
  for(i in 1:6){
    great=cov%>%filter(experiment %in% experiments[i])%>%mutate(counts=counts+noises)
    bruh=great$counts
    great=bruh[1:(2^power)]
    great2=bruh[(max-(2^power-1)):max]
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
makeset=function(acov, acov2, threshold){
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
    mutate(ratio=avgcontrol)%>%filter(!(is.na(position)))
  set2 <- set[,-1]
  rownames(set2) <- set[,1]
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))%>%mutate(sig=0)
  m=1
  for(i in 1:nrow(set)){
    if(set[[m,8]]>threshold | set[[m,9]]>threshold){
      set[[m, 11]]=1
    }
    if(set[[m,11]]==0){
      set=set[-m,]
      m=m-1
    }
    m=m+1
  }
  set=set%>%select(-sig)
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
makekappaset2=function(set){
  kappaset=set
  for(i in 1:nrow(set)){
    for(j in 2:7){
      if(set[[i,j]]>threshold){
        kappaset[[i,j]]=1
      }
      else{
        kappaset[[i,j]]=0
      }
    }
  }
  return(kappaset)
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
  meankappa=(le+le2)/2
  return(meankappa)
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
                      range.lower = 29,
                      range.upper = 36,
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
  set=makeset(acov, acov2, threshold)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  return(get_icc(set))
}
getzscorestats_l=function(transcript, size, threshold, ssize){
  cov <- get_coverage(ribo.object = pkm.ribo,
                      name        = transcript,
                      range.lower = 29,
                      range.upper = 36,
                      length      = TRUE,
                      alias       = TRUE,
                      compact = FALSE,
                      tidy=TRUE,
                      experiment  = get_experiments(pkm.ribo))
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=laggedtransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  if(nrow(acov2)==0){
    return(0)
  }
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2, threshold)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  return(get_icc(set))
}
wavelettransform=function(cov, threshold){
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  cov=wavetransform(cov)
  acov=cov%>%mutate(position=as.numeric(position))%>%group_by(experiment)%>%
    mutate(mean=counts/mean(counts))
  acov2=acov%>%filter(mean>threshold)
  if(nrow(acov2)==0){
    return(0)
  }
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2, threshold)
  if(nrow(set)==0){
    return (0)
  }
  kappaset=makekappaset2(set)
  return(list(acov, acov2, set, kappaset))
}

#for wavelet thresholded data, scores at each nucleotide position are calculated by dividing
#the count at that region by the global count average
getwavestats=function(transcript, threshold){
  cov <- get_coverage(ribo.object = pkm.ribo,
                      name        = transcript,
                      range.lower = 29,
                      range.upper = 36,
                      length      = TRUE,
                      alias       = TRUE,
                      compact = FALSE,
                      tidy=TRUE,
                      experiment  = get_experiments(pkm.ribo))
  ok=wavelettransform(cov, threshold)
  if(is.double(ok)){
    return(0)
  }
  set=ok[[3]]
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  if(nrow(set)>2){
    return(get_icc(set))
  }else{
    return(0)
  }
}

#applies p site adjustment to the coverage plot; shifts rl 29 and 30 backwards
getpcov=function(transcript){
  experiments=get_experiments(pkm.ribo)
  cov <- get_coverage(ribo.object = pkm.ribo,
                      name        = transcript,
                      range.lower = 29,
                      range.upper = 36,
                      length      = FALSE,
                      alias       = TRUE,
                      compact = FALSE,
                      tidy=TRUE,
                      experiment  = experiments[1:6])
  cov=cov%>%mutate(position=as.numeric(position))
  cov2=cov%>%filter(length%in%c(29:30))
  cov3=cov%>%filter(length%in%c(31:36))
  cov2=cov2%>%mutate(position=position-1)
  bind=cov2[c((nrow(cov2)-11):nrow(cov2)), ]
  bind=bind%>%mutate(position=position+1)
  cov2=cov2[c(13:nrow(cov2)),]
  acov=rbind(cov2, bind)
  cov=rbind(acov, cov3)
  cov=cov%>%group_by(experiment, position)%>%mutate(pcounts=sum(count))%>%filter(length==29)
  cov=cov%>%mutate(count=pcounts)%>%select(-length, -pcounts)
  return(cov)
}

wavestats_p=function(transcript, threshold){
  cov=getpcov(transcript)
  ok=wavelettransform(cov, threshold)
  if(is.double(ok)){
    return(0)
  }
  set=ok[[3]]
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  if(nrow(set)>2){
    return(c(get_icc(set)))
  }else{
    return(0)
  }
}
zscorestats_op=function(transcript, size, threshold, ssize){
  cov=getpcov(transcript)
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=zscoretransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  if(nrow(acov2)==0){
    return(0)
  }
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2, threshold)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  return(get_icc(set))
}
zscorestats_lp=function(transcript, size, threshold, ssize){
  cov=getpcov(transcripts)
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=laggedtransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  if(nrow(acov2)==0){
    return(0)
  }
  acov2=removepeaks(acov2, ssize)
  set=makeset(acov, acov2, threshold)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  return(get_icc(set))
}
differentialpeaks=function(set){
  set=set%>%mutate(statistic=ratio, p=ratio, significance=ratio)
  for(i in 1:nrow(set)){
    ouch=data.frame("id"=c(1, 2, 3, 1, 2, 3), 
                    "group"=c("KD", "KD", "KD", "control", "control", "control"), 
                    "score"=1:6)
    for(j in 1:6){
      ouch[j, 3]=set[i, j+1]
    }
    quack=ouch%>%t_test(score ~ group)%>%add_significance()
    set[i, 11]=quack[1, 6]
    set[i, 12]=quack[1, 8]
    set[i, 13]=quack[1, 9]
  }
  set=set[,c(1, 13, 11, 12, 8, 9, 2:7, 10)]
  gapdh=set%>%select(position, statistic, p, avgcontrol, avgKD)%>%filter(p<0.05)
  newlist=list(set, gapdh)
  return(newlist)
}

wavepeakstats=function(transcript, threshold){
  cov=getpcov(transcript)
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp2=wavelettransform(cov, threshold)
  if(is.double(temp2)){
    return(0)
  }
  set=temp2[[3]]
  positions=set$position
  count=length(positions)
  pcov=cov%>%filter(position %in% positions)
  set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
                 "control1"=1:count, "control2"=1:count, "control3"=1:count)
  set[,1]=positions
  for(i in 1:count){
    for(j in 2:7){
      set[i, j]=pcov[(i-1)*6+j-1, ncol(pcov)-1]
    }
  }
  set=set%>%mutate(avgKD=(KD1+KD2+KD3)/3, avgcontrol=(control1+control2+control3)/3)%>%
    mutate(ratio=log(avgKD/avgcontrol, 2))%>%filter(!(is.na(position)))%>%mutate(sig=0)
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  m=1
  for(i in 1:nrow(set)){
    if(set[[m,8]]>5 | set[[m,9]]>5){
      set[[m, 11]]=1
    }
    if(set[[m,11]]==0){
      set=set[-m,]
      m=m-1
    }
    m=m+1
  }
  set=set%>%select(-sig)
  set2 <- set[,-1]
  rownames(set2) <- set[,1]
  positions=set$position
  count=length(positions)
  pcov=cov%>%filter(position %in% positions)
  set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
                 "control1"=1:count, "control2"=1:count, "control3"=1:count)
  set[,1]=positions
  for(i in 1:count){
    for(j in 2:7){
      set[i, j]=pcov[(i-1)*6+j-1, ncol(pcov)]
    }
  }
  set=set%>%mutate(avgKD=(KD1+KD2+KD3)/3, avgcontrol=(control1+control2+control3)/3)%>%
    mutate(ratio=log(avgKD/avgcontrol, 2))%>%filter(!(is.na(position)))
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  kappaset=set
  ah=differentialpeaks(set)
  diffpeaks=ah[[2]]
  diffpositions=diffpeaks[,c(1:3)]
  set=set%>%mutate(sig=0)
  for(i in 1:nrow(set)){
    for(j in 2:7){
      if(set[[i,j]]>threshold){
        kappaset[[i,j]]=1
      }
      else{
        kappaset[[i,j]]=0
      }
    }
  }
  allpeaks=c(rep(0,9))
  diffpeaks=c(rep(0,9))
  kappaset=kappaset%>%select(-avgKD, -avgcontrol, -ratio)
  kappaset=kappaset%>%mutate(kd=0, control=0)
  
  if(!(nrow(kappaset)==0)){
    for(i in 1:nrow(kappaset)){
      for(j in 1:2){
        for(k in 1:3){
          if(kappaset[[i, (j-1)*3+k+1]]==1){
            kappaset[[i, 7+j]]=kappaset[[i, 7+j]]+1
          }
        }
      }
    }
    for(i in 1:nrow(kappaset)){
      if((kappaset[[i, 8]]+kappaset[[i, 9]])==1){
        allpeaks[[1]]=allpeaks[[1]]+1
      }
      else if((kappaset[[i, 8]]==1 & kappaset[[i, 9]]==1)){
        allpeaks[[2]]=allpeaks[[2]]+1
      }
      else if((kappaset[[i, 8]]==2 & kappaset[[i, 9]]==0) | (kappaset[[i, 8]]==0 & kappaset[[i, 9]]==2)){
        allpeaks[[3]]=allpeaks[[3]]+1
      }
      else if((kappaset[[i, 8]]==1 & kappaset[[i, 9]]==2) | (kappaset[[i, 8]]==2 & kappaset[[i, 9]]==1)){
        allpeaks[[4]]=allpeaks[[4]]+1
      }
      else if((kappaset[[i, 8]]==3 & kappaset[[i, 9]]==0) | (kappaset[[i, 8]]==0 & kappaset[[i, 9]]==3)){
        allpeaks[[5]]=allpeaks[[5]]+1
      }
      else if((kappaset[[i, 8]]==2 & kappaset[[i, 9]]==2)){
        allpeaks[[6]]=allpeaks[[6]]+1
      }
      else if((kappaset[[i, 8]]==3 & kappaset[[i, 9]]==1) | (kappaset[[i, 8]]==1 & kappaset[[i, 9]]==3)){
        allpeaks[[7]]=allpeaks[[7]]+1
      }
      else if((kappaset[[i, 8]]==3 & kappaset[[i, 9]]==2) | (kappaset[[i, 8]]==2 & kappaset[[i, 9]]==3)){
        allpeaks[[8]]=allpeaks[[8]]+1
      }
      else{
        allpeaks[[9]]=allpeaks[[9]]+1
      }
    }
  }
  
  ahh=diffpositions
  diffpositions=diffpositions$position
  diffkappaset=kappaset%>%filter(position %in% diffpositions)
  diffkappaset=diffkappaset%>%mutate(kd=0, control=0)
  
  if(!(length(diffpositions)==0)){
    for(i in 1:nrow(diffkappaset)){
      for(j in 1:2){
        for(k in 1:3){
          if(diffkappaset[[i, (j-1)*3+k+1]]==1){
            diffkappaset[[i, 7+j]]=diffkappaset[[i, 7+j]]+1
          }
        }
      }
    }
    for(i in 1:nrow(diffkappaset)){
      if((diffkappaset[[i, 8]]+diffkappaset[[i, 9]])==1){
        diffpeaks[[1]]=diffpeaks[[1]]+1
      }
      else if((diffkappaset[[i, 8]]==1 & diffkappaset[[i, 9]]==1)){
        diffpeaks[[2]]=diffpeaks[[2]]+1
      }
      else if((diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==0) | (diffkappaset[[i, 8]]==0 & diffkappaset[[i, 9]]==2)){
        diffpeaks[[3]]=diffpeaks[[3]]+1
      }
      else if((diffkappaset[[i, 8]]==1 & diffkappaset[[i, 9]]==2) | (diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==1)){
        diffpeaks[[4]]=diffpeaks[[4]]+1
      }
      else if((diffkappaset[[i, 8]]==3 & diffkappaset[[i, 9]]==0) | (diffkappaset[[i, 8]]==0 & diffkappaset[[i, 9]]==3)){
        diffpeaks[[5]]=diffpeaks[[5]]+1
      }
      else if((diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==2)){
        diffpeaks[[6]]=diffpeaks[[6]]+1
      }
      else if((diffkappaset[[i, 8]]==3 & diffkappaset[[i, 9]]==1) | (diffkappaset[[i, 8]]==1 & diffkappaset[[i, 9]]==3)){
        diffpeaks[[7]]=diffpeaks[[7]]+1
      }
      else if((diffkappaset[[i, 8]]==3 & diffkappaset[[i, 9]]==2) | (diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==3)){
        diffpeaks[[8]]=diffpeaks[[8]]+1
      }
      else{
        diffpeaks[[9]]=diffpeaks[[9]]+1
      }
    }
  }
  allpos=kappaset$position
  numpeaks=c(nrow(set), nrow(diffkappaset))
  newlist=list(allpeaks, diffpeaks, numpeaks, allpos, ahh)
  return(newlist)
}
zscorepeakstats=function(transcript, size, threshold, ssize){
  cov=getpcov(transcript)
  cov=cov%>%group_by(experiment)%>%mutate(counts=count/mean(count))
  temp=laggedtransform(cov, size, threshold, ssize)
  acov=temp[[1]]
  acov2=temp[[2]]
  if(nrow(acov2)==0){
    return(0)
  }
  set=makeset(acov, acov2, threshold)
  positions=set$position
  count=length(positions)
  pcov=cov%>%filter(position %in% positions)
  set=data.frame("position"=1:count, "KD1"=1:count, "KD2"=1:count, "KD3"=1:count, 
                 "control1"=1:count, "control2"=1:count, "control3"=1:count)
  set[,1]=positions
  for(i in 1:count){
    for(j in 2:7){
      set[i, j]=pcov[(i-1)*6+j-1, ncol(pcov)]
    }
  }
  set=set%>%mutate(avgKD=(KD1+KD2+KD3)/3, avgcontrol=(control1+control2+control3)/3)%>%
    mutate(ratio=log(avgKD/avgcontrol, 2))%>%filter(!(is.na(position)))
  set=set%>%filter(!is.na(ratio) & !is.infinite(ratio))
  set2 <- set[,-1]
  rownames(set2) <- set[,1]
  kappaset=set
  ah=differentialpeaks(set)
  diffpeaks=ah[[2]]
  diffpositions=diffpeaks[,c(1:3)]
  for(i in 1:nrow(set)){
    for(j in 2:7){
      if(set[[i,j]]>threshold){
        kappaset[[i,j]]=1
      }
      else{
        kappaset[[i,j]]=0
      }
    }
  }
  allpeaks=c(rep(0,9))
  diffpeaks=c(rep(0,9))
  kappaset=kappaset%>%select(-avgKD, -avgcontrol, -ratio)
  kappaset=kappaset%>%mutate(kd=0, control=0)
  
  if(!(nrow(kappaset)==0)){
    for(i in 1:nrow(kappaset)){
      for(j in 1:2){
        for(k in 1:3){
          if(kappaset[[i, (j-1)*3+k+1]]==1){
            kappaset[[i, 7+j]]=kappaset[[i, 7+j]]+1
          }
        }
      }
    }
    for(i in 1:nrow(kappaset)){
      if((kappaset[[i, 8]]+kappaset[[i, 9]])==1){
        allpeaks[[1]]=allpeaks[[1]]+1
      }
      else if((kappaset[[i, 8]]==1 & kappaset[[i, 9]]==1)){
        allpeaks[[2]]=allpeaks[[2]]+1
      }
      else if((kappaset[[i, 8]]==2 & kappaset[[i, 9]]==0) | (kappaset[[i, 8]]==0 & kappaset[[i, 9]]==2)){
        allpeaks[[3]]=allpeaks[[3]]+1
      }
      else if((kappaset[[i, 8]]==1 & kappaset[[i, 9]]==2) | (kappaset[[i, 8]]==2 & kappaset[[i, 9]]==1)){
        allpeaks[[4]]=allpeaks[[4]]+1
      }
      else if((kappaset[[i, 8]]==3 & kappaset[[i, 9]]==0) | (kappaset[[i, 8]]==0 & kappaset[[i, 9]]==3)){
        allpeaks[[5]]=allpeaks[[5]]+1
      }
      else if((kappaset[[i, 8]]==2 & kappaset[[i, 9]]==2)){
        allpeaks[[6]]=allpeaks[[6]]+1
      }
      else if((kappaset[[i, 8]]==3 & kappaset[[i, 9]]==1) | (kappaset[[i, 8]]==1 & kappaset[[i, 9]]==3)){
        allpeaks[[7]]=allpeaks[[7]]+1
      }
      else if((kappaset[[i, 8]]==3 & kappaset[[i, 9]]==2) | (kappaset[[i, 8]]==2 & kappaset[[i, 9]]==3)){
        allpeaks[[8]]=allpeaks[[8]]+1
      }
      else{
        allpeaks[[9]]=allpeaks[[9]]+1
      }
    }
  }
  
  ahh=diffpositions
  diffpositions=diffpositions$position
  diffkappaset=kappaset%>%filter(position %in% diffpositions)
  diffkappaset=diffkappaset%>%mutate(kd=0, control=0)
  
  if(!(length(diffpositions)==0)){
    for(i in 1:nrow(diffkappaset)){
      for(j in 1:2){
        for(k in 1:3){
          if(diffkappaset[[i, (j-1)*3+k+1]]==1){
            diffkappaset[[i, 7+j]]=diffkappaset[[i, 7+j]]+1
          }
        }
      }
    }
    for(i in 1:nrow(diffkappaset)){
      if((diffkappaset[[i, 8]]+diffkappaset[[i, 9]])==1){
        diffpeaks[[1]]=diffpeaks[[1]]+1
      }
      else if((diffkappaset[[i, 8]]==1 & diffkappaset[[i, 9]]==1)){
        diffpeaks[[2]]=diffpeaks[[2]]+1
      }
      else if((diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==0) | (diffkappaset[[i, 8]]==0 & diffkappaset[[i, 9]]==2)){
        diffpeaks[[3]]=diffpeaks[[3]]+1
      }
      else if((diffkappaset[[i, 8]]==1 & diffkappaset[[i, 9]]==2) | (diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==1)){
        diffpeaks[[4]]=diffpeaks[[4]]+1
      }
      else if((diffkappaset[[i, 8]]==3 & diffkappaset[[i, 9]]==0) | (diffkappaset[[i, 8]]==0 & diffkappaset[[i, 9]]==3)){
        diffpeaks[[5]]=diffpeaks[[5]]+1
      }
      else if((diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==2)){
        diffpeaks[[6]]=diffpeaks[[6]]+1
      }
      else if((diffkappaset[[i, 8]]==3 & diffkappaset[[i, 9]]==1) | (diffkappaset[[i, 8]]==1 & diffkappaset[[i, 9]]==3)){
        diffpeaks[[7]]=diffpeaks[[7]]+1
      }
      else if((diffkappaset[[i, 8]]==3 & diffkappaset[[i, 9]]==2) | (diffkappaset[[i, 8]]==2 & diffkappaset[[i, 9]]==3)){
        diffpeaks[[8]]=diffpeaks[[8]]+1
      }
      else{
        diffpeaks[[9]]=diffpeaks[[9]]+1
      }
    }
  }
  allpos=kappaset$position
  numpeaks=c(nrow(set), nrow(diffkappaset))
  newlist=list(allpeaks, diffpeaks, numpeaks, allpos, ahh)
  return(newlist)
}
#Code to run -- each function will return a 2 by 300 list with the ICC values for each
#algorithm on the 300 highest expressed transcripts in pkm.ribo. If the code throws any errors,
#please let me know and i can take a look. if it takes too long, no need to run the overlapping algo.

threshold=9
size=50
ssize=3

## CC ->  This had an error in one of the cores. This might be an edge-case. 
## It's encountered when using the top 350 transcripts. Ran fine with the top 283. 
## We might want to add some error handling to the code to isolate the problematic txn.  
wavestats=mcmapply(function(x,y){
  return(getwavestats(x,y))
}, x = transcripts, y = threshold, mc.cores=48)

# zscorestats_o=mcmapply(function(x,y,z,w){
#   return(getzscorestats_o(x,y,z,w))
# }, x = transcripts, y = size, z=threshold, w=ssize, 
# mc.cores=48)


## This ran fine with top 106; Encountered_issues with 283
zscorestats_l=mcmapply(function(x,y,z,w){
  return(getzscorestats_l(x,y,z,w))
}, x = transcripts, y = size, z=threshold, w=ssize,
mc.cores=48)

#new code for the P-site adjusted coverage plots
pwavestats=mcmapply(function(x,y){
  return(wavestats_p(x,y))
}, x = transcripts, y = threshold,
mc.cores=48)

# zscorestats_op=mcmapply(function(x,y,z,w){
#   return(zscorestats_op(x,y,z,w))
# }, x = transcripts, y = size, z=threshold, w=ssize,
# mc.cores=48)

zscorestats_lp=mcmapply(function(x,y,z,w){
  return(zscorestats_lp(x,y,z,w))
}, x = transcripts, y = size, z=threshold, w=ssize,
mc.cores=48)


# the peakstats functions return a list with five elements:
#1: a vector containing 9 elements, each one counting the # of pause sites detected by the algorithm in
#   a particular combination of trials. For example, the first element represents the # of pause sites that are only
#   detected in one trial, while the third element represents the # of pause sites detected in 2 trials for one condition
#   and 0 in the other. 
#2: The same as 1, but it only includes the distribution for differential pause sites.
#3: Total # of pause sites detected by the algorithm
#4: Total # of differential pause sites
#5: All nucleotide positions detected as peaks
#6: All nucleotide positions detected as differential peaks

wavepeakcounts=mcmapply(function(x,y){
  return(wavepeakstats(x,y))
}, x = transcripts, y = threshold, mc.cores=48)

zscorepeakcounts=mcmapply(function(x,y,z,w){
  return(zscorepeakstats(x,y,z,w))
}, x = transcripts, y = size, z=threshold, w=ssize, mc.cores=48)

