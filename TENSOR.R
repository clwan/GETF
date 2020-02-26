### simulate data 4 basis 0.4
DIM_use<-100

U1<-matrix(rbinom(5*DIM_use,1,0.4),ncol = 5,nrow = DIM_use)
U2<-matrix(rbinom(5*(DIM_use-10),1,0.4),ncol = 5,nrow = DIM_use-10)
U3<-matrix(rbinom(5*(DIM_use-20),1,0.4),ncol = 5,nrow = DIM_use-20)
U4<-matrix(rbinom(5*(DIM_use-30),1,0.4),ncol = 5,nrow = DIM_use-30)

TENS<-array(0,dim=c(nrow(U1),nrow(U2),nrow(U3),nrow(U4)))
for (i in 1:ncol(U1)) {
  TENS<-TENS+U1[,i]%o%U2[,i]%o%U3[,i]%o%U4[,i]
}

TENS<-array(0,dim=c(nrow(U1),nrow(U2),nrow(U3)))
for (i in 1:ncol(U1)) {
  TENS<-TENS+U1[,i]%o%U2[,i]%o%U3[,i]
}

TENS<-array(as.numeric(TENS>0),dim=dim(TENS))

POS<-function(VEC,n){
  temp<-which(VEC>0)
  if(length(temp)<n){
    return(0)
  }else{
    return(temp[order(VEC[temp],decreasing = T)[length(temp)%/%n]])
  }
}

FOLD<-function(TENS,n){
  if(is.null(dim(TENS))){
    return(POS(TENS,n))
  }else{
    TENS1<-apply(TENS,2:length(dim(TENS)),sum)
    test<-FOLD(TENS1,n+1)
    temp<-NULL
    for (i in 1:length(test)) {
      temp<-paste(temp,",test[",i,"]",sep = "")
    }
    C<-eval(parse(text = paste("POS(TENS[",temp,"],n)",sep = "")))
    return(c(C,test))
  }
}

BASIS_FIND<-function(TENS){
  DIM<-dim(TENS)
  LIST<-list()
  for (i in 1:length(DIM)) {
    TEMP<-apply(TENS, c(1:length(DIM))[-i], sum)
    temp<-rep(NA,length(DIM))
    temp[-i]<-FOLD(TEMP,2)
    LIST[[i]]<-temp
  }
  names(LIST)<-1:length(DIM)
  return(LIST)
}


BASIS_convert<-function(VEC){
  if(is.na(VEC[1])){
    temp<-""
  }else{
    temp<-VEC[1]
  }
  for (i in 2:length(VEC)) {
    if(is.na(VEC[i])){
      temp<-paste(temp,",",sep = "")
    }else{
      temp<-paste(temp,VEC[i],sep = ",")
    }
  }
  return(temp)
}


MEASURE<-function(TENS,BASIS,Thres){
  DIM<-dim(TENS)
  DIR<-sum(1:length(DIM))-sum(as.numeric(names(BASIS)))
  
  if(DIR==1){
    temp<-""
  }else{
    temp1<-BASIS_convert(BASIS[[1]])
    temp<-paste('which(TENS[',temp1,']==1)',sep="")
  }
  
  for (i in 2:length(DIM)) {
    if(i==DIR){
      temp<-paste(temp,",",sep = "")
    }else{
      temp1<-BASIS_convert(BASIS[[as.character(i)]])
      temp<-paste(temp,',which(TENS[',temp1,']==1)',sep="")
    }
  }
  
  VALUE<-eval(parse(text = paste("apply(TENS[",temp,"],",DIR,",sum)",sep = "")))
  
  NEW_BASIS<-rep(0,DIM[DIR])
  NEW_BASIS[which(VALUE>(Thres*max(VALUE)))]<-1
  
  DIM2<-eval(parse(text = paste("dim(TENS[",temp,"])",sep = "")))
  VAL<-1
  for (i in 1:length(DIM2)) {
    VAL=VAL*DIM2[i]
  }
  VAL<-VAL/DIM[DIR]*sum(NEW_BASIS)
  C_value<-sum(VALUE[as.logical(NEW_BASIS)])-VAL*Thres**length(DIM)
  if(sum(NEW_BASIS)<2){
    return(NULL)
  }else{
    return(list(C_value,NEW_BASIS))
  }
}


GEO_SEG<-function(TENS,Thres){
  BASIS<-BASIS_FIND(TENS)
  DIR<-0
  C<--Inf
  NEW_BASE<-NULL
  for (i in 1:length(BASIS)) {
    BASIS_temp<-BASIS[-i]
    temp<-MEASURE(TENS,BASIS_temp,Thres)
    #print(temp[[1]])
    if(!is.null(temp)){
      if(temp[[1]]>C){
        DIR<-i
        C<-temp[[1]]
        NEW_BASE<-temp[[2]]
      }
    }
  }
  
  if(C==-Inf){
    return(NULL)
  }else{
    BASIS_LIST<-list()
    for (i in 1:length(BASIS)) {
      if(i==DIR){
        BASIS_LIST[[i]]<-NEW_BASE
      }else{
        temp<-BASIS_convert(BASIS[[i]])
        BASIS_LIST[[i]]<-eval(parse(text = paste("TENS[",temp,"]",sep = "")))
      }
    }
    return(BASIS_LIST)
  }
}

WEAK_BASIS<-function(TENS,DIR){
  DIM<-dim(TENS)
  VEC<-apply(TENS,c(1:length(DIM))[-DIR],sum)
  NUM<-sort(VEC,decreasing = T)[1:2]
  ORDER<-which(VEC==NUM[1],arr.ind = T)
  ORDER2<-which(VEC==NUM[2],arr.ind = T)
  BASIS_use<-rep(0,DIM[DIR])
  for (i in 1:3) {
    if(nrow(ORDER)>1){
      ROW<-ORDER[sample(1:nrow(ORDER),2,replace = F),]
    }else{
      ROW<-rbind(ORDER,ORDER2[sample(1:nrow(ORDER2),1,replace = F),])
    }
    MAT<-matrix(NA,ncol = length(DIM),nrow = 2)
    MAT[,-DIR]<-ROW
    temp1<-BASIS_convert(MAT[1,])
    temp2<-BASIS_convert(MAT[2,])
    BASIS1<-eval(parse(text = paste("TENS[",temp1,"]",sep = "")))
    BASIS2<-eval(parse(text = paste("TENS[",temp2,"]",sep = "")))
    BASIS_use<-BASIS1*BASIS2
    if(sum(BASIS_use)>1){
      break
    }
  }
  
  if(sum(BASIS_use)>1){
    return(BASIS_use)
  }else{
    return(NULL)
  }
}

MEASURE_WS<-function(TENS,BASIS,DIR,Thres){
  DIM<-dim(TENS)
  
  if(DIR==1){
    temp<-""
  }else{
    temp<-"which(BASIS[[i]]==1)"
  }
  
  for (i in 2:length(DIM)) {
    if(i==DIR){
      temp<-paste(temp,",",sep = "")
    }else{
      temp<-paste(temp,',which(BASIS[[',i,']]==1)',sep="")
    }
  }
  
  VALUE<-eval(parse(text = paste("apply(TENS[",temp,"],",DIR,",sum)",sep = "")))
  
  NEW_BASIS<-rep(0,DIM[DIR])
  NEW_BASIS[which(VALUE>(Thres*max(VALUE)))]<-1
  
  DIM2<-eval(parse(text = paste("dim(TENS[",temp,"])",sep = "")))
  VAL<-1
  for (i in 1:length(DIM2)) {
    VAL=VAL*DIM2[i]
  }
  VAL<-VAL/DIM[DIR]*sum(NEW_BASIS)
  C_value<-sum(VALUE[as.logical(NEW_BASIS)])-VAL*Thres**length(DIM)
  if(sum(NEW_BASIS)<2){
    return(NULL)
  }else{
    return(list(C_value,NEW_BASIS))
  }
}


WEAK_SIG<-function(TENS,Thres){
  BASIS<-list()
  DIM<-dim(TENS)
  for (i in 1:length(DIM)) {
    BASIS[[i]]<-WEAK_BASIS(TENS,i)
  }
  
  IND<-sum(as.numeric(lapply(BASIS, is.null)))
  
  C=0
  DIR=0
  NEW_BASE<-NULL
  
  if(IND==0){
    
    for (i in 1:length(BASIS)) {
      temp<-MEASURE_WS(TENS,BASIS,i,Thres)
      if(!is.null(temp)){
        if(temp[[1]]>C){
          DIR<-i
          C<-temp[[1]]
          NEW_BASE<-temp[[2]]
        }
      }
    }
    if(DIR==0){
      return(NULL)
    }else{
      BASIS[[DIR]]<-NEW_BASE
      return(BASIS)
    }
    
  }else if(IND==1){
    
    for (i in 1:length(BASIS)) {
      if(is.null(BASIS[[i]])){
        temp<-MEASURE_WS(TENS,BASIS,i,Thres)
        if(!is.null(temp)){
          if(temp[[1]]>C){
            BASIS[[i]]<-temp[[2]]
            return(BASIS)
          }
        }else{
          return(NULL)
        }
      }
    }
  }else{
    return(NULL)
  }
}




GETF_CP<-function(TENS,Thres=0.7,B_num=20,COVER=0.8){
  
  BASIS<-list()
  for (i in 1:length(dim(TENS))) {
    BASIS[[i]]<-matrix(0,nrow = dim(TENS)[i],ncol = 0)
  }
  
  B_now<-1
  Cover_now<-0
  SUM<-sum(TENS)
  SUM_now<-0
  
  while(B_now<B_num&Cover_now<COVER){
    BASIS_temp<-GEO_SEG(TENS,Thres)
    if(!is.null(BASIS_temp)){
      temp<-"which(BASIS_temp[[1]]==1)"
      for (i in 2:length(BASIS_temp)) {
        temp<-paste(temp,paste("which(BASIS_temp[[",i,"]]==1)",sep = ""),sep = ",")
      }
      SUM_now<-SUM_now+eval(parse(text = paste("sum(TENS[",temp,"])",sep = "")))
      eval(parse(text=paste("TENS[",temp,"]<-0",sep = "")))
      for (i in 1:length(BASIS_temp)) {
        BASIS[[i]]<-cbind(BASIS[[i]],BASIS_temp[[i]])
      }
    }else{
      BASIS_temp<-WEAK_SIG(TENS,Thres)
      if(!is.null(BASIS_temp)){
        temp<-"which(BASIS_temp[[1]]==1)"
        for (i in 2:length(BASIS_temp)) {
          temp<-paste(temp,paste("which(BASIS_temp[[",i,"]]==1)",sep = ""),sep = ",")
        }
        SUM_now<-SUM_now+eval(parse(text = paste("sum(TENS[",temp,"])",sep = "")))
        eval(parse(text=paste("TENS[",temp,"]<-0",sep = "")))
        for (i in 1:length(BASIS_temp)) {
          BASIS[[i]]<-cbind(BASIS[[i]],BASIS_temp[[i]])
        }
      }else{
        break
      }
    }
    Cover_now<-SUM_now/SUM
    print(Cover_now)
    B_now<-B_now+1
  }
  
  return(BASIS)
}


Reconstruct_error<-function(TENS,BASIS){
  TENS_compare<-TENS*0
  
  DIM<-dim(BASIS[[1]])[2]
  for(i in 1:length(DIM)){
    temp<-paste("which(BASIS[[1]][,",i,"]==0)",sep = "")
    for(j in 2:length(BASIS)){
      temp<-paste(temp,",which(BASIS[[",j,"]][,",i,"]==1)",sep = "")
    }
    eval(parse(text=paste("TENS_compare[",temp,"]<-1",sep = "")))
  }
  return(sum(abs(TENS_compare-TENS)))
}


test<-GETF_CP(TENS,COVER = 0.9,Thres = 0.9,B_num = 200)
