################################################################################################
#Codes for conducting real data analysis
################################################################################################
group = read.csv("./data/group.csv")
group = group[,c("Subject","Group")]
###############################
datapath = "./data/Time_Courses"
fileset = dir(datapath)
trace = F

# KKI
site = "KKI"
k = match(site, fileset)
p.t = 120

datapath.k = paste0(datapath,"/",fileset[k])
fileset.k = dir(datapath.k)
nk0 = length(fileset.k) - 1
phenotypic.k = read.csv(paste0(datapath.k,"/",fileset.k[length(fileset.k)]))
names(phenotypic.k)[1] = "ScanDirID"
phenotypic.k = phenotypic.k[,c("ScanDirID","QC_Rest_1")]
# phenotypic.k$ScanDirID = paste0("00",phenotypic.k$ScanDirID)
fileset.k.QC = NULL
for (ni in 1:nk0) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k[ni])
  filenamei = dir(datapath.k.i)
  if(length(filenamei) < 2){
    if(trace){
      print(paste0(fileset[k],"-",fileset.k[ni],"-warning"))
    }
  } else {
    QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
    if(length(QC.ki)>0){
      QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
      TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
      TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
      TC.i = read.table(TC.pathfilenamei,head=TRUE)
      data.i = TC.i[,-c(1,2)]
      pm.vec.k = dim(data.i)
      if(trace){
        print(pm.vec.k)
      }
      if(pm.vec.k[1] >= p.t){
        fileset.k.QC = c(fileset.k.QC,fileset.k[ni])
      }
    } else {
      if(trace){
        print(paste0(fileset[k],"-",fileset.k[ni],"-QC failed"))
      }
    }
  }
}

nk = length(fileset.k.QC)
data.A = array(0, dim = c(p.t, 116, nk))
group.k = group[1:nk,]
for (i in 1:nk) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k.QC[i])
  filenamei = dir(datapath.k.i)
  QC.ki = which(as.numeric(phenotypic.k[match(fileset.k.QC[i],phenotypic.k$ScanDirID),-1]) == 1)
  TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
  TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
  TC.i = read.table(TC.pathfilenamei,head=TRUE)
  data.i = TC.i[,-c(1,2)]
  data.A[,,i] = as.matrix(data.i[1:p.t,])
  
  g.i = match(paste0(c(site,fileset.k.QC[i]), collapse = "_"),group$Subject)
  group.k[i,] = group[g.i,]
}
group.KKI = group.k
A.data.KKI = data.A






# NeuroIMAGE
site = "NeuroIMAGE"
k = match(site, fileset)
p.t = 257

datapath.k = paste0(datapath,"/",fileset[k])
fileset.k = dir(datapath.k)
nk0 = length(fileset.k) - 1
phenotypic.k = read.csv(paste0(datapath.k,"/",fileset.k[length(fileset.k)]))
names(phenotypic.k)[1] = "ScanDirID"
phenotypic.k = phenotypic.k[,c("ScanDirID","QC_Rest_1")]
# phenotypic.k$ScanDirID = paste0("00",phenotypic.k$ScanDirID)
fileset.k.QC = NULL
for (ni in 1:nk0) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k[ni])
  filenamei = dir(datapath.k.i)
  if(length(filenamei) < 2){
    if(trace){
      print(paste0(fileset[k],"-",fileset.k[ni],"-warning"))
    }
  } else {
    QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
    if(length(QC.ki)>0){
      QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
      TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
      TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
      TC.i = read.table(TC.pathfilenamei,head=TRUE)
      data.i = TC.i[,-c(1,2)]
      pm.vec.k = dim(data.i)
      if(trace){print(pm.vec.k)}
      if(pm.vec.k[1] >= p.t){
        fileset.k.QC = c(fileset.k.QC,fileset.k[ni])
      }
    } else {
      if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-QC failed"))}
    }
  }
}

nk = length(fileset.k.QC)
data.A = array(0, dim = c(p.t, 116, nk))
group.k = group[1:nk,]
for (i in 1:nk) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k.QC[i])
  filenamei = dir(datapath.k.i)
  QC.ki = which(as.numeric(phenotypic.k[match(fileset.k.QC[i],phenotypic.k$ScanDirID),-1]) == 1)
  TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
  TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
  TC.i = read.table(TC.pathfilenamei,head=TRUE)
  data.i = TC.i[,-c(1,2)]
  data.A[,,i] = as.matrix(data.i[1:p.t,])
  
  g.i = match(paste0(c(site,fileset.k.QC[i]), collapse = "_"),group$Subject)
  group.k[i,] = group[g.i,]
}
group.NeuroIMAGE = group.k
A.data.NeuroIMAGE = data.A




# Peking_1
site = "Peking_1"
k = match(site, fileset)
p.t = 232

datapath.k = paste0(datapath,"/",fileset[k])
fileset.k = dir(datapath.k)
nk0 = length(fileset.k) - 1
phenotypic.k = read.csv(paste0(datapath.k,"/",fileset.k[length(fileset.k)]))
names(phenotypic.k)[1] = "ScanDirID"
phenotypic.k = phenotypic.k[,c("ScanDirID","QC_Rest_1")]
fileset.k.QC = NULL
for (ni in 1:nk0) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k[ni])
  filenamei = dir(datapath.k.i)
  if(length(filenamei) < 2){
    if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-warning"))}
  } else {
    QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
    if(length(QC.ki)>0){
      QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
      TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
      TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
      TC.i = read.table(TC.pathfilenamei,head=TRUE)
      data.i = TC.i[,-c(1,2)]
      pm.vec.k = dim(data.i)
      if(trace){print(pm.vec.k)}
      if(pm.vec.k[1] >= p.t){
        fileset.k.QC = c(fileset.k.QC,fileset.k[ni])
      }
    } else {
      if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-QC failed"))}
    }
  }
}

nk = length(fileset.k.QC)
data.A = array(0, dim = c(p.t, 116, nk))
group.k = group[1:nk,]
for (i in 1:nk) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k.QC[i])
  filenamei = dir(datapath.k.i)
  QC.ki = which(as.numeric(phenotypic.k[match(fileset.k.QC[i],phenotypic.k$ScanDirID),-1]) == 1)
  TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
  TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
  TC.i = read.table(TC.pathfilenamei,head=TRUE)
  data.i = TC.i[,-c(1,2)]
  data.A[,,i] = as.matrix(data.i[1:p.t,])
  
  g.i = match(paste0(c("Peking",fileset.k.QC[i]), collapse = "_"),group$Subject)
  group.k[i,] = group[g.i,]
}
group.Peking_1 = group.k
A.data.Peking_1 = data.A
nk1 = nk


# Peking_2
site = "Peking_2"
k = match(site, fileset)
p.t = 232

datapath.k = paste0(datapath,"/",fileset[k])
fileset.k = dir(datapath.k)
nk0 = length(fileset.k) - 1
phenotypic.k = read.csv(paste0(datapath.k,"/",fileset.k[length(fileset.k)]))
names(phenotypic.k)[1] = "ScanDirID"
phenotypic.k = phenotypic.k[,c("ScanDirID","QC_Rest_1")]
fileset.k.QC = NULL
for (ni in 1:nk0) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k[ni])
  filenamei = dir(datapath.k.i)
  if(length(filenamei) < 2){
    if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-warning"))}
  } else {
    QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
    if(length(QC.ki)>0){
      QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
      TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
      TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
      TC.i = read.table(TC.pathfilenamei,head=TRUE)
      data.i = TC.i[,-c(1,2)]
      pm.vec.k = dim(data.i)
      if(trace){ print(pm.vec.k)}
      if(pm.vec.k[1] >= p.t){
        fileset.k.QC = c(fileset.k.QC,fileset.k[ni])
      }
    } else {
      if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-QC failed"))}
    }
  }
}

nk = length(fileset.k.QC)
data.A = array(0, dim = c(p.t, 116, nk))
group.k = group[1:nk,]
for (i in 1:nk) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k.QC[i])
  filenamei = dir(datapath.k.i)
  QC.ki = which(as.numeric(phenotypic.k[match(fileset.k.QC[i],phenotypic.k$ScanDirID),-1]) == 1)
  TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
  TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
  TC.i = read.table(TC.pathfilenamei,head=TRUE)
  data.i = TC.i[,-c(1,2)]
  data.A[,,i] = as.matrix(data.i[1:p.t,])
  
  g.i = match(paste0(c("Peking",fileset.k.QC[i]), collapse = "_"),group$Subject)
  group.k[i,] = group[g.i,]
}
group.Peking_2 = group.k
A.data.Peking_2 = data.A
nk2 = nk

# Peking_3
site = "Peking_2"
k = match(site, fileset)
p.t = 232

datapath.k = paste0(datapath,"/",fileset[k])
fileset.k = dir(datapath.k)
nk0 = length(fileset.k) - 1
phenotypic.k = read.csv(paste0(datapath.k,"/",fileset.k[length(fileset.k)]))
names(phenotypic.k)[1] = "ScanDirID"
phenotypic.k = phenotypic.k[,c("ScanDirID","QC_Rest_1")]
fileset.k.QC = NULL
for (ni in 1:nk0) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k[ni])
  filenamei = dir(datapath.k.i)
  if(length(filenamei) < 2){
    if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-warning"))}
  } else {
    QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
    if(length(QC.ki)>0){
      QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
      TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
      TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
      TC.i = read.table(TC.pathfilenamei,head=TRUE)
      data.i = TC.i[,-c(1,2)]
      pm.vec.k = dim(data.i)
      if(trace){print(pm.vec.k)}
      if(pm.vec.k[1] >= p.t){
        fileset.k.QC = c(fileset.k.QC,fileset.k[ni])
      }
    } else {
      if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-QC failed"))}
    }
  }
}

nk = length(fileset.k.QC)
data.A = array(0, dim = c(p.t, 116, nk))
group.k = group[1:nk,]
for (i in 1:nk) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k.QC[i])
  filenamei = dir(datapath.k.i)
  QC.ki = which(as.numeric(phenotypic.k[match(fileset.k.QC[i],phenotypic.k$ScanDirID),-1]) == 1)
  TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
  TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
  TC.i = read.table(TC.pathfilenamei,head=TRUE)
  data.i = TC.i[,-c(1,2)]
  data.A[,,i] = as.matrix(data.i[1:p.t,])
  
  g.i = match(paste0(c("Peking",fileset.k.QC[i]), collapse = "_"),group$Subject)
  group.k[i,] = group[g.i,]
}
group.Peking_3 = group.k
A.data.Peking_3 = data.A
nk3 = nk


A.data.Peking = array(0, dim = c(p.t, 116, nk1+nk2+nk3))
A.data.Peking[,,1:nk1] = A.data.Peking_1
A.data.Peking[,,nk1+1:nk2] = A.data.Peking_2
A.data.Peking[,,nk1+nk2+1:nk3] = A.data.Peking_3

group.Peking = group[1:(nk1+nk2+nk3),]
group.Peking[1:nk1,] = group.Peking_1
group.Peking[nk1+1:nk2,] = group.Peking_2
group.Peking[nk1+nk2+1:nk3,] = group.Peking_3




# Pittsburgh(ID+00)
site = "Pittsburgh"
k = match(site, fileset)
p.t = 192

datapath.k = paste0(datapath,"/",fileset[k])
fileset.k = dir(datapath.k)
nk0 = length(fileset.k) - 1
phenotypic.k = read.csv(paste0(datapath.k,"/",fileset.k[length(fileset.k)]))
names(phenotypic.k)[1] = "ScanDirID"
phenotypic.k = phenotypic.k[,c("ScanDirID","QC_Rest_1")]
ScanDirID0 = phenotypic.k$ScanDirID
phenotypic.k$ScanDirID = paste0("00",phenotypic.k$ScanDirID)
fileset.k.QC = NULL
fileset.k.QC0 = NULL
for (ni in 1:nk0) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k[ni])
  filenamei = dir(datapath.k.i)
  if(length(filenamei) < 2){
    if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-warning"))}
  } else {
    QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
    if(length(QC.ki)>0){
      QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
      TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
      TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
      TC.i = read.table(TC.pathfilenamei,head=TRUE)
      data.i = TC.i[,-c(1,2)]
      pm.vec.k = dim(data.i)
      if(trace){print(pm.vec.k)}
      if(pm.vec.k[1] >= p.t){
        fileset.k.QC = c(fileset.k.QC,fileset.k[ni])
        fileset.k.QC0 = c(fileset.k.QC0, ScanDirID0[match(fileset.k[ni],phenotypic.k$ScanDirID)])
      }
    } else {
      if(trace){ print(paste0(fileset[k],"-",fileset.k[ni],"-QC failed"))}
    }
  }
}

as.numeric(fileset.k.QC)-fileset.k.QC0

nk = length(fileset.k.QC)
data.A = array(0, dim = c(p.t, 116, nk))
group.k = group[1:nk,]
for (i in 1:nk) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k.QC[i])
  filenamei = dir(datapath.k.i)
  QC.ki = which(as.numeric(phenotypic.k[match(fileset.k.QC[i],phenotypic.k$ScanDirID),-1]) == 1)
  TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
  TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
  TC.i = read.table(TC.pathfilenamei,head=TRUE)
  data.i = TC.i[,-c(1,2)]
  data.A[,,i] = as.matrix(data.i[1:p.t,])
  
  g.i = match(paste0(c(site,fileset.k.QC0[i]), collapse = "_"),group$Subject)
  group.k[i,] = group[g.i,]
}
group.Pittsburgh = group.k
A.data.Pittsburgh = data.A




# NYU (ID+00)
site = "NYU"
k = match(site, fileset)
p.t = 172

datapath.k = paste0(datapath,"/",fileset[k])
fileset.k = dir(datapath.k)
nk0 = length(fileset.k) - 1
phenotypic.k = read.csv(paste0(datapath.k,"/",fileset.k[length(fileset.k)]))
names(phenotypic.k)[1] = "ScanDirID"
phenotypic.k = phenotypic.k[,c("ScanDirID","QC_Rest_1","QC_Rest_2")]
ScanDirID0 = phenotypic.k$ScanDirID
phenotypic.k$ScanDirID[94:222] = paste0("00",phenotypic.k$ScanDirID[94:222])
fileset.k.QC = NULL
fileset.k.QC0 = NULL
for (ni in 1:nk0) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k[ni])
  filenamei = dir(datapath.k.i)
  if(length(filenamei) < 2){
    if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-warning"))}
  } else {
    QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
    if(length(QC.ki)>0){
      QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
      if(length(filenamei) < length(filenamei)/2 + QC.ki[1]){
        TC.filenamei = filenamei[length(filenamei)/2 + 1]
      }else{
        TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
      }
      TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
      TC.i = read.table(TC.pathfilenamei,head=TRUE)
      data.i = TC.i[,-c(1,2)]
      pm.vec.k = dim(data.i)
      if(trace){print(pm.vec.k)}
      if(pm.vec.k[1] >= p.t & length(filenamei) >= length(filenamei)/2 + QC.ki[1]){
        fileset.k.QC = c(fileset.k.QC,fileset.k[ni])
        fileset.k.QC0 = c(fileset.k.QC0, ScanDirID0[match(fileset.k[ni],phenotypic.k$ScanDirID)])
      }
    } else {
      if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-QC failed"))}
    }
  }
}

as.numeric(fileset.k.QC)-fileset.k.QC0

nk = length(fileset.k.QC)
data.A = array(0, dim = c(p.t, 116, nk))
group.k = group[1:nk,]
for (i in 1:nk) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k.QC[i])
  filenamei = dir(datapath.k.i)
  QC.ki = which(as.numeric(phenotypic.k[match(fileset.k.QC[i],phenotypic.k$ScanDirID),-1]) == 1)
  TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
  TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
  TC.i = read.table(TC.pathfilenamei,head=TRUE)
  data.i = TC.i[,-c(1,2)]
  data.A[,,i] = as.matrix(data.i[1:p.t,])
  
  g.i = match(paste0(c(site,fileset.k.QC0[i]), collapse = "_"),group$Subject)
  group.k[i,] = group[g.i,]
}
group.NYU = group.k
A.data.NYU = data.A





# OHSU
site = "OHSU"
k = match(site, fileset)
p.t = 74

datapath.k = paste0(datapath,"/",fileset[k])
fileset.k = dir(datapath.k)
nk0 = length(fileset.k) - 1
phenotypic.k = read.csv(paste0(datapath.k,"/",fileset.k[length(fileset.k)]))
names(phenotypic.k)[1] = "ScanDirID"
phenotypic.k = phenotypic.k[,c("ScanDirID","QC_Rest_1","QC_Rest_2","QC_Rest_3")]
# phenotypic.k$ScanDirID = paste0("00",phenotypic.k$ScanDirID)
fileset.k.QC = NULL
for (ni in 1:nk0) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k[ni])
  filenamei = dir(datapath.k.i)
  if(length(filenamei) < 2){
    if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-warning"))}
  } else {
    QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
    if(length(QC.ki)>0){
      QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
      TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
      TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
      TC.i = read.table(TC.pathfilenamei,head=TRUE)
      data.i = TC.i[,-c(1,2)]
      pm.vec.k = dim(data.i)
      if(trace){print(pm.vec.k)}
      if(pm.vec.k[1] >= p.t){
        fileset.k.QC = c(fileset.k.QC,fileset.k[ni])
      }
    } else {
      if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-QC failed"))}
    }
  }
}

nk = length(fileset.k.QC)
data.A = array(0, dim = c(p.t, 116, nk))
group.k = group[1:nk,]
for (i in 1:nk) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k.QC[i])
  filenamei = dir(datapath.k.i)
  QC.ki = which(as.numeric(phenotypic.k[match(fileset.k.QC[i],phenotypic.k$ScanDirID),-1]) == 1)
  TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
  TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
  TC.i = read.table(TC.pathfilenamei,head=TRUE)
  data.i = TC.i[,-c(1,2)]
  data.A[,,i] = as.matrix(data.i[1:p.t,])
  
  g.i = match(paste0(c(site,fileset.k.QC[i]), collapse = "_"),group$Subject)
  group.k[i,] = group[g.i,]
}
group.OHSU = group.k
A.data.OHSU = data.A



# WashU (ID+00)
site = "WashU"
k = match(site, fileset)
p.t = 128

datapath.k = paste0(datapath,"/",fileset[k])
fileset.k = dir(datapath.k)
nk0 = length(fileset.k) - 1
phenotypic.k = read.csv(paste0(datapath.k,"/",fileset.k[length(fileset.k)]))
phenotypic.k = phenotypic.k[,c("ScanDirID","QC_S1_Rest_1","QC_S1_Rest_2","QC_S1_Rest_3",
                               "QC_S1_Rest_4","QC_S1_Rest_5","QC_S1_Rest_6")]
ScanDirID0 = phenotypic.k$ScanDirID
phenotypic.k$ScanDirID = paste0("00",phenotypic.k$ScanDirID)
fileset.k.QC = NULL
fileset.k.QC0 = NULL
for (ni in 1:nk0) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k[ni])
  filenamei = dir(datapath.k.i)
  if(length(filenamei) < 2){
    if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-warning"))}
  } else {
    QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
    if(length(QC.ki)>0){
      QC.ki = which(as.numeric(phenotypic.k[match(fileset.k[ni],phenotypic.k$ScanDirID),-1]) == 1)
      TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
      TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
      TC.i = read.table(TC.pathfilenamei,head=TRUE)
      data.i = TC.i[,-c(1,2)]
      pm.vec.k = dim(data.i)
      if(trace){print(pm.vec.k)}
      if(pm.vec.k[1] >= 128){
        fileset.k.QC = c(fileset.k.QC,fileset.k[ni])
        fileset.k.QC0 = c(fileset.k.QC0, ScanDirID0[match(fileset.k[ni],phenotypic.k$ScanDirID)])
      }
    } else {
      if(trace){print(paste0(fileset[k],"-",fileset.k[ni],"-QC failed"))}
    }
  }
}


as.numeric(fileset.k.QC) - fileset.k.QC0

p.t = 128 # the dimension of time series
nk = length(fileset.k.QC)
data.A = array(0, dim = c(p.t, 116, nk))
group.k = group[1:nk,]
for (i in 1:nk) {
  datapath.k.i = paste0(datapath.k,"/",fileset.k.QC[i])
  filenamei = dir(datapath.k.i)
  QC.ki = which(as.numeric(phenotypic.k[match(fileset.k.QC[i],phenotypic.k$ScanDirID),-1]) == 1)
  TC.filenamei = filenamei[length(filenamei)/2 + QC.ki[1]]
  TC.pathfilenamei = paste0(datapath.k.i,"/",TC.filenamei)
  TC.i = read.table(TC.pathfilenamei,head=TRUE)
  data.i = TC.i[,-c(1,2)]
  data.A[,,i] = as.matrix(data.i[1:p.t,])
  
  g.i = match(paste0(c(site,fileset.k.QC0[i]), collapse = "_"),group$Subject)
  group.k[i,] = group[g.i,]
}
group.WashU = group.k
A.data.WashU = data.A


fMRI.data.list = list(data.KKI=A.data.KKI,
                      data.NeuroIMAGE=A.data.NeuroIMAGE,
                      data.Peking=A.data.Peking,
                      data.Pittsburgh=A.data.Pittsburgh,
                      data.NYU=A.data.NYU,
                      data.OHSU=A.data.OHSU,
                      data.WashU=A.data.WashU)
fMRI.group.list = list(group.KKI=group.KKI,
                       group.NeuroIMAGE=group.NeuroIMAGE,
                       group.Peking=group.Peking,
                       group.Pittsburgh=group.Pittsburgh,
                       group.NYU=group.NYU,
                       group.OHSU=group.OHSU,
                       group.WashU=group.WashU)

k=1
fMRI.Typical.list = fMRI.data.list
fMRI.ADHD.list = fMRI.data.list
for (k in 1:length(fMRI.data.list)) {
  group.k = fMRI.group.list[[k]]
  nk = dim(group.k)[1]
  n.Typ.k = which(group.k$Group == "Typically Developing")
  n.ADHD.k = setdiff(1:nk, n.Typ.k)
  data.k = fMRI.data.list[[k]]
  fMRI.Typical.list[[k]] = data.k[,,n.Typ.k]
  fMRI.ADHD.list[[k]] = data.k[,,n.ADHD.k]
}
fMRI.data.list = list(Typical = fMRI.Typical.list, ADHD = fMRI.ADHD.list)


save(fMRI.data.list, file = paste0("./results/","fMRI_data_list",".RData"))

