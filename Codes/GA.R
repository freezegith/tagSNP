source("./Generate_Simulate_SNP_Data.R")
source("./datatransform.R")
source("./matrixToList.R")
source("./Predict.R")
source("./GeneticFunction.R")

LdList = matrixToList(Ld, 'value')
system("touch ./1.txt")
causalposname="/home/bioinformatics/Desktop/code/test/data/HTR2A650v3.txt"
causalpos=read.table(causalposname,header = T)
causalpos
N=1000
gene="HTR2A"
dieasegene=0
risk1=1.1
risk2=2.0
case=1000
control=1000

#generate SNP data

for (i in 1:5) {
  causal=causalpos[i,2]
  Generate_Simulate_SNP_Data(causal = causal, N = N, gene = gene, dieasegene = dieasegene,
                             risk1 = risk1, risk2 = risk2, case = case, control = control,
                             j = 1)
}


#transfer data
for (i in 1:5) {
  causal=causalpos[i,2]
  datatransform(causal = causal, N = N, gene = gene, dieasegene = dieasegene,
                risk1 = risk1, risk2 = risk2, case = case, control = control,
                j = 1)
}

genotypeName = ("/home/bioinformatics/Desktop/code/test/finalresult/46306571-0-1.1-2/genotype/1.txt")
genotype<-try(read.table(genotypeName, header = FALSE, 
                   stringsAsFactors = FALSE), silent = TRUE)
genotype
LdName=paste("/home/bioinformatics/Desktop/code/test/finalresult/46306571-0-1.1-2/RZ/1.txt")
Ld<-try(read.table(LdName, header = FALSE, 
                          stringsAsFactors = FALSE), silent = TRUE)
Ld

#p个染色体，k个tag SNP 
p = 50
tagSNPNumber = 10 
#SNPNumber表示SNP的总个数
SNPNumber = nrow(Ld)
#SNP数据
SNPData <- matrix(0, nrow = p, ncol =SNPNumber)
SNPData

#SNP数据初始化
for (i in 1:p) {
  for (j in 1:k) {
    ind <- sample(SNPNumber,1)
    f0 = ind
    if(SNPData[i,f0] == 1)
      j = j - 1
    else 
      SNPData[i,f0] = 1
  }
}


for (t in 1:50) {
  
  #求fitness
  fitness <- matrix(0, nrow = p, ncol = 2)
  
  for (i in 1:p) {
    for (j in 1:SNPNumber) {
      if(SNPData[i,j] == 1){
        for (k in 1:SNPNumber) {
          if(k != j){
            fitness[i,1] = fitness[i,1] + Ld[j,k]
            fitness[i,2] = i
          }
        }
      }
    }
  }
  
  #计算间距
  minBound = 0.8
  maxBound = 1.2
  dis = (maxBound - minBound) / (p - 1)
  
  #生成fitnessRank矩阵
  for (i in 1:p) {
    for(j in 1:(p - i)){
      if(length(fitness[j,1]) != 0 && length(fitness[j + 1, 1]) != 0){
        if(fitness[j, 1] > fitness[j + 1, 1]){
          temp = fitness[j, 1]
          fitness[j, 1] = fitness[j + 1, 1]
          fitness[j + 1, 1] = temp
          temp = fitness[j, 2]
          fitness[j, 2] = fitness[j + 1, 2]
          fitness[j + 1, 2] = temp
        }
      }
    }
  }
  
  fitnessRank <- matrix(0, nrow = p, ncol = 3)
  
  for (i in 1:p) {
    fitnessRank[i,1] = fitness[i,1]
    fitnessRank[i,2] = minBound + dis * (i - 1)
    fitnessRank[i,3] = fitness[i,2]
  }
  
  #选择操作
  fitnessSelection <- matrix(0, nrow = p, ncol = 3)
  j = 1
  for (i in 1:p) {
    if(fitnessRank[i,2] > 1){
      fitnessSelection[j,1] = fitnessRank[i,1]
      fitnessSelection[j,2] = fitnessRank[i,2] - 1
      fitnessSelection[j,3] = fitnessRank[i,3]
      j = j +  1
    }
  }
  i = 1
  while(j <= p){
    f1 = runif(1)
    if(fitnessRank[i,2] > f1){
      fitnessSelection[j,1] = fitnessRank[i,1]
      fitnessSelection[j,2] = fitnessRank[i,2]
      fitnessSelection[j,3] = fitnessRank[i,3]
      j = j + 1
      i = (i + 1) %% p
    }
  }
  
  SNPDatacol = ncol(SNPData)
  SNPDatarow = nrow(SNPData)
  SNPDataTemp <- matrix(0, nrow = SNPDatarow, ncol = SNPDatacol)
  for (i in 1:p) {
    for(j in 1:SNPDatacol)
    SNPDataTemp[i,j] = SNPData[fitnessSelection[i,3],j]
  }
  i
  SNPData = SNPDataTemp
  
  for (i in 1:nrow(fitnessSelection)) {
    cat(fitnessSelection[i,3])
    cat(" ")
  }
  
  
  
  #交叉
  crossRate = 0.8
  for (i in 1:(p - 1)) {
    f2 = runif(1)
    if(crossRate >= f2){
      ind <- sample(SNPNumber,1)
      r1 = ind
      ind <- sample(SNPNumber,1)
      r2 = ind
    }
    if(r1 - r2 >= 0){
      for (j in r2:r1) {
        temp = SNPData[i,j]
        SNPData[i + 1,j] = SNPData[i,j]
        SNPData[i,j] = temp
      }
    }
    if(r1 - r2 < 0){
      for (j in r1:r2) {
        temp = SNPData[i,j]
        SNPData[i + 1,j] = SNPData[i,j]
        SNPData[i,j] = temp
      }
    }
  }
  
  #变异
  mutationRate = 0.01
  for (i in 1:p) {
    for (j in 1:SNPNumber) {
      f3 = runif(1)
      if(f2 < mutationRate)
        SNPData[i,j] = 1 - SNPData[i,j]
    }
  }
  
  #修改
  for (i in 1:p) {
    k = 0
    for(j in 1:SNPNumber){
      if(SNPData[i,j] == 1)
        k = k+1
    }
    if(tagSNPNumber - k > 0) {
      differ = tagSNPNumber - k
      for(c in 1:differ){
        ind <- sample(SNPNumber,1)
        f3 = ind
        if(SNPData[i,f3] == 1)
          c = c - 1
        else 
          SNPData[i,f3] = 1
      }
    }
    if(tagSNPNumber - k < 0) {
      differ = k - tagSNPNumber
      for(c in 1:differ){
        ind <- sample(SNPNumber,1)
        f3 = ind
        if(SNPData[i,f3] == 0)
          c = c - 1
        else 
          SNPData[i,f3] = 0
      }
    }
  }
}

#求fitness
fitness <- matrix(0, nrow = p, ncol = 2)

for (i in 1:p) {
  for (j in 1:SNPNumber) {
    if(SNPData[i,j] == 1){
      for (k in 1:SNPNumber) {
        if(k != j){
          fitness[i,1] = fitness[i,1] + Ld[j,k]
          fitness[i,2] = i
        }
      }
    }
  }
}

#计算间距
minBound = 0.8
maxBound = 1.2
dis = (maxBound - minBound) / (p - 1)

#生成fitnessRank矩阵
for (i in 1:p) {
  for(j in 1:(p - i)){
    if(length(fitness[j,1]) != 0 && length(fitness[j + 1, 1]) != 0){
      if(fitness[j, 1] > fitness[j + 1, 1]){
        temp = fitness[j, 1]
        fitness[j, 1] = fitness[j + 1, 1]
        fitness[j + 1, 1] = temp
        temp = fitness[j, 2]
        fitness[j, 2] = fitness[j + 1, 2]
        fitness[j + 1, 2] = temp
      }
    }
  }
}

tagSNPSet<-vector(mode = "numeric", length = ncol(SNPData))
for (i in ncol(SNPData)) {
  tagSNPSet[i] = 0
}
for(i in 1:ncol(SNPData)){
  tagSNPSet[i] = SNPData[fitnessRank[p, 3], i]
}

 




