# 程序说明
# 当前版本： 分箱算法v8.0

# 版本号：分箱算法v8.0
# 更新时间：2016.8.5
# 作者：Crystal
# 更新内容：
# 1) 修改了value函数，解决了分割点0不出现问题。
# 2) 修改了DiscThree函数，把第二次分割的样本数限制条件改为总样本数*比例。

# 版本号：分箱算法v6.0
# 更新时间：2016.7.23
# 作者：Crystal
# 更新内容：
# 1) 修改了ClassM函数，解决table(x,y)，x,y must have the same length 问题。

# 版本号：分箱算法v5.0
# 更新时间：2016.7.23
# 作者：Crystal
# 更新内容：
# 1) 修复了ClassM中有关chi算法最后合并样本数量少的类别的bug。

# 版本号：分箱算法v4.0
# 更新时间：2016.7.13
# 作者：Crystal
# 更新内容：
# 1）在Disc，DiscThree，ClassM中加入了useNA参数（默认为FALSE）来决定是否支持NA，如果支持NA（useNA=TRUE）则按照按照之前说的有空值的分箱算法，否则会报错并停止运行。
# 2) 删除了characchange,修改了Attr，记得导入数据时用StringAsFactor = FALSE哦。
# 3）修改了Attr来防止在取出某个元素a进行判别时取到了NA。
# 4) 修改了ClassM函数，使得正常运行。
# 5) 增加了GetNAidx函数来找到为NULL，NA，“”的元素下标、
# 6）修改timeRecog考虑了NAidx的情况。

# 版本号：分箱算法v3.0
# 更新时间：2016.7.12
# 作者：Jacky
# 更新内容：
# 1）结合web版的逻辑，在主函数Disc中加入了typelist参数接口，变量类型统一由java后台判断后传值，R判断只作为一个验证
# 2）将函数中所有变量类型由汉字描述改成了英文描述。即“时间”改为“dateTime”，“分类”改成“string”，“连续”改成“double”
# 3）将timeRecog函数改成levi优化后的timeRecog函数，只支持4种时间格式，分别是：“2016/8/1 0:10”、“2016/8/2”、“1:16:00”、“1:16”，其它时间格式类型会在java后台转换成这四种类型之一
# 4）修正了Dis函数里面处理时间型变量部分中res<-DiscThree(.....)这一句代码中的method改成了methodNum，因为和Dis函数的参数赋值不对应，不修改会报错

library(lubridate)

# 综合算法 用于处理分类变量、连续变量、和时间变量
# Disc --------------------------------------------------------------------
Disc <- function(data, typelist, prop1=0.2, prop2=0.1, prop3=0.2, prop4=0.9, n=100,
                 methodNum='ent', methodFac='ent', three=FALSE, alpha=0.05,
                 sigma=0, beta=0.05, useNA=FALSE){
  
  #prop1, prop3指分箱后的样本数需要满足的最小比例
  #prop2, prop4指分箱后合格率差值需要满足最小比例
  #n表示当连续变量因子数过多的时候预分箱参数
  #methodNum连续变量默认分箱算法为ent
  #methodFac类别变量合并默认算法为ent
  #three默认分为两箱
  #alpha类别变量合并时chi算法的参数
  #sigma类别变量合并时ent算法的参数
  #beta限制当某类的样本数量小于总样本的beta的时候把该类合并到和其合格率最接近的类别
  #useNA表示是否支持空值，如果不支持就会报错
  
  p <- ncol(data)
  datanew <- data #复制数据
  typelistR <- c() #用于储存变量类型
  cutlist <- list() #分割点
  ratiolist <- list() #合格率列表
  numlist <- list() #样本数列表
  for (i in 1:(p-1)){
    typelistR[i] <- Attr(datanew[,i]) #找到所有变量的类别
    if(typelistR[i] != typelist[i])
      print(paste("warning:第",i,"列变量R判断的数据类型与java判断的数据类型不一致",sep=""))
  }
  idxTime <- which(typelist=='dateTime') #找到时间变量
  idxDisc <- which(typelist=='string') #找到类别变量
  idxCont <- which(typelist=='double') #找到连续变量
  valchange <- function(idx, res, time=FALSE){ #赋值函数 把discdata cutlist ratiolist numlist进行更新
    cuts <- res$cutlist
    if (time){ #如果是时间变量的化比较麻烦 要根据分割点下标找到原始数据中的时间
      for (i in 1:length(cuts)){
        cuts[[i]] <- as.character(data[cuts[[i]], idx[i]])
     }
    }
    cutlist[idx] <<- cuts #传递分割点列表
    discx <- res$discdata
    datanew[, idx] <<- discx[, -ncol(discx)] #得到离散后的数据
    ratiolist[idx] <<- res$ratiolist
    numlist[idx] <<- res$numlist
  }
  if (length(idxDisc) > 0) { #对分类变量进行处理
    dataDisc <- datanew[, c(idxDisc, p)]
    res <- ClassM(dataDisc, methodFac, alpha, sigma, beta, useNA)
    valchange(idxDisc, res) #更新
  }
  if (length(idxCont) > 0) { #对连续变量进行处理
    dataCont <- datanew[, c(idxCont, p)]
    res <- DiscThree(dataCont, prop1, prop2, prop3, prop4, n, methodNum, three, useNA)
    valchange(idxCont, res)
  }
  if (length(idxTime) > 0){ #对时间变量进行处理 转换为连续变量后分箱
    dataTemp <- matrix(sapply(as.matrix(datanew[, idxTime]), timeRecog),
                       ncol = length(idxTime)) #处理时间
    datanew[, idxTime] <- dataTemp
    dataTime <- datanew[, c(idxTime, p)]
    res <- DiscThree(dataTime, prop1, prop2, prop3, prop4, n, methodNum, three, useNA, fixed = TRUE) #对于时间变量的分割点取数据中已有的值
    valchange(idxTime, res, time = TRUE)
  }
  return (list(cutlist = cutlist, discdata = datanew, ratiolist = ratiolist, numlist = numlist))
}

# Attr --------------------------------------------------------------------
Attr <- function(attr){ #判断变量类型的函数
  cond1 <- length(grep('[年月日时分秒/:-]', attr)) > 0 #如果有年、月、日、-、/、:为时间
  NAidx <- GetNAidx(attr)
  if (length(NAidx)){
    temp <- attr[-NAidx] #取出一个
    a <- as.character(temp[1])
  }
  else{
    a <- as.character(attr[1]) #取出一个
  }
  cond2 <- length(unlist(strsplit(a, '[年月日时分秒/:-]'))) > 1 #如果只有单个的年、月、日、时、分、秒则还是看作类别变量
  if (class(attr) == 'factor' || class(attr) == 'character'){
    if (cond1 & cond2)
      result <- 'dateTime'
    else
      result <- 'string'
  }
  if (class(attr) == 'integer' || class(attr) == 'numeric'){
    result <- 'double'
  }
  return(result)
}

# timeRecog ---------------------------------------------------------------
timeRecog<- function(Dat){ #时间识别即转化为数字函数
  Sys.setlocale("LC_TIME", "C") #保证as.Date运行不出现NA
  Dat <- as.character(Dat)
  Fir <- NULL    #储存日期
  Sec <- NULL    #储存具体时间
  Dat_F <- NULL  #临时储存日期
  Dat_S <- NULL  #临时储存时间
  num1 <- 0      #储存日期的换算出的数字
  num2 <- 0      #储存具体时间换算出的数字
  
  NAidx <- GetNAidx(Dat)
  if (length(NAidx))
      return (Dat)  
  
  #条件判断
  cond_sf <- grepl('[:]', Dat)                    #搜索数据中是含有：的列
  thre <- '[ :/-]'                                #含有 ： /  - 的正则表达式
  splitDat <- unlist(strsplit(Dat, " "))          #把日期和时间分离

  #如果是年月日时
  if (length(splitDat) == 2){
    Dat_F <- unlist(strsplit(splitDat[1], thre))   #将含有 ： /  - 都去掉， Dat_F日期年月日
    Dat_S <- unlist(strsplit(splitDat[2], thre))   #将含有 ： /  - 都去掉， Dat_F日期年月日
  }

  else {
    if (cond_sf){
      Dat_S <- unlist(strsplit(splitDat, thre))
    }
    else {
      Dat_F <- unlist(strsplit(splitDat, thre))
    }
  }

  #标准化年月日
  if (!is.null(Dat_F)){
    temp <- paste(Dat_F, collapse = '-')
    if (length(Dat_F)==3){
      Fir <- ymd(temp)
    }
    else {
      Fir <- as.Date(temp, '%m-%d')
    }
  }

  #标准化时间
  if (!is.null(Dat_S)){
    temp <- paste(Dat_S, collapse = ':')
    if (length(Dat_S) == 3){
      Sec <- hms(temp)
    }
    else {
      Sec <- hm(temp)
    }
  }

  #计算年月日和时间的数值
  if (!is.null(Fir)){
    num1 <- yday(Fir) + (year(Fir) - 1990)*367
  }
  if (!is.null(Sec)){
    num2 <- round(period_to_seconds(Sec) / as.numeric(ddays(1)), 3)
  }
  result <- num1 + num2

  return(result)
}


# DiscThree ---------------------------------------------------------------
DiscThree <- function(data, prop1=0.2, prop2=0.1, prop3=0.2, prop4=0.9, n=100,
                      method='ent', three=FALSE, useNA=FALSE, fixed=FALSE){
  
  #prop1, prop3指分箱后的样本数需要满足的最小比例
  #prop2, prop4指分箱后合格率差值需要满足最小比例
  #n表示当连续变量因子数过多的时候预分箱参数
  #method默认分箱算法为ent three默认分为两箱
  #fixed表示分割点取列表中存在的数据 返回该数据下标
  
  p <- dim(data)[2] - 1 #连续变量属性个数
  N <- nrow(data) #样本数
  y <- data[, (p + 1)] #目标类别变量
  datanew <- data #储存离散后的数据
  ynew <- tarchange(y) #改变目标变量类型
  cutlist <- list() #储存分割点
  ratiolist <- list() #储存区间合格率列表
  numlist <- list() #储存区间样本数
  getminprop <- function(prep, diff=NULL){ #求分箱合格率，样本比例限制条件
    if (is.null(diff)) #diff为空的时候求的是第一次分割的限制条件
      minprop <- c(prep * prop1, (1 - prep) * prop1,  (1 - prep) * prop2)
    # prep * prop1 得到区间中合格率高样本数在总样本中的最小比例(第一次分割)
    # (1 - prep) * prop1 得到区间中合格率低样本数在总样本中的最小比例
    # (1 - prep) * prop2 得到区间合格率差值最小比例
    else  #diff不为空的时候求的就是第二次分割限制条件
      minprop <- c(prep * prop3, (1 - prep) * prop3,  diff * prop4)
    # prep * prop3 得到区间中合格率高样本数在总样本中的最小比例(第二次分割)
    # (1 - prep) * prop3 得到区间中合格率低样本数在总样本中的最小比例
    # (1 - prep) * prop4 得到区间合格率差值最小比例
    return (list(prep = prep, minprop = minprop))
  }
  cutIndex1 <- function(x0, xi, n){ #当因子数过多的时候进行等分预分箱
    len <- length(x0)
    if (length(xi) > n) {
      seqidx <- seq(1, len, by = round(len / n))
      cuttemp <- unique(x0[seqidx])
      cutidx <- xi[x0[xi] %in% cuttemp]
    }
    else {
      cutidx <- xi
    }
    xcut <- (x0[cutidx] + x0[cutidx+1]) / 2
    return(xcut)
  }
  cutIndex2 <- function(cuts, attr, tar, res, step = 1){ #分箱函数
    prep <- res$prep
    minprop <- res$minprop #找到合格率，样本数需要满足的限制条件
    p1 <- minprop[1]
    p2 <- minprop[2]
    p3 <- minprop[3]
    m <- length(cuts) #分割候选点个数
    n <- length(tar) #样本总数
    dlist <- c(min(attr), max(attr)) #初始化分割点列表
    initEnt <- Inf
    entropy <- initEnt #初始化熵值为很大的数
    bestc <- 0 #初始化c值
    ci <- NULL #初始化分割点
    ri <- prep #初始化合格率
    di <- dlist #初始化分割点列表
    diffi <- 0 #初始化合格率差值
    for (i in 1:m) { #从分割点候选列表中找到最佳分割点
      ct <- cuts[i]
      d <- unique(insert(ct, dlist))
      dff <- findInterval(attr, d, rightmost.closed = TRUE)
      res <- ratio(dff, tar)
      rate <- res$yprop
      diff <- res$diff
      if (step == 1){
        prop <- res$num / n
      }
      else if (step == 2){
        prop <- res$num / N
      }
      idx <- which.max(rate)
      if ((prop[idx] >= p1) & (prop[-idx] >= p2) & (diff >= p3)){
        if (method=='ent'){ #ent
          wx <- which(attr <= ct)
          wn <- length(wx) / n
          e1 <- wn * ent(tar[wx])
          e2 <- (1 - wn) * ent(tar[-wx])
          val <- e1 + e2
          if (val < entropy){
            entropy <- val
            ci <- ct
            ri <- rate
            diffi <- diff
          }
        }
        else{ #caim、modcaim、ameva、cacc、diff
          df <- table(dff, tar)
          if (method=='caim')
            c <- caimm(df)
          else if (method=='ameva')
            c <- amevam(df)
          else if (method=='cacc')
            c <- caccm(df)
          else if (method=='modcaim')
            c <- modcaimm(df)
          else if (method=='diff')
            c <- diff
          if (c > bestc){
            bestc <- c
            ci <- ct
            ri <- rate
            diffi <- diff
          }
        }
      }
    }
    di <- insert(ci, dlist)
    return(list(cut = di, rate = ri, diff = diffi)) #返回分割点列表和合格率
  }
  cutIndex3 <- function(attr, tar, xcut, minprop, part, cutpart){ #进行二次分箱的函数
    ca <- NULL #初始化分割点
    e <- 0 #初始化熵值
    if (any(cutpart==TRUE)) {
      attrNew <- attr[part]
      tarNew <- tar[part]
      xcutNew <- xcut[cutpart]
      res <- cutIndex2(xcutNew, attrNew, tarNew, minprop, step = 2)
      ci2 <- res$cut
      if (length(ci2) > 2)
        ca <- ci2[2]
      e <- ent(tarNew)
    }
    return(list(c = ca, e = e))
  }
  for(i in 1:p){
    attr <- datanew[, i]
    tar <- ynew
    NAidx <- GetNAidx(attr)
    if (length(NAidx)){
      if (!useNA){
        stop("数据含有空值")
        # stop("NA")
      }
      else{
        attr <- as.numeric(attr[-NAidx]) #把值不为NA的拿来分箱
        tar <- ynew[-NAidx]
      }
    }
    prep <- length(which(tar == 1)) / length(tar) #合格率
    minprop <- getminprop(prep) #找到minprop的几个限制条件
    od <- order(attr)
    x0 <- attr[od]
    y0 <- ynew[od]
    xi <- which(diff(x0) != 0) #找到第i个属性中相邻两值不同的样本d[i+1]-d[i] != 0
    if(!length(xi)){
      ci <- x0[1:2]
    }
    else {
      xcut <- cutIndex1(x0, xi, n) #等分预分箱
      res <- cutIndex2(xcut, attr, tar, minprop) #第一次分割
      ci <- res$cut
      ri <- res$rate
      diffi <- res$diff
      if ((length(ci) > 2) & three) { #如果继续分割的化
          cx <- ci[2] # 得到中间的分割点
          minprop <- getminprop(prep, diffi) #得到第二次分割限制条件
          resl <- cutIndex3(attr, tar, xcut, minprop, attr < cx, xcut < cx) #对A区间进行分割
          resr <- cutIndex3(attr, tar, xcut, minprop, attr > cx, xcut > cx) #对B区间进行分割
          ca <- c(resl$c, resr$c)
          ens <- c(resl$e, resr$e)
          idx <- which.max(ens)
          if (length(ca) < 2) #如果只有一个区间或没有区间可以分割
            ci <- insert(ca, ci)
          else #如果两个区间都可以分割则分割熵值大的区间
            ci <- insert(ca[idx], ci)
      }
    }
    xdisc <- findInterval(attr, ci, rightmost.closed = TRUE)
    if (length(NAidx))
      datanew[-NAidx, i] <- xdisc
    else 
      datanew[, i] <- xdisc
    if (fixed){ #对于时间变量设计
      for (k in 1:length(ci)){
         cnew <- max(x0[x0<=ci[k]]) #找到最接近且小于ci的数据cnew
         idx <- sort(which(data[, i]==cnew))[1] #找到最小的下标
         ci[k] <- idx #把下标idx赋值给ci[k+1]
      }
    }
    cutlist[[i]] <- unique(ci)
    res <- ratio(xdisc, tar)
    ratiolist[[i]] <- res$yprop
    numlist[[i]] <- res$num
  }
  return (list(cutlist = cutlist, discdata = datanew, ratiolist = ratiolist, numlist = numlist))
}

# measure -----------------------------------------------------------------
# 测度函数
modcaimm <- function(df){#修正的caim测度
  nr <- dim(df)[1] #区间个数
  nc <- dim(df)[2] #类别个数
  okr <- df[,2] #每区间ok的个数
  sumr <- apply(df, 1, sum) #每区间中所有类的总数
  return(sum(okr^2/sumr)/nr) #返回caim测度
}
caimm <- function(df){#caim测度
  nr <- dim(df)[1] #区间个数
  nc <- dim(df)[2] #类别个数
  maxr <- apply(df, 1, max) #每区间最多的类个数
  sumr <- apply(df, 1, sum) #每区间中所有类的总数
  return(sum(maxr^2/sumr)/nr) #返回caim测度
}
amevam <- function (tb){ #ameva测度
  nr <- dim(tb)[1]
  nc <- dim(tb)[2]
  den <- nr * (nc - 1)
  val = chiSq(tb)/den
  return(val)
}
caccm <- function (tb){#cacc测度
  n <- sum(tb)
  nr <- dim(tb)[1] + 1
  logn <- 0
  if (nr > 1)
    logn <- log(nr)
  yp <- chiSq(tb)
  val = sqrt(yp/(yp + n * logn))
  return(val)
}
chiSq <- function (tb){
  #计算列联表中的卡方值 卡方值可以用于度量多个群体或属性之间的相关性
  #卡方值越大则相关性越小
  tb <- tb + 1e-04
  e <- tb
  n <- sum(tb)
  p <- dim(tb)[1]
  q <- dim(tb)[2]
  mi <- numeric(p)
  ni <- numeric(q)
  for (i in 1:q) ni[i] <- sum(tb[, i])
  for (i in 1:p) {
    mi[i] <- sum(tb[i, ])
    e[i, ] <- mi[i] * ni/n
  }
  val <- sum((tb - e)^2/e)
  return(val)
}
entSq <- function(tb){ #类型变量合并时的评价标准
  m <- sum(tb)
  Hnew <- ent(tb[1, ] + tb[2, ], type = 1)
  Hold <- ent(tb, type = 2)
  f <- as.numeric(m * (Hnew - Hold))
  return(f)
}
ent <- function (y, type=0) {#计算信息熵的函数
  if (type == 0){ #当给的是array
    p <- prop.table(table(y))
    e <- -sum(p * mylog(p))
  }
  else if (type == 1){ #当给的是table的某一行
    num <- as.numeric(y)
    m <- sum(num)
    p <- num / m
    e <- -sum(p * mylog(p))
  }
  else if (type == 2){ #当给出的是整个table的时候算出整个table的熵值
    if (is.null(dim(y)))
      ma <- matrix(y, nrow = 1)
    else
      ma <- matrix(y, dim(y)[1])
    prop <- rowSums(ma) / sum(ma)
    entropy <- apply(ma, 1, ent, type = 1)
    e <- sum(prop * entropy)
  }
  return(e)
}

# ClassM --------------------------------------------------------------------
# 信息熵或chiMerge算法用于合并类别变量
ClassM <- function (data, type = 'ent', alpha = 0.05, sigma = 0, beta = 0.05, useNA=FALSE){
  
  # 用于分类变量的合并
  # type 为合并方法，分为ent和chi，默认为ent
  # alpha为chi算法中的参数 alpha越大则使得分类总数越少
  # sigma为ent算法的参数 sigma越大使得分类总数越少
  # beta限制当某类的样本数量小于总样本的某个比值的时候把该类合并到和其合格率最接近的类别
  
  p <- dim(data)[2] - 1 #数据中分类变量属性总个数
  y <- data[, (p+1)] #目标类别变量
  ynew <- tarchange(y) #目标变量类型转换
  datanew <- data #复制数据
  cutlist <- list() #用于储存分类变量合并后的分割点
  ratiolist <- list() #用于储存合格率列表
  numlist <- list() #用于储存样本列表
  for (i in 1:p) { #对每个分类变量进行处理
    x <- datanew[, i] #提取第i个变量的值
    tar <- ynew 
    NAidx <- GetNAidx(x)
    if (length(NAidx)){
       if (!useNA){
         stop("数据含有空值")
         # stop("NA")
       }  
       else{
         x <- x[-NAidx]
         tar <- ynew[-NAidx]
       }
    }
    cuts <- value(x, tar, type, alpha, sigma, beta) #得到第i个变量合并后的分割点列表
    discx <- merg(cuts, x) #得到处理之后的第i个变量的值
    if (length(NAidx))
      datanew[-NAidx,i] <- discx #将第i个属性中的分类变量进行合并 改变类别名称
    else
      datanew[, i] <- discx
    cutlist[[i]] <- cuts #将第i个变量分割点列表储存
    res <- ratio(discx, tar) #得到区间合格率和样本数
    ratiolist[[i]] <- res$yprop #合格率
    numlist[[i]] <- res$num #样本数
  }
  return(list(cutlist = cutlist,  discdata = datanew,
              ratiolist = ratiolist, numlist = numlist)) #返回所有属性的分割点列表和类别合并之后的数据
}

value <- function(x, y, type, alpha, sigma, prop){ #算法主函数
  n <- length(x) #样本个数
  s <- length(unique(y)) #属性中类别总个数
  if (s == 1){ #如果只有一类则全部合并在一起
    cutpoint <- paste(unique(x), collapse = '|')
    return (cutpoint)
  }
  tb <- as.array(table(x,y)) #得到目标函数值在不同类别中的分布
  NULLidx <- which((tb[, 1] + tb[, 2]) == 0)
  if (length(NULLidx))
    tb <- tb[-NULLidx, ]
  cutpoint <- rownames(tb) #初始化分割列表为类别总个数
  if (type == 'ent'){
    threshold <- n * log(2)
    repeat {
      m <- dim(tb)[1] #此时类别总个数
      if (is.null(m) || m < 2) #当此时类别数小于2时停止
        break
      entroOld <- ent(tb, type = 2) #此时合并之前的熵值
      testlist <- matrix(rep(threshold, m^2), nrow = m) #用于储存两个类别之间评价标准的矩阵
      for (i in 1 : (m - 1)){ #计算出不同类别两两之间的评价值并储存
        for (j in (i + 1) : m){
          d <- tb[c(i, j), ] #抽取第i，j类别在不同目标值中的分布并构造列联表
          testlist[i,j] = entSq(d) #计算出它们之间的评价标准
        }
      }
      k<- min(testlist) #找到最小的评价指标
      idx <- which(testlist==k, arr.ind = TRUE)[1,] #找到矩阵中等于最小准则值的位置
      i <- idx[1]
      j <- idx[2] #i，j即为将要合并的类别
      tbNew <- tb
      tbNew[i, ] <- tbNew[i, ] + tbNew[j, ] #将i，j两个类进行合并 并将合并后的值赋予类别i
      tbNew <- tbNew[-j, ] #除掉类别j的所有值
      entroNew <- ent(tbNew, type = 2)
      entrodiff <- entroNew - entroOld #熵值差
      entrolev <- ifelse((entroOld==0)&(entrodiff==0), 0, entrodiff / entroOld)
      if (entrolev > sigma) #当整体的熵值开始增大的时候（有一个松弛到sigma，即合并之后熵值与原来熵值差占原来熵值比例为sigma的时候）停止合并
        break
      tb <- tbNew
      cutpoint[i] <- paste(cutpoint[i], '|', cutpoint[j]) #把类别i名字改变 例如A改为A|B
      cutpoint <- cutpoint[-j] #因为i和j类别已经合并则在分割点列表中去除j类的名字
    }
  }
  else if (type == 'chi'){
    threshold <- qchisq(1 - alpha, s - 1)
    repeat {
      m <- dim(tb)[1] #此时类别总个数
      if (is.null(m) || m < 2) #当此时类别数小于2或目标类别个数小于2时停止
        break
      testlist <- matrix(rep(threshold + 1, m^2), nrow = m) #用于储存两个类别之间的卡方值的矩阵
      for (i in 1 : (m - 1)){ #计算出不同类别两两之间的卡方值并储存
        for (j in (i + 1) : m){
          d <- tb[c(i, j), ] #抽取第i，j类别在不同目标值中的分布 构造列联表
          testlist[i,j] = chiSq(d) #计算出它们之间的卡方值
        }
      }
      k<- min(testlist) #找到最小的卡方值
      if (k > threshold)
        break #当最小的卡方值大于阈值的时候停止合并过程
      if (k <= threshold){ #当最小的卡方值小于等于阈值的时候
        idx <- which(testlist==k, arr.ind = TRUE)[1,] #找到矩阵中等于最小卡方值的位置
        i <- idx[1]
        j <- idx[2] #i，j即为将要合并的类别
        tb[i, ] <- tb[i, ] + tb[j, ] #将i，j两个类进行合并 并将合并后的值赋予类别i
        cutpoint[i] <- paste(cutpoint[i], '|', cutpoint[j]) #把类别i名字改变 例如A改为A|B
        cutpoint <- cutpoint[-j] #因为i和j类别已经合并则在分割点列表中去除j类的名字
        tb <- tb[-j,] #除掉类别j的所有值
      }
    }
  }
  if (!is.null(dim(tb))){
    sumtb <- margin.table(tb, margin = 1)
    combIdx <- which(sumtb < prop * n)
    num <- length(combIdx)
    if (num > 0){ #如果有类别它的样本数小于prop * 总样本数那么合并到合格率最接近的那一类
       proptb <- prop.table(tb, margin = 1)
       preplist <- proptb[, 2]
       preplistc <- proptb[-combIdx, 2]
       prepidx <- 1:nrow(tb)
       prepidxc <- prepidx[-combIdx]
       for (i in combIdx){
           prep <- preplist[i]
           diff <- abs(preplistc - prep)
           idx <- prepidxc[which.min(diff)]
           cutpoint[idx] <- paste(cutpoint[idx], '|', cutpoint[i])
           cutpoint[i] <- NA
       }
       cutpoint <- cutpoint[which(!is.na(cutpoint))]
    }
  }
  return(cutpoint) #返回合并后的分割点列表
}

merg <- function(cuts, x){ #对数据合并后的结果 改变类别的名字
  len <- length(cuts) #合并后的总类数
  data <- as.character(x)
  for (i in 1:len){
    x <- cuts[i] #第i个类别
    y <- trimws(unlist(strsplit(x, '[|]')), which = c('both')) #合并后的分割类别提取出来 如A|B提取为 A B
    if (length(y) >= 2){ #如果类别需要合并
      for (z in y){
        idx <- which(data == z)
        data[idx] <- rep(x, length(idx))
      }
    }
  }
  return(data) #返回改变后的数据
}

# other functions ---------------------------------------------------------
# 插入函数
insert <- function (x, l){ #把元素x插入列表l中并保持表中元素的顺序为升序
  p <- length(l)
  w <- which(l > x)
  len <- length(w)
  if (len == p)
    return(c(x, l))
  if (len == 0)
    return(c(l, x))
  return(c(l[1:(w[1]-1)], x, l[w[1]:p]))
}

#返回区间合格率和区间之间合格率的差值
ratio <- function(discx, y){
  tb <- table(discx, y)
  proptb <- prop.table(tb, margin = 1)
  if (ncol(proptb) > 1){
    yprop <- as.numeric(proptb[, 2]) #合格率
    diff <- abs(diff(yprop)) #合格率差值
  }
  else{
    yprop <- rep(0,2)
    diff <- 0
  }
  num <- as.numeric(table(discx)) #区间样本数
  return (list(yprop = yprop, diff = diff, num = num))
}

tarchange <- function(y){#改变目标变量的值 good==1 bad==0
  ynew <- rep(0, length(y))
  if (length(grep('ok', y, ignore.case = TRUE)) != 0 ){
    idx <- grepl('ok', y, ignore.case = TRUE) & !grepl('bad', y, ignore.case = TRUE)
    ynew[idx] <-1
  }
  else if(length(grep('good', y, ignore.case = TRUE)) != 0){
    idx <- grepl('good', y, ignore.case = TRUE) & !grepl('bad', y, ignore.case = TRUE)
    ynew[idx] <-1
  }
  else stop('no target')
  return(ynew)
}

mylog <- function (x){ #计算log函数
  x[which(x <= 1e-10)] <- 1 #当x足够小的时候定义log(x) = 0
  return(log(x))
}

GetNAidx <- function(x){
  res <- which(sapply(x, function(x) is.na(x) || 
                 is.null(x) || x == "") == TRUE) #找到值为空的坐标
  return (res)  
}



