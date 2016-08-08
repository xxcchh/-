#-----------------------------------------------------------------------------------------------
# 更新时间：2016.8.7
# 作者：Crystal
# 更新内容：
# 1）删除step8，把所有规则放在一起比较相似性。
#-----------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------
# 时间：2016年8月4日8:34:20
# 目的：step5生成规则V4.0
#-----------------------------------------------------------------------------------------------

#主函数(rulesmining)： 
#输入参数接口说明：
#data_pre_handle:step4执行得到的结果，包含numlist、ratiolist、cutlist、discdata四个list的大list
#GainsR_value：step4执行得到的结果
#appCols：用户指定的必须出现在规则中的变量,以向量（数组）的形式传值
#conCols：用户选择的可控变量，以向量（数组）的形式传值
#n：控制选取重要性最大的变量个数，此接口需开放给后台
#mu：设置supp时的参数，此接口需开放给后台
#eta：判断规则相似的参数，此接口需开放给后台
#alpha：显著性水平，此接口需开放给后台
#minNum：第一阶段至少得到的规则数量，此接口需开放给后台
#maxLhs：第三阶段最大前项限制，此接口需开放给后台
#appointed：是否最后得到的规则中必须有指定变量的那几个变量，指定变量数量最多为3，由web前端做限制。默认TRUE(是)，此接口需开放给后台
#theta:step8规则置信度阈值系数，此接口需开放给后台

##########
#输出参数接口说明：
#rulesFinal：最终用于呈现的规则
#ruleslist：每一步得到的具体规则，后台不需调用
#rulesnum：每一步的规则数量，后台只需调用rulesnum[[8]](即得到规则系数量)
#rulesnumPerseries：最后一步得到的每个规则系的规则数量
#puritylossAndCov：最后一步得到的规则系中每条规则每个前项对应的合格率损失和覆盖率
#---------------------------------------------------------------

library(R6)
library(utils)
library(rlist)
library(arules)

# rulesmining -------------------------------------------------------------
rulesmining <- function(data_pre_handle, GainsR_value, appCols, conCols, n=10, mu=0.2, eta=0.5,
                        alpha=0.95, minNum=150, maxLhs=15, appointed=TRUE, theta=0.9){
  
  #GainsR_value为变量重要性排序
  #data_pre_handle 为连续数据离散化之后得到的list结果，包含有discdata ratiolist cutlist numlist
  #appCols为用户指定的必须出现在规则中的变量，当前端不指定变量时需默认传值NULL，且appointed对应需传值FALSE
  #conCols为可控变量，用户必须选择，不能为NULL
  #n控制选取重要性最大的变量个数
  #mu为设置supp时的参数
  #eta为判断规则相似的参数
  #alpha为显著性水平
  #minNum为第一阶段至少得到的规则数量
  #maxLhs为第三阶段最大前项限制
  #appointed为是否最后得到的规则中必须有指定变量的那几个变量，指定变量数量最多为3，由web前端做限制。默认TRUE(是)
  #theta为第三阶段增加前项时候supp的限制阈值，即当supp > rate*0.2*theta时该前项才能被添加入规则
  
  data <- ChoosedData(data_pre_handle,appCols,conCols,GainsR_value,n)
  p <- ncol(data)
  N <- nrow(data)
  y <- data[, p]
  ynew <- tarchange(y)
  target <- colnames(data)[p]
  value <- as.character(unique(y[ynew==1]))
  rate <- length(which(ynew == 1)) / N
  supp <- rate * mu #设置为先验概率*mu
  conf <- rate #设置为先验概率
  minlen <- 2
  maxlen <- 6
  lift <- 1.01
  ruleslist <- list() #记录每一步产生的规则
  rulesnum <- list() #记录每一步产生的规则数目
  ######################################
  # 阶段一: Apriori算法挖掘规则
  res <- GLFX(data, value, supp, conf, minlen, maxlen, lift, alpha, minNum, 
              appCols, conCols, appointed)
  rulesseries <- res$rules
  # 提前结束
  if (is.null(rulesseries)){
    # warning("规则数为0了")
    warning("0")
    return(res)
  }
  #step 1-4
  ruleslist[c(1:4)] <- res$ruleslist
  #step 5
  ruleslist[[5]] <- rulesseries
  rulesnum[c(1:4)] <- res$rules_num
  rulesnum[[5]] <-  sum(unlist(lapply(rulesseries, function(x) length(x))))
  ######################################
  # 阶段二: 筛选出有差异的规则集合
  interestFilter <- function(rules){ #对规则系中的规则计算val并筛选出每个系中val最大的规则
    len <- length(rules) 
    rulesNew <- c() 
    for (i in 1:len){
      r <- rules[[i]]
      inter <- Interest(r@quality$supp, r@quality$conf, rate, mu)
      idx <- which.max(inter)
      res <- r[idx]
      res@quality$interest <- inter[idx]
      if (i == 1){
        rulesNew <- res
      }
      else{
        rulesNew <- union(rulesNew, res)
      }
    }
    return(rulesNew)
  }
  subsetFilter <- function(rules){ #找到仍然为某集合子集的规则并且并入到一起作为规则系选出val最大的规则
    ridx <- 1:length(rules) 
    m <- is.subset(rules, proper = TRUE)
    msum <- rowSums(m)
    subidx <- as.numeric(which(msum >= 1))
    for (i in subidx){
      idx <- as.numeric(which(m[i, ] == TRUE))
      if ((length(idx) == 0) | (ridx[i] == 0))
        next
      r <- rules[i]
      rs <- c(r, rules[idx]) # 子集和其母集的集合
      val <- rs@quality$interest
      repidxl <- which(val[1] < val[-1])
      repidxr <- which(val[1] > val[-1])
      condl <- length(repidxl)
      condr <- length(repidxr)
      if (condr){ # 如果存在母集比子集val小那么删除该母集
        ridx[idx[repidxr]] <- 0
      }
      if (condl){ # 如果子集比母集val小那么删除该子集
        ridx[i] <- 0 
      }
    }
    rulesNew <- rules[ridx[ridx != 0]]
    return(rulesNew)
  }
  #step 6
  rules <- interestFilter(rulesseries) #从每个规则系中挑选规则
  rules <- subsetFilter(rules) #如果仍然有超集存在则形成规则筛选规则
  ruleslist[[6]] <- rules
  n <- length(rules) #此时规则个数
  rulesnum[[6]] <- n
#   lhslen <- c() #规则的前项个数储存列表
#   for (i in 1:n){
#     r <- rules[i]
#     lhslen[i] <- length(r@lhs@data@i)
#   }
#   groupinfo <- grouping(lhslen) #根据前项长度将规则分组信息
#   groupindx <- c(0, attr(groupinfo, 'ends'))
#   sameJudge <- function(rules){ #判断前项长度相同的规则之间是否相似并且保留val最大的规则
#     appColsNum <- length(appCols)
#     reserveIdx <- rep(1, n)
#     for (k in 1:(length(groupindx) - 1)){
#       idx <- groupinfo[(groupindx[k] + 1) : groupindx[k+1]]
#       ntotal <- lhslen[idx[1]]
#       r <- rules[idx]
#       l <- length(r)
#       if (l == 1)
#         next
#       m <- matrix(rep(0, l^2), nrow = l)
#       keepidx <- rep(1, l)
#       for (i in 1:(l-1)){
#         for (j in (i+1):l){
#           ilhs <- r[i]@lhs@data@i 
#           jlhs <- r[j]@lhs@data@i 
#           nsame <- length(intersect(ilhs, jlhs))
#           m[i, j] <- (nsame-appColsNum) / (ntotal-appColsNum)
#         }
#       }
#       sameidx <- which(m >= eta, arr.ind = TRUE)
#       keepidx[unique(c(sameidx))] <- 0
#       lidx <- sameidx[, 1]
#       ridx <- sameidx[, 2]
#       lval <- r[lidx]@quality$interest
#       rval <- r[ridx]@quality$interest
#       diff <- (rval - lval) > 0
#       keepidx[c(ridx[diff], lidx[!diff])] <- 1
#       keepidx[c(ridx[!diff], lidx[diff])] <- 0
#       reserveIdx[idx[keepidx == 0]] <- 0
#     }
#     rulesNew <- rules[which(reserveIdx != 0)]
#     return(rulesNew)
#   }
    sameJudge <- function(r){ 
      n <- length(r)
      if (n == 1)
        return(r)
      appColsNum <- length(appCols)
      idx <- 1:n
      reserveIdx <- rep(1, n)
      m <- matrix(rep(0, n^2), nrow = n)
      for (i in 1:(n-1)){
        for (j in (i+1):n){
          ilhs <- r[i]@lhs@data@i + 1
          jlhs <- r[j]@lhs@data@i + 1
          nsame <- length(intersect(ilhs, jlhs))
          ntotal <- length(union(ilhs, jlhs))
          m[i, j] <- (nsame-appColsNum) / (ntotal-appColsNum)
        }
      }
      sameidx <- which(m >= eta, arr.ind = TRUE)
      if (length(sameidx) == 0)
        return(r)
      reserveIdx[unique(c(sameidx))] <- 0
      lidx <- sameidx[, 1]
      ridx <- sameidx[, 2]
      lval <- r[lidx]@quality$interest
      rval <- r[ridx]@quality$interest
      diff <- (rval - lval) > 0
      reserveIdx[c(ridx[diff], lidx[!diff])] <- 1
      reserveIdx[c(ridx[!diff], lidx[diff])] <- 0
      rulesNew <- r[which(reserveIdx != 0)]
      return(rulesNew)
    }
  #step 7
  rules <- sameJudge(rules)
  ruleslist[[7]] <- rules
  rulesnum[[7]] <- length(rules)
  ######################################
  # 阶段三: 为规则添加前项
  #step 8
  datatemp <- data_pre_handle$discdata
  datatemp <- datatemp[, c(which(colnames(datatemp) %in% 
                                   c(appCols, conCols)), ncol(datatemp))]
  res <- seqcovermining(datatemp, appCols, rules, rate, maxLhs, mu, theta)
#   rules <- res$rules 
#   ruleslist[[8]] <- rules
#   rulesnum[[8]] <- sum(unlist(lapply(rules, function(x) nrow(x))))
  # 删减前项
  #step 9
  rules <- res$prules
  rulesFinal <- lapply(rules, function(x) {x$supp <- x$supp / rate; return(x)}) # 最终呈现的规则
#   ruleslist[[9]] <- rules
#   rulesnum[[9]] <- sum(unlist(lapply(rules, function(x) nrow(x))))
  ruleslist[[8]] <- rules
  rulesnum[[8]] <- sum(unlist(lapply(rules, function(x) nrow(x))))
  # step 10
  # 所有规则所有前项的purityloss和coverage
  df <- res$pdf
  
  rulesnumPerseries <- list()
  rulesnumPerseries <- lapply(rulesFinal,nrow)
  result <- list(rules = rulesFinal, ruleslist = ruleslist, 
                 rulesnum = rulesnum, rulesnumPerseries = rulesnumPerseries, puritylossAndCov = df)
  return(result)
}

# Interest ----------------------------------------------------------------
Interest <- function(supp, conf, lambda, mu){
  # lambda为合格率
  val <- ((conf - lambda) / lambda) * ((supp - lambda*mu) / (lambda*mu))
  return(val)
}

# ChoosedData ----------------------------------------------------------------
#参数说明：data_pre_handle，分箱后的数据(包括4个list)；appCol，用户指定的一定要出现的变量；
#          conCol，可控变量；GainsR_value，重要性排序结果；n，纳入关联分析的变量数
ChoosedData=function(data_pre_handle,appCols,conCols,GainsR_value,n) {
  
  tempData <- data_pre_handle$discdata
  appNum <- length(appCols)
  conNum <- length(conCols)
  totalNum <- appNum + conNum
  colNum <- dim(tempData)[2]
  tarCol <- colnames(tempData)[colNum]
  
  if(appNum > n ) {
    idx <- which((as.character(GainsR_value[, 1])) %in% appCols)
    idx <- idx[1:n]
    colsChoice <- as.character(GainsR_value[idx, 1])
    data <- tempData[, c(colsChoice, tarCol)]
  }
  else if(appNum == n) {
    data <- tempData[, c(appCols, tarCol)]
  }
  else {
    if(totalNum > n) {
      tempNum <- n - appNum
      idx <- which((as.character(GainsR_value[, 1])) %in% conCols)
      idx <- idx[1:tempNum]
      colsChoice <- as.character(GainsR_value[idx, 1])
      data <- tempData[, c(appCols, colsChoice, tarCol)]
    }
    else {
      data <- tempData[, c(appCols, conCols, tarCol)]
    }
  }
  
  return(data)
}

# GLFX --------------------------------------------------------------------
######求关联规则######
GLFX <- function(data, value, sup_def, conf_def, minlen_def,
                 maxlen_def, lift_def, alpha, minNum, 
                 appCols, conCols, appointed) {
  cond <- sapply(data, function(x) length(levels(x)) == 1) #将水平数为1的变量删除（水平数为1对关联分析没有意义，反倒会增加计算量导致内存溢出）
  idx <- which(cond == TRUE)
  if (length(idx) != 0) {
    data <- data[, -idx]
  }
  if (length(value) != 1) {
    value <- paste(value, collapse="&")
  }
  p <- ncol(data)
  y <- data[, p]
  target <- colnames(data)[p]
  temp <- paste(target, "=", value, sep="")
  data <- factorchange(data)
  data_trans <- as(data, "transactions")#将数据转换成稀疏矩阵
  rules <- list() #储存每一步的规则
  num <- c() #储存每一步的规则数目
  
  ######################################
  
  ######################################
  #rules_get_step1(调用关联规则算法)
  rules[[1]] <- apriori(data_trans, parameter = list(support=sup_def, confidence=conf_def, minlen=minlen_def, maxlen=maxlen_def),
                        appearance = list(rhs=c(temp), default="lhs"), control=list(verbose=F))
  if(appointed == TRUE) {
    rules[[1]] <- getAppRules(rules[[1]], appCols)
  }
  num[1] <- length(rules[[1]])
  ######################################
  
  ######################################
  #rules_get_step2(按提升度筛选)
  rules[[2]] <- subset(rules[[1]], subset = lift >= lift_def)
  num[2] <- length(rules[[2]])
  ######################################
  
  ######################################
  #rules_get_step3(用矩阵法去除冗余规则)
  if (nrow(rules[[2]]@quality) >=2 ) {
    #     tempsort <- sort(rules[[2]], by="lift", decreaing=TRUE)#根据关联规则结果中的提升度进行降序排序
    #     redundant <- is.redundant(tempsort, measure='lift')
    #     rules[[3]] <- tempsort[!redundant] #从规则中去除这些列
    tempsort <- sort(rules[[2]], by="lift", decreaing=TRUE)#根据关联规则结果中的提升度进行降序排序
    tempm <- is.subset(tempsort, tempsort, proper=TRUE)#生成一个关联规则的子集矩阵
    tempm[lower.tri(tempm, diag=TRUE)] <- NA#将矩阵对角线以下的元素置为空
    redundant <- colSums(tempm, na.rm=TRUE) >= 1#将子集矩阵中每列元素和大于等于1的列找出来
    rules[[3]] <- tempsort[!redundant]#从规则中去除这些列
  }
  else {
    rules[[3]] <- rules[[2]]
  }
  num[3] <- length(rules[[3]])
  ######################################
  
  ######################################
  #rules_get_step4(用多重统计显著性检验得到重要规则)
  temp <- sort(rules[[3]][is.significant(rules[[3]], data_trans)],
               by='confidence', decreasing = TRUE)
  if (length(temp) > minNum) {
    rules[[4]] <- temp[1:minNum]
    num[4] <- minNum
  }
  else {
    n <- length(temp)
    rules[[4]] <- temp
    num[4] <- n
  }
  ######################################
  # 判断终止
  if (num[4] == 0){
    idx <- which(num > 0)
    res <- list('ruleslist'= rules[idx], 'rulesnum' = num[idx])
    return(res)
  }
  ######################################
  
  ######################################
  
  ######################################
  #rules_get_step5(划分规则)
  rules_series <- list()                    #规则系列表
  rules_raw <- list()
  if(num[4] == 1){
    rules_raw[[1]] <- rules[[4]]
    rules_series <- as(rules[[4]], "data.frame")
  }
  else{
    booleanm <- is.subset(rules[[4]])
    superidx <- as.vector(which(rowSums(booleanm) == 1))
    super_rules <- rules[[4]][superidx]
    sub_rules <- rules[[4]][-superidx]
    len <- length(superidx)
    for(i in 1:len){
      rulesidx <- is.subset(sub_rules, super_rules[i]) #规则系的索引
      ruleslist <- c(super_rules[i], sub_rules[rulesidx]) #一个超集的规则系
      rules_raw[[i]] <- ruleslist
      rulesdf <- as(ruleslist, "data.frame")
      rules_series[[i]] <- rulesdf
    }
  }
  ##############################################
  
  ##############################################
  rules_BP <- as(rules[[4]], "data.frame")
  result <- list('rules' = rules_raw, 'ruleslist' = rules,
                 'rules_BP' = rules_BP, 'rules_num' = num,
                 'rules_series' = rules_series)
  return(result)
}

# getAppRules -------------------------------------------------------------
getAppRules <- function(rules, appCols){
  tempCols <- sub(" ", "", paste(appCols, "="))
  ruleslist <- sapply(tempCols, function(x) subset(rules, lhs %pin% x))
  res <- rules
  n <- length(ruleslist)
  for (i in 1:n){
    res <- intersect(res, ruleslist[[i]])
  }
  return(res)
}

# seqcovermining -------------------------------------------------------------
seqcovermining <- function(data, appCols, rules, rate, maxLhs, mu, theta){
  rulestree <- RulesTree$new(data, appCols, rules, rate, maxLhs, mu, theta)
  # rules <- rulestree$rules
  prules <- rulestree$prules
  df <- rulestree$pdf
  return(list(prules = prules, pdf = df))
}

# rnode ------------------------------------------------------------------
rnode <- R6Class("rnode", #节点
                 lock_class = TRUE,
                 public = list(
                   #变量名称和生长到该节点时规则覆盖的样本数和其中good的样本数
                   name = "",
                   supp = NULL,
                   conf = NULL,
                   #初始化
                   initialize = function(parent, item, supp, conf, isroot=FALSE){
                     self$name <- as.character(item)
                     self$supp <- supp
                     self$conf <- conf
                     private$root <- isroot
                     private$p_parent <- parent
                     private$len <- length(unlist(strsplit(self$name, split = ',')))
                     if (!is.null(parent)){
                       private$len <- private$len + parent$GetLen()
                     }
                   },
                   #增加子节点
                   AddChild = function(child){
                     if (!missing(child)) {
                       private$p_children[[child$name]] <- child
                       child$parent <- self
                     }
                   },
                   # 删除子节点
                   DeleteChild = function(item){
                     child <- self$FindChild(item)
                     if (!is.null(child)){
                       private$p_children[[child$name]] <- NULL
                     }
                     else{
                       warning('no such child')
                     }
                   },
                   #判断是否为叶节点
                   IsLeaf = function(){
                     if (is.null(private$p_children))
                       return(TRUE)
                     else
                       return(FALSE)
                   },
                   #判断是否为根节点
                   IsRoot = function(){
                     private$root
                   },
                   #返回高度
                   GetHeight = function(){
                     height <- 1
                     if (self$IsLeaf())
                       return(height)
                     else{
                       for (child in private$p_children){
                         h <- child$GetHeight()
                         if (h >= height)
                           height <- height + h
                       }
                     }
                     return(height)
                   },
                   #得到前项长度
                   GetLen = function(){
                     if (!is.null(self))
                       return (private$len)
                     else
                       return(0)
                   },
                   #找到所有子节点
                   GetAllChildren = function(){
                     return (private$p_children)
                   },
                   #找到所有节点
                   GetAllNodesList = function(){
                     if (self$IsLeaf()){
                       return (self)
                     }
                     else{
                       val <- list(root = self, children = private$p_children)
                       return (val)
                     }
                   },
                   #根据名字找到子节点
                   FindChild = function(item){
                     itemname <- as.character(item)
                     for (c in private$p_children){
                       if (itemname == c$name)
                         return(c)
                     }
                     return (NULL)
                   },
                   #得到节点名称
                   GetRawName = function(){
                     val <- gsub("=.*", "", self$name)
                     return (val)
                   },
                   #增加属性
                   AddItem = function(fact, itemlist, N, lambda, maxLhs, mu, theta){
                     tb <- self$GetTable(fact, itemlist, N, lambda, mu, theta)
                     idx <- c() #记录增加属性的坐标
                     if (max(tb$fit) != -Inf){
                       idx <- which.max(tb$fit)
                       item <- rownames(tb[idx, ])
                       supp <- tb[idx, 'supp']
                       conf <- tb[idx, 'conf']
                       c <- rnode$new(self, item, supp, conf)
                       self$AddChild(c)
                       fact <- fact[which(fact[, idx] == itemlist[, idx]), ]
                     }
                     else{
                       c <- self
                     }
                     minfidx <- which(colnames(itemlist) %in% rownames(tb[which(tb$fit == -Inf), ]))
                     fact[, c(minfidx, idx)] <- NULL
                     itemlist[, c(minfidx, idx)] <- NULL
                     c$Grow(fact, itemlist, N, lambda, maxLhs, mu, theta)
                   },
                   Grow = function(fact, itemlist, N, lambda, maxLhs, mu, theta){
                     if ((ncol(itemlist) > 0) & (self$GetLen() < maxLhs))
                       self$AddItem(fact, itemlist, N, lambda, maxLhs, mu, theta)
                   },
                   # 构造表格
                   GetTable = function(fact, l, N, lambda, mu, theta){
                     n <- ncol(l)
                     name <- colnames(l)
                     df <- data.frame(conf = rep(0,n), supp = rep(0,n), fit = rep(0, n))
                     p_supp <- self$supp
                     p_conf <- self$conf
                     for (i in 1:n){
                       temp <- fact[which(fact[, i] == l[, i]), (n+1)]
                       count <- length(temp)
                       gcount <- length(which(temp==1))
                       if (count > 0)
                         conf <- gcount / count
                       else
                         conf <- 0
                       supp <- gcount / N
                       df[i, 'supp'] <- supp
                       df[i, 'conf'] <- conf
                       #################################
                       # 变量添加的判断条件
                       cond1 <- Interest(supp, conf, lambda, mu) -
                         Interest(p_supp, p_conf, lambda, mu)
                       cond2 <- supp - lambda*mu*theta
                       # Interest(supp = , conf = , lambda = )
                       # 增加规则项之后的interest的差值
                       if (cond1 <= 0 || cond2 <= 0)
                         df[i, 'fit'] <- -Inf
                       else
                         df[i, 'fit'] <- cond1
                     }
                     rownames(df) <- name
                     return(df)
                   },
                   # 删减属性
                   DeleteItem = function(N, lambda, mu){
                     c <- self
                     itemlist <- list()
                     dflist <- list() #储存每条规则的每个前项purity_loss和coverage信息的list
                     i <- 1
                     while (length(c$children) > 0){
                       res <- c$GetFrame(N, lambda, mu)
                       item <- res$item
                       info <- res$info
                       df <- res$df
                       itemlist[[item]] <- info
                       dflist[[i]] <- df
                       i <- i + 1
                       c$DeleteChild(item)
                     }
                     return(list(itemlist = itemlist, dflist = dflist))
                   },
                   # 找到删除后使得interest最大的某个节点
                   GetFrame = function(N, lambda, mu){
                     children <- self$GetAllChildren()
                     n <- length(children)
                     item <- NULL
                     maxint <- 0
                     suppi <- 0
                     confi <- 0
                     df <- data.frame() #储存purity_loss和coverage
                     confraw <- length(list.common(lapply(children, function(x) x$conf))) / 
                       length(list.common(lapply(children, function(x) x$supp)))
                     for (i in 1:n){
                       c <- children[[i]]
                       if (n > 1){
                         resichildren <- children[-i]
                         temps <- length(list.common(lapply(resichildren, function(x) x$supp)))
                         tempc <- length(list.common(lapply(resichildren, function(x) x$conf)))
                         supp <- tempc / N
                         conf <- tempc / temps
                       }
                       else{
                         temps <- length(c$supp)
                         tempc <- length(c$conf)
                         supp <- tempc / N
                         conf <- tempc / temps
                       }
                       int <- Interest(supp, conf, lambda, mu)
                       if (int >= maxint){
                         maxint <- int
                         item <- c$name
                         suppi <- supp
                         confi <- conf
                       }
                       # 计算purity_loss和coverage
                       name <- c$GetRawName()
                       df[i, 'name'] <- name
                       df[i, 'purityloss'] <- ifelse(n > 1, confraw - conf, confraw)
                       df[i, 'coverage'] <- length(c$supp) / N
                     }
                     res <- list(supp = suppi, conf = confi, interest = maxint)
                     return(list(item = item, info = res, df = df))
                   },
                   # 删减前项最优规则
                   DeleteItemNew = function(idx, N, lambda, mu){
                     children <- self$GetAllChildren()
                     appchildren <- children[idx]
                     itemlist <- list() #储存不同长度规则前项的选择信息
                     infolist <- list() #储存supp，conf，interest信息
                     dflist <- list() #储存每条规则的每个前项purity_loss和coverage信息的list
                     n <- length(children)
                     k <- length(idx)
                     for (m in n:1){
                       # 如果有appCols
                       if (k > 0){
                         if (m > k)
                           res <- self$GetFrameNew(children[-idx], m-k, n-k, N, lambda, mu, appchildren)
                         else
                           res <- self$GetFrameNew(children[idx], m, k, N, lambda, mu)
                       }
                       else{
                         res <- self$GetFrameNew(children, m, n, N, lambda, mu)
                       }
                       i <- n - m + 1
                       itemlist[[i]] <- res$item
                       infolist[[i]] <- res$info
                       dflist[[i]] <- res$df
                     }
                     return(list(itemlist = itemlist, infolist = infolist, dflist = dflist))
                   },
                   # 在删减前项的过程中寻找全局最优的组合
                   GetFrameNew = function(children, m, n, N, lambda, mu, appchildren = NULL){
                     item <- NULL
                     maxint <- 0
                     suppi <- 0
                     confi <- 0
                     df <- data.frame() #储存purity_loss和coverage
                     numlist <- combn(n, m, simplify = FALSE)
                     for (k in numlist){
                       choosechildren <- children[k]
                       if (!is.null(appchildren))
                         choosechildren <- c(appchildren, choosechildren)
                       temps <- length(list.common(lapply(choosechildren, function(x) x$supp)))
                       tempc <- length(list.common(lapply(choosechildren, function(x) x$conf)))
                       supp <- tempc / N
                       conf <- tempc / temps
                       int <- Interest(supp, conf, lambda, mu)
                       
                       if (int >= maxint){
                         maxint <- int
                         item <- k
                         suppi <- supp
                         confi <- conf
                       }
                       
                     }
                     p <- c(appchildren, children[item])
                     mn <- length(p)
                     for (i in 1:mn){
                       c <- p[[i]]
                       if (mn > 1){
                         resic <- p[-i]
                         temps <- length(list.common(lapply(resic, function(x) x$supp)))
                         tempc <- length(list.common(lapply(resic, function(x) x$conf)))
                         conf <- tempc / temps
                       }
                       # 计算purity_loss和coverage
                       name <- c$GetRawName()
                       df[i, 'name'] <- name
                       df[i, 'purityloss'] <- ifelse(mn > 1, confi - conf, confi)
                       df[i, 'coverage'] <- length(c$supp) / N
                     }
                     item <- unlist(sapply(p, function(x) x$name))
                     res <- list(supp = suppi, conf = confi, interest = maxint)
                     return(list(item = item, info = res, df = df))
                   }
                 ),
                 active = list(
                   parent = function(val){
                     if (missing(val)) return(private$p_parent)
                     private$p_parent <- val
                   },
                   children = function(val){
                     if (missing(val)) return(private$p_children)
                     private$children <- val
                   }
                 ),
                 private = list(
                   p_children = NULL,
                   p_parent = NULL,
                   len = 0,
                   root = FALSE
                 )
)

# treebuilder -------------------------------------------------------------
RulesTree  <- R6Class('RulesTree', #产生规则树
                      lock_class = TRUE,
                      public = list(
                        items = NULL,#储存属性合格率高的区间列表
                        tree = NULL,#储存tree
                        rules = NULL,#储存规则
                        ptree = NULL, #储存Paralleltree
                        prules = NULL, #储存逐步删除前项的规则系
                        pdf = NULL, #储存每条规则的所有前项的purityloss和coverage
                        initialize = function(data, appCols, rules, prep, maxLhs, mu, theta){
                          res <- ratiochange(data)
                          self$items <- res$itemlist
                          private$data <- res$data
                          private$prep <- prep
                          private$maxLhs <- maxLhs
                          private$theta <- theta
                          private$mu <- mu
                          private$ratio <- res$ratio
                          private$initRules <- rules
                          private$itemLabels <- itemLabels(rules)
                          private$N <- dim(private$data)[1]
                          # self$tree <- self$GrowTree()
                          # self$rules <- self$GetRules()
                          self$rules <- self$GetRawRules()
                          self$ptree <- self$GrowParallelTree()
                          res <- self$GetHierarchyRules(appCols)
                          self$prules <- res$df
                          self$pdf <- res$pdf
                        },
                        #生成树
                        GrowTree = function(){
                          rules <- private$initRules
                          n <- length(rules)
                          tree <- rnode$new(NULL, "Rulestree", 0, 0, isroot = TRUE)
                          for (i in 1:n){
                            r <- rules[i]
                            node <- self$GetNode(r)
                            tree$AddChild(node)
                          }
                          return(tree)
                        },
                        #生成平行树
                        GrowParallelTree = function(){
                          rules <- self$rules
                          data <- private$data
                          n <- length(rules)
                          ptree <- list()
                          for (i in 1:n){
                            tree <- rnode$new(NULL, 'PRulestree', 0, 0, isroot = TRUE)
                            r <- rules[[i]]
                            itemlist <- unlist(strsplit(r$rule, ','))
                            res <- strsplit(itemlist, '=')
                            attr <- as.character(unlist(lapply(res, function(x) x[1]))) #元素的属性
                            val <- as.character(unlist(lapply(res, function(x) x[2]))) #元素的属性值
                            m <- length(itemlist)
                            for (j in 1:m){
                              c <- itemlist[j]
                              supp <- which(data[, attr[j]] == val[j])
                              conf <- which(data[, attr[j]] == val[j] & data[, ncol(data)] == 1)
                              node <- rnode$new(NULL, c, supp, conf)
                              tree$AddChild(node)
                            }
                            ptree[[i]] <- tree
                          }
                          return(ptree)
                        },
                        # 得到节点
                        GetNode = function(r){
                          res <- self$GetRuleInfo(r)
                          item <- res$item
                          fact <- res$fact
                          itemlist <- res$itemlist
                          supp <- r@quality$support
                          conf <- r@quality$confidence
                          node <- rnode$new(NULL, item, supp, conf)
                          node$Grow(fact, itemlist, private$N, private$prep, private$maxLhs, private$mu, private$theta)
                          return(node)
                        },
                        # 从某规则中得到需要考虑的剩下属性信息
                        GetRuleInfo = function(rule){
                          elem <- private$itemLabels[rule@lhs@data@i + 1]
                          n <- length(elem)
                          res <- strsplit(elem, split = '=')
                          attr <- as.character(unlist(lapply(res, function(x) x[1]))) #元素的属性
                          val <- as.character(unlist(lapply(res, function(x) x[2]))) #元素的属性值
                          idxlist <- list()
                          for (k in 1:n){
                            idxlist[[k]] <- which(private$data[, attr[k]] == val[k])
                          }
                          idx <- list.common(idxlist)
                          colsdata <- which(!colnames(private$data) %in% attr)
                          colsitem <- which(!colnames(self$items) %in% attr)
                          fact <- private$data[idx, colsdata]
                          itemlist <- self$items[, colsitem]
                          name <- paste(elem, collapse = ",")
                          return(list(item = name, fact = fact, itemlist = itemlist))
                        },
                        # 得到规则
                        GetRules = function(){
                          tr <- self$tree
                          n <- length(private$initRules)
                          rules <- list()
                          for (i in 1:n){
                            node <- tr$GetAllChildren()[[i]]
                            name <- node$name
                            df <- data.frame()
                            len <- node$GetHeight()
                            if (len > 1){
                              for (j in 1:(len - 1)){
                                node <- node$GetAllChildren()[[1]]
                                name <- paste0(name, ',', node$name, '=', self$items[, node$name])
                              }
                            }
                            df[1, 'rule'] <- name
                            df[1, 'supp'] <- node$supp
                            df[1, 'conf'] <- node$conf
                            df[1, 'interest'] <- Interest(node$supp, node$conf, private$prep, private$mu)
                            rules[[i]] <- df
                          }
                          return(rules)
                        },
                        # 对原始规则进行操作
                        GetRawRules = function(){
                          n <- length(private$initRules)
                          rules <- list()
                          for (i in 1:n){
                            df <- data.frame()
                            rule <- private$initRules[i]
                            elem <- private$itemLabels[rule@lhs@data@i + 1]
                            df[1, 'rule'] <- paste(elem, collapse = ',')
                            df[1, 'supp'] <- rule@quality$supp
                            df[1, 'conf'] <- rule@quality$conf
                            df[1, 'interest'] <- rule@quality$interest
                            rules[[i]] <- df
                          }
                          return(rules)
                        },
                        #逐步删减前项得到规则
                        GetParallelRules = function(){
                          tr <- self$ptree
                          rules <- self$rules
                          namelist <- list()
                          df <- list()
                          pdf <- list()
                          n <- length(tr)
                          for (i in 1:n){
                            ri <- rules[[i]]
                            tri <- tr[[i]]
                            dfi <- ri
                            res <- tri$DeleteItem(private$N, private$prep, private$mu)
                            lhslist <- res$itemlist
                            pdfi <- res$dflist
                            itemi <- unlist(strsplit(ri$rule, ','))
                            namelist <- names(lhslist)
                            m <- length(lhslist)
                            if (m > 1){
                              for (j in 1:(m-1)){
                                namej <- namelist[j]
                                idx <- which(itemi %in% namej)
                                itemi <- itemi[-idx]
                                rulej <- paste(itemi, collapse = ',')
                                lhsj <- lhslist[[j]]
                                idxj <- j + 1
                                dfi[idxj, 'rule'] <- rulej
                                dfi[idxj, 'supp'] <- lhsj$supp
                                dfi[idxj, 'conf'] <- lhsj$conf
                                dfi[idxj, 'interest'] <- lhsj$interest
                              }
                            }
                            df[[i]] <- dfi
                            pdf[[i]] <- pdfi
                          }
                          return(list(df = df, pdf = pdf))
                        },
                        #删减前项得到全局最优解
                        GetHierarchyRules = function(appCols){
                          tr <- self$ptree
                          rules <- self$rules
                          appColsIdx <- self$GetAppIdx(appCols)
                          df <- list()
                          pdf <- list()
                          n <- length(tr)
                          for (i in 1:n){
                            dfi <- data.frame()
                            tri <- tr[[i]]
                            appColsIdxi <- appColsIdx[[i]]
                            res <- tri$DeleteItemNew(appColsIdxi, private$N, private$prep, private$mu)
                            itemlisti <- res$itemlist
                            infolisti <- res$infolist
                            pdfi <- res$dflist
                            m <- length(itemlisti)
                            for (j in 1:m){
                              rulej <- itemlisti[[j]]
                              rulej <- paste(rulej, collapse = ',')
                              lhsj <- infolisti[[j]]
                              idxj <- j 
                              dfi[idxj, 'rule'] <- rulej
                              dfi[idxj, 'supp'] <- lhsj$supp
                              dfi[idxj, 'conf'] <- lhsj$conf
                              dfi[idxj, 'interest'] <- lhsj$interest
                            }
                            df[[i]] <- dfi
                            pdf[[i]] <- pdfi
                          }
                          return(list(df = df, pdf = pdf))
                        },
                        # 得到appCols的坐标
                        GetAppIdx = function(appCols){
                          tr <- self$ptree
                          appColsIdx <- list()
                          n <- length(tr)
                          for (i in 1:n){
                            tri <- tr[[i]] 
                            childi <- tri$GetAllChildren()
                            choicei <- as.character(sapply(childi, function(x) x$GetRawName()))
                            appColsIdxi <- which(appCols %in% choicei)
                            appColsIdx[[i]] <- appColsIdxi
                          }
                          return (appColsIdx)
                        }
                      ),
                      private = list(
                        data = NULL,
                        ratio = NULL,
                        N = NULL,
                        initRules = NULL,
                        itemLabels = NULL,
                        prep = 0,
                        maxLhs = 0,
                        theta = 0,
                        mu = 0
                      )
)

# ratiochange -------------------------------------------------------------------
ratiochange <- function(data){
  p <- dim(data)[2] - 1
  datanew <- data
  datanew[, (p+1)] <- tarchange(datanew[, (p+1)])
  itemlist <- list()
  ratio <- list()
  for (i in 1:p){
    tb <- table(datanew[, i], datanew[, (p+1)])
    proptb <- prop.table(tb, margin = 1)
    if (ncol(proptb) > 1){
      if (nrow(proptb) > 1)
        item <- names(which.max(proptb[, 2]))
      else
        item <- rownames(proptb)
      yprop <- as.numeric(proptb[, 2])
    }
    else{
      item <- rownames(proptb)[1]
      yprop <- rep(0, nrow(proptb))
    }
    itemlist[[i]] <- item
    ratio[[i]] <- yprop
  }
  itemlist <- data.frame(itemlist, stringsAsFactors = FALSE)
  colnames(itemlist) <- colnames(datanew[, -(p+1)])
  return(list(itemlist = itemlist, ratio = ratio, data = datanew))
}

# factorchange ------------------------------------------------------------
factorchange <- function(data){
  cols <- colnames(data)
  data[, cols] <- lapply(data[, cols], factor)
  return(data)
}

# tarchange ---------------------------------------------------------------
tarchange <- function(y) {#good==1 bad==0
  ynew <- rep(0, length(y))
  if (length(grep('ok', y, ignore.case = TRUE)) != 0 ) {
    idx <- grepl('ok', y, ignore.case = TRUE) & !grepl('bad', y, ignore.case = TRUE)
    ynew[idx] <-1
  }
  else if(length(grep('good', y, ignore.case = TRUE)) != 0) {
    idx <- grepl('good', y, ignore.case = TRUE) & !grepl('bad', y, ignore.case = TRUE)
    ynew[idx] <-1
  }
  else stop('no target')
  return(ynew)
}

# StratSample -------------------------------------------------------------
StratSample <- function(data, rate=0.2){ #分层抽样1 <- data[which(y == 1), ]
  group0 <- data[which(y == 0), ]
  n1 <- nrow(group1)
  n0 <- nrow(group0)
  idx1 <- sample(1:n1, round(n1*rate), replace = FALSE)
  idx0 <- sample(1:n0, round(n0*rate), replace = FALSE)
  train1 <- group1[idx1, ]
  train0 <- group0[idx0, ]
  train <- rbind(train1, train0)
  train <- train[order(as.integer(rownames(train))), ]
  return(train)
}

