paritiallypaired <- function(Sample1, Sample2, alternative = "two.sided", mu = 0, conf = 0.95){
  newdata <-data.frame(Sample1, Sample2)
  
  #data 1 (n1)
  data1 <- na.omit(newdata)
  data1
  
  #data 2 (n2)
  data2__1<-newdata[is.na(newdata$Sample2),]
  data2 <-data2__1[is.na(data2__1$Sample1)==FALSE,]
  data2
  
  #data 3 (n3)
  data3__1<-newdata[is.na(newdata$Sample1),]
  data3 <-data3__1[is.na(data3__1$Sample2)==FALSE,]
  data3
  
  #data 4 (n4)
  data4_1 <- newdata[is.na(newdata$Sample1),]
  data4 <- data4_1[is.na(data4_1$Sample2),]
  data4
  
  n = nrow(newdata)
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  n3 <- nrow(data3)
  n4 <- nrow(data4)

  result <- 0
  
  if(n4 == n){ #n1,n2,n3 ==0
    # no data
    print("Invalid")
    return
  }else if (n2 == n){#n1=0, n2>0, n3=0
    #If data is normal -> one sample t-test
    if (shapiro.test(data2$Sample1)$p.value>0.05){
      result <- t.test(data2$Sample1, alternative = alternative , mu = mu, conf.level = conf )

    }
    #if data is not normal -> Wilcoxon signed rank test
    else{
      result <-wilcox.test(data3$Sample2,alternative = alternative , mu = mu, conf.level = conf)
    }
  }else if(n3 ==n){ # n1=0, n2=0, n3>0
  #If data is normal -> one sample t-test
    if (shapiro.test(data3$Sample2)$p.value>0.05){
    result <-t.test(data3$Sample2,alternative = alternative , mu = mu, conf.level = conf)
    }
  #if data is not normal -> Wilcoxon signed rank test
    else{
    result <-wilcox.test(data3$Sample2,alternative = alternative , mu = mu, conf.level = conf)
    }
  }else if(n2 + n3 == n){
    # To check two data's normality -> Shapiro wilk
    # if data is normal -> variance test
    if (shapiro.test(data2$Sample1)$p.value > 0.05 & shapiro.test(data3$Sample2)$p.value > 0.05){
      #if data is normal 
      if(var.test(data2$Sample1,data3$Sample2)$p.value >0.05){
        # equal variance -> pooled 2 sample t test
        result <- t.test(data2$Sample1,data3$Sample2,alternative = alternative , mu = mu, conf.level = conf)
      }
      else{
        result <- t.test(data2$Sample1,data3$Sample2,alternative = alternative , mu = mu, conf.level = conf)
      }
    }else{# if data is not normal -> two sample wilcoxon test
      result <- wilcox.test(data2$Sample1,data3$Sample2,alternative = alternative , mu = mu, conf.level = conf)
    }
  }else if (n1 == n){
    # To check normality -> Shapiro Wilk test
    diff <-data1$Sample1 - data1$Sample2
    #if data is normal -> two sample t-test
    if (shapiro.test(diff)$p.value > 0.05){
      result<- t.test(data1$Sample1,data1$Sample2, paired = T,alternative = alternative , mu = mu, conf.level = conf)
    }
    # if data is not normal -> wilcox test
    else{
      result<- wilcox.test(data1$Sample1,data1$Sample2, paired =T,alternative = alternative , mu = mu, conf.level = conf)
    }
  }else if (n1 != 0 & n2 !=0 & n3 !=0){ # n1,n2,n3 >0
    
    # variables
    D_ <- mean(data1$Sample1 - data1$Sample2)
    T_ <- mean(data2$Sample1)
    N_ <- mean(data3$Sample2)
    SD <- sd(data1$Sample1 - data1$Sample2)
    ST <- sd(data2$Sample1)
    SN <- sd(data3$Sample2)
    data23 <- c(data2$Sample1,data3$Sample2)
    nh <- 1/mean(1/data23)
    
    #Kim's et al.s modified t statistic
    t3num <- ((n1*D_)+(nh*(T_-N_)))
    t3den <- (sqrt(n1*SD^2)+nh^2*((ST^2/n2)+(SN^2/n3)))
    t3 <- t3num/t3den
    if (alternative == "greater"){
      pvalue <- 1 - pnorm(t3, mean = 0, sd = 1);
      z <- qnorm(conf)
    }else if (alternative == "less"){
      pvalue <- pnorm(t3,mean= 0, sd=1)
      z <- qnorm(conf)
    }else if(alternative == "two.sided"){
      pvalue <- 2 * (1- pnorm(abs(t3), mean= 0, sd=1))
      z <- qnorm(1- conf)
    }else{
      print("Invalid")
      return ()
    }
    meandifference <- n1 * (D_-mu) + nh*((T_-N_))
    upper <- meandifference + abs(z) * t3den
    lower <- meandifference - abs(z) * t3den
    conf.interval <- c(lower,upper)
    
    result <- c(pvalue, conf.interval)
    return(result)
  }else{ # n1 > 0, n2>0, n3=0 / n1>0, n2 =0, n3>0
    data1_1 <- data1$Sample1
    data2_1 <- data2$Sample1
    data12_1 <- c(data1_1,data2_1)
    T_2 <- mean(data12_1)
    data1_2 <- data1$Sample2
    data3_2 <- data3$Sample2
    data13_2 <- c(data1_2,data3_2)
    N_2 <- mean(data13_2)
    ST_2 <- sd(data12_1)
    SN_2 <- sd(data13_2)
    STN <- cov(data1$Sample1, data1$Sample2)
    
    numerator_zcorr <-(T_2 - N_2)
    denominator_zcorr <- (sqrt(ST_2^2/(n1+n2) + SN_2^2/(n1+n3) - (2*n1*STN)/((n1+n2)*(n1+n3))))
    zcorr<- numerator_zcorr/denominator_zcorr
    
    if (alternative == "greater"){
      pvalue <- 1 - pnorm(zcorr, mean = 0, sd = 1);
      z <- qnorm(conf)
    }else if (alternative == "less"){
      pvalue <- pnorm(zcorr,mean= 0, sd=1)
      z <- qnorm(conf)
    }else if(alternative == "two.sided"){
      pvalue <- 2 * (1- pnorm(abs(zcorr), mean= 0, sd=1))
      z <- qnorm(1- conf)
    }else{
      print("Invalid")
      return()
    }
    upper <- (numerator_zcorr-mu) + abs(z) * denominator_zcorr
    lower <- (numerator_zcorr-mu) - abs(z) * denominator_zcorr
    conf.interval<- c(lower, upper)
    
    result <- c(pvalue, conf.interval)
    return(result)
  }
  
  return(result)
}











