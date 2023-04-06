rm(list = ls())
#example using MH sampling method, in this case π(x) = exp(-x), P(x prime| x) = 0.5exp(-0.5x)
set.seed(1)
n <- 323700 
burnin <- 0#burn-in 


chain <- rep(NA,323700)
x <- 0.1 #First sample chain is staring(default) value
chain[1] <- x

for (i in 2:n) {
  
  proposed_x <- rexp(1,rate = 0.5)
  aprob <- min(1, (dexp(proposed_x, rate = 1)/dexp(x, rate = 1)*dexp(x, rate = 0.5)/dexp(proposed_x, rate = 0.5))) 
  #acceptance rate, where:
  
  u <- runif(1)
  if (u < aprob) 
    x <- proposed_x
  chain[i] <- x
}
if (burnin == 0){
  chain=as.matrix(chain)
}
if (burnin != 0){
  chain=as.matrix(chain[-(1:burnin)]) #discord the burn-in part
}


summary(chain)
var(chain)

library("ggplot2")
chaind <-as.data.frame(chain)
p <- ggplot(chaind, aes(x = chain))
p+ geom_density(color = "black", fill = "gray")

library("mcmcse")
estvssamp(chain)

par(mfrow=c(2,1)) 
plot(ts(chain), xlab="Chain", ylab="Draws")
abline(h = mean(chain), lwd="2", col="red" )
acf(chain)



set.seed(2)
n <- 323700 
burnin <- 0#burn-in 


chain2 <- rep(NA,323700)
x <- 0.1 #First sample chain is staring(default) value
chain2[1] <- x

for (i in 2:n) {
  
  proposed_x <- rexp(1,rate = 5)
  aprob <- min(1, (dexp(proposed_x, rate = 1)/dexp(x, rate = 1)*dexp(x, rate = 5)/dexp(proposed_x, rate = 5))) 
  #acceptance rate, where:
  
  u <- runif(1)
  if (u < aprob) 
    x <- proposed_x
  chain2[i] <- x
}

if (burnin == 0){
  chain2=as.matrix(chain2)
}
if (burnin != 0){
  chain2=as.matrix(chain2[-(1:burnin)]) #discord the burn-in part
}

# rho2 <- stats::acf(x = chain2, lag.max = 323699, plot = FALSE)$acf
# sum(rho2)
# [1] 0.5

chiand2 <- as.data.frame(chain2)
library("ggplot2")
p <- ggplot(chaind2, aes(x = chain2))
p+ geom_density(color = "black", fill = "gray")

library("mcmcse")
estvssamp(chain2)

par(mfrow=c(2,1)) 
plot(ts(chain2), xlab="Chain", ylab="Draws")
abline(h = mean(chain2), lwd="2", col="red" )
acf(chain2)


#ESS
#ESS by package "mcmcse"
multiESS(chaind)
#16181.71
minESS(1)
#6146 


#ESS function by myself

#Minnimum ESS for different p, alpha and epsilon 
min.ESS <- function(p, alpha = 0.05, epsilon = 0.05){
  ess <- (2^(2/p)*pi/(p*gamma(1/2))^2/p)*(qchisq(1-alpha,p)/epsilon^2)
  return(ess)
}
min.ESS(2, epsilon = 0.05)
# 153658.4
min.ESS(3)


#Calculate mcse using method "batchmeans"
mcse_batchmeans <- function(x, size = min(100, floor(length(x)/5)), trim = 0, progress = FALSE){
  
  # Set number of batches and calculate batch means
  n <- floor(length(x) / size)
  batch.means <- matrix(NA, n, 1)
  for(i in 1:n){
    batch.means[i] <- mean(x[(1+size*(i-1)):(size*i)])
  }
  
  # Trim the batch means and estimate the standard deviation
  trimmed.batch.means <- batch.means[((1+trim):(n-trim))]
  s <- sd(trimmed.batch.means) / sqrt(length(trimmed.batch.means))
  
  # If progress argument is TRUE, show progress bar
  if(progress){
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    on.exit(close(pb))
    setTxtProgressBar(pb, 1)
  }
  
  return(s)
}

# revised version
get.ESS <- function(x, ESS = 3) {
  # Center the data
  x <- x - mean(x)
  # Compute the length of the vector
  M <- length(x)
  # Compute the autocorrelations if necessary
  if (ESS %in% c(1, 4)) {
    rho <- stats::acf(x = x, lag.max = M - 2, plot = FALSE)$acf
  }
  # Compute the variance and autocovariance
  lambda.sq <- var(x)
  if(ESS == 1){           # file:///C:/Users/14817/Desktop/Convergence%20for%20MCMC/Convergence%20diagnostics%20for%20MCMC.pdf
    E <- M / (1 + 2 * sum(rho))
  }
  else if (ESS == 2) { # file:///C:/Users/14817/Desktop/Convergence%20for%20MCMC/mcmcse_vignette.pdf
    sigma.sq <- mcse.multi(x)$cov
  } 
  else if (ESS == 3) {  # n*lambda/sigma = length(chain)*var(chain)/(mcse_batchmeans(chain))^2*length(chain) 
    sigma.sq <- (mcse_batchmeans(x))^2 * M # it is supposed to have the same value in "ESS == 2", but it seems like there are some errors
  } 
  else if(ESS == 4) { # https://stats.stackexchange.com/questions/429470/what-is-the-correct-effective-sample-size-ess-calculation
    sigma.sq <- lambda.sq + 2 * sum(rho)
  }
  # Compute the effective sample size
  E <- M * lambda.sq / sigma.sq
  return(E)
}



get.ESS(chain)

get.ESS(chain,1)
#161849.3 the most closed value to the paper.
get.ESS(chain,2)
#154088.2
get.ESS(chain,3)
#155330.8
get.ESS(chain,4)
#162173


effectiveSize(chain)
#155823.8


get.ESS(chain2,1)
#161849 /// 4.783389
get.ESS(chain2,2)
#84.15389
get.ESS(chain2,3)
#3557.791
multiESS(chaind2)
#84.15389


library("coda")
effectiveSize(chain2)
#225




rho <- stats::acf(x = chain2, lag.max = 327000-2, plot = FALSE)$acf
sum(rho)
#0.5


# Geweke

# In the function, I calculate "Geweke's statistic", which is used to test the convergence of the Markov chain. The calculation method of 
# this statistic is to take the first a*n samples and the last b*n samples of the chain as two independent sample sets, calculate the 
# standardized value of the mean difference of the two sample sets, and then use the two-sided t distribution Tests whether this 
# normalized value is significantly different from zero.
# In the code, I provided two ways to calculate the sample variance, one is the var() function, and the other is the power spectral 
# density estimate calculated by the spectrum0.ar() function. Users should choose which method to use for calculation through the method 
# parameter.

# 在函数中，我计算了 "Geweke's statistic" 统计量，用于检验马尔科夫链的收敛性。这个统计量的计算方法是将链的前 a*n 个样本和后 b*n 个样本分
# 别作为两个独立的样本集，计算两个样本集的均值差的标准化值，然后用双侧 t 分布检验这个标准化值是否显著不等于 0。
# 在代码中，我提供了两种计算样本方差的方法，一种是 var() 函数，另一种是 spectrum0.ar() 函数计算的功率谱密度估计值。使用者应通过 method
# 参数来选择使用哪种方法进行计算。


# This is an example function that only outputs the statistics in the geweke test once, that is to say, regardless of batch.
# By default we calculate the difference between first a = 10% part and last b = 50% part.

# 这是一个示例函数，他只输出一次Geweke检验中的统计量，也就是说不分batch多次计算.默认情况下我们计算前10%和后50%的差异.

geweke_toy <- function(chain, a = 0.1, b = 0.5, method = "coda") {
  n <- length(chain)
  A <- floor(a * n)
  B <- floor(b * n) 
  # Calculating Geweke’s statistic
  if (method == "normal"){
    ssa <- var(chain[1:A])
    ssb <- var(chain[B:n-1])
    z <- (mean(chain[1:A]) - mean(chain[B:n-1])) / 
      sqrt(ssa / length(chain[1:A]) + ssb / length(chain[B:n-1]))
  }
  
  else if (method == "coda"){
    # This ssa and ssb are the another calculating method of variance (The spectral density at frequency zero) which are using in "coda"
    # Specifically, the spectrum0.ar function performs an autoregressive fit on time series data and computes the variance of the 
    # fitted residuals. This variance can be used as the variance estimate required in the calculation of the Geweke statistic.
    
    # 此处的ssa和ssb是在“coda”中使用的另一种方差计算方法（零频率处的频谱密度）
    # 具体而言，spectrum0.ar函数会对时间序列数据进行自回归拟合，并计算拟合残差的方差。这个方差可以作为Geweke统计量的计算中所需的方差估计。
    ssa <- spectrum0.ar(chain[1:A])$spec
    ssb <- spectrum0.ar(chain[B:n-1])$spec
    z <- (mean(chain[1:A]) - mean(chain[B:n-1])) / 
      sqrt(ssa / length(chain[1:A]) + ssb / length(chain[B:n-1]))
  }
  
  p_value <- 2 * pmin(pnorm(abs(z), lower.tail = FALSE), pnorm(abs(z), lower.tail = TRUE))
  result <- list(z, p_value)
  names(result) <- c("Geweke’s statistic", "p_value")
  return(result)
}

geweke_toy(chain,method = "coda")   
# $`Geweke’s statistic`
# [1] 0.9614776
# $p_value
# [1] 0.3363121

geweke_toy(chain2,method = "coda")
# $`Geweke’s statistic`
# [1] 0.7455853
# $p_value
# [1] 0.455918



# Inside the function, the input time series matrix is first converted to a numeric matrix. Then, for each column, fit it to a linear 
# regression model and check if its residuals are 0. If the residual is 0, it means that the column can be regarded as a constant 
# sequence, its spectral density is 0, and the order of the autoregressive model is also 0; otherwise, the AR model is used to fit the 
# column, and its spectral density and autoregressive model are calculated. The order of the regression model, saved into the 
# corresponding vector. Finally, the spectral density vector and order vector are returned as a list as the result.

# spectrum0.ar会首先将输入的时间序列矩阵转换为一个数值矩阵。然后，对于每一列，将其拟合为一个线性回归模型并检查其残差是否为0。如果残差为
# 0，说明该列数据中所有值都相同（常数序列），其谱密度为0，自回归模型的阶数也为0；否则，使用AR模型拟合该列，计算出其谱密度和自回归模型的
# 阶数，并保存到相应的向量中。最后，将谱密度向量和阶数向量作为结果以列表的形式返回。

spectrum0.ar <- function (x) 
{
  x <- as.matrix(x)
  v0 <- order <- numeric(ncol(x))
  names(v0) <- names(order) <- colnames(x)
  z <- 1:nrow(x)
  for (i in 1:ncol(x)) {
    lm.out <- lm(x[, i] ~ z)
    if (identical(all.equal(sd(residuals(lm.out)), 0), TRUE)) {
      v0[i] <- 0
      order[i] <- 0
    }
    else {
      ar.out <- ar(x[, i], aic = TRUE)
      v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
      order[i] <- ar.out$order
    }
  }
  return(list(spec = v0, order = order))
}

# note:
# 谱估计是一种用于估计信号的频率成分的方法。在时域中，信号的波形随时间变化，而在频域中，信号可以表示成不同频率成分的叠加。谱估计的目的就
# 是从信号的样本数据中估计其在不同频率上的成分大小。

# 在上述代码中，使用自回归模型来估计信号的谱。自回归模型是一种用于建模时间序列的模型，它假设当前时刻的数据值是前几个时刻的数据值的线性组
# 合。在这个模型中，自回归系数表示不同时刻的数据之间的相关性，而谱估计值可以表示为自回归系数的函数。

# 具体来说，ar.out$var.pred表示使用自回归模型预测当前时刻的数据时的方差，而sum(ar.out$ar)表示自回归系数之和，其值通常在-1和1之间。因此，
# ar.out$var.pred/(1 - sum(ar.out$ar))^2这个公式可以将自回归系数的影响从谱估计中去除，得到更准确的谱估计值。

# 关于谱估计值，它是一种频率域上的幅度值，可以用于表示信号在不同频率上的能量分布情况。在谱估计中，通常使用周期图法或基于模型的方法（如自
# 回归模型）来估计谱估计值。

# 值得注意的是，这里使用的是一种特殊的谱估计方法，即零频率处的谱估计方法，它对应于自回归模型中的方差估计。在这个函数中，使用谱估计值来
# 估计方差是合理的，因为Geweke统计量的计算需要用到两个时间窗口中的方差估计值。

# 为了计算Geweke统计量，需要分别估计这两个时间窗口中的方差，并计算它们的加权平均值。在这个函数中，使用了谱估计值来计算这两个时间窗口中的
# 方差估计值。具体而言，使用了spectrum0.ar函数来计算自回归模型在零频率处的谱估计值，然后将它们分别用于估计前a%和后b%的方差。最后，通过对
# 这两个方差的加权平均值来计算Geweke统计量。这是因为零频率处的谱估计值可以反映时间序列数据中的整体方差，而且可以将自回归系数的影响从谱估
# 计中去除，从而更为准确地估计总体方差。因此，将零频率处的谱估计值作为方差估计值用于Geweke统计量的计算是合理的。

# 自回归谱估计方法（如上面代码中的 spectrum0.ar 函数）可以通过考虑数据中的自相关结构来估计谱密度函数，这种方法在信号处理和时间序列分析中
# 经常使用。这种方法的优点是可以充分利用时间序列数据中的自相关结构，因此可以更准确地估计信号的频率成分。

# 在收敛效果不好的情况下，样本方差可能会不准确，因为它没有考虑到时间序列数据中的自相关结构。因此，使用零频率的谱估计值作为方差估计值可能
# 会更好，因为它能够在某种程度上考虑到时间序列数据的自相关结构，特别是在自相关系数较大的情况下。

# 但是，自回归谱估计并不总是比直接计算样本方差更准确。在某些情况下，直接计算样本方差可能是更准确的，例如，当数据中存在离群值或异常值时--
# 自回归谱估计方法假设时间序列是平稳的，并且在该假设下进行谱估计。然而，存在离群值或异常值时，时间序列的平稳性可能会受到影响。在这种情况
# 下，使用自回归谱估计方法估计的谱密度可能会失真，因此可能不可靠。更具体地说，离群值或异常值可能导致时间序列的方差或自协方差结构发生变化，
# 从而违反了平稳性假设。这可能会导致自回归模型的拟合效果变差，从而导致自回归谱估计方法估计的谱密度失真。

# 因此，在选择方差估计方法时需要考虑数据的特点和应用场景。




# A more complex approach

# When calculating the Geweke statistic, the purpose of dividing the sampling sequence into multiple sub sequences is to increase the
# sample size and improve the accuracy of the statistical test. If only the entire sampling sequence is used for the Geweke test, May 
# lead to inaccurate test results, especially if sampling autocorrelation is high. By dividing the sampling sequence into multiple sub 
# sequences, we can obtain more samples and reduce the influence of sampling autocorrelation on the Geweke test. Simultaneously, using 
# multiple sub sequences for the Geweke test can also provide more information in order to identify possible problems or biases in MCMC 
# sampling.

# 在计算Geweke统计量时，将采样序列分成多个子序列是为了增加样本量，提高统计检验的准确性。如果只使用整个采样序列进行Geweke检验，可能会导致
# 检验结果不准确，特别是在采样自相关性很高的情况下。通过将采样序列分成多个子序列，我们可以获得更多的样本数量，并减小采样自相关性对Geweke
# 检验的影响。同时，使用多个子序列进行Geweke检验也可以提供更多的信息，以便识别MCMC采样中可能存在的问题或偏差。

# 于是我编写了一个使用多个batches的新函数，该函数的输入参数包括：

# x：一个包含MCMC采样过程中的所有状态变量的向量。
# a：表示前期采样占总采样长度的比例，默认为0.1。
# b：表示后期采样占总采样长度的比例，默认为0.5。
# num_batches：表示将采样过程分成的批次数量，默认为20。
# method：表示计算Geweke统计量时使用的方法，默认为"normal"。可选的值包括"normal"和"coda"。
# 如果method为"normal"，则使用常规方法计算Geweke统计量。如果method为"coda"，则使用coda软件包提供的spectrum0.ar函数计算Geweke统计量。

geweke <- function(x, a = 0.1, b = 0.5, num_batches = 20, method = "normal"){
  n <- length(x)
  A <- floor(a * n)
  B <- floor(b * n)
  batch_size <- floor((n - B) / num_batches)
  
  # Initialize vector of results
  z <- numeric(num_batches)
  p_value <- numeric(num_batches)
  
  if (method == "normal"){
    for (i in 1:num_batches) {
      # Calculate start and end indices for current batch
      start <- B + (i - 1) * batch_size
      end <- start + batch_size
      
      # Calculate Geweke’s statistic for current batch
      ssa <- var(x[1:A])
      ssb <- var(x[start:end])
      z[i] <- (mean(x[1:A]) - mean(x[start:end])) / sqrt(ssa / length(x[1:A]) + ssb / length(x[start:end]))
      
      # Calculate p-value for current batch
      p_value[i] <- 2 * pmin(pnorm(abs(z[i]), lower.tail = FALSE), pnorm(abs(z[i]), lower.tail = TRUE))
    }
  }
  else if (method == "coda"){
    for (i in 1:num_batches) {
      # Calculate start and end indices for current batch
      start <- B + (i - 1) * batch_size
      end <- start + batch_size
      
      # Calculate Geweke’s statistic for current batch
      ssa <- spectrum0.ar(x[1:A])$spec
      ssb <- spectrum0.ar(x[start:end])$spec
      z[i] <- (mean(x[1:A]) - mean(x[start:end])) / sqrt(ssa / length(x[1:A]) + ssb / length(x[start:end]))
      
      # Calculate p-value for current batch
      p_value[i] <- 2 * pmin(pnorm(abs(z[i]), lower.tail = FALSE), pnorm(abs(z[i]), lower.tail = TRUE))
    }
  }
  # Combine results from all batches
  result <- list(z = mean(z), p_value = mean(p_value))
  names(result) <- c("Geweke’s statistic", "p_value")
  return(result)
}


geweke(as.matrix(chain))
# $`Geweke’s statistic`
# [1] 0.4778486
# $p_value
# [1] 0.4073557

geweke(chain,method = "coda")
# $`Geweke’s statistic`
# [1] 0.3312266
# $p_value
# [1] 0.5306944


geweke(chain2)
# $`Geweke’s statistic`
# [1] -95.21099
# $p_value
# [1] 0.00768044

geweke(chain2,method = "coda")
# $`Geweke’s statistic`
# [1] -5.065428
# $p_value
# [1] 0.154139

# Wouldn't this violate the original intention of geweke? By default, geweke calculates the variance for the last 50% of the samples.

# The original intention of the Geweke method is to test whether the last 50% of the samples come from the same distribution, rather 
# than to test whether each batch comes from the same distribution. In this function, in order to achieve batch calculation, the sample 
# is divided into multiple small batches. Therefore, what this function implements is not the standard Geweke method, but a 
# modification of the Geweke method.


# So we should choose as many batches as possible when there are enough independent and identically distributed samples in each batch, 
# right?

# Generally speaking, we hope that there are enough independent and identically distributed samples in each batch, and we also want to 
# choose as many batches as possible to more accurately test whether the samples come from the same distribution.
# When using the Geweke method in practice, some basic analysis is usually performed on the data set to determine the nature and 
# distribution of the sample data, and then select the appropriate batch size and batch number based on this information.
# It should be noted that the selection of batch size and number of batches may have an impact on the test results, so when using the 
# Geweke method, it is necessary to adjust it according to the specific situation and conduct multiple experiments to verify the 
# reliability of the results.


library(coda)
geweke.diag(chain)
#var1 
#0.663 
geweke.diag(chain2)
#var1 
#0.01928
#-3.865 


# A version that can handle multi-dimensional situations XD.
geweke_multi <- function(x, a = 0.1, b = 0.5, num_batches = 20, method = "normal") {
  n <- dim(x)[1]
  m <- dim(x)[2]
  results <- lapply(1:m, function(i) {
    geweke(x[,i], a, b, num_batches, method)
  })
  names_vec <- c("Geweke’s statistic", "p_value")
  names(results) <- paste0("Chain", 1:m)
  for (i in 1:m) {
    names(results[[i]]) <- names_vec
  }
  return(results)
}

set.seed(123)
x <- MASS::mvrnorm(n = 1000, mu = c(0,0), Sigma = matrix(c(1,0.5,0.5,2), nrow = 2))
# Calculate Geweke's test for each dimension
geweke_multi(x)


# Geweke's conclusions may be inaccurate if the Markov chain takes a long time to warm up. Geweke's approach is based on the assumption 
# that the Markov chain is already stable. If the Markov chain has not yet converged to a stable distribution, Geweke's results may be 
# affected by the warm-up period, leading to erroneous conclusions. Therefore, samples from the warm-up period usually need to be 
# discarded before applying the Geweke method. 



# Heidelberger and Welch 
# The Heidelberger and Welch test (Heidelberger and Welch 1981, 1983) consists of two parts: a stationary portion test and a half-width
# test. The stationarity test assesses the stationarity of a Markov chain by testing the hypothesis that the chain comes from a covariance 
# stationary process. The half-width test checks whether the Markov chain sample size is adequate to estimate the mean values accurately.

Heidelberger <- function(x, eps = 0.1, pvalue = 0.05) {
  
  # Convert input to matrix 
  x <- as.matrix(x)
  
  # Initialize the Heidelberger and Welch object
  HW.mat0 <- matrix(0, ncol = 6, nrow = ncol(x))
  dimnames(HW.mat0) <- 
    list(colnames(x), c("Stationarity test", "Starting point", "pvalue", "Halfwidth test", "Mean", "Halfwidth"))
  
  HW.mat <- HW.mat0
  
  # Test each variable
  for (j in 1:ncol(x)) {
    # Generate a set of starting positions from the start position to halfway through the iterated sequence
    start.vec <- seq(from = 1, to = floor(nrow(x)/2), by = floor(nrow(x)/10))
    # Get the iteration sequence of the current variable
    Y <- x[, j, drop = TRUE]
    n1 <- length(Y)
    
    # Schruben's test, by iteratively reducing the length of the sequence, checks whether the sequence is convergent
    
    # Spectral estimates for the second half
    S0 <- spectrum0.ar(window(Y, start=nrow(x)/2))$spec
    converged <- FALSE
    for (i in seq_along(start.vec)) {
      Y <- window(Y, start = start.vec[i])
      n <- length(Y)
      ybar <- mean(Y)
      B <- cumsum(Y) - ybar * (1:n)
      Bsq <- (B * B)/(n * S0)
      I <- sum(Bsq)/n  # Actually this is Bn(s)_square
      
      if (converged <- !is.na(I) && pcramer(I) < 1 - pvalue)
        break
    }
    
    # If it converges, recalculate S0, calculate the half-width and whether it passes the half-width test
    S0ci <- spectrum0.ar(Y)$spec
    halfwidth <- 1.96 * sqrt(S0ci/n)
    passed.hw <- !is.na(halfwidth) & (abs(halfwidth/ybar) <= eps)
    
    # Returns NA if not converged or cannot be computed
    if (!converged || is.na(I) || is.na(halfwidth)) {
      nstart <- NA
      passed.hw <- NA
      halfwidth <- NA
      ybar <- NA
    } else {
      nstart <- start(Y)[1]
    }
    
    # Convert test results to text format
    converged <- ifelse(converged, "passed", "failed")
    passed.hw <- ifelse(passed.hw, "passed", "failed")
    # Store the test results of the current variable in HW.mat
    # Use the 'format' function to control the number of decimal places for pvalue, Mean, and Halfwidth
    HW.mat[j, ] <- c(converged, nstart, format(1 - pcramer(I), digits = 3), 
                     passed.hw, format(ybar, digits = 4), format(halfwidth, digits = 3))
  }
  
  # Return the matrix without the class attribute
  return(HW.mat)
}

# The pcramer function computes the probability of exceeding the given 
# value of a test statistic, based on the Cramer-Lundberg approximation.

pcramer <- function (q, eps = 1e-05) 
{
  # Set the logarithmic epsilon value
  log.eps <- log(eps)
  # Initialize a matrix of zeros to store the results
  y <- matrix(0, nrow = 4, ncol = length(q))
  for (k in 0:3) {
    # Compute the Cramer-Lundberg approximation for the given q
    z <- gamma(k + 0.5) * sqrt(4 * k + 1)/(gamma(k + 1) * pi^(3/2) * sqrt(q))
    # Compute the value of the test statistic
    u <- (4 * k + 1)^2/(16 * q)
    # If the value of the test statistic is too large (i.e., the 
    # exponent is too negative), set the result to zero to avoid underflow.
    # Otherwise, compute the value of the Bessel function of the second kind
    # with nu = 1/4 and multiply it with z and exp(-u).
    y[k + 1, ] <- ifelse(u > -log.eps, 0, z * exp(-u) * besselK(x = u, nu = 1/4))
  }
  # Return the sum of the rows of y
  return(apply(y, 2, sum))
}


Heidelberger(chain)
library(coda)
heidel.diag(chain)



# Gelman-Rubin
# Uses parallel chains with dispersed initial values to test whether they all converge to the same target distribution. Failure 
# could indicate the presence of a multi-mode posterior distribution (different chains converge to different local modes) or the need
# to run a longer chain (burn-in is yet to be completed).One-sided test based on a variance ratio test statistic. Large bRc values 
# indicate rejection.

# There are two main problems with the toy model. The first small problem is that the improved version of Brooks and Gelman (1997) was
# not considered when calculating the PSRF statistics (considering the possible existence of multicollinearity between samples, the 
# effective freedom is considered at the end of the calculation. degree). The second problem is that the function only considers chains
# with the same target distribution, that is, even if multiple Markov chains are input, they are considered to have the same target 
# distribution.

gelman_toy <- function(x,plot=TRUE){
  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.numeric(x)) stop("x must be numeric")
  if(dim(x)[1]<5) stop("each chain must have at least 5 iterations")
  if(dim(x)[2]<2) stop("there must be at least 2 chains")
  m <- dim(x)[2]
  n <- dim(x)[1]
  
  # between-chain variance B:
  
  theta_bar_m <- colMeans(x) 
  theta_bar <- (1/m)*sum(theta_bar_m)
  b <- rep(NA,dim(x)[2])
  for (i in 1:m){
    b[i] <- (theta_bar_m[i] - theta_bar) ^ 2
  }
  B <- (n/(m-1)) * sum(b)
  
  # within-chain variance W:
  s_m_square <- apply(x,2,var) 
  W <- (1/m) * sum(s_m_square)
  
  # The estimate of the variance V_hat:
  V_hat <- ((n-1)/n) * W + ((m+1)/n*m) * B
  
  # A refined version of PSRF R_hat:
  R_hat <- sqrt(varhat/W) # In the "real refined " version, R_hat will be recalculated. 
   
  if(plot){
    plot(R_hat,type="b",ylab="R hat",xlab="Parameter")
    abline(h=c(1.1,1.2),lty=c(2,3),col=c("red","orange"))
    legend("topright",lty=c(2,3),col=c("red","orange"),
           legend=c("R hat = 1.1","R hat = 1.2"))
  }
  cat("Potential scale reduction factors:\n")
  return(R_hat)
}

x <- cbind(chain,chain2)
gelman(x)

# Next, we consider ways to solve the above problems, In Brooks and Gelman's (1998) modified version of PSRF, df.V represents the 
# effective number of degrees of freedom for the model parameters, calculated as: (2 * V^2)/var.V. At the same time, when calculating
# the variance, we need to calculate the covariance between different variables (because the target distribution may be different at 
# this time), so we need to use cov.wb to calculate the covariance and use if to judge whether it is multivariate.
# Specific calculation method is in P24 of:
# file:///C:/Users/14817/Desktop/Convergence%20for%20MCMC/SAS_introbayes.pdf

# Please note that this function has only one variable by default, and its name is "x". If you need to put chain sample matrices of 
# different target distributions into this function as input, you need to input the number of variables and it is best to name them 
# respectively.

# Warning: When the autoburnin parameter is set to TRUE by default, the gelman.diag function automatically removes the first half 
# of the MCMC chain (the burn-in period). In practical experience, I have repeatedly verified that when autoburnin is set to 
# TRUE, the function tends to produce non-convergent results (the potential scale reduction factors increase, as do the point est.
# and Upper C.I.). Investigation is needed to determine whether setting autoburnin to TRUE by default is appropriate.

gelman <- function (x, nvar = 1, xnames = c("x"), confidence = 0.95, autoburnin = TRUE, multivariate = TRUE) 
{
  if (dim(x)[1] < 2) 
    stop("You need at least two chains")
  Niter <- dim(x)[1]
  Nchain <- dim(x)[2]
  Nvar <- nvar
  if (autoburnin && start(x)[1] < end(x)[1]/2) 
    x <- window(x, start = end(x)[1]/2 + 1)
  S2 <- array(sapply(as.data.frame(x), var, simplify = TRUE), dim = c(Nvar, Nvar, Nchain)) # calculating the variance in each chain.
  W <- apply(S2, c(1, 2), mean) # within-chain variance W 
  xbar <- matrix(apply(x, 2, mean), nrow = Nvar, ncol = Nchain) # calculating the mean in each chain.
  B <- Niter * var(t(xbar)) # between-chain variance B 
  if (Nvar > 1 && multivariate) { # Check if multivariate case
    if (is.R()) { 
      CW <- chol(W)
      emax <- eigen(backsolve(CW, t(backsolve(CW, B, transpose = TRUE)), 
                              transpose = TRUE), symmetric = TRUE, only.values = TRUE)$values[1] # Finding the largest eigenvalue
    }
    else { # If using non-R implementation
      emax <- eigen(qr.solve(W, B), symmetric = FALSE, 
                    only.values = TRUE)$values
    }
    mpsrf <- sqrt((1 - 1/Niter) + (1 + 1/Nvar) * emax/Niter) # Calculate multivariate potential scale reduction factor (MPSRF)
  }
  else mpsrf <- NULL # If not multivariate, set MPSRF to NULL
  
  # Get the diagonal of the within-chain variance
  w <- diag(W)
  # Get the diagonal of the between-chain variance
  b <- diag(B)
  # Get the variance of each chain
  s2 <- matrix(apply(S2, 3, diag), nrow = Nvar, ncol = Nchain) 
  # Compute the estimate of the marginal posterior mean for each variable
  muhat <- apply(xbar, 1, mean) 
  # Calculate the within-chain variance
  var.w <- apply(s2, 1, var)/Nchain 
  # Calculate the between-chain variance
  var.b <- (2 * b^2)/(Nchain - 1) 
  # Compute the covariance between the within-chain and between-chain variance estimates
  cov.wb <- (Niter/Nchain) * diag(var(t(s2), t(xbar^2)) - 2 * 
                                    muhat * var(t(s2), t(xbar))) 
  # Calculate the estimate of the marginal posterior variance for each variable
  V <- (Niter - 1) * w/Niter + (1 + 1/Nchain) * b/Niter 
  # Calculate the variance of the estimate of the marginal posterior variance for each variable
  var.V <- ((Niter - 1)^2 * var.w + (1 + 1/Nchain)^2 * var.b + 
              2 * (Niter - 1) * (1 + 1/Nchain) * cov.wb)/Niter^2 
  # Calculate the effective number of independent samples
  df.V <- (2 * V^2)/var.V 
  # Calculate the adjusted degrees of freedom, as suggested by Brooks and Gelman (1997)
  df.adj <- (df.V + 3)/(df.V + 1) # R_hat_adj = sqrt((d+3)*V_hat/(d-1)*W)
  # Between-chain degrees of freedom
  B.df <- Nchain - 1 
  # Within-chain degrees of freedom
  W.df <- (2 * w^2)/var.w 
  # Estimate of the fixed variance component
  R2.fixed <- (Niter - 1)/Niter 
  # Estimate of the random variance component
  R2.random <- (1 + 1/Nchain) * (1/Niter) * (b/w) 
  # Total estimate of the variance component
  R2.estimate <- R2.fixed + R2.random 
  # Upper bound of the variance component
  R2.upper <- R2.fixed + qf((1 + confidence)/2, B.df, W.df) * R2.random
  # Compute the potential scale reduction factor (PSRF)
  psrf <- cbind(sqrt(df.adj * R2.estimate), sqrt(df.adj * R2.upper))
  dimnames(psrf) <- list(xnames, c("Point est.", "Upper C.I."))
  out <- list(psrf = psrf, mpsrf = mpsrf)
  return(out)
}

  


x <- cbind(chain,chain2)
gelman(x)

library(coda)
x <- mcmc.list(as.mcmc(chain),as.mcmc(chain2))
gelman.diag(x)
gelman.plot(x)

# A large PSRF indicates that the between-chain variance is substantially greater than the within-chain variance, so that longer simulation
# is needed. If the PSRF is close to 1, you can conclude that each of the M chains has stabilized, and they are likely to have reached the 
# target distribution.
# It is best to choose different initial values for all M chains. The initial values should be as dispersed from each other as possible so
# that the Markov chains can fully explore different parts of the distribution before they converge to the target. Similar initial values 
# can be risky because all of the chains can get stuck in a local maximum; that is something this convergence test cannot detect. If you do
# not supply initial values forall the different chains, the procedures generate them for you.
