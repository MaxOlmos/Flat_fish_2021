x <- rlnorm(10,log(10),10)
hist(x)

order <- order(x,decreasing=T)
cumsum <- cumsum(x[order])/sum(x)
threshold <- max(x[order[cumsum>0.95]])
plot(x[order],cumsum,type="l")
layer_t = ifelse( x > threshold[,u], 1, 0 )