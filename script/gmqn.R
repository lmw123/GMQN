#authour: Mengwei Li
#email: limengwei@big.ac.cn


gmqn <- function(m, um, probe,
                 t1.green.ref.mean=c(361.0165,9646.7525),
                 t1.green.ref.sd=c(278.9551,5109.423),
                 t1.red.ref.mean=c(696.0015,12036.3621),
                 t1.red.ref.sd=c(519.3141,7002.863)){
  # Gaussian Mixture Quantile Normalization (GMQN) implement
  #
  # Args:
  #   m: A dataframe stores methylation signals with probes in the rows, samples in the columns
  #   um: A dataframe stores unmethylation signals with probes in the rows, samples in the columns
  #   probe: the probe id for m and um
  #   t1.green.ref.mean: the reference means of two-state Gaussian mixture distribution for type I green probe
  #   t1.green.ref.sd: the reference standard errors of two-state Gaussian mixture distribution for type I green probe
  #   t1.red.ref.mean: the reference means of two-state Gaussian mixture distribution for type I red probe
  #   t1.red.ref.sd: the reference standard errors of two-state Gaussian mixture distribution for type I red probe
  #
  # Returns:
  #   Normalized beta values

  require(mclust)
  
  
  print("----------Fitting Gaussian mixture model for probe of type 1 red----------")
  
  t1.red.index = match(t1.red, probe)
  t1.red.index = t1.red.index[which(!is.na(t1.red.index))]
  
  print(length(t1.red.index))
  
  t1.red.signal = c(m[t1.red.index], um[t1.red.index])
  
  t1.red.model = Mclust(t1.red.signal, G=2)
  t1.red.mean = t1.red.model$parameters$mean
  t1.red.sd = sqrt(t1.red.model$parameters$variance$sigmasq)

  
  print("----------Normalizing probe of type 1 red----------")
  
  t1.red.signal[which(t1.red.model$classification == 1)] = qnorm(pnorm(t1.red.signal[which(t1.red.model$classification == 1)],
                                                            t1.red.mean[1], t1.red.sd[1]), t1.red.ref.mean[1], t1.red.ref.sd[1])
  print(length(which(t1.red.model$classification == 1)))
  
  t1.red.signal[which(t1.red.model$classification == 2)] = qnorm(pnorm(t1.red.signal[which(t1.red.model$classification == 2)],
                                                            t1.red.mean[2], t1.red.sd[2]), t1.red.ref.mean[2], t1.red.ref.sd[2])
  m[t1.red.index] = t1.red.signal[1:(length(t1.red.signal)/2)]
  um[t1.red.index] = t1.red.signal[(length(t1.red.signal)/2+1):length(t1.red.signal)]


  print("----------Fitting Gaussian mixture model for probe of type 1 green----------")
  
  t1.green.index= match(t1.green,probe)
  t1.green.index = t1.green.index[which(!is.na(t1.green.index))]
  
  t1.green.signal = c(m[t1.green.index],um[t1.green.index])
  
  t1.green.model = Mclust(t1.green.signal, G=2)
  t1.green.mean = t1.green.model$parameters$mean
  t1.green.sd = sqrt(t1.green.model$parameters$variance$sigmasq)
  
  
  print("----------Normalizing probe of type 1 green----------")
  
  t1.green.signal[which(t1.green.model$classification == 1)] = qnorm(pnorm(t1.green.signal[which(t1.green.model$classification == 1)],
                                                            t1.green.mean[1], t1.green.sd[1]), t1.green.ref.mean[1], t1.green.ref.sd[1])
  t1.green.signal[which(t1.green.model$classification == 2)] = qnorm(pnorm(t1.green.signal[which(t1.green.model$classification == 2)],
                                                            t1.green.mean[2], t1.green.sd[2]), t1.green.ref.mean[2], t1.green.ref.sd[2])

  m[t1.green.index] = t1.green.signal[1:(length(t1.green.signal)/2)]
  um[t1.green.index] = t1.green.signal[(length(t1.green.signal)/2+1):length(t1.green.signal)]

  
  print("----------Detecting p value----------")
  
  t2.index = match(t2, probe)
  t2.index = t2.index[which(!is.na(t2.index))]
  
  pIR <- apply(cbind(1 - pnorm(m[t1.red.index], t1.red.mean[1], t1.red.sd[1]), 
                     1 - pnorm(um[t1.red.index], t1.red.mean[1], t1.red.sd[1])),
               1, min)
  pIG <- apply(cbind(1 - pnorm(m[t1.green.index], t1.green.mean[1], t1.green.sd[1]), 
                     1 - pnorm(um[t1.green.index], t1.green.mean[1], t1.green.sd[1])),
               1, min)
  pII <- apply(cbind(1 - pnorm(um[t2.index], t1.green.mean[1], t1.green.sd[1]), 
                     1 - pnorm(m[t2.index], t1.red.mean[1], t1.red.sd[1])),
               1, min)

  p=cbind(c(pIR, pIG , pII),
          c(probe[t1.red.index], probe[t1.green.index], probe[t2.index]))
  row.names(p) = c(probe[t1.red.index], probe[t1.green.index], probe[t2.index])

  m[which(m <= 0)] = min(m[which(m > 0)])
  um[which(um <= 0)] = min(um[which(um > 0)])

  normalized.signal = data.frame(cbind(round(as.numeric(m), 0),
                         round(as.numeric(um), 0),
                         as.numeric(p[probe, 1])))

  row.names(normalized.signal) = probe
  names(normalized.signal) = c("m", "um", "p")
  
  
  print("----------Probe type normalization----------")
  
  
  normalized.signal.pass = normalized.signal[which(normalized.signal$p <= 0 & !is.na(normalized.signal$p)),]
  
  normalized.signal.pass[which(normalized.signal.pass$m == Inf), "m"] = max(normalized.signal.pass$m[normalized.signal.pass$m != Inf])
  normalized.signal.pass[which(normalized.signal.pass$um == Inf), "um"] = max(normalized.signal.pass$um[normalized.signal.pass$um != Inf])
  
  beta = normalized.signal.pass$m / (normalized.signal.pass$m + normalized.signal.pass$um)
  type = probe.type[row.names(normalized.signal.pass),2]
  
  type[which(type=='I')] = 1
  type[which(type=='II')] = 2
  
  type = as.numeric(type)
  
  beta.adjust = BMIQ(beta, type, plots = F)
  beta.adjust = beta.adjust$nbeta
  
  normalized.signal[row.names(normalized.signal.pass),'beta'] = beta.adjust
  
  return(normalized.signal)
}


