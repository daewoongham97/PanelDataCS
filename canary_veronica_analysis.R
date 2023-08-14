df <- read.csv(file="~/Downloads/veronica_df.csv")


day1 = df[df$utc_date == "20230719", ]

day2 = df[df$utc_date == "20230720", ]

get_CS = function(metric, day1, day2, eta = 0.77, alpha = 0.05) {
  ## ith hours jth minute
  betahat = s = vector()
  for (i in 0:23) {
    ith_hour = day1[day1$utc_hour == i, ]
    for (j in 1:60) {
      ith_jth_minute = ith_hour[ith_hour$minute_of_hour == (j -1), ]
      y = as.numeric(ith_jth_minute[, metric])
      trt = ith_jth_minute$test_cell_nbr
      trt = trt - 1
      lm_mod = lm(y ~ trt)
      hat_beta = coef(lm_mod)[2]
      var = as.numeric(diag(vcovHC(lm_mod, type = "HC")))[2]
      betahat = c(betahat, hat_beta)
      s = c(s, var)
    }
    print(paste0("Day 1: Hour ", i, "Done"))
  }
  
  ## ith hours jth minute
  for (i in 0:23) {
    ith_hour = day2[day2$utc_hour == i, ]
    for (j in 1:60) {
      ith_jth_minute = ith_hour[ith_hour$minute_of_hour == (j -1), ]
      y = as.numeric(ith_jth_minute[, metric])
      trt = ith_jth_minute$test_cell_nbr
      trt = trt - 1
      lm_mod = lm(y ~ trt)
      hat_beta = coef(lm_mod)[2]
      var = as.numeric(diag(vcovHC(lm_mod, type = "HC")))[2]
      betahat = c(betahat, hat_beta)
      s = c(s, var)
    }
    print(paste0("Day 2: Hour ", i, "Done"))
  }
  
  eta = 0.77; alpha = 0.05
  t = length(s)
  S = cumsum(s)
  center = cumsum(betahat)/(1:t)
  radius = sqrt((S*eta^2 + 1)/eta^2 * log((S*eta^2 + 1)/alpha^2))*(1/(1:t))
  lower = center - radius
  upper = center + radius
  result_df = data.frame(center, lower, upper)
  return(result_df)
}

get_plot = function(result_df, ylabel, ylimits, title = "Canary Sequential Tests") {
  lower = result_df$lower
  upper = result_df$upper
  plotting_df = data.frame(time = 1:length(lower), lower = lower,
                           upper = upper)
  
  s = 5
  w = 50
  s2 = 5
  a1 = 15
  a2 = 20
  
  Q_CS = ggplot(data = plotting_df, aes(x = time)) + geom_ribbon(aes(ymin = lower, ymax = upper, fill = "blue"), alpha = 0.5) +
    geom_hline(yintercept = 0, col = "black", linetype = "dashed", size = 1) + ylim(c(ylimits[1],ylimits[2])) + 
    xlab("Minute") + ggtitle(title) + ylab(ylabel) + scale_x_continuous(breaks=c( 0, 1000, 2000)) + 
    theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position= "none", legend.title=element_blank(),legend.key=element_blank(), legend.text=element_text(size=a1),legend.key.size = unit(4,"line"),plot.title = element_text(size = a2, face = "bold"), axis.title.x = element_text(vjust=-0.5))
  
  return(Q_CS)
}


# network bytes:
# network interface controller controlling the physical number of bits coming out of the server
# Server and device and sending content the stream of bits that stream of bits is exactly that network bytes 
# not just content and video thumbnail
# flow of data out of OCA 
# Aggregate all the bytes netfli from all the netflix members getting from the server 

# network bytes
df$network_GB = as.numeric(df$network_bytes_per_sec_out)/10^9

network_bytes = get_CS(metric = "network_GB", day1 = day1, day2 = day2)

which(network_bytes$upper <0)[1] # 128 th minute

network_bytes_plot = get_plot(network_bytes, "Network GB", ylimits = c(-0.5, 0.5))

network_bytes_plot


# network bytes
network_bytes2 = get_CS(metric = "network_bytes_per_sec_out", day1 = day1, day2 = day2)

which(network_bytes2$upper <0)[1] # 1595th minute

network_bytes2_plot = get_plot(network_bytes2, "Network Bytes", ylimits = c(-10^9, 10^9))

network_bytes2_plot


# busy_cpu
# percent how busy is your cpu from a 0-100 scale
# Cost effficient to see how busy your cpu 
# don't want to see that be too large 

cpu_df = get_CS(metric = "cpu_busy", day1 = day1, day2 = day2)

which(cpu_df$upper <0)[1] # 34 th minute

cpu_busy_plot = get_plot(cpu_df, "Busy CPU", ylimits = c(-0.5, 0.5), title = "")

cpu_busy_plot

# ngix
ngix_df = get_CS(metric = "nginx_http_bytes_total_occ", day1 = day1, day2 = day2)

which(ngix_df$upper <0)[1] # never

ngix_plot = get_plot(ngix_df, "Ngix Total Occ", ylimits = c(-5000000, 5000000))

ngix_plot

library(ggpubr)


pdf(file = "~/Downloads/canary_test.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 4.5)

ggarrange(network_bytes_plot, cpu_busy_plot, nrow = 1)


dev.off()








