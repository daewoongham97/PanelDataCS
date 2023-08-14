df = read.csv("~/Downloads/omakase_500000.csv")

###
t = 60
# decision to focus on accounts that at least have 60 days of allocation.
in_accs =  unique(df$account_id[df$days_since_allocation == 60])

trt_accs = unique(df$account_id[df$test_cell_nbr == 3])
ctrl_accs = unique(df$account_id[df$test_cell_nbr == 1])

trt_accs = trt_accs[trt_accs %in% in_accs]
ctrl_accs = ctrl_accs[ctrl_accs %in% in_accs]


trt_accs = sort(trt_accs)
ctrl_accs = sort(ctrl_accs)

y_trts = matrix(NA, nrow = length(trt_accs), ncol = t)
y_ctrls = matrix(NA, nrow = length(ctrl_accs), ncol = t)

for (i in 1:t) {
  
  in_df = df[df$days_since_allocation == as.character(i), ]
  
  in_df = in_df[order(in_df$account_id), ]
  
  ## added treatment
  new_y_trt = rep(0, length(trt_accs))
  
  trt_accs_considered = in_df$account_id[in_df$account_id %in% trt_accs]
  
  new_y_trt[trt_accs %in% trt_accs_considered] = in_df$qualified_play[in_df$account_id %in% trt_accs]
  
  y_trts[, i] = new_y_trt
  
  ## added control
  new_y_ctrl = rep(0, length(ctrl_accs))
  
  ctrl_accs_considered = in_df$account_id[in_df$account_id %in% ctrl_accs]
  
  new_y_ctrl[ctrl_accs %in% ctrl_accs_considered] = in_df$qualified_play[in_df$account_id %in% ctrl_accs]
  
  y_ctrls[, i] = new_y_ctrl
}

S1 = cov(y_trts)
S0 = cov(y_ctrls)


hatbeta = sigma = added_trts = added_ctrls = vector()
n1 = length(trt_accs)
n0 = length(ctrl_accs)
i = 1
in_df = df[df$days_since_allocation == as.character(i), ]
trt = in_df$qualified_play[in_df$test_cell_nbr == 3]
ctrl = in_df$qualified_play[in_df$test_cell_nbr == 1]
added_trts[i] = sum(!(trt_accs %in% in_df$account_id[in_df$test_cell_nbr == 3]))

added_ctrls[i] = sum(!(ctrl_accs %in% in_df$account_id[in_df$test_cell_nbr == 1]))
new_trt = c(trt, rep(0, sum(!(trt_accs %in% in_df$account_id[in_df$test_cell_nbr == 3]))))

new_ctrl = c(ctrl, rep(0, sum(!(ctrl_accs %in% in_df$account_id[in_df$test_cell_nbr == 1]))))

hatbeta[i] = mean(new_trt) - mean(new_ctrl)
sigma[i] = S1[1, 1]/n1 + S0[1, 1]/n0
unaccounted_var = sigma

for (i in 2:t) {
  in_df = df[df$days_since_allocation == as.character(i), ]
  trt = in_df$qualified_play[in_df$test_cell_nbr == 3]
  ctrl = in_df$qualified_play[in_df$test_cell_nbr == 1]
  added_trts[i] = sum(!(trt_accs %in% in_df$account_id[in_df$test_cell_nbr == 3]))
  
  added_ctrls[i] = sum(!(ctrl_accs %in% in_df$account_id[in_df$test_cell_nbr == 1]))
  new_trt = c(trt, rep(0, sum(!(trt_accs %in% in_df$account_id[in_df$test_cell_nbr == 3]))))
  
  new_ctrl = c(ctrl, rep(0, sum(!(ctrl_accs %in% in_df$account_id[in_df$test_cell_nbr == 1]))))
  
  hatbeta[i] = mean(new_trt) - mean(new_ctrl)
  sigma_t1 = as.numeric((S1[i, i] - matrix(S1[i, (1:(i-1))], ncol = i - 1)%*% solve(S1[1:(i-1), 1:(i-1)]) %*% matrix(S1[(1:(i-1)), i], nrow = (i-1))))/n1
  sigma_t0 = as.numeric((S0[i, i] - matrix(S0[i, (1:(i-1))], ncol = i - 1)%*% solve(S0[1:(i-1), 1:(i-1)]) %*% matrix(S0[(1:(i-1)), i], nrow = (i-1))))/n0
  
  sigma[i] = (sigma_t1 + sigma_t0)
  unaccounted_var[i] = S1[i, i]/n1 + S0[i,i]/n0
}

library(ggplot2)
eta = 0.77; alpha = 0.05
S = cumsum(sigma)
center = cumsum(hatbeta)/(1:t)
radius = sqrt((S*eta^2 + 1)/eta^2 * log((S*eta^2 + 1)/alpha^2))*(1/(1:t))
lower = center - radius
upper = center + radius

plotting_df = data.frame(time = 1:length(lower), lower = lower,
                         upper = upper)

S = cumsum(unaccounted_var)
radius = sqrt((S*eta^2 + 1)/eta^2 * log((S*eta^2 + 1)/alpha^2))*(1/(1:t))
lower = center - radius
upper = center + radius


plotting_df2 = data.frame(time = 1:length(lower), lower = lower,
                          upper = upper)

final_plot_df = rbind(plotting_df, plotting_df2)


final_plot_df$group = factor(rep(c("Accounted Var", "Unaccounted Var"), each = t))


# plotting configs
s = 5
w = 50
s2 = 5
a1 = 15
a2 = 20

Q_CS = ggplot(data = final_plot_df, aes(x = time)) + geom_ribbon(aes(ymin = lower, ymax = upper, col = group), fill = "white", size = 1, alpha = 0.15, position = position_jitter(width = 0.0, height = 0.05)) +
  scale_color_manual(name = "", values = c("Accounted Var" = "blue", "Unaccounted Var" = "red")) + 
  geom_hline(yintercept = 0, col = "black", linetype = "dashed", size = 1) + 
  xlab("Day") + ggtitle("Recommendation Algorithm Experiment: n = 1500") + ylab("ATE: Qualified Play Days") +
  theme(axis.text=element_text(size=a1), axis.title=element_text(size=a2,face="bold"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), 
        legend.position= "top", legend.text=element_text(size=a1),legend.key.size = unit(1,"line"),plot.title = element_text(size = a2, face = "bold"), axis.title.x = element_text(vjust=-0.5))


pdf(file = "~/Downloads/omakase_analysis_fig.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 5)

Q_CS

dev.off()


which(plotting_df$upper <= 0)[1]


