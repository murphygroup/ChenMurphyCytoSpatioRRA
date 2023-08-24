library(spatstat)
library(spatstat.utils)
library(spatstat.data)
library(ggplot2)
library(dplyr)
library(permute)
library(data.table)
# script_dir = '/home/hrchen/Documents/Research/CancerECT/scripts'
script_dir = '/home/haoranch/projects/HuBMAP/ppm/script'
setwd(script_dir)
source('mppm.ppp.R')
source('mppm.quad.ppp.R')
source('bt.frame.ppp.R')
source('quadscheme.ppp.R')
source('default.dummy.ppp.R')
source('default.n.tiling.ppp.R')
source('mpl.engine.ppp.R')
source('mpl.prepare.ppp.R')
source('evalInteraction.ppp.R')
source('quadBlockSizes.ppp.R')
source('evalInterEngine.ppp.R')
source('evaluate.ppp.R')
source('evalPairPotential.ppp.R')

# 
# get_shuffle_all = function(data){
#   pattern_ID_list = unique(data$pattern_ID)
#   pattern_ID_list = sample(pattern_ID_list, length(pattern_ID_list))
#   pattern_ID_list_num = length(pattern_ID_list)
#   first_half = c()
#   second_half = c()
#   
#   for (i in 1:pattern_ID_list_num){
#     pattern_current = data[which(data$pattern_ID == pattern_ID_list[i]),]
#     # print(c('pattern_current', length(pattern_current$x)))
#     # print(pattern_ID_list[i])
#     if (i < floor(pattern_ID_list_num / 2)){
#       if (length(first_half) == 0){
#         first_half = pattern_current
#       } else{
#         first_half = rbindlist(list(first_half, pattern_current))
#       }
#     } else{
#       if (length(second_half) == 0){
#         second_half = pattern_current
#       } else{
#         second_half = rbindlist(list(second_half, pattern_current))
#       } 
#     }
#     
#   }
#   first_half = as.data.frame(first_half)
#   second_half = as.data.frame(second_half)
#   first_half = first_half[, 1-ncol(first_half)]
#   second_half = second_half[, 1-ncol(second_half)]
#   return(list(test1 = first_half, test2 = second_half))
# }
# 
# get_shuffle_training_validation = function(training_validation){
#   
#   pattern_ID_list = unique(training_validation$pattern_ID)
#   pattern_ID_list = sample(pattern_ID_list, length(pattern_ID_list))
#   pattern_ID_list_num = length(pattern_ID_list)
#   training = c()
#   validation = c()
#   
#   for (i in 1:pattern_ID_list_num){
#     pattern_current = training_validation[which(training_validation$pattern_ID == pattern_ID_list[i]),]
#     # print(c('pattern_current', length(pattern_current$x)))
#     # print(pattern_ID_list[i])
#     if (i < floor(pattern_ID_list_num / 2)){
#       if (length(training) == 0){
#         training = pattern_current
#       } else{
#         training = rbindlist(list(training, pattern_current))
#       }
#     } else{
#       if (length(validation) == 0){
#         validation = pattern_current
#       } else{
#         validation = rbindlist(list(validation, pattern_current))
#       } 
#     }
#     
#   }
#   training = as.data.frame(training)
#   validation = as.data.frame(validation)
#   training = training[, 1-ncol(training)]
#   validation = validation[, 1-ncol(validation)]
#   Quad_training = Quad_all
#   Quad_training$moadf = training
#   # training
#   # model_training = mppm.fit.pp3(Y ~ marks, data=H, moadf=pattern_training, interaction=MultiStrauss(radii=r_matrix))
#   
#   model_train = mppm.fit.pp3(Data=Quad_training, interaction=MultiStrauss(radii=r_matrix))
#   return(list(train = training, val = validation, model = model_train))
# }
# 
# 
# get_shuffle_training_validation_pointwise = function(training_validation){
#   pattern_ID_list = unique(training_validation$pattern_ID)
#   pattern_ID_list = sample(pattern_ID_list, length(pattern_ID_list))
#   print(pattern_ID_list)
#   pattern_ID_list_num = length(pattern_ID_list)
#   
#   for (i in 1:pattern_ID_list_num){
#     pattern_current = training_validation[which(training_validation$pattern_ID == pattern_ID_list[i]),]
#     # print(c('pattern_current', length(pattern_current$x)))
#     # print(pattern_ID_list[i])
#     if (i == 1){
#       training_validation_shuf = pattern_current
#     } else{
#       training_validation_shuf = rbind(training_validation_shuf, pattern_current)
#     }
#     
#   }
#   
#   total_num = dim(training_validation_shuf)[1]
#   training = training_validation_shuf[1:floor(total_num/2), 1-ncol(training_validation_shuf)]
#   validation = training_validation_shuf[floor(total_num/2):total_num, 1-ncol(training_validation_shuf)]
#   Quad_training = Quad_all
#   Quad_training$moadf = training
#   # training
#   # model_training = mppm.fit.pp3(Y ~ marks, data=H, moadf=pattern_training, interaction=MultiStrauss(radii=r_matrix))
#   
#   model_train = mppm.fit.pp3(Data=Quad_training, interaction=MultiStrauss(radii=r_matrix))
#   return(list(train = training, val = validation, model = model_train))
# }
# 
# # moadf = Quad_training$moadf
# # fmla = Quad_training$fmla
# # glmmsubset <- .mpl.SUBSET <- moadf$.mpl.SUBSET
# # .mpl.W <- moadf$.mpl.W
# # weights = .mpl.W * caseweight
# # caseweight <- moadf$caseweight
# # fitter <- "glm"
# # ctrl <- do.call(glm.control, resolve.defaults(gcontrol, list(maxit = 50)))
# # FIT <- glm(fmla, family = quasi(link = "log", variance = "mu"), 
# #            weights = .mpl.W * caseweight, data = fit_data, subset = (.mpl.SUBSET == "TRUE"), control = ctrl)
# 
# get_maxlogpl = function(){
#   
#   # if (test_set == 'training'){
#   #   maxlogpl_training = model_training$maxlogpl
#   #   devi = model_training$Fit$FIT$deviance
#   #   SUBSET <- model_training$Fit$FIT$data$.mpl.SUBSET
#   #   Z <- (model_training$Fit$FIT$data$.mpl.Y != 0)
#   #   # maxlogpl_training = model_training$maxlogpl
#   #   return(c(maxlogpl_training, sum(Z & SUBSET), devi))
#   #   # return(c(maxlogpl_training, length(model_training$Fit$FIT$data$.mpl.Y), devi))
#   # } else{
#     fmla = model_training$Fit$fmla
#     family = model_training$Fit$FIT$family
#     # fit_data = pattern
#     .mpl.W = fit_data$.mpl.W
#     caseweight = fit_data$caseweight
#     # weights = pattern_all$.mpl.W
#     # data = pattern_all[,1:ncol(data)-1]
#     gcontrol = list()
#     coef <- model_training$Fit$FIT$coefficients
#     # FIT <- glm(fmla, family = family, weights = .mpl.W, 
#     #            data = test1, subset = .mpl.SUBSET, control = gcontrol, 
#     #            model = FALSE)
#     ctrl <- do.call(glm.control, resolve.defaults(gcontrol, 
#                                                   list(maxit = 50)))
#     devi_list <- glm.prep(fmla, family = quasi(link = "log", variance = "mu"), weights = .mpl.W * caseweight, data = fit_data, start = coef, subset = (fit_data$.mpl.SUBSET == "TRUE"), control = ctrl)
#     
#     # devi = glm.prep(fmla, family, weights = weights,
#     #                 data = data, subset = data$.mpl.SUBSET, start = coef, control = gcontrol,
#     #                 model = FALSE)
#     W <- with(fit_data, .mpl.W * caseweight)
#     SUBSET <- fit_data$.mpl.SUBSET
#     Z <- (fit_data$.mpl.Y != 0)
#     devi = sum(devi_list[Z & SUBSET])
#     # devi = sum(devi_list)
#     maxlogpl_all = -(devi/2 + sum(log(W[Z & SUBSET])) + sum(Z & SUBSET))
#     # maxlogpl_all = -(devi/2 + sum(log(W[Z & SUBSET])) + sum(Z & SUBSET))
#     return(c(maxlogpl_all, sum(Z & SUBSET), devi))
#     # return(c(maxlogpl_all, length(fit_data$.mpl.Y), devi))
#   # }
# }
# 
# get_dataframe = function(data, data_type){
#   data_list = list(shuffle_num = 1:length(data), maxlogpl=data, type=rep(data_type, length(data)))
#   data_dataframe = as.data.frame(do.call(cbind, data_list))
#   data_dataframe$maxlogpl = as.numeric(as.character(data_dataframe$maxlogpl))
#   data_dataframe$shuffle_num = as.numeric(as.character(data_dataframe$shuffle_num))
#   return(data_dataframe)
# }
# 
# get_relative_error = function(training, test, data_type){
#   error = abs((training - test) / training) * 100
#   data_list = list(shuffle_num = 1:length(error), relative_error=error, type=rep(data_type, length(error)))
#   data_dataframe = as.data.frame(do.call(cbind, data_list))
#   data_dataframe$relative_error = as.numeric(as.character(data_dataframe$relative_error))
#   data_dataframe$shuffle_num = as.numeric(as.character(data_dataframe$shuffle_num))
#   return(data_dataframe)
# }
# 
# get_figures = function(){  
#   # training_maxlogpl_dataframe = get_dataframe(training_maxlogpl, 'training')
#   # validation_maxlogpl_dataframe = get_dataframe(validation_maxlogpl, 'validation')
#   # test_maxlogpl_dataframe = get_dataframe(test_maxlogpl, 'test')
#   # final_dataframe = rbind(training_maxlogpl_dataframe, validation_maxlogpl_dataframe, test_maxlogpl_dataframe)
#   # p = ggplot(data=final_dataframe, aes(x=shuffle_num, y=maxlogpl, group=type, color=type)) + geom_line()
#   # p = p + xlab("Shuffle Number")
#   # p = p + ylab("Pseudo-likelihood")
#   # # p = p + ylim(c(2.4, 2.65))
#   # ggsave('/data2/CancerECT/fig/normal_training_radius_5_cluster_0_2_original_maxlogpl.png', p, width = 7, height = 5)  
#   # 
#   # training_maxlogpl_dataframe = get_dataframe(training_maxlogpl / training_num, 'training')
#   # validation_maxlogpl_dataframe = get_dataframe(validation_maxlogpl / validation_num, 'validation')
#   # test_maxlogpl_dataframe = get_dataframe(test_maxlogpl / test_num, 'test')
#   # final_dataframe = rbind(training_maxlogpl_dataframe, validation_maxlogpl_dataframe, test_maxlogpl_dataframe)
#   # p = ggplot(data=final_dataframe, aes(x=shuffle_num, y=maxlogpl, group=type, color=type)) + geom_line()
#   # p = p + xlab("Shuffle Number")
#   # p = p + ylab("Normalized Pseudo-likelihood")
#   # # p = p + ylim(c(2.4, 2.65))
#   # p
#   # ggsave('/data2/CancerECT/fig/normal_training_radius_5_cluster_0_2_normalized_maxlogpl.png', p, width = 7, height = 5)  
#   # 
#   # training_maxlogpl_dataframe = get_dataframe(training_num, 'training')
#   # validation_maxlogpl_dataframe = get_dataframe(validation_num, 'validation')
#   # test_maxlogpl_dataframe = get_dataframe(test_num, 'test')
#   # final_dataframe = rbind(training_maxlogpl_dataframe, validation_maxlogpl_dataframe, test_maxlogpl_dataframe)
#   # 
#   # p = ggplot(data=final_dataframe, aes(x=shuffle_num, y=maxlogpl, group=type, color=type)) + geom_line()
#   # p = p + xlab("Shuffle Number")
#   # p = p + ylab("Number of Particles")
#   # # p = p + ylim(c(2.4, 2.65))
#   # p 
#   # ggsave('/data2/CancerECT/fig/normal_training_radius_5_cluster_0_2_num.png', p, width = 7, height = 5)  
#   # 
#   # training_validation_error = get_relative_error(training_maxlogpl / training_num, validation_maxlogpl / validation_num, 'training vs validation')
#   # training_all_error = get_relative_error(training_maxlogpl / training_num, test_maxlogpl / test_num, 'training vs test')
#   # final_dataframe = rbind(training_validation_error, training_all_error)
#   # 
#   # p = ggplot(data=final_dataframe, aes(x=shuffle_num, y=relative_error, group=type, color=type)) + geom_line()
#   # p = p + xlab("Shuffle Number")
#   # p = p + ylab("Percent Difference")
#   # # p = p + ylim(c(2.4, 2.65))
#   # p 
#   # ggsave('/data2/CancerECT/fig/normal_training_radius_5_cluster_0_2_normalized_maxlogpl_error.png', p, width = 7, height = 5)  
#   # 
#   training_devi_dataframe = get_dataframe(training_devi / training_num, 'training')
#   validation_devi_dataframe = get_dataframe(validation_devi / validation_num, 'validation')
#   test1_devi_dataframe = get_dataframe(test1_devi / test1_num, 'test1')
#   test2_devi_dataframe = get_dataframe(test2_devi / test2_num, 'test2')
#   final_dataframe = rbind(training_devi_dataframe, validation_devi_dataframe, test1_devi_dataframe, test2_devi_dataframe)
#   
#   p = ggplot(data=final_dataframe, aes(x=shuffle_num, y=maxlogpl, group=type, color=type)) + geom_line()
#   p = p + xlab("Shuffle Number")
#   p = p + ylab("Normalized Deviance")
#   # p = p + ylim(c(2.8, 5.4))
#   p 
#   ggsave('/data3/CancerECT/all_batch/fig/normal_training_radius_5_cluster_0_2_normalized_deviance.png', p, width = 7, height = 5)  
#   
#   training_devi_dataframe = get_dataframe(training_devi, 'training')
#   validation_devi_dataframe = get_dataframe(validation_devi, 'validation')
#   test1_devi_dataframe = get_dataframe(test1_devi, 'test1')
#   test2_devi_dataframe = get_dataframe(test2_devi, 'test2')
#   final_dataframe = rbind(training_devi_dataframe, validation_devi_dataframe, test1_devi_dataframe, test2_devi_dataframe)
#   
#   p = ggplot(data=final_dataframe, aes(x=shuffle_num, y=maxlogpl, group=type, color=type)) + geom_line()
#   p = p + xlab("Shuffle Number")
#   p = p + ylab("Deviance")
#   # p = p + ylim(c(2.4, 2.65))
#   p 
#   ggsave('/data3/CancerECT/all_batch/fig/normal_training_radius_5_cluster_0_2_deviance.png', p, width = 7, height = 5)  
#   
#   
#   training_validation_error = get_relative_error(training_devi / training_num, validation_devi / validation_num, 'training vs validation')
#   training_all1_error = get_relative_error(training_devi / training_num, test1_devi / test1_num, 'training vs test1')
#   training_all2_error = get_relative_error(training_devi / training_num, test2_devi / test1_num, 'training vs test2')
#   final_dataframe = rbind(training_validation_error, training_all1_error, training_all2_error)
#   
#   p = ggplot(data=final_dataframe, aes(x=shuffle_num, y=relative_error, group=type, color=type)) + geom_line()
#   p = p + xlab("Shuffle Number")
#   p = p + ylab("Percent Difference")
#   # p = p + ylim(c(-2, 90))
#   p 
#   ggsave('/data3/CancerECT/all_batch/fig/normal_training_radius_5_cluster_0_2_normalized_deviance_error.png', p, width = 7, height = 5)  
#   
# }
# 
get_avg_particle_num = function(tomograms){
  particle_num_total = 0
  tomogram_num = length(tomograms)
  for (i in 1:tomogram_num){
    particle_num_total = particle_num_total + tomograms[[i]]$n
  }
  return(floor(particle_num_total / tomogram_num))
}
# 

print('Generating quadrature scheme...')

data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
args = commandArgs(TRUE)
tissue = args[1]
intensity_type = args[2]
n = as.numeric(args[3])
r = as.numeric(args[4])
img_idx = args[6]
if (tissue == 'LI' | tissue == 'SI'){
  TMC = 'Stanford'
} else{
  TMC = 'Florida'
}
print(paste(tissue, intensity_type, n, r, sep = ' '))
load(file.path(data_dir, TMC, tissue, paste('cell_list_', n, '_', intensity_type, '_across3_', img_idx, '.Rda', sep = '')))
#load(file.path(data_dir, TMC, tissue, 'img_name_list_all.Rda'))
P = cell_list
print(paste('image num = ', length(P), sep = ''))
avg_particle_num = get_avg_particle_num(P)
interaction = list()
for (i in 1:length(P)) {
  num = length(unique(P[[i]]$marks))
  ir_matrix = matrix(rep(r, num*num), nrow = num)
  hr = as.numeric(args[5])
  hr_matrix = matrix(rep(hr, num*num), nrow = num)
  interaction[[i]] = MultiStraussHard(iradii = ir_matrix, hradii = hr_matrix)
}
H = hyperframe(Y = P)
I = hyperframe(Interaction = interaction)
Quad_all = mppm.quad.ppp(Y ~ marks, data=H, interaction = I, average_number = avg_particle_num, run_quad = F)
filename = file.path(data_dir, TMC, tissue, paste('quad_', n, '_', r, '_', hr, '_', intensity_type, '_across3_', img_idx ,'_self_dummy_grid_eps_50.Rda', sep = ''))
save(Quad_all, file=filename)
gc()
print(warnings())


  

