library(spatstat)
data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
# data_dir = '/data2/PPM/HUBMAP_DATA_new'


args = commandArgs(TRUE)
tissue = args[1]
intensity_type = args[2]
num = as.numeric(args[3])
r = as.numeric(args[4])
hr = as.numeric(args[5])

# tissue = 'LN'
# intensity_type = 'total'
# num = 5
# r = 150
# hr = 1

if (tissue == 'LI' | tissue == 'SI'){
  TMC = 'Stanford'
} else{
  TMC = 'Florida'
}


# load(file = file.path(data_dir, TMC, tissue, paste('quad_', num, '_50_', hr, '_', intensity_type, '_across3.Rda', sep = '')))
# Quad_all$moadf = NULL
# Quad_all_all = Quad_all
# rm(Quad_all)
# gc()
for (i in 1:30){
  print(i)
  quad_image_file = file.path(data_dir, TMC, tissue, paste('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_', i, '_self_dummy_grid_eps_50.Rda', sep = ''))
  if (file.exists(quad_image_file)){
    load(file = quad_image_file)
    Quad_all$moadf$pattern_ID = i
    if (i == 1){
      Quad_all_all = Quad_all
    } else{
      Quad_all_all$moadf = rbind(Quad_all_all$moadf, Quad_all$moadf)
      Quad_all_all$Info$rownames = c(Quad_all_all$Info$rownames, as.character(i))
      rownames(Quad_all$Inter$interaction) = as.character(i)
      Quad_all_all$Inter$interaction = rbind.hyperframe(Quad_all_all$Inter$interaction, Quad_all$Inter$interaction)
      Quad_all_all$npat = Quad_all_all$npat + 1
      rownames(Quad_all$data) = as.character(i)
      Quad_all_all$data = rbind.hyperframe(Quad_all_all$data, Quad_all$data)
      levels(Quad_all_all$data$id) = c(levels(Quad_all_all$data$id), as.character(i))
      Quad_all_all$data$id[i] = i
      Quad_all_all$Y[[i]] = Quad_all$Y[[1]]
      levels(Quad_all_all$datadf$id) = c(levels(Quad_all_all$datadf$id), as.character(i))
      Quad_all_all$datadf = rbind(Quad_all_all$datadf, i)
    }
    rm(Quad_all)
    gc()
  }
}

save(Quad_all_all, file =file.path(data_dir, TMC, tissue, paste('quad_', num, '_', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_50.Rda', sep = '')))
