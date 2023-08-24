data_dir = '/home/haoranch/projects/HuBMAP/ppm/HUBMAP_DATA_new'
args = commandArgs(TRUE)
tissue = args[1]
intensity_type = args[2]
n = as.numeric(args[3])
r = as.numeric(args[4])
hr = as.numeric(args[5])

if (tissue == 'LI' | tissue == 'SI'){
  TMC = 'Stanford'
} else{
  TMC = 'Florida'
}
for (radii in seq(100, r, 100)){
  load(file.path(data_dir, TMC, tissue, paste('quad_', n, '_', radii, '_1_total_across3_self_dummy_grid_eps_50.Rda', sep = '')))
  current_quad = Quad_all_all$moadf
  print(colnames(current_quad))
  rm(Quad_all_all)
  gc()
  for (c in (7:(dim(current_quad)[2]-2))){
    colnames(current_quad)[c] = paste(colnames(current_quad)[c], radii, sep = 'x')
  }
  if (radii == 100){
    last_current_quad = current_quad[,7:(dim(current_quad)[2]-2)]
    concatenated_quad = current_quad[,1:(dim(current_quad)[2]-2)]
    pattern_ID = current_quad[,(dim(current_quad)[2]-1)]
    caseweight = current_quad[,dim(current_quad)[2]]
  } else{
    current_quad = current_quad[,7:(dim(current_quad)[2]-2)]
    concatenated_quad = cbind(concatenated_quad, current_quad-last_current_quad)
    last_current_quad = current_quad
  }
}
rm(current_quad)
gc()
concatenated_quad = cbind(concatenated_quad, pattern_ID)
concatenated_quad = cbind(concatenated_quad, caseweight)
load(file.path(data_dir, TMC, tissue, paste('quad_', n, '_', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_50.Rda', sep = '')))
Quad_all_all$moadf = concatenated_quad
save(Quad_all_all, file = file.path(data_dir, TMC, tissue, paste('quad_', n, '_100-', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_50.Rda', sep = '')))
print(colnames(concatenated_quad))
#for (radii in seq(100, r, 100)){
#file.remove(file.path(data_dir, TMC, tissue, paste('quad_', n, '_', r, '_', hr, '_', intensity_type, '_across3_self_dummy_grid_eps_50.Rda', sep = '')))
#}
