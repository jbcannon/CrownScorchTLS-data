reflectance = function(int) -25 + 4.577804e-4*int 
intensity = function(refl) 54601.111 + 2184.451*refl

#Add reflectance if it is missing
add_reflectance = function(las) {
  cols = colnames(las@data)
  if('Reflectance' %in% cols) return(las)
  ref = reflectance(las$Intensity)
  las = lidR::add_lasattribute(las, ref, 'Reflectance', 'Reflectance')
  return(las)  
}

SORWithCloudCompare = function(
    infile, outfile=NULL, 
    path_to_cloudcompare = 'C:\\Program Files\\CloudCompareStereo\\CloudCompare.exe') {

  if(is.null(outfile)) outfile = gsub('.las', '_SOR.las', infile)
  if(file.exists(outfile)) {
    warning(outfile, ' already exists')
    return(outfile)
  }
  in1_tmp = tempfile(fileext = '.las')
  out1_tmp = gsub('.las', '_SOR.las', in1_tmp)
  file.copy(infile, in1_tmp)
  cmd = paste0('"', path_to_cloudcompare, '"', " -SILENT -o -GLOBAL_SHIFT AUTO ",
               in1_tmp, " -C_EXPORT_FMT LAS -NO_TIMESTAMP -SOR 6 1")
  cmd = gsub('/', '\\\\', cmd)
  shell(cmd)
  file.copy(out1_tmp, outfile)
  unlink(in1_tmp)
  unlink(out1_tmp)
  return(outfile)
}

# Plot both trees to make sure they are registered, separately, then together
plotRegistered = function(las0, las1) {
  p1 = plot(las0, color="Intensity")
  p2 = plot(las1, color="Intensity")
  p3 = plot(las0, color="Intensity")
  plot(las1, add=p3)
  readline('type n to continue')
}

# return density histogram with defined breaks
get_histogram = function(las, breaks) {
  ref = las$Reflectance
  ref = ref[ref < max(breaks)]
  ref = ref[ref > min(breaks)]
  histogram = hist(ref, breaks=breaks, plot = FALSE)
  histogram = data.frame(intensity = histogram$mids,
                         density = histogram$density)
  return(histogram)
}


tree_locs_from_treeLS = function(boles){
  map = try(TreeLS::treeMap(boles))
  if(class(map)=='try-error') return(NULL)
  boles = TreeLS::treePoints(boles, map, method=TreeLS::trp.crop())
  boles = TreeLS::stemPoints(boles)
  dbh_algo = TreeLS::shapeFit(shape='cylinder', algorithm = 'ransac', n=100, 
                              n_best = 5, inliers=.95, z_dev=10)
  tree_locs = TreeLS::tlsInventory(boles, hp = .99, d_method = dbh_algo)
  tree_locs =   dplyr::select(dplyr::mutate(tree_locs, X=PX,Y=PY,Z=PZ),
                              TreeID, X, Y, Z, Radius, Error)
  return(tree_locs)
}

remove_stem = function(las) {
  las$Z = las$Z - quantile(las$Z,0.001)
  las = TreeLS::stemPoints(las)
  las = lidR::filter_poi(las, !Stem)
  las = lidR::filter_poi(las, Z>1)
  return(las)
}

alignWithCloudCompare = function(
    fn_reference, fn_data, fn_output = NULL,
    path_to_cloudcompare = 'C:\\Program Files\\CloudCompare\\CloudCompare.exe',
    errorDifference = '0.00001', rotation=c('NONE','XYZ')) {
  
  stopifnot(file.exists(fn_reference))
  stopifnot(file.exists(fn_data))
  stopifnot(file.exists(path_to_cloudcompare))
  rotation = match.arg(rotation)
  if(is.null(fn_output)) fn_output = gsub('.las', '_REG.las', fn_data)
  if(file.exists(fn_output)) {
    warning('`fn_output` already exists')
    return(fn_output)
  }
  
  # Create local temp copies of needed files
  ref_fileext = paste0('.',tools::file_ext(fn_reference))
  reference_basename = tempfile(fileext = ref_fileext)
  data_fileext = paste0('.',tools::file_ext(fn_data))
  data_basename = tempfile(fileext = data_fileext)
  file.copy(fn_reference, reference_basename)
  file.copy(fn_data, data_basename)
  
  # Iterative registration
  cmd = paste0('"', path_to_cloudcompare, '"', " -SILENT -o -GLOBAL_SHIFT AUTO ",
               data_basename, " -o -GLOBAL_SHIFT FIRST ", 
               reference_basename, " -C_EXPORT_FMT LAS -NO_TIMESTAMP -ICP -FARTHEST_REMOVAL -ROT ", rotation, " -MIN_ERROR_DIFF ", errorDifference, " -RANDOM_SAMPLING_LIMIT 100000")
  cmd = gsub('/', '\\\\', cmd)
  shell(cmd)
  
  # write good files, and delete temporary files
  file.copy(gsub(data_fileext, paste0('_REGISTERED', '.las'), data_basename), fn_output)
  patt = gsub('.las|.laz', '', paste0(basename(data_basename), '|', basename(reference_basename)))
  unlink(list.files(dirname(data_basename), pattern=patt, full.names=TRUE))
  return(fn_output)
}

make_clean_boles = function(las_norm, int.threshold=41000, vert.threshold=0.8) {
  # Clean up las
  las = lidR::classify_noise(las_norm, lidR::ivf(0.25, 3))
  las = lidR::filter_poi(las, Classification != lidR::LASNOISE)
  #normalize las
  boles = lidR::filter_poi(las, Intensity>int.threshold, Z>0.33)
  boles = lidR::add_lasattribute(boles, round(spanner::eigen_metrics(boles)$Verticality,2), 'verticality', 'verticality')
  boles = lidR::filter_poi(boles, verticality > 0.8)
  return(boles)
}