# Estimate crown scorch from lidar intensity returns

# Script to read in all manually segmented trees,
# create histogram of intensity scores
# break into segments and use Random Forests to predict
# scorch. Visualize histograms from different scorch patterns

library(tidyverse)
library(randomForest)
library(betareg)
source('helper_fxns.R')

# Search and clean segmented tree list
tree_dir = 'data/manual-clip-trees/'
histogram_dir = 'data/histograms-manual/'
scorch_data = "data/mortality_scorch_survey - data.csv"
histograms_file = 'data/compiledHistogram-manual.csv'
#histogramBreaks = seq(10000,55000,500)
#pretty(reflectance(histogramBreaks), 100)
histogramBreaks = seq(-20,0, by = 0.2)


update = FALSE

# Create directories if they don't exist
if(!dir.exists(tree_dir)) dir.create(tree_dir, recursive = TRUE)
if(!dir.exists(histogram_dir)) dir.create(histogram_dir, recursive = TRUE)

# From clipped trees, compile database of all files
fn = list.files(tree_dir, full.names = TRUE)
fn = lapply(fn, function(i) {
  treeid = gsub('.laz|.las|_pre|_post', '', basename(i))
  path = i
  tile = substr(treeid, 0, 4)
  phase = ifelse(grepl('pre', basename(i)), 'pre', 'post')
  df = data.frame(treeid = treeid, tile=tile, phase=phase, path = i)
  return(df)
}); fn = do.call(rbind, fn)

# tabular scorch data
scorch_ocular = read_csv(scorch_data)
s = scorch_ocular$`% SCORCH` 
s = ifelse(s == '<1', 0.5, s)
s = ifelse(s == '1-5', mean(c(1,5)), s)
s = ifelse(s == '5-10', mean(c(5,10)), s)
s = as.numeric(s)
scorch_ocular$`% SCORCH` = s
if(exists('s')) rm(s)
scorch_ocular = mutate(scorch_ocular, UNQ = paste0(TILE, '-', TAG_ID))

# Generate histograms of scanned trees
make_plot = FALSE
for(t in (unique(fn$treeid))) {
  tree_info = fn %>% filter(treeid == t)
  # ensure they have a match
  if(nrow(tree_info) < 2) {
    warning('no match found for', t, 'skipping\n')
    next
  }
  
  # ensure it hasn't already been completed
  histogram_csv_fn = paste0(histogram_dir, '/', t, '.csv')
  if(file.exists(histogram_csv_fn)) {
    cat('histogram', t, 'already complete\n')
    next
  }
  
  # Make histogram for individual trees
  tree_info = tree_info %>% arrange(desc(phase))
  his = lapply(tree_info$path, function(fn) {
    las = lidR::readLAS(fn)
    las = add_reflectance(las)
    las = remove_stem(las)
    his = get_histogram(las, breaks=histogramBreaks)
    return(his)
  })
  
  # combine histograms and get difference in intensity
  h0 = his[[1]]
  h0$dens1 = his[[2]]$density
  colnames(h0) = c('intensity', 'dens0', 'dens1')
  h0$delta = h0$dens1 - h0$dens0
  h0 = h0[, c('intensity', 'dens0', 'dens1', 'delta')]
  h0$treeid = t
  h0$method='manual clip'
  
  # plot and write to csv
  if(make_plot) h0 %>% ggplot(aes(y = delta, x = intensity)) + geom_line() + theme_bw()
  write_csv(h0, histogram_csv_fn)
  cat('tree', t, 'histogram complete\n')
}

# Compile all histograms from individual lidar scans
if(update) {
  hist_fns = list.files(histogram_dir, pattern='.csv$', full.names=TRUE)
  compiledHistograms = do.call(rbind, lapply(hist_fns, function(f) {
    t = gsub('.csv', '', basename(f))
    x = read_csv(f, show_col_types=FALSE) 
    x$treeid = t
    return(x)
  }))
  write_csv(compiledHistograms, histograms_file)  
}
compiledHistograms = read_csv(histograms_file)

# Generate predictors for Random Forests using histogram_breaks
out = lapply(unique(compiledHistograms$treeid), function(id) {
  treeHistogram = filter(compiledHistograms, treeid==id)
  tmp = t(treeHistogram[c('dens1', 'delta')])
  tmp = c(tmp[1,], tmp[2,])
  vars = c('intensity', 'delta')
  colnames = c(paste0(vars[1], '_', round(treeHistogram$intensity,1)),
               paste0(vars[2], '_', round(treeHistogram$intensity,1)))
  tmp = data.frame(values=tmp, names=colnames)
  tmp = t(tmp)
  colnames(tmp) = tmp[2,]
  tmp = as.data.frame(t(tmp[1,]))
  scorch = as.numeric(filter(scorch_ocular, UNQ == id)$`% SCORCH`)
  if(length(scorch)==0) scorch = NA
  dbh = as.numeric(filter(scorch_ocular, UNQ == id)$`DBH_CM`)
  if(length(dbh)==0) dbh = NA
  scorchHeight = filter(scorch_ocular, UNQ == id)$`SCORCH \nHEIGHT`
  scorchHeight = ifelse(scorchHeight == 'X', NA, as.numeric(scorchHeight))
  if(length(scorchHeight)==0) scorchHeight=NA
  db = data.frame(scorch=scorch, sorchHt = scorchHeight)
  out = cbind(db, tmp)
  out = cbind(data.frame(id=id, dbh = dbh), out)
  return(out)  
}); out = do.call(rbind,out)
rownames(out)=NULL
all_data = out
for(c in 2:ncol(all_data)) {
  all_data[,c] = as.numeric(trimws(all_data[,c]))
}

if(update) {
  
  #intialize
  if(any(all_data$scorch>1)) all_data$scorch = all_data$scorch/100
  seed = 781699
  set.seed(seed)
  
  # Random forest model (Intensity)
  data.int = select(all_data, contains(c('scorch', 'intensity')))
  impvars.int = Boruta::Boruta(data.int[,-1], data.int[,1], maxRuns=250)
  impvars.int = names(impvars.int$finalDecision)[impvars.int$finalDecision == 'Confirmed']
  RF_scorch_int = randomForest::randomForest(y = data.int[,1],
                                             x = data.int[,impvars.int],
                                             ntree=250,
                                             keep.forest=TRUE,
                                             mtry=5,
                                             corr.bias=TRUE); print(RF_scorch_int)
  #save that data for output
  df.rf.int = data.frame(id = all_data$id,
                         dbh = all_data$dbh,
                         observed = data.int[,'scorch'], 
                         predicted = predict(RF_scorch_int),
                         model = 'Random forests',
                         response = 'Intensity')
  
  # Beta regression model (Intensity)
  data.int$sc = data.int$scorch
  data.int$sc[data.int$scorch == 0] = 0.001
  data.int$sc[data.int$scorch == 1] = 0.999
  f_int = paste0('sc ~ ', paste0('`', impvars.int, '`', collapse = ' + '))
  train_rows = sample(1:nrow(data.int), 0.7*ceiling(nrow(data.int)))
  br_mod_int = betareg(f_int, data=data.int[train_rows,])
  
  df.br.int = data.frame(id = all_data$id[-train_rows],
                         dbh = all_data$dbh[-train_rows],
                         observed = data.int[-train_rows,1], 
                         predicted = predict(br_mod_int, data.int[-train_rows,-1]),
                         model = 'Beta regression',
                         response = 'Intensity')
  
  # Random forest model (delta Intensity)
  data.delt = select(all_data, contains(c('scorch', 'delta')))
  impvars.delt = Boruta::Boruta(data.delt[,-1], data.delt[,1], maxRuns=1000)
  impvars.delt = names(impvars.delt$finalDecision)[impvars.delt$finalDecision == 'Confirmed']
  set.seed(seed)
  RF_scorch_delta = randomForest::randomForest(y = data.delt[,1],
                                             x = data.delt[,impvars.delt],
                                             ntree=250,
                                             keep.forest=TRUE,
                                             mtry=5,
                                             corr.bias=TRUE); print(RF_scorch_delta)
  df.rf.delt = data.frame(id = all_data$id,
                          dbh = all_data$dbh,
                          observed = data.delt[,1], 
                          predicted = predict(RF_scorch_delta),
                          model = 'Random forests',
                          response = 'Delta intensity')
  
  
  # Beta regression with  Delta intensity
  data.delt$sc = data.int$sc
  f_delt = formula(paste0('sc ~ ', paste0('`', impvars.delt, '`', collapse=' + ')))
  br_mod_delt = betareg(f_delt, data=data.delt[train_rows,])
  df.br.delt = data.frame(id = all_data$id[-train_rows],
                          dbh = all_data$dbh[-train_rows],
                          observed = data.delt[-train_rows,1],
                          predicted = predict(br_mod_delt, data.delt[-train_rows,-1]),
                          model = 'Beta regression',
                          response = 'Delta intensity')
  
  
  # Combine and save prediction outputs
  df = rbind(df.rf.int, df.br.int, df.rf.delt, df.br.delt)
  
  
  write_csv(df, 'data/model-predictions.csv')
  # Output random forest models
  save(RF_scorch_int, RF_scorch_delta, lm_mod_int, lm_mod_delt, file='data/RF_models.Rdata')
}
 
load('data/RF_models.Rdata')
model_predictions = read_csv('data/model-predictions.csv')

model_outcomes = model_predictions %>% group_by(response, model) %>%
    mutate(resid = observed-predicted) %>%
    summarize(mean_diff = mean(resid), 
              rmse = sqrt(mean(resid^2)),
              nrmse = sqrt(mean(resid^2/mean(observed)))*100,
              r2 = summary(lm(observed~predicted))$r.squared,
              smdape = median(200*abs(observed-predicted)/(abs(observed)+abs(predicted))),
              p = summary(lm(predicted~observed))[[4]][2,4])
model_outcomes

write_csv(model_outcomes, 'figs/model-outcomes.csv')

model_outcomes = model_outcomes %>% mutate(
                             lab.r2 = paste0('R^2 == ', round(r2,3)),
                             lab.rmse = paste0('RMSE == ', round(rmse,1)),
                             lab.nrmse = paste0('nRMSE == ', round(nrmse,1),"*\'%\'"),
                             lab.smdape = paste0('SMdAPE == ', round(smdape,1), "*\'%\'"))

fig.allmods = model_predictions %>% ggplot(aes(x = observed, y = predicted)) + 
    geom_point(alpha = 0.6) +
    facet_grid(rows = vars(model), cols = vars(response), scales='free_y') +
    #lims(y = c(0,100), x = c(0,100)) +
    geom_smooth(method = 'lm', se=FALSE, color='dodgerblue3') +
    geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
    theme_bw() + 
    labs(x = expression(Scorch[observed]), y = expression(Scorch[predicted])) +
    geom_text(data = model_outcomes, size = 3, aes(x=0,y=.99, label = lab.r2), parse=TRUE, hjust=0) +
    geom_text(data = model_outcomes, size = 3, aes(x=0,y=.92, label = lab.rmse), parse=TRUE, hjust=0) +
    geom_text(data = model_outcomes, size = 3, aes(x=0,y=.85, label = lab.nrmse), parse=TRUE, hjust=0) +
    geom_text(data = model_outcomes, size = 3, aes(x=0,y=.78, label = lab.smdape), parse=TRUE, hjust=0)

   
fig.allmods
# Model Output figures
ggsave('figs/modelFits.jpg', fig.allmods, width = 8, height = 8, dpi = 600)
ggsave('figs/modelFits.pdf', fig.allmods, width = 8, height = 8, dpi = 600)
  
# Compare estimtes to arkin et al.
model_predictions %>% filter(model == 'Random forests') %>%
  filter(response == 'Delta intensity') %>%
  mutate(difference = predicted-observed,
         abs_diff = abs(difference),
         diff_lt_10per = abs_diff<=0.10) %>% 
  group_by(dbh>50) %>% summarize(sum(diff_lt_10per)/length(diff_lt_10per))

# Variable importance plots
ImportancePlot_int = randomForest::varImpPlot(RF_scorch_int) # -7.7.3 to -6.3 db look really important
ImportancePlot_delt = randomForest::varImpPlot(RF_scorch_delta) # -6.7 to -5.9 db look really important
impdf = rbind(ImportancePlot_int, ImportancePlot_delt)
impdf = data.frame(intensity=rownames(impdf), purity = as.numeric(impdf))
impdf = impdf %>% mutate(model = ifelse(grepl('intensity', intensity), 'Intensity', 'Delta Intensity'))
impdf$intensity = as.numeric(gsub('delta_|intensity_','', impdf$intensity))
impdf = impdf %>% mutate(model = factor(model, 
                                  levels = c('Intensity', 'Delta Intensity'),
                                  labels = c(expression(RF~Model:Intensity), expression(RF~Model:Delta~Intensity))))
fig.varImportance = ggplot(data=impdf, aes(y=purity, x=intensity)) + 
    facet_wrap(~model, labeller = label_parsed) + 
    geom_line() + geom_point(alpha=0.7) +
    theme_bw() + labs(x='Return intensity (dB)', y='Increase in node purity')
  fig.varImportance
  ggsave('figs/varImportanceplot.jpg', width=8, height=4, dpi=600)
  ggsave('figs/varImportanceplot.pdf', width=8, height=4, dpi=600)

##### Example histograms
donetrees = unique(compiledHistograms$treeid)
select_trees = 
  c('M-04-15549' , 'L-05-14669',# 0% scorch
    'E-08-9269', 'B-04-4286',   #  ~50% scorch
    'D-03-10867', 'C-04-11029') # 95% scorch
scorchdb = data.frame(treeid=select_trees, 
                      scorch = as.factor(c('0%','0%','50%','50%','95%','95%')),
                      tree = c('A','B','A','B','A','B'))

fig.exampleHistograms = compiledHistograms %>% 
  filter(treeid %in% select_trees) %>%
  pivot_longer(c(dens0,dens1)) %>%
  left_join(scorchdb) %>% 
  ggplot(aes(x=intensity, color=name, y=value)) +
  geom_line(size=1) + 
  facet_grid(rows = vars(scorch), cols=vars(tree)) +
  theme_bw() + 
  labs(x='Return intensity (dB)', y='Density') +
  scale_color_manual(labels = c('pre-burn', 'post-burn'), values=c('darkgreen', 'orange')) +
  theme(strip.text.x = element_blank()) +
  theme(legend.position= c(0.2, 0.75),
        legend.background = element_blank(),
        legend.title= element_blank()) +
  geom_text(data=scorchdb %>% mutate(name='dens0'), aes(x=-20,y=0.0, label=treeid), color = grey(0.3), hjust=0)
fig.exampleHistograms
ggsave('figs/example_histograms.pdf', width=6, height=8)
ggsave('figs/example_histograms.jpg', width=6, height=8, dpi=600)

# All histograms combined
histdf = scorch_ocular %>% mutate(treeid = UNQ) %>% select(treeid, `% SCORCH`) %>%
  right_join(compiledHistograms) %>%
  mutate( scorch_cls = cut(`% SCORCH`, seq(0,100,25))) %>%
  na.omit


a = ggplot(histdf, aes(x=intensity, y=dens0)) + 
  geom_line(aes(group=treeid), color='grey', alpha=0.5, size=1) +
  geom_smooth(span=0.1, se=FALSE, aes(color=scorch_cls)) + theme_bw() +
  labs(y='Density', x = 'Return intensity (dB)') +
  lims(y=c(0,0.2)) + labs(color = '% Scorch') +
  ggtitle('Pre-burn return intensity') +
  theme(legend.position = c(0.22,0.70), legend.background = element_blank())

b = ggplot(histdf, aes(x=intensity, y=dens1)) + 
  geom_line(aes(group=treeid), color='grey', alpha=0.5, size=1) +
  geom_smooth(span=0.1, se=FALSE, aes(color=scorch_cls)) + theme_bw() +
  labs(y='Density', x = 'Return intensity (db)') +
  lims(y=c(0,0.2)) +
  ggtitle('Post-burn return intensity') +
  theme(legend.position = 'none')

c = ggplot(histdf, aes(x=intensity, y=delta)) + 
  geom_line(aes(group=treeid), color='grey', alpha=0.5, size=1) +
  geom_smooth(span=0.1, se=FALSE, aes(color=scorch_cls)) + theme_bw() +
  geom_hline(yintercept=0, linetype='dashed') +
  labs(y='Density', x = 'Return intensity (dB)') +
  ggtitle(expression(Delta~density[postburn-preburn])) +
  theme(legend.position = 'none') 

fig.combinedHistograms = ggpubr::ggarrange(a,b,c, nrow=1, labels=LETTERS)
fig.combinedHistograms
ggsave('figs/combinedHistograms.jpg', width=8, units='in', height=3, dpi=600)
ggsave('figs/combinedHistograms.pdf', width=8, units='in', height=3, dpi=600)



mod_df = data.frame(mod = c('RF_scorch_int', 'RF_scorch_delta', 'lm_mod_int', 'lm_mod_delt'),
                    response = c('Intensity', 'Delta intensity', 'Intensity', 'Delta intensity'),
                    model = c('Random forest', 'Random forest', 'Linear model', 'Linear model'))



assessment = model_predictions
assessment$dbh_cls = cut(assessment$dbh, breaks=seq(0,60,10))

fig.sizeAssess = ggplot(assessment %>% filter(model=='Random forests'), aes(x=observed, y=predicted)) +
  geom_point(alpha=0.6, size=1) + 
  facet_grid(rows=vars(model, response), cols=vars(dbh_cls)) +
  geom_smooth(method = 'lm', se=FALSE, color='dodgerblue3') +
  coord_fixed() + 
  geom_abline(slope=1, intercept=0, linetype='dashed') +
  theme_bw() + 
  scale_x_continuous(breaks=c(0,0.5,1),limits=0:1) +
  scale_y_continuous(breaks=c(0,0.5,1), limits=0:1) +
  #geom_text(data = mod.assessment, aes(x = 0, y = 102, label = r2.lab), parse=TRUE, hjust=0, size=2.5) + 
  #geom_text(data = mod.assessment, aes(x = 0, y = 92, label = rmse.lab), parse=TRUE, hjust=0, size = 2.5) +
  labs(y = expression(Scorch[predicted]), x = expression(Scorch[observed]))
fig.sizeAssess
ggsave('figs/sizeAssess.jpg', width=8, height=3, dpi=600)
ggsave('figs/sizeAssess.pdf', width=8, height=3, dpi=600)


mod.assessment = assessment %>% group_by(model, response, dbh_cls) %>%
  summarize(mean_diff = mean(predicted-observed), 
            rmse = sqrt(mean((predicted-observed)^2)),
            nrmse = sqrt(mean((predicted-observed)^2/mean(observed)))*100,
            r2 = summary(lm(observed~predicted))$r.squared,
            smdape = median(200*abs(observed-predicted)/(abs(observed)+abs(predicted))))


mod.assessment = mod.assessment %>% mutate(mod = paste0(model, response),
                                           dbh = as.numeric(gsub(']', '', gsub("(.*?)*,", "", dbh_cls)))-5)

mod.assessment %>% filter(model == 'Random forests', dbh_cls == '(50,60]')
mod.assessment %>% filter(model == 'Random forests', dbh_cls == '(0,10]')

a = mod.assessment %>% filter(model == 'Random forests') %>% ungroup %>%
        ggplot(aes(x = dbh_cls, y = r2, color = response, group=response)) +
        geom_point() + geom_line() +
        theme_bw() + labs(y = expression(LOOCV~R^2), x = 'Size class (dbh, cm)') +
        lims(y = c(0,1)) + labs(color=element_blank()) +
        scale_color_brewer(palette='Dark2') +
        theme(legend.position = c(0.25, 0.85), 
        legend.text = element_text(size=8),
        legend.background = element_blank(),
        legend.title = element_text(size=8))
b = mod.assessment %>% filter(model == 'Random forests') %>% ungroup %>% 
  ggplot(aes(x = dbh_cls, y = mean_diff, color = response, group=response)) +
  geom_point() + geom_line() +
  theme_bw() + labs(y = 'Mean difference', x = 'Size class (dbh, cm)') +
  lims(y = c(-0.2,0.2)) + 
  scale_color_brewer(palette='Dark2') +
  theme(legend.position = 'none') +
  geom_hline(yintercept = 0, linetype='dashed')
c = mod.assessment %>% filter(model == 'Random forests') %>% ungroup %>% 
  ggplot(aes(x = dbh_cls, y = rmse, color = response, group=response)) +
  geom_point() + geom_line() +
  theme_bw() + labs(y = 'RMSE', x = 'Size class (dbh, cm)') +
  lims(y = c(0,.25)) + 
  scale_color_brewer(palette='Dark2') +
  theme(legend.position = 'none')
d = mod.assessment %>% filter(model == 'Random forests') %>% ungroup %>% 
  ggplot(aes(x = dbh_cls, y = nrmse, color = response, group=response)) +
  geom_point() + geom_line() +
  theme_bw() + labs(y = 'nRMSE (%)', x = 'Size class (dbh, cm)') +
  lims(y = c(0,100)) + 
  scale_color_brewer(palette='Dark2') +
  theme(legend.position = 'none')
e = mod.assessment %>% filter(model == 'Random forests') %>% ungroup %>% 
  ggplot(aes(x = dbh_cls, y = smdape, color = response, group=response)) +
  geom_point() + geom_line() +
  theme_bw() + labs(y = 'SMdAPE (%)', x = 'Size class (dbh, cm)') +
  lims(y = c(0,100)) + 
  scale_color_brewer(palette='Dark2') +
  theme(legend.position = 'none')

figModassess = ggpubr::ggarrange(
              ggpubr::ggarrange(a,b,labels=LETTERS, nrow=1),
              ggpubr::ggarrange(c,d, e,labels=LETTERS[3:5],nrow=1), nrow=2)
              
figModassess
ggsave('figs/modAssess_by_size.jpg', width=8, height=6, dpi=600)
ggsave('figs/modAssess_by_size.pdf', width=8, height=6, dpi=600)

trees_used = assessment %>% 
  filter(model == 'Random forests', response== 'Intensity') %>%
  select(dbh,observed)


stratification = trees_used %>%
  mutate(dbh_cls = cut(dbh, breaks = seq(0,60,10)),
         sch_cls = cut(observed, breaks = c(-1, seq(0,1,0.1)))) %>%
  group_by(dbh_cls, sch_cls) %>%
  summarize(n = length(observed)) %>%
  pivot_wider(names_from = `dbh_cls`, values_from = n)
write_csv(stratification, 'figs/tree-stratification.csv')


a = trees_used %>% ggplot(aes(x=dbh)) + geom_histogram(color='black') + 
  theme_bw() + labs(y = 'Count', x = 'DBH (cm)')
b = trees_used %>% ggplot(aes(x=observed)) + geom_histogram(color='black') +
  theme_bw() + labs(y = 'Count', x = expression(Scorch[observed]*"(%)"))
c = trees_used %>% ggplot(aes(x=dbh, y=observed)) + geom_point(alpha=0.5) +
  theme_bw() + labs(x = 'DBH (cm)', y = expression(Scorch[observed]*"(%)")) +
  geom_smooth(method = 'lm')
tree_selection.fig = ggpubr::ggarrange(a,b,c, labels=LETTERS, nrow=1)
tree_selection.fig
ggsave('figs/tree-stratification.jpg', height=3, width=9, dpi = 600)
ggsave('figs/tree-stratification.pdf', height=3, width=9, dpi = 600)

trees_used %>% lm(formula = observed ~ dbh) %>% summary()
nrow(trees_used)
range(trees_used$dbh)
range(trees_used$observed)

