# 总体分析：主模型、补充检验、环境风险及保护优先级
local({
# 在 RStudio 打开本脚本，点击 Source；也可以在仓库根目录用 Rscript 运行。
paths <- unlist(lapply(sys.frames(), function(f) f$ofile), use.names = FALSE)
args <- grep("^--file=", commandArgs(FALSE), value = TRUE)
paths <- c(rev(paths), sub("^--file=", "", args))
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  paths <- c(paths, tryCatch(rstudioapi::getSourceEditorContext()$path,
                            error = function(e) ""))
}
paths <- paths[!is.na(paths) & nzchar(paths)]
dirs <- dirname(normalizePath(gsub("~+~", " ", paths, fixed = TRUE), mustWork = FALSE))
candidates <- unique(c(dirs, dirname(dirs), getwd(), dirname(getwd()),
                       getOption("chapter1.project_folder")))
has_data <- function(p) all(file.exists(file.path(p, "data",
                                      c("Final_data.csv", "Supporting_data.rds"))))
found <- candidates[vapply(candidates, has_data, logical(1))]
if (!length(found) && interactive()) {
  message("请选择仓库 data 文件夹中的 Final_data.csv。")
  selected <- file.choose()
  candidate <- dirname(dirname(normalizePath(selected)))
  if (has_data(candidate)) found <- candidate
}
if (!length(found)) stop("未找到 data/Final_data.csv 和 data/Supporting_data.rds，请完整下载仓库。")
root <- normalizePath(found[1])
data_dir <- file.path(root, "data")
output_dir <- file.path(root, "outputs")
dir.create(output_dir, showWarnings = FALSE)
options(chapter1.project_folder = root)
old_wd <- setwd(output_dir)
on.exit(setwd(old_wd), add = TRUE)
# Mac 使用自带 Quartz；避免未安装 XQuartz 时的 Cairo 报错。
if (Sys.info()[["sysname"]] == "Darwin") options(bitmapType = "quartz")
support <- readRDS(file.path(data_dir, "Supporting_data.rds"))
suppressPackageStartupMessages({
    library(metafor)
    library(readr)
    library(dplyr)
    library(stringr)
    library(tidyr)
    library(tibble)
    library(terra)
    library(sf)
})

# 1. 从唯一的主数据表取出 Evidence 和已纳入的 Meta。
data <- read_csv(file.path(data_dir, "Final_data.csv"), show_col_types = FALSE, guess_max = Inf)
common <- c("reference_no", "pollinator_species", "publish_year")
meta <- data %>% filter(meta_included == 1) %>%
select(all_of(common), starts_with("meta_")) %>%
rename_with(~sub("^meta_", "", .x), starts_with("meta_")) %>% select(-included)
stopifnot(!anyDuplicated(meta$effect_id), !anyNA(meta$effect_id))
stopifnot(setequal(meta$effect_id, support$meta_order))
meta <- meta[match(support$meta_order, meta$effect_id), ]
evidence <- data %>% select(all_of(common), start_year, end_year, country, region,
  latitude, longitude, pollinator_order, pollinator_family, study_method, trend,
  drivers_reported, continent, starts_with("evidence_")) %>%
rename_with(~sub("^evidence_", "", .x), starts_with("evidence_")) %>%
rename(effect_id = record_id)
time_data <- left_join(meta, support$time_audit, by = "effect_id")
time_data$annual_yi <- with(time_data, ifelse(annual_status == "paired_cumulative",
  yi / audit_duration, ifelse(annual_status == "native_annual", yi, NA_real_)))
time_data$annual_vi <- with(time_data, ifelse(annual_status == "paired_cumulative" &
  variance_type == "source_sampling_or_posterior_variance", vi / audit_duration^2, vi))
time_data$annual_vi_all_scaled <- with(time_data,
  ifelse(annual_status == "paired_cumulative", vi / audit_duration^2, vi))
stopifnot(nrow(meta) == 251, nrow(evidence) == 595,
  all(is.finite(meta$yi)), all(is.finite(meta$vi) & meta$vi > 0))
cat("Evidence:", nrow(evidence), "rows; Meta:", nrow(meta), "effects /",
  n_distinct(meta$reference_no), "sources\n")

# 2. 主模型、分组结果及权重敏感性。
rand<-list(~1|reference_no,~1|country,~1|pollinator_family)
clean<-function(d){
  # D419 is a regional projection, not a single-country sample. Keep original country
  # missing in the deliverable; use one explicitly named regional random-effect group.
  regional<-is.na(d$country)|d$country==''|d$country=='NA'
  stopifnot(all(d$reference_no[regional]=='D419'))
  d$country[regional]<-'Regional_D419_Southeast_Europe'
  d$study_method_clean<-case_when(grepl('historical',d$study_method,ignore.case=TRUE)&grepl('field',d$study_method,
    ignore.case=TRUE)~'Historical + Field',grepl('historical',d$study_method,ignore.case=TRUE)~'Historical',
    grepl('field',d$study_method,ignore.case=TRUE)~'Field',TRUE~d$study_method)
  d$continent<-str_to_title(d$continent)
  d$continent[d$continent%in%c('South America','Central-South America')|d$country%in%c('Costa Rica','Panama',
    'Dominican Republic','Brazil','Chile','Ecuador')]<-'Central-South America'
  d$continent[d$country%in%c('Canada','Canade','USA','USA, Canada','Canada/USA','Mexico')]<-'North America'
  stopifnot(all(is.finite(d$yi)),all(is.finite(d$vi)),all(d$vi>0))
  d
}
new<-clean(meta)
ans<-list()
logs<-list()
coeff<-list()
membership<-list()
models<-list()
run<-function(id,d,mods=~1,random=rand,v='vi'){
  cat(id,nrow(d),'\n')
  flush.console()
  warn<-character()
  d$.v<-d[[v]]
  fitfun<-function(control=list())withCallingHandlers(rma.mv(yi,.v,mods=mods,random=random,data=d,method='REML',
    control=control),warning=function(w){warn<<-c(warn,conditionMessage(w))
      invokeRestart('muffleWarning')})
  fit<-tryCatch(fitfun(),error=function(e)e)
  if(inherits(fit,'error'))fit<-tryCatch(fitfun(list(optimizer='optim',optmethod='BFGS')),error=function(e)e)
  membership[[id]]<<-data.frame(model=id,effect_id=d$effect_id,reference_no=d$reference_no)
  if(inherits(fit,'error')){logs[[id]]<<-data.frame(model=id,status='failed',message=conditionMessage(fit))
    return(NULL)}
  logs[[id]]<<-data.frame(model=id,status=if(length(warn))'warning' else 'completed',message=paste(unique(warn),
    collapse='; '))
  b<-as.numeric(fit$b[1])
  se<-fit$se[1]
  pi<-if(length(fit$b)==1)predict(fit,pi.type='t') else NULL
  vv<-c(fit$sigma2,rep(NA_real_,3))[1:3]
  ans[[id]]<<-data.frame(model=id,k=nrow(d),publications=n_distinct(d$reference_no),estimate=b,SE=se,ci_lb=fit$ci.lb[1],
    ci_ub=fit$ci.ub[1],p=fit$pval[1],pct=100*expm1(b),pct_lb=100*expm1(fit$ci.lb[1]),pct_ub=100*expm1(fit$ci.ub[1]),
    pi_lb=if(is.null(pi))NA else pi$pi.lb,pi_ub=if(is.null(pi))NA else pi$pi.ub,Q=fit$QE,Q_df=fit$k-fit$p,
    Q_p=fit$QEp,I2_proxy=100*sum(fit$sigma2)/(sum(fit$sigma2)+mean(d$.v)),variance_publication=vv[1],variance_country=vv[2],
    variance_family=vv[3],QM=fit$QM,QM_df=fit$QMdf[1],QM_p=fit$QMp)
  coeff[[id]]<<-data.frame(model=id,term=rownames(fit$b),estimate=as.numeric(fit$b),SE=fit$se,ci_lb=fit$ci.lb,
    ci_ub=fit$ci.ub,p=fit$pval)
  models[[id]]<<-fit
  fit
}
for(s in 'updated'){
  d<-new
  main<-run(paste0(s,'_overall'),d)
  for(v in c('pollinator_order','continent','study_method_clean')){
    z<-d[!is.na(d[[v]])&d[[v]]!='',]
    run(paste(s,'QM',v,sep='_'),z,mods=reformulate(v))
    for(g in unique(z[[v]])){q<-z[z[[v]]==g,]
      if(nrow(q)>=2)run(paste(s,'sub',v,g,sep='_'),q)}
  }
  long<-bind_rows(lapply(seq_len(nrow(d)),function(i){drv<-trimws(strsplit(d$drivers_reported[i],',')[[1]])
      drv<-unique(drv[drv%in%c('Multiple','Habitat','Climate','Urbanization','Pesticides','Pathogens','Invasives')])
      if(!length(drv))return(NULL)
      q<-d[rep(i,length(drv)),]
      q$driver<-drv
      q}))
  run(paste0(s,'_QM_driver'),long,mods=~driver)
  for(g in unique(long$driver))if(sum(long$driver==g)>1)run(paste(s,'sub_driver',g,sep='_'),long[long$driver==g,
    ])
  run(paste0(s,'_random_publication'),d,random=list(~1|reference_no))
  run(paste0(s,'_random_publication_country'),d,random=list(~1|reference_no,~1|country))
  d$vi_cap<-1/pmin(1/d$vi,quantile(1/d$vi,.95))
  d$vi_equal<-1
  run(paste0(s,'_weight_cap95'),d,v='vi_cap')
  run(paste0(s,'_weight_equal'),d,v='vi_equal')
  eg<-lm(residuals(main)~sqrt(d$vi))
  egger_residual<-coef(summary(eg))
  ag<-bind_rows(lapply(split(d,d$reference_no),function(x){if(nrow(x)==1)return(data.frame(reference_no=x$reference_no,
      yi=x$yi,vi=x$vi))
      m<-rma(yi,vi,data=x,method='REML')
      data.frame(reference_no=x$reference_no[1],yi=as.numeric(m$b),vi=m$se^2)}))
  st<-rma(yi,vi,data=ag,method='REML')
  tf<-trimfill(st)
  et<-regtest(st,model='rma',predictor='sei')
  diagnostics<-list(study=st,trimfill=tf,aggregate=ag,main=main)
  bias_results<-data.frame(model=c('study_level','trimfill'),k=c(st$k,tf$k),estimate=c(st$b,tf$b),ci_lb=c(st$ci.lb,
    tf$ci.lb),ci_ub=c(st$ci.ub,tf$ci.ub),p=c(st$pval,tf$pval),missing=c(0,tf$k0),egger_z=et$zval,egger_p=et$pval)
  inf<-influence(st)
  influence_results<-cbind(reference_no=ag$reference_no,as.data.frame(inf$inf))
  leave_one_out<-cbind(reference_no=ag$reference_no,as.data.frame(leave1out(st)))
}
# 按响应和提取方式建立敏感性分析子集。
txt<-tolower(paste(new$response,new$unit,new$time_definition))
extraction_txt<-tolower(paste(new$extraction_method,new$effect_source_notes))
flags<-data.frame(effect_id=new$effect_id,reference_no=new$reference_no,
  acute=new$reference_no=='H165',habitat=grepl('suitability|suitable area',txt),
  digitized=(!is.na(new$digitized)&new$digitized==1)|grepl('digitiz|digitis|graph estimate|read from.*fig|figure-derived',
  extraction_txt),
  assemblage=grepl('spp\\.|assemblage|complex',new$pollinator_species,ignore.case=TRUE),
  spatial=grepl('spatial comparison|canopy.*ground|patch comparison',txt),
  composition=grepl('composition|relative frequency|relative abundance|proportion of all',txt),
  zero_corrected=grepl('0\\.5 correction|0\\.5 continuity|a0\\.5|continuity correction',tolower(new$effect_source_notes)),

  proxy_variance=abs(new$vi-1/new$sample_size)<1e-10)
for(nm in c('acute','habitat','digitized','assemblage','spatial','composition','zero_corrected')){
  ex<-flags[[nm]]
  if(any(ex))run(paste0('updated_exclude_',nm),new[!ex,])
}
primary <- list(results = bind_rows(ans), coefficients = bind_rows(coeff),
  membership = bind_rows(membership), status = bind_rows(logs), models = models,
  flags = flags, diagnostics = diagnostics, bias = bias_results,
  egger_residual = egger_residual, influence = influence_results,
  leave_one_out = leave_one_out)
driver_data <- long

# 3. 年化与响应定义的补充检验。
d <- time_data
d$country[is.na(d$country)|d$country==''|d$country=='NA']<-'Regional_D419_Southeast_Europe'
ans<-list()
coefs<-list()
members<-list()
logs<-list()
models<-list()
run<-function(id,x,response='yi',variance='vi',mods=~1){
  cat(id,nrow(x),'\n')
  flush.console()
  warnings<-character()
  x$.y<-x[[response]]
  x$.v<-x[[variance]]
  stopifnot(nrow(x)>3,all(is.finite(x$.y)),all(is.finite(x$.v)&x$.v>0))
  fitfun<-function(control=list())withCallingHandlers(rma.mv(.y,.v,mods=mods,random=list(~1|reference_no,
    ~1|country,~1|pollinator_family),data=x,method='REML',control=control),warning=function(w){warnings<<-c(warnings,
      conditionMessage(w))
      invokeRestart('muffleWarning')})
  fit<-tryCatch(fitfun(),error=function(e)e)
  if(inherits(fit,'error'))fit<-fitfun(list(optimizer='optim',optmethod='BFGS'))
  b<-as.numeric(fit$b[1])
  v<-c(fit$sigma2,rep(NA_real_,3))[1:3]
  ans[[id]]<<-data.frame(model=id,k=nrow(x),sources=n_distinct(x$reference_no),estimate=b,SE=fit$se[1],
    CI_low=fit$ci.lb[1],CI_high=fit$ci.ub[1],z=fit$zval[1],p=fit$pval[1],percentage_change=100*expm1(b),pct_low=100*expm1(fit$ci.lb[1]),
    pct_high=100*expm1(fit$ci.ub[1]),Q=fit$QE,Q_df=fit$k-fit$p,Q_p=fit$QEp,I2_proxy=100*sum(fit$sigma2)/(sum(fit$sigma2)+mean(x$.v)),
    variance_publication=v[1],variance_country=v[2],variance_family=v[3],QM=fit$QM,QM_df=fit$QMdf[1],QM_p=fit$QMp,
    AIC=AIC(fit),response=response,variance=variance)
  coefs[[id]]<<-data.frame(model=id,term=rownames(fit$b),estimate=as.numeric(fit$b),SE=fit$se,p=fit$pval,
    CI_low=fit$ci.lb,CI_high=fit$ci.ub)
  members[[id]]<<-data.frame(model=id,effect_id=x$effect_id,reference_no=x$reference_no)
  logs[[id]]<<-data.frame(model=id,status=if(length(warnings))'warning' else 'completed',message=paste(unique(warnings),
    collapse='; '))
  models[[id]]<<-fit
  invisible(fit)
}
p<-d[d$annual_status=='paired_cumulative',]
annual<-d[d$annual_status%in%c('paired_cumulative','native_annual'),]
stopifnot(max(abs(p$annual_yi*p$audit_duration-p$yi))<1e-10)
run('S1_complete',d)
run('S1_paired_cumulative',p)
run('S1_paired_annual',p,'annual_yi','annual_vi')
run('S1_expanded_annual',annual,'annual_yi','annual_vi')
run('S1_QM_cumulative_duration',p,mods=~log(audit_duration))
run('S1_QM_paired_annual_duration',p,'annual_yi','annual_vi',~log(audit_duration))
run('S1_QM_expanded_annual_duration',annual,'annual_yi','annual_vi',~log(audit_duration))
for(cut in c(1,5)){
  run(paste0('S1_cumulative_gt',cut),p[p$audit_duration>cut,])
  run(paste0('S1_paired_annual_gt',cut),p[p$audit_duration>cut,],'annual_yi','annual_vi')
  run(paste0('S1_expanded_annual_gt',cut),annual[annual$audit_duration>cut,],'annual_yi','annual_vi')
}
run('S1_no_acute',d[d$acute_class!='confirmed_acute',])
run('S1_no_acute_or_abrupt',d[!d$acute_class%in%c('confirmed_acute','abrupt_step'),])
run('S1_annual_all_scaled_variances',p,'annual_yi','annual_vi_all_scaled')
run('S1_exact_cumulative',p[p$approximate==0,])
run('S1_exact_annual',p[p$approximate==0,],'annual_yi','annual_vi')
# Response/extraction rules apply to corrected responses, not superseded encodings.
txt<-tolower(paste(d$response,d$unit,d$time_definition,d$effect_source_notes))
f<-data.frame(effect_id=d$effect_id,reference_no=d$reference_no,
  habitat=d$reference_no%in%c('H67','H1341','D419'),acute=d$acute_class=='confirmed_acute',abrupt=d$acute_class=='abrupt_step',

  assemblage=is.na(d$pollinator_species)|d$pollinator_species=='NA'|grepl('spp\\.|assemblage',d$pollinator_species,
  ignore.case=TRUE),
  species_complex=d$reference_no=='H13',
  digitized=(!is.na(d$digitized)&d$digitized==1)|grepl('digitiz|digitis|bar.height|read from.*fig|figure.derived|graph estimate',
  txt),
  inferred_N=d$reference_no%in%c('H165','C30','C1221','L818','D50','D1682'),
  N_scope=d$reference_no%in%c('H165','D13','C8'),
  transformed=d$reference_no=='H23',native_annual=d$reference_no%in%c('H13','H137'),
  non_temporal=d$reference_no%in%c('C103','C1822','L1254','H2732','NS0126','H80','H181','D419'),
  composition=d$reference_no%in%c('H185','NS513'),genetic=d$reference_no=='C1472',mixed_activity=d$reference_no=='NS410')
f[is.na(f)]<-FALSE
ex<-list(M01=f$habitat,M02=f$acute,M03=f$acute|f$abrupt,M04=f$assemblage,M05=f$assemblage|f$species_complex,
  M06=f$N_scope,M07=f$N_scope|f$inferred_N,M08=f$digitized,M09=f$digitized|f$transformed,M10=f$habitat|f$acute|f$assemblage|f$species_complex|f$N_scope|f$digitized,
  M11=f$transformed,M12=f$transformed|f$native_annual,M13=f$non_temporal,M14=f$habitat|f$genetic|f$composition|f$mixed_activity|d$reference_no=='NS0172',
  M15=f$habitat|f$genetic|f$composition|f$mixed_activity|f$transformed|f$native_annual|d$reference_no=='NS0172')
run('S2_M00',d)
for(id in names(ex))run(paste0('S2_',id),d[!ex[[id]],])
pre<-bind_rows(meta,support$cycle_meta)
pre<-pre[match(support$pre_cycle_order,pre$effect_id),]
pre$country[is.na(pre$country)]<-'Regional_D419_Southeast_Europe'
run('S2_including_cycle_comparisons',pre)
z<-d[!grepl('^N',d$reference_no),]
run('S4_legacy_sources_only',z)
run('S4_no_original_supplementary_trap_activity',d[is.na(d$trap_or_activity)|d$trap_or_activity!=1,])
run('S4_no_Fox_report',d[d$reference_no!='NEW_FOX2015',])
supplement <- list(results = bind_rows(ans), coefficients = bind_rows(coefs),
  membership = bind_rows(members), status = bind_rows(logs), flags = f, models = models)

# 4. 空间风险和保护优先级；保留已采用的300次文献bootstrap。
layers <- lapply(support$rasters, terra::unwrap)
risk <- local({
    set.seed(20260513)
    n_boot <- 300
    aggregate_factor <- 20  # 0.05° * 20 = 1° analysis grid
    high_risk_quantile <- 0.75

    # ------------------------------------------------------------
    # 1) Helpers
    # ------------------------------------------------------------
    clean_order <- function(x) {
      x <- str_to_lower(str_squish(as.character(x)))
      case_when(
      x == "coleoptera"  ~ "Coleoptera",
      x == "diptera"     ~ "Diptera",
      x == "hymenoptera" ~ "Hymenoptera",
      x == "lepidoptera" ~ "Lepidoptera",
      TRUE ~ str_to_title(x)
    )
  }

    clean_method <- function(x) {
      x <- str_to_lower(str_squish(as.character(x)))
      case_when(
      str_detect(x, "field") & str_detect(x, "historical") ~ "Historical + Field",
      str_detect(x, "historical") ~ "Historical",
      str_detect(x, "field") ~ "Field",
      TRUE ~ NA_character_
    )
  }

    std_z <- function(x) {
      x <- as.numeric(x)
      sx <- sd(x, na.rm = TRUE)
      if (!is.finite(sx) || sx == 0) return(rep(NA_real_, length(x)))
      as.numeric(scale(x))
  }

    norm01_vec <- function(x) {
      rng <- range(x, na.rm = TRUE)
      if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(NA_real_, length(x)))
      (x - rng[1]) / diff(rng)
  }

    extract_raster_at_points <- function(raster_path, pts_ll, out_name) {
      r <- raster_path
      if (is.na(crs(r)) || crs(r) == "") crs(r) <- "EPSG:4326"

      pts <- vect(pts_ll, geom = c("longitude", "latitude"), crs = "EPSG:4326")
      if (!same.crs(pts, r)) pts <- project(pts, crs(r))

      val <- terra::extract(r, pts, method = "bilinear")
      tibble(effect_id = pts_ll$effect_id, !!out_name := as.numeric(val[[2]]))
  }

    fit_fast_weights <- function(dat) {
      dat <- dat %>%
      filter(
      is.finite(yi), is.finite(vi), vi > 0,
      is.finite(pesticide_z), is.finite(habitat_z), is.finite(climate_z),
      !is.na(pollinator_order), !is.na(study_method_clean)
    )

      if (nrow(dat) < 40 || n_distinct(dat$reference_no) < 10) {
        return(c(Pesticides = NA_real_, Habitat = NA_real_, Climate = NA_real_))
    }

      fit <- try(
      lm(
      yi ~ pesticide_z + habitat_z + climate_z + pollinator_order + study_method_clean,
      data = dat,
      weights = 1 / vi
    ),
      silent = TRUE
    )

      if (inherits(fit, "try-error")) {
        return(c(Pesticides = NA_real_, Habitat = NA_real_, Climate = NA_real_))
    }

      b <- coef(fit)[c("pesticide_z", "habitat_z", "climate_z")]
      names(b) <- c("Pesticides", "Habitat", "Climate")

      # In lnRR, more negative means stronger decline. Convert coefficients to risk direction.
      # If a coefficient is positive in a bootstrap sample, that layer does not increase
      # predicted decline risk in that sample, so its risk weight is set to zero.
      risk_b <- pmax(-as.numeric(b), 0)

      # If all signs are positive in a bootstrap sample, fall back to absolute size to avoid
      # losing the sample. This is rare but keeps the analysis numerically stable.
      if (sum(risk_b, na.rm = TRUE) == 0) risk_b <- abs(as.numeric(b))

      w <- risk_b / sum(risk_b, na.rm = TRUE)
      names(w) <- c("Pesticides", "Habitat", "Climate")
      w
  }


    # ------------------------------------------------------------
    # 2) Build effect-size exposure table for bootstrap weights
    # ------------------------------------------------------------
    meta <- meta %>%
    mutate(
    effect_id = as.character(effect_id),
    reference_no = as.character(reference_no),
    pollinator_order = clean_order(pollinator_order),
    study_method_clean = clean_method(study_method),
    latitude = suppressWarnings(as.numeric(latitude)),
    longitude = suppressWarnings(as.numeric(longitude)),
    yi = suppressWarnings(as.numeric(yi)),
    vi = suppressWarnings(as.numeric(vi))
  ) %>%
    filter(
    is.finite(yi), is.finite(vi), vi > 0,
    is.finite(latitude), is.finite(longitude),
    latitude >= -90, latitude <= 90,
    longitude >= -180, longitude <= 180
  )

    pts_ll <- meta %>% select(effect_id, longitude, latitude)

    exposure <- pts_ll %>%
    left_join(extract_raster_at_points(layers$pest, pts_ll, "pesticide"), by = "effect_id") %>%
    left_join(extract_raster_at_points(layers$habitat, pts_ll, "habitat"), by = "effect_id") %>%
    left_join(extract_raster_at_points(layers$climate, pts_ll, "climate"), by = "effect_id")

    dat <- meta %>%
    left_join(exposure %>% select(-longitude, -latitude), by = "effect_id") %>%
    mutate(
    pesticide_z = std_z(pesticide),
    habitat_z = std_z(habitat),
    climate_z = std_z(climate)
  )


    # Bootstrap at publication level
    pubs <- unique(dat$reference_no)
    weights_boot <- replicate(n_boot, {
      sampled_pubs <- sample(pubs, length(pubs), replace = TRUE)
      boot_dat <- bind_rows(lapply(sampled_pubs, function(id) dat %>% filter(reference_no == id)))
      fit_fast_weights(boot_dat)
  })

    weights_boot <- t(weights_boot) %>% as.data.frame()
    weights_boot <- weights_boot %>%
    mutate(iter = row_number()) %>%
    filter(if_all(c(Pesticides, Habitat, Climate), is.finite))

    weights_summary <- weights_boot %>%
    summarise(
    n_boot_valid = n(),
    pesticide_mean = mean(Pesticides),
    pesticide_sd = sd(Pesticides),
    habitat_mean = mean(Habitat),
    habitat_sd = sd(Habitat),
    climate_mean = mean(Climate),
    climate_sd = sd(Climate)
  )


    cat("\n== Bootstrap weight summary ==\n")
    print(weights_summary)

    # ------------------------------------------------------------
    # 3) Make 1° pressure grid and Kew/WCVP orchid mask
    # ------------------------------------------------------------
    r_pest <- aggregate(layers$pest, fact = aggregate_factor, fun = mean, na.rm = TRUE)
    r_hab  <- aggregate(layers$habitat, fact = aggregate_factor, fun = mean, na.rm = TRUE)
    r_clim <- aggregate(layers$climate, fact = aggregate_factor, fun = mean, na.rm = TRUE)

    names(r_pest) <- "Pesticides"
    names(r_hab) <- "Habitat"
    names(r_clim) <- "Climate"

    pressure <- c(r_pest, r_hab, r_clim)

    grid_df <- as.data.frame(pressure, xy = TRUE, na.rm = FALSE) %>%
    rename(lon = x, lat = y)

    kew <- support$kew %>%
    st_make_valid() %>%
    st_transform(4326) %>%
    mutate(
    orchid_hotspot = log1p(as.numeric(total)) * as.numeric(pc_threat_pred),
    orchid_present_kew = ifelse(is.finite(orchid_hotspot) & orchid_hotspot > 0, 1L, NA_integer_)
  ) %>%
    select(orchid_hotspot, orchid_present_kew)

    grid_sf <- st_as_sf(grid_df, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
    grid_kew <- st_join(grid_sf, kew, join = st_within, left = TRUE) %>%
    st_drop_geometry()

    grid_mat <- grid_kew %>%
    mutate(
    Pesticides = norm01_vec(Pesticides),
    Habitat = norm01_vec(Habitat),
    Climate = norm01_vec(Climate),
    in_orchid_region = !is.na(orchid_present_kew)
  )

    # ------------------------------------------------------------
    # 4) Predict risk and confidence/uncertainty per grid cell
    # ------------------------------------------------------------
    pressure_mat <- as.matrix(grid_mat[, c("Pesticides", "Habitat", "Climate")])
    valid_cell <- grid_mat$in_orchid_region & rowSums(is.na(pressure_mat)) == 0

    risk_store <- matrix(NA_real_, nrow = nrow(pressure_mat), ncol = nrow(weights_boot))

    for (i in seq_len(nrow(weights_boot))) {
      w <- as.numeric(weights_boot[i, c("Pesticides", "Habitat", "Climate")])
      risk_i <- pressure_mat %*% w
      risk_i <- as.numeric(risk_i)
      risk_i[!valid_cell] <- NA_real_
      risk_i <- norm01_vec(risk_i)
      risk_store[, i] <- risk_i
  }

    risk_mean <- rowMeans(risk_store, na.rm = TRUE)
    risk_sd <- apply(risk_store, 1, sd, na.rm = TRUE)

    high_flags <- apply(risk_store, 2, function(x) {
      thr <- as.numeric(quantile(x, high_risk_quantile, na.rm = TRUE))
      as.numeric(x >= thr)
  })
    confidence_high_risk <- rowMeans(high_flags, na.rm = TRUE)

    out_grid <- grid_mat %>%
    transmute(
    lon, lat,
    predicted_risk = risk_mean,
    risk_uncertainty_sd = risk_sd,
    confidence_high_risk = confidence_high_risk,
    orchid_hotspot
  )


    # Convert back to rasters for GeoTIFF outputs
    template <- r_pest
    risk_r <- rast(out_grid[, c("lon", "lat", "predicted_risk")], type = "xyz", crs = "EPSG:4326")
    sd_r   <- rast(out_grid[, c("lon", "lat", "risk_uncertainty_sd")], type = "xyz", crs = "EPSG:4326")
    conf_r <- rast(out_grid[, c("lon", "lat", "confidence_high_risk")], type = "xyz", crs = "EPSG:4326")

    risk_r <- resample(risk_r, template, method = "near")
    sd_r <- resample(sd_r, template, method = "near")
    conf_r <- resample(conf_r, template, method = "near")



    stopifnot(nrow(weights_boot)>0,all(is.finite(as.matrix(weights_boot[,c("Pesticides","Habitat","Climate")]))))
    list(exposure = dat, weights_boot = weights_boot, weights_summary = weights_summary,
    grid = out_grid, risk_r = risk_r, confidence_r = conf_r)

})
priority <- local({
    priority_levels <- c("Very Low", "Low", "Medium", "High", "Very High")
    cols_priority <- c(
    "Very Low" = "#bdbdbd",
    "Low" = "#cbddf0",
    "Medium" = "#9ecae1",
    "High" = "#fc9272",
    "Very High" = "#de2d26"
  )

    # ------------------------------------------------------------
    # 1) Old-style 2-degree grid
    # ------------------------------------------------------------
    bbox <- st_as_sfc(
    st_bbox(c(xmin = -180, xmax = 180, ymin = -60, ymax = 85), crs = 4326)
  )

    grid2 <- st_sf(geometry = st_make_grid(bbox, cellsize = c(2, 2), square = TRUE)) %>%
    mutate(
    index = row_number(),
    centroid = st_centroid(geometry),
    lon_bin = st_coordinates(centroid)[, 1],
    lat_bin = st_coordinates(centroid)[, 2]
  ) %>%
    st_set_crs(4326) %>%
    select(index, lon_bin, lat_bin, geometry)

    # ------------------------------------------------------------
    # 2) Kew/WCVP orchid hotspot layer
    # ------------------------------------------------------------
    kew <- support$kew %>%
    st_make_valid() %>%
    st_transform(4326) %>%
    mutate(
    ml_index = log1p(as.numeric(total)) * as.numeric(pc_threat_pred),
    kew_total = as.numeric(total),
    kew_threat = as.numeric(pc_threat_pred),
    tdwg_name = as.character(LEVEL3_NAM),
    tdwg_code = as.character(area_code_l3)
  )

    grid_cent <- st_as_sf(
    grid2 %>% st_drop_geometry() %>% select(index, lon_bin, lat_bin),
    coords = c("lon_bin", "lat_bin"),
    crs = 4326
  )

    grid_join <- st_join(grid_cent, kew, join = st_within, left = TRUE) %>%
    st_drop_geometry()

    hot_grid_ml <- grid2 %>%
    st_drop_geometry() %>%
    left_join(
    grid_join %>%
    select(index, tdwg_name, tdwg_code, kew_total, kew_threat, ml_index),
    by = "index"
  ) %>%
    transmute(
    index, lon_bin, lat_bin,
    tdwg_name, tdwg_code,
    kew_total = coalesce(kew_total, 0),
    kew_threat,
    hot_ml = ifelse(is.na(ml_index), 0, ml_index)
  )

    # ------------------------------------------------------------
    # 3) Evidence layer from MainData.xlsx
    # ------------------------------------------------------------
    main_dat <- evidence %>%
    mutate(latitude=as.numeric(latitude),longitude=as.numeric(longitude),reference_no=as.character(reference_no)) %>%
    filter(is.finite(latitude),is.finite(longitude),between(latitude,-90,90),between(longitude,-180,180))
    evidence_grid <- st_as_sf(main_dat, coords = c("longitude", "latitude"), crs = 4326) %>%
    st_join(grid2, join = st_intersects) %>%
    st_drop_geometry() %>%
    distinct(reference_no, index, .keep_all = TRUE) %>%
    group_by(index) %>%
    summarise(evidence = n(), .groups = "drop")

    # ------------------------------------------------------------
    # 4) Research mismatch layer using old code logic
    # ------------------------------------------------------------
    world <- support$world %>%
    st_transform(4326) %>%
    st_make_valid()
    coast <- NULL
    world_union <- st_union(world)

    mismatch_data <- grid2 %>%
    left_join(hot_grid_ml %>% select(index, hot_ml), by = "index") %>%
    left_join(evidence_grid, by = "index") %>%
    mutate(
    hot = coalesce(hot_ml, 0),
    evidence = coalesce(evidence, 0)
  ) %>%
    filter(!(hot == 0 & evidence == 0)) %>%
    filter(lengths(st_intersects(., world_union)) > 0) %>%
    mutate(
    z_hot = scale(log1p(hot))[, 1],
    z_evi = scale(log1p(evidence))[, 1],
    M = z_hot - z_evi,
    bin = cut(
    M,
    breaks = c(-Inf, -2, -1, 1, 2, Inf),
    labels = c("<= -2σ (over-studied)", "-2σ~ -1σ", "-1σ~ +1σ",
    "+1σ~ +2σ", ">= +2σ (under-studied)"),
    include.lowest = TRUE,
    right = TRUE,
    ordered_result = TRUE
  )
  )

    gap_df <- mismatch_data %>%
    st_drop_geometry() %>%
    select(index, M, bin) %>%
    mutate(
    gap3 = case_when(
    bin %in% c(">= +2σ (under-studied)", "+1σ~ +2σ") ~ "Under",
    bin %in% c("<= -2σ (over-studied)", "-2σ~ -1σ") ~ "Over",
    TRUE ~ "Balanced"
  )
  )

    # ------------------------------------------------------------
    # 5) Risk layers: predicted vs confidence, aggregated to 2 degrees
    # ------------------------------------------------------------
    extract_grid_mean <- function(path, name) {
      r <- path
      x <- terra::extract(r, vect(grid2), fun = mean, na.rm = TRUE)
      names(x)[1:2] <- c("index", name)
      x[[name]][is.na(x[[name]])] <- 0
      x
  }

    risk_df <- extract_grid_mean(risk$risk_r, "predicted_risk") %>%
    left_join(extract_grid_mean(risk$confidence_r, "confidence_risk"), by = "index")

    # ------------------------------------------------------------
    # 6) Priority logic copied from old PriorityMap.r
    # ------------------------------------------------------------
    build_priority <- function(risk_col, risk_layer_label) {
      hot_vals <- hot_grid_ml$hot_ml[hot_grid_ml$hot_ml > 0]
      hot_q <- quantile(hot_vals, probs = c(0.5, 0.75), na.rm = TRUE)
      names(hot_q) <- c("med", "high")

      risk_vals <- risk_df[[risk_col]][risk_df[[risk_col]] > 0]
      risk_q <- quantile(risk_vals, probs = c(1/3, 2/3), na.rm = TRUE)
      names(risk_q) <- c("low_high", "med_high")

      grid2 %>%
      left_join(hot_grid_ml, by = "index") %>%
      left_join(evidence_grid, by = "index") %>%
      left_join(gap_df %>% select(index, gap3), by = "index") %>%
      left_join(risk_df %>% select(index, all_of(risk_col)), by = "index") %>%
      mutate(
      hot = coalesce(hot_ml, 0),
      evid = coalesce(evidence, 0),
      risk = coalesce(.data[[risk_col]], 0),
      B_class = case_when(
      hot >= hot_q["high"] ~ "High",
      hot >= hot_q["med"] ~ "Medium",
      hot > 0 ~ "Low",
      TRUE ~ "None"
    ),
      R_class = case_when(
      risk >= risk_q["med_high"] ~ "High",
      risk >= risk_q["low_high"] ~ "Medium",
      risk > 0 ~ "Low",
      TRUE ~ "None"
    ),
      gap3 = recode(gap3, .missing = "Under"),
      BasePriority = case_when(
      B_class == "High"   & R_class == "High"   ~ "Very High",
      B_class == "High"   & R_class == "Medium" ~ "High",
      B_class == "High"   & R_class == "Low"    ~ "Medium",
      B_class == "High"   & R_class == "None"   ~ "Low",
      B_class == "Medium" & R_class == "High"   ~ "High",
      B_class == "Medium" & R_class == "Medium" ~ "Medium",
      B_class == "Medium" & R_class == "Low"    ~ "Low",
      B_class == "Medium" & R_class == "None"   ~ "Low",
      B_class == "Low"    & R_class == "High"   ~ "Medium",
      B_class == "Low"    & R_class == "Medium" ~ "Low",
      B_class == "Low"    & R_class == "Low"    ~ "Very Low",
      B_class == "Low"    & R_class == "None"   ~ "Very Low",
      B_class == "None"   & R_class == "High"   ~ "Low",
      B_class == "None"   & R_class == "Medium" ~ "Low",
      B_class == "None"   & R_class == "Low"    ~ "Very Low",
      B_class == "None"   & R_class == "None"   ~ "Very Low",
      TRUE ~ "Very Low"
    ),
      BaseRank = case_when(
      BasePriority == "Very Low" ~ 1L,
      BasePriority == "Low" ~ 2L,
      BasePriority == "Medium" ~ 3L,
      BasePriority == "High" ~ 4L,
      BasePriority == "Very High" ~ 5L,
      TRUE ~ 1L
    ),
      AdjRank = case_when(
      gap3 == "Under" ~ pmin(BaseRank + 1L, 5L),
      gap3 == "Over" ~ pmax(BaseRank - 1L, 1L),
      TRUE ~ BaseRank
    ),
      AdjRank = case_when(
      gap3 == "Under" & risk > 0 & AdjRank < 2L ~ 2L,
      TRUE ~ AdjRank
    ),
      Priority_chr = case_when(
      AdjRank == 1L ~ "Very Low",
      AdjRank == 2L ~ "Low",
      AdjRank == 3L ~ "Medium",
      AdjRank == 4L ~ "High",
      AdjRank == 5L ~ "Very High",
      TRUE ~ "Very Low"
    ),
      Priority_chr = ifelse(hot == 0 & risk == 0 & evid == 0, NA_character_, Priority_chr),
      Priority = factor(Priority_chr, levels = priority_levels),
      risk_layer = risk_layer_label
    )
  }

    normalise_priority_coords <- function(x) {
      x %>%
      mutate(
      lon_bin = coalesce(.data[["lon_bin.x"]], .data[["lon_bin.y"]]),
      lat_bin = coalesce(.data[["lat_bin.x"]], .data[["lat_bin.y"]])
    ) %>%
      select(-any_of(c("lon_bin.x", "lat_bin.x", "lon_bin.y", "lat_bin.y")))
  }

    priority_predicted <- build_priority("predicted_risk", "Predicted environmental risk") %>%
    normalise_priority_coords()
    priority_confidence <- build_priority("confidence_risk", "Confidence high-risk layer") %>%
    normalise_priority_coords()

    priority_all <- bind_rows(
    priority_predicted %>% st_drop_geometry(),
    priority_confidence %>% st_drop_geometry()
  )


    priority_summary <- priority_all %>%
    filter(!is.na(Priority)) %>%
    count(risk_layer, Priority, name = "n_cells") %>%
    group_by(risk_layer) %>%
    mutate(percent = 100 * n_cells / sum(n_cells)) %>%
    ungroup()



    mask <- support$mask
    z<-st_drop_geometry(priority_confidence)
    z$keep<-mask$strict_presence_decision_kingsley[match(z$index,mask$index)]=='keep'
    z$rank_masked<-ifelse(z$keep,z$AdjRank,NA)

    global_grid <- z
    lin <- support$lineage
    risk_df<-risk_df[,c("index","confidence_risk")]
    build_group_grid<-function(group_name){z<-lin[lin$group==group_name,]
      z<-z[!duplicated(z$index),]
      z[,c('index','lon_bin','lat_bin','richness','hot_ml')]}
    build_gap <- function(group_grid) {
      combo <- grid2 %>%
      left_join(group_grid %>% select(index, hot_ml), by = "index") %>%
      left_join(evidence_grid, by = "index") %>%
      mutate(
      hot = replace_na(hot_ml, 0),
      evidence = replace_na(evidence, 0)
    ) %>%
      filter(!(hot == 0 & evidence == 0)) %>%
      filter(lengths(st_intersects(., world_union)) > 0)

      if (nrow(combo) == 0 || length(unique(combo$hot)) < 2 || length(unique(combo$evidence)) < 2) {
        return(combo %>% st_drop_geometry() %>% transmute(index, M = NA_real_, gap3 = "Balanced"))
    }

      combo %>%
      mutate(
      z_hot = as.numeric(scale(log1p(hot))),
      z_evi = as.numeric(scale(log1p(evidence))),
      M = z_hot - z_evi,
      bin = cut(
      M,
      breaks = c(-Inf, -2, -1, 1, 2, Inf),
      labels = c("<= -2sigma (over-studied)", "-2sigma to -1sigma",
      "-1sigma to +1sigma", "+1sigma to +2sigma",
      ">= +2sigma (under-studied)"),
      include.lowest = TRUE,
      right = TRUE,
      ordered_result = TRUE
    ),
      gap3 = case_when(
      bin %in% c(">= +2sigma (under-studied)", "+1sigma to +2sigma") ~ "Under",
      bin %in% c("<= -2sigma (over-studied)", "-2sigma to -1sigma") ~ "Over",
      TRUE ~ "Balanced"
    )
    ) %>%
      st_drop_geometry() %>%
      select(index, M, gap3)
  }

    build_priority <- function(group_name, group_type_label) {
      group_grid <- build_group_grid(group_name)
      gap_df <- build_gap(group_grid)

      hot_vals <- group_grid$hot_ml[group_grid$hot_ml > 0]
      hot_q <- quantile(hot_vals, probs = c(0.5, 0.75), na.rm = TRUE)
      names(hot_q) <- c("med", "high")

      risk_vals <- risk_df$confidence_risk[risk_df$confidence_risk > 0]
      risk_q <- quantile(risk_vals, probs = c(1/3, 2/3), na.rm = TRUE)
      names(risk_q) <- c("low_high", "med_high")

      grid2 %>%
      left_join(group_grid, by = "index") %>%
      left_join(evidence_grid, by = "index") %>%
      left_join(gap_df %>% select(index, gap3), by = "index") %>%
      left_join(risk_df, by = "index") %>%
      mutate(
      group = group_name,
      group_type = group_type_label,
      hot = replace_na(hot_ml, 0),
      richness = replace_na(richness, 0),
      evid = replace_na(evidence, 0),
      risk = replace_na(confidence_risk, 0),
      B_class = case_when(
      hot >= hot_q["high"] ~ "High",
      hot >= hot_q["med"] ~ "Medium",
      hot > 0 ~ "Low",
      TRUE ~ "None"
    ),
      R_class = case_when(
      risk >= risk_q["med_high"] ~ "High",
      risk >= risk_q["low_high"] ~ "Medium",
      risk > 0 ~ "Low",
      TRUE ~ "None"
    ),
      gap3 = recode(gap3, .missing = "Under"),
      BasePriority = case_when(
      B_class == "High"   & R_class == "High"   ~ "Very High",
      B_class == "High"   & R_class == "Medium" ~ "High",
      B_class == "High"   & R_class == "Low"    ~ "Medium",
      B_class == "High"   & R_class == "None"   ~ "Low",
      B_class == "Medium" & R_class == "High"   ~ "High",
      B_class == "Medium" & R_class == "Medium" ~ "Medium",
      B_class == "Medium" & R_class == "Low"    ~ "Low",
      B_class == "Medium" & R_class == "None"   ~ "Low",
      B_class == "Low"    & R_class == "High"   ~ "Medium",
      B_class == "Low"    & R_class == "Medium" ~ "Low",
      B_class == "Low"    & R_class == "Low"    ~ "Very Low",
      B_class == "Low"    & R_class == "None"   ~ "Very Low",
      B_class == "None"   & R_class == "High"   ~ "Low",
      B_class == "None"   & R_class == "Medium" ~ "Low",
      B_class == "None"   & R_class == "Low"    ~ "Very Low",
      B_class == "None"   & R_class == "None"   ~ "Very Low",
      TRUE ~ "Very Low"
    ),
      BaseRank = case_when(
      BasePriority == "Very Low" ~ 1L,
      BasePriority == "Low" ~ 2L,
      BasePriority == "Medium" ~ 3L,
      BasePriority == "High" ~ 4L,
      BasePriority == "Very High" ~ 5L,
      TRUE ~ 1L
    ),
      AdjRank = case_when(
      gap3 == "Under" ~ pmin(BaseRank + 1L, 5L),
      gap3 == "Over" ~ pmax(BaseRank - 1L, 1L),
      TRUE ~ BaseRank
    ),
      AdjRank = case_when(
      gap3 == "Under" & risk > 0 & AdjRank < 2L ~ 2L,
      TRUE ~ AdjRank
    ),
      Priority_chr = case_when(
      AdjRank == 1L ~ "Very Low",
      AdjRank == 2L ~ "Low",
      AdjRank == 3L ~ "Medium",
      AdjRank == 4L ~ "High",
      AdjRank == 5L ~ "Very High",
      TRUE ~ "Very Low"
    ),
      # Lineage-specific maps should only colour cells where the focal
      # subfamily/genus has WCVP native distribution support. Risk/evidence
      # layers can modify priority inside the lineage range, but must not
      # create priority cells outside that range.
      Priority_chr = ifelse(richness <= 0, NA_character_, Priority_chr),
      Priority = factor(Priority_chr, levels = priority_levels)
    ) %>%
      mutate(
      lon_bin = coalesce(.data[["lon_bin.x"]], .data[["lon_bin.y"]]),
      lat_bin = coalesce(.data[["lat_bin.x"]], .data[["lat_bin.y"]])
    ) %>%
      select(-any_of(c("lon_bin.x", "lat_bin.x", "lon_bin.y", "lat_bin.y")))
  }
    lineages<-bind_rows(lapply(unique(lin$group),function(g){x<-st_drop_geometry(build_priority(g,lin$group_type[match(g,
      lin$group)]))
      m<-lin[lin$group==g,]
      x$final_keep<-m$final_keep[match(x$index,m$index)]
      x$rank_masked<-ifelse(x$final_keep,x$AdjRank,NA)
      x}))
    list(global_masked = global_grid, lineages = lineages, all_classes = priority_all,
    summary = priority_summary, evidence_grid = evidence_grid,
    mismatch_grid = st_drop_geometry(mismatch_data))

})

# 5. 所有结果存到一个文件；各张图直接读取它。
model_ok <- risk$exposure %>% filter(is.finite(yi), is.finite(vi), vi > 0,
  is.finite(pesticide_z), is.finite(habitat_z), is.finite(climate_z),
  !is.na(pollinator_order), !is.na(study_method_clean))
regions <- priority$global_masked %>% left_join(
  st_drop_geometry(support$kew) %>% select(area_code_l3, LEVEL1_COD),
  by = c("tdwg_code" = "area_code_l3")) %>%
mutate(continent = case_when(LEVEL1_COD == 1 ~ "Europe", LEVEL1_COD == 2 ~ "Africa",
  LEVEL1_COD %in% c(3, 4) ~ "Asia", LEVEL1_COD %in% c(5, 6) ~ "Oceania",
  LEVEL1_COD == 7 ~ "North America", LEVEL1_COD == 8 ~ "South America", TRUE ~ "Other")) %>%
filter(is.finite(rank_masked)) %>% group_by(continent) %>%
summarise(cells = n(), very_high = sum(rank_masked == 5),
  high_or_very_high = sum(rank_masked >= 4), .groups = "drop")
risk$risk_r <- risk$confidence_r <- NULL
results <- list(meta = meta, evidence = evidence, driver_data = driver_data,
  time_data = time_data, primary = primary, supplement = supplement,
  risk = risk, priority = priority, regions = regions, risk_membership = model_ok,
  input_md5 = unname(tools::md5sum(file.path(data_dir, "Final_data.csv"))),
  supporting_md5 = unname(tools::md5sum(file.path(data_dir, "Supporting_data.rds"))),
  session = capture.output(sessionInfo()))
# 各模型的数值、成员和诊断均已保存；不重复保存几十份大矩阵。
results$primary$models <- results$supplement$models <- NULL
detach_formula_environment <- function(x) {
  if (inherits(x, "formula")) { environment(x) <- baseenv()
    return(x) }
  if (is.list(x)) {
    a <- attributes(x)
    y <- unclass(x)
    for (i in seq_along(y)) y[i] <- list(detach_formula_environment(y[[i]]))
    attributes(y) <- a
    x <- y
  }
  x
}
results$primary$diagnostics <- detach_formula_environment(results$primary$diagnostics)
saveRDS(results, "Results.rds", compress = "gzip")
summary <- bind_rows(
  primary$results %>% select(model, k, publications, estimate, SE, ci_lb, ci_ub, p, QM, QM_p),
  supplement$results %>% transmute(model, k, publications = sources, estimate, SE,
  ci_lb = CI_low, ci_ub = CI_high, p, QM, QM_p))
write_csv(summary, "Results.csv")
print(summary %>% filter(model %in% c("updated_overall", "updated_QM_pollinator_order")))
print(primary$bias)
cat("分析完成：Results.csv 为结果摘要，Results.rds 供各图脚本读取。\n")
})
