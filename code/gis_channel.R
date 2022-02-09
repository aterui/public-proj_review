
# setup -------------------------------------------------------------------

rm(list = ls())
pacman::p_load(sf,
               raster,
               whitebox,
               mapview,
               tidyverse,
               foreach)

source(here::here("review/code/function_arc2d8.R"))

# extraction --------------------------------------------------------------

## watershed polygon layer
sf_wsd <- st_read(dsn = here::here("review/data_raw/gis"), 
                  layer = "epsg4326_watershed") %>% 
  dplyr::select(NULL) %>% 
  st_transform(crs = 32654) %>% 
  mutate(wsd_id = seq_len(nrow(.)),
         area = units::set_units(st_area(.), "km^2"))

## raster layer
upa <- raster(here::here("review/data_raw/gis/epsg4326_upa.tif"))
fdir <- raster(here::here("review/data_raw/gis/epsg4326_dir.tif"))
fdir <- arc2d8(fdir)

## export to temporary directory
file_upa <- paste0(tempdir(), "\\upa_tmp.tif")
file_dir <- paste0(tempdir(), "\\dir_tmp.tif")

writeRaster(upa,
            filename = file_upa,
            overwrite = TRUE)

writeRaster(fdir,
            filename = file_dir,
            overwrite = TRUE)

## objects for watershed analysis
streamgrid <- paste0(tempdir(), "\\streamgrid.tif")
channel <- paste0(tempdir(), "\\channel.shp")
a_t <- seq(1, 15, by = 2)

df_chl <- foreach(i = seq_len(length(a_t)),
                  .combine = bind_rows) %do% {
                    
                    wbt_extract_streams(flow_accum = file_upa,
                                        output = streamgrid,
                                        threshold = a_t[i])
                    
                    wbt_raster_streams_to_vector(streams = streamgrid,
                                                 d8_pntr = file_dir,
                                                 output = channel)
                    
                    sf_chl <- st_read(dsn = tempdir(), layer = "channel") %>% 
                      st_set_crs(4326) %>% 
                      st_transform(32654) %>% 
                      st_join(sf_wsd) %>% 
                      mutate(length = units::set_units(st_length(.), "km")) %>% 
                      as_tibble() %>% 
                      mutate(length = as.numeric(length),
                             x = length / mean(length)) %>% 
                      drop_na(wsd_id) %>% 
                      group_by(wsd_id) %>% 
                      filter(n() > 2) %>% 
                      group_by(wsd_id) %>% 
                      summarize(rate = fitdistrplus::fitdist(length, "exp")$estimate,
                                mean = 1 / rate,
                                tl = sum(length),
                                p_branch = pexp(1, rate),
                                n_branch = n(),
                                area = unique(area)) %>% 
                      mutate(a_t = a_t[i])
                    
                    return(sf_chl)
                  }

save(df_chl, file = here::here("review/data_fmt/df_chl.RData"))

# analysis ----------------------------------------------------------------

## load river network data
load(file = here::here("review/data_fmt/df_chl.RData"))

df_chl_filtered <- df_chl %>% 
  filter(area > units::set_units(300, "km^2")) %>% 
  group_by(wsd_id) %>%
  filter(n() == n_distinct(.$a_t)) %>% 
  ungroup() %>% 
  mutate(wsd_id = factor(wsd_id),
         scl_rate = rate * sqrt(a_t),
         p_r = n_branch / tl)
  
fit <- lme4::lmer(log(rate) ~ 
                    log(a_t) + 
                    (1 + log(a_t)|wsd_id),
                  df_chl_filtered)

z <- coef(fit)$wsd_id %>% 
  mutate(wsd_id = as.factor(rownames(.))) %>% 
  as_tibble() %>% 
  dplyr::select(k = `(Intercept)`,
                z = `log(a_t)`,
                wsd_id) %>% 
  mutate(k = exp(k))

df_m <- df_chl_filtered %>% 
  left_join(z, by = "wsd_id") %>% 
  mutate(scl_rate = rate / (a_t)^z)


# figure ------------------------------------------------------------------

df_m %>% 
  ggplot(aes(x = a_t,
             y = scl_rate,
             color = wsd_id,
             alpha = 0.1)) +
  geom_point() + 
  geom_line() +
  geom_hline(aes(yintercept = k,
                 color = wsd_id)) +
  facet_wrap(facets = ~wsd_id, ncol = 5) +
  theme(legend.position = "none")
  
df_m %>% 
  filter(a_t == 1) %>% 
  select(k, rate) %>% 
  cor(method = "spearman")
