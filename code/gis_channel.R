
# setup -------------------------------------------------------------------

rm(list = ls())
pacman::p_load(sf,
               raster,
               whitebox,
               mapview,
               tidyverse,
               foreach)

source(here::here("code/function_arc2d8.R"))


# stream extraction -------------------------------------------------------

## watershed polygon layer
sf_wsd <- st_read(dsn = here::here("data_raw/gis"), 
                  layer = "epsg4326_watershed") %>% 
  dplyr::select(NULL) %>% 
  st_transform(crs = 32654) %>% 
  mutate(wsd_id = seq_len(nrow(.)), # watershed id
         area = units::set_units(st_area(.), "km^2"))

## raster layer - upstream watershed area (upa) and flow direction (fdir)
upa <- raster(here::here("data_raw/gis/epsg4326_upa.tif"))
fdir <- raster(here::here("data_raw/gis/epsg4326_dir.tif"))
fdir <- arc2d8(fdir) # convert flow direction from Arc scheme to D8 scheme

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
a_t <- seq(1, 10, by = 1) # threshold values for stream extraction (unit km^2)

df_chl <- foreach(i = seq_len(length(a_t)),
                  .combine = bind_rows) %do% {
                    
                    # stream grid extraction
                    wbt_extract_streams(flow_accum = file_upa,
                                        output = streamgrid,
                                        threshold = a_t[i])
                    
                    # grid to vector
                    wbt_raster_streams_to_vector(streams = streamgrid,
                                                 d8_pntr = file_dir,
                                                 output = channel)
                    
                    sf_chl <- st_read(dsn = tempdir(), layer = "channel") %>% 
                      st_set_crs(4326) %>% 
                      st_transform(32654) %>% 
                      st_join(sf_wsd) %>% # spatial intersect with watershed polygons
                      mutate(length = units::set_units(st_length(.), "km")) %>% 
                      as_tibble() %>% 
                      mutate(length = as.numeric(length)) %>% 
                      drop_na(wsd_id) %>% 
                      group_by(wsd_id) %>% 
                      filter(n() > 2) %>% # select watersheds more than 2 links
                      group_by(wsd_id) %>% 
                      summarize(rate = fitdistrplus::fitdist(length, "exp")$estimate, # rate parameter
                                mean = 1 / rate, # mean link length
                                tl = sum(length), # total river length
                                p_branch = pexp(1, rate), # branching probability per 1km river distance
                                n_branch = n(), # number of links
                                area = unique(area) # total watershed area A (unit km^2)
                                ) %>% 
                      mutate(a_t = a_t[i])
                    
                    return(sf_chl)
                  }

## save data
save(df_chl, file = here::here("data_fmt/df_chl.RData"))


# analysis ----------------------------------------------------------------

## in the following analysis, watersheds > 300km^2 were selected
## small watersheds may have an issue of estimating the rate parameter (i.e., the number of liks is limited)

## load river network data
load(file = here::here("data_fmt/df_chl.RData"))

df_chl_filtered <- df_chl %>% 
  filter(area > units::set_units(300, "km^2")) %>% # select watersheds > 300 km^2
  group_by(wsd_id) %>%
  filter(n() == n_distinct(.$a_t)) %>% 
  ungroup() %>% 
  mutate(wsd_id = factor(wsd_id),
         p_r = n_branch / tl) # branching ration ~ rate parameter

## power law scaling between rate and A_T  
fit <- lme4::lmer(log(rate) ~ 
                    log(a_t) + 
                    (1 + log(a_t) | wsd_id),
                  df_chl_filtered)

z <- coef(fit)$wsd_id %>% 
  mutate(wsd_id = as.factor(rownames(.))) %>% 
  as_tibble() %>% 
  dplyr::select(ln_k = `(Intercept)`, # scaling constant
                z = `log(a_t)`, # scaling exponent
                wsd_id) %>% 
  mutate(k = exp(ln_k),
         mu_z = mean(z)) %>% 
  arrange(desc(k)) %>% 
  mutate(rank = seq_len(nrow(.)))

df_m <- df_chl_filtered %>% 
  left_join(z, by = "wsd_id") %>% 
  mutate(scl_rate = rate / (a_t)^z)

df_1km <- df_m %>% 
  filter(a_t == 1)

cor(df_1km$k, df_1km$scl_rate, method = "spearman")

# figure ------------------------------------------------------------------

df_m %>% 
  ggplot(aes(x = a_t,
             y = scl_rate,
             color = k,
             alpha = 0.05)) +
  geom_point() + 
  geom_line() +
  geom_hline(aes(yintercept = k,
                 color = k)) +
  facet_wrap(facets = ~rank, ncol = 6) +
  labs(x = expression(A[T]),
       y = "Rescaled rate parameter") +
  scale_color_viridis_c() +
  guides(alpha = "none") +
  theme_bw()

