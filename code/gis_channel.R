
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
                    
                    df_length <- st_read(dsn = tempdir(), layer = "channel") %>% 
                      st_set_crs(4326) %>% # wgs84
                      st_transform(32654) %>% # utm54n
                      st_join(sf_wsd) %>% # spatial intersect with watershed polygons
                      mutate(length = units::set_units(st_length(.), "km")) %>% 
                      as_tibble() %>% 
                      mutate(length = as.numeric(length)) %>% 
                      drop_na(wsd_id) %>% 
                      group_by(wsd_id) %>% 
                      filter(n() > 2) %>%  # select watersheds more than 2 links
                      ungroup()
                    
                    save(df_length,
                         file = here::here(paste0("data_fmt/df_length_wsd", i, ".RData")))
                    
                    df_summary <- df_length %>% 
                      group_by(wsd_id) %>% 
                      summarize(rate = fitdistrplus::fitdist(length, "exp")$estimate, # rate parameter
                                mean = 1 / rate, # mean link length
                                tl = sum(length), # total river length
                                p_branch = pexp(1, rate), # branching probability per 1km river distance
                                n_branch = n(), # number of links
                                p_r = n_branch / tl, # branching ration ~ rate parameter
                                area = unique(area) # total watershed area A (unit km^2)
                                ) %>% 
                      mutate(a_t = a_t[i])
                    
                    return(df_summary)
                  }

## save data
save(df_chl, file = here::here("data_fmt/df_chl.RData"))

