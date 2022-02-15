
# setup -------------------------------------------------------------------

rm(list = ls())
pacman::p_load(tidyverse,
               sf,
               foreach)

# analysis ----------------------------------------------------------------

## in the following analysis, watersheds > 300km^2 were selected
## small watersheds could have an issue of estimating the rate parameter
## i.e., the number of links is limited

## network properties ####
load(file = here::here("data_fmt/df_chl.RData"))

df_chl_filtered <- df_chl %>% 
  filter(area > units::set_units(300, "km^2")) %>% # watersheds > 300km^2 were selected
  group_by(wsd_id) %>%
  filter(n() == n_distinct(.$a_t)) %>% 
  ungroup() %>% 
  mutate(wsd_id = factor(wsd_id)) 

## power law scaling between rate and A_T  
fit <- glmmTMB::glmmTMB(log(rate) ~ 
                          log(a_t) + 
                          (1 + log(a_t) | wsd_id),
                        df_chl_filtered)

z <- coef(fit)[[1]]$wsd_id %>% 
  mutate(wsd_id = as.factor(rownames(.))) %>% 
  as_tibble() %>% 
  dplyr::select(ln_k = `(Intercept)`, # ln scaling constant
                z = `log(a_t)`, # scaling exponent
                wsd_id) %>%
  mutate(k = exp(ln_k)) %>% # scaling constant
  arrange(desc(k)) %>% 
  mutate(rank = seq_len(nrow(.)))

df_m <- df_chl_filtered %>% 
  left_join(z, by = "wsd_id") %>% 
  mutate(scl_rate = rate / (a_t)^z)


## link length ####

rdata <- list.files(path = here::here("data_fmt"),
                    full.names = TRUE) %>% 
  as_tibble() %>% 
  filter(str_detect(.$value, "df_length")) %>% 
  pull()

df_link <- foreach(i = seq_len(length(rdata)),
                   .combine = bind_rows) %do% {
                  
                     load(file = rdata[i])
                     df_length <- mutate(df_length,
                                         a_t = str_extract(rdata[i], pattern = "\\d{1,}"),
                                         a_t = as.numeric(a_t))
                     
                     return(df_length)
                   }


df_freq <- df_link %>% 
  mutate(wsd_id = factor(wsd_id)) %>% 
  filter(wsd_id %in% unique(z$wsd_id)) %>%
  group_by(wsd_id, a_t) %>% 
  summarize(class = as.numeric(names(table(ceiling(length)))),
            freq = c(table(ceiling(length))),
            cum_freq = cumsum(freq),
            prop = freq / n(),
            cum_prop = cum_freq / n()) %>% 
  ungroup() %>% 
  left_join(df_m, by = c("wsd_id", "a_t")) %>% 
  mutate(cum_prob = pexp(class, rate = rate))
  

# figure ------------------------------------------------------------------

g1 <- df_freq %>% 
  ggplot(aes(x = class,
             y = cum_prop,
             color = factor(a_t))) +
  geom_point(alpha = 0.5,
             size = 0.2) + 
  geom_line(aes(x = class,
                y = cum_prob),
            alpha = 0.5) + 
  facet_wrap(facets = ~rank,
             ncol = 6,
             scales = "free") +
  labs(y = "Cumulative proportion or probability",
       x = "Length distance class (km)",
       color = expression(A[T]~"("*km^2*")")) +
  theme_bw()

g2 <- df_m %>% 
  ggplot(aes(x = a_t,
             y = scl_rate,
             color = k,
             alpha = 0.05)) +
  geom_point() + 
  geom_line() +
  geom_hline(aes(yintercept = k,
                 color = k)) +
  facet_wrap(facets = ~rank, ncol = 6) +
  labs(x = expression(A[T]~"("*km^{-2}*")"),
       y = expression("Rescaled rate parameter ("*lambda*"')")) +
  scale_color_viridis_c() +
  guides(alpha = "none") +
  theme_bw()
