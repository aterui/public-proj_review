
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

df_chl <- df_chl %>% 
  group_by(wsd_id) %>%
  filter(n() == n_distinct(.$a_t)) %>% 
  ungroup() %>% 
  mutate(wsd_id = factor(wsd_id)) 

## power law scaling between rate and A_T  
fit <- glmmTMB::glmmTMB(log(rate) ~ 
                          log(a_t) + 
                          (1 + log(a_t) | wsd_id),
                        df_chl)

z <- coef(fit)[[1]]$wsd_id %>% 
  mutate(wsd_id = as.factor(rownames(.))) %>% 
  as_tibble() %>% 
  dplyr::select(ln_k = `(Intercept)`, # ln scaling constant
                z = `log(a_t)`, # scaling exponent
                wsd_id) %>%
  mutate(k = exp(ln_k)) %>% # scaling constant
  arrange(desc(k)) %>% 
  mutate(rank = seq_len(nrow(.)))

df_m <- df_chl %>% 
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
                     df_length <- df_length %>% 
                       mutate(wsd_id = factor(wsd_id),
                              a_t = str_extract(rdata[i], pattern = "\\d{1,}"),
                              a_t = as.numeric(a_t)) %>% 
                       select(-geometry)
                     
                     return(df_length)
                   }


df_freq <- df_link %>% 
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
  

# Horton's ratio ----------------------------------------------------------

df_horton <- df_link %>%
  group_by(wsd_id) %>% 
  filter(a_t == 1) %>% 
  group_by(wsd_id, tributary, order) %>% 
  summarize(link = sum(length)) %>% 
  group_by(wsd_id, order) %>% 
  summarize(ave_length = mean(link),
            n_link = n())

# length ratio
fit_rl <- glmmTMB::glmmTMB(log(ave_length) ~ order + (1 + order | wsd_id),
                           df_horton)

rl <- coef(fit_rl)[[1]]$wsd_id %>% 
  mutate(wsd_id = as.factor(rownames(.))) %>% 
  as_tibble() %>% 
  dplyr::select(ln_rl = `order`, # scaling exponent
                wsd_id) %>%
  mutate(rl = exp(ln_rl))

# branching ratio
fit_rb <- glmmTMB::glmmTMB(log(n_link) ~ order + (1 + order | wsd_id),
                           df_horton)

rb <- coef(fit_rb)[[1]]$wsd_id %>% 
  mutate(wsd_id = as.factor(rownames(.))) %>% 
  as_tibble() %>% 
  dplyr::select(ln_rb = `order`, # scaling exponent
                wsd_id) %>%
  mutate(rb = exp(-ln_rb))

df_r <- df_m %>% 
  left_join(rl, by = "wsd_id") %>% 
  left_join(rb, by = "wsd_id") %>% 
  filter(a_t == 1) %>% 
  mutate(beta = 1 + ln_rl / ln_rb)

# figure ------------------------------------------------------------------

g1 <- df_freq %>% 
  filter(area >= units::set_units(500, "km^2")) %>% 
  mutate(rank = n_distinct(.$wsd_id) + 1 - as.numeric(factor(k))) %>% 
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
  filter(area >= units::set_units(500, "km^2")) %>% 
  mutate(rank = n_distinct(.$wsd_id) + 1 - as.numeric(factor(k))) %>% 
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

g3 <- df_r %>% 
  ggplot(aes(y = rb,
             x = k,
             alpha = 0.1)) +
  geom_point() +
  geom_smooth(method = "lm", color = grey(0.1)) +
  guides(alpha = "none") +
  labs(y = expression(R[B]),
       x = expression("Scaling constant k (="~lambda*"')")) +
  theme_bw()
