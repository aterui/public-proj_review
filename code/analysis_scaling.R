
# setup -------------------------------------------------------------------

rm(list = ls())
pacman::p_load(tidyverse)

# analysis ----------------------------------------------------------------

## in the following analysis, watersheds > 300km^2 were selected
## small watersheds could have an issue of estimating the rate parameter
## i.e., the number of links is limited

## load river network data
load(file = here::here("data_fmt/df_chl.RData"))

df_chl_filtered <- df_chl %>% 
  filter(area > units::set_units(300, "km^2")) %>% # select watersheds > 300 km^2
  group_by(wsd_id) %>%
  filter(n() == n_distinct(.$a_t)) %>% 
  ungroup() %>% 
  mutate(wsd_id = factor(wsd_id)) 

## power law scaling between rate and A_T  
fit <- lme4::lmer(log(rate) ~ 
                    log(a_t) + 
                    (1 + log(a_t) | wsd_id),
                  df_chl_filtered)

z <- coef(fit)$wsd_id %>% 
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


# figure ------------------------------------------------------------------

g <- df_m %>% 
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

print(g)
