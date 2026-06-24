#### SETUP ---------------------------------------------------------------------
# setup age classes, fertility, mortality
age_classes = c(seq(4, 180, 4),seq(240,1200,by=60))/12
bin_width = diff(c(0, age_classes))

#### DEFINE WAIFW MATRICES TO TEST ---------------------------------------------
## series from flat to peaky to diagonal to polymod
waifw = get_waifws(age_classes, background = 0.001, rescale = TRUE)

# plot WAIFW matrices
waifw_labs = c("flat", "peak age 5", "peak age 10", "diagonal", "POLYMOD")
names(waifw_labs) = 1:5

lapply(waifw, melt) %>%
  bind_rows(.id = "waifw_id") %>%
  mutate(value_scaled = value, .by = c("waifw_id")) %>%
  ggplot(aes(x = Var1, y = Var2, fill = value_scaled)) + 
  geom_tile() + 
  facet_wrap(vars(waifw_id), labeller = labeller(waifw_id = waifw_labs)) +
  scale_fill_viridis_c() + 
  scale_x_continuous(expand = c(0,0),
                     # breaks = which(age_classes_jcm %in% c(2, 4, 6, 8, 10, 30, 50, 70)),
                     # labels = c(2, 4, 6, 8, 10, 30, 50, 70),
                     name = "age (years)") +
  scale_y_continuous(expand = c(0,0),
                     # breaks = which(age_classes_jcm %in% c(2, 4, 6, 8, 10, 30, 50, 70)),
                     # labels = c(2, 4, 6, 8, 10, 30, 50, 70),
                     name = "age (years)")

#### GET FULL WAIFW ------------------------------------------------------------
age_classes_full = c(seq(4,900,by=4))/12
waifw_full = get_waifws(age_classes_full, background = 0.001, rescale = TRUE)

