getwd()
rm(list = ls())
# load penguins
library(readxl)
sheet_path = './r4ds/data/penguins.xlsx'
tab_names = excel_sheets(sheet_path)

library(plyr)
penguins = lapply(tab_names, function(x){
    read_excel(path = sheet_path, sheet = x)}) %>% 
    rbind.fill() %>% as_tibble()

str(penguins)
# change cols types

penguins %<>% 
    mutate_at(3:6,as.integer) -> penguins

str(penguins)
library(ggplot2)
penguins %>% ggplot(
    aes(
    x = flipper_length_mm,
    y = body_mass_g)
    )+
    geom_point(
       mapping = aes(color = species)
    )+
    geom_smooth(method = "lm")

library(ggthemes)
ggplot(
    data = penguins,
    mapping = aes(x = flipper_length_mm, y = body_mass_g)
) +
    geom_point(aes(color = species, shape = species,
                   ), na.rm=TRUE) +
    geom_smooth(method = "lm") +
    labs(
        title = "Body mass and flipper length",
        subtitle = "Dimensions for Adelie, Chinstrap, and Gentoo Penguins",
        x = "Flipper length (mm)", y = "Body mass (g)",
        color = "Species", shape = "Species"
    ) +
    scale_color_colorblind()


penguins %>% ggplot(
    aes(
        x = flipper_length_mm,
        y = body_mass_g)
)+
    geom_point(
        mapping = aes(color = bill_depth_mm),
        na.rm = TRUE
    )+
    geom_smooth()

ggplot(
    data = penguins,
    mapping = aes(x = flipper_length_mm, y = body_mass_g, color = island)
) +
    geom_point() +
    geom_smooth(se = FALSE)

ggplot(
    data = penguins,
    mapping = aes(x = flipper_length_mm, y = body_mass_g)
) +
    geom_point() +
    geom_smooth()

ggplot() +
    geom_point(
        data = penguins,
        mapping = aes(x = flipper_length_mm, y = body_mass_g)
    ) +
    geom_smooth(
        data = penguins,
        mapping = aes(x = flipper_length_mm, y = body_mass_g)
    )
