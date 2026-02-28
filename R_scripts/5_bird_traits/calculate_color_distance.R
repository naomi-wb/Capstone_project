## This script calculates a color distance
## for every species
## and saves the color data out as cleaned data
## ready for processing
## it also calculates a 'brightness' measure for every color
library(readr)
library(dplyr)
library(tidyr)
library(ggcorrplot)

color <- read_csv("Data/color data/Dale_2015_data_allpoints.csv")

# getting a brown from a satellite image of a scrubland in mexico
# a<-earthtones::get_earthtones(latitude = 24.425069,longitude = -105.645420, zoom = 9, number_of_colors = 5)
# I changed this to a more standard 'brown' color - I think 'chocolate' technically
brown <- data.frame(R = 102, B = 68, G = 0)
brown.lab <- data.frame(convertColor(brown, from = "sRGB", "Lab"))
color.lab <-
  convertColor(data.frame(
    r = color$red,
    g = color$green,
    b = color$blue
  ),
  from = "sRGB",
  "Lab")
color.complete <- cbind(color, color.lab)

if ("a.x" %in% names(color.complete)) {
  rename(color.complete, a = a.x)
}

color.complete %>%
  mutate(brown_L = brown.lab$L) %>%
  mutate(brown_a = brown.lab$a) %>%
  mutate(brown_b = brown.lab$b) %>%
  mutate(distance_from_brown = sqrt(((L - brown_L) ^ 2) + ((a.x - brown_a) ^ 2) +
                                 ((b - brown_b) ^ 2))) -> color_dists

color_dists %>%
  group_by(TipLabel) %>%
  summarise(max_distance = max(distance_from_brown),
            med_distance = median(distance_from_brown),
            mean_distance = mean(distance_from_brown)) %>%
  distinct() %>%
  left_join(., color_dists %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, distance_from_brown) %>%
              rename(lower_breast_distance=distance_from_brown) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(lower_breast_distance=mean(lower_breast_distance)), by="TipLabel") %>%
  left_join(., color_dists %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, distance_from_brown) %>%
              rename(upper_breast_distance=distance_from_brown) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(upper_breast_distance=mean(upper_breast_distance)), by="TipLabel") %>%
  ungroup() %>%
  mutate(mean_breast_distance=(lower_breast_distance+upper_breast_distance)/2) -> both_sexes_combined

color_dists %>%
  dplyr::filter(gender==1) %>%
  group_by(TipLabel) %>%
  summarise(m.max_distance = max(distance_from_brown),
            m.med_distance = median(distance_from_brown),
            m.mean_distance = mean(distance_from_brown)) %>%
  distinct() %>%
  left_join(., color_dists %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, distance_from_brown) %>%
              rename(lower_breast_distance=distance_from_brown) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(m.lower_breast_distance=mean(lower_breast_distance)), by="TipLabel") %>%
  left_join(., color_dists %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, distance_from_brown) %>%
              rename(upper_breast_distance=distance_from_brown) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(m.upper_breast_distance=mean(upper_breast_distance)), by="TipLabel") %>%
  ungroup() %>%
  mutate(m.mean_breast_distance=(m.lower_breast_distance+m.upper_breast_distance)/2) -> males

color_dists %>%
  dplyr::filter(gender==0) %>%
  group_by(TipLabel) %>%
  summarise(f.max_distance = max(distance_from_brown),
            f.med_distance = median(distance_from_brown),
            f.mean_distance = mean(distance_from_brown)) %>%
  distinct() %>%
  left_join(., color_dists %>%
              dplyr::filter(gender==0) %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, distance_from_brown) %>%
              rename(lower_breast_distance=distance_from_brown) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(f.lower_breast_distance=mean(lower_breast_distance)), by="TipLabel") %>%
  left_join(., color_dists %>%
              dplyr::filter(gender==0) %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, distance_from_brown) %>%
              rename(upper_breast_distance=distance_from_brown) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(f.upper_breast_distance=mean(upper_breast_distance)), by="TipLabel") %>%
  ungroup() %>%
  mutate(f.mean_breast_distance=(f.lower_breast_distance+f.upper_breast_distance)/2) -> females

color_differences <- both_sexes_combined %>%
  left_join(., males, by="TipLabel") %>%
  left_join(., females, by="TipLabel")

color_differences %>%
  dplyr::select(-TipLabel) %>%
  cor(.) %>%
  ggcorrplot(., lab=TRUE)

# A lot of really strong correlation among many variables describing the distance from brown
# for males, females, and regardless of sex
# the max distance (regardless of sex) appears to be one of the better predictors of intercept/slopes from 
# detectability modelling




############ distance from brown doesn't seem to be working well
############ when correlating that with the detectability of a species
############ so will try a different approach, using 'brightness'
############ https://stackoverflow.com/questions/596216/formula-to-determine-brightness-of-rgb-color
color_brightness <- color %>%
  mutate(brightness_1 = 0.2126*red + 0.7152*green + 0.0722*blue) %>%
  mutate(brightness_2 = 0.299*red + 0.587*green + 0.114*blue) %>%
  mutate(brightness_3 = 0.299*(red^2)+0.587*(green^2)+0.114*(blue^2))

color_brightness %>%
  group_by(TipLabel) %>%
  summarise(max_brightness_1 = max(brightness_1),
            med_brightness_1 = median(brightness_1),
            mean_brightness_1 = mean(brightness_1)) %>%
  distinct() %>%
  left_join(., color_brightness %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, brightness_1) %>%
              rename(lower_breast_brightness_1=brightness_1) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(lower_breast_brightness_1=mean(lower_breast_brightness_1)), by="TipLabel") %>%
  left_join(., color_brightness %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, brightness_1) %>%
              rename(upper_breast_brightness_1=brightness_1) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(upper_breast_brightness_1=mean(upper_breast_brightness_1)), by="TipLabel") %>%
  ungroup() %>%
  mutate(mean_breast_brightness_1=(lower_breast_brightness_1+upper_breast_brightness_1)/2) -> brightness_1_both_sexes_combined

color_brightness %>%
  dplyr::filter(gender==1) %>%
  group_by(TipLabel) %>%
  summarise(m.max_brightness_1 = max(brightness_1),
            m.med_brightness_1 = median(brightness_1),
            m.mean_brightness_1 = mean(brightness_1)) %>%
  distinct() %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, brightness_1) %>%
              rename(lower_breast_brightness_1=brightness_1) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(m.lower_breast_brightness_1=mean(lower_breast_brightness_1)), by="TipLabel") %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, brightness_1) %>%
              rename(upper_breast_brightness_1=brightness_1) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(m.upper_breast_brightness_1=mean(upper_breast_brightness_1)), by="TipLabel") %>%
  ungroup() %>%
  mutate(m.mean_breast_brightness_1=(m.lower_breast_brightness_1+m.upper_breast_brightness_1)/2) -> brightness_1_males

color_brightness %>%
  dplyr::filter(gender==0) %>%
  group_by(TipLabel) %>%
  summarise(f.max_brightness_1 = max(brightness_1),
            f.med_brightness_1 = median(brightness_1),
            f.mean_brightness_1 = mean(brightness_1)) %>%
  distinct() %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, brightness_1) %>%
              rename(lower_breast_brightness_1=brightness_1) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(f.lower_breast_brightness_1=mean(lower_breast_brightness_1)), by="TipLabel") %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, brightness_1) %>%
              rename(upper_breast_brightness_1=brightness_1) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(f.upper_breast_brightness_1=mean(upper_breast_brightness_1)), by="TipLabel") %>%
  ungroup() %>%
  mutate(f.mean_breast_brightness_1=(f.lower_breast_brightness_1+f.upper_breast_brightness_1)/2) -> brightness_1_females


color_brightness_1 <- brightness_1_both_sexes_combined %>%
  left_join(., brightness_1_males, by="TipLabel") %>%
  left_join(., brightness_1_females, by="TipLabel")

# repeat the above but for 'brightness 2'
color_brightness %>%
  group_by(TipLabel) %>%
  summarise(max_brightness_2 = max(brightness_2),
            med_brightness_2 = median(brightness_2),
            mean_brightness_2 = mean(brightness_2)) %>%
  distinct() %>%
  left_join(., color_brightness %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, brightness_2) %>%
              rename(lower_breast_brightness_2=brightness_2) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(lower_breast_brightness_2=mean(lower_breast_brightness_2)), by="TipLabel") %>%
  left_join(., color_brightness %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, brightness_2) %>%
              rename(upper_breast_brightness_2=brightness_2) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(upper_breast_brightness_2=mean(upper_breast_brightness_2)), by="TipLabel") %>%
  ungroup() %>%
  mutate(mean_breast_brightness_2=(lower_breast_brightness_2+upper_breast_brightness_2)/2) -> brightness_2_both_sexes_combined

color_brightness %>%
  dplyr::filter(gender==1) %>%
  group_by(TipLabel) %>%
  summarise(m.max_brightness_2 = max(brightness_2),
            m.med_brightness_2 = median(brightness_2),
            m.mean_brightness_2 = mean(brightness_2)) %>%
  distinct() %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, brightness_2) %>%
              rename(lower_breast_brightness_2=brightness_2) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(m.lower_breast_brightness_2=mean(lower_breast_brightness_2)), by="TipLabel") %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, brightness_2) %>%
              rename(upper_breast_brightness_2=brightness_2) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(m.upper_breast_brightness_2=mean(upper_breast_brightness_2)), by="TipLabel") %>%
  ungroup() %>%
  mutate(m.mean_breast_brightness_2=(m.lower_breast_brightness_2+m.upper_breast_brightness_2)/2) -> brightness_2_males

color_brightness %>%
  dplyr::filter(gender==0) %>%
  group_by(TipLabel) %>%
  summarise(f.max_brightness_2 = max(brightness_2),
            f.med_brightness_2 = median(brightness_2),
            f.mean_brightness_2 = mean(brightness_2)) %>%
  distinct() %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, brightness_2) %>%
              rename(lower_breast_brightness_2=brightness_2) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(f.lower_breast_brightness_2=mean(lower_breast_brightness_2)), by="TipLabel") %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, brightness_2) %>%
              rename(upper_breast_brightness_2=brightness_2) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(f.upper_breast_brightness_2=mean(upper_breast_brightness_2)), by="TipLabel") %>%
  ungroup() %>%
  mutate(f.mean_breast_brightness_2=(f.lower_breast_brightness_2+f.upper_breast_brightness_2)/2) -> brightness_2_females


color_brightness_2 <- brightness_2_both_sexes_combined %>%
  left_join(., brightness_2_males, by="TipLabel") %>%
  left_join(., brightness_2_females, by="TipLabel")

# repeat the above but for 'brightness 3'
color_brightness %>%
  group_by(TipLabel) %>%
  summarise(max_brightness_3 = max(brightness_3),
            med_brightness_3 = median(brightness_3),
            mean_brightness_3 = mean(brightness_3)) %>%
  distinct() %>%
  left_join(., color_brightness %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, brightness_3) %>%
              rename(lower_breast_brightness_3=brightness_3) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(lower_breast_brightness_3=mean(lower_breast_brightness_3)), by="TipLabel") %>%
  left_join(., color_brightness %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, brightness_3) %>%
              rename(upper_breast_brightness_3=brightness_3) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(upper_breast_brightness_3=mean(upper_breast_brightness_3)), by="TipLabel") %>%
  ungroup() %>%
  mutate(mean_breast_brightness_3=(lower_breast_brightness_3+upper_breast_brightness_3)/2) -> brightness_3_both_sexes_combined

color_brightness %>%
  dplyr::filter(gender==1) %>%
  group_by(TipLabel) %>%
  summarise(m.max_brightness_3 = max(brightness_3),
            m.med_brightness_3 = median(brightness_3),
            m.mean_brightness_3 = mean(brightness_3)) %>%
  distinct() %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, brightness_3) %>%
              rename(lower_breast_brightness_3=brightness_3) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(m.lower_breast_brightness_3=mean(lower_breast_brightness_3)), by="TipLabel") %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, brightness_3) %>%
              rename(upper_breast_brightness_3=brightness_3) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(m.upper_breast_brightness_3=mean(upper_breast_brightness_3)), by="TipLabel") %>%
  ungroup() %>%
  mutate(m.mean_breast_brightness_3=(m.lower_breast_brightness_3+m.upper_breast_brightness_3)/2) -> brightness_3_males

color_brightness %>%
  dplyr::filter(gender==0) %>%
  group_by(TipLabel) %>%
  summarise(f.max_brightness_3 = max(brightness_3),
            f.med_brightness_3 = median(brightness_3),
            f.mean_brightness_3 = mean(brightness_3)) %>%
  distinct() %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="lower_breast") %>%
              dplyr::select(TipLabel, brightness_3) %>%
              rename(lower_breast_brightness_3=brightness_3) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(f.lower_breast_brightness_3=mean(lower_breast_brightness_3)), by="TipLabel") %>%
  left_join(., color_brightness %>%
              dplyr::filter(gender==1) %>%
              dplyr::filter(patch=="upper_breast") %>%
              dplyr::select(TipLabel, brightness_3) %>%
              rename(upper_breast_brightness_3=brightness_3) %>%
              distinct() %>%
              group_by(TipLabel) %>%
              summarize(f.upper_breast_brightness_3=mean(upper_breast_brightness_3)), by="TipLabel") %>%
  ungroup() %>%
  mutate(f.mean_breast_brightness_3=(f.lower_breast_brightness_3+f.upper_breast_brightness_3)/2) -> brightness_3_females


color_brightness_3 <- brightness_3_both_sexes_combined %>%
  left_join(., brightness_3_males, by="TipLabel") %>%
  left_join(., brightness_3_females, by="TipLabel")

# put the three brightnesses together
brightness_all <- color_brightness_1 %>%
  left_join(color_brightness_2) %>%
  left_join(color_brightness_3)

brightness_all %>%
  dplyr::select(-TipLabel) %>%
  cor(.) %>%
  ggcorrplot(.)

# all pretty correlated, among the three measures of 'brightness'
# so will just go with 'brightness_1' for now
brightness_summary <- color_brightness %>% 
  group_by(TipLabel, gender) %>% 
  summarize(total_brightness=sum(brightness_1)) %>% 
  pivot_wider(names_from=gender, values_from=total_brightness) %>% 
  rename(female_total_brightness=`0`) %>% 
  rename(male_total_brightness=`1`) %>% 
  mutate(median_total_brightness=(female_total_brightness+male_total_brightness)/2) %>%
  left_join(., color_brightness %>%
              group_by(TipLabel) %>%
              summarize(max_brightness=max(brightness_1)), by="TipLabel")

distance_summary <- color_dists %>%
  group_by(TipLabel, gender) %>%
  summarize(total_distance=sum(distance_from_brown)) %>%
  pivot_wider(names_from=gender, values_from=total_distance) %>%
  rename(female_total_distance_from_brown=`0`) %>%
  rename(male_total_distance_from_brown=`1`) %>%
  mutate(median_total_distance_from_brown=(female_total_distance_from_brown+male_total_distance_from_brown)/2) %>%
  left_join(., color_dists %>%
              group_by(TipLabel) %>%
              summarize(max_distance=max(distance_from_brown)), by="TipLabel")

color_summary <- brightness_summary %>%
  left_join(distance_summary)

# plot correlations
color_summary %>%
  ungroup() %>%
  dplyr::select(-TipLabel) %>%
  cor(.) %>%
  ggcorrplot(., lab=TRUE)

# lots of correlation among these variables
# so let's just pick the median and remove male and female
# and then can remove other variables later on down the track
color_summary.2 <- color_summary %>%
  dplyr::select(TipLabel, max_distance, max_brightness, 
                median_total_distance_from_brown, median_total_brightness)

# plot again
color_summary.2 %>%
  ungroup() %>%
  dplyr::select(-TipLabel) %>%
  cor(.) %>%
  ggcorrplot(., lab=TRUE)

# will write these our for now
# write out the summary file
write_csv(color_summary.2,
          "Data/color data/color_data_by_species.csv")
