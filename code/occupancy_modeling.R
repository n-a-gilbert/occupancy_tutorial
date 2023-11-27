library(tidyverse)
library(unmarked)
library(spOccupancy)
library(sf)
library(janitor)
library(here)

setwd(here::here("data"))

final <- readr::read_csv("ybch_data_v01.csv") |> 
  dplyr::mutate(route = as.numeric(factor(route)))

dplyr::glimpse(final)

#--------------------------------------------------------------
# simplest implementation: unmarked
#--------------------------------------------------------------

# formatting data for going into models is the hardest part!

# y is the detection/nondetection history
# 1 if Yellow-breasted Chat was detected, 0 if not
# 173 rows becacuse 173 sites were visited
# 3 columns because there were 3 repeated visits
y <- final |> 
  dplyr::select(y_1, y_2, y_3) |> 
  base::as.matrix() |> 
  base::unname()

# site covariates
# % canopy cover within 600 m of point
# route ID to be used as a random effect
siteCovs <- final |> 
  dplyr::select(canopy, route) |> 
  base::as.data.frame()

# detection covariates
# variables unique for each site-visit combo
# time (measured as minutes after sunrise)
# number of cars passing during survey
# wind speed, 
# temperature
obsCovs <- list(
  date = final |> dplyr::select(ordinal_date_1:ordinal_date_3) |> base::as.matrix() |> base::unname(),
  time = final |> dplyr::select(min_sun_1:min_sun_3) |> base::as.matrix() |> base::unname(),
  cars = final |> dplyr::select(cars_1:cars_3) |> base::as.matrix() |> base::unname(),
  wind = final |> dplyr::select(wind_1:wind_3) |> base::as.matrix() |> base::unname(), 
  temp = final |> dplyr::select(temp_1:temp_3) |> base::as.matrix() |> base::unname())

# now you have to package all that into a special "unmarked frame" for unmarked
umf <- unmarked::unmarkedFrameOccu( y = y, 
                                    siteCovs = siteCovs, 
                                    obsCovs = obsCovs )

# simplest occupancy model
# implemented in unmarked, maximum likelihood
# date, time, time^2, number of cars, wind speed, and temp used to predict detection
# canopy and canopy^2 used to predict occupancy
m1 <- unmarked::occu(~date + time + I(time^2) + cars + wind + temp ~ canopy + I(canopy^2), data = umf)

# strong evidence for quadratic effect of canopy cover on occupancy
# strong evidence of effect of date on detection (lower detection later in season), 
# evidence of time effect (highest around sunrise/ soon after, dropping later in the morning)
# evidence of wind effect (lower detection when windy)
summary(m1)

# marginal effect plots of the predictors
# unmarked has nice built in function for visualizing these
# depending on the version of unmarked you have installed, this may or may not work
unmarked::plotEffects(m1, type = "state", covariate = "canopy")
unmarked::plotEffects(m1, type = "det", covariate = "date")
unmarked::plotEffects(m1, type = "det", covariate = "time")
unmarked::plotEffects(m1, type = "det", covariate = "cars")
unmarked::plotEffects(m1, type = "det", covariate = "wind")
unmarked::plotEffects(m1, type = "det", covariate = "temp")

#---------------------------------------------------------------------
# spOccupancy time!
# this is a big improvement over unmarked because you can fit
# random effects or spatial models
#---------------------------------------------------------------------

# visualize the locations of surveys
# 173 point-count locations were conducted along roadsides
# 15 routes; each route had ~12 points, separated by ~1km
# potentially important to account for non-independence between points?
final |> 
  ggplot2::ggplot(aes(x = long, y = lat, color = factor(route))) + 
  ggplot2::geom_point() + 
  ggplot2::coord_fixed(1.3)

# package up data for spOccupancy
ybch <- list(
  y = y,
  occ.covs = final |> dplyr::select(canopy, route), 
  det.covs = obsCovs, 
  coords = final |> 
    dplyr::select(long, lat) |> 
    sf::st_as_sf(coords = c("long", "lat"), 
                 crs = 4326, 
                 agr = "constant") |> 
    sf::st_transform(crs = 26916) |> 
    sf::st_coordinates()
)

# model 2: put random intercept for route in the occupancy linear predictor
m2 <- spOccupancy::PGOcc(
  occ.formula = ~canopy + I(canopy^2) + (1 | route),
  det.formula = ~date + time + I(time^2) + cars + wind + temp,
  data = ybch,
  n.samples = 1000,
  n.chains = 3)

summary(m2)

# grab parameters
# unfortunately, not built in functionality to visualize marginal effects (yet)
params <- m2$beta.samples |> 
  base::cbind( tibble::as_tibble(m2$beta.star.samples)) |> 
  tibble::as_tibble() |> 
  janitor::clean_names() |> 
  tidyr::pivot_longer(route_1:route_15, names_to = "route", values_to = "int_re") 

all <- merge( params, tibble::tibble(canopy = seq(from = min(final$canopy), 
                                                  to = max(final$canopy),
                                                  by = 0.5)),
              by = NULL) |> 
  tibble::as_tibble() |> 
  dplyr::mutate(lpsi = intercept + canopy.x*canopy.y + i_canopy_2 * canopy.y * canopy.y + int_re) |> 
  dplyr::mutate(psi = plogis(lpsi)) |> 
  dplyr::group_by(route, canopy.y) |> 
  dplyr::summarise(mean = mean(psi), 
                   l95 = quantile(psi, c(0.025)), 
                   u95 = quantile(psi, c(0.975)))

all |> 
  ggplot2::ggplot(aes(x = canopy.y, y = mean, color = route)) + 
  ggplot2::geom_line(size = 2)

# spatial occupancy model 
m3 <- spPGOcc(
  occ.formula = ~canopy + I(canopy^2),
  det.formula = ~date + time + I(time^2) + cars + wind + temp,
  data = ybch,
  n.batch = 400,
  batch.length = 25,
  n.chains = 3)

summary(m3)

post <- m3$beta.samples |> 
  tibble::as_tibble() |> 
  janitor::clean_names() |> 
  dplyr::slice(1:2000)

merge(post, data.frame(can = seq(from = min(ybch$occ.covs$canopy), to = max(ybch$occ.covs$canopy), 
                                 by = 0.1)),
      by = NULL) |> 
  tibble::as_tibble() |> 
  dplyr::mutate(lpsi = intercept + canopy*can + i_canopy_2 * can * can ) |> 
  dplyr::mutate(psi = plogis(lpsi)) |> 
  dplyr::group_by(can) |> 
  dplyr::summarise(l95 = quantile(psi , c(0.025)), 
                   mean = mean(psi), 
                   u95 = quantile(psi, c(0.975))) |> 
  
  ggplot2::ggplot(aes(x = can, y = mean)) +
  ggplot2::geom_ribbon(aes(ymin = l95, ymax = u95), color = NA, alpha = 0.2) + 
  ggplot2::geom_line(size = 1.5)

alpha_post <- m3$alpha.samples |> 
  tibble::as_tibble() |> 
  janitor::clean_names() |> 
  dplyr::slice(1:2000)

merge(alpha_post, 
      data.frame( date = seq(from = min(ybch$det.covs$date), 
                             to = max(ybch$det.covs$date),
                             by = 0.25)),
      by = NULL) |> 
  tibble::as_tibble() |> 
  dplyr::mutate(lp = intercept + date.x * date.y  ) |> 
  dplyr::mutate(p = plogis(lp)) |> 
  dplyr::group_by(date.y) |> 
  dplyr::summarise(l95 = quantile(p , c(0.025)), 
                   mean = mean(p), 
                   u95 = quantile(p, c(0.975))) |> 
  
  ggplot2::ggplot(aes(x = date.y, y = mean)) +
  ggplot2::geom_ribbon(aes(ymin = l95, ymax = u95), color = NA, alpha = 0.2) + 
  ggplot2::geom_line(size = 1.5) +
  ggplot2::labs(x = "Date", 
                y = "Predicted detection probability")

merge(alpha_post, 
      data.frame( time = seq(from = min(ybch$det.covs$time), 
                             to = max(ybch$det.covs$time),
                             by = 0.25)),
      by = NULL) |> 
  tibble::as_tibble() |> 
  dplyr::mutate(lp = intercept + time.x * time.y + i_time_2 * time.y * time.y  ) |> 
  dplyr::mutate(p = plogis(lp)) |> 
  dplyr::group_by(time.y) |> 
  dplyr::summarise(l95 = quantile(p , c(0.025)), 
                   mean = mean(p), 
                   u95 = quantile(p, c(0.975))) |> 
  
  ggplot2::ggplot(aes(x = time.y, y = mean)) +
  ggplot2::geom_ribbon(aes(ymin = l95, ymax = u95), color = NA, alpha = 0.2) + 
  ggplot2::geom_line(size = 1.5) +
  ggplot2::labs(x = "Minutes after sunrise", 
                y = "Predicted detection probability")

merge(alpha_post, 
      data.frame( cars = seq(from = min(ybch$det.covs$cars), 
                             to = max(ybch$det.covs$cars),
                             by = 0.25)),
      by = NULL) |> 
  tibble::as_tibble() |> 
  dplyr::mutate(lp = intercept + cars.x * cars.y  ) |> 
  dplyr::mutate(p = plogis(lp)) |> 
  dplyr::group_by(cars.y) |> 
  dplyr::summarise(l95 = quantile(p , c(0.025)), 
                   mean = mean(p), 
                   u95 = quantile(p, c(0.975))) |> 
  
  ggplot2::ggplot(aes(x = cars.y, y = mean)) +
  ggplot2::geom_ribbon(aes(ymin = l95, ymax = u95), color = NA, alpha = 0.2) + 
  ggplot2::geom_line(size = 1.5) +
  ggplot2::labs(x = "Number of cars", 
                y = "Predicted detection probability")

merge(alpha_post, 
      data.frame( wind = seq(from = min(ybch$det.covs$wind), 
                             to = max(ybch$det.covs$wind),
                             by = 0.25)),
      by = NULL) |> 
  tibble::as_tibble() |> 
  dplyr::mutate(lp = intercept + wind.x * wind.y  ) |> 
  dplyr::mutate(p = plogis(lp)) |> 
  dplyr::group_by(wind.y) |> 
  dplyr::summarise(l95 = quantile(p , c(0.025)), 
                   mean = mean(p), 
                   u95 = quantile(p, c(0.975))) |> 
  
  ggplot2::ggplot(aes(x = wind.y, y = mean)) +
  ggplot2::geom_ribbon(aes(ymin = l95, ymax = u95), color = NA, alpha = 0.2) + 
  ggplot2::geom_line(size = 1.5) +
  ggplot2::labs(x = "Wind speed", 
                y = "Predicted detection probability")

merge(alpha_post, 
      data.frame( temp = seq(from = min(ybch$det.covs$temp), 
                             to = max(ybch$det.covs$temp),
                             by = 0.25)),
      by = NULL) |> 
  tibble::as_tibble() |> 
  dplyr::mutate(lp = intercept + temp.x * temp.y  ) |> 
  dplyr::mutate(p = plogis(lp)) |> 
  dplyr::group_by(temp.y) |> 
  dplyr::summarise(l95 = quantile(p , c(0.025)), 
                   mean = mean(p), 
                   u95 = quantile(p, c(0.975))) |> 
  
  ggplot2::ggplot(aes(x = temp.y, y = mean)) +
  ggplot2::geom_ribbon(aes(ymin = l95, ymax = u95), color = NA, alpha = 0.2) + 
  ggplot2::geom_line(size = 1.5) +
  ggplot2::labs(x = "Temperature", 
                y = "Predicted detection probability")