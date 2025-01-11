library(tidyverse)

# Load the data
GRooT <- read_csv("GRooTFullVersion.csv") |> 
  glimpse() 
names(GRooT)
#getting fred vs groot counts

FRED <- read_csv("FRED3_Filtered_Database_20240222-211455_KRC_04_08.csv") |>
  glimpse()

names(FRED)

FRED_counts <- FRED |>
  count(`Belowground part`)
print(FRED_counts)

GRooT_counts <- GRooT |>
  count(`belowgroundEntities`)
print(GRooT_counts)


#Filter groot data to only include root data on fine roots

GRooT_filtered <- GRooT |>
  filter(belowgroundEntities == "FR") |>
  select(decimalLongitud, decimalLatitude, genusTNRS, speciesTNRS, 
         belowgroundEntities, traitName, traitValue, 
         mycorrhizalAssociationTypeFungalRoot, growthForm)  # Adding mycorrhizal type and growthForm columns

names(GRooT_filtered)

# Remove points with no lat and long
GRooT_parsed <- GRooT_filtered |>
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitud) & !is.na(traitValue)) # 37,580 obs-----if just filter by FR: 53,543 obs

# Filter to Bergmann variables: N, D, SRL, RTD
GRooT_bergmann <- GRooT_parsed |>
  filter(traitName %in% c("Root_N_concentration", "Mean_Root_diameter", "Root_tissue_density", "Specific_root_length")) # 19,796 obs-----if just filter by FR: 27,672 obs

GRooT_bergmann <- GRooT_bergmann |>
  filter(!is.na(traitValue))

#Load and glimpse the soilgrids x GRooT unique points data
GROOT_SG <- read_csv("GRooT_uniqueplots_soilgrids.csv") |> glimpse()

#Select relevant columns from GROOT_SG to keep for joining to groot
GROOT_SG <- GROOT_SG |>
  select(soc, bdod, lat, lon) |> glimpse()

#Rename columns to match for the join
GRooT_parsed <- GRooT_parsed |>
  rename(
    lon = decimalLongitud,
    lat = decimalLatitude
  )

# Perform the left join. Here we are joining by latitude and longitude (location)
GRooT_and_SG_joined <- left_join(GRooT_parsed, GROOT_SG, by = c("lat", "lon")) |> glimpse()


#List unique root traits frmo joined data
unique(GRooT_and_SG_joined$traitName)



#Filter to Bergmann variables again if needed
GRooT_and_SG_joined_BergmannVar <- GRooT_and_SG_joined |>
  filter(traitName %in% c("Root_N_concentration", "Mean_Root_diameter", "Root_tissue_density", "Specific_root_length", "soc")) # 27,672 obs

#Summarize by taking the first value
GRooT_and_SG_joined_wide <- GRooT_and_SG_joined_BergmannVar |>
  pivot_wider(names_from = traitName, values_from = traitValue, values_fn = function(x) x[1])

#Print reshaped data
print(GRooT_and_SG_joined_wide)


#Combine genus and species columns into a new column
GRooT_and_SG_joined_wide <- GRooT_and_SG_joined_wide |>
  mutate(genus_species = paste(`genusTNRS`, `speciesTNRS`)) |> glimpse()



names(GRooT_and_SG_joined_wide)


library(plyr)
spp_GRooT_and_SG_joined_wide <- ddply(GRooT_and_SG_joined_wide , c("genus_species", "mycorrhizalAssociationTypeFungalRoot", "growthForm"), summarize,
                                      Root_N_avg = mean(`Root_N_concentration`, na.rm = TRUE),
                                      Root_N_count = sum(!is.na(`Root_N_concentration`)),
                                      Root_N_sd = sd(`Root_N_concentration`, na.rm = TRUE),
                                      Root_N_se = Root_N_sd / sqrt(Root_N_count),
                                      Root_D_avg = mean(`Mean_Root_diameter`, na.rm = TRUE),
                                      Root_D_count = sum(!is.na(`Mean_Root_diameter`)),
                                      Root_D_sd = sd(`Mean_Root_diameter`, na.rm = TRUE),
                                      Root_D_se = Root_D_sd / sqrt(Root_D_count),
                                      Root_RTD_avg = mean(`Root_tissue_density`, na.rm = TRUE),
                                      Root_RTD_count = sum(!is.na(`Root_tissue_density`)),
                                      Root_RTD_sd = sd(`Root_tissue_density`, na.rm = TRUE),
                                      Root_RTD_se = Root_RTD_sd / sqrt(Root_RTD_count),
                                      SRL_avg = mean(`Specific_root_length`, na.rm = TRUE),
                                      SRL_count = sum(!is.na(`Specific_root_length`)),
                                      SRL_sd = sd(`Specific_root_length`, na.rm = TRUE),
                                      SRL_se = SRL_sd / sqrt(SRL_count),
                                      soc_avg = mean(soc, na.rm = TRUE),
                                      soc_count = sum(!is.na(soc)),
                                      soc_sd = sd(soc, na.rm = TRUE),
                                      soc_se = soc_sd / sqrt(soc_count))
names(spp_GRooT_and_SG_joined_wide)


#now creating a subset df with only trait and soc avg columns
spp_GRooT_and_SG_joined_wide_avgs <- spp_GRooT_and_SG_joined_wide |>
  select(genus_species,
         Root_N_avg,
         Root_D_avg, 
         Root_RTD_avg,
         SRL_avg,
         soc_avg)

#now we want to remove NAs
spp_GRooT_and_SG_joined_wide_avgs <- na.omit(spp_GRooT_and_SG_joined_wide_avgs) |>
  glimpse()       #375 observations of 6 variables

####################################################################################################

#here we are doing a break-point analysis (rank-based transformation) to determine mineral versus organic cutoff
#unload plyr bc conflicts with dplyr
detach("package:plyr", unload = TRUE)


#Load libraries
library(segmented)

# Sort the DataFrame by SOC values and add a rank column
sorted <- spp_GRooT_and_SG_joined_wide |>
  arrange(desc(soc_avg)) |>
  mutate(rank = row_number())


# Clean the data by removing rows with missing values
clean <- sorted |>
  filter(!is.na(soc_avg), !is.na(rank))


# Fit the initial linear model for rank ~ soc_avg
lm_model <- lm(rank ~ soc_avg, data = clean)

# Fit the segmented model (using SOC as the predictor, not rank)
segmented_model <- segmented(lm_model, seg.Z = ~ rank, npsi = 3)

# Get the estimated breakpoints
breakpoints <- segmented_model$psi[, 2]
print(breakpoints)  # Print the estimated breakpoints
#psi1.rank psi2.rank psi3.rank 
#55.0500  207.7386  975.5443

# Predict using the segmented model
predicted_values <- predict(segmented_model)

# Plot the segmented regression with the original data and multiple breakpoints
ggplot(clean, aes(x = soc_avg, y = rank)) +
  geom_point() +
  geom_line(aes(y = predicted_values), color = "blue", size = 0.8) + # Use predicted values correctly
  
  # Add vertical lines for each breakpoint
  geom_vline(xintercept = breakpoints, linetype = "dashed", color = "red", size = 1.0) +
  
  # Add error bars for SOC
  geom_errorbarh(aes(xmin = soc_avg - soc_sd, xmax = soc_avg + soc_sd), height = 0.2) +
  
  # Add labels for breakpoints
  annotate("text", x = breakpoints[1], y = 1250, label = paste("Breakpoint 1: ", round(breakpoints[1], 2), " mg/g SOC"), hjust = -0.05, color = "black") +
  annotate("text", x = breakpoints[2], y = 1150, label = paste("Breakpoint 2: ", round(breakpoints[2], 2), " mg/g SOC"), hjust = -0.05, color = "black") +
  annotate("text", x = breakpoints[3], y = 1050, label = paste("Breakpoint 3: ", round(breakpoints[3], 2), " mg/g SOC"), hjust = -0.05, color = "black") +
  
  # Add plot labels and theme
  labs(
    title = "Breakpoint Analysis of avg SOC",
    x = "avg SOC",
    y = "Rank"
  ) +
  theme_minimal()


####################################################################################################

#now we want to create histograms for each of the top 4 FR traits, then one for all of them pooled together
#i also created a custom theme to implement for each for clarity and consistency

custom_theme <- theme(
  panel.background = element_rect(fill = "white"),
  panel.grid.major = element_blank(),  # Remove major grid lines
  panel.grid.minor = element_blank(),  # Remove minor grid lines
  axis.line = element_line(color = "black")  # Add black axis lines
)

names(spp_GRooT_and_SG_joined_wide_avgs)


#for N soc histogram
# Ensure that dplyr is loaded
library(dplyr)

# Select the columns and remove NA values
# Explicitly call dplyr's select function
spp_GRooT_and_SG_joined_wide_N <- dplyr::select(spp_GRooT_and_SG_joined_wide_avgs, Root_N_avg, soc_avg) |>
  na.omit()

# Check the data
glimpse(spp_GRooT_and_SG_joined_wide_N)


# ggplot histograms with custom theme and specified color scheme
rootN_soc_hist_gg <- ggplot(data = spp_GRooT_and_SG_joined_wide_N, aes(x = soc_avg)) +
  geom_histogram(color = "black", fill = "#51127c") +
  labs(title = "Root N", x = " ", y = " ") +
  custom_theme

print(rootN_soc_hist_gg)

#for RTD soc histogram
spp_GRooT_and_SG_joined_wide_D <- dplyr::select(spp_GRooT_and_SG_joined_wide_avgs, Root_D_avg, soc_avg) |>
  na.omit()
  glimpse(spp_GRooT_and_SG_joined_wide_D)


rootD_soc_hist_gg <- ggplot(data = spp_GRooT_and_SG_joined_wide_D, aes(x = soc_avg)) +
  geom_histogram(color = "black", fill = "#721f81") +
  labs(title = "Root D", x = " ", y = " ") +
  custom_theme
print(rootD_soc_hist_gg)


#for SRL soc histogram
spp_GRooT_and_SG_joined_wide_SRL <- dplyr::select(spp_GRooT_and_SG_joined_wide_avgs, SRL_avg , soc_avg) |>
  na.omit() |>
  glimpse()

rootSRL_soc_hist_gg <- ggplot(data = spp_GRooT_and_SG_joined_wide_SRL, aes(x = soc_avg)) +
  geom_histogram(color = "black", fill = "#f1605d") +
  labs(title = "SRL", x = " ", y = " ") +
  custom_theme
print(rootSRL_soc_hist_gg)


#for RTD soc histogram
spp_GRooT_and_SG_joined_wide_RTD <- dplyr::select(spp_GRooT_and_SG_joined_wide_avgs, Root_RTD_avg , soc_avg) |>
  na.omit() |>
  glimpse()

rootRTD_soc_hist_gg <- ggplot(data = spp_GRooT_and_SG_joined_wide_RTD, aes(x = soc_avg)) +
  geom_histogram(color = "black", fill = "#b73779") +
  labs(title = "RTD", x = " ", y = " ") +
  custom_theme
print(rootRTD_soc_hist_gg)

library(patchwork)
library(grid)

# Combine the plots using the patchwork package in a 2x2 grid
plots_combined <- (rootN_soc_hist_gg + rootD_soc_hist_gg) / 
  (rootRTD_soc_hist_gg + rootSRL_soc_hist_gg) +
  plot_layout(guides = 'collect')

# Create central labels
central_x_label <- wrap_elements(grid::textGrob("SOC", gp = grid::gpar(fontsize = 14, fontface = "bold")))
central_y_label <- wrap_elements(grid::textGrob(" ", rot = 90, gp = grid::gpar(fontsize = 14, fontface = "bold")))

# Adjust the layout
traitSOC_patchwork_plots <- (
  central_y_label +
    plots_combined + 
    plot_spacer()
) / (
  plot_spacer() + 
    central_x_label + 
    plot_spacer()
) +
  plot_layout(widths = c(1, 20), heights = c(20, 1, 1)) +
  plot_annotation(
    title = ' ',
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )
  )

#print the patchwork plot
print(traitSOC_patchwork_plots)

###################################################################################################

#now creating the histogram for all traits together (not averaged)
GRooT_and_SG_joined_wide_main4 <-GRooT_and_SG_joined_wide |>
  dplyr::select(soc, Root_tissue_density, Root_N_concentration, Specific_root_length, Mean_Root_diameter) |>
  glimpse()


# Ensure the soc column is numeric and remove NA values
GRooT_and_SG_joined_wide_main4 <- GRooT_and_SG_joined_wide_main4 |>
  mutate(soc = as.numeric(soc)) |>
  filter(!is.na(soc)) |>
  glimpse()

# Create a histogram of the non-NA values in the soc column
groot_soc_HIST <- ggplot(GRooT_and_SG_joined_wide_main4, aes(x = soc)) +
  geom_histogram(binwidth = 30, fill = "#fcfdbf", color = "black") +
  labs(title = "SOC distribution within top 4 GRooT fine-root traits",
       x = "SOC (mg/g)",
       y = "Frequency") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_rect(fill = "white"),  # White background
    axis.line = element_line(color = "black"),  # Black axis lines
    axis.title = element_text(color = "black"),  # Black axis titles
    axis.text = element_text(color = "black"),  # Black axis text
    axis.text.x = element_text(angle = 90, hjust = 1)  # Adjust x-axis text
  )

# Print the histogram
print(groot_soc_HIST)

###############################################################################################

#now i want to know number of observations categorized as mineral/organic and get a percentage rep for each


#GRooT FR Traits by organic v mineral soils >205 because that is our previously established cutoff point

OvM_trait_grootSG_wide_avgs <- spp_GRooT_and_SG_joined_wide_avgs |>
  mutate(soil_type = if_else(soc_avg > 207.7386, "organic", "mineral"))


OvM_trait_grootSG_wide_avgs <- OvM_trait_grootSG_wide_avgs |>
  filter(!is.na(soil_type))

#final DataFrame
glimpse(OvM_trait_grootSG_wide_avgs)

# Count the total number of values in the soil_type column (excluding NA values)
total_soil_type_count <- OvM_trait_grootSG_wide_avgs |>
  summarize(total = sum(!is.na(soil_type)))


print(total_soil_type_count) #375 total non NA values

# Count how many are 'organic' excluding NA values
organic_count <- OvM_trait_grootSG_wide_avgs |>
  filter(soil_type == "organic") |>
  summarize(total = sum(!is.na(soil_type)))

print(organic_count)   #91

#percent representation: 
91/325 #28.0% organic

# Count how many are 'mineral' excluding NA values
mineral_count <- OvM_trait_grootSG_wide_avgs |>
  filter(soil_type == "mineral") |>
  summarize(total = sum(!is.na(soil_type)))

print(mineral_count)   #284

#percent representation: 
284/325 #87.38% mineral 

#####################################################################################################
#now we will get numbers of observations for each of top 4 FR traits
OvM_GRooT_and_SG_joined_wide_main4 <- GRooT_and_SG_joined_wide_main4 |>
  mutate(soil_type = if_else(soc > 207.7386, "organic", "mineral"))

# now i want to group by trait
names(OvM_GRooT_and_SG_joined_wide_main4)

# Count organic observations for each root trait
trait_organic_counts <- OvM_GRooT_and_SG_joined_wide_main4 |>
  filter(soil_type == "organic") |>
  summarize(
    Root_N_content = sum(!is.na(`Root_N_concentration`)),
    Root_diameter = sum(!is.na(`Mean_Root_diameter`)),
    RTD = sum(!is.na(`Root_tissue_density`)),
    SRL = sum(!is.na(`Specific_root_length`))
  )

print(trait_organic_counts)


# Count mineral observations for each root trait
trait_mineral_counts <- OvM_GRooT_and_SG_joined_wide_main4 |>
  filter(soil_type == "mineral") |>
  summarize(
    Root_N_content = sum(!is.na(`Root_N_concentration`)),
    Root_diameter = sum(!is.na(`Mean_Root_diameter`)),
    RTD = sum(!is.na(`Root_tissue_density`)),
    SRL = sum(!is.na(`Specific_root_length`))
  )


print(trait_mineral_counts)

####################################################################################################

#now we want to do t-tests (while accounting for unequal variance) to see differences in traits mineral v organic

# Filter the data for mineral and organic soils
mineral_soils <- OvM_GRooT_and_SG_joined_wide_main4 |>
  filter(soil_type == "mineral")

organic_soils <- OvM_GRooT_and_SG_joined_wide_main4 |>
  filter(soil_type == "organic")


# Perform t-tests for each trait, here we use welch two sample t test to account for unequal variance!
t_test_Root_N_concentration <- t.test(mineral_soils$Root_N_concentration, organic_soils$Root_N_concentration, na.rm = TRUE)
t_test_Mean_Root_diameter <- t.test(mineral_soils$Mean_Root_diameter, organic_soils$Mean_Root_diameter, na.rm = TRUE)
t_test_Specific_root_length <- t.test(mineral_soils$Specific_root_length, organic_soils$Specific_root_length, na.rm = TRUE)
t_test_Root_tissue_density <- t.test(mineral_soils$Root_tissue_density, organic_soils$Root_tissue_density, na.rm = TRUE)

# Print the results
print(t_test_Root_N_concentration)



print(t_test_Mean_Root_diameter)


print(t_test_Specific_root_length)


print(t_test_Root_tissue_density)

############################################################################################
#set cutoff for mineral vs organic
############################################################################################

OvM_spp_GRooT_and_SG_joined_wide_avgs <- spp_GRooT_and_SG_joined_wide_avgs |>
  mutate(soil_type = if_else(soc_avg > 207.7386, "organic", "mineral")) 

OvM_spp_GRooT_and_SG_joined_wide_avgs <- OvM_spp_GRooT_and_SG_joined_wide_avgs |>
  filter(!is.na(soil_type))

glimpse(OvM_spp_GRooT_and_SG_joined_wide_avgs)

############################################################################################
#time to make PCAs by color of PFTs for mineral+organic, mineral, and organic soil data
############################################################################################

names(OvM_spp_GRooT_and_SG_joined_wide_avgs)
names(spp_GRooT_and_SG_joined_wide)

#growthFrom to dataset
OvM_spp_GRooT_and_SG_joined_wide_avgs <- OvM_spp_GRooT_and_SG_joined_wide_avgs |>
  left_join(
    spp_GRooT_and_SG_joined_wide |> dplyr::select(genus_species, growthForm),
    by = "genus_species"
  )

glimpse(OvM_spp_GRooT_and_SG_joined_wide_avgs)

#select relevant columns, including growthForm and genus_species
all_soils_pca_data <- na.omit(OvM_spp_GRooT_and_SG_joined_wide_avgs[, c("genus_species", "growthForm", "Root_N_avg", "Root_D_avg", "Root_RTD_avg", "SRL_avg")])

names(all_soils_pca_data)
str(all_soils_pca_data)

# Store genus_species and growthForm separately for use in the plot
species <- all_soils_pca_data$genus_species
growthForm <- all_soils_pca_data$growthForm

# Remove these columns from PCA input
all_soils_pca_data$genus_species <- NULL
all_soils_pca_data$growthForm <- NULL

# Perform PCA on the selected root traits
all_soils_pca <- prcomp(all_soils_pca_data, center = TRUE, scale. = TRUE)
summary(all_soils_pca)

# Extract PCA scores and add genus_species and growthForm back
pca_scores <- as.data.frame(all_soils_pca$x)
pca_scores$genus_species <- species
pca_scores$growthForm <- growthForm


# Filter out "fern" and "subshrub" from pca_scores and remove NA values
pca_scores <- na.omit(pca_scores[!(pca_scores$growthForm %in% c("fern", "subshrub")), ])

library(paletteer)
# Define a colorblind-friendly palette using paletteer
color_palette <- paletteer_d("ggthemes::Tableau_10")



#######here we are loading in our personalized ggplot2 function#######
ggbiplot2<-function (pcobj, choices = 1:2, scale = 1, pc.biplot = TRUE, 
                     obs.scale = 1 - scale, var.scale = scale, groups = NULL, 
                     ellipse = FALSE, ellipse.prob = 0.95, labels = NULL, labels.size = 6, 
                     alpha = 0.5, var.axes = TRUE, circle = FALSE, circle.prob = 0.69, 
                     varname.size = 12, varname.adjust = 1.5, varname.abbrev = FALSE,
                     color = "black", size_arrow=2,# <- add new arguments to the function,
                     linetype="solid",alpha_arrow=1)
{
  library(ggplot2)
  library(plyr)
  library(scales)
  library(grid)
  stopifnot(length(choices) == 2)
  if (inherits(pcobj, "prcomp")) {
    nobs.factor <- sqrt(nrow(pcobj$x) - 1)
    d <- pcobj$sdev
    u <- sweep(pcobj$x, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$rotation
  }
  else if (inherits(pcobj, "princomp")) {
    nobs.factor <- sqrt(pcobj$n.obs)
    d <- pcobj$sdev
    u <- sweep(pcobj$scores, 2, 1/(d * nobs.factor), FUN = "*")
    v <- pcobj$loadings
  }
  else if (inherits(pcobj, "PCA")) {
    nobs.factor <- sqrt(nrow(pcobj$call$X))
    d <- unlist(sqrt(pcobj$eig)[1])
    u <- sweep(pcobj$ind$coord, 2, 1/(d * nobs.factor), FUN = "*")
    v <- sweep(pcobj$var$coord, 2, sqrt(pcobj$eig[1:ncol(pcobj$var$coord), 
                                                  1]), FUN = "/")
  }
  else if (inherits(pcobj, "lda")) {
    nobs.factor <- sqrt(pcobj$N)
    d <- pcobj$svd
    u <- predict(pcobj)$x/nobs.factor
    v <- pcobj$scaling
    d.total <- sum(d^2)
  }
  else {
    stop("Expected a object of class prcomp, princomp, PCA, or lda")
  }
  choices <- pmin(choices, ncol(u))
  df.u <- as.data.frame(sweep(u[, choices], 2, d[choices]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^var.scale, FUN = "*")
  df.v <- as.data.frame(v[, choices])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  if (pc.biplot) {
    df.u <- df.u * nobs.factor
  }
  r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  if (obs.scale == 0) {
    u.axis.labs <- paste("standardized PC", choices, sep = "")
  }
  else {
    u.axis.labs <- paste("PC", choices, sep = "")
  }
  u.axis.labs <- paste(u.axis.labs, sprintf("(%0.1f%% explained var.)", 
                                            100 * pcobj$sdev[choices]^2/sum(pcobj$sdev^2)))
  if (!is.null(labels)) {
    df.u$labels <- labels
  }
  if (!is.null(groups)) {
    df.u$groups <- groups
  }
  if (varname.abbrev) {
    df.v$varname <- abbreviate(rownames(v))
  }
  else {
    df.v$varname <- rownames(v)
  }
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  g <- ggplot(data = df.u, aes(x = xvar, y = yvar)) + xlab(u.axis.labs[1]) + 
    ylab(u.axis.labs[2]) + coord_equal()
  if (var.axes) {
    if (circle) {
      theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, 
                                                length = 50))
      circle <- data.frame(xvar = r * cos(theta), yvar = r * 
                             sin(theta))
      g <- g + geom_path(data = circle, color = muted("white"), 
                         size = 1/2, alpha = 1/3)
    }
    g <- g + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                           xend = xvar, yend = yvar), arrow = arrow(length = unit(1/2, 
                                                                                                  "picas")), color = color, linetype=linetype,alpha=alpha_arrow,size=size_arrow)
  }
  if (!is.null(df.u$labels)) {
    if (!is.null(df.u$groups)) {
      g <- g + geom_text(aes(label = labels, color = groups), 
                         size = labels.size)
    }
    else {
      g <- g + geom_text(aes(label = labels), size = labels.size)
    }
  }
  else {
    if (!is.null(df.u$groups)) {
      g <- g + geom_point(aes(color = groups), alpha = alpha)
    }
    else {
      g <- g + geom_point(alpha = alpha)
    }
  }
  if (!is.null(df.u$groups) && ellipse) {
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- cbind(cos(theta), sin(theta))
    ell <- ddply(df.u, "groups", function(x) {
      if (nrow(x) <= 2) {
        return(NULL)
      }
      sigma <- var(cbind(x$xvar, x$yvar))
      mu <- c(mean(x$xvar), mean(x$yvar))
      ed <- sqrt(qchisq(ellipse.prob, df = 2))
      data.frame(sweep(circle %*% chol(sigma) * ed, 2, 
                       mu, FUN = "+"), groups = x$groups[1])
    })
    names(ell)[1:2] <- c("xvar", "yvar")
    g <- g + geom_path(data = ell, aes(color = groups, group = groups))
  }
  if (var.axes) {
    g <- g + geom_text(data = df.v, aes(label = varname, 
                                        x = xvar, y = yvar, angle = angle, hjust = hjust), 
                       color = "black", size = varname.size)
  }
  return(g)
}

library(ggrepel)

# Define the offset factor and arrow labels as provided
offset_factor <- 3.5  # Adjust this factor as needed to position labels beyond arrow tips
arrow_labels <- data.frame(
  PC1 = c(all_soils_pca$rotation["Root_RTD_avg", "PC1"],
          all_soils_pca$rotation["SRL_avg", "PC1"],
          all_soils_pca$rotation["Root_D_avg", "PC1"],
          all_soils_pca$rotation["Root_N_avg", "PC1"]) * offset_factor,
  PC2 = c(all_soils_pca$rotation["Root_RTD_avg", "PC2"],
          all_soils_pca$rotation["SRL_avg", "PC2"],
          all_soils_pca$rotation["Root_D_avg", "PC2"],
          all_soils_pca$rotation["Root_N_avg", "PC2"]) * offset_factor,
  label = c("RTD", "SRL", "D", "N")  # New labels
)

# PCA plot with colorblind-friendly palette
gALL <- ggbiplot2(all_soils_pca, 
                  scale = 1, 
                  obs.scale = 1, 
                  var.scale = 1, 
                  varname.size = 0, 
                  ellipse = TRUE, 
                  ellipse.prob = 0.95, 
                  circle = FALSE, 
                  color = "black", 
                  linetype = "solid", 
                  size_arrow = 1) +
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = growthForm), size = 2, alpha = 0.6) + # Points colored based on growthForm
  scale_color_manual(values = as.vector(color_palette)) +  # Use the colorblind-friendly color palette
  
  # Customize the theme
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.title.x = element_text(size = 12, vjust = -1),
        axis.title.y = element_text(size = 12, vjust = 1),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "right") +  # Show legend
  scale_x_continuous(limits = c(-6, 6)) +  # Adjust axis limits as necessary
  scale_y_continuous(limits = c(-6, 6)) +  # Adjust axis limits as necessary
  
  # Add the arrow labels using ggrepel
  geom_text_repel(data = arrow_labels, aes(x = PC1, y = PC2, label = label), size = 4, color = "black") 

# Print the plot
print(gALL)


############################################################################################

#moving on to doing mineral PCA
############################################################################################
# Step 1: Filter for mineral soils
mineral_soils_data <- OvM_spp_GRooT_and_SG_joined_wide_avgs[OvM_spp_GRooT_and_SG_joined_wide_avgs$soil_type == "mineral", ]  # Replace with the correct soil type column

# Step 2: Select relevant columns, including growthForm and genus_species
mineral_soils_data <- na.omit(mineral_soils_data[, c("genus_species", "growthForm", "Root_N_avg", "Root_D_avg", "Root_RTD_avg", "SRL_avg")])

# Store genus_species and growthForm separately for use in the plot
species <- mineral_soils_data$genus_species
growthForm <- mineral_soils_data$growthForm

# Remove these columns from PCA input
mineral_soils_data$genus_species <- NULL
mineral_soils_data$growthForm <- NULL

# Step 3: Perform PCA on the selected root traits
mineral_soils_data_pca <- prcomp(mineral_soils_data, center = TRUE, scale. = TRUE)
summary(mineral_soils_data_pca)





# Extract PCA scores and add genus_species and growthForm back
pca_scores <- as.data.frame(mineral_soils_data_pca$x)
pca_scores$genus_species <- species
pca_scores$growthForm <- growthForm

# Filter out "fern" and "subshrub" from pca_scores and remove NA values
pca_scores <- na.omit(pca_scores[!(pca_scores$growthForm %in% c("fern", "subshrub")), ])

library(Hmisc)
MIN_spp_avg_bergmannSG_wide_pca<-as.matrix(mineral_soils_data[,c("Root_N_avg", "Root_D_avg", "Root_RTD_avg", "SRL_avg")])
rcorr(MIN_spp_avg_bergmannSG_wide_pca,type=c("pearson"))# produces p values associated with correlated variables




################
#   2.  Run PCA on above data; calculated significance of variables
##############
# Step 1: Filter for mineral soils
mineral_soils_data <- OvM_spp_GRooT_and_SG_joined_wide_avgs[OvM_spp_GRooT_and_SG_joined_wide_avgs$soil_type == "mineral", ]  # Replace with the correct soil type column

# Step 2: Select relevant columns, including growthForm and genus_species
mineral_soils_data <- na.omit(mineral_soils_data[, c("genus_species", "growthForm", "Root_N_avg", "Root_D_avg", "Root_RTD_avg", "SRL_avg")])

# Store genus_species and growthForm separately for use in the plot
species <- mineral_soils_data$genus_species
growthForm <- mineral_soils_data$growthForm

# Remove these columns from PCA input
mineral_soils_data$genus_species <- NULL
mineral_soils_data$growthForm <- NULL

# Step 3: Perform PCA on the selected root traits
mineral_soils_data_pca <- prcomp(mineral_soils_data, center = TRUE, scale. = TRUE)
summary(mineral_soils_data_pca)

# Extract PCA scores and add genus_species and growthForm back
pca_scores <- as.data.frame(mineral_soils_data_pca$x)
pca_scores$genus_species <- species
pca_scores$growthForm <- growthForm

# Filter out "fern" and "subshrub" from pca_scores and remove NA values
pca_scores <- na.omit(pca_scores[!(pca_scores$growthForm %in% c("fern", "subshrub")), ])

library(Hmisc)
MIN_spp_avg_bergmannSG_wide_pca <- as.matrix(mineral_soils_data[, c("Root_N_avg", "Root_D_avg", "Root_RTD_avg", "SRL_avg")])
rcorr(MIN_spp_avg_bergmannSG_wide_pca, type = c("pearson"))  # Produces p-values associated with correlated variables

# Generate PCA for the mineral soil data using the complete dataset
MIN_spp_avg_bergmannSG_wide_pca <- OvM_spp_GRooT_and_SG_joined_wide_avgs %>%
  filter(soil_type == "mineral")

# Prepare the data for PCA (excluding genus_species)
MIN_spp_avg_bergmannSG_wide_pca <- na.omit(OvM_spp_GRooT_and_SG_joined_wide_avgs[, c("genus_species", "Root_N_avg", "Root_D_avg", "Root_RTD_avg", "SRL_avg")])

# Extract genus_species column and store it
MINspp_avg_bergmann_spp <- MIN_spp_avg_bergmannSG_wide_pca$genus_species

# Remove genus_species column for PCA
MIN_spp_avg_bergmannSG_wide_pca$genus_species <- NULL

# Perform PCA
MIN_SPP_avg_bergmann_wide_PCA <- prcomp(MIN_spp_avg_bergmannSG_wide_pca, center = TRUE, scale. = TRUE, retx = TRUE)
summary(MIN_SPP_avg_bergmann_wide_PCA)

# Extract PCA scores and combine with genus_species and soc_avg
MINpca_scores <- as.data.frame(MIN_SPP_avg_bergmann_wide_PCA$x)
MINpca_scores$genus_species <- MINspp_avg_bergmann_spp
MINpca_scores$soc_avg <- MIN_spp_avg_bergmannSG_wide_pca$soc_avg

# Plot PCA result
plot(MIN_SPP_avg_bergmann_wide_PCA, type = "l")
biplot(MIN_SPP_avg_bergmann_wide_PCA)

# Significance calculation for the variables in PCA loadings
sigpca2 <- function(x, permutations = 1000) {
  pcnull <- princomp(x, cor = TRUE)
  res <- pcnull$loadings
  out <- matrix(0, nrow = nrow(res), ncol = ncol(res))
  N <- nrow(x)
  for (i in 1:permutations) {
    pc <- princomp(x[sample(N, replace = TRUE), ], cor = TRUE)
    pred <- predict(pc, newdata = x)
    r <- cor(pcnull$scores, pred)
    k <- apply(abs(r), 2, which.max)
    reve <- sign(diag(r[k, ]))
    sol <- pc$loadings[, k]
    sol <- sweep(sol, 2, reve, "*")
    out <- out + ifelse(res > 0, sol <= 0, sol >= 0)
  }
  out / permutations
}

output <- sigpca2(MIN_spp_avg_bergmannSG_wide_pca, permutations = 10000)
output

# Generate the PCA biplot with custom labels and visualization
arrow_labels_mineral <- data.frame(
  PC1 = c(mineral_soils_data_pca$rotation["Root_RTD_avg", "PC1"],
          mineral_soils_data_pca$rotation["SRL_avg", "PC1"],
          mineral_soils_data_pca$rotation["Root_D_avg", "PC1"],
          mineral_soils_data_pca$rotation["Root_N_avg", "PC1"]) * 0.5,
  PC2 = c(mineral_soils_data_pca$rotation["Root_RTD_avg", "PC2"],
          mineral_soils_data_pca$rotation["SRL_avg", "PC2"],
          mineral_soils_data_pca$rotation["Root_D_avg", "PC2"],
          mineral_soils_data_pca$rotation["Root_N_avg", "PC2"]) * 0.5,
  label = c("RTD", "SRL", "D", "N")  # New labels for mineral soils
)

# Final PCA plot using ggplot2
gMIN <- ggbiplot2(mineral_soils_data_pca,
                  scale = 1,
                  obs.scale = 1,
                  var.scale = 1,
                  varname.size = 0,
                  ellipse = TRUE,
                  ellipse.prob = 0.95,
                  circle = FALSE,
                  color = "black",
                  linetype = "solid",
                  size_arrow = 1) +
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = growthForm), size = 2, alpha = 0.6) +  # Points colored based on growthForm
  scale_color_manual(values = color_palette ) +  # Adjust the color palette to your needs
  
  # Customize the theme
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.title.x = element_text(size = 12, vjust = -1),
        axis.title.y = element_text(size = 12, vjust = 1),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "right") +  # Show legend
  scale_x_continuous(limits = c(-6, 6)) +  # Adjust axis limits as necessary
  scale_y_continuous(limits = c(-6, 6)) +  # Adjust axis limits as necessary
  
  # Add the arrow labels using ggrepel
  geom_text_repel(data = arrow_labels_mineral, aes(x = PC1, y = PC2, label = label), size = 4, color = "black")

# Print the PCA plot
print(gMIN)



########################################
#now doing organic soils pca
#######################################
organic_soils_data <- OvM_spp_GRooT_and_SG_joined_wide_avgs[OvM_spp_GRooT_and_SG_joined_wide_avgs$soil_type == "organic", ]  # Replace with the correct soil type column

# Step 2: Select relevant numeric columns for PCA (excluding genus_species and growthForm)
organic_soils_data_pca <- na.omit(organic_soils_data[, c("Root_N_avg", "Root_D_avg", "Root_RTD_avg", "SRL_avg")])

# Step 3: Perform PCA on the numeric data
organic_soils_data_pca_results <- prcomp(organic_soils_data_pca, center = TRUE, scale. = TRUE)

# Display PCA summary
summary(organic_soils_data_pca_results)

names(pca_scores)

dim(orgpca_scores)
names(orgpca_scores)

# Step 4: Extract PCA scores and add genus_species and growthForm back
orgpca_scores <- as.data.frame(organic_soils_data_pca_results$x)
organic_soils_data_pca_results$growthForm <- growthForm

# Filter out "fern" and "subshrub" from PCA scores and remove NA values
orgpca_scores <- na.omit(pca_scores[!(pca_scores$growthForm %in% c("fern", "subshrub")), ])

# Step 5: Correlation analysis for root traits
ORG_spp_avg_bergmannSG_wide_pca <- as.matrix(organic_soils_data[, c("Root_N_avg", "Root_D_avg", "Root_RTD_avg", "SRL_avg")])
rcorr(ORG_spp_avg_bergmannSG_wide_pca, type = c("pearson"))  # Correlation matrix


# Calculate significance of variables in loadings
sigpca2 <- function (x, permutations = 1000) {
  pcnull <- princomp(x, cor = TRUE)
  res <- pcnull$loadings
  out <- matrix(0, nrow = nrow(res), ncol = ncol(res))
  N <- nrow(x)
  for (i in 1:permutations) {
    pc <- princomp(x[sample(N, replace = TRUE), ], cor = TRUE)
    pred <- predict(pc, newdata = x)
    r <- cor(pcnull$scores, pred)
    k <- apply(abs(r), 2, which.max)
    reve <- sign(diag(r[k, ]))
    sol <- pc$loadings[, k]
    sol <- sweep(sol, 2, reve, "*")
    out <- out + ifelse(res > 0, sol <= 0, sol >= 0)
  }
  out / permutations
}

output <- sigpca2(organic_soils_data_pca, permutations = 10000)
output

# Step 6: Plot PCA results
arrow_labels <- data.frame(
  PC1 = c(organic_soils_data_pca_results$rotation["Root_RTD_avg", "PC1"],
          organic_soils_data_pca_results$rotation["SRL_avg", "PC1"],
          organic_soils_data_pca_results$rotation["Root_D_avg", "PC1"],
          organic_soils_data_pca_results$rotation["Root_N_avg", "PC1"]) * offset_factor,
  PC2 = c(organic_soils_data_pca_results$rotation["Root_RTD_avg", "PC2"],
          organic_soils_data_pca_results$rotation["SRL_avg", "PC2"],
          organic_soils_data_pca_results$rotation["Root_D_avg", "PC2"],
          organic_soils_data_pca_results$rotation["Root_N_avg", "PC2"]) * offset_factor,
  label = c("RTD", "SRL", "D", "N")  # Labels for root traits
)

names(organic_soils_data_pca_results)


ORGpca_scores <- merge(orgpca_scores, OvM_spp_GRooT_and_SG_joined_wide_avgs[, c("genus_species", "growthForm")], by = "genus_species", all.x = TRUE)

ORGpca_scores$growthForm <- as.factor(orgpca_scores$growthForm)

# Filter out rows where growthForm is "fern" or "subshrub"
orgpca_scores <- orgpca_scores %>%
  filter(!growthForm %in% c("fern", "subshrub"))



gORG<- ggbiplot2(organic_soils_data_pca_results,
                 scale = 1, 
                 obs.scale = 1, 
                 var.scale = 1, 
                 varname.size = 0, 
                 ellipse = TRUE, 
                 ellipse.prob = 0.95, 
                 circle = FALSE, 
                 color = "black", 
                 linetype = "solid", 
                 size_arrow = 1) 
gORG


gORG <- gORG +  
  geom_point(data = orgpca_scores, aes(x = PC1, y = PC2, color = growthForm), size = 2, alpha = 0.6) + # Points colored based on growthForm
  scale_color_manual(values = as.vector(color_palette)) +  
  
  # Customize the theme
  theme_bw() +
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.ticks = element_line(colour = "black", linewidth = 0.5),
        axis.title.x = element_text(size = 12, vjust = -1),
        axis.title.y = element_text(size = 12, vjust = 1),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "right") +  # Show legend
  scale_x_continuous(limits = c(-6, 6)) +  # Adjust axis limits as necessary
  scale_y_continuous(limits = c(-6, 6)) +  
  # Add the arrow labels using ggrepel
  geom_text_repel(data = arrow_labels, aes(x = PC1, y = PC2, label = label), size = 4, color = "black")
print(gORG)


names(orgpca_scores)
view(orgpca_scores)

########run this then re-run pca plot codes to ensure consistent colors
# Get unique growthForm values from all datasets
unique_growth_forms <- unique(c(orgpca_scores$growthForm, MINpca_scores$growthForm, pca_scores$growthForm))


# Create a color palette (using RColorBrewer or manually defined colors)
library(RColorBrewer)
color_palette <- setNames(brewer.pal(length(unique_growth_forms), "Set3"), unique_growth_forms)

# Preview the palette
print(color_palette)


#############################################################################
#here i am making tables eigenvectors and importance of compnenets for each PCA
# Function to extract eigenvector and importance tables
extract_pca_results <- function(pca_result) {
  # Eigenvectors (rotation matrix)
  eigenvectors <- as.data.frame(pca_result$rotation)
  colnames(eigenvectors) <- paste0("PC", 1:ncol(eigenvectors))
  
  # Importance of components
  importance <- summary(pca_result)$importance
  importance_df <- as.data.frame(t(importance))
  colnames(importance_df) <- c("Standard Deviation", "Proportion of Variance", "Cumulative Proportion")
  
  list(eigenvectors = eigenvectors, importance = importance_df)
}

# Extract results for all soils
all_soils_pca_results <- extract_pca_results(all_soils_pca)
all_soils_eigenvectors <- all_soils_pca_results$eigenvectors
all_soils_importance <- all_soils_pca_results$importance

# Extract results for mineral soils
mineral_soils_pca_results <- extract_pca_results(mineral_soils_data_pca)
mineral_soils_eigenvectors <- mineral_soils_pca_results$eigenvectors
mineral_soils_importance <- mineral_soils_pca_results$importance

# Assuming you also performed PCA for organic soils
organic_soils_pca_results <- extract_pca_results(organic_soils_data_pca_results)
organic_soils_eigenvectors <- organic_soils_data_pca_results$eigenvectors
organic_soils_importance <- organic_soils_pca_results$importance

# Display the tables
list(
  All_Soils = list(Eigenvectors = all_soils_eigenvectors, Importance = all_soils_importance),
  Mineral_Soils = list(Eigenvectors = mineral_soils_eigenvectors, Importance = mineral_soils_importance),
  Organic_Soils = list(Eigenvectors = organic_soils_eigenvectors, Importance = organic_soils_importance)
)





##############################################################################
#now doing three-way anovas to see which factor impacts traits signifciantly 
##############################################################################
#doing aov for trait data 3 way to show overall 'what matters': soil type, mfat, pft....?

names(spp_GRooT_and_SG_joined_wide)
spp_GRooT_and_SG_joined_wide <- spp_GRooT_and_SG_joined_wide|>
  mutate(soil_type = if_else(soc_avg > 207.7386, "organic", "mineral")) 


ggplot(data=spp_GRooT_and_SG_joined_wide, aes(x= sqrt(Root_N_avg)))+geom_histogram()
rootN_aov<-aov(data = spp_GRooT_and_SG_joined_wide, sqrt(Root_N_avg) ~ growthForm * mycorrhizalAssociationTypeFungalRoot * soil_type)
summary(rootN_aov)

plot(rootN_aov)




ggplot(data=spp_GRooT_and_SG_joined_wide,aes(x=log(Root_RTD_avg)))+geom_histogram()
RTD_aov<-aov(data = spp_GRooT_and_SG_joined_wide, log(Root_RTD_avg) ~ growthForm * mycorrhizalAssociationTypeFungalRoot* soil_type)
summary(RTD_aov)

plot(RTD_aov)


ggplot(data=spp_GRooT_and_SG_joined_wide,aes(x=log(Root_D_avg)))+geom_histogram()
D_aov<-aov(data = spp_GRooT_and_SG_joined_wide, log(Root_D_avg) ~ growthForm * mycorrhizalAssociationTypeFungalRoot* soil_type)
summary(D_aov)

plot(D_aov)  #also try sqrt as well??



ggplot(data=spp_GRooT_and_SG_joined_wide,aes(x= log(SRL_avg)))+geom_histogram()
SRL_aov<-aov(data = spp_GRooT_and_SG_joined_wide, log(SRL_avg) ~ growthForm * mycorrhizalAssociationTypeFungalRoot* soil_type)
summary(SRL_aov)

plot(SRL_aov)


################################################
#now doing box plots for pft and then for MFAT
################################################

#PFT
 
names(GRooT_and_SG_joined_wide)
#add soiltype
GRooT_and_SG_joined_wide<- GRooT_and_SG_joined_wide |>
  mutate(soil_type = if_else(soc > 207.7386, "organic", "mineral")) 

library(dplyr)

filtered_pft <- GRooT_and_SG_joined_wide |>
  dplyr::select(growthForm, Specific_root_length , Root_N_concentration, Mean_Root_diameter, Root_tissue_density, soil_type)


filtered_pft <- GRooT_and_SG_joined_wide |> filter(!is.na(`growthForm`)) |>
                                                     filter(!is.na(`soil_type`))
glimpse(filtered_pft)

# Convert growthForm to a factor with specified levels
filtered_pft$growthForm <- factor(filtered_pft$growthForm, 
                                  levels = c("graminoid", "herb", "herb/shrub", "shrub", "shrub/tree", "tree"))


names(filtered_pft)

rootN_pft<- ggplot(filtered_pft, aes(x = growthForm, y = Root_N_concentration , fill = growthForm)) +  
  stat_summary(fun = "mean", geom = "bar", color = "black", width = 0.7) + 
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.3, color = "black") +
  labs(title = "",      
       x = "PFT",       
       y = "Root N ") + 
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),       
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_viridis_d(option = "magma", name = "PFT")


print(rootN_pft)

###################
#now let's split this data into mineral versus organic using cutoff from breakpoint
names(filtered_pft)
filtered_pft <- filtered_pft |>
  filter(!is.na(soc)) |>
  mutate(soil_type = ifelse(soc > 205, "organic", "mineral"))

# Count the number of datapoints for each PFT category for each soil type

pft_counts <- filtered_pft |>
  dplyr::group_by(growthForm, soil_type) |>
 dplyr:: summarize(
    Root_N_count = sum(!is.na(Root_N_concentration)),
    Root_D_count = sum(!is.na(Mean_Root_diameter)),
    Root_RTD_count = sum(!is.na(Root_tissue_density)),
    SRL_count = sum(!is.na(Specific_root_length)),
    SOC_count = sum(!is.na(soc)))

print(pft_counts)


# Root N concentration by PFT, faceted by soil type and i want to display counts
# Join the counts to the original data

filtered_pft <- filtered_pft |>
  filter(!growthForm %in% c("subshrub", "fern") & !is.na(growthForm))


view(filtered_pft)


# Root N by PFT, boxplot
rootN_pft_facet <- ggplot(filtered_pft, aes(x = growthForm, y = Root_N_concentration, fill = soil_type)) +  
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) + 
  labs(title = "", x = "PFT", y = "Root N") + 
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),       
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("mineral" = "#7C4D79", "organic" = "#EEC9E5"), 
                    name = "Soil Type") 
print(rootN_pft_facet)

pft_counts <- filtered_pft |>
  dplyr::group_by(growthForm, soil_type) |>
  dplyr::summarize(
    Root_N_count = sum(!is.na(Root_N_concentration)),
    Root_D_count = sum(!is.na(Mean_Root_diameter)),
    Root_RTD_count = sum(!is.na(Root_tissue_density)),
    SRL_count = sum(!is.na(Specific_root_length)),
    SOC_count = sum(!is.na(soc)))


# Root D by PFT, boxplot
rootD_pft_facet <- ggplot(filtered_pft, aes(x = growthForm, y = Mean_Root_diameter, fill = soil_type)) +  
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) + 
  labs(title = "", x = "PFT", y = "Root D") + 
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_manual(values = c("mineral" = "#7C4D79", "organic" = "#EEC9E5"), 
                    name = "Soil Type") 
print(rootD_pft_facet)

# RTD by PFT, boxplot
RTD_pft_facet <- ggplot(filtered_pft, aes(x = growthForm, y = Root_tissue_density, fill = soil_type)) +  
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) + 
  labs(title = "", x = "PFT", y = "RTD") + 
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_manual(values = c("mineral" = "#7C4D79", "organic" = "#EEC9E5"), 
                    name = "Soil Type") 
print(RTD_pft_facet)

# SRL by PFT, boxplot
SRL_pft_facet <- ggplot(filtered_pft, aes(x = growthForm, y = Specific_root_length, fill = soil_type)) +  
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) + 
  labs(title = "", x = "PFT", y = "SRL") + 
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"), 
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_manual(values = c("mineral" = "#7C4D79", "organic" = "#EEC9E5"), 
                    name = "Soil Type") 
print(SRL_pft_facet)

# Combine the plots using the patchwork package
SW_patchwork_plots <- (rootN_pft_facet + rootD_pft_facet + RTD_pft_facet + SRL_pft_facet +
                         plot_layout(guides = 'collect') +  # Collect legends
                         plot_annotation(title = 'Most commonly measured fine-root traits among PFTs by soil type'))
print(SW_patchwork_plots)



# Root N concentration
table_Root_N <- table(filtered_pft$growthForm[!is.na(filtered_pft$Root_N_concentration)], 
                      filtered_pft$soil_type[!is.na(filtered_pft$Root_N_concentration)])



# Root Diameter
table_Root_D <- table(filtered_pft$growthForm[!is.na(filtered_pft$Mean_Root_diameter)], 
                      filtered_pft$soil_type[!is.na(filtered_pft$Mean_Root_diameter)])



# Root Tissue Density
table_Root_RTD <- table(filtered_pft$growthForm[!is.na(filtered_pft$Root_tissue_density)], 
                        filtered_pft$soil_type[!is.na(filtered_pft$Root_tissue_density)])


# Specific Root Length
table_SRL <- table(filtered_pft$growthForm[!is.na(filtered_pft$Specific_root_length)], 
                   filtered_pft$soil_type[!is.na(filtered_pft$Specific_root_length)])



# SOC
table_SOC <- table(filtered_pft$growthForm[!is.na(filtered_pft$soc)], 
                   filtered_pft$soil_type[!is.na(filtered_pft$soc)])






#now for anovas and tukeys for traits and pft

#doing anovas
#rootN

# ANOVA
modelNpft <- aov(Root_N_concentration ~ growthForm * soil_type, data = filtered_pft)
summary(modelNpft)


# Post hoc Tukey HSD test
tukey_Npft <- TukeyHSD(modelNpft, "growthForm:soil_type")
print(tukey_Npft)

# Use emmeans to get estimated marginal means and letters
library(emmeans)
library(multcompView)
library(multcomp)
emmeans_Npft <- emmeans(modelNpft, ~ growthForm * soil_type)
cld_Npft <- cld(emmeans_Npft, Letters = letters)


emmeans_interactionpft <- emmeans(modelNpft, ~ growthForm * soil_type)

# View the pairwise comparisons

print(emmeans_interactionpft)



# Run Tukey-adjusted comparisons
tukey_test <- pairs(emmeans_interactionpft, adjust = "tukey")

# View pairwise comparison results
print(tukey_test)

cld_results <- cld(emmeans_interactionpft, Letters = letters)
# View the letter results
print(cld_results)

# Pairwise comparison for organic vs mineral within each mycorrhizal type
pairwise_comparison <- contrast(emmeans_interactionpft, method = "pairwise", by = "growthForm")





#now for diameter
modelDpft <- aov(`Mean_Root_diameter` ~ growthForm * soil_type, data = filtered_pft)
summary(modelDpft)

# Tukey HSD post hoc
tukey_diameter <- TukeyHSD(modelDpft, "growthForm:soil_type")
print(tukey_diameter)

# Estimated marginal means and CLD for Mean_Root_diameter
emmeans_diameter <- emmeans(modelDpft, ~ growthForm * soil_type)
cld_diameter <- cld(emmeans_diameter, Letters = letters)
print(cld_diameter)



# Pairwise comparisons by soil type within each growthForm
pairwise_comparison_diameter <- contrast(emmeans_diameter, method = "pairwise", by = "growthForm")
print(pairwise_comparison_diameter)









# ANOVA for Root_tissue_density
model_RTD <- aov(Root_tissue_density ~ growthForm * soil_type, data = filtered_pft)
summary(model_RTD)

# Tukey HSD post hoc test
tukey_RTD <- TukeyHSD(model_RTD, "growthForm:soil_type")
print(tukey_RTD)

# Calculate emmeans and compact letter display (CLD) for Root_tissue_density
emmeans_RTD <- emmeans(model_RTD, ~ growthForm * soil_type)
cld_RTD <- cld(emmeans_RTD, Letters = letters)
print(cld_RTD)




# Pairwise comparison for organic vs mineral within each growthForm for Root_tissue_density
pairwise_comparison_RTD <- contrast(emmeans_RTD, method = "pairwise", by = "growthForm")
print(pairwise_comparison_RTD)






# ANOVA for Specific_root_length
model_SRL <- aov(Specific_root_length ~ growthForm * soil_type, data = filtered_pft)
summary(model_SRL)
# 


# Tukey HSD post hoc test
tukey_SRL <- TukeyHSD(model_SRL, "growthForm:soil_type")
print(tukey_SRL)

# Calculate emmeans and compact letter display (CLD) for Specific_root_length
emmeans_SRL <- emmeans(model_SRL, ~ growthForm * soil_type)
cld_SRL <- cld(emmeans_SRL, Letters = letters)
print(cld_SRL)


# Pairwise comparison for organic vs mineral within each growthForm for Specific_root_length
pairwise_comparison_SRL <- contrast(emmeans_SRL, method = "pairwise", by = "growthForm")
print(pairwise_comparison_SRL)


#####################################################################################################
#####now doing this for MFATs

mf_grootbergmannSG_wide<- GRooT_and_SG_joined_wide |>
  dplyr::select(mycorrhizalAssociationTypeFungalRoot, Specific_root_length , Root_N_concentration, Mean_Root_diameter, Root_tissue_density, soil_type)



 



mf_grootbergmannSG_wide_rm.NA <-mf_grootbergmannSG_wide |>
  filter(!is.na(`mycorrhizalAssociationTypeFungalRoot`) & !is.na(`soil_type`)) |>
  view()

#10 observations of 6 variables


names(mf_grootbergmannSG_wide_rm.NA)

# Convert  to a factor with specified levels
mf_grootbergmannSG_wide_rm.NA$mycorrhizalAssociationTypeFungalRoot <- factor(mf_grootbergmannSG_wide_rm.NA$mycorrhizalAssociationTypeFungalRoot, 
                                                                                  levels = c("AM", "EcM", "EcM-AM", "ErM", "NM", "NM-AM", "OM" ,"species-specific: AM or rarely EcM-AM or AM", "uncertain"))


# Filter out the rows where mycorrhizalAssociationType is 'uncertain' or 'species-specific'
mf_grootbergmannSG_wide_rm.NA <- mf_grootbergmannSG_wide_rm.NA[!(mf_grootbergmannSG_wide_rm.NA$mycorrhizalAssociationTypeFungalRoot %in% 
                                                         c("uncertain", "species-specific: AM or rarely EcM-AM or AM", "OM", "EcM-AM", "NM-AM", "ErM")), ]


rootN_mf <- ggplot(mf_grootbergmannSG_wide_rm.NA, aes(x = mycorrhizalAssociationTypeFungalRoot, y = Root_N_concentration, fill = soil_type)) +  
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) + 
  labs(title = "",      
       x = " ",       
       y = "Root N") + 
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),       
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("mineral" = "#7C4D79", "organic" = "#EEC9E5"), 
                    name = "Soil Type") 

print(rootN_mf)





rootD_mf <- ggplot(mf_grootbergmannSG_wide_rm.NA, aes(x = mycorrhizalAssociationTypeFungalRoot, y = Mean_Root_diameter, fill = soil_type)) +  
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) + 
  labs(title = "",      
       x = " ",       
       y = "Root D") + 
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),       
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("mineral" = "#7C4D79", "organic" = "#EEC9E5"), 
                    name = "Soil Type") 

print(rootD_mf)





rootRTD_mf <- ggplot(mf_grootbergmannSG_wide_rm.NA, aes(x = mycorrhizalAssociationTypeFungalRoot, y = Root_tissue_density, fill = soil_type)) +  
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) + 
  labs(title = "",      
       x = " ",       
       y = "RTD") + 
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),       
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("mineral" = "#7C4D79", "organic" = "#EEC9E5"), 
                    name = "Soil Type") 

print(rootRTD_mf)




rootSRL_mf <- ggplot(mf_grootbergmannSG_wide_rm.NA, aes(x = mycorrhizalAssociationTypeFungalRoot, y = Specific_root_length, fill = soil_type)) +  
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.7) + 
  labs(title = "",      
       x = " ",       
       y = "SRL") + 
  theme_minimal() +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = "white"),       
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("mineral" = "#7C4D79", "organic" = "#EEC9E5"), 
                    name = "Soil Type") 

print(rootSRL_mf)




library(patchwork)
library(plotly)

# Combine the plots using the patchwork package
SW_patchwork_plots <- (rootN_mf + rootD_mf + rootRTD_mf + rootSRL_mf  +
                         plot_layout(guides = 'collect') +  # Collect legends
                         plot_annotation(title = 'Most commonly measured fine-root traits among mycorrhizal association types'))

# Print the patchwork plot

print(SW_patchwork_plots)

# Summarize the count of observations by mycorrhizal type and soil type
mfa_counts <- mf_grootbergmannSG_wide_rm.NA |>
  dplyr::group_by(mycorrhizalAssociationTypeFungalRoot , soil_type) |>
 dplyr::summarize(
    Root_N_count = sum(!is.na(Root_N_concentration)),
    Root_D_count = sum(!is.na(Mean_Root_diameter)),
    Root_RTD_count = sum(!is.na(Root_tissue_density)),
    SRL_count = sum(!is.na(Specific_root_length)),
    SOC_coutn = sum(!is.na(soil_type)))




#doing anovas

#rootN

names(mf_grootbergmannSG_wide_rm.NA)
modelN <- aov(Root_N_concentration ~ mycorrhizalAssociationTypeFungalRoot * soil_type, data = GRooT_and_SG_joined_wide)

summary(modelN)





#Post hoc comparisons for interaction effects
emmeans_interaction <- emmeans(modelN, ~ mycorrhizalAssociationTypeFungalRoot * soil_type)

# View the pairwise comparisons

print(emmeans_interaction)



# Run Tukey-adjusted comparisons
tukey_test <- pairs(emmeans_interaction, adjust = "tukey")

# View pairwise comparison results
print(tukey_test)

cld_results <- cld(emmeans_interaction, Letters = letters)
# View the letter results
print(cld_results)



# Pairwise comparison for organic vs mineral within each mycorrhizal type
pairwise_comparison <- contrast(emmeans_interaction, method = "pairwise", by = "mycorrhizalAssociationTypeFungalRoot")



#now for diameter
modelD <- aov(`Mean_Root_diameter` ~ mycorrhizalAssociationTypeFungalRoot * soil_type, data = GRooT_and_SG_joined_wide)
summary(modelD)


emmeans_interaction_D <- emmeans(modelD, ~ mycorrhizalAssociationTypeFungalRoot * soil_type)


# View pairwise comparisons
print(emmeans_interaction_D)


# Run Tukey-adjusted comparisons
tukey_test_D <- pairs(emmeans_interaction_D, adjust = "tukey")
print(tukey_test_D)

# View letter results
cld_results_D <- cld(emmeans_interaction_D, Letters = letters)
print(cld_results_D)

# Pairwise comparison for organic vs mineral within each mycorrhizal type

pairwise_comparison_D <- contrast(emmeans_interaction_D, method = "pairwise", by = "mycorrhizalAssociationTypeFungalRoot")
print(pairwise_comparison_D)





#now for RTD
modelRTD <- aov(`Root_tissue_density` ~ mycorrhizalAssociationTypeFungalRoot * soil_type , data = GRooT_and_SG_joined_wide)
summary(modelRTD)



# Post hoc comparisons for interaction effects
emmeans_interaction_RTD <- emmeans(modelRTD, ~ mycorrhizalAssociationTypeFungalRoot * soil_type)

# View pairwise comparisons
print(emmeans_interaction_RTD)


# Run Tukey-adjusted comparisons
tukey_test_RTD <- pairs(emmeans_interaction_RTD, adjust = "tukey")

print(tukey_test_RTD)


# View letter results
cld_results_RTD <- cld(emmeans_interaction_RTD, Letters = letters)
print(cld_results_RTD)


# Pairwise comparison for organic vs mineral within each mycorrhizal type
pairwise_comparison_RTD <- contrast(emmeans_interaction_RTD, method = "pairwise", by = "mycorrhizalAssociationTypeFungalRoot")
print(pairwise_comparison_RTD)







#now for SRL
modelSRL <- aov(`Specific_root_length` ~ mycorrhizalAssociationTypeFungalRoot * soil_type, data = GRooT_and_SG_joined_wide)
summary(modelSRL)

# Run the ANOVA for Specific Root Length
modelSRL <- aov(Specific_root_length ~ mycorrhizalAssociationTypeFungalRoot * soil_type, data = GRooT_and_SG_joined_wide)
summary(modelSRL)


# Post hoc comparisons for interaction effects
emmeans_interaction_SRL <- emmeans(modelSRL, ~ mycorrhizalAssociationTypeFungalRoot * soil_type)

# View pairwise comparisons
print(emmeans_interaction_SRL)


# Run Tukey-adjusted comparisons
tukey_test_SRL <- pairs(emmeans_interaction_SRL, adjust = "tukey")
print(tukey_test_SRL)

# View letter results
cld_results_SRL <- cld(emmeans_interaction_SRL, Letters = letters)
print(cld_results_SRL)


# Pairwise comparison for organic vs mineral within each mycorrhizal type
pairwise_comparison_SRL <- contrast(emmeans_interaction_SRL, method = "pairwise", by = "mycorrhizalAssociationTypeFungalRoot")
print(pairwise_comparison_SRL)


##############################################################################################
#now moving on to figures/tables in the supplemental materials
##############################################################################################

#first, re-run code lines 1 to 48!!! and use GRooT_and_SG_joined df.

names(GRooT_and_SG_joined)
unique(GRooT_and_SG_joined$traitName)

GRooT_and_SG_joined <- GRooT_and_SG_joined |>
  pivot_wider(names_from = traitName, values_from = traitValue, values_fn = function(x) x[1])


groottrait_counts <- GRooT_and_SG_joined |>
  plyr::summarize(Root_production = sum(!is.na(`Root_production`)),
            Root_Ca = sum(!is.na(`Root_Ca_concentration`)),
            Root_K = sum(!is.na(`Root_K_concentration`)),
            Root_Mg = sum(!is.na(`Root_Mg_concentration`)),
            Root_P = sum(!is.na(`Root_P_concentration`)),
            Root_length_density_vol = sum(!is.na(`Root_length_density_volume`)),
            Root_C_N_ratio = sum(!is.na(`Root_C_N_ratio`)),
            Specific_root_resp = sum(!is.na(`Specific_root_respiration`)),
            Root_turnover_rate = sum(!is.na(`Root_turnover_rate`)),
            Root_N= sum(!is.na(`Root_N_concentration`)),
            Root_mass_density = sum(!is.na(`Root_mass_density`)),
            Root_C= sum(!is.na(`Root_C_concentration`)),
            Root_Mn = sum(!is.na(`Root_Mn_concentration`)),
            Fine_root_mass_leaf_mass_ratio = sum(!is.na(`Fine_root_mass_leaf_mass_ratio`)),
            Root_mass_fraction = sum(!is.na(`Root_mass_fraction`)),
            Root_D = sum(!is.na(`Mean_Root_diameter`)),
            RTD = sum(!is.na(`Root_tissue_density`)),
            Specific_root_area = sum(!is.na(`Specific_root_area`)),
            SRL = sum(!is.na(`Specific_root_length`)),
            Root_stele_diameter = sum(!is.na(`Root_stele_diameter`)),
            Root_stele_fraction = sum(!is.na(`Root_stele_fraction`)),
            Root_branching_density = sum(!is.na(`Root_branching_density`)),
            Root_lignin_concentration = sum(!is.na(`Root_lignin_concentration`)),
            Root_total_structural_carb_concentration = sum(!is.na(`Root_total_structural_carbohydrate_concentration`)),
            Root_litter_mass_loss_rate = sum(!is.na(`Root_litter_mass_loss_rate`)),
            Root_lifespan_median = sum(!is.na(`Root_lifespan_median`)),
            Root_vessel_diameter = sum(!is.na(`Root_vessel_diameter`)),
            Root_mycorrhizal_col = sum(!is.na(`Root_mycorrhizal colonization`)),
            Root_branching_ratio = sum(!is.na(`Root_branching_ratio`)),
            Root_dry_matter_content = sum(!is.na(`Root_dry_matter_content`)),
            Root_cortex_thickness = sum(!is.na(`Root_cortex_thickness`)),
            Coarse_root_fine_root_mass_ratio = sum(!is.na(`Coarse_root_fine_root_mass_ratio`)),
            Rooting_depth = sum(!is.na(`Rooting_depth`)),
            Root_lifespan_mean = sum(!is.na(`Root_lifespan_mean`)),
            Net_nitrogen_uptake_rate = sum(!is.na(`Net_nitrogen_uptake_rate`)),
            Root_xylem_vessel_number = sum(!is.na(`Root_xylem_vessel_number`)),
            Root_N_P_ratio = sum(!is.na(`Root_N_P_ratio`))    )

print(groottrait_counts)



#Convert to data frame
groottrait_table <- data.frame(
  Trait = names(groottrait_counts),
  Count = as.numeric(groottrait_counts)
)


# Print as arranged table

groottrait_table <- groottrait_table |>
  dplyr::arrange(desc(Count))

#########################################################################################
#########################################################################################
#here is the t test figures code for figure s2

# Custom color scheme
custom_colors <- c("#51127c", "#b73779")
names(OvM_GRooT_and_SG_joined_wide_main4)

# Calculate the means and standard errors
soil_summary <- OvM_GRooT_and_SG_joined_wide_main4 |>
  dplyr::group_by(soil_type) |>
  dplyr::summarize(
    avg_Root_N_concentration = mean(Root_N_concentration, na.rm = TRUE),
    se_Root_N_concentration = sd(Root_N_concentration, na.rm = TRUE) / sqrt(n()),
    avg_Mean_Root_diameter = mean(Mean_Root_diameter, na.rm = TRUE),
    se_Mean_Root_diameter = sd(Mean_Root_diameter, na.rm = TRUE) / sqrt(n()),
    avg_Specific_root_length = mean(Specific_root_length, na.rm = TRUE),
    se_Specific_root_length = sd(Specific_root_length, na.rm = TRUE) / sqrt(n()),
    avg_Root_tissue_density = mean(Root_tissue_density, na.rm = TRUE),
    se_Root_tissue_density = sd(Root_tissue_density, na.rm = TRUE) / sqrt(n())
  )

# Function to create the plots with custom theme, color, and standard error bars
plot_trait <- function(data, trait_mean, trait_se, trait_name) {
  ggplot(data, aes(x = soil_type, y = get(trait_mean), fill = soil_type)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    geom_errorbar(aes(ymin = get(trait_mean) - get(trait_se), ymax = get(trait_mean) + get(trait_se)),
                  position = position_dodge(width = 0.7), width = 0.25) +
    scale_fill_manual(values = custom_colors) +
    labs(x = " ", y = "Average", fill = "Soil Type",
         title = paste(trait_name)) +
    theme_minimal(base_size = 15) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      legend.position = "none"  # Suppress legend in individual plots
    )
}

# Create individual plots with standard error bars
plot_Root_N_concentration <- plot_trait(soil_summary, "avg_Root_N_concentration", "se_Root_N_concentration", "Root N Concentration")
plot_Mean_Root_diameter <- plot_trait(soil_summary, "avg_Mean_Root_diameter", "se_Mean_Root_diameter", "Mean Root Diameter")
plot_Specific_root_length <- plot_trait(soil_summary, "avg_Specific_root_length", "se_Specific_root_length", "Specific Root Length")
plot_Root_tissue_density <- plot_trait(soil_summary, "avg_Root_tissue_density", "se_Root_tissue_density", "Root Tissue Density")

# Combine plots into a single figure using patchwork
combined_plot <- (plot_Root_N_concentration + plot_Mean_Root_diameter) / (plot_Specific_root_length + plot_Root_tissue_density) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# Display the combined plot
print(combined_plot)


anova_Root_N_concentration <- aov(Root_N_concentration ~ soil_type, data = OvM_GRooT_and_SG_joined_wide_main4)
summary(anova_Root_N_concentration)



anova_Mean_Root_diameter <- aov(Mean_Root_diameter ~ soil_type, data = OvM_GRooT_and_SG_joined_wide_main4)
summary(anova_Mean_Root_diameter)

anova_Specific_root_length <- aov(Specific_root_length ~ soil_type, data = OvM_GRooT_and_SG_joined_wide_main4)
summary(anova_Specific_root_length)

anova_Root_tissue_density <- aov(Root_tissue_density ~ soil_type, data = OvM_GRooT_and_SG_joined_wide_main4)
summary(anova_Root_tissue_density)



##############################################
#counts for MFAT by soil type and PFT
##############################################
names(spp_GRooT_and_SG_joined_wide)

pft_counts <- spp_GRooT_and_SG_joined_wide |>
  filter(!growthForm %in% c("fern", NA, "subshrub"),
         !is.na(soil_type) # Exclude rows with NA in soil_type
  )|>
  dplyr::group_by(growthForm, soil_type) |>
  dplyr::summarize(
    Root_N_count = sum(Root_N_count, na.rm = TRUE),
    Root_D_count = sum(Root_D_count, na.rm = TRUE),
    Root_RTD_count = sum(Root_RTD_count, na.rm = TRUE),
    SRL_count = sum(SRL_count, na.rm = TRUE),
    SOC_count = sum(soc_count, na.rm = TRUE),
    .groups = "drop" # Optional: remove grouping after summarizing
  )

print(pft_counts)




mfat_counts <- spp_GRooT_and_SG_joined_wide |>
  filter(!mycorrhizalAssociationTypeFungalRoot %in% c("EcM-AM", "NA", "NM-AM", "ErM", "OM", "uncertain", "species-specific: AM or rarely EcM-AM or AM"),
         !is.na(soil_type),
         !is.na(mycorrhizalAssociationTypeFungalRoot)
  )|>
  dplyr::group_by(mycorrhizalAssociationTypeFungalRoot, soil_type) |>
  dplyr::summarize(
    Root_N_count = sum(Root_N_count, na.rm = TRUE),
    Root_D_count = sum(Root_D_count, na.rm = TRUE),
    Root_RTD_count = sum(Root_RTD_count, na.rm = TRUE),
    SRL_count = sum(SRL_count, na.rm = TRUE),
    SOC_count = sum(soc_count, na.rm = TRUE),
    .groups = "drop" 
  )

print(mfat_counts)


#####stats for pft and mfat boxplots--- for tables s6 and s8


# ANOVA
modelNpft <- aov(Root_N_concentration ~ growthForm * soil_type, data = filtered_pft)
summary(modelNpft)


# Post hoc Tukey HSD test
tukey_Npft <- TukeyHSD(modelNpft, "growthForm:soil_type")
print(tukey_Npft)

# Use emmeans to get estimated marginal means and letters
library(emmeans)
emmeans_Npft <- emmeans(modelNpft, ~ growthForm * soil_type)
cld_Npft <- cld(emmeans_Npft, Letters = letters)


emmeans_interactionpft <- emmeans(modelNpft, ~ growthForm * soil_type)

# View the pairwise comparisons

print(emmeans_interactionpft)

library(multcompView)
library(multcomp)
# Run Tukey-adjusted comparisons
tukey_test <- pairs(emmeans_interactionpft, adjust = "tukey")

# View pairwise comparison results
print(tukey_test)

cld_results <- cld(emmeans_interactionpft, Letters = letters)
# View the letter results
print(cld_results)


# Pairwise comparison for organic vs mineral within each mycorrhizal type
pairwise_comparison <- contrast(emmeans_interactionpft, method = "pairwise", by = "growthForm")




#now for diameter
modelDpft <- aov(`Mean_Root_diameter` ~ growthForm * soil_type, data = filtered_pft)
summary(modelDpft)

# Tukey HSD post hoc
tukey_diameter <- TukeyHSD(model_diameter, "growthForm:soil_type")
print(tukey_diameter)

# Estimated marginal means and CLD for Mean_Root_diameter
emmeans_diameter <- emmeans(model_diameter, ~ growthForm * soil_type)
cld_diameter <- cld(emmeans_diameter, Letters = letters)
print(cld_diameter)

# Pairwise comparisons by soil type within each growthForm
pairwise_comparison_diameter <- contrast(emmeans_diameter, method = "pairwise", by = "growthForm")
print(pairwise_comparison_diameter)









# ANOVA for Root_tissue_density
model_RTD <- aov(Root_tissue_density ~ growthForm * soil_type, data = filtered_pft)
summary(model_RTD)

# Tukey HSD post hoc test
tukey_RTD <- TukeyHSD(model_RTD, "growthForm:soil_type")
print(tukey_RTD)

# Load emmeans package
library(emmeans)

# Calculate emmeans and compact letter display (CLD) for Root_tissue_density
emmeans_RTD <- emmeans(model_RTD, ~ growthForm * soil_type)
cld_RTD <- cld(emmeans_RTD, Letters = letters)
print(cld_RTD)



# Pairwise comparison for organic vs mineral within each growthForm for Root_tissue_density
pairwise_comparison_RTD <- contrast(emmeans_RTD, method = "pairwise", by = "growthForm")
print(pairwise_comparison_RTD)




# ANOVA for Specific_root_length
model_SRL <- aov(Specific_root_length ~ growthForm * soil_type, data = filtered_pft)
summary(model_SRL)


# Tukey HSD post hoc test
tukey_SRL <- TukeyHSD(model_SRL, "growthForm:soil_type")
print(tukey_SRL)

# Calculate emmeans and compact letter display (CLD) for Specific_root_length
emmeans_SRL <- emmeans(model_SRL, ~ growthForm * soil_type)
cld_SRL <- cld(emmeans_SRL, Letters = letters)
print(cld_SRL)


# Pairwise comparison for organic vs mineral within each growthForm for Specific_root_length
pairwise_comparison_SRL <- contrast(emmeans_SRL, method = "pairwise", by = "growthForm")
print(pairwise_comparison_SRL)






