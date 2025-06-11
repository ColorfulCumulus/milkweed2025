#### MORPHOLOGY AND VOCS NMDS ####

#### Load ####
# clear environment
rm(list = ls())

# Packages
library(tidyverse)
library(vegan)
library(GGally)
library(corrplot)
library(ggfortify)
library(tidygraph)
library(ggraph)
library(ggpubr)
library(patchwork)
library(systemfonts)
library(flextable) 
library(igraph)
library(ggforce)
library(graphlayouts)
library(reshape2)
library(grid)
library(gridExtra)
library(cowplot)
library(officer)

# Set directory
# setwd(<filepath>)


# Set color-blind palette
cb <- c("#000000", # 1 - black
        "#E69F00", # 2 - orange
        "#56B4E9", # 3 - light blue
        "#009E73", # 4 - green
        "#F0E442", # 5 - yellow
        "#0072B2", # 6 - dark blue
        "#D55E00", # 7 - red
        "#CC79A7") # 8 - pink

# Common plotting settings 
# Set colors 
cols <- c("control" = cb[4],
          "damage" = cb[7])

# Set labels
labs <- c("Control", 
          "Herbivory")

# Load data 
chems <- read.csv("data/2019_all_traits.csv", 
                     header = TRUE, 
                     check.names = FALSE)

# replace all NAs with 0
chems[is.na(chems)] <- 0

names(chems) <- make.unique(names(chems)) # so next line can work

#### Data Formatting #####
# as_tibble
chems <- as_tibble(chems)

# Declare as.factor vars and drop unused vars
chems <- chems %>% 
  mutate(across(c(pid, treatment), 
                as.factor)) %>% 
  select(!c(filename))

#extracting retention indices
ret_index <- chems[1,] %>%
  select(c(where(~is.numeric(.x) &&
                   sum(.) > 0))) 
chems <- chems[-1,]

#extracting dry weights 
dry_weights <- chems$dry_weight

# selecting morphological traits
morphotraits <- chems %>% 
  select(c(where(~is.factor(.x)),
           where(~is.numeric(.x) &&
                   sum(.) < 1100))) 

# Rename morphological vars
morpho_names <- c("Plant ID",
                  "treatment",
                  "Flower number", 
                  "Flower unopened", 
                  "Flower diameter", 
                  "Hood height", 
                  "Hood width", 
                  "Dry weight")

names(morphotraits) <- morpho_names


# removing morphological traits
chems <- chems %>% 
  select(c(where(~is.numeric(.x) &&
                   sum(.) > 1100)))

# standardizing chemicals by dry weight 
chems_w <- chems %>% 
  select(c(where(~is.numeric(.x) &&
                   sum(.) > 1100))) %>%
  mutate(./dry_weights)

# now you have two just chemical datasets:
names(chems)
names(chems_w)

# to create a full dataset again: 
  # chemicals + morphotraits, unweighted
  traits19a <- cbind(chems, morphotraits)

  # chemicals + morphotraits, weighted by dry weight
  traits19w <- cbind(chems_w, morphotraits)

# Summaries for Chemicals  ####
# get means & se values for each treatment type
treatments <- traits19w$treatment

# summarizing means by treatment
means_e<- cbind(treatments, chems_w) %>%
  group_by(treatments) %>%
  summarize_all(list(average = mean))
means_e <- means_e %>% 
  pivot_longer(!treatments, names_to = "chemical",
              values_to = "concentration")

# first getting standard deviation
std_errors <- cbind(treatments, chems_w) %>%
  group_by(treatments) %>%
  summarize_all(list(se = sd))

# dividing std devs by sqrt(n)
std_errors <- std_errors %>%
  select(where(~is.numeric(.x))) %>%
  mutate(./sqrt((traits19w %>% count(treatment))$n))
std_errors$treatment <- c("control", "damage")
std_errors <- std_errors %>% 
  pivot_longer(!treatment, names_to = "chemical",
               values_to = "concentration")

# Putting it all together
emissions_output <-bind_cols(
  names(chems_w),
  t(ret_index[1,]),
  means_e[c(1:55),3], 
  std_errors[c(1:55),3],
  means_e[c(56:110),3], 
  std_errors[c(56:110),3])

# exporting to .docx
names(emissions_output) <- c("Compound Name", "RI",  "Control Mean", "Control SE",
                             "Herbivory Mean", "Herbivory SE")
emissions_output <- as.data.frame(t(emissions_output))
colnames(emissions_output) <- emissions_output[1,]

#appending m/z names
{lookup <- c(`m/z 69, 83, 95, 109, 120, 127` = "RI 1014",
            `m/z 69, 83, 95, 109, 120, 137` = "RI 1077",
            `m/z 69, 83, 95, 109, 135, 150` ="RI 1095",
            `m/z 79, 81, 91, 93, 109` = "RI 1143",
            `m/z 79, 95, 107, 123, 135, 151, 166` = "RI 1179", 
            `m/z 91, 93, 107, 135 (RI 1219)` = "RI 1219", 
            `m/z 77, 81, 93, 121` = "RI 1233",
            `m/z 68, 79, 95, 97` = "RI 1274",
            `m/z 91, 92, 176` = "RI 1310",
            `m/z 81, 91, 95, 109, 123, 124` = "RI 1349",
            `m/z 93, 107, 123, 135 (RI 1441)` = "RI 1441",
            `m/z 93, 107, 123, 135 (RI 1623)` = "RI 1623",
            `m/z 91, 93, 107, 135 (RI 1673)` = "RI 1673",
            `m/z 93, 107, 121, 135` = "RI 1678",
            `m/z 77, 91, 123, 151` = "RI 1694",
            `m/z 91, 97, 184` = "RI 1739",
            `m/z 91, 109, 123, 135` = "RI 1742",
            `m/z 77, 91, 151, 77, 199` = "RI 1751",
            `m/z 91, 109, 137, 151` = "RI 1758",
            `m/z 91, 109, 151, 199` = "RI 1768",
            `m/z 81, 91, 107, 121` = "RI 1851",
            `m/z 91, 121, 151, 242` = "RI 1884",
            `m/z 91, 95, 121, 184` = "RI 1895",
            `m/z 135, 191, 226` = "RI 1915",
            `m/z 91, 95, 199, 242` = "RI 1955")
}

# rename
emissions_output <- rename(emissions_output, 
                    all_of(lookup))

emissions_output <- as.data.frame(t(emissions_output))

emissions_output[,1] <- rownames(emissions_output)
rownames(emissions_output) <- NULL

for(n in 2:6){
emissions_output[,n] <- as.numeric(emissions_output[,n])
}
# exporting happens at line 380. 

####  NMDS #####

# Standardize variables
# 2019
# Hellinger's standardization
voctraits19 <- chems_w %>%
  decostand(method = "hellinger")

# Use range for morphology
morphotraits19 <- morphotraits %>% 
  select(where(~is.numeric(.x))) %>% 
  decostand(method = "range")

# Replace raw VOC's
s_traits19 <- 
  bind_cols(voctraits19, morphotraits19)

# Get distance matrix
dist19 <- s_traits19 %>% 
  select(where(~is.numeric(.x))) %>% 
  vegdist(method = "bray")


# Get factors
fact_traits19 <- morphotraits %>% 
  select(c(`Plant ID`, treatment))

# Run NMDS
nmds19 <- metaMDS(dist19, 
                  distance = "bray",
                  autotransform = FALSE,
                  k = 3)


###### NMDS Plot ########
# Get scores
nmds19_scores <- as_tibble(scores(nmds19))

# Add treatment and pid
nmds19_scores <- nmds19_scores %>% 
  mutate(pid = fact_traits19$`Plant ID`,
         treatment = fact_traits19$treatment)

# Format data for hull plots 
# Hull values - control
grp19_con <- 
  nmds19_scores[nmds19_scores$treatment == "control", ][chull(nmds19_scores[
    nmds19_scores$treatment == "control", c("NMDS1", "NMDS2")]), ] 

# Hull values - damage
grp19_dam <- 
  nmds19_scores[nmds19_scores$treatment == "damage", ][chull(nmds19_scores[
    nmds19_scores$treatment == "damage", c("NMDS1", "NMDS2")]), ]  

# Combine
hull19_scores <- rbind(grp19_con, 
                       grp19_dam)  

# Plot
plot_nmds19 <- ggplot(nmds19_scores) + 
  geom_point(aes(x = NMDS1, # points
                 y = NMDS2, 
                 color = treatment, 
                 shape = treatment), 
             size = 3) + 
  geom_polygon(data = hull19_scores, # shaded area
               aes(x = NMDS1,
                   y = NMDS2, 
                   fill = treatment, 
                   group = treatment),
               alpha = 0.3) + 
  scale_color_manual(name = "Treatment",
                     values = cols, 
                     labels = labs) +
  scale_fill_manual(name = "Treatment",
                    values = cols,
                    labels = labs) +
  scale_shape_manual(name = "Treatment",
                     values = c(18,19),
                     labels = labs) +
  geom_text(aes(label = paste("Stress = ", 
                              round(nmds19$stress, 
                                    digits = 3)), 
                x = 0.18, 
                y = 0.22)) + 
  theme_classic(base_size = 18) +
  theme(legend.position = "bottom")

# Save
ggsave(plot = plot_nmds19,
       filename = "figures/Figure_2.pdf",
       width = 12,
       height = 10, 
       units = "cm")

##### Distance Tests #####
(permanova19 <- adonis2(dist19 ~ treatment,
                       data = fact_traits19, 
                       permutations = 10000, 
                       method = "bray"))

(anosim19 <-  anosim(dist19, 
                    fact_traits19$treatment, 
                    permutations = 10000))

(disptest19 <- anova(betadisper(dist19, 
                               fact_traits19$treatment)))

simper19 <- simper(comm = select(s_traits19, 
                                 where(~is.numeric(.x))), 
                   group = fact_traits19$treatment,
                   permutations = 10000)
summary(simper19)

##### SIMPER Table ####
# Store output
simper_output <- summary(simper19)

# Keep only results
simper_output <- simper_output$damage_control

# Row names as a column
simper_output$`Trait/Compound Name` <- rownames(simper_output)
rownames(simper_output) <- NULL

# Make trait/compound names first column
simper_output <- simper_output %>% 
  relocate(`Trait/Compound Name`, 
           .before = average)

names(simper_output)[2:ncol(simper_output)] <- 
  str_to_sentence(names(simper_output[2:ncol(simper_output)]))

names(simper_output)[5] <- "Av-cont"
names(simper_output)[6] <- "Av-herb"


# Convert to flextable, export

#specifying page properties
page_properties <- prop_section(
  page_size = page_size(
    orient = "landscape",
    width = 8.5, height = 11
  ),
  page_margins = page_mar()
)

# custom code to make autofit work
fit_flex_to_page <- function(ft, pgwidth = x){
  
  ft_out <- ft %>% autofit()
  
  ft_out <- width(ft_out, width = dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
  return(ft_out)
}

#exporting SIMPER TABLE
simper_ft <- as_flextable(simper_output, 
                          max_row = nrow(simper_output), 
                          show_coltype = F) %>% 
  colformat_double(digits = 3) %>% 
  theme_vanilla() %>% 
  bold(part = "body", # bold only significant 
       ~P < 0.05) %>%
  fit_flex_to_page(pgwidth = 9) %>%
  fontsize(size = 12) %>%
  fontsize(size = 12, part = "header") %>%
  save_as_docx(path = "tables/Table_S2.docx", pr_section = page_properties) #save


# using SIMPER values to bold chems for Table S1
pvalsS1 <- data.frame(simper_output$`Trait/Compound Name`,
                      simper_output$P)
pvalsS1 <- pvalsS1[-(1:6),] #removing morphotraits
pvalsS1<- pvalsS1[match(names(chems_w), pvalsS1[,1]), ]

# creating base flex table
emissions_ft <- as_flextable(emissions_output, 
                             max_row = nrow(emissions_output), 
                             show_coltype = F) %>% 
  fontsize(size =12)%>%
  fontsize(size = 12, part = "header") %>%
  colformat_double(digits = 0, big.mark = "") %>% 
  colformat_int(big.mark = "") %>%
  theme_vanilla() %>% 
  bold(part = "body", # bold only significant 
       ~pvalsS1[,2] < 0.05) %>%
  fit_flex_to_page(pgwidth = 9) 


# exporting emissions table
emissions_ft %>%
 save_as_docx(path = "tables/Table_S1.docx", pr_section = page_properties) 

# creating a graph

#first, selecting top 10 compounds only
tops <- simper_output[c(7:20),]

#then, specifying what isn't significant in simper
tops$sig <- ifelse(tops$P<0.05, "Y", "N")
tops<- tops %>%
  arrange(desc(sig))
tops <- tops[,1]

#then, grab the mean ± se values
simper_tops <- means_e
simper_tops$se <- std_errors$concentration

#removing "average"
simper_tops$chemical <- gsub("_average", "", simper_tops$chemical)

#then, select only the top 10 compounds
cc <- data.frame(matrix(0, nrow = 2*length(tops),
                       ncol = ncol(simper_tops)))
colnames(cc) <- colnames(simper_tops)

for(k in 1:length(tops)){
bb <- simper_tops %>% 
        filter(chemical == tops[k])

for(j in 0:(nrow(bb)-1)){
cc[((2*k-1)+j),] <- bb[(j+1),]
}
}

# for some reason control = 1, treatment = 2. fixing
cc$treatments <- gsub(1, "control", cc$treatments)
cc$treatments <- gsub(2, "treatment", cc$treatments)

## plot
plot_list <- list()

for(i in 1:length(tops)){
plot_list[[i]] <- local({
  i <- i
  ggplot(cc[(2*i-1):(2*i),], aes(x=treatments, y=concentration, group=1)) +
      geom_errorbar(width=.1, 
      aes(ymin=concentration-se, 
          ymax=concentration+se)) +
      geom_point(shape=21, size=5, color = cols, fill= cols) +
      labs(x = substr(tops[i], 1, 21),
           y ="") +
      #labs(y = expression(paste("VOC Emissions"~ "("*ng%.% h^-1*")"))) +
      theme_classic() +
      theme(axis.title = element_text(size = 12, face = "bold")) 
    })
}      
# stitching it all together
allsigs <- ggarrange(plot_list[[1]], plot_list[[2]], plot_list[[3]], 
                     plot_list[[4]], plot_list[[5]], plot_list[[6]], 
                     plot_list[[7]], plot_list[[8]], plot_list[[9]], 
                     plot_list[[10]], plot_list[[11]], plot_list[[12]],
                     plot_list[[13]], plot_list[[14]], 
                     labels = c("**", "***", "*", "**", "***", "**", "***",
                                "", "", "", "", "", "", ""),
                     label.x = 0.5,
                     font.label = list(size = 20, color = "black", 
                                       face = "bold", family = NULL),
                     ncol = 5, nrow = 3)
#exporting
ggsave(plot = allsigs,
       filename = "figures/Figure_3.pdf",
       width = 25,
       height = 15, 
       units = "cm")



#### Correlation Structure ####
# Stratify by treatment
# Control
control_traits19 <- traits19w %>% 
  filter(treatment == "control") %>% 
  select(!c(`Plant ID`, 
            treatment, 
            `Flower unopened`))

# Treatment
treatment_traits19 <- traits19w %>% 
  filter(treatment == "damage") %>% 
  select(!c(`Plant ID`,
            treatment, 
            `Flower unopened`))

##### Correlation Tests ####
# Compute correlation matrices for control and treatment
control_cor_matrix19 <- cor(control_traits19, 
                            use = "pairwise.complete.obs")
treatment_cor_matrix19 <- cor(treatment_traits19, 
                              use = "pairwise.complete.obs")

# get p-values for corrs
control_pmat19 <- cor.mtest(control_traits19)$p
treatment_pmat19 <- cor.mtest(treatment_traits19)$p

# Get difference
cor_diff_matrix19 <- abs(treatment_cor_matrix19) - abs(control_cor_matrix19)

# Clone matrix
cor_diff_dir <- cor_diff_matrix19

# Abbreviate row names to the first 10 letters, change greek letters
corr_rownames <- substr(rownames(cor_diff_matrix19), 1, 12)
corr_rownames <- gsub("α-", "alpha-", corr_rownames)
corr_rownames <- gsub("β-", "beta-", corr_rownames)
rownames(cor_diff_dir) <- corr_rownames

# Abbreviate column names to the first 10 letters, change greekletters
corr_colnames <- substr(colnames(cor_diff_matrix19), 1, 12)
corr_colnames <- gsub("α-", "alpha-", corr_colnames)
corr_colnames <- gsub("β-", "beta-", corr_colnames)
colnames(cor_diff_dir)  <- corr_colnames


for(i in 1:nrow(control_cor_matrix19)){
  for(j in 1:ncol(control_cor_matrix19)){
    
    # 1 if r increases, -1 if r decreases
    if(treatment_cor_matrix19[i,j] > control_cor_matrix19[i,j]){
      cor_diff_dir[i,j] <- 1
    } else {
      cor_diff_dir[i,j] <- -1
    }
     
  }
}


# Change to 
# cor_diff_matrix19 <- control_cor_matrix19 - treatment_cor_matrix19

# Sample sizes
control_n19 <- nrow(control_traits19)
treatment_n19 <- nrow(treatment_traits19)

# Calculate Z stats
control_z19 <- 0.5 * log((1 + control_cor_matrix19) / 
                           (1 - control_cor_matrix19))
treatment_z19 <- 0.5 * log((1 + treatment_cor_matrix19) / 
                             (1 - treatment_cor_matrix19))

# Fishers' Z-test
test_statistic19 <- (control_z19 - treatment_z19) / 
  sqrt(1 / (control_n19 - 3) + 1 / (treatment_n19 - 3))

# Get p-values
p_value_matrix19 <- 2 * (1 - pnorm(abs(test_statistic19)))
colnames(p_value_matrix19)  <- corr_colnames
rownames(p_value_matrix19) <- corr_rownames

# Print the matrices
# print(control_cor_matrix19)
# print(treatment_cor_matrix19)
# print(cor_diff_matrix19)
# print(cor_diff_dir)
# print(p_value_matrix19)


# Scale find directions of changes 
# scaled_diff_matrix19 <- 2 * (cor_diff_matrix19 - min(cor_diff_matrix19)) / (max(cor_diff_matrix19) - min(cor_diff_matrix19)) - 1



##### Heatmaps Plots ####
# ### Control Traits
# jpeg(file = "figures/control_traits_heatmap.jpeg", 
#      height = 2000, 
#      width = 2000)
# 
# 
# # Heatmap
# corrplot(control_cor_matrix19 , 
#          p.mat = control_pmat19, 
#          sig.level = alpha, 
#          insig = "pch", 
#          pch.col = "red",
#          tl.col = "black",
#          method = "color", 
#          addCoef.col = "black", 
#          diag = FALSE, 
#          type = "upper")
# 
# dev.off()

### Damage Traits
# jpeg(file = "figures/damage_traits_heatmap.jpeg", 
#      height = 2000, 
#      width = 2000)
# 
# # Heatmap
# corrplot(treatment_cor_matrix19, 
#          p.mat = treatment_pmat19, 
#          sig.level = alpha, 
#          insig = "pch", 
#          pch.col = 'red',
#          tl.col = "black",
#          method = "color", 
#          addCoef.col = "black", 
#          diag = FALSE, 
#          type = "upper")
# 
# dev.off()

### Differences 
pdf(file = "figures/Figure_S3.pdf",
     height = unit(12, "cm"),
     width = unit(12, "cm"))

# significance value
alpha <- 0.05

# Heatmap
 corrplot(cor_diff_dir,
         p.mat = p_value_matrix19,
         sig.level = alpha,
         col = c(cb[7], cb[6]),
         method = "square",
         is.corr = FALSE,
         insig = "blank",
         tl.col = "black",
         diag = FALSE,
         type = "upper",
         tl.cex = .8,
         cl.cex = 1)

# End
dev.off()

###### Combined heatmap #### 
# Create new matrix that includes control trait correlations in lower triangle and,
# damage trait correlation in upper.
# begin w/ control corrs in lower triangle
contrast_cor_matrix <- control_cor_matrix19
contrast_cor_matrix[upper.tri(contrast_cor_matrix)] <- 0

# Insert damage
contrast_cor_matrix[upper.tri(contrast_cor_matrix)] <- 
  treatment_cor_matrix19[upper.tri(treatment_cor_matrix19)]

colnames(contrast_cor_matrix) <- corr_colnames
rownames(contrast_cor_matrix) <- corr_rownames

# Create new matrix that includes p values for lower and upper triangles
# begin w/ control corrs in lower triangle
contrast_p_matrix <- control_pmat19

# add damage p-values
contrast_p_matrix[upper.tri(contrast_p_matrix)] <- 
  treatment_pmat19[upper.tri(treatment_pmat19)]

# assign names - not going to be on final plot
  # but removes annoying error in ggplot code
colnames(contrast_p_matrix) <- corr_colnames
rownames(contrast_p_matrix) <- corr_rownames

# Create Heatmap
pdf(file = "figures/Figure_4.pdf", 
     height = unit(12, "cm"), 
     width = unit(12, "cm"))
par(xpd = TRUE)

# significance value
alpha <- 0.05

corrplot(contrast_cor_matrix, 
         p.mat = contrast_p_matrix, 
         sig.level = alpha, 
         insig = "blank", 
         tl.col = "black",
         tl.cex = .8,
         cl.cex = 1,
         cl.pos = "b",
         method = "square", 
         order = "original", 
         mar = c(0,4.5,4.5,0),
         diag = FALSE)

# Add triangle labels
# Damage 
text(x = 28, 
     y = 90, 
     labels = "Trait-Pair Correlations \n(Damage)", 
     cex = 1.5, 
     col = "black")

# Control
text(x = -30, 
     y = 30, 
     labels = "Trait-Pair Correlations \n(Control)", 
     cex = 1.5, 
     col = "black",
     srt = 90)

# Add separating line through diagonal
segments(60, 1, 1, 60, 
         col = cb[7], 
         lwd = 5)

# Segments highlighting morphological traits
segments(-9, 5, -9, 1, 
         col = cb[8], 
         lwd = 5)

segments(56, 70, 60, 70, 
         col = cb[8], 
         lwd = 5)

# Segments highlighting morphological traits
segments(0.25, 5, 0.25, 1, 
         col = cb[8], 
         lwd = 5)
segments(56, 60.75, 60, 60.75, 
         col = cb[8], 
         lwd = 5)

# End
dev.off()

##### Histogram Plots ###### 
# Data Format: Create new corr. matrix, w/ only significant correlations
# 
# Control #
control_graph19 <- matrix(data = NA,
                          ncol = ncol(control_cor_matrix19),
                          nrow = nrow(control_cor_matrix19),
                          dimnames = dimnames(control_cor_matrix19))

# Create binary graph
control_binary19 <- control_graph19

# Only keep significant values
for(i in 1:nrow(control_cor_matrix19)){
  for(j in 1:ncol(control_cor_matrix19)){
    if(!is.nan(control_pmat19[i,j]) & 
       control_pmat19[i,j] <= 0.05){
      control_graph19[i,j] <- control_cor_matrix19[i,j]
      control_binary19[i,j] <- 1
    } else {
      control_graph19[i,j] <- 0
      control_binary19[i,j] <- 0
    }
  }
}


# Damage
treatment_graph19 <- matrix(data = NA,
                            ncol = ncol(treatment_cor_matrix19),
                            nrow = nrow(treatment_cor_matrix19),
                            dimnames = dimnames(treatment_cor_matrix19))

# Create binary graph
treatment_binary19 <- treatment_graph19

# Only keep significant values
for(i in 1:nrow(treatment_cor_matrix19)){
  for(j in 1:ncol(treatment_cor_matrix19)){
    if(!is.nan(treatment_pmat19[i,j]) & 
       treatment_pmat19[i,j] <= 0.05){
      treatment_graph19[i,j] <- treatment_cor_matrix19[i,j]
      treatment_binary19[i,j] <- 1
    } else {
      treatment_graph19[i,j] <- 0
      treatment_binary19[i,j] <- 0
    }
  }
}

# Create long form for histograms
control_corrs <- as_tibble(melt(control_cor_matrix19))
treatment_corrs <- as_tibble(melt(treatment_cor_matrix19))

#get mean Pearson's R for each. Use this for x-interctepts of lines in plot
rmean_control<- mean(control_corrs$value)
rmean_treatment<- mean(treatment_corrs$value)

# load legend for plot
#common_legend <- readRDS("figures/common_legend.RDS")

# Create Histogram
histo <- ggplot() +
  geom_histogram(data = control_corrs, 
                 aes(x = value), 
                 fill = cb[4], 
                 alpha = 0.5, 
                 bins = 40) + 
  geom_histogram(data = treatment_corrs,
                 aes(x = value), 
                 fill = cb[7],
                 alpha = 0.5, 
                 bins = 40) +
  geom_vline(xintercept = mean(control_corrs$value), 
             color = cb[4],
             linetype = 2, 
             size = 1.5) +
  geom_vline(xintercept = mean(treatment_corrs$value), 
             color = cb[7],
             linetype = 2, 
             size = 1.5) +
  labs(x = expression(paste("Trait-Pair Correlations",
                            "\n (Pearson's ", 
                            italic("r"),
                            ")")), 
       y = "Frequency") +
  lims(x = c(-1,1)) +
  theme_classic()

#histo <- histo / common_legend + 
 # plot_layout(heights = unit(c(7,1), 
  #                           c('cm', 'null')))

# Save
ggsave(plot = histo,
       filename = "figures/Figure_S2.pdf",
       width = 12,
       height = 10, 
       units = "cm")


#### Emission Rates #####
# Subset only VOCs
vocs <- cbind(treatments, chems_w) %>%
  mutate(total = rowSums(pick(where(is.numeric))))


###### Mann Whitney Test ####
# Subset totals
control_totals <- vocs %>% 
  filter(treatments == "control") %>% 
  select(total) %>% 
  pull() # as vector
herbivory_totals <- vocs %>% 
  filter(treatments == "damage") %>% 
  select(total) %>% 
  pull() # as vector


# Mann-Whitney test
wilcox.test(control_totals, # control
            herbivory_totals) # damage

# Calculate means
mean(control_totals) # 58357.51
mean(herbivory_totals) # 71633.8

###### Emissions Boxplot #####
emissions_plot <- ggplot(data = vocs) +
  geom_boxplot(aes(x = treatments, 
                   y = total, 
                   fill = treatments),
               color = cb[1]) + 
  theme_classic(base_size = 15) + 
  labs(x = "Treatments", 
       y = expression(paste("VOC Emissions"~ "("*ng%.% g^-1 %.% h^-1*")"))) +
  scale_fill_manual(name = "Treatments",
                     values = cols, 
                     labels = labs) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("Control", 
                              "Herbivory")) +
  annotate("text",
           x = 1.5, y = 90000, 
           label = "NS", 
           size = 8)


# Save
ggsave(plot = emissions_plot,
       filename = "figures/Figure_S1.pdf",
       width = 8,
       height = 10, 
       units = "cm")

#### Changes in Floral Morphology #####

###### Mann Whitney Tests ####
# Mann-Whitney test
# Flower diameter
wilcox.test(`Flower diameter` ~ treatment, 
            data = morphotraits) # significant

# Hood width
wilcox.test(`Hood width` ~ treatment, 
            data = morphotraits) # significant

# Hood height
wilcox.test(`Hood height` ~ treatment, 
            data = morphotraits) # significant 

# Flower number
wilcox.test(`Flower number` ~ treatment, 
            data = morphotraits, 
            exact = F) #  exact = F because of ties. no effect

# Dry weight
wilcox.test(`Dry weight` ~ treatment, 
            data = morphotraits, 
            exact = F) # exact = F because of ties. no effect


###### Morphology Plots ####
# Set repeating plotting elements
labels <- scale_x_discrete(labels = c("Control", "Herbivory"))
theme_size <- theme_classic(base_size = 16)
modified_legend <- scale_color_manual(name = "Treatment",
                                      values = cols, 
                                      labels = labs) 
x_lab <- xlab("Treatment")
no_legend <- theme(legend.position = "none")

# Flower diameter
diam_plot <- ggplot(morphotraits, 
                     aes(x = treatment)) +
  geom_jitter(aes(y = `Flower diameter`, 
                  color = treatment), 
              height = 0, # jitter mods
              width = .15) +
  labs(y = "Flower Diameter (mm)", 
       tag = "a") +
  labels + 
  theme_size +
  modified_legend +
  x_lab +
  no_legend +
  annotate("text", # significance stars
           x = 1.5, 
           y = 9.3,
           label = "*", 
           size = 8)


# Hood width
width_plot <- ggplot(morphotraits, 
                     aes(x = treatment)) +
  geom_jitter(aes(y = `Hood width`, 
                  color = treatment), 
              height = 0, # jitter mods
              width = .15) +
  labs(y = "Hood Width (mm)", 
       tag = "b") + 
  labels + 
  theme_size +
  modified_legend +
  x_lab +
  no_legend +
  annotate("text", 
           x = 1.5, 
           y = 2.8,
           label = "***", 
           size = 8)

# Hood height
height_plot <- ggplot(morphotraits, 
                    aes(x = treatment)) +
  geom_jitter(aes(y = `Hood height`, 
                   color = treatment), 
               height = 0, # jitter mods
               width = .15) +
  labs(y = "Hood Height (mm)", 
       tag = "c") +
  labels + 
  theme_size +
  modified_legend +
  x_lab +
  no_legend +
  annotate("text", 
           x = 1.5,
           y = 6.37,
           label = "***", 
           size = 8)

# Flower Number
count_plot <- ggplot(morphotraits, 
                     aes(x = treatment)) +
  geom_jitter(aes(y = `Flower number`, 
                  color = treatment), 
              height = 0, # jitter mods
              width = .15) + 
  labs(y = "Flower No. \nper Inflorescence", 
       tag = "d") +
  labels + 
  theme_size +
  modified_legend +
  x_lab +
  no_legend

# Dry weight
weight_plot <- ggplot(morphotraits, 
                      aes(x = treatment)) +
  geom_jitter(aes(y = `Dry weight`, 
                  color = treatment), 
              height = 0, # jitter mods
              width = .15) +
  labs(y = "Inflorescence \nDry Weight (g)", 
       tag = "e") +
  labels + 
  theme_size +
  modified_legend +
  x_lab +
  no_legend


morpho_plots <- diam_plot + width_plot + height_plot + # row 1
  count_plot + weight_plot + guide_area() + # row 2
  plot_layout(guides = "collect") &
  theme(legend.position = "top") & 
  labs(x = "Treatment") & 
  theme_classic()

# Save 
ggsave(plot = morpho_plots,
       filename = "figures/Figure_1.pdf",
       width = 16,
       height = 9, 
       units = "cm")

# Unused code ####
if(FALSE){
  # remove columns w/ all zeros 
  # this line of code asks: 
  #(1) is the data a factor?
  #(2) is the data numeric, and does it have less than 15 zeros?
  # if YES to 1 or 2, the data is kept.
  
  chems <- chems %>% 
    select(c(where(~is.factor(.x)), 
             where(~is.numeric(.x) && 
                     sum(. == 0) < 15))) 
# Update names
oldnames <- names(traits19)
for(i in 4:ncol(traits19)){
  if(colnames(traits19)[i] %in% updated_names$old_name){
    # check matches
    print(paste("name found", 
                match(colnames(traits19)[i], 
                      updated_names$old_name))) 
    
    # insert new name on 
    pos <- match(colnames(traits19)[i], 
                 updated_names$old_name)
    colnames(traits19)[i] <- updated_names$new_name[pos]
  }
}

# make column names.unique
names(traits19) <- make.unique(names(traits19))

# load new name and info
updated_names <- read.csv("data/updated_names.csv", 
                          header = TRUE)

# Rename chemicals
if(FALSE){ # change to FALSE if names are updated.
  # even if TRUE, will not change names that do not exist
  lookup <- c(`Propylbenzene` = "Benzene, propyl-",
              `B-Myrcene` = "beta-Myrcene" ,
              `m/z 69, 83, 95, 109, 120, 127` = "Unknown 1014",
              `trans-B-Ocimene` = "(E)-beta Ocimene" ,
              `1,2,3,4a,8a-hexahydronaphthalene` = "Naphthalene, 1,2,3,4,4a,8a-hexahydro-",
              `m/z 69, 83, 95, 109, 120, 137` = "Unknown 1077",
              `2,5-dimethyl-1,6-heptadiene` = "1,6-Heptadiene, 2,5-dimethyl-", 
              `m/z 69, 83, 95, 109, 135, 150` = "Unknown 1095",
              `m/z 79, 81, 91, 93, 109` = "Unknown 1143",
              `Pentylbenzene` = "Benzene, pentyl-",
              `m/z 79, 95, 107, 123, 135, 151, 166` = "Unknown 1179",
              `m/z 91, 93, 107, 135 (RT 1219)` = "Unknown 1219", 
              `m/z 77, 81, 93, 121` = "Unknown 1233",
              `m/z 68, 79, 95, 97` = "Unknown 1274",
              `Heptylbenzene (1)` = "Benzene, heptyl-", 
              `Heptylbenzene (2)` =  "Benzene, heptyl-.1", 
              `m/z 81, 91, 95, 109, 123, 124` = "Unknown 1349",
              `Heptylbenzene (3)` = "Benzene, heptyl-.2",
              `Benzyl isovalerate` =
                "Butanoic acid, 3-methyl-, phenylmethyl ester",
              `trans-Isoeugenol` = "Phenol, 2-methoxy-4-(1-propenyl)-, (E)-",
              `m/z 93, 107, 123, 135 (RT 1441)` = "Unknown 1441",
              `cis-Isoeugenol` = "Isoeugenol",
              `a-Farnesene` = "alpha-Farnesene" ,
              `m/z 93, 107, 123, 135 (RT 1623)` = "Unknown 1623",
              `m/z 91, 93, 107, 135 (RT 1673)` = "Unknown 1673",
              `m/z 93, 107, 121, 135` ="Unknown 1678",
              `m/z 77, 91, 123, 151` = "Unknown 1694",
              `Methyl phenylacetate` =
                "Phenylacetic acid, (-)-menthyl ester",
              `m/z 91, 97, 184` = "Unknown 1739",
              `m/z 91, 109, 123, 135` = "Unknown 1742",
              `m/z 77, 91, 151, 77, 199` = "Unknown 1751",
              `m/z 91, 109, 137, 151` = "Unknown 1758",
              `m/z 91, 109, 151, 199` = "Unknown 1768",
              `m/z 81, 91, 107, 121` = "Unknown 1851",
              `m/z 91, 121, 151, 242` = "Unknown 1884",
              `m/z 91, 95, 121, 184` = "Unknown 1895",
              `m/z 135, 191, 226` = "Unknown 1915",
              `m/z 91, 95, 199, 242` = "Unknown 1955")
  
  # just chemicals, unweighted
  chems <- rename(chems, 
                  all_of(lookup))
  # just chemicals, weighted
  chems_w <- rename(chems_w,
                    all_of(lookup))
}
}

