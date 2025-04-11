library(SingleCellExperiment)
library(ggplot2)
library(tidyverse)
library(Seurat)

#### Function: change the first letter to capital
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
#### Load All MOB Slices in Matrix 
countMat = readRDS("~/data/Collaboration/AlexF/MOB/AllSection_NormCount.rds")


######## IRIS domain and 3D Coordiante
##### IRIS Results and 3D Coordinate
spatial_countMat_list = readRDS("~/data/Collaboration/AlexF/MOB/spatial_countMat_list.RDS")
spatial_location_list = readRDS("~/data/Collaboration/AlexF/MOB/spatial_location_list.RDS")
# spatial_countMat_list = spatial_countMat_list[106:126]
# spatial_location_list = spatial_location_list[106:126]
spatial_countMat_list = spatial_countMat_list[37:146]
spatial_location_list = spatial_location_list[37:146]

numCluster = 10
IRIS_Domain = read.csv(paste0("~/data/Collaboration/AlexF/MOB/IRIS_Leiden_classifyDV_Analysis/DomainLabel/Fleischmann_leiden_clusters.csv.gz"))
coord3d = read.csv("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE/symmetrized_spots_correct.csv")

##### Combine IRIS Results and Visium 3D Coordinate and D-V layers
coord3d$slice = unlist(lapply(coord3d$folder,function(x){strsplit(x,split = "_")[[1]][1]}))
coord3d$slice =  unlist(lapply(coord3d$slice,function(x){
  temp = strsplit(x,split = "-")[[1]]
  return(paste0("Slide",temp[2],"_",temp[3]))}))

### for Leiden Label File
IRIS_Domain.slide =  unlist(lapply(IRIS_Domain$zone,function(x){strsplit(x,split = "_ob")[[1]][1]}))
IRIS_Domain.section =  unlist(lapply(IRIS_Domain$zone,function(x){strsplit(x,split = "_ob")[[1]][2]}))
IRIS_Domain$Slice = paste0(IRIS_Domain.slide,"_",IRIS_Domain.section)

coord3d  = coord3d %>%
  filter(slice %in% names(spatial_countMat_list))

section_slide = coord3d %>%
  dplyr::select(slice,section) %>%
  unique() %>%
  dplyr::filter(slice %in% names(spatial_countMat_list)) %>%
  dplyr::arrange(section)

section_s = section_slide$section
visium_s = section_slide$slice

coord3d = coord3d %>%
  filter(section %in% section_s)
IRIS_Domain = IRIS_Domain %>%
  filter(Slice %in% visium_s)

# IRIS_Domain$slice_barcode = paste0(IRIS_Domain$Slice,"_",IRIS_Domain$spotName)
IRIS_Domain.spotName = paste0(unlist(lapply(IRIS_Domain$barcode,function(x){strsplit(x,split = "-")[[1]][1]})),"-",
                              (unlist(lapply(IRIS_Domain$barcode,function(x){strsplit(x,split = "-")[[1]][2]}))) )
IRIS_Domain$spotName = IRIS_Domain.spotName
IRIS_Domain$slice_barcode = paste0(IRIS_Domain$Slice,"_",IRIS_Domain$spotName)

coord3d$slice_barcode = paste0(coord3d$slice,"_",coord3d$barcode,"-1")
# IRIS_Domain$iris_domain = paste0("IRIS Domain",IRIS_Domain$IRIS_domain)
IRIS_Domain$iris_domain = paste0("Leiden Domain",IRIS_Domain$leiden_1_2)

coord3d_IRIS = merge(coord3d,
                     IRIS_Domain,
                     by.x = c("slice_barcode"),
                     by.y = c("slice_barcode"))




######## Find Gene Plane Deep granule cell
# countMat = readRDS("~/data/Collaboration/AlexF/MOB/AllSection_NormCount.rds")
slice_list = visium_s
domain_algorithm = "Leiden_110Section"
select_domain = c("Leiden Domain1")
select_domain_name = "Domain1"
gene_list = rownames(countMat)
# gene_list = unique(granularMarkerTop200$gene)
# gene_list = intersect(gene_list,rownames(countMat))

DEA_rotateAxis = NULL
coord3d_IRIS_slice_domain = coord3d_IRIS %>%
  filter(iris_domain %in% select_domain)

coord3d_IRIS_slice_domain_left = coord3d_IRIS_slice_domain %>%
  filter(x.x < 0)
coord3d_IRIS_slice_domain_right = coord3d_IRIS_slice_domain %>%
  filter(x.x >= 0)

glomer_axis_angle_left = round(atan(6.7446/1.3679) * 180 / pi,0)
glomer_axis_angle_right = round(atan(-6.7446/1.3679) * 180 / pi,0) + 180

coord3d_IRIS_slice_domain_left$x.center = coord3d_IRIS_slice_domain_left$x.x - mean(coord3d_IRIS_slice_domain_left$x.x)
coord3d_IRIS_slice_domain_left$y.center = coord3d_IRIS_slice_domain_left$y.x - mean(coord3d_IRIS_slice_domain_left$y.x)
coord3d_IRIS_slice_domain_right$x.center = coord3d_IRIS_slice_domain_right$x.x - mean(coord3d_IRIS_slice_domain_right$x.x)
coord3d_IRIS_slice_domain_right$y.center = coord3d_IRIS_slice_domain_right$y.x - mean(coord3d_IRIS_slice_domain_right$y.x)
coord3d_IRIS_slice_domain_left$x.center.rotate = coord3d_IRIS_slice_domain_left$x.center
coord3d_IRIS_slice_domain_left$y.center.rotate = coord3d_IRIS_slice_domain_left$y.center

ggplot(coord3d_IRIS_slice_domain_left,aes(x = x.center,y = y.center)) + geom_point()
ggplot(coord3d_IRIS_slice_domain_left,aes(x = x.center.rotate,y = y.center.rotate)) + geom_point()

coord3d_IRIS_slice_domain_right$x.center.rotate = coord3d_IRIS_slice_domain_right$x.center
coord3d_IRIS_slice_domain_right$y.center.rotate = coord3d_IRIS_slice_domain_right$y.center

ggplot(coord3d_IRIS_slice_domain_right,aes(x = x.center,y = y.center)) + geom_point()
ggplot(coord3d_IRIS_slice_domain_right,aes(x = x.center.rotate,y = y.center.rotate)) + geom_point()

coord3d_IRIS_slice_domain_right$x.center.flipped = -coord3d_IRIS_slice_domain_right$x.center
coord3d_IRIS_slice_domain_right$y.center.flipped = coord3d_IRIS_slice_domain_right$y.center


ggplot(coord3d_IRIS_slice_domain_right,aes(x = x.center.flipped,y = y.center.flipped)) + geom_point()
ggplot(coord3d_IRIS_slice_domain_left,aes(x = x.center,y = y.center)) + geom_point()


coord3d_IRIS_slice_domain_left$x.center.flipped = coord3d_IRIS_slice_domain_left$x.center
coord3d_IRIS_slice_domain_left$y.center.flipped = coord3d_IRIS_slice_domain_left$y.center
coord3d_IRIS_slice_domain_left$side = "left"
coord3d_IRIS_slice_domain_right$side = "right"
coord3d_IRIS_slice_domain_reflected = rbind(coord3d_IRIS_slice_domain_left,coord3d_IRIS_slice_domain_right)
  
glomer_axis_rad_left  <- glomer_axis_angle_left  * pi / 180
slope1 <- sin(glomer_axis_rad_left) / cos(glomer_axis_rad_left)
p = ggplot(coord3d_IRIS_slice_domain_reflected,aes(x = x.center.flipped,y = y.center.flipped,color = side)) + 
  geom_point() +
  geom_abline(color = "blue",slope = slope1,size = 2,linetype = "dashed") +
  labs(color = "MOB")

pdf(paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit/Left_Right_Overlaid_Map/",
           domain_algorithm,"_",
           select_domain_name,".pdf"),width = 4,height = 4)
p
dev.off()


#### Angel of axis
coord3d_IRIS_slice_domain_left = coord3d_IRIS_slice_domain_reflected
coord3d_IRIS_slice_domain_left$x.center.rotate = coord3d_IRIS_slice_domain_left$x.center.flipped
coord3d_IRIS_slice_domain_left$y.center.rotate = coord3d_IRIS_slice_domain_left$y.center.flipped


num_angles <- 10
left_angles = seq( glomer_axis_angle_left,glomer_axis_angle_left - 180, length.out = num_angles + 1)[-(num_angles + 1)] # Exclude the last point to avoid duplicating 0 degrees
right_angles = seq( glomer_axis_angle_right, glomer_axis_angle_right - 180,length.out = num_angles + 1)[-(num_angles + 1)] # Exclude the last point to avoid duplicating 0 degrees
left_angles_nomial = seq(0,180,length.out = num_angles + 1)[-(num_angles + 1)]
right_angles_nomial = seq(0,180,length.out = num_angles + 1)[-(num_angles + 1)]
####### Left OB: Partition a region into two parts using a rotating axis
results <- list()
for (angle in left_angles) {
  # Convert angle to radians for computation
  angle_rad <- angle * pi / 180

  # Handle vertical line case (90 and 270 degrees)
  if (cos(angle_rad) == 0) {
    coord3d_IRIS_slice_domain_left$cluster <- ifelse(coord3d_IRIS_slice_domain_left$x.center.rotate > 0, "Cluster 1", "Cluster 2")
  } else {
    # Compute the slope of the line based on the angle
    slope <- sin(angle_rad) / cos(angle_rad)

    # Compute the y-value of the line for each x-value
    coord3d_IRIS_slice_domain_left$line_y <- slope * coord3d_IRIS_slice_domain_left$x.center.rotate

    # Label points based on their position relative to the line
    coord3d_IRIS_slice_domain_left$cluster <- ifelse(coord3d_IRIS_slice_domain_left$y.center.rotate > coord3d_IRIS_slice_domain_left$line_y, "Cluster 1", "Cluster 2")
  }

  # Store the results for this angle
  results[[paste0("Angle_", angle)]] <- data.frame(x = coord3d_IRIS_slice_domain_left$x.center.rotate,
                                                   y = coord3d_IRIS_slice_domain_left$y.center.rotate,
                                                   cluster = coord3d_IRIS_slice_domain_left$cluster,
                                                   slice_barcode = coord3d_IRIS_slice_domain_left$slice_barcode,
                                                   angle = angle)
}
final_left <- do.call(rbind, results)


glomer_axis_rad_left  <- glomer_axis_angle_left  * pi / 180
slope1 <- sin(glomer_axis_rad_left) / cos(glomer_axis_rad_left)
p = ggplot(final_left, aes(x = x, y = y, color = cluster)) +
  geom_point(size = 0.5) +
  geom_abline(color = "blue",slope = slope1,size = 2,linetype = "dashed") +
  labs(title = "2D Clustering by a Rotating Line",
       x = "X Coordinate",
       y = "Y Coordinate") +
  facet_wrap(~ angle,nrow = 4) +
  theme_minimal()

pdf(paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit/Left_Right_Overlaid_RotatingAxes_Map/",
           domain_algorithm,"_",
           select_domain_name,".pdf"),width = 8,height = 10)
p
dev.off()
# 
# # ####### Right OB: Partition a region into two parts using a rotating axis
# # results <- list()
# # for (angle in right_angles) {
# #   # Convert angle to radians for computation
# #   angle_rad <- angle * pi / 180
# # 
# #   # Handle vertical line case (90 and 270 degrees)
# #   if (cos(angle_rad) == 0) {
# #     coord3d_IRIS_slice_domain_right$cluster <- ifelse(coord3d_IRIS_slice_domain_right$x.center.rotate > 0, "Cluster 1", "Cluster 2")
# #   } else {
# #     # Compute the slope of the line based on the angle
# #     slope <- sin(angle_rad) / cos(angle_rad)
# # 
# #     # Compute the y-value of the line for each x-value
# #     coord3d_IRIS_slice_domain_right$line_y <- slope * coord3d_IRIS_slice_domain_right$x.center.rotate
# # 
# #     # Label points based on their position relative to the line
# #     coord3d_IRIS_slice_domain_right$cluster <- ifelse(coord3d_IRIS_slice_domain_right$y.center.rotate > coord3d_IRIS_slice_domain_right$line_y, "Cluster 1", "Cluster 2")
# #   }
# # 
# #   # Store the results for this angle
# #   results[[paste0("Angle_", angle)]] <- data.frame(x = coord3d_IRIS_slice_domain_right$x.center.rotate,
# #                                                    y = coord3d_IRIS_slice_domain_right$y.center.rotate,
# #                                                    cluster = coord3d_IRIS_slice_domain_right$cluster,
# #                                                    slice_barcode = coord3d_IRIS_slice_domain_right$slice_barcode,
# #                                                    angle = angle)
# # }
# # final_right <- do.call(rbind, results)
# 
# # ggplot(final_right, aes(x = x, y = y, color = cluster)) +
# #   geom_point(size =  0.5) +
# #   labs(title = "2D Clustering by a Rotating Line",
# #        x = "X Coordinate",
# #        y = "Y Coordinate") +
# #   facet_wrap(~ angle,nrow = 4) +
# #   theme_minimal()
# 
# ####### Left OB: Differential Expression
# DEA_left = NULL
# for(a in 1:length(left_angles)){
#   meta.left = final_left %>% filter(angle == left_angles[a])
#   count.left = countMat[gene_list,meta.left$slice_barcode]
#   mob_seurat = CreateSeuratObject(count = count.left)
#   mob_seurat@assays$RNA$data =  count.left
#   Idents(mob_seurat) = as.factor(meta.left$cluster)
#   dea = FindAllMarkers(mob_seurat,logfc.threshold = 0)
# 
#   if(nrow(dea) == 0){
#     next
#   }
#   dea = dea %>% filter(!duplicated(gene))
#   dea$angle = left_angles[a]
#   dea$semi = "left"
# 
#   # DEA_Empty = data.frame(gene = gene_list,
#   #                        angle = angles[a],
#   #                        slice = slice_list[s],
#   #                        semi = "left")
#   # dea = merge(dea,DEA_Empty,
#   #             by.x = c("gene","angle","slice","semi"),
#   #             by.y = c("gene","angle","slice","semi"),
#   #             all.y = T)
# 
#   DEA_left = rbind(DEA_left,dea)
# }
# 
# 
# DEA_rotateAxis = DEA_left
# 
# saveRDS(DEA_rotateAxis,paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit/Left_Right_Overlaid_DEdataframe/",
#                               domain_algorithm,"_",select_domain_name,
#                               "_DEgene_DEaxes.rds"))

DEA_rotateAxis = readRDS(paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit/Left_Right_Overlaid_DEdataframe/",
                                                              domain_algorithm,"_",select_domain_name,
                                                               "_DEgene_DEaxes.rds"))

##### Filter DE Cases by different P-value threshold
DEA_rotateAxis_left_pval000001 = DEA_rotateAxis %>% filter(semi == "left" & p_val_adj <= 0.05)
DEA_rotateAxis_right_pval000001 = DEA_rotateAxis %>% filter(semi == "right" & p_val_adj <= 0.05)


### Left 
DEA_rotateAxis_left_pval000001_optForEachGene = DEA_rotateAxis_left_pval000001 %>%
  group_by(gene) %>%
  top_n(1,abs(avg_log2FC)) %>%
  ungroup()

### Right 
DEA_rotateAxis_right_pval000001_optForEachGene = DEA_rotateAxis_right_pval000001 %>%
  group_by(gene) %>%
  top_n(1,abs(avg_log2FC)) %>%
  ungroup()

## Top 100 Genes with highest logfold change
DEA_rotateAxis_left_pval000001_optForEachGene_topCase = DEA_rotateAxis_left_pval000001_optForEachGene %>% top_n(50,abs(avg_log2FC))
DEA_rotateAxis_right_pval000001_optForEachGene_topCase = DEA_rotateAxis_right_pval000001_optForEachGene %>% top_n(50,abs(avg_log2FC))


saveRDS(DEA_rotateAxis_left_pval000001_optForEachGene_topCase,paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit/Left_Right_Overlaid_DEdataframe_Top50//",
                              domain_algorithm,"_",select_domain_name,
                              "_DEgene_DEaxes_Top50Bylog2FC.rds"))

#### Leftp
##### Given Panel Equation from Mihaly  6.7446x + 1.3679y + z  + 8280.1816 = 0
library(ggplot2)
library(ggforce)
library(ggrepel)
library(dplyr)
df <- DEA_rotateAxis_left_pval000001_optForEachGene_topCase %>%
  select(gene,p_val_adj,avg_log2FC,angle)

df$angle_mihaly = - (df$angle - glomer_axis_angle_left )

colnames(df) = c("gene","pval","log2FC","angle","angle_mihaly")
df = df %>%
  mutate(pvalFactor = case_when(pval < 0.0001 ~ "p<0.0001",
                                pval >= 0.0001 & pval < 0.0005 ~ "0.0001<p<0.0005",
                                pval >= 0.0005 ~ "p>0.0005"),
         log2FCabs = abs(log2FC)
         )
df$pvalFactor = factor(df$pvalFactor,levels = c("p>0.0005","0.0001<p<0.0005","p<0.0001"))

# Convert log2FC to polar coordinates
df <- df %>%
  mutate(
    radius = log2FCabs,  # Shift log2FC to be positive for better visualization
    x = case_when(log2FC > 0 ~ radius * cos(angle * pi / 180),
                  log2FC < 0 ~ radius * cos((180+angle) * pi / 180)),
    y = case_when(log2FC > 0 ~ radius * sin(angle * pi / 180),
                  log2FC < 0 ~ radius * sin((180+angle) * pi / 180))
  )


# angles_apdx <- c(left_angles,162 + c(1:10) * 18)# Exclude the last point to avoid duplicating 0 degrees
angles_apdx = c(left_angles,left_angles[length(left_angles)] - c(1:10) * 18)
radius_circle = max(df$radius)



angle_labels <- data.frame(
  angles_apdx = angles_apdx ,  # Angle labels from 0° to 330° in steps of 30°
  angle_label = c(left_angles_nomial,left_angles_nomial),
  radius_circle = radius_circle + 0.2  # Place the labels slightly outside the circle
) %>%
  mutate(
    x = radius_circle * cos(angles_apdx * pi / 180),
    y = radius_circle * sin(angles_apdx * pi / 180),
    label = paste0(angle_label, "°")  # Format labels
  )

# df$gene_show = df$gene
# df$gene_show[which(!df$gene_show %in% c("Clca3a1", "Klk6", "Gal3st1", "Cldn11"))] = NA
df$gene_show = NA

glomer_axis_rad_left  <- glomer_axis_angle_left  * pi / 180
slope <- sin(glomer_axis_rad_left) / cos(glomer_axis_rad_left)
# Create the plot
p = ggplot(df) +
  # Circle outline
  geom_circle(aes(x0 = 0, y0 = 0, r = radius_circle), color = "black",size = 2) +
  geom_text(data = angle_labels, aes(x = x, y = y, label = label), size =14) +
  # Dashed lines for each angle
  geom_segment(data = data.frame(angles_apdx = angles_apdx),aes(x = 0, y = 0, xend = cos(angles_apdx * pi / 180) * radius_circle ,
                   yend = sin(angles_apdx * pi / 180) * radius_circle), linetype = "dashed", color = "black") +

  # Points for genes
  geom_point(aes(x = x, y = y, size = pvalFactor)) +
  # Labels for genes
  # geom_text(aes(x = x, y = y, label = gene), hjust = 0, vjust = 0, size = 10) +
  geom_abline(slope = slope,size = 2, color = "red") +
  # geom_text_repel(aes(x = x, y = y, label = gene_show),
  #                  size = 12,
  #                  direction = "both",
  #                  max.overlaps = Inf,
  #                  force = 2,
  #                  nudge_x = 2,
  #                  nudge_y = ,
  #                 box.padding = 0.5) +

  # Theme adjustments
  labs(size = "Adjusted P-value") +
  scale_x_continuous(limits=c(-2.5,2.5)) + 
  scale_y_continuous(limits=c(-2.5,2.5)) + 
  # coord_fixed() +
  # theme_void() +
  scale_size_manual(values = c(3,6,12)) +  # Adjust point sizes
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30)) +
  guides(
    size = guide_legend(nrow = 3)
  )

pdf(paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit//Left_Right_Overlaid_radiusmap//",
           domain_algorithm,"_",
           select_domain_name,"_",
           "Top50_DEgeneBylog2FC.pdf"
           ),
    width = 11,height =12.5)
p
dev.off()


## Top 50 Genes with highest logfold change Gene Expression map
countMat = readRDS("~/data/Collaboration/AlexF/MOB/AllSection_NormCount.rds")


#### Left OB: Differential Expression Signature ######
### Top 1 Axis
de_axis = DEA_rotateAxis_left_pval000001_optForEachGene_topCase
ang = -47
ang_de_axis = de_axis %>% filter(angle == ang)
ang_de_axis_pos = ang_de_axis %>% filter(avg_log2FC > 0)
ang_de_axis_neg = ang_de_axis %>% filter(avg_log2FC < 0)
##### Gene with Positive Log2FC
gene = ang_de_axis_pos$gene
meta.left = final_left %>% filter(angle == ang)
count = countMat[gene,meta.left$slice_barcode, drop=FALSE]
spatial = coord3d_IRIS[coord3d_IRIS$slice_barcode %in% meta.left$slice_barcode,]
rownames(spatial) = spatial$slice_barcode
spatial = spatial[colnames(count),]
all(names(count) == rownames(spatial))
spatial$center.x = spatial$x.x - mean(spatial$x.x)
spatial$center.y = spatial$y.x - mean(spatial$y.x)
spatial$expr = as.numeric(colMeans(count))

glomer_axis_rad_left  <- glomer_axis_angle_left  * pi / 180
slope1 <- sin(glomer_axis_rad_left) / cos(glomer_axis_rad_left)
ang_rad  <- (ang)  * pi / 180
slope2 <- sin(ang_rad) / cos(ang_rad)

##### All Different gene for Angel
gene = ang_de_axis$gene
meta.left = final_left %>% filter(angle == ang)
count = countMat[gene,meta.left$slice_barcode,drop=FALSE]
spatial = coord3d_IRIS_slice_domain_left[coord3d_IRIS_slice_domain_left$slice_barcode %in% meta.left$slice_barcode,]
rownames(spatial) = spatial$slice_barcode
spatial = spatial[colnames(count),]
all(names(count) == rownames(spatial))

spatial$center.x = spatial$x.center.flipped
spatial$center.y = spatial$y.center.flipped

spatial$expr = as.numeric(colMeans(count))

glomer_axis_rad_left  <- glomer_axis_angle_left  * pi / 180
slope1 <- sin(glomer_axis_rad_left) / cos(glomer_axis_rad_left)
ang_rad  <- (ang)  * pi / 180
slope2 <- sin(ang_rad) / cos(ang_rad)
p.all =
  spatial %>%
  # filter(section %in% 107:127) %>%
  ggplot(aes(x = center.x, y = center.y)) +
  geom_point(aes(color = expr),size = 3) +
  labs(color = "Signature") + 
  # scale_colour_gradientn(colours = c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d80032")) +
  # scale_colour_gradientn(colours = c("#ebd4cb","#da9f93","#b6465f","#890620","#2c0703")) +
  geom_abline(color = "blue",slope = slope1,size = 2,linetype = "dashed") +
  geom_abline(color = "green",slope = slope2,size = 2,linetype = "dashed") +
  # viridis::scale_colour_viridis(option = "D") + 
  # scale_colour_gradientn(colors = rev(c("#970005","#ed0101","#ffffff","#0c44ac","#000052"))) +
  labs(title = paste0("Angle=",left_angles_nomial[which(left_angles == ang)],"")) +
  # facet_wrap(~ slice,nrow = 3,scales = "free") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 30),
        strip.text = element_text(size = 30),
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 20))


pdf(paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit//Left_Right_Overlaid_exprSignature///",
           domain_algorithm,"_",select_domain_name,
           "_left_angle",
           left_angles_nomial[which(left_angles == ang)],
           "_all.pdf"),width = 12,height = 12)
print(cowplot::plot_grid(p.all,nrow = 1))
dev.off()

library(RColorBrewer)
colormap = brewer.pal(9, "YlOrRd")
p.all =
  spatial %>%
  # filter(section %in% 107:127) %>%
  ggplot(aes(x = center.x, y = center.y)) +
  geom_point(aes(color = expr),size = 3) +
  labs(color = "Signature") + 
  scale_colour_gradientn(colours = colormap) +
  # scale_colour_gradientn(colours = c("#ebd4cb","#da9f93","#b6465f","#890620","#2c0703")) +
  geom_abline(color = "blue",slope = slope1,size = 2,linetype = "dashed") +
  geom_abline(color = "green",slope = slope2,size = 2,linetype = "dashed") +
  # viridis::scale_colour_viridis(option = "D") + 
  # scale_colour_gradientn(colors = rev(c("#970005","#ed0101","#ffffff","#0c44ac","#000052"))) +
  labs(title = paste0("Angle=",left_angles_nomial[which(left_angles == ang)],"")) +
  # facet_wrap(~ slice,nrow = 3,scales = "free") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 30),
        strip.text = element_text(size = 30),
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 20))




pdf(paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit//Left_Right_Overlaid_exprSignature///",
           domain_algorithm,"_",select_domain_name,
           "_left_angle",
           left_angles_nomial[which(left_angles == ang)],
           "_all_YlOrRd.pdf"),width = 12,height = 12)
print(cowplot::plot_grid(p.all,nrow = 1))
dev.off()

### Top 2 Axis
de_axis = DEA_rotateAxis_left_pval000001_optForEachGene_topCase
ang = -65
ang_de_axis = de_axis %>% filter(angle == ang)
ang_de_axis_pos = ang_de_axis %>% filter(avg_log2FC > 0)
ang_de_axis_neg = ang_de_axis %>% filter(avg_log2FC < 0)
##### Gene with Positive Log2FC
gene = ang_de_axis_pos$gene
meta.left = final_left %>% filter(angle == ang)
count = countMat[gene,meta.left$slice_barcode, drop=FALSE]
spatial = coord3d_IRIS[coord3d_IRIS$slice_barcode %in% meta.left$slice_barcode,]
rownames(spatial) = spatial$slice_barcode
spatial = spatial[colnames(count),]
all(names(count) == rownames(spatial))
spatial$center.x = spatial$x.x - mean(spatial$x.x)
spatial$center.y = spatial$y.x - mean(spatial$y.x)
spatial$expr = as.numeric(colMeans(count))

glomer_axis_rad_left  <- glomer_axis_angle_left  * pi / 180
slope1 <- sin(glomer_axis_rad_left) / cos(glomer_axis_rad_left)
ang_rad  <- (ang)  * pi / 180
slope2 <- sin(ang_rad) / cos(ang_rad)

##### All Different gene for Angel
gene = ang_de_axis$gene
meta.left = final_left %>% filter(angle == ang)
count = countMat[gene,meta.left$slice_barcode,drop=FALSE]
spatial = coord3d_IRIS_slice_domain_left[coord3d_IRIS_slice_domain_left$slice_barcode %in% meta.left$slice_barcode,]
rownames(spatial) = spatial$slice_barcode
spatial = spatial[colnames(count),]
all(names(count) == rownames(spatial))

spatial$center.x = spatial$x.center.flipped
spatial$center.y = spatial$y.center.flipped

spatial$expr = as.numeric(colMeans(count))

glomer_axis_rad_left  <- glomer_axis_angle_left  * pi / 180
slope1 <- sin(glomer_axis_rad_left) / cos(glomer_axis_rad_left)
ang_rad  <- (ang)  * pi / 180
slope2 <- sin(ang_rad) / cos(ang_rad)
p.all =
  spatial %>%
  # filter(section %in% 107:127) %>%
  ggplot(aes(x = center.x, y = center.y)) +
  geom_point(aes(color = expr),size = 3) +
  labs(color = "Signature") + 
  # scale_colour_gradientn(colours = c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d80032")) +
  # scale_colour_gradientn(colours = c("#ebd4cb","#da9f93","#b6465f","#890620","#2c0703")) +
  geom_abline(color = "blue",slope = slope1,size = 1.5,linetype = "dashed") +
  geom_abline(color = "green",slope = slope2,size = 2,linetype = "dashed") +
  viridis::scale_colour_viridis(option = "D") + 
  # scale_colour_gradientn(colors = rev(c("#970005","#ed0101","#ffffff","#0c44ac","#000052"))) +
  labs(title = paste0("Angle=",left_angles_nomial[which(left_angles == ang)],"")) +
  # facet_wrap(~ slice,nrow = 3,scales = "free") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 30),
        strip.text = element_text(size = 30),
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 20))




pdf(paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit//Left_Right_Overlaid_exprSignature///",
           domain_algorithm,"_",select_domain_name,
           "_left_angle",
           left_angles_nomial[which(left_angles == ang)],
           "_all.pdf"),width = 12,height = 12)
print(cowplot::plot_grid(p.all,nrow = 1))
dev.off()



library(RColorBrewer)
colormap = brewer.pal(9, "YlOrRd")
p.all =
  spatial %>%
  # filter(section %in% 107:127) %>%
  ggplot(aes(x = center.x, y = center.y)) +
  geom_point(aes(color = expr),size = 3) +
  labs(color = "Signature") + 
  scale_colour_gradientn(colours = colormap) +
  # scale_colour_gradientn(colours = c("#ebd4cb","#da9f93","#b6465f","#890620","#2c0703")) +
  geom_abline(color = "blue",slope = slope1,size = 1.5,linetype = "dashed") +
  geom_abline(color = "green",slope = slope2,size = 2,linetype = "dashed") +
  # viridis::scale_colour_viridis(option = "D") + 
  # scale_colour_gradientn(colors = rev(c("#970005","#ed0101","#ffffff","#0c44ac","#000052"))) +
  labs(title = paste0("Angle=",left_angles_nomial[which(left_angles == ang)],"")) +
  # facet_wrap(~ slice,nrow = 3,scales = "free") +
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 30),
        strip.text = element_text(size = 30),
        legend.position = "right",
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 20))




pdf(paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit//Left_Right_Overlaid_exprSignature///",
           domain_algorithm,"_",select_domain_name,
           "_left_angle",
           left_angles_nomial[which(left_angles == ang)],
           "_all_YlOrRd.pdf"),width = 12,height = 12)
print(cowplot::plot_grid(p.all,nrow = 1))
dev.off()



####### DE Gene Expression ########
#### By Log Fc
library(RColorBrewer)
colormap = brewer.pal(9, "YlOrRd")
DE_Top15 = DEA_rotateAxis_left_pval000001_optForEachGene_topCase %>%
  top_n(10,abs(avg_log2FC))
for(g in 1:nrow(DE_Top15)){
  gene = DE_Top15$gene[g]
  ang = DE_Top15$angle[g]
  meta.left = final_left %>% filter(angle == ang)
  count = countMat[gene,meta.left$slice_barcode, drop=FALSE]
  spatial = coord3d_IRIS_slice_domain_left[coord3d_IRIS_slice_domain_left$slice_barcode %in% meta.left$slice_barcode,]
  rownames(spatial) = spatial$slice_barcode
  spatial = spatial[colnames(count),]
  all(names(count) == rownames(spatial))
  spatial$center.x = spatial$x.center.flipped
  spatial$center.y = spatial$y.center.flipped
  spatial$expr = as.numeric(colMeans(count))

  glomer_axis_rad_left  <- glomer_axis_angle_left  * pi / 180
  slope1 <- sin(glomer_axis_rad_left) / cos(glomer_axis_rad_left)
  ang_rad  <- (ang)  * pi / 180
  slope2 <- sin(ang_rad) / cos(ang_rad)
  p.pos =
    spatial %>%
    # filter(section %in% 107:127) %>%
    ggplot(aes(x = center.x, y = center.y)) +
    geom_point(aes(color = expr),size = 3) +
    # scale_colour_gradientn(colours = c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d80032")) +
    # scale_colour_gradientn(colours = c("#ebd4cb","#da9f93","#b6465f","#890620","#2c0703")) +
    geom_abline(color = "blue",slope = slope1,size = 1.5,linetype = "dashed") +
    geom_abline(color = "green",slope = slope2,size = 1.5,linetype = "dashed") +
    # viridis::scale_colour_viridis(option = "D") + 
    scale_colour_gradientn(colours = colormap) +
    labs(title = paste0(gene,",Angle=",left_angles_nomial[which(left_angles == ang)],"")) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 30),
          strip.text = element_text(size = 30),
          legend.position = "right",
          legend.title = element_text(size = 30),
          legend.text = element_text(size = 20))
  
  pdf(paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit//Left_Right_Overlaid_DEgeneExpr_byLFC/",
             domain_algorithm,"_",select_domain_name,
             "_",
             gene,"_Angle_",left_angles_nomial[which(left_angles == ang)],
             "_YlOrRd.pdf"),width = 12,height = 12)
  print(cowplot::plot_grid(p.pos,nrow = 1))
  dev.off()
}


#### By pVALUE
DE_Top15 = DEA_rotateAxis_left_pval000001_optForEachGene_topCase %>%
  top_n(10,-p_val_adj)
for(g in 1:nrow(DE_Top15)){
  gene = DE_Top15$gene[g]
  ang = DE_Top15$angle[g]
  meta.left = final_left %>% filter(angle == ang)
  count = countMat[gene,meta.left$slice_barcode, drop=FALSE]
  spatial = coord3d_IRIS_slice_domain_left[coord3d_IRIS_slice_domain_left$slice_barcode %in% meta.left$slice_barcode,]
  rownames(spatial) = spatial$slice_barcode
  spatial = spatial[colnames(count),]
  all(names(count) == rownames(spatial))
  spatial$center.x = spatial$x.center.flipped
  spatial$center.y = spatial$y.center.flipped
  spatial$expr = as.numeric(colMeans(count))
  
  glomer_axis_rad_left  <- glomer_axis_angle_left  * pi / 180
  slope1 <- sin(glomer_axis_rad_left) / cos(glomer_axis_rad_left)
  ang_rad  <- (ang)  * pi / 180
  slope2 <- sin(ang_rad) / cos(ang_rad)
  p.pos =
    spatial %>%
    # filter(section %in% 107:127) %>%
    ggplot(aes(x = center.x, y = center.y)) +
    geom_point(aes(color = expr),size = 3) +
    # scale_colour_gradientn(colours = c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d80032")) +
    # scale_colour_gradientn(colours = c("#ebd4cb","#da9f93","#b6465f","#890620","#2c0703")) +
    geom_abline(color = "blue",slope = slope1,size = 1.5,linetype = "dashed") +
    geom_abline(color = "green",slope = slope2,size = 1.5,linetype = "dashed") +
    # viridis::scale_colour_viridis(option = "D") + 
    scale_colour_gradientn(colours = colormap) +
    labs(title = paste0(gene,",Angle=",left_angles_nomial[which(left_angles == ang)],"")) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size = 30),
          strip.text = element_text(size = 30),
          legend.position = "right",
          legend.title = element_text(size = 30),
          legend.text = element_text(size = 20))
  
  pdf(paste0("~/data/Collaboration/AlexF/MOB/IRIS_identifyAxis_optimizeDE_Submit//Left_Right_Overlaid_DEgeneExpr_byPVal/",
             domain_algorithm,"_",select_domain_name,
             "_",
             gene,"_Angle_",left_angles_nomial[which(left_angles == ang)],
             "_YlOrRd.pdf"),width = 12,height = 12)
  print(cowplot::plot_grid(p.pos,nrow = 1))
  dev.off()
}


