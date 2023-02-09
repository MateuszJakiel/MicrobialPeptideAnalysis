#Load dependencies
library(ggfortify)
library(ggplot2)
library(umap)
library(tidymodels)
library(tidyr)
library(dplyr)
library(MASS)
library(readJDX)
library(caret)
library(GGally)
library(rpart)
library(rpart.plot)
library(dbscan)										                       
library(rmcfs)

#Load and clean data
qpID <- readRDS("qpID.rds")
data <- read.csv('complete_database.csv')
data <- data[-1]
data <- data[-5]
t_data <- data[69:170]
data <- data[1:68,]

#Perform PCA reduction
pca_data <- data[-14]
pca_data <- pca_data[vapply(pca_data, function(x) length(unique(x)) > 1, logical(1L))]
pca <- prcomp(pca_data, scale = TRUE)
pca_data$Receptor.nr <- as.factor(pca_data$Receptor.nr)
autoplot(pca, data = pca_data, colour = 'Receptor.nr', label = TRUE)

receptor <- data$Receptor.nr
receptor <- as.data.frame(receptor[-29])
receptor <- as.data.frame(receptor[-33,])

be <- data$beg_end
be <- as.data.frame(be[-29])
be <- as.data.frame(be[-33,])

pca_data <- pca_data[-34,]
pca_data <- pca_data[-29,]

pca_data$Receptor.nr <- as.numeric(pca_data$Receptor.nr)
pca_data <- pca_data[vapply(pca_data, function(x) length(unique(x)) > 1, logical(1L))]
pca <- prcomp(pca_data, scale = TRUE)
pca_data$Receptor.nr <- as.factor(pca_data$Receptor.nr)
autoplot(pca, data = pca_data, label = TRUE) 

#Perform UMAP reduction
umap_data <- pca_data
custom_config <- umap.defaults
custom_config$n_neighbors <- 5

manifolds_n5 <- umap(pca_data, config = custom_config)
umap_plot_n5 <- as.data.frame(cbind(manifolds_n5$layout, umap_data$Receptor.nr))
umap_plot_n5$V3 <- as.factor(umap_plot_n5$V3)

n5 <- ggplot(
  umap_plot_n5,
  aes(
    x = V1,
    y = V2)) +
  geom_point() + labs(title = "UMAP with kNN set to n = 5", x = "UMAP-X", y = "UMAP-Y")

custom_config <- umap.defaults
custom_config$n_neighbors <- 15

manifolds_n15 <- umap(umap_data, config = custom_config)
umap_plot_n15 <- as.data.frame(cbind(manifolds_n15$layout, umap_data$Receptor.nr))
umap_plot_n15$V3 <- as.factor(umap_plot_n15$V3)

n15 <- ggplot(
  umap_plot_n15,
  aes(
    x = V1,
    y = V2)) +
  geom_point() + labs(title = "UMAP with kNN set to n = 15", x = "UMAP-X", y = "UMAP-Y")

custom_config <- umap.defaults
custom_config$n_neighbors <- 20

manifolds_n20 <- umap(umap_data, config = custom_config)
umap_plot_n20 <- as.data.frame(cbind(manifolds_n20$layout, umap_data$Receptor.nr))
umap_plot_n20$V3 <- as.factor(umap_plot_n20$V3)

n20 <- ggplot(
  umap_plot_n20,
  aes(
    x = V1,
    y = V2)) +
  geom_point() + labs(title = "UMAP with kNN set to n = 20", x = "UMAP-X", y = "UMAP-Y")

list_umap <- list(n5, n15, n20)

cowplot::plot_grid(plotlist = list_umap, ncol = 3, nrow = 1)

#Perform k-means clustering
kclusts_pca <- 
  tibble(k = 1:9) %>%
  mutate(
    kclust = map(k, ~kmeans(pca_data, .x)),
    tidied = map(kclust, tidy),
    glanced = map(kclust, glance),
    augmented = map(kclust, augment, pca_data)
  )

clusterings_pca <- 
  kclusts_pca %>%
  unnest(cols = c(glanced))

ggplot(clusterings_pca, aes(k, tot.withinss)) +
  geom_line() +
  geom_point() + labs(x = "Number of clusters", y = "Total within-cluster sum of squares")
k_means_pca <- kmeans(pca_data, centers = 2, iter.max = 20)

rot <- as.data.frame(pca$x[,1:2])

table(k_means_pca$cluster)

pca_cluster <- ggplot(rot, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = as.factor(k_means_pca$cluster))) + labs(x = "PC1", y = "PC2", color = "Clusters")


umap_cluster <- ggplot(
  umap_plot,
  aes(
    x = V1,
    y = V2)) +
  geom_point(aes(col = as.factor(k_means_pca$cluster))) + labs(color = "Clusters")

k_means_pca2 <- kmeans(pca_data, centers = 5, iter.max = 20)

table(k_means_pca2$cluster)

pca_cluster_k5 <- ggplot(rot, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = as.factor(k_means_pca2$cluster))) + labs(x = "PC1", y = "PC2", color = "Clusters")


umap_cluster_k5 <- ggplot(
  umap_plot,
  aes(
    x = V1,
    y = V2)) +
  geom_point(aes(col = as.factor(k_means_pca2$cluster))) + labs(x = "UMAP-X", y = "UMAP-Y", color = "Clusters")

rot <- as.data.frame(pca[["x"]])
rownames(rot) <- qpID

k_means_pca2 = kmeans(rot[1:2], centers = 2, nstart = 50)
k2_pca <- factoextra::fviz_cluster(k_means_pca2, data = rot[1:2], geom = "text", show.clust.cent = FALSE)

k_means_umap_2 = kmeans(umap_plot_n5[1:2], centers = 2, nstart = 50)
k2_umap <- factoextra::fviz_cluster(k_means_umap_2, data = umap_plot_n5[1:2], geom = "text", show.clust.cent = FALSE)

k_means_pca_5 = kmeans(rot[1:2], centers = 5, nstart = 50)
k5_pca <- factoextra::fviz_cluster(k_means_pca_5, data = rot[1:2], geom = "text", show.clust.cent = FALSE)

k_means_umap_5 = kmeans(umap_plot_n5[1:2], centers = 5, nstart = 50)
k5_umap <- factoextra::fviz_cluster(k_means_umap_5, data = umap_plot_n5[1:2], geom = "text", show.clust.cent = FALSE)

k2_pca <- k2_pca + labs(title = "PCA with k-means clustering, k = 2", x = "PC1 (7.59%)", 
                        y = "PC2 (5.56%)", color = "Clusters") + guides(fill = "none")
k5_pca <- k5_pca + labs(title = "PCA with k-means clustering, k = 5", x = "PC1 (7.59%)", 
                        y = "PC2 (5.56%)", color = "Clusters") + guides(fill = "none")
k2_umap <- k2_umap + labs(title = "UMAP with k-means clustering, k = 2", x = "UMAP-X", 
                          y = "UMAP-Y", color = "Clusters") + guides(fill = "none")
k5_umap <- k5_umap + labs(title = "UMAP with k-means clustering, k = 5", x = "UMAP-X", 
                          y = "UMAP-Y", color = "Clusters") + guides(fill = "none")


k_means_list <- list(k2_pca, k2_umap, k5_pca, k5_umap)

cowplot::plot_grid(plotlist = k_means_list, ncol = 2, nrow = 2)


k_means_comparison <- data.frame(1:66)
k_means_comparison$PCA_k5 <- k_means_pca_5$cluster
k_means_comparison$UMAP_k5 <- k_means_umap_5$cluster
k_means_comparison <- k_means_comparison %>% group_by(PCA_k5) %>% count(UMAP_k5)
k_means_comparison$UMAP_k5 <- as.factor(k_means_comparison$UMAP_k5)

k_means_comparison_plot <- ggplot(data=k_means_comparison, aes(x=PCA_k5, y=n, fill = UMAP_k5)) +
  geom_bar(stat="identity", position=position_dodge()) 
k_means_comparison_plot +  labs(x = "PCA k-means clusters, k = 5", y = "Number of common samples", fill = "UMAP clusters")

#Density-based spatial clustering of applications with noise - DBSCAN
dbscan_clustering5 <- hdbscan(umap_plot_n5[1:2], minPts = 3)
dbscan_clustering15 <- hdbscan(umap_plot_n15[1:2], minPts = 3)
dbscan_clustering20 <- hdbscan(umap_plot_n20[1:2], minPts = 3)

umap_plot_n5$groups <- as.factor(dbscan_clustering$cluster)
umap_dbscan <- ggplot(
  umap_plot_n5,
  aes(
    x = V1,
    y = V2, color = groups)) +
  geom_point() + 
  labs(x = "UMAP-X", y = "UMAP-Y", color = "Clusters", title = "UMAP with DBSCAN") +
  stat_ellipse(geom = "polygon",
               aes(fill = groups), 
               alpha = 0.25) +
  guides(fill = "none")

pca_plot_input$groups <- as.factor(dbscan_clustering$cluster)
pca_dbscan <- ggplot(pca_plot_input, aes(x = PC1, y = PC2, color = groups)) +
  geom_point() + labs(x = "PC1 (7.59%)", y = "PC2 (5.56%)", color = "Cluster", title = "PCA with DBSCAN") + 
  stat_ellipse(geom = "polygon",
               aes(fill = groups), 
               alpha = 0.25) +
  guides(fill = "none")

list3 <- list(pca_dbscan, umap_dbscan)
cowplot::plot_grid(plotlist = list3, nrow = 1, ncol = 2)


#Perform LDA
grouping_input <- readRDS("grouping_input.rds")
grouping_input$qpID <- qpID
smp_size <- floor(0.95 * nrow(grouping_input))
set.seed(121)
train_ind <- sample(seq_len(nrow(grouping_input)), size = smp_size)

train <- grouping_input[train_ind, ]
test <- grouping_input[-train_ind, ]
test_groups <- test$groups

#Perform regression tree
tree_test_super <- rpart(groups~., train, method ="class")     predictions_tree <- predict(tree_test_super, test, type = "class")
paste("The prediction accuracy of the grouping is", mean(predictions_tree==test$groups))

test <- subset(test, select=-c(groups))


LDA_results <- lda(groups~., data = train)

predictions_grouping <- LDA_results %>% predict(test)
paste("The prediction accuracy of the grouping is", mean(predictions_super$class==test_groups))

values <- predict(LDA_results)

#Perform Monte Carlo Feature Selection
MCFS <- mcfs(formula = groups~., train, 
             threadsNumber = 12, mode = 1, 
             finalRuleset = FALSE, splitSetSize    = 1000)
group_attributes <- MCFS$RI$attribute[1:MCFS$cutoff_value]


#Plot PCA and UMAP
pca_clean <- ggplot(rot, aes(x = PC1, y = PC2)) +
  geom_text(aes(label = as.character(qpID))) + labs(x = "PC1 (7.59%)", y = "PC2 (5.56%)", text = "QuorumPeps ID", title = "PCA")

umap_clean <- ggplot(
  umap_plot_n5,
  aes(
    x = V1,
    y = V2)) +
  geom_text(aes(label = as.character(qpID))) + labs(x = "UMAP-X", y = "UMAP-Y", text = "QuorumPeps ID", title = "UMAP")

list <- list(pca_clean, umap_clean)
cowplot::plot_grid(plotlist = list, nrow = 1, ncol = 2)

#Perform the manual grouping according to self created principles
best <- data.frame(pca_data$Viability==1 & pca_data$Differentiation==1 &
                     pca_data$Inflammation==0 & pca_data$GDF15==0)
worst <- data.frame(pca_data$Viability==0 & pca_data$Differentiation==0 &
                      pca_data$Inflammation==1 & pca_data$GDF15==1)

V1 <- data.frame(pca_data$Viability==1 & pca_data$Differentiation==0 &
                   pca_data$Inflammation==0 & pca_data$GDF15==0)

D1 <- data.frame(pca_data$Viability==0 & pca_data$Differentiation==1 &
                   pca_data$Inflammation==0 & pca_data$GDF15==0)

meh1 <- data.frame(pca_data$Viability==1 & pca_data$Differentiation==0)

meh2 <- data.frame(pca_data$Viability==0 & pca_data$Differentiation==1)

bad1 <- data.frame(pca_data$Inflammation==1 & pca_data$GDF15==0)

bad2 <- data.frame(pca_data$Inflammation==0 & pca_data$GDF15==1)

nothing1 <-  data.frame(pca_data$Viability==1 & pca_data$Differentiation==1 &
                          pca_data$Inflammation==1 & pca_data$GDF15==1)

nothing2 <-  data.frame(pca_data$Viability==0 & pca_data$Differentiation==0 &
                          pca_data$Inflammation==0 & pca_data$GDF15==0)


grouping <- bind_cols(best, worst, V1, D1, meh1, meh2, bad1, bad2, nothing1, nothing2)
names(grouping)[1] <- "best"
names(grouping)[2] <- "worst"
names(grouping)[3] <- "V1"
names(grouping)[4] <- "D1"
names(grouping)[5] <- "meh1"
names(grouping)[6] <- "meh2"
names(grouping)[7] <- "bad1"
names(grouping)[8] <- "bad2"
names(grouping)[9] <- "nothing1"
names(grouping)[10] <- "nothing2"


grouping$best <- as.numeric(grouping$best)
grouping$worst <- as.numeric(grouping$worst)
grouping$V1 <- as.numeric(grouping$V1)
grouping$D1 <- as.numeric(grouping$D1)
grouping$med1 <- as.numeric(grouping$meh1)
grouping$med2 <- as.numeric(grouping$meh2)
grouping$bad1 <- as.numeric(grouping$bad1)
grouping$bad2 <- as.numeric(grouping$bad2)
grouping$nothing1 <- as.numeric(grouping$nothing1)
grouping$nothing2 <- as.numeric(grouping$nothing2)


grouping$positive <- grouping$V1 + grouping$D1
grouping$negative <- grouping$bad1 + grouping$bad2
grouping$medium <- grouping$med1 + grouping$med2
grouping$null_case <- grouping$nothing1 + grouping$nothing2

grouping <- grouping[-(3:12)]

grouping <- as.data.frame(grouping)

grouping$medium[21] <- 0
grouping$medium[26] <- 0
grouping$medium[29] <- 0
grouping$medium[34] <- 0
grouping$medium[36] <- 0
grouping$medium[37] <- 0
grouping$medium[39] <- 0
grouping$medium[51] <- 0
grouping$medium[2] <- 0
grouping$medium[3] <- 0
grouping$medium[15] <- 0
grouping$medium[22] <- 0
grouping$medium[27] <- 0
grouping$medium[31] <- 0
grouping$medium[35] <- 0
grouping$medium[38] <- 0
grouping$medium[42] <- 0
grouping$medium[44] <- 0
grouping$medium[49] <- 0
grouping$medium[53] <- 0
grouping$medium[54] <- 0
grouping$medium[56] <- 0


names <- colnames(grouping)

for (x in 1:ncol(grouping)){
  for (i in 1:nrow(grouping)) {
    if (grouping[i,x] == 1){
      grouping[i,x] <- names[x]
      next
    } else {
      grouping[i,x] <- NA
    }
  } 
}

groups <- grouping %>% unite("groups", na.rm = TRUE, remove = FALSE)
groups <- groups[1]

pca_data$groups <- groups$groups
grouping_input <- pca_data[4:1561]
rm(list=ls()) 

#Visualization of manual grouping on PCA and UMAP
grouping_input$groups <- factor(grouping_input$groups, levels = c("best", "positive", "medium", "negative", "worst" ,"null_case"))

umap_grouping <- ggplot(
  umap_plot_n5,
  aes(
    x = V1,
    y = V2)) +
  geom_point(aes(col = grouping_input$groups)) + 
  labs(x = "UMAP-X", y = "UMAP-Y", color = "Groups", title = "UMAP with manual grouping")

pca_grouping <- ggplot(rot, aes(x = PC1, y = PC2)) +
  geom_point(aes(col = grouping_input$groups)) + 
  labs(x = "PC1 (7.59%)", y = "PC2 (5.56%)", color = "Groups", title = "PCA with manual grouping") + theme(legend.position="none")

list3 <- list(pca_grouping, umap_grouping)
cowplot::plot_grid(plotlist = list3, nrow = 1, ncol = 2)

#Comparison of k-means clustering and DBSCAN
dbscan_kmeans_comp <- data.frame(1:66)
dbscan_kmeans_comp$dbscan <- as.factor(dbscan_clustering$cluster)
dbscan_kmeans_comp$kmeans <- k_means_umap_5$cluster
dbscan_kmeans_comp <- dbscan_kmeans_comp %>% group_by(dbscan) %>% count(kmeans)
dbscan_kmeans_comp$kmeans <- as.factor(dbscan_kmeans_comp$kmeans)

ggplot(data=dbscan_kmeans_comp, aes(x=dbscan, y=n, fill = kmeans)) +
  geom_bar(stat="identity", position=position_dodge()) + labs(x = "DBSCAN clusters", y = " Number of common samples", fill = "K-means clusters")

#Comparison of DBSCAN and manual grouping
clusterings <- data.frame(1:66)
clusterings$grouping <- grouping_input$groups
clusterings$grouping <- factor(clusterings$grouping, levels = c("best", "positive", "medium", "negative", "worst" ,"null_case"))
clusterings$dbscan <- dbscan_clustering5$cluster
clusterings <- clusterings %>% group_by(grouping) %>% count(dbscan)
clusterings$dbscan <- as.factor(clusterings$dbscan)

ggplot(data=clusterings, aes(x=grouping, y=n, fill = dbscan)) +
  geom_bar(stat="identity", position=position_dodge()) + labs(x = "Manual grouping classes", y = " Number of samples in DBSCAN")

#Plotting LDA
grouping_input$groups <- factor(grouping_input$groups, 
                                levels = c("best", "positive", "medium", "negative", "worst" ,"null_case"))

lda_plot <- ggpairs(values$x, legend = 1, aes(colour = train$groups, labels = train$groups), upper = list(continuous = wrap(ggally_cor, size = 0, alpha = 0))) + labs(fill = "Groups")

#Plot the regression tree
rpart.plot(test, box.palette = "Greens", type = 4, clip.right.labs = TRUE, branch = .3, under = TRUE, tweak = 1.2)
