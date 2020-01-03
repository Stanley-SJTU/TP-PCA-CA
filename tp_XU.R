# Loading
library("readxl")
library("ggplot2")

# xls files
raw_df <- read_excel("TP4_covC1234_DS19_20.xlsx")
bar_df <- subset(raw_df, TYPE != "?") # remove "?" from TYPE

for (i in colnames(bar_df)){
  
  bar_df[[i]] <- ifelse((bar_df[[i]] == 0), mean(bar_df[[i]], na.rm=TRUE), bar_df[[i]]) # fill "0" with mean
  
}

hiver_mean = colMeans(subset(bar_df, SAISON == "hiver")[, 2:15]) # mean of every composé in "hiver"
ete_mean = colMeans(subset(bar_df, SAISON == "été")[, 2:15])

plot(seq(1, 14, 1), ete_mean)
lines(seq(1, 14, 1), ete_mean, col="blue")
lines(seq(1, 14, 1), hiver_mean, col="green") # much higher gas emission in summer

BF2_mean = colMeans(subset(bar_df, Campagne == "BF2")[, 2:15])
BF3_mean = colMeans(subset(bar_df, Campagne == "BF3")[, 2:15])
CA1_mean = colMeans(subset(bar_df, Campagne == "CA1")[, 2:15])
CA2_mean = colMeans(subset(bar_df, Campagne == "CA2")[, 2:15])
CA3_mean = colMeans(subset(bar_df, Campagne == "CA3")[, 2:15])
CA4_mean = colMeans(subset(bar_df, Campagne == "CA4")[, 2:15])


plot(seq(1, 14, 1), CA2_mean, ylab="Concentration")
lines(seq(1, 14, 1), BF2_mean, col="blue")
lines(seq(1, 14, 1), BF3_mean, col="red")
lines(seq(1, 14, 1), CA1_mean, col="black")
lines(seq(1, 14, 1), CA2_mean, col="green")
lines(seq(1, 14, 1), CA3_mean, col="yellow")
lines(seq(1, 14, 1), CA4_mean, col="pink")
par(xpd=TRUE)
legend("topright", inset=c(-0.1,-0.1), c("BF2","BF3", "CA1", "CA2", "CA3", "CA4"), 
       fill=c("blue","red", "black", "green", "yellow", "pink"))

BF_mean = colMeans(rbind(BF2_mean, BF3_mean))
CA_mean = colMeans(rbind(CA1_mean, CA2_mean, CA3_mean, CA4_mean))

plot(seq(1, 14, 1), CA_mean, ylab="Concentration")
lines(seq(1, 14, 1), BF_mean, col="blue")
lines(seq(1, 14, 1), CA_mean, col="green")

############################## ACP

library(factoextra)

res.pca <- prcomp(bar_df[, 2:15], scale = TRUE)
fviz_eig(res.pca) # visualize the eigenvalue
fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

groups <- as.factor(bar_df$SAISON[1:138])
fviz_pca_ind(res.pca,
             col.ind = groups, # color by groups
             palette = c("red",  "green"),
             legend.title = "Groups",
             repel = TRUE
)

groups <- as.factor(bar_df$Campagne[1:138])
fviz_pca_ind(res.pca,
             geom="point",
             col.ind = groups, # color by groups
             palette = c("red",  "pink", "blue", "brown", "black", "gray"),
             legend.title = "Groups",
             repel = TRUE
)

# Results for Variables
res.var <- get_pca_var(res.pca)
# Results for individuals
res.ind <- get_pca_ind(res.pca)
res.ind.quality = res.ind$cos2           # Quality of representation

########################## AFD

library(MASS)
library(ggplot2)

# AFD sur SAISON

bar_lda = bar_df[, 2:17]
bar_lda[, 15] = NULL 
res.lda = lda(SAISON~., data=bar_lda)

bar_lda_coord = data.matrix(bar_lda[, 1:14]) %*% res.lda$scaling
plot(bar_lda_coord)

saison_df_lda = data.frame(saison = bar_df$SAISON, coord = bar_lda_coord)
plot(saison_df_lda)

var(bar_lda_coord) # variance totale
var(subset(saison_df_lda, saison == "hiver")[, 2]) # variance hiver
var(subset(saison_df_lda, saison == "été")[, 2]) # variance été

# AFD sur TYPE

bar_type_lda = bar_df[, 2:16]
res_type.lda = lda(TYPE~., data=bar_type_lda)

bar_type_lda_coord = data.matrix(bar_type_lda[, 1:14]) %*% res_type.lda$scaling
plot(bar_type_lda_coord[, 1], bar_type_lda_coord[, 2])

type_df_lda = data.frame(type = bar_df$TYPE, coord = bar_type_lda_coord[, 2])
plot(type_df_lda)

# AFD sur BF/CA

bar_ca_lda = bar_df[, 2:18]
bar_ca_lda[, 16] = NULL
bar_ca_lda[, 15] = NULL
bar_ca_lda$Campagne[bar_ca_lda$Campagne == "BF2" | bar_ca_lda$Campagne == "BF3"] <- "BF"  # replace BF2 or BF3 with BF
bar_ca_lda$Campagne[bar_ca_lda$Campagne == "CA1" | bar_ca_lda$Campagne == "CA2" | 
                      bar_ca_lda$Campagne == "CA3" | bar_ca_lda$Campagne == "CA4"] <- "CA"
res_ca.lda = lda(Campagne~., data=bar_ca_lda)

bar_ca_lda_coord = data.matrix(bar_ca_lda[, 1:14]) %*% res_ca.lda$scaling  # new coordinates
plot(bar_ca_lda_coord)

ca_df_lda = data.frame(camp = bar_ca_lda$Campagne, coord = bar_ca_lda_coord)
plot(ca_df_lda)

var(bar_ca_lda_coord) # variance totale
var(subset(ca_df_lda, camp == "BF")[, 2]) # variance BF
var(subset(ca_df_lda, camp == "CA")[, 2]) # variance CA

# Comparison ACP and AFD

bar_pca_coord = res.ind$coord[, 1]
saison_df_pca = data.frame(saison = bar_df$SAISON, coord = bar_pca_coord)

var(bar_pca_coord) # variance totale
var(subset(saison_df_pca, saison == "hiver")[, 2]) # variance hiver
var(subset(saison_df_pca, saison == "été")[, 2]) # variance été

##########################3 AFMD

library(FactoMineR)
bar_famd = bar_df[, 2:19]
res.famd <- FAMD(bar_famd, graph = FALSE)

eig.val <- get_eigenvalue(res.famd)
fviz_screeplot(res.famd)

# variables quantitatives
fviz_famd_var(res.famd, "quanti.var", repel = TRUE,
              col.var = "black")

# variables qualitatives
fviz_famd_var(res.famd, "quali.var", repel = TRUE, col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# variables
# Plot of variables
fviz_famd_var(res.famd, repel = TRUE)
