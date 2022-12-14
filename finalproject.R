# Final Project

##Initial Tests

###setting up environment
setwd("/personal/ksdixo23/BI377/FinalProject")
load("Howell R data.Rdata", verbose = TRUE)

###data cleaning
####raw data (only coords)
raw_data <- Howell_dat$dat.raw
####clean data (only coords)
clean_data <- Howell_dat$dat.super$coords
####list data object for analysis
toepads <- list()
toepads$coords <- Howell_dat$dat.super$coords
toepads$landmark.number <- 123
toepads$specimen.number <- 99
toepads$metadata <- Howell_dat$dat.field
####subset by landmarks
subset_toepads <- subsetgmm(toepads, landmarks = c(1:19)) #go back for semi-landmarks
landmark.plot(subset_toepads)
####create links matrix
links <- matrix(c(1,3, 3,5, 5,18, 18,16, 16,14, 14,12, 12,10, 10,7, 7,9, 9,8, 8,11, 11,13, 13,15, 15,17, 17,19, 19,6, 6,4, 4,2, 10,11, 12,13, 14,15, 16,17, 18,19),
                ncol = 2, byrow = TRUE)
landmark.plot(subset_toepads, links = links)

###gpa analysis
gpa <- align.procrustes(subset_toepads, outlier.analysis = TRUE)

###pca/shape space
pca <- gm.prcomp(gpa$gdf$coords)
shape.space(pca) #PC1 = 28.07%, PC2 = 17.48%
####by habitat
shape.space(pca, group = toepads$metadata$Habitat_Category)
shape.space(pca, group = toepads$metadata$Habitat_Category, convex.hulls = TRUE)
###by municipalities
shape.space(pca, group = toepads$metadata$Municipality)
shape.space(pca, group = toepads$metadata$Municipality, convex.hulls = TRUE)

###regression (Csize + urban/forest or muncipilities)
####null model
i <- 1e4-1
size.model <- procD.lm(coords~log(Csize), data = gpa.fw$gdf, iter = i)
anova(size.model)
####habitat_category model
gpa$gdf$habitat_category <- toepads$metadata$Habitat_Category
habitat.model <- procD.lm(coords~log(Csize) + habitat_category, data = gpa$gdf, iter = i)
anova(habitat.model)
####municipalities model
gpa$gdf$municipality <- toepads$metadata$Municipality
municipality.model <- procD.lm(coords~log(Csize) + municipality, data = gpa$gdf, iter = i)
anova(municipality.model)

###modularity test
modularity.hypothesis1 <- rep("toe", gpa$landmark.number)
modularity.hypothesis1[10:19] <- "lamelle"
m1 <- modularity.test(gpa$gdf$coords, modularity.hypothesis1, CI=TRUE, iter=99)
m1

#############################################################

##Changing Number of Landmarks Test 

##original landmarks plus some sliding landmarks (Kayla)

###data cleaning
####subset by different number of landmarks
subset_toepads2 <- subsetgmm(toepads, landmarks = c(1:19, 22, 25, 30, 33, 41, 47, 56, 62, 71, 73, 78, 81, 86, 89, 94, 97, 102, 105, 110, 113, 118, 121))
subset_toepads2$landmark.number <- 41
landmark.plot(subset_toepads2)
####create links matrix
links2 <- matrix(c(1,20, 20,21, 21,3, 3,22, 22,23, 23,5, 5,18, 18,16, 16,14, 14,12, 12,10, 10,24, 24,7, 7,25, 25,9, 9,26, 26,8, 8,27, 27,11, 11,13, 13,15, 15,17, 17,19, 19,6, 6,28, 28,29, 29,4, 4,30, 30,31, 31,2, 10,32, 32,33, 33,11, 12,34, 34,35, 35,13, 14,36, 36,37, 37,15, 16,38, 38,39, 39,17, 18,40, 40,41, 41,19),
                ncol = 2, byrow = TRUE)
landmark.plot(subset_toepads2, links = links2)

###gpa analysis
gpa2 <- align.procrustes(subset_toepads2, outlier.analysis = TRUE)

###pca/shape space
pca2 <- gm.prcomp(gpa2$gdf$coords)
shape.space(pca2) #PC1 = 23.38%, PC2 = 18.68%
####by habitat
shape.space(pca2, group = toepads$metadata$Habitat_Category)
shape.space(pca2, group = toepads$metadata$Habitat_Category, convex.hulls = TRUE,
            include.legend = TRUE)
###by municipalities
shape.space(pca2, group = toepads$metadata$Municipality)
shape.space(pca2, group = toepads$metadata$Municipality, convex.hulls = TRUE, 
            include.legend = TRUE)

###regression (Csize + urban/forest or muncipilities)
####null model
i <- 1e4-1
size.model <- procD.lm(coords~log(Csize), data = gpa2$gdf, iter = i)
anova(size.model)
####habitat_category model
gpa2$gdf$habitat_category <- toepads$metadata$Habitat_Category
habitat.model <- procD.lm(coords~log(Csize) + habitat_category, data = gpa2$gdf, iter = i)
anova(habitat.model)
####municipalities model
gpa2$gdf$municipality <- toepads$metadata$Municipality
municipality.model <- procD.lm(coords~log(Csize) + municipality, data = gpa2$gdf, iter = i)
anova(municipality.model)

###modularity test
modularity.hypothesis1 <- rep("toe", gpa2$landmark.number)
modularity.hypothesis1[c(10:19, 32:41)] <- "lamellae"
m1 <- modularity.test(gpa2$gdf$coords, modularity.hypothesis1, CI=TRUE, iter=99)
m1

##only landmarks for top section of toepads (Ryan)

###data cleaning
####subset by different number of landmarks
subset_toepads3 <- subsetgmm(toepads, landmarks = c(3,30,33,35,18,14,10,41,7,45,48,51,9,52,55,58,8,62,11,15,19,70,73,4,118,121,102,105,86,89))
subset_toepads3$landmark.number <- 30
landmark.plot(subset_toepads3)
####create links matrix
links3 <- matrix(c(1,2, 2,3, 3,4, 4,5, 5,6, 6,7, 7,8, 8,9, 9,10, 10,11, 11,12, 12,13, 13,14, 14,15, 15,16, 16,17, 17,18, 18,19, 19,20, 20,21, 21,22, 22,23, 23,24, 5,25, 25,26, 26,21, 6,27, 27,28, 28,20, 7,29, 29,30, 30,19),
                  ncol = 2, byrow = T)
landmark.plot(subset_toepads3, links = links3)

###gpa analysis
gpa3 <- align.procrustes(subset_toepads3, outlier.analysis = TRUE)

###pca/shape space
pca3 <- gm.prcomp(gpa3$gdf$coords)
shape.space(pca3) #PC1 = 31.55%, PC2 = 20.04%
####by habitat
shape.space(pca3, group = toepads$metadata$Habitat_Category)
shape.space(pca3, group = toepads$metadata$Habitat_Category, convex.hulls = TRUE,
            include.legend = TRUE)
###by municipalities
shape.space(pca3, group = toepads$metadata$Municipality)
shape.space(pca3, group = toepads$metadata$Municipality, convex.hulls = TRUE, 
            include.legend = TRUE)

###regression (Csize + urban/forest or muncipilities)
####null model
i <- 1e4-1
size.model <- procD.lm(coords~log(Csize), data = gpa3$gdf, iter = i)
anova(size.model)
####habitat_category model
gpa3$gdf$habitat_category <- toepads$metadata$Habitat_Category
habitat.model <- procD.lm(coords~log(Csize) + habitat_category, data = gpa3$gdf, iter = i)
anova(habitat.model)
####municipalities model
gpa3$gdf$municipality <- toepads$metadata$Municipality
municipality.model <- procD.lm(coords~log(Csize) + municipality, data = gpa3$gdf, iter = i)
anova(municipality.model)

###modularity test
modularity.hypothesis2 <- rep("toe", gpa3$landmark.number)
modularity.hypothesis2[c(5,6,7,21,20,19,25,27,29,26,28,30)] <- "lamellae"
m2 <- modularity.test(gpa3$gdf$coords, modularity.hypothesis2, CI=TRUE, iter=99)
m2


##only landmarks for lamellae (Lexi)

###data cleaning
####subset by different number of landmarks
subset_toepads4 <- subsetgmm(toepads, landmarks = c(10:19, 118,110,102,94,86,121,113,105,97,89))
subset_toepads4$landmark.number <- 20
landmark.plot(subset_toepads4)
####create links matrix
links4 <- matrix(c(1,15, 15,20, 20,2, 3,14, 14,19, 19,4, 5,13, 13,18, 18,6, 7,12, 12,17, 17,8, 9,11, 11,16, 16,10, 1,3, 3,5, 5,7, 7,9, 2,4, 4,6, 6,8, 8,10),
                              ncol = 2, byrow = TRUE)
landmark.plot(subset_toepads4, links = links4)

###gpa analysis
gpa4 <- align.procrustes(subset_toepads4, outlier.analysis = TRUE)

###pca/shape space
pca4 <- gm.prcomp(gpa4$gdf$coords)
shape.space(pca4) #PC1 = 42.05%, PC2 = 29.01%
####by habitat
shape.space(pca4, group = toepads$metadata$Habitat_Category)
shape.space(pca4, group = toepads$metadata$Habitat_Category, convex.hulls = TRUE,
            include.legend = TRUE)
###by municipalities
shape.space(pca4, group = toepads$metadata$Municipality)
shape.space(pca4, group = toepads$metadata$Municipality, convex.hulls = TRUE, 
            include.legend = TRUE)

###regression (Csize + urban/forest or muncipilities)
####null model
i <- 1e4-1
size.model <- procD.lm(coords~log(Csize), data = gpa4$gdf, iter = i)
anova(size.model)
####habitat_category model
gpa4$gdf$habitat_category <- toepads$metadata$Habitat_Category
habitat.model <- procD.lm(coords~log(Csize) + habitat_category, data = gpa4$gdf, iter = i)
anova(habitat.model)
####municipalities model
gpa4$gdf$municipality <- toepads$metadata$Municipality
municipality.model <- procD.lm(coords~log(Csize) + municipality, data = gpa4$gdf, iter = i)
anova(municipality.model)

###modularity test
modularity.hypothesis3 <- rep("lamellae", gpa4$landmark.number)
modularity.hypothesis3[c(1, 3, 5, 7, 9, 2, 4, 6, 8, 10)] <- "edge"
m3 <- modularity.test(gpa4$gdf$coords, modularity.hypothesis3, CI=TRUE, iter=99)
m3
