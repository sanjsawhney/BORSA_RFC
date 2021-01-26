# BORSA_RFC
All code used in Sawhney &amp; Ransom et al. (2021), including the Random Forest Classifier for BORSA identification.

## Fig. S1 - Lineage, MLST Barplots=
library(ggplot2)
library(RColorBrewer)

#Read in barplot file
df_core_lineage <- read.csv('r_barplot_coregenome_lineage.csv',
                            sep=",",
                            header = T)

sapply(df_core_lineage, class)
df_core_lineage$Lineage <- as.factor(df_core_lineage$Lineage)

#Make barplot
barplot_core_lineage <- ggplot(data=df_core_lineage, aes(x=Lineage,y=Isolate_Count,fill=factor(Class,levels=c("MRSA","BORSA","MSSA"))))+
  geom_bar(position="stack", stat="identity", width=0.75, color="black")+
  scale_y_continuous(name="Number of Isolates", limits=c(0,50), expand=c(0,0))+
  scale_fill_manual(values=c("#006666","#66FFFF","#6699FF"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), axis.ticks.x=element_blank())+
  xlab("Lineage")+
  ggtitle("Distribution by Lineage")

barplot_core_lineage+labs(fill="Class")

#Read in barplot file
df_core_mlst <- read.csv('r_barplot_NICU_metadata.csv',
                         sep=",",
                         header = T)

#Make barplot
barplot_core_mlst <- ggplot(data=df_core_mlst, aes(x=MLST,y=Isolate_Count,fill=factor(Class,levels=c("MRSA","BORSA","MSSA"))))+
  geom_bar(position="stack", stat="identity", width=0.75, color="black")+
  scale_y_continuous(name="Number of Isolates", limits=c(0,21), expand=c(0,0))+
  scale_fill_manual(values=c("#006666","#66FFFF","#6699FF"))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), axis.ticks.x=element_blank())+
  #xlim("1.1","1.2","1.3","4","5","6")+
  xlab("Lineage")+
  ggtitle("Distribution by Sequence Type")

barplot_core_mlst+labs(fill="Class")

# Fig. S2 - ACE
require(phytools)
library(ape)
#Read in tree
borsa.tree<-read.tree("6I5emvN5wOmcmugzZ77Uew_newick.txt")

#Read in metadata linking isolate to Class
metadata_nicu<-read.csv ("r_NICU_metadata.csv",
                         header=T,
                         stringsAsFactors = F,
                         comment.char = "",
                         row.names=1,
                         quote = "")
#Set class metadata
borsa.class<-setNames(metadata_nicu[,1],rownames(metadata_nicu))
#Set colors for each class
cols<-setNames(c("#66FFFF","black","#6699FF"),levels(borsa.class))
#Estimate ancestral states under ER model. ACE = ancestral character estimation
classER<-ace(borsa.class,borsa.tree,model="ER",type="discrete")
#Plot tree file
plotTree(borsa.tree,type="fan",fsize=0.7, ftype="i")
#Add in 
nodelabels(node=1:borsa.tree$Nnode+Ntip(borsa.tree),
           pie=classER$lik.anc,piecol=cols,cex=0.3)

# Fig. S3 - S. aureus Pangenome size
#!/usr/bin/env Rscript
#Take the output files from the Roary pan genome pipeline and create nice plots.
library(ggplot2)
conserved = colMeans(read.table("number_of_conserved_genes.Rtab"))
total = colMeans(read.table("number_of_genes_in_pan_genome.Rtab"))

genes = data.frame( genes_to_genomes = c(conserved,total),
                    genomes = c(c(1:length(conserved)),c(1:length(conserved))),
                    Key = c(rep("Conserved genes",length(conserved)), rep("Total genes",length(total))) )

ggplot(data = genes, aes(x = genomes, y = genes_to_genomes, group = Key, linetype=Key)) +geom_line()+
  theme_classic() +
  ylim(c(1,max(total)))+
  xlim(c(1,length(total)))+
  xlab("No. of genomes") +
  ylab("No. of genes")+ theme_bw(base_size = 16) + theme(legend.position="bottom")

# Fig. 2B - Accessory_Genome_PCoA

##For inclusion of MRSA isolates: replace accessorygenome_df1.csv and accessorygenome_df2.csv with accessorygenome_df1_mrsa.csv and accessorygenome_df2_mrsa.csv

#load vegan for jaccard distance, ape for pcoa, ggplot2 for visualization
library(ape)
library(vegan)
library(ggplot2)

#Read in accessorygenome.csv (.Rtab file)
df_accessorygenome<-read.csv('accessorygenome_df1.csv',
                             sep=",",
                             header = T,
                             row.names = 1)

#Calculate jaccard distance. Set to matrix
jaccard_accessorygenome<-vegdist(df_accessorygenome, method='jaccard',binary=TRUE)
jaccard_accessorygenome<-as.matrix(jaccard_accessorygenome)

#Make PCoA (correction is to account for negative eigenvalues) - can use "lingoes" or "cailliez"
#Look through pcoa_accessorygenome_corr on first use to gain an understanding of what $values, $vectors, and $vectors.cor mean
#pcoa_accessorygenome<-pcoa(jaccard_accessorygenome)
pcoa_accessorygenome_corr<-pcoa(jaccard_accessorygenome, correction = "lingoes")

#Get PCoA vectors to plot ordination as axis x axis. Set to data frame
pcoavectors_accessorygenome_corr<-pcoa_accessorygenome_corr$vectors.cor
pcoavectors_accessorygenome_corr<-as.data.frame(pcoavectors_accessorygenome_corr)

#Get % variance captured by each axis. Rel_corr_eig = Relative eigenvalues following correction method. Sums to 1.0000
rel_eigen_accessorygenome_corr<-pcoa_accessorygenome_corr$values$Rel_corr_eig
rel_eigen_accessorygenome_corr

#Scree plot
#barplot(pca_accessorygenome$values$Relative_eig[1:10])
barplot(pcoa_accessorygenome_corr$values$Rel_corr_eig[1:10])
#biplot.pcoa(pca_accessorygenome_corr, df_accessorygenome)

#Read in metadata.csv
accessorygenome_df2<-read.csv('accessorygenome_df2.csv',
                              sep=",",
                              header = T,
                              row.names = 1)

#plot first two principal coordinate axes
ggplot(
  pcoavectors_accessorygenome_corr
  ,aes(x=Axis.1, y=Axis.2))+
  geom_point(color = "black", size = 2.5, show.legend = TRUE)+
  scale_color_manual(values=c(MSSA="#6699FF",BORSA="#66FFFF",MRSA="#CC99FF"))+
  geom_point(aes(color=accessorygenome_df2$Class))+
  theme_bw()+
  theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black"), panel.grid.major = element_blank(), legend.title=element_blank(), plot.title = element_text(hjust = 0.5))+
  ggtitle("Accessory Gene Content")+xlab("PCA1 (24.5%)")+ylab("PCA2 (8.6%)")+
  coord_fixed((rel_eigen_accessorygenome_corr[2]/rel_eigen_accessorygenome_corr[1]))


#Adonis test for significance in clustering (permANOVA)
adonis(jaccard_accessorygenome~accessorygenome_df2$Class,pcoavectors_accessorygenome_corr[,c(1,2)],permutations=1000)


# Fig. 2C - ARG_Heatmap
library(reshape2)

resistance_df<-read.csv('r_resistance_genes.csv',
                        sep=",",
                        header = T)

pres.abs_df<-dcast(resistance_df, Isolate~resistance_df$Gene.symbol, length)
write.csv(pres.abs_df, file='resistance_presence_absence.csv')

anno_colors = list(
  Class = c(MSSA="#6699FF",BORSA="#66FFFF",MRSA="#006666")
)

library(dendsort)
library(pheatmap)

heat_pres.abs_df<-read.csv('resistance_presence_absence.csv', head = T, row.names=1, check.names=F)
heat_pres.abs_df
metadata_nicu<-read.csv ("r_NICU_metadata.csv",
                         header=T,
                         stringsAsFactors = F,
                         comment.char = "",
                         row.names=1,
                         quote = "")

metadata_drug<-read.csv ("r_drug_metadata.csv",
                         header=T,
                         stringsAsFactors = F,
                         comment.char = "",
                         row.names=1,
                         quote = "")


anno_colors = list(
  Class = c(BORSA="#66FFFF",MSSA="#6699FF",MRSA="#336666"),
  #Location = c(Anterior_Nares_NICU="#660000",Blood_Adult_Child ="#CC9999", Blood_NICU="#FFFFFF"),
  MLST = c(ST1="red",ST5="yellow",ST8="#996600",ST15="#FF33FF",ST30="#006600",ST45="purple",ST59="#66CCCC",ST72="#FFCCFF",ST88="#FF9999",ST97="orange",ST398="#0000FF",Other="#CCCCCC"))

#set callback function for sorting dendrogram
callback = function(hc, ...) {
  dendsort(hc)
}

drows=dist(heat_pres.abs_df, method = 'binary')

pheatmap(t(heat_pres.abs_df),
         cluster_rows = F,
         cluster_cols = T,
         #clustering_distance_cols=drows,
         #clustering_method = 'complete',
         clustering_callback = callback,
         cellheight = 10,
         cellwidth = 8,
         gaps_row = c(4,6,7,8,9,10,16,17,19,20),
         color=c('white','#666666'),
         annotation_colors = anno_colors,
         annotation_col = metadata_nicu,
         annotation_row = metadata_drug,
         border_color = 'black',
         fontsize_col = 6,
         fontsize_row = 7,
         angle_col =270,
         legend = F)

# Fig. 3A - Protein_Mutation
library(reshape2)
library(dendsort)
library(pheatmap)

#Read in presence-absence (0s,1s) data frame
mutation_pres.abs_df<-read.csv('Mutation_metadata.csv', head = T, row.names=1, check.names=F)

#Read in metadata linking isolate to Class
metadata_nicu<-read.csv ("r_NICU_metadata.csv",
                         header=T,
                         stringsAsFactors = F,
                         comment.char = "",
                         row.names=1,
                         quote = "")

#Read in metadata linking mutation to Gene
metadata_mutation<-read.csv ("r_mutation_metadata_expanded.csv",
                             header=T,
                             stringsAsFactors = F,
                             comment.char = "",
                             row.names=1,
                             quote = "")

anno_colors = list(
  Class = c(BORSA="#66FFFF",MSSA="#6699FF"),
  Protein = c(PBP1 = "#FF99FF",PBP2 = "#99FFCC",PBP3 = "#FF9966",PBP4 = "#CC66FF",GdpP = "#CC9966"),
  #Location = c(Anterior_Nares_NICU="#660000",Blood_Adult_Child ="#CC9999", Blood_NICU="#FFFFFF"),
  #MLST = c(ST1="red",ST5="yellow",ST8="#996600",ST15="#FF33FF",ST30="#006600",ST45="purple",ST59="#66CCCC",ST72="#FFCCFF",ST88="#FF9999",ST97="orange",ST398="#0000FF",Other="#CCCCCC")
  BORSA_Only = c(Yes = "#666600", No = "#FFFFFF"),
  Unreported = c(Yes = "#999900", No = "#FFFFFF")
)


#set callback function for sorting dendrogram
callback = function(hc, ...) {
  dendsort(hc)
}

#Non-callback/dendsort based clustering scheme (uses clustering_distance_cols and clustering_method)
#drows=dist(mutation_pres.abs_df, method = 'binary')

pheatmap(t(mutation_pres.abs_df),
         cluster_rows = F,
         cluster_cols = T,
         # clustering_distance_cols=drows,clustering_method = 'ward.D',
         clustering_callback = callback,
         cellheight = 7,
         cellwidth = 7,
         gaps_row = c(24,45,65,82),
         #gaps_row = c(31,60,86,106),
         color=c('white','#666666'),
         annotation_colors = anno_colors,
         annotation_col = metadata_nicu,
         annotation_row = metadata_mutation,
         border_color = 'black',
         fontsize_col = 4,
         fontsize_row = 5,
         angle_col =270,
         legend = F)
# Fig. 3B - Venn Diagram (mutations)
library(eulerr)
venn_mutations<-c("BORSA (Previously reported)"=73,"BORSA (This study)"=15,"MSSA (This study)"=26,"BORSA (Previously reported)&BORSA (This study)"=1,"BORSA (Previously reported)&MSSA (This study)"=5,"BORSA (This study)&MSSA (This study)"=26,"BORSA (Previously reported)&BORSA (This study)&MSSA (This study)"=19)
venn_colors<-c("gray","#66FFFF","#6699FF")
#plot(venn(venn_mutations), fills = venn_colors, legend = list(side="right"))
plot(euler(venn_mutations, shape="ellipse"),quantities=TRUE, fills = venn_colors, legend = list(side="right"))

# pvalue
#Fisher's exact test with BH correction for mutation association by BORSA or MSSA status

#_AAC refers to the combined datasets of our work and Argudín et al. (2018)
alltables_AAC <- read.csv('r_pvalue_AAC.csv',
                            sep=",",
                            header = T,
                            row.names = 1)
pvalues_AAC<-apply(alltables_AAC,1,function(x) fisher.test(matrix(x,nr=2))$p.value)
pvalues_AAC
p.adjust(pvalues_AAC, "BH")
#gdpP_premature_stop p=0.001729039
#GdpP I52V p=0.025245965

#_noAAC refers to our dataset alone
alltables_noAAC <- read.csv('r_pvalue_noAAC.csv',
                      sep=",",
                      header = T,
                      row.names = 1)
pvalues_noAAC<-apply(alltables_noAAC,1,function(x) fisher.test(matrix(x,nr=2))$p.value)
p.adjust(pvalues_noAAC, "BH")

job1 = matrix(c(3,5,0,6,1,2,0,1,4,3,5,3,6,6,11,2,14,3,4,10,0,1),11,2)
job1 = matrix(c(3,3,5,0,6,1,2,0,1,4,3,5,8,2,6,5,11,2,10,3,4,7,0,1),12,2)
job1 = matrix(c(7,3,5,0,6,2,1,4,5,13,2,6,5,11,10,4,7,1),9,2)
job2 = matrix(c(2,4,20,1,6,10,7,20,5,11),5,2)

job1
job2
fisher.test(job2,simulate.p.value = TRUE, B=1e7)


# Fig. 4 - Random Forest Classifier
#---STEP 1: FEATURE ELIMINATION (caret)---
##Recursive feature elimination using caret
library(ggplot2)
library(caret)
library(vegan)


#Read in metadata.csv
mutation_all_df<-read.csv('RFC_metadata_noClass.csv',
                          sep=",",
                          header = T,
                          row.names = 1)

#Change mutation 0s and 1s from integers to factors
mutation_all_df[sapply(mutation_all_df, is.integer)] <- lapply(mutation_all_df[sapply(mutation_all_df, is.integer)], as.numeric)
transform(mutation_all_df,clav_effect=as.numeric(clav_effect))
mutation_all_df$MLST=as.numeric(mutation_all_df$MLST)
sapply(mutation_all_df, class)
#calculate correlation matrix
correlationMatrix <- cor(mutation_all_df[,1:128])
#summarize the correlation matrix
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.75)
highlyCorrelated = sort(highlyCorrelated)
#Remove highly correlated variables from dataset. Left with 72 uncorrelated features
reduced_data = mutation_all_df[,-c(highlyCorrelated)]
#Change mutation 0s and 1s from integers to factors
reduced_data[sapply(reduced_data, is.numeric)] <- lapply(reduced_data[sapply(reduced_data, is.numeric)], as.factor)
reduced_data$clav_effect=as.numeric(reduced_data$clav_effect)
write.csv(reduced_data,"reduced_data.csv")
##Add  Class column to reduced_data.csv
reduced_data<-read.csv('reduced_data_corr_withClass.csv',
                       sep=",",
                       header = T,
                       row.names = 1)
#Define the control using a random forest selection function
feature_elimination_accuracy=data.frame(matrix(NA, nrow = 72, ncol = 2))
for (i in 1:25) {
  control <- rfeControl(functions=rfFuncs, method="cv", number=10)
  ## run the RFE algorithm
  #results <- rfe(reduced_data[,2:129], reduced_data[,1], sizes=c(1:128), rfeControl=control)
  results <- rfe(x=reduced_data[,2:73], y=reduced_data[,1], sizes=c(1:72), rfeControl=control)
  feature_elimination_accuracy[,i]=results$results$Accuracy
}
## summarize the results, list the minimum number of features required, and plot the results
print(results)
predictors(results)
write.csv(feature_elimination_accuracy,"feature_elimination_acc_corr.csv")

#---STEP 2: RANDOM FOREST CLASSIFIER---
#Read in metadata.csv. Originally ran on the set of 72 uncorrelated features. (not final model)
##Used the importance function (see below) to determine the 6 most important features for Mean Decrease in Accuracy.
##Re-ran on a new metadata.csv file of just these 6 features (final, sparse model)
mutation_df<-read.csv('RFC_metadata_6_corr_USE_THIS.csv',
                      sep=",",
                      header = T,
                      row.names = 1)

#Change mutation 0s and 1s from integers to factors
mutation_df[sapply(mutation_df, is.integer)] <- lapply(mutation_df[sapply(mutation_df, is.integer)], as.factor)
transform(mutation_df,clav_effect=as.numeric(clav_effect))
#Check class of each column
sapply(mutation_df,class)

#Define empty variables to be filled for each iteration of the RFC, over 100 iterations
importance_plot=data.frame(matrix(NA, nrow = 6, ncol = 100))
#Variable for area under the curve
rfc_auc=c()
##(a1, a2) <--> (pp1, pp2) are (sensitivity, specificity coordinates) for the ROC curve
a1=c()
b1=c()
c1=c()
d1=c()
e1=c()
f1=c()
g1=c()
h1=c()
i1=c()
j1=c()
k1=c()
l1=c()
m1=c()
n1=c()
o1=c()
p1=c()
q1=c()
r1=c()
s1=c()
t1=c()
u1=c()
v1=c()
w1=c()
x1=c()
y1=c()
aa1=c()
bb1=c()
cc1=c()
dd1=c()
ee1=c()
ff1=c()
gg1=c()
hh1=c()
ii1=c()
jj1=c()
kk1=c()
ll1=c()
mm1=c()
nn1=c()
oo1=c()
pp1=c()

a2=c()
b2=c()
c2=c()
d2=c()
e2=c()
f2=c()
g2=c()
h2=c()
i2=c()
j2=c()
k2=c()
l2=c()
m2=c()
n2=c()
o2=c()
p2=c()
q2=c()
r2=c()
s2=c()
t2=c()
u2=c()
v2=c()
w2=c()
x2=c()
y2=c()
aa2=c()
bb2=c()
cc2=c()
dd2=c()
ee2=c()
ff2=c()
gg2=c()
hh2=c()
ii2=c()
jj2=c()
kk2=c()
ll2=c()
mm2=c()
nn2=c()
oo2=c()
pp2=c()

##set.seed(8), mtry=3 - AUC: 0.8862 +/- 0.0080 (6_corr)
set.seed(8)
## Create for loop to iterate the RFC 100x
for (i in 1:100) {
  
  input_BORSA <- mutation_df[which(mutation_df$Class == 'BORSA'),]
  input_MSSA <- mutation_df[which(mutation_df$Class == 'MSSA'),]
  
  #set.seed(146) (Fig. 4C)
  #Split isolates into 2 sets that account for class biases (1 to train, 1 to test)
  input_BORSA_training_rows <- sample(1:nrow(input_BORSA), 0.7*nrow(input_BORSA))
  input_MSSA_training_rows <- sample(1:nrow(input_MSSA), 0.7*nrow(input_MSSA))
  
  #Make training dataset
  training_BORSA <- input_BORSA[input_BORSA_training_rows,]
  training_MSSA <- input_MSSA[input_MSSA_training_rows,]
  trainingData <- rbind(training_BORSA,training_MSSA)
  
  #Make test dataset
  testing_BORSA <- input_BORSA[-input_BORSA_training_rows,]
  testing_MSSA <- input_MSSA[-input_MSSA_training_rows,]
  testingData <- rbind(testing_BORSA,testing_MSSA)
  
  #Run randomForest classifier on trainingData
  library(randomForest)
  mutation_rf <- randomForest(Class~., data=trainingData,ntree=5000, proximity=F, importance=T, mtry=3)

  ##Determine most important features for Mean Decrease in Accuracy (Fig. 4B)
  impo<-importance(mutation_rf)
  importance_plot[,i]=impo[,3]
  
  #Run randomForest classifier on testData - Remove hashtags only when determining confusion matrix accuracy (Fig. 4C). Is not iterated.
  # rf_prediction_class <- predict(mutation_rf, newdata=testingData, type = "class")
  # table(observed=testingData[,1],predicted=rf_prediction_class)
  #mean(rf_prediction_class == testingData$Class)

  #Make ROC Curve of VALIDATION dataset (pROC) for reach iteration of RFC (Fig. 4A)
  library(pROC)
  rf_prediction_prob=predict(mutation_rf,type="prob",newdata = testingData)
  testing.pROC<-roc(testingData$Class, rf_prediction_prob[,2] )
  plot(testing.pROC,main="ROC Curve", col=c("light gray"), lwd=0.3, add=TRUE)
  
  #Gather sensitivity and specificity coordinates for each RFC iteration. This will be used for the mean ROC Curve over all 100 iterations.
  a1[i]=testing.pROC$sensitivities[1]
  b1[i]=testing.pROC$sensitivities[2]
  c1[i]=testing.pROC$sensitivities[3]
  d1[i]=testing.pROC$sensitivities[4]
  e1[i]=testing.pROC$sensitivities[5]
  f1[i]=testing.pROC$sensitivities[6]
  g1[i]=testing.pROC$sensitivities[7]
  h1[i]=testing.pROC$sensitivities[8]
  i1[i]=testing.pROC$sensitivities[9]
  j1[i]=testing.pROC$sensitivities[10]
  k1[i]=testing.pROC$sensitivities[11]
  l1[i]=testing.pROC$sensitivities[12]
  m1[i]=testing.pROC$sensitivities[13]
  n1[i]=testing.pROC$sensitivities[14]
  o1[i]=testing.pROC$sensitivities[15]
  p1[i]=testing.pROC$sensitivities[16]
  q1[i]=testing.pROC$sensitivities[17]
  r1[i]=testing.pROC$sensitivities[18]
  s1[i]=testing.pROC$sensitivities[19]
  t1[i]=testing.pROC$sensitivities[20]
  u1[i]=testing.pROC$sensitivities[21]
  v1[i]=testing.pROC$sensitivities[22]
  w1[i]=testing.pROC$sensitivities[23]
  x1[i]=testing.pROC$sensitivities[24]
  y1[i]=testing.pROC$sensitivities[25]
  aa1[i]=testing.pROC$sensitivities[26]
  bb1[i]=testing.pROC$sensitivities[27]
  cc1[i]=testing.pROC$sensitivities[28]
  dd1[i]=testing.pROC$sensitivities[29]
  ee1[i]=testing.pROC$sensitivities[30]
  ff1[i]=testing.pROC$sensitivities[31]
  gg1[i]=testing.pROC$sensitivities[32]
  hh1[i]=testing.pROC$sensitivities[33]
  ii1[i]=testing.pROC$sensitivities[34]
  jj1[i]=testing.pROC$sensitivities[35]
  kk1[i]=testing.pROC$sensitivities[36]
  ll1[i]=testing.pROC$sensitivities[37]
  mm1[i]=testing.pROC$sensitivities[38]
  nn1[i]=testing.pROC$sensitivities[39]
  oo1[i]=testing.pROC$sensitivities[40]
  pp1[i]=testing.pROC$sensitivities[41]
  
  a2[i]=testing.pROC$specificities[1]
  b2[i]=testing.pROC$specificities[2]
  c2[i]=testing.pROC$specificities[3]
  d2[i]=testing.pROC$specificities[4]
  e2[i]=testing.pROC$specificities[5]
  f2[i]=testing.pROC$specificities[6]
  g2[i]=testing.pROC$specificities[7]
  h2[i]=testing.pROC$specificities[8]
  i2[i]=testing.pROC$specificities[9]
  j2[i]=testing.pROC$specificities[10]
  k2[i]=testing.pROC$specificities[11]
  l2[i]=testing.pROC$specificities[12]
  m2[i]=testing.pROC$specificities[13]
  n2[i]=testing.pROC$specificities[14]
  o2[i]=testing.pROC$specificities[15]
  p2[i]=testing.pROC$specificities[16]
  q2[i]=testing.pROC$specificities[17]
  r2[i]=testing.pROC$specificities[18]
  s2[i]=testing.pROC$specificities[19]
  t2[i]=testing.pROC$specificities[20]
  u2[i]=testing.pROC$specificities[21]
  v2[i]=testing.pROC$specificities[22]
  w2[i]=testing.pROC$specificities[23]
  x2[i]=testing.pROC$specificities[24]
  y2[i]=testing.pROC$specificities[25]
  aa2[i]=testing.pROC$specificities[26]
  bb2[i]=testing.pROC$specificities[27]
  cc2[i]=testing.pROC$specificities[28]
  dd2[i]=testing.pROC$specificities[29]
  ee2[i]=testing.pROC$specificities[30]
  ff2[i]=testing.pROC$specificities[31]
  gg2[i]=testing.pROC$specificities[32]
  hh2[i]=testing.pROC$specificities[33]
  ii2[i]=testing.pROC$specificities[34]
  jj2[i]=testing.pROC$specificities[35]
  kk2[i]=testing.pROC$specificities[36]
  ll2[i]=testing.pROC$specificities[37]
  mm2[i]=testing.pROC$specificities[38]
  nn2[i]=testing.pROC$specificities[39]
  oo2[i]=testing.pROC$specificities[40]
  pp2[i]=testing.pROC$specificities[41]
  
  #Determine AUC of validation dataset (ROCR)
  library(ROCR)
  rf_prediction_prob=predict(mutation_rf,type="prob",newdata = testingData)
  testing.ROCR=prediction(rf_prediction_prob[,2],testingData$Class)
  auc.ROCR=performance(testing.ROCR,measure='auc')
  rfc_auc[i]=auc.ROCR@y.values

#End for loop
}

#Get mean and confidence intervals of AUC over 100 iterations
rfc_auc <- as.numeric(rfc_auc)
t.test(rfc_auc)

#Set NAs in sensitivity/specificity coordinates to 0 and 1 respectively
q1[is.na(q1)] <- 0
r1[is.na(r1)] <- 0
s1[is.na(s1)] <- 0
t1[is.na(t1)] <- 0
u1[is.na(u1)] <- 0
v1[is.na(v1)] <- 0
w1[is.na(w1)] <- 0
x1[is.na(x1)] <- 0
y1[is.na(y1)] <- 0
aa1[is.na(aa1)] <- 0
bb1[is.na(bb1)] <- 0
cc1[is.na(cc1)] <- 0
dd1[is.na(dd1)] <- 0
ee1[is.na(ee1)] <- 0
ff1[is.na(ff1)] <- 0
gg1[is.na(gg1)] <- 0
hh1[is.na(hh1)] <- 0
ii1[is.na(ii1)] <- 0
jj1[is.na(jj1)] <- 0
kk1[is.na(kk1)] <- 0
ll1[is.na(ll1)] <- 0
mm1[is.na(mm1)] <- 0
nn1[is.na(nn1)] <- 0
oo1[is.na(oo1)] <- 0
pp1[is.na(pp1)] <- 0

q2[is.na(q2)] <- 1
r2[is.na(r2)] <- 1
s2[is.na(s2)] <- 1
t2[is.na(t2)] <- 1
u2[is.na(u2)] <- 1
v2[is.na(v2)] <- 1
w2[is.na(w2)] <- 1
x2[is.na(x2)] <- 1
y2[is.na(y2)] <- 1
aa2[is.na(aa2)] <- 1
bb2[is.na(bb2)] <- 1
cc2[is.na(cc2)] <- 1
dd2[is.na(dd2)] <- 1
ee2[is.na(ee2)] <- 1
ff2[is.na(ff2)] <- 1
gg2[is.na(gg2)] <- 1
hh2[is.na(hh2)] <- 1
ii2[is.na(ii2)] <- 1
jj2[is.na(jj2)] <- 1
kk2[is.na(kk2)] <- 1
ll2[is.na(ll2)] <- 1
mm2[is.na(mm2)] <- 1
nn2[is.na(nn2)] <- 1
oo2[is.na(oo2)] <- 1
pp2[is.na(pp2)] <- 1

#Find average sensitivity and specificity coordinates to plot average ROC curve over 100 iterations
mean_sensitivity=c(mean(a1), mean(b1), mean(c1), mean(d1), mean(e1), mean(f1), mean(g1), mean(h1), mean(i1), mean(j1), mean(k1), mean(l1), mean(m1), mean(n1), mean(o1), mean(p1), mean(q1), mean(r1), mean(s1), mean(t1), mean(u1), mean(v1), mean(w1), mean(x1), mean(y1), mean(aa1), mean(bb1), mean(cc1), mean(dd1), mean(ee1), mean(ff1), mean(gg1), mean(hh1), mean(ii1), mean(jj1))
#, mean(kk1), mean(ll1), mean(mm1), mean(nn1), mean(oo1), mean(pp1))
mean_specificity=c(mean(a2), mean(b2), mean(c2), mean(d2), mean(e2), mean(f2), mean(g2), mean(h2), mean(i2), mean(j2), mean(k2), mean(l2), mean(m2), mean(n2), mean(o2), mean(p2), mean(q2), mean(r2), mean(s2), mean(t2), mean(u2), mean(v2), mean(w2), mean(x2), mean(y2),  mean(aa2), mean(bb2), mean(cc2), mean(dd2), mean(ee2), mean(ff2), mean(gg2), mean(hh2), mean(ii2), mean(jj2))
#, mean(kk2), mean(ll2), mean(mm2), mean(nn2), mean(oo2), mean(pp2))

#Plot average ROC curve over 100 iterations (Fig. 4A)
segments(1.2,-0.2,-0.2,1.2)
lines(x=mean_specificity, y=mean_sensitivity, col = "red", lwd=2.5)
points(x=mean_specificity, y=mean_sensitivity, pch=16, col="red", cex=0.75)
legend(0.2,0.1,title="AUC = 0.81 ± 0.01", legend=c("ROC Curve","Luck"), col=c("red","black"),lty=1:1, lwd=2.5:0.5, cex=0.8)


#---STEP 3: IMPORTANCE PLOT (Fig. 4B)---
library(ggplot2)
write.csv(importance_plot,"importance_plot_100x_6_corr.csv")
#Add SEM column to importance_plot.csv
imp_df<-read.csv('importance_plot_100x_6_corr.csv',
                 sep=",",
                 header = T)

#Order mutation by csv instead of alphabetical
imp_df$Feature<-factor(imp_df$Feature, levels=imp_df$Feature[order(imp_df$Average)])

ggplot(imp_df, aes(x=Feature, y=Average_Accuracy_Decrease))+
  geom_point(stat='identity')+
  coord_flip()+
  geom_errorbar(aes(ymin=Average_Accuracy_Decrease-SE, ymax=Average_Accuracy_Decrease+SE),width=0.2,position=position_dodge(0.05))+
  theme_bw()+
  scale_y_continuous(limits = c(0,100))+
  ylab("Mean Decrease in Accuracy")+xlab("Feature")


#EXTRA CODE

#---Tuning RFC---
##Determine optimal number of variables to test at each node
a=c()
for (i in 1:72) {
  model3<-randomForest(Class~., data=trainingData, ntree=5000, proximity=F, importance=T, mtry=i)
  pred<-predict(model3, testingData, type="class")
  a[i]=mean(pred==testingData$Class)
}
a

##Determine optimal number of trees in random forest
val = c(10,50,100,200,500,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000)
optimal_tree=c()
for (i in val) {
  model3<-randomForest(Class~., data=trainingData, ntree=i, proximity=F, importance=T, mtry=1)
  pred<-predict(model3, testingData)
  optimal_tree[i]=mean(pred==testingData$Class)
}
optimal_tree <- na.omit(optimal_tree)
optimal_tree
