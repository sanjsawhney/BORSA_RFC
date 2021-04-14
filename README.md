# BORSA_RFC
Code relevant to the Random Forest Classifier for BORSA identification in Sawhney &amp; Ransom et al. (2021).

## Random Forest Classifier
### STEP 1: FEATURE ELIMINATION (caret)
<pre>
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
  ##run the RFE algorithm
  #results <- rfe(reduced_data[,2:129], reduced_data[,1], sizes=c(1:128), rfeControl=control)
  results <- rfe(x=reduced_data[,2:73], y=reduced_data[,1], sizes=c(1:72), rfeControl=control)
  feature_elimination_accuracy[,i]=results$results$Accuracy
}

##summarize the results, list the minimum number of features required, and plot the results
print(results)
predictors(results)
write.csv(feature_elimination_accuracy,"feature_elimination_acc_corr.csv")
</pre>
### STEP 2: RANDOM FOREST CLASSIFIER
<pre>
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
#### Create for loop to iterate the RFC 100x
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
  #rf_prediction_class <- predict(mutation_rf, newdata=testingData, type = "class")
  #table(observed=testingData[,1],predicted=rf_prediction_class)
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
</pre>

### STEP 3: IMPORTANCE PLOT
<pre>
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
</pre>

#### EXTRA CODE - Tuning RFC
<pre>
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
</pre>
