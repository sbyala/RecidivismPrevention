library(ROCR)
library(caret)
library(gains)
library(ggplot2)
library(lattice)
library(corrplot)
#Logistic Regression Model
#load data
mydata <- read.csv("C:/Users/shikh/OneDrive/Desktop/Assignments/BDA/arrest data -1.csv")
mydata

#converting integers into factors
mydata$age_fac <- as.factor(mydata$age_fac)
mydata$level <- as.factor(mydata$level)


# Remove missing values
mydata[is.na(mydata)] <- 0

# Extract numeric columns
num_cols <- sapply(mydata, is.numeric)
mydata_num <- mydata[, num_cols]

# Calculate correlation matrix
cor_matrix <- cor(mydata_num)

# Display correlation matrix
print(cor_matrix)
# Display correlation matrix as simple plot
corrplot(cor_matrix)


# Set the seed for reproducibility
set.seed(111)

#random selection of 60% data
trainindex<-sample(c(1:dim(mydata)[1]),dim(mydata)[1]*0.6) 
trainindex

#Split data into 60% and 40% 
train.df <-mydata[trainindex,] 
valid.df<-mydata[-trainindex,] 

#Logistic Regression
glm.out=glm(arrest ~ fin + race + age_fac+ wexp + prio + level + emp1 + emp2 + emp3 + emp4 + emp5 + emp6 + emp7 + emp8 + emp9 + emp10 + emp11 + emp12 + emp13 + emp14 + emp15 + emp16 + emp17 + emp18 + emp19 + emp20 + emp21 + emp22 + emp23 + emp24 + emp25 + emp26 + emp27 + emp28 + emp29 + emp30 + emp31 + emp32 + emp33 + emp34 + emp35 + emp36 + emp37 + emp38 + emp39 + emp40 + emp41 + emp42 + emp43 + emp44 + emp45 + emp46 + emp47 + emp48 + emp49 + emp50 + emp51 + emp52, family = binomial(logit), data = train.df)
summary(glm.out)


# draw the confusion matrix on validation set
pred <- predict(glm.out,valid.df,type='response')            
valid.df <- cbind(valid.df,pred)
pred2 <- ifelse(valid.df$pred > 0.5,1,0)
valid.df <- cbind(valid.df,pred2)
confusionMatrix(data = factor(valid.df$pred2),reference = factor(valid.df$arrest), positive = '1')
library(gains)
# create the lift table
gain1 <- gains(valid.df$arrest,valid.df$pred,groups=10)
gain1
valid.df$pred
sort(valid.df$pred,decreasing = TRUE)[103]

# draw the lift curve
plot(c(0, gain1$cume.pct.of.total*sum(valid.df$arrest)) ~ c(0, gain1$cume.obs),type='l',
     xlab = 'number of arrests',ylab = 'cumulative',col="pink", lty=1, main = "Lift chart")


# load necessary packages
library(ggplot2)

# create data frame with lift data
lift_data <- data.frame(cume.obs = gain1$cume.obs, cume.pct.of.total = gain1$cume.pct.of.total)

# create lift curve using ggplot2
ggplot(lift_data, aes(x = cume.obs, y = cume.pct.of.total)) +
  geom_line(color = "pink", size = 1.2, linetype = "solid") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray") +
  scale_x_continuous(expand = c(0,0), limits = c(0, max(gain1$cume.obs))) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.1)) +
  labs(title = "Lift Chart",
       x = "Number of Arrests",
       y = "Cumulative Percentage of Total") +
  theme_minimal()

# draw the decile lift curve
barplot(gain1$mean.resp / mean(valid.df$arrest), names.arg = gain1$depth,xlab = 'percentile',
        ylab= 'mean response',col="pink", main = "Decile chart")



# create a vector of colors for the bars
colors <- c("pink", "deeppink", "hotpink", "pink", "deeppink", "hotpink",
            "pink", "deeppink", "hotpink", "pink")

# draw the decile lift curve with colored bars
barplot(gain1$mean.resp / mean(valid.df$arrest), names.arg = gain1$depth, xlab = 'percentile',
        ylab= 'mean response', col=colors, main = "Decile chart",
        border = NA, ylim = c(0, max(gain1$mean.resp / mean(valid.df$arrest)) * 1.2))

# add a horizontal line at y=1 to indicate baseline
abline(h = 1, col = "gray", lty = 2)

# add text labels above each bar showing the lift value
text(1:length(gain1$depth), (gain1$mean.resp / mean(valid.df$arrest)) + 0.1, 
     round(gain1$mean.resp / mean(valid.df$arrest), 2), pos = 3)




# generate predictions using threshold of 0.0014124
pred <- ifelse(predict(glm.out, valid.df, type = 'response') > 0.0014124, 1, 0)

# create confusion matrix with the new threshold
confusionMatrix(data = factor(pred), reference = factor(valid.df$arrest), positive = '1')


if (!requireNamespace("pROC", quietly = TRUE)) {
  install.packages("pROC")
}
library(ROCR)

pred_prob <- predict(glm.out, valid.df, type = 'response')
pred_obj <- prediction(pred_prob, valid.df$arrest)

# ROC curve
roc_obj <- performance(pred_obj, measure = "tpr", x.measure = "fpr")
plot(roc_obj, main = "ROC Curve",col="pink", lty=1)


# compute AUC
auc_obj <- performance(pred_obj, measure = "auc")
auc_val <- auc_obj@y.values[[1]]

# plot AUC curve
plot(roc_obj, main = paste("AUC Curve (AUC =", round(auc_val, 3), ")"),
     col = "pink", lty = 1)

# draw AUC line
lines(x = c(0, 1), y = c(0, 1), col = "gray")
text(x = 0.5, y = 0.5, labels = paste("AUC =", round(auc_val, 3)), cex = 1.5, col = "black")


## SURVIVAL ANALYSIS
#Exponential
library(SurvRegCensCov)
library("survminer")
#load data
mydata <- read.csv("C:/Users/shikh/OneDrive/Desktop/Assignments/BDA/arrest data -1.csv")
mydata
library(survival)

names(mydata)

attach(mydata)
#FOR EXP MODEL
exp.model <- survreg(Surv(week,arrest) ~ fin, data=mydata, dist='exponential')
exp.model
summary(exp.model)
exp(-(summary(exp.model)$table[,1]))

#FOR WEIBULL
wb.model <- survreg(Surv(week,arrest) ~ fin, data=mydata, dist='weibull')
wb.model
summary(wb.model)
ConvertWeibull(wb.model)

"To check if the two models are statistically different"
#Log Rank Test is the used to check if there is a significant difference between the two models
" H0 : survival in two groups is same
  Ha : surv not"
survdiff(Surv(week,arrest) ~ fin + race + mar + wexp + paro + prio + educ)
#Based on the result that the p value is less than 5%, we reject the null hypo, meaning that it is statisically signifant and it is imp that survival depends of whether one is over40
