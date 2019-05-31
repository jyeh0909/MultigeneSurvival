#' @title Multiple Gene Survival Analysis
#' @description
#' This function used LASSO/Elastic net in cox regreesion to determine the important genes and used Kmeans (K=2) to divide them into 2 groups in survival.
#' The output includes K-M plot, a list of important genes selected from LASSO/Elastic in Cox model.
#' @details
#' 1. For each gene, the gene expression will be standardized to be mean equal to 0 and standard deviation equal to 1. If we don’t do so, we couldn’t compare the coefficient estimates across genes.
#'
#' 2. Use LASSO/Elastic method for variables selection in the Cox model. The predictor is gene and the outcome is survival. For those coefficient converging to 0 in Cox model are not related to survival, while the remaining nonzero coefficients are chosen for the next step.
#'
#' 3.	To use LASSO/Elastic, we have to determine tuning parameter λ. The cross validation for 100 times was used to determine the most appropriate λ so that the cross validation error is minimized.
#'
#' 4.	After choosing important genes, i.e. coefficients are not equal to 0, we will use standardized gene expression to separate subjects into two groups based on k-means (k=2 here).
#'
#' 5.	Plot the Kaplan-Meier curve from the k-means cluster result.
#'
#' @param dataset A dataset. This dataset contains only gene expression and survival information.
#' @param time A character. The time variable name in dataset.
#' @param event A character. The event variable name in dataset.
#' @param title A character. Usually it's data set name.
#' @param alpha A number. alpha = 1 is LASSO, 0 < alpha < 1 is Elastic net.
#'
#' @examples
#' \dontrun{
#' data(chin_expr_surv)
#' MultigeneSurvival(dataset=chin_expr_surv,time="time",event="event",title="Chin",alpha=1)
#' }
#'
#' @export
#' @import survival
#' @import glmnet
#' @import ggplot2
#' @import survminer

MultigeneSurvival <- function(dataset,time,event,title,alpha){

  a <- which(colnames(dataset) == time)
  b <- which(colnames(dataset) == event)

  colnames(dataset)[a] <- "time"
  colnames(dataset)[b] <- "event"

  gene <- dataset[,-c(a,b)]
  surv <- dataset[,c(a,b)]

  ### Standardized with mean = 0 and SD = 1
  gene1 <- scale(gene,scale = T,center = T)
  y <- Surv(surv$time, surv$event)

  alpha <- alpha

  ### Using 10 times 10 folds Cross validation procedure to choose the tuning parameter when cross validation error is minimized.
  lambda <- NULL

  for(j in 1:10){
    cvfit <- cv.glmnet(gene1, y, family="cox", alpha = alpha, nlambda = 1000)
    lambda[j] <- cvfit$lambda.min
  }

  ### Average the 10 lambdas to be the new tuning parameter for LASSO in cox regression
  lambda1 <- mean(lambda)

  ### LASSO in cox regression
  cvfit1 <- glmnet(gene1, y, family="cox", alpha = alpha, lambda = lambda1)

  ### A list of all coefficients corresponding to each gene
  coef.min = data.frame(as.matrix(coef(cvfit1, s = lambda1)))

  ### Select the gene when its coefficient after LASSO is not zero
  gene.min = rownames(coef.min)[which(coef.min != 0)]

  ### combine the gene selected and survival information to be a new dataset
  dataused1 <- cbind(gene1[,which(colnames(gene1) %in% gene.min)],surv)

  ### do kmeans based on these selected gene expression
  kmeans <- kmeans(dataused1[,-c((ncol(dataused1)-1),ncol(dataused1))], 2 , iter.max = 10000000)
  dataused1$cluster <- kmeans$cluster

  ### Do survival curves based on kmeans cluster
  fit <- survfit(Surv(time,event) ~ cluster, data=dataused1)
  title <- paste0("Dataset = ",title," / ","Alpha = ", alpha)

  ### KM plot
  g <- ggsurvplot(
        fit,
        data=dataused1,
        risk.table=TRUE,
        pval=TRUE,
        conf.int=TRUE,
        ggtheme=theme_minimal(),
        risk.table.y.text.col=TRUE,
        risk.table.y.text=FALSE,
        title = title)

 return(list(g,gene.min))
}
