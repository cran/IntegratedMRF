% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IntegratedMRF.R
\name{Combination}
\alias{Combination}
\title{Weights for combination of predictions from different data subtypes using Least Square Regression based on various error estimation techniques}
\usage{
Combination(finalX, finalY_train, Cell, finalY_train_cell, n_tree, m_feature,
  min_leaf, Confidence_Level)
}
\arguments{
\item{finalX}{List of Matrices where each matrix represent a specific data subtype (such as genomic characterizations for
drug sensitivity prediction). Each subtype can have different types of features. For example, if there are three subtypes containing
 100, 200 and 250 features respectively,  finalX will be a list containing 3 matrices of sizes M x 100, M x 200 and M x 250
 where M is the number of Samples.}

\item{finalY_train}{A M x T matrix of output features for training samples, where M is number of samples and T is the number of output features.
The dataset is assumed to contain no missing values. If there are missing values, an imputation method should be applied before using the function.
A function 'Imputation' is included within the package.}

\item{Cell}{It contains a list of samples (the samples can be represented either numerically by indices or by names) for each data subtype.
For the example of 3 data subtypes, it will be a list containing 3 arrays where each array contains the sample information for each data subtype.}

\item{finalY_train_cell}{Sample names of output features for training samples}

\item{n_tree}{Number of trees in the forest, which must be positive integer}

\item{m_feature}{Number of randomly selected features considered for a split in each regression tree node, Valid Input is a positive integer, which is less than N (which is equal to number of input features for the smallest genomic characterization)}

\item{min_leaf}{Minimum number of samples in the leaf node, which must be positive integer and less than or equal to M (number of training samples)}

\item{Confidence_Level}{Confidence level for calculation of confidence interval (User Defined), which must be between 0 and 100}
}
\value{
List with the following components:
\item{BSP_coeff}{Combination weights using Bootstrap Error Estimation Model, where index is in list format.
If the number of genomic characterizations or subtypes of dataset is 5, there will be 2^5-1=31 list of weights}
\item{Nfold_coeff}{Combination weights using N fold cross validation Error Estimation Model, where index is in list format.
If the number of genomic characterizations or subtypes of dataset is 5, there will be 2^5-1=31 list of weights}
\item{BSP632plus_coeff}{Combination weights using 0.632+ Bootstrap Error Estimation Model, where index is in list format.
If the number of genomic characterizations or subtypes of dataset is 5, there will be 2^5-1=31 list of weights}
\item{LOO_coeff}{Combination weights using Leave-One-Out Error Estimation Model, where index is in list format.
If the number of genomic characterizations or subtypes of dataset is 5, there will be 2^5-1=31 list of weights}
\item{Error}{Matrix of Mean Absolute Error, Mean Square Error and correlation between actual and predicted responses for integrated model based
on Bootstrap, N fold cross validation, 0.632+ Bootstrap and Leave-one-out error estimation sampling techniques for the integrated model
containing all the data subtypes}
\item{Confidence Interval}{Low and High confidence interval for a user defined confidence level for the drug using Jackknife-After-Bootstrap Approach in a list}
\item{BSP_error_all_mae}{Bootstrap Mean Absolute Errors (MAE) for all combinations of the dataset subtypes. Size C x R, where C is the number of
combinations and R is the number of output responses. C is in decreasing order, which means first value is combination of all subtypes
and next ones are in decreasing order. For example, if a dataset has 3 subtypes, then C is equal to 2^3-1=7.  The ordering of C is the combination of
subtypes [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3] }
\item{Nfold_error_all_mae}{N fold cross validation Mean Absolute Errors (MAE) for all combinations of the dataset subtypes. Size C x R, where C is the number of
combinations and R is the number of output responses. C is in decreasing order, which means first value is combination of all subtypes
and next ones are in decreasing order. For example, if a dataset has 3 subtypes, then C is equal to 2^3-1=7.  The ordering of C is the combination of
subtypes [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3] }
\item{BSP632plus_error_all_mae}{0.632+ Bootstrap Mean Absolute Errors (MAE) for all combinations of the dataset subtypes. Size C x R, where C is the number of
combinations and R is the number of output responses. C is in decreasing order, which means first value is combination of all subtypes
and next ones are in decreasing order. For example, if a dataset has 3 subtypes, then C is equal to 2^3-1=7.  The ordering of C is the combination of
subtypes [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3] }
\item{LOO_error_all_mae}{Leave One Out Mean Absolute Errors (MAE) for all combinations of the dataset subtypes. Size C x R, where C is the number of
combinations and R is the number of output responses. C is in decreasing order, which means first value is combination of all subtypes
and next ones are in decreasing order. For example, if a dataset has 3 subtypes, then C is equal to 2^3-1=7.  The ordering of C is the combination of
subtypes [1 2 3], [1 2], [1 3], [2 3], [1], [2], [3] }
The function also returns figures of different error estimations in .tiff format
}
\description{
Calculates combination weights for different subtypes of dataset combinations to generate integrated Random Forest (RF) or Multivariate Random Forest (MRF) model based on different error estimation models such as Bootstrap, 0.632+ Bootstrap, N-fold cross validation or Leave one out.
}
\details{
The function takes all the subtypes of dataset in matrix format and its corresponding sample information.
For calculation purpose, we have considered the data of the samples that are common in all the subtypes and output training responses.
For example, consider a dataset of 3 sub-types with different number of samples and features, with indices of samples in subtype 1, 2, 3  and output feature matrix
is 1:10, 3:15, 5:16 and 5:11 respectively. Thus, features of sample index 5:10 (common to all subtypes and output feature matrix) of all subtypes and output feature
matrix will be selected and considered for all calculations.

For M x N dataset, N number of bootstrap sampling sets are considered. For each bootstrap sampling set and each subtype, a Random Forest (RF)
or, Multivariate Random Forest (MRF) model is generated, which is used for calculating the prediction performance for out-of-bag samples.
The prediction performance for each subtype of the dataset is based on the averaging over different bootstrap training sets.
The combination weights (regression coefficients) for each combination of subtypes are generated using least Square Regression from the
individual subtype predictions and used later to calculate mean absolute error, mean square error and correlation coefficient between
predicted and actual values.

For N-fold cross validation error estimation with M cell lines, N models are generated for each subtype of dataset, where for each partition (M/N)*(N-1) cell
lines are used for training and the remaining cell lines are used to
estimate errors and combination weights for different data subtype combinations.

In 0.632 Bootstrap error estimation, bootstrap and re-substitution error estimates are combined based on
0.632xBootstrap Error + 0.368xRe-substitution Error. While 0.632+ Bootstrap error estimation considers the overfitting of re-substitution error
with no information error rate \eqn{\gamma}. An estimate of \eqn{\gamma} is obtained by permuting the responses \eqn{y[i]} and predictors \eqn{x[j]}.
\deqn{\gamma=sum(sum(error(x[j],y[i]),j=1,m),i=1,m)/m^2}
The relative overfitting rate is defined as \eqn{R=(Bootstrap Error-Resubstitution Error)/(\gamma-Resubstitution Error)} and weight distribution
between bootstrap error and Re-substitution Error is defined as \eqn{w=0.632/(1-0.368*R)}. So, 0.632+ Bootstrap error is equal to
\eqn{(1-w)*Bootstrap Error+w*Resubstitution Error}.
These prediction results are then used to compute the errors and combination weights for different data subtype combinations.

Confidence Interval has been calculated using Jackkniffe-After-Bootstrap Approach and prediction result of bootstrap error estimation.

For leave-one-out error estimation using M cell lines, M models are generated for each subtype of dataset, which are then used to
calculate the errors and combination weights for different data subtype combinations.
}
\examples{
library(IntegratedMRF)
data(Dream_Dataset)
Tree=1
Feature=1
Leaf=5
Confidence=80
finalX=Dream_Dataset[[1]]
Cell=Dream_Dataset[[2]]
Y_train_Dream=Dream_Dataset[[3]]
Y_train_cell=Dream_Dataset[[4]]
Y_test=Dream_Dataset[[5]]
Y_test_cell=Dream_Dataset[[6]]
Drug=c(1,2,3)
Y_train_Drug=matrix(Y_train_Dream[,Drug],ncol=length(Drug))
Result=Combination(finalX,Y_train_Drug,Cell,Y_train_cell,Tree,Feature,Leaf,Confidence)

}
\references{
Wan, Qian, and Ranadip Pal. "An ensemble based top performing approach for NCI-DREAM drug sensitivity prediction challenge." PloS one 9.6 (2014): e101183.

Rahman, Raziur, John Otridge, and Ranadip Pal. "IntegratedMRF: random forest-based framework for integrating prediction from different data types." Bioinformatics (Oxford, England) (2017).

Efron, Bradley, and Robert Tibshirani. "Improvements on cross-validation: the 632+ bootstrap method." Journal of the American Statistical Association 92.438 (1997): 548-560.
}

