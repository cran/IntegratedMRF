% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CombPredict.R
\name{CombPredict}
\alias{CombPredict}
\title{Integrated Prediction of Testing samples using Combination Weights from integrated RF or MRF model}
\usage{
CombPredict(finalX, finalY_train, Cell, finalY_train_cell, finalY_test_cell,
  n_tree, m_feature, min_leaf, Coeff)
}
\arguments{
\item{finalX}{List of Matrices where each matrix represents a specific data subtype (such as genomic characterizations for
drug sensitivity prediction). Each subtype can have different types of features. For example, if there are three subtypes containing
 100, 200 and 250 features respectively,  finalX will be a list containing 3 matrices of sizes M x 100, M x 200 and M x 250
 where M is the number of Samples.}

\item{finalY_train}{A M x T matrix of output features for training samples, where M is number of samples and
T is the number of output features. The dataset is assumed to contain no missing values. If there are missing values, an imputation method
should be applied before using the function. A function 'Imputation' is included within the package.}

\item{Cell}{It contains a list of samples (the samples can be represented either numerically by indices or by names) for each data subtype.
For the example of 3 data subtypes, it will be a list containing 3 arrays where each array contains the sample information for each data subtype.}

\item{finalY_train_cell}{Cell lines of output features for training samples}

\item{finalY_test_cell}{Cell lines of output features for testing samples}

\item{n_tree}{number of trees in the forest, which must be positive integer}

\item{m_feature}{Number of randomly selected features considered for a split in each regression tree node, which must be positive integer}

\item{min_leaf}{minimum number of samples in the leaf node, which must be positive integer and less than or equal to M (number of training samples)}

\item{Coeff}{Combination Weights. The user can supply the weights based on either Bootstrap, Re-substitution, 0.632Bootstrap or Leave-one-out
error estimation approaches.}
}
\value{
Final Prediction of testing samples based on provided testing sample names.
}
\description{
Generates Random Forest or Multivariate Random Forest model for each subtype of dataset and predicts testing samples using the generated models.
Subsequently, the prediction for different subtypes of dataset are combined using the Combination weights generated from 'Combination' function.
}
\details{
Input matrix and output response of training samples have been used to build Random Forest or Multivariate Random Forest model for each subtype of
a dataset. These models are used to calculate prediction of
testing samples for each subtype separately. Subsequently Combination Weights (different errors have different combination weights
and the user should select the one to be used) are used to integrate the predictions from data subtypes.
Note that the combination weights are linear regression coefficients generated using the training samples.

The specific set of combination weights to be used for testing samples will depend on the number of data subtypes available
for the testing samples. Note that not all subtype information maybe available for all samples.
As an example with three data subtypes, a testing sample with all subtype data available will use
the combination weights corresponding to Serial [1 2 3] where if subtype 3 is not available, the function will
using the combination weights corresponding to Serial [1 2].
}
\examples{
library(IntegratedMRF)
data(Dream_Dataset)
Tree=1
Feature=1
Leaf=10
Confidence=80
finalX=Dream_Dataset[[1]]
Cell=Dream_Dataset[[2]]
Y_train_Dream=Dream_Dataset[[3]]
Y_train_cell=Dream_Dataset[[4]]
Y_test=Dream_Dataset[[5]]
Y_test_cell=Dream_Dataset[[6]]
Drug=1
Y_train_Drug=matrix(Y_train_Dream[,Drug],ncol=length(Drug))
Result=Combination(finalX,Y_train_Drug,Cell,Y_train_cell,Tree,Feature,Leaf,Confidence)

CombPredict(finalX,Y_train_Drug,Cell,Y_train_cell,Y_test_cell,Tree,Feature,Leaf,Result[[1]])
}

