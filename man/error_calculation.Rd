% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/error_calculation.R
\name{error_calculation}
\alias{error_calculation}
\title{Error calculation for integrated model}
\usage{
error_calculation(final_pred, final_actual)
}
\arguments{
\item{final_pred}{A n x p matrix of predicted features, where n is the number of samples and p is the number of data subtypes with prediction}

\item{final_actual}{A n x 1 vector of original output responses}
}
\value{
List with the following components: 
\item{Integrated Prediction}{Integrated Prediction based on combining predictions from data subtypes using Least Square Regression}
\item{error_mae}{Mean Absolute Error between Integrated Prediction and Original Output Responses}
\item{error_mse}{Mean Square Error between Integrated Prediction and Original Output Responses}
\item{error_corr}{Pearson Correlation Coefficient between Integrated Prediction and Original Output Responses}
}
\description{
Combines Prediction from different data subtypes through Least Square Regression and computes Mean Absolute Error, 
Mean Square Error and Pearson Correlation Coefficient between Integrated Prediction and Original Output feature.
}
\details{
If final_pred is a vector, it refers to the prediction result for one subtype of dataset and this function will return 
Mean Absolute Error, Mean Square Error and Pearson Correlation Coefficient between predicted and Original Output response. 
If final_pred is a matrix containing prediction results for more than one subtype of dataset, Least Square 
Regression will be used to calculate the weights for combining the predictions and generate an integrated prediction of size n x 1. 
Subsequently, Mean Absolute Error, Mean Square Error and Pearson Correlation Coefficient between 
Integrated Prediction and Original Output responses are calculated.
}
\seealso{
\code{lsei}
}

