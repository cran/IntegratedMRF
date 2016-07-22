#' Prediction for testing samples using specific combination weights from integrated RF or MRF model 
#' 
#' Generates Random Forest (One Output Feature) or Multivariate Random Forest (More than One Output Feature) 
#' model for each subtype of dataset and predicts testing samples using these models. The predictions are 
#' combined using the specific combination weights provided by the user. For the input combination weights, 
#' the testing cell lines should have the subtype data corresponding to the non-zero weight subtypes. 
#'  
#' @param finalX List of Matrices where each matrix represent a specific data subtype (such as genomic characterizations for 
#' drug sensitivity prediction). Each subtype can have different types of features. For example, if there are three subtypes containing
#'  100, 200 and 250 features respectively,  finalX will be a list containing 3 matrices of sizes M x 100, M x 200 and M x 250 
#'  where M is the number of Samples.
#' @param finalY_train A M x T matrix of output features for training samples, where M is number of samples and T is the number of output features. 
#' The dataset is assumed to contain no missing values. If there are missing values, an imputation method should be applied before using the function. 
#' A function 'Imputation' is included within the package.
#' @param Cell It contains a list of samples (the samples can be represented either numerically by indices or by names) for each data subtype. 
#' For the example of 3 data subtypes, it will be a list containing 3 arrays where each array contains the sample information for each data subtype.
#' @param finalY_train_cell Sample names of output features for training samples
#' @param finalY_test_cell Sample names of output features for testing samples (All these testing samples
#' must have features for each subtype of dataset)
#' @param n_tree Number of trees in the forest, which must be a positive integer
#' @param m_feature Number of randomly selected features considered for a split in each regression tree node, which must be a positive integer
#' @param min_leaf Minimum number of samples in the leaf node, which must be a positive integer less than or equal to M (number of training samples)  
#' @param Coeff Combination Weights (user defined or some combination weights generated using the 'Combination' function). 
#' The size must be C, which is equal to the number of subtypes of dataset given in finalX.
#' 
#' @return Final Prediction of testing samples based on provided testing sample names
#' @examples
#' library(IntegratedMRF)
#' data(Dream_Dataset)
#' Tree=1
#' Feature=1
#' Leaf=10
#' Confidence=80
#' finalX=Dream_Dataset[[1]]
#' Cell=Dream_Dataset[[2]]
#' Y_train_Dream=Dream_Dataset[[3]]
#' Y_train_cell=Dream_Dataset[[4]]
#' Y_test=Dream_Dataset[[5]]
#' Y_test_cell=Dream_Dataset[[6]]
#' Drug=c(1,2,3)
#' Y_train_Drug=matrix(Y_train_Dream[,Drug],ncol=length(Drug))
#' Result=Combination(finalX,Y_train_Drug,Cell,Y_train_cell,Tree,Feature,Leaf,Confidence)
#' 
#' CombPredictSpecific(finalX,Y_train_Drug,Cell,Y_train_cell,Y_test_cell,Tree,
#'         Feature,Leaf,runif(length(Cell)*1))
#' @details
#' Input feature matrix and output feature matrix have been used to generate Random Forest (One Output Feature) 
#' or Multivariate Random Forest (More than One Output Feature) model for each subtype of dataset separately. 
#' The prediction of testing samples using each subtype trained model is generated. The predictions are combined 
#' using the specific combination weights provided by the user. For the input combination weights, the testing cell lines 
#' should have the subtype data corresponding to the non-zero weight subtypes. For instance, if combination weights are 
#' [0.6 0.3 0 0.1], subtypes 1, 2 and 4 needs to be present for the testing samples. Furthermore, all the features 
#' should be present for the required subtypes for the testing samples.
#' @importFrom caTools combs
#' @export
CombPredictSpecific <- function(finalX,finalY_train,Cell,finalY_train_cell,finalY_test_cell,n_tree,m_feature,min_leaf,Coeff){
  Serial=NULL
  #   library(caTools)
  for (p in length(Cell):1){
    nk=combs(1:length(Cell),p)
    sk=length(Serial)
    for (q in 1:dim(nk)[1]){
      Serial[[sk+q]]=nk[q, ]
    }
  }
  ##
  Common_cell_dataset=NULL
  for (q in 1:length(Cell)){
    Common_cell_dataset[[q]]=intersect(finalY_train_cell,Cell[[q]])
  }
  Variable_number=ncol(finalY_train)
  if (ncol(finalY_train)==1){
    Command=1
  }else if (ncol(finalY_train)>1){
    Command=2
  }
  
  final_dataset=NULL
  finalY_dataset=NULL
  for (q in 1:length(Cell)){
    Cell_ind=match(Common_cell_dataset[[q]],Cell[[q]])
    final_dataset[[q]]=finalX[[q]][Cell_ind, ]
    final_dataset[[q]]=matrix(as.numeric(final_dataset[[q]]),nrow = dim(final_dataset[[q]])[1], ncol = dim(final_dataset[[q]])[2])
    
    Cell_ind_Y=match(Common_cell_dataset[[q]],finalY_train_cell)
    finalY_dataset[[q]]=matrix(finalY_train[Cell_ind_Y,],ncol=Variable_number)
    if (class(min_leaf)=="character" || min_leaf%%1!=0 || min_leaf<1 || min_leaf>nrow(finalY_dataset[[q]])) stop('Minimum leaf number can not be fractional or negative integer or string or greater than number of samples')  
  }
  if (class(n_tree)=="character" || n_tree%%1!=0 || n_tree<1) stop('Number of trees in the forest can not be fractional or negative integer or string')
  if (class(m_feature)=="character" || m_feature%%1!=0 || m_feature<1) stop('Number of randomly selected features considered for a split can not be fractional or negative integer or string')
  
  Y_hat_test=NULL#matrix(rep(0,length(Cell)*length(finalY_test_cell)),nrow=length(finalY_test_cell))
  final_test=NULL
  Common_cell_test=NULL
  for (q in 1:length(Cell)){
    Common_cell_test=intersect(finalY_test_cell,Cell[[q]])
    Cell_ind=match(Common_cell_test,Cell[[q]])
    final_test[[q]]=finalX[[q]][Cell_ind, ]
    final_test[[q]]=matrix(as.numeric(final_test[[q]]),nrow = dim(final_test[[q]])[1], ncol = dim(final_test[[q]])[2])
    
    Y_hat_test[[q]] = build_forest_predict(final_dataset[[q]],finalY_dataset[[q]], n_tree, m_feature, min_leaf, final_test[[q]])
  }
  
  Coeff2=rep( list(NULL), length(Serial) )   
  for (q in 1:length(Serial)){
    Temp2=Coeff[Serial[q][[1]]]
    Coeff2[[q]]=Temp2/sum(Temp2)
  }
  
  Common_cell=NULL
  Final_test=NULL
  Drug_sensitivity_cell_test2=finalY_test_cell
  for (q in 1:length(Serial)){
    D=NULL
    D=Drug_sensitivity_cell_test2
    for (check in 1:length(Serial[[q]])){
      D=intersect(D,Cell[Serial[[q]][check]][[1]])
    }
    Common_cell[[q]]=D
    Drug_sensitivity_cell_test2=setdiff(Drug_sensitivity_cell_test2,D)
    
    F_test=matrix(rep(0,length(Common_cell[[q]])*Variable_number),ncol=Variable_number)
    match_ind=NULL
    W=Serial[[q]]
    for (RR in 1:Variable_number){
      for (R in 1:length(W)){
        match_ind[[R]]=match(Common_cell[[q]],Common_cell_test[[W[R]]])
        F_test[,RR]=F_test[,RR]+Coeff2[[q]][R]*matrix(Y_hat_test[[W[R]]][match_ind[[R]],RR],ncol=1)
      }
    }
    
    Final_test[[q]]=F_test
  }
  Final_result=NULL
  for (q in 1:length(Serial)){
    if (length(Common_cell[[q]])>0){
      Final_result=rbind(Final_result,cbind(matrix(Common_cell[[q]],ncol=1),matrix(Final_test[[q]],ncol=Variable_number)))
    }
  }
  match_ind=match(finalY_test_cell,Final_result[,1])
  Final_prediction=Final_result[match_ind,]
  return(Final_prediction)
}