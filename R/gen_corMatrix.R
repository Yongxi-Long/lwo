#' Generate correlation matrix
#'
#' @param n_visits Numerical value specifying the number of visits of the repeated measurements, which
#' determines the dimension of the correlation matrix
#' @param rho Correlation coefficient
#' @param corstr A character string specifying the correlation structure. Options include `"independence"|"ind"`, `"exchangeable"|"cs"|"compound symmetry"|"exch"`, `"ar1"`
#'
#' @return A matrix of dimension n_visit * n_visit
#'
#' @examples
#' gen_corMatrix(n_visits=3,rho=0.6,corstr="ar1")
#' @export
gen_corMatrix <- function(n_visits,rho,corstr)
{
  if(corstr %in% c("independence","ind"))
  {
    diag(n_visits)
  } else if (corstr == "ar1")
  {
    exponent <- abs(matrix(1:n_visits - 1, nrow = n_visits, ncol = n_visits, byrow = TRUE) -
                      (1:n_visits - 1))
    rho^exponent
  } else if (corstr %in% c("exchangeable","cs","compound symmetry","exch"))
  {
    mat <- matrix(rho, nrow = n_visits, ncol = n_visits)
    diag(mat) <- 1
    mat
  } else
  {
    stop("unknown correlation structure")
  }
}
