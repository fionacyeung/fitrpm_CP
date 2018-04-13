#' This is the version of Gale-Shapley stable matching algorithm (translated from Menzel's Matlab code?)
#' 
#' This code allows the self-matched option
#' 
#' @param U The utility matrix for the women's side. Each row is a woman, each column is a man.
#' The matrix entry (i,j) is the utility that woman \code{i} gains from pairing with man \code{j}. 
#' In other words, the utility is computed from woman \code{i}'s perspective.
#' @param V The utility matrix for the men's side. Each row is a man, each column is a woman.
#' The matrix entry (i,j) is the utility that man \code{i} gains from pairing with woman \code{j}. 
#' In other words, the utility is computed from man \code{i}'s perspective.
#' @return This function returns the following matrix: 
#' \item{mu}{The matching matrix, where 1 represents a pairing, 0 otherwise. 
#' Each row is a woman, each column is a man. The order of the rows is the same as the 
#' rows in \code{U}. The order of the columns is the same as the rows in \code{V}.}
#' @seealso fitrpm_R_CP
#' @references Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @keywords models
#' @examples
#' 
Gale_Shapley <- function(U,V){

# first argument proposing side
# proposing side: individuals correspond to rows

nw <- nrow(U)
nm <- ncol(U)
U_temp <- U
nmax <- 10*nw*nm

for (i in 1:nmax){

 # check dimensions!
 # no need to worry about ties
    Prop <- sweep(U_temp, 2, apply(rbind(U_temp,0),2,max),"==")
    Rej <- (Prop*V < sweep(Prop,1,apply(cbind(Prop*V,0),1,max),"*"))
    U_temp[Rej&Prop] = -1

    if (sum(sum(Rej))==0){
        break
    }
}

mu <- sweep(U_temp, 2, apply(rbind(U_temp,0),2,max),"==")

#mu <- cbind(1:nrow(mu),apply(mu,1,function(x){match(1,x,nomatch=0)}))
#mu[mu[,2]>0,]
1*mu
}
