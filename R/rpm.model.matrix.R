#' Creates a model matrix to estimate the parameters of a Revealed Preference Matchings Model
#' 
#' \code{\link{rpm.model.matrix}} assumes a bipartite network (i.e. two-sided matching market)
#' It creates a model matrix according to the formula passed in.
#' 
#' @param model.terms For the details on the possible \code{<model terms>}, see
#' \code{\link{rpm-terms}}.
#' @param model.terms.coef.names the covariates used to construct the model matrix.
#' They are used in conjunction with the model terms. 
#' @param Xall the unique types of b1
#' @param Zall the unique types of b2
#' @return A list consists of the following elements:
#' \item{X}{the model matrix for b1.}
#' \item{Z}{the model matrix for b2.}
#' \item{Xnames}{the names of the covariates for b1.} 
#' \item{Znames}{the names of the covariates for b2.}
#' @seealso fitrpm
#' @references Menzel, Konrad (2015).
#' \emph{Large Matching Markets as Two-Sided Demand Systems}
#' Econometrica, Vol. 83, No. 3 (May, 2015), 897-941.
#' @keywords models
#' @examples
#' 
#' #
#' # load the Florentine marriage data matrix
#' #
#' data(flo)
#' #
#' # attach the sociomatrix for the Florentine marriage data
#' # This is not yet a network object.
#' #
#' flo
#' #
#' # Create a network object out of the adjacency matrix
#' #
#' flomarriage <- network(flo,directed=FALSE)
#' flomarriage
#' #
#' # print out the sociomatrix for the Florentine marriage data
#' #
#' flomarriage[,]
#' #
#' # create a vector indicating the wealth of each family (in thousands of lira) 
#' # and add it as a covariate to the network object
#' #
#' flomarriage %v% "wealth" <- c(10,36,27,146,55,44,20,8,42,103,48,49,10,48,32,3)
#' flomarriage
#' #
#' # create a plot of the social network
#' #
#' plot(flomarriage)
#' #
#' # now make the vertex size proportional to their wealth
#' #
#' plot(flomarriage, vertex.cex="wealth", main="Marriage Ties")
#' #
#' # Use 'data(package = "ergm")' to list the data sets in a
#' #
#' data(package="ergm")
#' #
#' # Load a network object of the Florentine data
#' #
#' data(florentine)
#' #
#' # Fit a model where the propensity to form ties between
#' # families depends on the absolute difference in wealth
#' #
#' gest <- ergm(flomarriage ~ edges + absdiff("wealth"))
#' summary(gest)
#' #
#' # add terms for the propensity to form 2-stars and triangles
#' # of families 
#' #
#' gest <- ergm(flomarriage ~ kstar(1:2) + absdiff("wealth") + triangle)
#' summary(gest)
#' 
#' # import synthetic network that looks like a molecule
#' data(molecule)
#' # Add a attribute to it to mimic the atomic type
#' molecule %v% "atomic type" <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3)
#' #
#' # create a plot of the social network
#' # colored by atomic type
#' #
#' plot(molecule, vertex.col="atomic type",vertex.cex=3)
#' 
#' # measure tendency to match within each atomic type
#' gest <- ergm(molecule ~ edges + kstar(2) + triangle + nodematch("atomic type"),
#'   MCMCsamplesize=10000)
#' summary(gest)
#' 
#' # compare it to differential homophily by atomic type
#' gest <- ergm(molecule ~ edges + kstar(2) + triangle
#'                               + nodematch("atomic type",diff=TRUE),
#'   MCMCsamplesize=10000)
#' summary(gest)
#'
#' @export rpm.model.matrix
#' 
rpm.model.matrix = function(model.terms.names, model.terms.coef.names, Xall, Zall)		
{
    # assumes Xall and Zall are the unique types
    
    # no. of members on each side
    nb1 = nrow(Xall) 
    nb2 = nrow(Zall)

    ncov = length(model.terms.names)
    
    # returned variables
    X = matrix(1, nrow=nb1, ncol=nb2)
    Z = matrix(1, nrow=nb2, ncol=nb1)
    Xnames = "b1.intercept"
    Znames = "b2.intercept"
    
    for(i in 1:ncov)
    {
        switch(model.terms.names[i],
               b1cov = {
                   for(j in 1 : length(model.terms.coef.names[[i]]))
                   {	
                       attrname = model.terms.coef.names[[i]][j]
                       
                       # every row has the same value -- group 1's raw values for this attribute
                       X = abind(X,
                                 matrix(rep(Xall[, attrname], times = nb2), ncol = nb2, byrow = F),
                                 along = 3)
                       Xnames = c(Xnames, paste("b1cov", attrname, sep="."))
                   }
               },
               b2cov = {
                   for(j in 1 : length(model.terms.coef.names[[i]]))
                   {	
                       attrname = model.terms.coef.names[[i]][j]
                       
                       # every row has the same value -- group 2's raw values for this attribute
                       Z = abind(Z,
                                 matrix(rep(Zall[, attrname], times = nb1), ncol = nb1, byrow = F),
                                 along = 3)
                       Znames = c(Znames, paste("b2cov", attrname, sep="."))
                   }
               },
               b1factor = {
                   for(j in 1 : length(model.terms.coef.names[[i]]))
                   {	
                       attrname = model.terms.coef.names[[i]][j]
                       
                       # get all the factor levels possible for both sides
                       temp = sort(unique(c(Xall[,attrname], Zall[,attrname])))
                       
                       for(k in temp)
                       {
                           # every row has the same value --  1 if group 1's value for this attribute is same as k; 0 otherwise
                           X = abind(X, 
                                     matrix(rep((1 * (Xall[, attrname] == k)), times = nb2), ncol = nb2, byrow = F),
                                     along = 3)
                           Xnames = c(Xnames, paste("b1factor", attrname, k, sep="."))
                       }
                   }
               },
               b2factor = {
                   for(j in 1 : length(model.terms.coef.names[[i]]))
                   {			
                       attrname = model.terms.coef.names[[i]][j]
                       
                       # get all the factor levels possible for both sides
                       temp = sort(unique(c(Xall[,attrname], Zall[,attrname])))
                       
                       for(k in temp)
                       {
                           # every row has the same value --  1 if group 2's value for this attribute is same as k; 0 otherwise
                           Z = abind(Z, 
                                     matrix(rep((1 * (Zall[, attrname] == k)), times = nb1), ncol = nb1, byrow = F),
                                     along = 3)
                           Znames = c(Znames, paste("b2factor", attrname, k, sep="."))
                       }
                   }
               },
               b1nodematch = {
                   parts = model.terms.coef.names[[i]]
                   
                   # only the matches for the main effect 
                   if(length(parts) == 1)
                   {
                       attrname = parts
                       
                       # get all the factor levels possible for both sides
                       u = sort(unique(c(Xall[,attrname], Zall[,attrname])))
                       
                       # matching on each level
                       for (k in u)
                       {
                           # only matching on kth level
                           b1 = 1 * (Xall[,attrname] == k)
                           b2 = 1 * (Zall[,attrname] == k)
                           
                           X = abind(X,
                                      1*outer(b1, b2, "*"),
                                      along = 3) 
                           Xnames = c(Xnames, paste("b1nodematch", attrname, k, sep="."))
                       }
                   }
                   else # the interaction between matching attributes
                   {
                       unique.combo = unique(rbind(Xall[,parts], Zall[,parts]))
                       
                       # matching on each combination of level
                       for (j in 1:nrow(unique.combo))
                       {
                           b1 = 1 * apply(Xall[,parts], 1, function(x) identical(x, unique.combo[j,]))
                           b2 = 1 * apply(Zall[,parts], 1, function(x) identical(x, unique.combo[j,]))
                           
                           X = abind(X,
                                      1*outer(b1, b2, "*"),
                                      along = 3) 
                           Xnames = c(Xnames, paste("b1nodematch", paste(parts, unique.combo[j,], sep=".", collapse = "."), sep="."))
                       }
                   }
               },
               b2nodematch = {
                   parts = model.terms.coef.names[[i]]
                   
                   # only the matches for the main effect 
                   if(length(parts) == 1)
                   {
                       attrname = parts
                       
                       # get all the factor levels possible for both sides
                       u = sort(unique(c(Xall[,attrname], Zall[,attrname])))
                       
                       # matching on each level
                       for (k in u)
                       {
                           # only matching on kth level
                           b1 = 1 * (Xall[,attrname] == k)
                           b2 = 1 * (Zall[,attrname] == k)
                           
                           Z = abind(Z,
                                      1*outer(b2, b1, "*"),
                                      along = 3) 
                           Znames = c(Znames, paste("b2nodematch", attrname, k, sep="."))
                       }
                   }
                   else # the interaction between matching attributes
                   {
                       unique.combo = unique(rbind(Xall[,parts], Zall[,parts]))
                       
                       # matching on each combination of level
                       for (j in 1:nrow(unique.combo))
                       {
                           b1 = 1 * apply(Xall[,parts], 1, function(x) identical(x, unique.combo[j,]))
                           b2 = 1 * apply(Zall[,parts], 1, function(x) identical(x, unique.combo[j,]))
                           
                           Z = abind(Z,
                                      1*outer(b2, b1, "*"),
                                      along = 3) 
                           Znames = c(Znames, paste("b2nodematch", paste(parts, unique.combo[j,], sep=".", collapse = "."), sep="."))
                       }
                   }
               },
               b1homophily = {
                   parts = model.terms.coef.names[[i]]
                   
                   # only the matches for the main effect 
                   if(length(parts) == 1)
                   {
                       attrname = parts
                    
                       # matching values on attribute
                       X = abind(X,
                                  1*outer(Xall[, attrname], Zall[, attrname], "=="),
                                  along = 3) 
                       Xnames = c(Xnames, paste("b1homophily", attrname, sep="."))
                   }
                   else # matching on more than a single attribute
                   {
                       unique.combo = unique(rbind(Xall[,parts], Zall[,parts]))
                       
                       # assign unique group membership for group 1 
                       Xgroupmem = rep(NA,nrow(Xall))
                       for(i in 1:nrow(unique.combo)){
                           Xgroupmem[apply(Xall[,parts], 1, function(x) identical(x, unique.combo[i,]))] = i
                       }
                       
                       # assign unique group membership for group 2 
                       Zgroupmem = rep(NA,nrow(Zall))
                       for(i in 1:nrow(unique.combo)){
                           Zgroupmem[apply(Zall[,parts], 1, function(x) identical(x, unique.combo[i,]))] = i
                       }
                       
                       # matching values on attribute
                       X = abind(X,
                                  1*outer(Xgroupmem, Zgroupmem, "=="),
                                  along = 3) 
                       Xnames = c(Xnames, paste("b1homophily", paste(parts, sep=".", collapse = "."), sep="."))
                   }
               },
               b2homophily = {
                   parts = model.terms.coef.names[[i]]
                   
                   # only the matches for the main effect 
                   if(length(parts) == 1)
                   {
                       attrname = parts
                       
                       # matching values on attribute
                       Z = abind(Z,
                                  1*outer(Zall[, attrname], Xall[, attrname], "=="),
                                  along = 3) 
                       Znames = c(Znames, paste("b2homophily", attrname, sep="."))
                       
                   }
                   else # matching on more than a single attribute
                   {
                       unique.combo = unique(rbind(Xall[,parts], Zall[,parts]))
                       
                       # assign unique group membership for group 1 
                       Xgroupmem = rep(NA,nrow(Xall))
                       for(i in 1:nrow(unique.combo)){
                           Xgroupmem[apply(Xall[,parts], 1, function(x) identical(x, unique.combo[i,]))] = i
                       }
                       
                       # assign unique group membership for group 2 
                       Zgroupmem = rep(NA,nrow(Zall))
                       for(i in 1:nrow(unique.combo)){
                           Zgroupmem[apply(Zall[,parts], 1, function(x) identical(x, unique.combo[i,]))] = i
                       }
                       
                       # matching values on attribute
                       Z = abind(Z,
                                  1*outer(Zgroupmem, Xgroupmem, "=="),
                                  along = 3) 
                       Znames = c(Znames, paste("b2homophily", paste(parts, sep=".", collapse = "."), sep="."))
                   }
               },
               b1smallerthan = {
                   for(j in 1 : length(model.terms.coef.names[[i]]))
                   {			
                       attrname = model.terms.coef.names[[i]][j]
                   
                       X = abind(X,
                                 1*outer(Xall[,attrname], Zall[,attrname],">"),  # a liking for group 2 of lower level
                                 along = 3)
                       Xnames = c(Xnames, paste("b1smallerthan", attrname, sep="."))
                   }
               },
               b1greaterthan = {
                   for(j in 1 : length(model.terms.coef.names[[i]]))
                   {			
                       attrname = model.terms.coef.names[[i]][j]
                       
                       X = abind(X,
                                 1*outer(Xall[,attrname], Zall[,attrname],"<"),  # a liking for group 2 of higher level
                                 along = 3)
                       Xnames = c(Xnames, paste("b1greaterthan", attrname, sep="."))
                   }
               },
               b2smallerthan = {
                   for(j in 1 : length(model.terms.coef.names[[i]]))
                   {			
                       attrname = model.terms.coef.names[[i]][j]
                       
                       Z = abind(Z,
                                 1*outer(Zall[,attrname], Xall[,attrname],">"),  # a liking for group 1 of lower level
                                 along = 3)
                       Znames = c(Znames, paste("b2smallerthan", attrname, sep="."))
                   }
               },
               b2greaterthan = {
                   for(j in 1 : length(model.terms.coef.names[[i]]))
                   {			
                       attrname = model.terms.coef.names[[i]][j]
                       
                       Z = abind(Z,
                                 1*outer(Zall[,attrname], Xall[,attrname],"<"),  # a liking for group 1 of higher level
                                 along = 3)
                       Znames = c(Znames, paste("b2greaterthan", attrname, sep="."))
                   }
               },
               b1absdiff = {
                   
                   attrname = model.terms.coef.names[[i]][1]
                   difflev = model.terms.coef.names[[i]][2]
                       
                   X = abind(X,
                             1*outer(Xall[,attrname], Zall[,attrname], function(x,z){1*(abs(x-z)==difflev)}),  
                             along = 3)
                   Xnames = c(Xnames, paste("b1absdiff", attrname, difflev, sep="."))
               },
               b2absdiff = {
                   
                   attrname = model.terms.coef.names[[i]][1]
                   difflev = model.terms.coef.names[[i]][2]
                   
                   Z = abind(Z,
                             1*outer(Zall[,attrname], Xall[,attrname], function(z,x){1*(abs(x-z)==difflev)}),  
                             along = 3)
                   Znames = c(Znames, paste("b2absdiff", attrname, difflev, sep="."))
               },
               # default
               {
                   print(paste0("Error: unrecognized ERGM term: ", model.terms.names[i]))
               }
               
               )
    }

    return(list(X=X, Z=Z, Xnames=Xnames, Znames=Znames))
}
