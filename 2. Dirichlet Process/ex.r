  1 #
  2 # Dirichlet Process Mixture Model example.
  3 # Last modified on Apr 19, 2018 by Yuelin Li from an example given by
  4 # Tamara Broderick, whose original source code is available at:
  5 # https://raw.githubusercontent.com/tbroderick/bnp_tutorial
  6 #      /2016mlss/ex7_sampler.R
  7 ###### Below are headers in the original source code file ######
  8 # A Gibbs sampler for a CRP Gaussian mixture model
  9 # Algorithm 3 in Neal 2000
 10 # Copyright (C) 2015, Tamara Broderick
 11 # www.tamarabroderick.com
 12 #
 13 # This program is free software: you can redistribute it and/or modify
 14 # it under the terms of the GNU General Public License as published by
 15 # the Free Software Foundation, either version 3 of the License, or
 16 # (at your option) any later version.
 17 #
 18 # This program is distributed in the hope that it will be useful,
 19 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 20 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 21 # GNU General Public License for more details.
 22 #
 23 # You should have received a copy of the GNU General Public License
 24 # along with this program. If not, see <-ttp://www.gnu.org/licenses/>.
 25 # useful for working with multivariate normal distributions
 26 #
 27 ##### Below are modified by YL. Variable names are changed from
 28 # the original to match the best we can the notation in the tutorial.
 29 crp_gibbs <- function(data,alpha=0.01,mu0,sigma0,sigma_y,c_init,maxIters=1000)
 30 {
 31 # data: an N by D matrix of data points
 32 # A small alpha encourages fewer clusters
 33 #    alpha <- 0.01
 34 # sigma_y: measurement error of data y, assumed known, same for all clusters.
 35 #    sigma_y <- diag(data_dim) * 1
 36 # mu0, sigma0: prior mean and variance around unknown mean mu
 37 #    mu0 <- matrix(rep(0, data_dim), ncol = data_dim, byrow = TRUE)
 38 #    sigma0 <- diag(data_dim) * 3^2
 39 # c_init: initial assignments of data points to clusters
 40 #
 41 require(mvtnorm)
 42 # dimension of the data points
 43 data_dim <- ncol(data)
 44 # number of data points
 45 N <- nrow(data)
 46 #####
 47 # Priors
 48 #####
 49 # prior precision on the unknown mean mu.
 50 tau0 <- solve(sigma0) # prior precision on $mu$, inverse of prior covariance
 51 # cluster-specific precision, assumed known, all clusters are assumed to
 52 # share identical measurement error of y ~ N(mu, sigma_y).
 53 tau_y <- solve(sigma_y)
 54 # initialize the CRP Gibbs sampler
 55 z <- c_init # initial cluster membership assignments
 56 n_k <- as.vector(table(z)) # initial data counts at each cluster, Eq (4).
 57 Nclust <- length(n_k) # initial number of clusters
 58 ##
 59 # Chinese Restaurant Process (CRP) Gibbs sampler begins
 60 ##
 61 res <- matrix(NA, nrow = N, ncol = maxIters) # cluster membership storage
 62 pb <- txtProgressBar(min = 0, max = maxIters, style = 3)
 63 for(iter in 1:maxIters) { # maxIters also prevents endless loops
 64 for(n in 1:N) { # one data point (customer) at a time
 65 # when nth customer enters the Chinese restaurant, we need to first
 66 # un-assign his/her initial cluster membership, then use Eq (12)
 67 # [already occupied table] and (13) [new table] to calculate the
 68 # updated probability of table assignment.
 69 c_i <- z[n] # what is the nth persons table assignment?
 70 n_k[c_i] <- n_k[c_i] - 1 # remove the nth person from table
 71 # if the table becomes empty when the nth person is removed,
 72 # then that table/cluster is removed.
 73 if( n_k[c_i] == 0 )
 74 {
 75   n_k[c_i] <- n_k[Nclust]  # last cluster to replace this empty cluster
 76   loc_z <- ( z == Nclust ) # who are in the last cluster?
 77   z[loc_z] <- c_i          # move them up to fill just emptied cluster
 78   n_k <- n_k[ -Nclust ]    # take out the last cluster, now empty
 79   Nclust <- Nclust - 1     # decrease total number of clusters by 1
 80 }
 81 z[n] <- −1 # ensures z[n] doesnt get counted as a cluster #####
 82 ##
 83 # Now we are ready to update table assignment by Eqs (12) and (13).
 84 ##
 85 # log probabilities for the clusters, add previously unoccupied table
 86 logp <- rep( NA, Nclust + 1 )
 87 # loop over already occupied tables 1:J and calculate pr as per Eq (13).
 88 for( c_i in 1:Nclust ) {
 89   tau_p <- tau0 + n_k[c_i] * tau_y # cluster precision as per Eq (4)
 90   sig_p <- solve(tau_p) # cluster variance, inverse of tau_c
 91   # find all of the points in this cluster
 92   loc_z <- which(z == c_i)
 93   # sum all the points in this cluster
 94   if(length(loc_z) > 1) {
 95     sum_data <- colSums(data[z == c_i, ]) }
 96   else {
 97     sum_data <- data[z == c_i, ]
 98   }
 99 #
100 # We need to use the predictive distribution of each already
101 # occupied table to predict the next customer sitting there.
102 #
103 # Use Eq (4) to find the conditional posterior distribution for
104 # the cluster means, (y * n_k * tau_j + mu0 * s0) /tau_p,and
105 # then use the predictive distribution of y_j in Eq (11) to
106 # predict new data value c_i from c-i.
107 #
108 mean_p <- sig_p %*% (tau_y %*% sum_data + tau0 %*% t(mu0))
109 logp[c_i] <- log(n_k[c_i]) +
110 dmvnorm(data[n,], mean = mean_p, sigma = sig_p + sigma_y, log = TRUE) }
111 #
112 # We are done looping over already occupied tables. Next, we use
113 # Eq (12) to calcualte the log probability of a previously
114 # unoccupied, 'new' table. Essentially,it is the prior predicitive
115 # distribution of the DP.
116 #
117 logp[ Nclust+1 ] <- log(alpha) +
118     dmvnorm(data[n,], mean = mu0, sigma = sigma0 + sigma_y, log = TRUE)
119 # transform unnormalized log probabilities into probabilities
120 max_logp <- max(logp)
121 logp <- logp - max_logp
122 loc_probs <- exp(logp)
123 loc_probs <- loc_probs / sum(loc_probs)
124 # draw a sample of which cluster this new customer should belong to
125 newz <- sample(1:(Nclust+1), 1, replace = TRUE, prob = loc_probs)
126 # spawn a new cluster if necessary
127 if(newz == Nclust + 1) {
128   n_k <- c(n_k, 0)
129   Nclust <- Nclust + 1
130   }
131 z[n] <- newz
132 n_k[newz] <- n_k[newz] + 1 # update the cluster n_k
133 }
134 setTxtProgressBar(pb, iter) # update text progress bar after each iter
135 res[, iter] <- z     # cluster membership of N observations
136 }
137 close(pb)            # close text progress bar
138 invisible(res)       # return results, N by maxIters matrix
139 }
