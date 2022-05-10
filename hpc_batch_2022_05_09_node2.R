
##
## Artemis batch file run: Date 09/05/2022
##

library(data.table)
library(pbapply)

## generate data
##
gen_data_arm <- function(n=c(90),
                     rand_assignment=c(1,1),
                     rand_pr=c(0.5,0.5)){
  ##
  ## n is the sample in one interim 
  ##
  dat <- data.table::CJ(
    id = 1:n,
    id_arm = NA_integer_
  )
  if(is.null(rand_assignment)){
    n_arm <- length(rand_assignment)
    rand_pr <- rand_assignment/sum(rand_assignment)
  }
  else{
    n_arm <- length(rand_pr)
    rand_pr <- rand_pr
  }
  dat$id_arm <- sample(1:n_arm, size=nrow(dat), replace=TRUE, prob=rand_pr)
  dat
}
#dat <- gen_data_arm(n=90,rand_assignment=c(1,1))
#table(dat$id_arm)
##
trial_betabino <- function(n,
                           rand_assignment=c(1,1,1),
                           rand_pr=c(0.3,0.3,0.3),
                           n_theta=c(0.75,0.75,0.75),
                           delta=0.15,
                           iter=1000,
                           prior_a=1,prior_b=1){
  ##
  y <- gen_data_arm(n=n,rand_assignment=rand_assignment,rand_pr=rand_pr)
  ##
  n_sample_arm <- table(y$id_arm)
  y[, evt := rbinom(n_sample_arm[id_arm], 1, n_theta[id_arm])]
  d_res <- merge(
    data.table(id_arm = 1:length(n_sample_arm)),
    y[, .(n = .N, y = sum(evt)), keyby = id_arm],
    by = "id_arm", all.x = T
  )
  post <- do.call(cbind, lapply(1:d_res[, .N],
                                function(k){
                                  rbeta(iter, 
                                        prior_a + d_res[id_arm ==k, y], 
                                        prior_b + d_res[id_arm ==k, n] - d_res[id_arm ==k, y])
                                }))
  decision_pr <- c()
  for(i in 1:(length(n_sample_arm)-1)){
    decision_pr[i] <- mean(post[,i+1]>post[,1]-delta)
  }
  d_res$decision_pr <- c(0,decision_pr)
  d_res$est_theta <- apply(post,2,mean)
  d_res$est_theta_sd <- apply(post,2,sd)
  d_res
}
##

## run model for different operating characteristics

n_theta <- cbind(seq(0.60,0.90,by=0.05), seq(0.60,0.90,by=0.05))
delta <- c(0.1,0.15,0.2)
res_mat <- NULL
ds_fnc <- function(x){
  c(mean(x>0.900),mean(x>0.950),mean(x>0.975),mean(x<0.025),mean(x<0.05),mean(x<0.10))
}
nSim <- 5000
id_scenario <- 0
for(i in 1:length(delta)){
  for(j in 1:nrow(n_theta)){
    for(jj in 1:nrow(n_theta)){
        out <- pblapply(1:nSim, function(x) trial_betabino(n=120,rand_assignment=c(1,1),rand_pr=c(0.5,0.5),
                                                           n_theta=c(n_theta[j,1],n_theta[jj,2]),
                                                           delta=delta[i]))
        est_theta <- apply(sapply(out, function(x) x$est_theta),1,mean)
        est_theta_sd <- apply(sapply(out, function(x) x$est_theta_sd),1,mean)
        decision_pr <- sapply(out, function(x) x$decision_pr)
        decision_pr <- t(apply(decision_pr,1,ds_fnc))
		id_scenario <- id_scenario + 1
        res_mat <- rbind(res_mat,cbind(1:2,id_scenario,60,delta[i],c(n_theta[j,1],n_theta[jj,2]),est_theta,est_theta_sd,decision_pr))
    }
  }
}
dimnames(res_mat)[[2]] <- c("arms","id_scenario","sample_size","delta","true_theta","est_theta","est_theta_sd",
                            "pr_noninferior_0900","pr_noninferior_0950","pr_noninferior_0975",
							"pr_inferior_025","pr_inferior_050","pr_inferior_010")
res_mat <- data.frame(res_mat)
##
write.csv(res_mat,file=paste0("result_60_arm2_prior1_",Sys.Date(),".csv"),row.names=FALSE)
##
