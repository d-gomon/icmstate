
scenario <- c(3)
n <- 500
n_obs <- 6
N <- 500
methods <- c("binomial", "msm", "poisson") #Change to have any of c("poisson", "multinomial", "msm")
RNG <- c("", "RNG2")


#Create file names to load:
load_names <- expand.grid(scenario, n, n_obs, N, methods, RNG)
#Because we have only ran n <- 500 with N <- 500 we need an extra line:
#which_remove <- which(load_names[, 1] == 4 & load_names[, 5] == "poisson")
#load_names <- load_names[-which_remove,]
colnames(load_names) <- c("scenario", "n", "n_obs", "N", "method", "RNG")
load_names[, 5] <- as.character(load_names[, 5])
load_names[, 6] <- as.character(load_names[, 6])
load_names <- load_names[order(load_names$scenario, load_names$method, load_names$RNG), ]
names_load <- c()
var_names <- c()
var_names_new <- c()
for(i in 1:nrow(load_names)){
  x <- load_names[i,]
  names_load <- c(names_load, paste0("sc", x[1], "n", x[2], "obs", trimws(x[3]), "N", x[4], x[5], x[6], ".Rdata"))
  var_names <- c(var_names, paste0("sc", x[1], "n", x[2], "obs", trimws(x[3]), "N", x[4], x[5], x[6]))
}

load_names2 <- load_names
load_names2[which(load_names[,5] == "binomial"), 5] <- "multinomial"
for(i in 1:nrow(load_names2)){
  x <- load_names2[i,]
  var_names_new <- c(var_names_new, paste0("sc", x[1], "n", x[2], "obs", trimws(x[3]), "N", 1000, x[5]))  
}



for(i in 1:(length(names_load)/2)){
  load(names_load[(i-1)*2 + 1])
  load(names_load[(i-1)*2 + 2])
  assign(var_names_new[(i-1)*2 + 1], c(get(var_names[(i-1)*2 + 1]), get(var_names[(i-1)*2 + 2])))
  save(list = var_names_new[(i-1)*2 + 1], file = paste0("newfiles/", var_names_new[(i-1)*2 + 1], ".Rdata"))
  rm(list = c(var_names_new[(i-1)*2+1], var_names[(i-1)*2+1], var_names[(i-1)*2+2]))
}






