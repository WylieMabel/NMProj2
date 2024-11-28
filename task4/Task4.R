setwd("C:\\Users\\spisa\\Desktop\\Assignment2")

library(RSiena)
library(sna)
library(parallel)
library(igraph)
library(RColorBrewer)
library(ggraph)

plotNet = function(net, attributes, alcohol, distance, index)
{
  n = network(net, directed = TRUE)
  n %v% "gender" = attributes$gender
  n %v% "alcohol" = alcohol[,index]
  n %e% "distance" = distance
  
  layout = create_layout(n, layout = "auto")
  ggraph(layout) + 
    geom_edge_link0(
      arrow = arrow(
        angle = 10,
        length = unit(2, "mm"),
        type = "closed"),
      aes(edge_alpha = distance)
    ) + 
    scale_edge_alpha(range = c(1, 0.1), name = "Distance") + 
    geom_node_point(
      size = 3, 
      aes(
        fill = as.factor(alcohol),
        shape = as.factor(gender))
    ) +
    scale_fill_discrete(
      name = "Alcohol consumption",
      labels = c(
        "1" = "none",
        "2" = "once or twice a year",
        "3" = "once a month",
        "4" = "once a week",
        "5" = "more than once a week")
    ) +
    guides(fill = guide_legend(
      override.aes = list(shape = 21),
      labels = c(
        "1" = "none",
        "2" = "once or twice a year",
        "3" = "once a month",
        "4" = "once a week",
        "5" = "more than once a week")
    )) +
    scale_shape_manual (
      name = "Gender",
      values = c("1" = 22, "2" = 21),
      labels = c("1" = "boy", "2" = "girl")
    )
}

#' Adapted from section 6.3.3 RSiena manual
#' Run a \code{siena07} algorithm until convergence and export results
#' @param max_iter maximum iteration steps before giving up
#' @param name name of the model, used to make backups
#' @param prevAns previous model to use as a start
fitModel = function(alg, data, effects, max_iter, name = alg$projname, prevAns = NULL, ...)
{
  ############FILE MANAGMENT
  name = paste(name,".bin",sep="")
  if(!dir.exists("task4\\fits")) 
    dir.create(path = "task4\\fits", recursive = TRUE)
  if(!dir.exists("task4\\results")) 
    dir.create(path = "task4\\results", recursive = TRUE)
  results_file = paste(alg$projname,".txt",sep="")
  results_path = paste("task4\\results\\",results_file,sep="")
  ##########################
  
  iter = 0
  print("Begin fitting.")
  out = siena07(alg, data = data, effects = effects, prevAns = prevAns, ...)
  repeat 
  {
    save(out, file = file.path("task4\\fits",name))
    if(file.exists(results_file))
    {
      file.copy(from = results_file, to = results_path, overwrite = TRUE)
      file.remove(results_file)
    }
    iter = iter + 1
    conv = out$tconv.max
    print(paste("Iteration: ", iter, ". Convergence: ", conv, sep = ""))
    if ((conv < 0.25) && (sum(abs(out$tconv)<0.1) == length(out$tconv))) 
    {
      print("Convergence succeded.")
      return(out)
    }
    if (conv > 10) 
    {
      print("Covergence failed.")
      return(NULL)
    }
    if (iter > max_iter) 
    {
      print("Timed out.")
      return(NULL)
    }
    out = siena07(alg, data = data, effects = effects, prevAns = out, ...)
  }
}

gof = function(alg, model)
{
  name = alg$projname
  plot_path = paste("task4\\plots\\gof\\",name,sep="")
  if(!dir.exists(plot_path)) 
    dir.create(path = plot_path, recursive = TRUE)
  # use parallel computation
  cl <- makeCluster(4) # with 4 workers
  # Indegree distribution
  gofId <- sienaGOF(
    model,
    verbose = FALSE,
    varName = "net", IndegreeDistribution,
    cluster = cl
  )
  
  # Outdegree distribution
  gofOd <- sienaGOF(
    model,
    verbose = FALSE,
    varName = "net", OutdegreeDistribution,
    cluster = cl
  )
  
  # Triad census
  gofTC <- sienaGOF(
    model,
    verbose = FALSE,
    varName = "net", TriadCensus,
    cluster = cl
  )
  
  # Geodesic distance
  GeodesicDistribution <- function(
    i, data, sims, period, groupName,
    varName, levls = c(1:5, Inf), cumulative = TRUE) {
    x <- networkExtraction(i, data, sims, period, groupName, varName)
    require(sna)
    a <- sna::geodist(symmetrize(x))$gdist
    if (cumulative) {
      gdi <- sapply(levls, function(i) {
        sum(a <= i)
      })
    }
    else {
      gdi <- sapply(levls, function(i) {
        sum(a == i)
      })
    }
    names(gdi) <- as.character(levls)
    gdi
  }
  gofGD <- sienaGOF(
    model,
    verbose = FALSE,
    varName = "net", GeodesicDistribution
  )
  gofBeh <- sienaGOF(
    model,
    verbose = FALSE,
    varName = "alcohol", BehaviorDistribution,
    cluster = cl
  )
  stopCluster(cl)
  
  directories = list(
    paste(plot_path,"indegree_distribution.png",sep="\\"),
    paste(plot_path,"outdegree_distribution.png",sep="\\"),
    paste(plot_path,"triad_census.png",sep="\\"),
    paste(plot_path,"geodist_distribution.png",sep="\\"),
    paste(plot_path,"alcohol_behaviour.png",sep="\\")
  )
  settings = list(FALSE,FALSE,TRUE,FALSE,FALSE)
  gofs = list(gofId,gofOd,gofTC,gofGD,gofBeh)
  
  return(list(directories,settings,gofs))
}

if(!dir.exists("task4")) dir.create("task4")

net1 <- as.matrix(read.csv("Glasgow/f1.csv", header = FALSE))
net2 <- as.matrix(read.csv("Glasgow/f2.csv", header = FALSE))
net3 <- as.matrix(read.csv("Glasgow/f3.csv", header = FALSE))

attributes = as.data.frame(read.csv("Glasgow/demographic.csv", header = TRUE))
alcohol_matrix = as.matrix(read.csv("Glasgow/alcohol.csv", header = TRUE))
distance_matrix = as.matrix(read.csv("Glasgow/logdistance.csv", header = FALSE))

# plotNet(net1,attributes,alcohol_matrix,distance_matrix,1)
# plotNet(net2,attributes,alcohol_matrix,distance_matrix,2)
# plotNet(net3,attributes,alcohol_matrix,distance_matrix,3)

net = sienaDependent(array(c(net1,net2,net3), dim = c(129,129,3)))
age = coCovar(attributes$age)
gender = coCovar(attributes$gender)
distance = coDyadCovar(distance_matrix)
alcohol = sienaDependent(array(alcohol_matrix, dim = c(129,1,3)))

data = sienaDataCreate(net, gender, age, distance, alcohol)
#' Task 4.2.1: Jacard index is in the following report, in the last section. But 
#' here is a copy.
#' Period 1-->2 0.304
#' Period 2-->3 0.351
print01Report(data, modelname = "task4_data")
file.copy(from = "task4_data.txt", to = "task4\\task4_data.txt", overwrite = TRUE)
file.remove("task4_data.txt")

# Already includes density and reciprocity effects.
eff = getEffects(data)
{
  eff = includeEffects(eff, inPopSqrt)
  eff = includeEffects(eff, X, interaction1 = "distance")
  eff = includeEffects(eff, indeg, name = "alcohol", interaction1 = "net")
}

#' Using default values, see section 6.2 of Siena Manual. Increasing nsub and
#' n3 at the same time might improve convergence.
alg0 = sienaAlgorithmCreate(projname = "model0", seed = 37, nsub = 4, n3 = 1000)

# load("task4\\fits\\model0.bin")
# model0 = out # load restores the object with its original name: out (see fitModel), this changes it.

model0 = fitModel(alg0, data = data, effects = eff, max_iter = 30,
                  useCluster = TRUE, nbrNodes = 4,   # Parallel computing args
                  returnDeps = TRUE)  # Return simulated networks for gof

gof0 = gof(alg0,model0)
png(filename = gof0[[1]][[1]])
plot(gof0[[3]][[1]], center = gof0[[2]][[1]], scale = gof0[[2]][[1]])
dev.off()
png(filename = gof0[[1]][[2]])
plot(gof0[[3]][[2]], center = gof0[[2]][[2]], scale = gof0[[2]][[2]])
dev.off()
png(filename = gof0[[1]][[3]])
plot(gof0[[3]][[3]], center = gof0[[2]][[3]], scale = gof0[[2]][[3]])
dev.off()
png(filename = gof0[[1]][[4]])
plot(gof0[[3]][[4]], center = gof0[[2]][[4]], scale = gof0[[2]][[4]])
dev.off()
png(filename = gof0[[1]][[5]])
plot(gof0[[3]][[5]], center = gof0[[2]][[5]], scale = gof0[[2]][[5]])
dev.off()


{
  eff = includeEffects(eff, egoX, altX, sameX, interaction1 = "gender")
  eff = includeEffects(eff, avSim, name = "alcohol", interaction1 = "net")
}

alg1 = sienaAlgorithmCreate(projname = "model1", seed = 37, nsub = 4, n3 = 1000)
model1 = fitModel(alg1, data = data, effects = eff, max_iter = 30, prevAns = model0,
                  useCluster = TRUE, nbrNodes = 4, returnDeps = TRUE)
gof1 = gof(alg1,model1)
png(filename = gof1[[1]][[1]])
plot(gof1[[3]][[1]], center = gof1[[2]][[1]], scale = gof1[[2]][[1]])
dev.off()
png(filename = gof1[[1]][[2]])
plot(gof1[[3]][[2]], center = gof1[[2]][[2]], scale = gof1[[2]][[2]])
dev.off()
png(filename = gof1[[1]][[3]])
plot(gof1[[3]][[3]], center = gof1[[2]][[3]], scale = gof1[[2]][[3]])
dev.off()
png(filename = gof1[[1]][[4]])
plot(gof1[[3]][[4]], center = gof1[[2]][[4]], scale = gof1[[2]][[4]])
dev.off()
png(filename = gof1[[1]][[5]])
plot(gof1[[3]][[5]], center = gof1[[2]][[5]], scale = gof1[[2]][[5]])
dev.off()

{
  eff = includeEffects(eff, transTrip, cycle3)
  eff = includeEffects(eff,outActSqrt,outPopSqrt)
}

alg2 = sienaAlgorithmCreate(projname = "model2", seed = 37, nsub = 4, n3 = 1000)
model2 = fitModel(alg2, data = data, effects = eff, max_iter = 30, prevAns = model1,
                  useCluster = TRUE, nbrNodes = 4, returnDeps = TRUE)
gof2 = gof(alg2,model2)
png(filename = gof2[[1]][[1]])
plot(gof2[[3]][[1]], center = gof2[[2]][[1]], scale = gof2[[2]][[1]])
dev.off()
png(filename = gof2[[1]][[2]])
plot(gof2[[3]][[2]], center = gof2[[2]][[2]], scale = gof2[[2]][[2]])
dev.off()
png(filename = gof2[[1]][[3]])
plot(gof2[[3]][[3]], center = gof2[[2]][[3]], scale = gof2[[2]][[3]])
dev.off()
png(filename = gof2[[1]][[4]])
plot(gof2[[3]][[4]], center = gof2[[2]][[4]], scale = gof2[[2]][[4]])
dev.off()
png(filename = gof2[[1]][[5]])
plot(gof2[[3]][[5]], center = gof2[[2]][[5]], scale = gof2[[2]][[5]])
dev.off()

{
  eff = includeEffects(eff, outIso)
  eff = includeEffects(eff, name = "net", altX, egoX, simX, interaction1 = "alcohol")
}

alg3 = sienaAlgorithmCreate(projname = "model3", seed = 37, nsub = 4, n3 = 1000)
model3 = fitModel(alg3, data = data, effects = eff, max_iter = 30, prevAns = model2,
                  useCluster = TRUE, nbrNodes = 4, returnDeps = TRUE)
gof3 = gof(alg3,model3)
png(filename = gof3[[1]][[1]])
plot(gof3[[3]][[1]], center = gof3[[2]][[1]], scale = gof3[[2]][[1]])
dev.off()
png(filename = gof3[[1]][[2]])
plot(gof3[[3]][[2]], center = gof3[[2]][[2]], scale = gof3[[2]][[2]])
dev.off()
png(filename = gof3[[1]][[3]])
plot(gof3[[3]][[3]], center = gof3[[2]][[3]], scale = gof3[[2]][[3]])
dev.off()
png(filename = gof3[[1]][[4]])
plot(gof3[[3]][[4]], center = gof3[[2]][[4]], scale = gof3[[2]][[4]])
dev.off()
png(filename = gof3[[1]][[5]])
plot(gof3[[3]][[5]], center = gof3[[2]][[5]], scale = gof3[[2]][[5]])
dev.off()

siena.table(model3, type = "html", sig = TRUE)
  





