####################################################################################
#                      Script for benchmark data sets                              #
#                                                                                  #
# In: Parsimonious Ultrametric Gaussian Mixture Models. Submitted for publication. #                                                      #
####################################################################################

# If necessary
# install.packages("MASS")
library(MASS)
# If necessary
# install.packages("mclust")
library(mclust)
# If necessary
# install.packages("pgmm")
library(pgmm)
# If necessary
# install.packages("HDclassif")
library(HDclassif)
# If necessary
# install.packages("IMIFA")
library(IMIFA)
# If necessary
# install.packages("remotes")
remotes::install_github("PUGMM-authors/PUGMM")
library(PUGMM)

################################ WINE 13 #######################################
data(wine, package = "HDclassif")
x <- scale(wine[, -1])
lab.wine <- wine[, 1]
G.max <- length(unique(lab.wine)) + 2
if (dim(x)[2] < 5){m.max <- dim(x)[2]} else {m.max <- 5}

pugmm.wine <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.wine <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.wine <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.wine <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.wine$G
pugmm.wine$m
pugmm.wine$model.name
pugmm.wine$pm.free
adjustedRandIndex(lab.wine, pugmm.wine$label)
# GPCMs
mclust.wine$G
mclust.wine$modelName
nMclustParams(mclust.wine$modelName, mclust.wine$d, mclust.wine$G)
adjustedRandIndex(lab.wine, mclust.wine$classification)
# PGMMs
pgmm.wine$g
pgmm.wine$q
pgmm.wine$model
PGMM_dfree(pgmm.wine$q, mclust.wine$d, pgmm.wine$g, pgmm.wine$model)
adjustedRandIndex(lab.wine, pgmm.wine$map)
# HDDC
hddc.wine$K
hddc.wine$d
hddc.wine$model
hddc.wine$complexity
adjustedRandIndex(lab.wine, hddc.wine$class)

################################ WINE 27 #######################################
data(wine, package = "pgmm")
x <- scale(wine[, -1])
lab.wine27 <- wine[, 1]
G.max <- length(unique(lab.wine27)) + 2
if (dim(x)[2] < 5){m.max <- dim(x)[2]} else {m.max <- 5}

pugmm.wine27 <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.wine27 <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.wine27 <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.wine27 <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.wine27$G
pugmm.wine27$m
pugmm.wine27$model.name
pugmm.wine27$pm.free
adjustedRandIndex(lab.wine27, pugmm.wine27$label)
# GPCMs
mclust.wine27$G
mclust.wine27$modelName
nMclustParams(mclust.wine27$modelName, mclust.wine27$d, mclust.wine27$G)
adjustedRandIndex(lab.wine27, mclust.wine27$classification)
# PGMMs
pgmm.wine27$g
pgmm.wine27$q
pgmm.wine27$model
PGMM_dfree(pgmm.wine27$q, mclust.wine27$d, pgmm.wine27$g, pgmm.wine27$model)
adjustedRandIndex(lab.wine27, pgmm.wine27$map)
# HDDC
hddc.wine27$K
hddc.wine27$d
hddc.wine27$model
hddc.wine27$complexity
ari.hddc <- adjustedRandIndex(lab.wine27, hddc.wine27$class)
table(lab.wine27, hddc.wine27$class)
hddc.wine27$loglik
hddc.wine27$BIC

# Gopt
# PGMMs
pgmm.wine27.Gopt <- pgmmEM(x, rG = length(unique(lab.wine27)), rq = 1:m.max, relax = TRUE)
pgmm.wine27.Gopt$q
pgmm.wine27.Gopt$model
PGMM_dfree(pgmm.wine27.Gopt$q, mclust.wine27$d, pgmm.wine27.Gopt$g, pgmm.wine27.Gopt$model)  #UUCU
adjustedRandIndex(lab.wine27, pgmm.wine27.Gopt$map)
# HDDC
hddc.wine27.Gopt <- hddc(x, K = length(unique(lab.wine27)), model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))
hddc.wine27.Gopt$model
hddc.wine27.Gopt$complexity
adjustedRandIndex(lab.wine27, hddc.wine27.Gopt$class)

################################ THYROID #######################################
data(thyroid, package = "mclust")
thyroid <- as.data.frame(thyroid)
lab.thyroid <- as.numeric(thyroid[, 1])
x <- scale(thyroid[, -1])
G.max <- length(unique(lab.thyroid)) + 2
if (dim(x)[2] < 5){m.max <- dim(x)[2]} else {m.max <- 5}

pugmm.thyroid <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.thyroid <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.thyroid <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.thyroid <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.thyroid$G
pugmm.thyroid$m
pugmm.thyroid$model.name
pugmm.thyroid$pm.free
adjustedRandIndex(lab.thyroid, pugmm.thyroid$label)
# GPCMs
mclust.thyroid$G
mclust.thyroid$modelName
nMclustParams(mclust.thyroid$modelName, mclust.thyroid$d, mclust.thyroid$G)
adjustedRandIndex(lab.thyroid, mclust.thyroid$classification)
# PGMMs
pgmm.thyroid$g
pgmm.thyroid$q
pgmm.thyroid$model
PGMM_dfree(pgmm.thyroid$q, mclust.thyroid$d, pgmm.thyroid$g, pgmm.thyroid$model)
adjustedRandIndex(lab.thyroid, pgmm.thyroid$map)
# HDDC
hddc.thyroid$K
hddc.thyroid$d
hddc.thyroid$model
hddc.thyroid$complexity
adjustedRandIndex(lab.thyroid, hddc.thyroid$class)

################################ SOBAR #########################################
x <- read.csv("DIRECTORY_NAME/sobar.csv", sep = " ")
lab.sobar <- x[, 20]
G.max <- length(unique(lab.sobar)) + 2
if (dim(x)[2] < 5){m.max <- dim(x)[2]} else {m.max <- 5}

pugmm.sobar <- pugmm(scale(x[, -20]), G = 1:G.max, m = 1:m.max)
mclust.sobar <- Mclust(scale(x[, -20]), G = 1:G.max, control = emControl(itmax = 500))
pgmm.sobar <- pgmmEM(scale(x[, -20]), rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.sobar <- hddc(scale(x[, -20]), K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.sobar$G
pugmm.sobar$m
pugmm.sobar$model.name
pugmm.sobar$pm.free
adjustedRandIndex(lab.sobar, pugmm.sobar$label)
# GPCMs
mclust.sobar$G
mclust.sobar$modelName
nMclustParams(mclust.sobar$modelName, mclust.sobar$d, mclust.sobar$G)
adjustedRandIndex(lab.sobar, mclust.sobar$classification)
# PGMMs
pgmm.sobar$g
pgmm.sobar$q
pgmm.sobar$model
PGMM_dfree(pgmm.sobar$q, mclust.sobar$d, pgmm.sobar$g, pgmm.sobar$model) #CUCU
adjustedRandIndex(lab.sobar, pgmm.sobar$map)
# HDDC
hddc.sobar$K
hddc.sobar$d
hddc.sobar$model
hddc.sobar$complexity
adjustedRandIndex(lab.sobar, hddc.sobar$class)

# Gopt
# GPMCs
mclust.sobar.Gopt <- Mclust(scale(x[, -20]), G = length(unique(lab.sobar)), control = emControl(itmax = 500))
mclust.sobar.Gopt$modelName
nMclustParams(mclust.sobar.Gopt$modelName, mclust.sobar.Gopt$d, mclust.sobar.Gopt$G)
adjustedRandIndex(lab.sobar, mclust.sobar.Gopt$classification)
# HDDC
hddc.sobar.Gopt <- hddc(scale(x[, -20]), K = length(unique(lab.sobar)), model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))
hddc.sobar.Gopt$model
hddc.sobar.Gopt$complexity
adjustedRandIndex(lab.sobar, hddc.sobar.Gopt$class)

################################################################################
################################ KIDNEY ########################################
data(ckd, package = "teigen")
x <- scale(ckd[, -1])
lab.kidney <- as.numeric(ckd[, 1])
G.max <- length(unique(lab.kidney)) + 2
if (dim(x)[2] < 5){m.max <- dim(x)[2]} else {m.max <- 5}

pugmm.kidney <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.kidney <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.kidney <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.kidney <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.kidney$G
pugmm.kidney$m
pugmm.kidney$model.name
pugmm.kidney$pm.free
adjustedRandIndex(lab.kidney, pugmm.kidney$label)
# GPCMs
mclust.kidney$G
mclust.kidney$modelName
nMclustParams(mclust.kidney$modelName, mclust.kidney$d, mclust.kidney$G)
adjustedRandIndex(lab.kidney, mclust.kidney$classification)
# PGMMs
pgmm.kidney$g
pgmm.kidney$q
pgmm.kidney$model
PGMM_dfree(pgmm.kidney$q, mclust.kidney$d, pgmm.kidney$g, pgmm.kidney$model) #CUU
adjustedRandIndex(lab.kidney, pgmm.kidney$map)
# HDDC
hddc.kidney$K
hddc.kidney$d
hddc.kidney$model
hddc.kidney$complexity
adjustedRandIndex(lab.kidney, hddc.kidney$class)

# Gopt
# PUGMMs
pugmm.kidney.Gopt <- pugmm(x, G = length(unique(lab.kidney)), m = 1:m.max)
pugmm.kidney.Gopt$m
pugmm.kidney.Gopt$model.name
pugmm.kidney.Gopt$pm.free
adjustedRandIndex(lab.kidney, pugmm.kidney.Gopt$label)
#GPCMs
mclust.kidney.Gopt <- Mclust(x, G = length(unique(lab.kidney)), control = emControl(itmax = 500))
mclust.kidney.Gopt$modelName
nMclustParams(mclust.kidney.Gopt$modelName, mclust.kidney.Gopt$d, mclust.kidney.Gopt$G)
adjustedRandIndex(lab.kidney, mclust.kidney.Gopt$classification)
# PGMMs
pgmm.kidney.Gopt <- pgmmEM(x, rG = length(unique(lab.kidney)), rq = 1:m.max, relax = TRUE)
pgmm.kidney.Gopt$q
pgmm.kidney.Gopt$model
PGMM_dfree(pgmm.kidney.Gopt$q, mclust.kidney.Gopt$d, pgmm.kidney.Gopt$g, pgmm.kidney.Gopt$model) #UUCU
adjustedRandIndex(lab.kidney, pgmm.kidney.Gopt$map)
# HDDC
hddc.kidney.Gopt <- hddc(x, K = length(unique(lab.kidney)), model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))
hddc.kidney.Gopt$model
hddc.kidney.Gopt$complexity
adjustedRandIndex(lab.kidney, hddc.kidney.Gopt$class)

################################ ECONOMICS #####################################
data(Economics, package = "datasetsICR")
x <- scale(Economics[, -13])
lab.economics <- as.numeric(Economics[, 13])
G.max <- length(unique(lab.economics)) + 2
if (dim(x)[2] < 5){m.max <- dim(x)[2]} else {m.max <- 5}

pugmm.economics <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.economics <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.economics <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.economics <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.economics$G
pugmm.economics$m
pugmm.economics$model.name
pugmm.economics$pm.free
adjustedRandIndex(lab.economics, pugmm.economics$label)
# GPCMs
mclust.economics$G
mclust.economics$modelName
nMclustParams(mclust.economics$modelName, mclust.economics$d, mclust.economics$G)
adjustedRandIndex(lab.economics, mclust.economics$classification)
# PGMMs
pgmm.economics$g
pgmm.economics$q
pgmm.economics$model
PGMM_dfree(pgmm.economics$q, mclust.economics$d, pgmm.economics$g, pgmm.economics$model) #UUU
adjustedRandIndex(lab.economics, pgmm.economics$map)
# HDDC
hddc.economics$K
hddc.economics$d
hddc.economics$model
hddc.economics$complexity
adjustedRandIndex(lab.economics, hddc.economics$class)

# Gopt
# PUGMMs
pugmm.economics.Gopt <- pugmm(x, G = length(unique(lab.economics)), m = 1:m.max)
pugmm.economics.Gopt$m
pugmm.economics.Gopt$model.name
pugmm.economics.Gopt$pm.free
adjustedRandIndex(lab.economics, pugmm.economics.Gopt$label)
# GPCMs
mclust.economics.Gopt <- Mclust(x, G = length(unique(lab.economics)), control = emControl(itmax = 500))
mclust.economics.Gopt$modelName
nMclustParams(mclust.economics.Gopt$modelName, mclust.economics.Gopt$d, mclust.economics.Gopt$G)
adjustedRandIndex(lab.economics, mclust.economics.Gopt$classification)
# HDDC
hddc.economics.Gopt <- hddc(x, K = length(unique(lab.economics)), model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))
hddc.economics.Gopt$model
hddc.economics.Gopt$complexity
adjustedRandIndex(lab.economics, hddc.economics.Gopt$class)

################################ TETRAGONULA ###################################
x <- read.csv2("DIRECTORY_NAME/Tetragonula.txt")
x <- x[, 2:6]
lab.tetragonula <- x[, 1]
G.max <- length(unique(lab.tetragonula)) + 2
if ((dim(x)[2] - 1) < 5){m.max <- dim(x)[2] - 1} else {m.max <- 5}

pugmm.tetragonula <- pugmm(x[, -1], G = 1:G.max, m = 1:m.max)
mclust.tetragonula <- Mclust(x[, -1], G = 1:G.max, control = emControl(itmax = 500))
pgmm.tetragonula <- pgmmEM(x[, -1], rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.tetragonula <- hddc(x[, -1], K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.tetragonula$G
pugmm.tetragonula$m
pugmm.tetragonula$model.name
pugmm.tetragonula$pm.free
adjustedRandIndex(lab.tetragonula, pugmm.tetragonula$label)
# GPCMs
mclust.tetragonula$G
mclust.tetragonula$modelName
nMclustParams(mclust.tetragonula$modelName, mclust.tetragonula$d, mclust.tetragonula$G)
adjustedRandIndex(lab.tetragonula, mclust.tetragonula$classification)
# PGMMs
# NULL (errors occur)
# HDDC
hddc.tetragonula$K
hddc.tetragonula$d
hddc.tetragonula$model
hddc.tetragonula$complexity
adjustedRandIndex(lab.tetragonula, hddc.tetragonula$class)
table(x[, 1], hddc.tetragonula$class)

# Gopt
# PUGMMs
pugmm.tetragonula.Gopt <- pugmm(x[, -1], G = length(unique(lab.tetragonula)), m = 1:m.max)
pugmm.tetragonula.Gopt$m
pugmm.tetragonula.Gopt$model.name
pugmm.tetragonula.Gopt$pm.free
adjustedRandIndex(lab.tetragonula, pugmm.tetragonula.Gopt$label)
# GPCMs
mclust.tetragonula.Gopt <- Mclust(x[, -1], G = length(unique(lab.tetragonula)), control = emControl(itmax = 500))
mclust.tetragonula.Gopt$modelName
nMclustParams(mclust.tetragonula.Gopt$modelName, mclust.tetragonula.Gopt$d, mclust.tetragonula.Gopt$G)
adjustedRandIndex(lab.tetragonula, mclust.tetragonula.Gopt$classification)
# PGMMs
pgmm.tetragonula.Gopt <- pgmmEM(x[, -1], rG = length(unique(lab.tetragonula)), rq = 1:m.max, relax = TRUE)
pgmm.tetragonula.Gopt$q
pgmm.tetragonula.Gopt$model
PGMM_dfree(pgmm.tetragonula.Gopt$q, mclust.tetragonula.Gopt$d, pgmm.tetragonula.Gopt$g, pgmm.tetragonula.Gopt$model)  #CUU
adjustedRandIndex(lab.tetragonula, pgmm.tetragonula.Gopt$map)
# HDDC
hddc.tetragonula.Gopt <- hddc(x[, -1], K = length(unique(lab.tetragonula)), model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))
hddc.tetragonula.Gopt$model
hddc.tetragonula.Gopt$complexity
adjustedRandIndex(lab.tetragonula, hddc.tetragonula.Gopt$class)

################################ DIABETES ######################################
data(diabetes, package = "mclust")
x <- scale(diabetes[, -1])
lab.diabetes <- as.numeric(diabetes[, 1])
G.max <- length(unique(lab.diabetes)) + 2
if (dim(x)[2] < 5){m.max <- dim(x)[2]} else {m.max <- 5}

pugmm.diabetes <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.diabetes <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.diabetes <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.diabetes <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.diabetes$G
pugmm.diabetes$m
pugmm.diabetes$model.name
pugmm.diabetes$pm.free
adjustedRandIndex(lab.diabetes, pugmm.diabetes$label)
# GPCMs
mclust.diabetes$G
mclust.diabetes$modelName
nMclustParams(mclust.diabetes$modelName, mclust.diabetes$d, mclust.diabetes$G)
adjustedRandIndex(lab.diabetes, mclust.diabetes$classification)
# PGMMs
pgmm.diabetes$g
pgmm.diabetes$q
pgmm.diabetes$model
PGMM_dfree(pgmm.diabetes$q, mclust.diabetes$d, pgmm.diabetes$g, pgmm.diabetes$model) #UCU
adjustedRandIndex(lab.diabetes, pgmm.diabetes$map)
# HDDC
hddc.diabetes$K
hddc.diabetes$d
hddc.diabetes$model
hddc.diabetes$complexity
adjustedRandIndex(lab.diabetes, hddc.diabetes$class)

# Gopt
# PUGMMs
pugmm.diabetes.Gopt <- pugmm(x, G = length(unique(lab.diabetes)), m = 1:m.max)
pugmm.diabetes.Gopt$m
pugmm.diabetes.Gopt$model.name
pugmm.diabetes.Gopt$pm.free
adjustedRandIndex(lab.diabetes, pugmm.diabetes.Gopt$label)

################################ BANKNOTES #####################################
data(banknote, package = "mclust")
x <- scale(banknote[, -1])
lab.banknotes <- as.numeric(banknote[, 1])
G.max <- length(unique(lab.banknotes)) + 2
if (dim(x)[2] < 5){m.max <- dim(x)[2]} else {m.max <- 5}

pugmm.banknotes <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.banknotes <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.banknotes <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.banknotes <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.banknotes$G
pugmm.banknotes$m
pugmm.banknotes$model.name
pugmm.banknotes$pm.free
adjustedRandIndex(lab.banknotes, pugmm.banknotes$label)
# GPCMs
mclust.banknotes$G
mclust.banknotes$modelName
nMclustParams(mclust.banknotes$modelName, mclust.banknotes$d, mclust.banknotes$G)
adjustedRandIndex(lab.banknotes, mclust.banknotes$classification)
# PGMMs
pgmm.banknotes$g
pgmm.banknotes$q
pgmm.banknotes$model
PGMM_dfree(pgmm.banknotes$q, mclust.banknotes$d, pgmm.banknotes$g, pgmm.banknotes$model) #CCUU
adjustedRandIndex(lab.banknotes, pgmm.banknotes$map)
# HDDC
hddc.banknotes$K
hddc.banknotes$d
hddc.banknotes$model
hddc.banknotes$complexity
adjustedRandIndex(lab.banknotes, hddc.banknotes$class)

# Gopt
# PUGMMs
pugmm.banknotes.Gopt <- pugmm(x, G = length(unique(lab.banknotes)), m = 1:m.max)
pugmm.banknotes.Gopt$m
pugmm.banknotes.Gopt$model.name
pugmm.banknotes.Gopt$pm.free
adjustedRandIndex(lab.banknotes, pugmm.banknotes.Gopt$label)
# GPCMs
mclust.banknotes.Gopt <- Mclust(x, G = length(unique(lab.banknotes)), control = emControl(itmax = 500))
mclust.banknotes.Gopt$modelName
nMclustParams(mclust.banknotes.Gopt$modelName, mclust.banknotes.Gopt$d, mclust.banknotes.Gopt$G)
adjustedRandIndex(lab.banknotes, mclust.banknotes.Gopt$classification)
# PGMMs
pgmm.banknotes.Gopt <- pgmmEM(x, rG = length(unique(lab.banknotes)), rq = 1:m.max, relax = TRUE)
pgmm.banknotes.Gopt$q
pgmm.banknotes.Gopt$model
PGMM_dfree(pgmm.banknotes.Gopt$q, mclust.banknotes.Gopt$d, pgmm.banknotes.Gopt$g, pgmm.banknotes.Gopt$model)  #UUCU
adjustedRandIndex(lab.banknotes, pgmm.banknotes.Gopt$map)
# HDDC
hddc.banknotes.Gopt <- hddc(x, K = length(unique(lab.banknotes)), model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))
hddc.banknotes.Gopt$model
hddc.banknotes.Gopt$complexity
adjustedRandIndex(lab.banknotes, hddc.banknotes.Gopt$class)

################################ CERAMIC #######################################
x <- read.csv("DIRECTORY_NAME/Ceramic.txt", header = FALSE)
lab.ceramic <- x[, 2]
G.max <- length(unique(lab.ceramic)) + 2
if ((dim(x)[2] -2) < 5){m.max <- dim(x)[2] - 2} else {m.max <- 5}

pugmm.ceramic <- pugmm(scale(x[, -c(1,2)]), G = 1:G.max, m = 1:m.max)
mclust.ceramic <- Mclust(scale(x[, -c(1,2)]), G = 1:G.max, control = emControl(itmax = 500))
pgmm.ceramic <- pgmmEM(scale(x[, -c(1,2)]), rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.ceramic <- hddc(scale(x[, -c(1,2)]), K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.ceramic$G
pugmm.ceramic$m
pugmm.ceramic$model.name
pugmm.ceramic$pm.free
adjustedRandIndex(lab.ceramic, pugmm.ceramic$label)
# GPCMs
mclust.ceramic$G
mclust.ceramic$modelName
nMclustParams(mclust.ceramic$modelName, mclust.ceramic$d, mclust.ceramic$G)
adjustedRandIndex(lab.ceramic, mclust.ceramic$classification)
# PGMMs
pgmm.ceramic$g
pgmm.ceramic$q
pgmm.ceramic$model
PGMM_dfree(pgmm.ceramic$q, mclust.ceramic$d, pgmm.ceramic$g, pgmm.ceramic$model) #UUU
adjustedRandIndex(lab.ceramic, pgmm.ceramic$map)
# HDDC
hddc.ceramic$K
hddc.ceramic$d
hddc.ceramic$model
hddc.ceramic$complexity
adjustedRandIndex(lab.ceramic, hddc.ceramic$class)

# Gopt
# PUGMMs
pugmm.ceramic.Gopt <- pugmm(scale(x[, -c(1,2)]), G = length(unique(lab.ceramic)), m = 1:m.max)
pugmm.ceramic.Gopt$m
pugmm.ceramic.Gopt$model.name
pugmm.ceramic.Gopt$pm.free
adjustedRandIndex(lab.ceramic, pugmm.ceramic.Gopt$label)

################################ AIS ###########################################
data(ais, package = "sn")
x <- scale(ais[, -c(1,2)])
lab.ais <- as.numeric(ais[, 1])
G.max <- length(unique(lab.ais)) + 2
if (dim(x)[2] < 5){m.max <- dim(x)[2]} else {m.max <- 5}

pugmm.ais <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.ais <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.ais <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.ais <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.ais$G
pugmm.ais$m
pugmm.ais$model.name
pugmm.ais$pm.free
adjustedRandIndex(lab.ais, pugmm.ais$label)
# GPCMs
mclust.ais$G
mclust.ais$modelName
nMclustParams(mclust.ais$modelName, mclust.ais$d, mclust.ais$G)
adjustedRandIndex(lab.ais, mclust.ais$classification)
# PGMMs
pgmm.ais$g
pgmm.ais$q
pgmm.ais$model
PGMM_dfree(pgmm.ais$q, mclust.ais$d, pgmm.ais$g, pgmm.ais$model) #UCU
adjustedRandIndex(lab.ais, pgmm.ais$map)
# HDDC
hddc.ais$K
hddc.ais$d
hddc.ais$model
hddc.ais$complexity
adjustedRandIndex(lab.ais, hddc.ais$class)

# Gopts
# PUGMMs
pugmm.ais.Gopt <- pugmm(x, G = length(unique(lab.ais)), m = 1:m.max)
pugmm.ais.Gopt$m
pugmm.ais.Gopt$model.name
pugmm.ais.Gopt$pm.free
adjustedRandIndex(lab.ais, pugmm.ais.Gopt$label)
# GPCMs
mclust.ais.Gopt <- Mclust(x, G = length(unique(lab.ais)), control = emControl(itmax = 500))
mclust.ais.Gopt$modelName
nMclustParams(mclust.ais.Gopt$modelName, mclust.ais.Gopt$d, mclust.ais.Gopt$G)
adjustedRandIndex(lab.ais, mclust.ais.Gopt$classification)
# PGMMs
pgmm.ais.Gopt <- pgmmEM(x, rG = length(unique(lab.ais)), rq = 1:m.max, relax = TRUE)
pgmm.ais.Gopt$q
pgmm.ais.Gopt$model
PGMM_dfree(pgmm.ais.Gopt$q, mclust.ais.Gopt$d, pgmm.ais.Gopt$g, pgmm.ais.Gopt$model) #UUCU
adjustedRandIndex(lab.ais, pgmm.ais.Gopt$map)
# HDDC
hddc.ais.Gopt  <- hddc(x, K = length(unique(lab.ais)), model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))
hddc.ais.Gopt$model
hddc.ais.Gopt$complexity
adjustedRandIndex(lab.ais, hddc.ais.Gopt$class)

################################ COFFEE ########################################
data(coffee, package = "pgmm")
x <- scale(coffee[, -c(1, 2)])
lab.coffee <- coffee[, 1]
G.max <- length(unique(lab.coffee)) + 2
if (dim(x)[2] < 5){m.max <- dim(x)[2]} else {m.max <- 5}

pugmm.coffee <- pugmm(x, G = 1:G.max, m = 1:m.max)
mclust.coffee <- Mclust(x, G = 1:G.max, control = emControl(itmax = 500))
pgmm.coffee <- pgmmEM(x, rG = 1:G.max, rq = 1:m.max, relax = TRUE)
hddc.coffee <- hddc(x, K = 1:G.max, model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))

# Results
# PUGMMs
pugmm.coffee$G
pugmm.coffee$m
pugmm.coffee$model.name
pugmm.coffee$pm.free
adjustedRandIndex(lab.coffee, pugmm.coffee$label)
# GPCMs
mclust.coffee$G
mclust.coffee$modelName
nMclustParams(mclust.coffee$modelName, mclust.coffee$d, mclust.coffee$G)
adjustedRandIndex(lab.coffee, mclust.coffee$classification)
# PGMMs
pgmm.coffee$g
pgmm.coffee$q
pgmm.coffee$model
PGMM_dfree(pgmm.coffee$q, mclust.coffee$d, pgmm.coffee$g, pgmm.coffee$model) #CCUU
adjustedRandIndex(lab.coffee, pgmm.coffee$map)
# HDDC
hddc.coffee$K
hddc.coffee$d
hddc.coffee$model
hddc.coffee$complexity
adjustedRandIndex(lab.coffee, hddc.coffee$class)

# Gopt
# PUGMMs
pugmm.coffee.Gopt <- pugmm(x, G = length(unique(lab.coffee)), m = 1:m.max)
pugmm.coffee.Gopt$m
pugmm.coffee.Gopt$model.name
pugmm.coffee.Gopt$pm.free
adjustedRandIndex(lab.coffee, pugmm.coffee.Gopt$label)
# GPCMs
mclust.coffee.Gopt <- Mclust(x, G = length(unique(lab.coffee)), control = emControl(itmax = 500))
mclust.coffee.Gopt$modelName
nMclustParams(mclust.coffee.Gopt$modelName, mclust.coffee.Gopt$d, mclust.coffee.Gopt$G)
adjustedRandIndex(lab.coffee, mclust.coffee.Gopt$classification)
# HDDC
hddc.coffee.Gopt <- hddc(x, K = length(unique(lab.coffee)), model = "ALL", itermax = 500, mc.cores = 8, kmeans.control = list(nstart = 100))
hddc.coffee.Gopt$model
hddc.coffee.Gopt$complexity
