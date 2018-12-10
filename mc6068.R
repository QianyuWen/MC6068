devtools::create("MC6068")
devtools::document()

library(jpeg)
library(gplots)
figure1 <- readJPEG("mona.jpg")
mona1 <- figure1[,,1]
heatmap.2(mona1, Rowv = NA, 
          Colv = NA, key = F, scale = "none", 
          col = grey.colors,
          trace = "none",
          labRow = "", labCol = ""
)
n <- length(mona1)
set.seed(32607)
naindex <- sample(n,0.2*n)
mona1_missing <- mona1
mona1_missing[naindex] <- NA
heatmap.2(mona1_missing, Rowv = NA, 
          Colv = NA, key = F, scale = "none", 
          col = grey.colors,
          trace = "none",
          labRow = "", labCol = "")
mona <- as.data.frame(mona1_missing)
devtools::use_data(mona, overwrite = T)

devtools::build()

install.packages(file.path("C:/Users/shpwen/Desktop/MC6068_0.0.0.9000.tar.gz")
                 ,repos=NULL,type="source")
library(MC6068)
remove.packages("MC6068")

help(package="MC6068")
data(mona)

cm <- AdmmMC(x=as.matrix(mona), silence = T)
heatmap.2(x=cm, Rowv = NA, 
          Colv = NA, key = F, scale = "none", 
          col = grey.colors,
          trace = "none",
          labRow = "", labCol = "")
