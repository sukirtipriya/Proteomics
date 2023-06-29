#Load the package

library(DEqMS)

# Download and Read the input protein table

url <- "https://ftp.ebi.ac.uk/pride-archive/2016/06/PXD004163/Yan_miR_Protein_table.flatprottable.txt"
download.file(url, destfile = "./miR_Proteintable.txt",method = "auto")

df.prot = read.table("miR_Proteintable.txt",stringsAsFactors = FALSE,
                     header = TRUE, quote = "", comment.char = "",sep = "\t")

#Extract quant data columns for DEqMS

TMT_columns = seq(15,33,2)
dat = df.prot[df.prot$miR.FASP_q.value<0.01,TMT_columns]
rownames(dat) = df.prot[df.prot$miR.FASP_q.value<0.01,]$Protein.accession

dat.log = log2(dat)
#remove rows with NAs
dat.log = na.omit(dat.log)

#Use boxplot to check if the samples have medians centered. if not, do median centering.

boxplot(dat.log,las=2,main="TMT10plex data PXD004163")

#Make design table.

cond = as.factor(c("ctrl","miR191","miR372","miR519","ctrl",
"miR372","miR519","ctrl","miR191","miR372"))

# The function model.matrix is used to generate the design matrix
design = model.matrix(~0+cond) # 0 means no intercept for the linear model
colnames(design) = gsub("cond","",colnames(design))

#Make contrasts
x <- c("miR372-ctrl","miR519-ctrl","miR191-ctrl",
       "miR372-miR519","miR372-miR191","miR519-miR191")
contrast =  makeContrasts(contrasts=x,levels=design)
fit1 <- lmFit(dat.log, design)
fit2 <- contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)

#DEqMS analysis

library(matrixStats)
count_columns = seq(16,34,2)
psm.count.table = data.frame(count = rowMins(
  as.matrix(df.prot[,count_columns])), row.names =  df.prot$Protein.accession)
fit3$count = psm.count.table[rownames(fit3$coefficients),"count"]
fit4 = spectraCounteBayes(fit3)

#Visualize the fit curve - variance dependence on quantified PSM
VarianceBoxplot(fit4,n=30,main="TMT10plex dataset PXD004163",xlab="PSM count")

VarianceScatterplot(fit4,main="TMT10plex dataset PXD004163")

#Extract the results as a data frame and save it.
DEqMS.results = outputResult(fit4,coef_col = 1)
#if you are not sure which coef_col refers to the specific contrast,type
head(fit4$coefficients)

head(DEqMS.results)

# Save it into a tabular text file
write.table(DEqMS.results,"DEqMS.results.miR372-ctrl.txt",sep = "\t",
            row.names = F,quote=F)
