#code i don't think i need
#Cmatrix.diag = Cmatrix #store a copy before removing diag
#diag(Cmatrix) = rep(0,47)
#cited/citing
#cited = rowSums(Cmatrix)
#citing = colSums(Cmatrix)
#cite.ratio = cited/citing
#normalized
#Cmatrix.norm = t(t(Cmatrix)/citing) #entries are percentage of j's (columnn's) citations to i (row)
#Cnet.norm = as.network(Cmatrix.norm, directed=T, matrix.type="a", ignore.eval=F,  names.eval="citations")
#------


#Code from Carlen and Handcock's "Discussion on 'Statistical Modelling of Citation Exchange Between Statistics
#Journals’ by Cristiano Varin, Manuela Cattelan and David Firth"


# Pre-Processing -------

library(latentnet)
library(ergm.count)
library(stats)

download.file("http://cristianovarin.weebly.com/uploads/1/5/1/5/15156956/varin_cattelan_firth_supplement.zip", "varin_cattelan_firth_supplement.zip")
unz("varin_cattelan_firth_supplement.zip", "Varin_Cattelan_Firth_supplement")
Cmatrix <- as.matrix(read.csv("Varin_Cattelan_Firth_supplement/Data/cross-citation-matrix.csv", row.names = 1))

#as valued net (#see Modeling valued networks with statnet paper)
#removes self-citations; Stigler method also excluded self-citations
Cnet = as.network(Cmatrix, directed=T, matrix.type="a", ignore.eval=F,  names.eval="citations")

# Pre-Processing 2: Code from from Varin_Cattelan_Firth_supplement.R published with article -------

## From Section 1 & 2

journal.abbr <- rownames(Cmatrix)
Tmatrix <- Cmatrix + t(Cmatrix)
diag(Tmatrix) <- diag(Cmatrix)
journals.cluster <- hclust(d = as.dist(1 - cor(Tmatrix)))

## From Section 3: Quasi-Stigler Model
library(BradleyTerry2)

Cdata <- countsToBinomial(Cmatrix)
fit <- BTm(outcome = cbind(win1, win2), player1 = player1, player2 = player2, data = Cdata)

npairs <- NROW(Cdata)
njournals <- nlevels(Cdata$player1)
phi <- sum(residuals(fit, "pearson")^2) / (npairs - (njournals - 1))

### "3.1 Journal residuals"
journal.res <- rep(NA, njournals)
res <- residuals(fit, type = "pearson")
coefs <- c(0, coef(fit)) # 0 is the coefficient of the first journal
for(i in 1:njournals){
  A <- which(Cdata$player1 == journal.abbr[i])
  B <- which(Cdata$player2 == journal.abbr[i])
  y <- c(res[A], -res[B])
  x <- c(-coefs[Cdata$player2[A]], -coefs[Cdata$player1[B]])
  journal.res[i] <- sum(y * x) / sqrt(phi * sum(x ^ 2))
}
names(journal.res) <- journal.abbr

library(qvcalc)
cov.matrix <- matrix(0, nrow = njournals, ncol = njournals)
cov.matrix[-1, -1] <- vcov(fit)
qse <- qvcalc(phi * cov.matrix , estimates = c(0, coef(fit)),
              labels = journal.abbr)

export.scores <- qse$qvframe$estimate
export.scores <- export.scores - mean(export.scores)
names(export.scores) <- journal.abbr

sort.id <- sort(export.scores, decreasing = TRUE, index.return = TRUE)$ix
fit.table <- data.frame(quasi = export.scores[sort.id], qse = qse$qvframe$quasiSE[sort.id])

## Some post-processing/alphabetizing (not from author supplement)

rownames(fit.table)
rownames(fit.table)[c(1,6,20)] = c("JRSS.B", "JRSS.A", "JRSS.C")
match(rownames(fit.table),Cnet%v%"vertex.names")
fit.table2 = fit.table[order(match(rownames(fit.table),Cnet%v%"vertex.names")),]

# Network model 1: Sender-receiver model with citation counts ~ Poisson -------

latent.srp1 = ergmm(Cnet~ sender(base=0) + receiver(base=0) - 1, response = "citations",
                    family="Poisson.log", control = control.ergmm(pilot.runs=1), seed=111)

# Correlation of output with Stigler (0.95)

cor(latent.srp1$mcmc.mle$beta[48:94] - latent.srp1$mcmc.mle$beta[1:47], fit.table2$quasi)


# Network model 2: Two-dimensional latent-space model  -------

latent.srp2 = ergmm(Cnet~euclidean(d=2) + sender(base=0) + receiver(base=0) - 1, response = "citations", 
                    family="Poisson.log", seed=123, control=ergmm.control(interval=200, sample.size=10000, burnin=100000)) 

# Correlation of output with Stigler (0.99)

cor(latent.srp2$mcmc.mle$beta[48:94] - latent.srp2$mcmc.mle$beta[1:47], fit.table2$quasi)

# Network Plot 1 (Fig. 12a) -------

plot(latent.srp2, labels = T, cex=.7, label.cex=.5,
     plot.vars=F, vertex.col = cutree(journals.cluster, h = 0.6)+3, 
     edge.col=0, print.formula=F, main = NA, vertex.border = 0,
     vertex.cex = vc2, label.pos=3, suppress.axes=T,  xlab=NA, ylab=NA)

# Network Plot 2 (Fig. 12b) -------

# cloud/uncertainty plot

p = latent.srp2[["sample"]][["Z"]]
m = matrix(0,0,3)
for (i in 1:47) {
  p1 = cbind(p[,i,][,1],p[,i,][,2], rep((cutree(journals.cluster, h = 0.6)+3)[i],10000))
  m = rbind(m,p1)
}

#null plot:
plot.ergmm(latent.srp2, xlab = NA, ylab = NA, vertex.col=0, edge.col=0, 
           plot.vars=F, suppress.axes=T, print.formula=F, vertex.border=0,
           xlim = c(-4,4), ylim = c(-3,5), main = NA, pie=F)

#add points
s = sample(1:470000)
for (i in 1:47000) {
  points(m[s[i],1], m[s[i],2], col = m[s[i],3], pch=".")
}

legend("topleft", col = 4:11, pch = 16, legend=c("review", "general", "theory/method", "apply/health", "computational","eco/envir", "JSS", "StataJ"), cex=.8, box.col=0)

