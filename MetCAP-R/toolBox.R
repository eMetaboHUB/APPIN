# This file is part of [MetCAP-R]
# 
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software. You can use, 
# modify and/or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
# 
# As a counterpart to the access to the source code and rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty and the software's author, the holder of the
# economic rights, and the successive licensors have only limited
# liability.
# 
# See the LICENSE file for more details.

## This script contains several functions used for the database creation and matching Shiny app.
## Author: jfmartin



###### Function to plot the best correlations between the variables of 2 matrix A and B 
### Suppress missing data by pairs
plotrsup = function(A, B, rmin, method="pearson")
{
  for(i in 1:dim(A)[2]) {
    if(is.numeric(A[,i]))
      for(j in 1:dim(B)[2]) {
        if(is.numeric(B[,j])) 
        {
          AnNa <- A[!is.na(A[,i]) & !is.na(B[,j]),i]
          BnNa <- B[!is.na(A[,i]) & !is.na(B[,j]),j]
          r = cor(AnNa,BnNa, method = method)
          if(abs(r) >= rmin && r<0.99999)
          {
            plot(AnNa, BnNa, xlab = dimnames(A)[[2]][i], ylab = dimnames(B)[[2]][j])
            abline(lsfit(AnNa, BnNa))
            title(paste("r=", as.character(round(r, 2))))
          }
        }
      }
  }
}

###### Funtion to write in a dataframe the best correlations between the variables of 2 matrix A and B 
### Suppress missing data by pairs
write.rsup = function(A, B, rmin)
{
  resr <- array(0,c(dim(A)[2]*dim(B)[2],3)); it <- 1
  for(i in 1:dim(A)[2]) {
    if(is.numeric(A[,i]))
      for(j in 1:dim(B)[2]) {
        if(is.numeric(B[,j])) 
        {
          AnNa=A[!is.na(A[,i]) & !is.na(B[,j]),i]
          BnNa=B[!is.na(A[,i]) & !is.na(B[,j]),j]
          r = cor(AnNa,BnNa)
          if(abs(r) >= rmin && r<0.99999)
          {
            resr[it,1] <- dimnames(A)[[2]][i]
            resr[it,2] <- dimnames(B)[[2]][j]
            resr[it,3] <- r
            it <- it +1
          }
        }
      }
  }
  resr <- resr[resr[,1]!=0,]
  return(resr)
}

###### Plot A,B  with different symbols depending the s factor
plotby = function(s, A, B, lx ="", ly ="", reg=T ,regby=F, leg=T ,rmin)
{
if(!is.factor(s))  stop("Le premier argument doit etre un FACTEUR")
if(!is.numeric(A)) stop("Le deuxieme argument doit etre NUMERIQUE")
if(!is.numeric(B)) stop("Le troisieme argument doit etre NUMERIQUE")

   AnNa=A[!is.na(A)& !is.na(B)]
   BnNa=B[!is.na(A)& !is.na(B)]
   snNa=s[!is.na(A)& !is.na(B)]
   D <- data.frame(f=snNa, x=AnNa, y=BnNa)
   r <- cor (D$x, y=D$y)
   if(abs(r) >= rmin && r<0.99999) {

      plot (D$x , D$y , type="n", xlab=lx , ylab=ly )
      if ( reg==T ) abline(lsfit(D$x , D$y), lty=1)
      title (paste("r=", as.character(round(r, 2))))
      nl <- length(levels(D$f))
      cl=c(1,17,16,2,3,4); sy=c(1,16,17,3,5,12)
      for (i in 1:nl) {
        sD <- D[ D$f==levels(D$f)[i],]
        points (sD$x , sD$y, pch=sy[i], col=cl[i])
        if ( regby==T ) abline(lsfit( sD$x , sD$y), lty=i+1)
      }
      if (leg==T) {
         py <- max(D$y)
         if (r >= 0) px <- min(D$x)
         if (r <  0) px <- min(D$x)+(max(D$x)-min(D$x))/2
         legend(px, py , levels(D$f),pch=sy[1:nl])
      }
   }
}


###### Plot the boxplots of the variables from a dataframe containing observations at +/- nsd Standard deviation
plot_outlier = function (ids, by, nsd)
{
for (i in 1:dim(ids)[2])
  if (is.numeric(ids[,i])) 
  {
    mini=min(ids[!is.na(ids[,i]),i])
    maxi=max(ids[!is.na(ids[,i]),i])
    moy=mean(ids[!is.na(ids[,i]),i])
    ect=sd(ids[!is.na(ids[,i]),i])  
    if (mini < moy- nsd*ect | maxi > moy+nsd*ect) boxplot (ids[,i]~ by, ylab=dimnames(ids)[[2]][i]) 
  }
}

acpEllipse=function(ids,fact,firstvar,lastvar=dim(ids)[2],scaling="uv",meth="svd",plotloading=T,tc=0.66) {
	# Make a PCA on ids knowing column 1 contains the individual indification
	# tc  is the tall of character used for the plots
	# Column 2 contains the factor defining individuals colors 
	# The next columns 3:dim(ids)[2] are the quantitative data used for the PCA
	# The PCA is realized with the scaling "uv","none","pareto"
	# The PCA algorithm is : nipals (if missing data) or svdImpute or Bayesian or svd or probalistic
	# The PCA loading is plotted if plotloading==T (default value)
  require(ade4)
	classe=factor(ids[[fact]])
	idions=dimnames(ids[,firstvar:lastvar])[2]
	colour=1:length(levels(classe))
	ions=as.matrix(ids[,c(firstvar:lastvar)])
	# Choice of scaling : "uv","none","pareto"
	object=prep(ions, scale=scaling, center=TRUE)
	# ALGO: nipals,svdImpute, Bayesian, svd, probalistic=F
	result <- pca(object, center=F, method=meth, nPcs=2)
	# ADE4 : individuals ellispoids representation for each class
	s.class(result@scores, classe, cpoint = 1,xax=1,yax=2,col=colour,sub=sprintf("Scores - PCs %sx%s",1,2), possub="bottomright")
	if (plotloading==T) 
		{s.label(result@loadings,cpoint = 0,boxes=F,clabel=tc, xax=1,yax=2,sub="Loadings",possub="bottomright")}
}

acplight=function(ids, scaling="uv") {
	# Make a PCA on ids knowing the column 1 contains the individuals identifications
	# The column 2:nf contains the factors defining the individuals colors
	for (i in 1:3) {
		# i=1 
		idss=data.frame(ids[ids[,i+1]!="",])
		classe=as.factor(idss[[i+1]])
		idsample=as.character(idss[[1]])
		colour=1:length(levels(classe))
		ions=as.matrix(idss[,5:dim(idss)[2]])
		# Choice of scaling : "uv","none","pareto"
		object=prep(ions, scale=scaling, center=TRUE)
		# ALGO: nipals,svdImpute, Bayesian, svd, probalistic=F
		result <- pca(object, center=F, method="svd", nPcs=2)
		# ADE4 : individuals ellispoids representation for each class
		s.class(result@scores, classe, cpoint = 1,xax=1,yax=2,col=colour,sub=sprintf("Scores - PCs %sx%s",1,2), possub="bottomright")
		#s.label(result@loadings,label = ions, cpoint = 0, clabel=0.4, xax=1,yax=2,sub="Loadings",possub="bottomright")
		}
}

# Program to plot ions with highlighting of factors 

plotInterCI = function (ids,indf1,indf2,firstvar) {
   # Function to plot variables (ions) with two cross factors. If kinetics use it as factor 2 
   # As entry, ids input dataframe, indf1 and indf2 index in the 2 factors dataframe
   # ids=x;i=25;indf1=6;indf2=10;indSubject=4
   lastvar=dim(ids)[2]
   for (i in firstvar:lastvar) { 
      
      par(las=2)
      xp=ids[!is.na(ids[[i]]),c(indf1,indf2,i)];
      xp=data.frame(xp,interaction(xp[[1]],xp[[2]],sep="."));
      dimnames(xp)[[2]][4]="Inter";
      levInter=c(1:length(levels(xp[[4]])));nl1=length(levels(xp[[1]]));nl2=length(levels(xp[[2]]));
      levInter=array(data=levInter,dim=c(nl1,nl2))
      connection=vector("list",nl1)
      for (i in 1:nl1) {connection[[i]]=levInter[i,]}
      plotmeans(xp[[3]]~ xp[[4]],
                connect=connection,
                ccol=c(1:nl1), barwidth=1.5, pch=16,
                main=dimnames(xp)[[2]][3],xlab=" ",	
                ylab="Intensity"
      );
      #	plot.design(xp[[3]]~xp[[1]]+xp[[2]]+xp[[4]],fun=mean,main=dimnames(xp)[[2]][3],xlab=" ",ylab=" ");
      # 	boxplot(xp[[3]]~xp[[4]],ylab="Intensity");
   }
}


BoxHistDesign = function (ids, fact, firstvar,lastvar=dim(ids)[2]) {
# ids input dataset; fact type vector containing index of studied factors; firstvar 1st quantitative variable
# Function that for a dataset ids plot a histogramm of each quantitative variable with MINIMUM 2 FACTORS 
# Plot a boxplot by factor fact[i]
# Plot the plo.design of factors
# If the was a kinetic, put the index in first in the vector index "fact"
#	lastvar=dim(ids)[2]
	nbfact=length(fact); nbinter=0.5*(nbfact^2)-0.5*nbfact; gco=c(1,2,3,3,4);gli=c(3,3,4,6,6)
	for (i in firstvar:lastvar) {
		if (nbfact < 6) { par(mfrow=c(gli[nbfact],gco[nbfact]))} else { par(mfrow=c(nbfact-1,nbinter-1))}
		# cat("var ",i," -> ",names(ids)[i],"\n")
		# Suppression of missing data 
		ods=ids[!is.na(ids[[i]]),]; lab=names(ods)[i]
		xpf=ods[,fact]
		xpv=ods[,i]
#		xpi=array(data="inter",nrow=dim(ods)[1],ncol=nbinter)	
		xb=data.frame(xpf,xpv); indvar=dim(xb)[2];lab=names(ods)[i];names(xb)[indvar]=lab
		hist(xb[[indvar]],main=lab,xlab=" ",ylab=" ")
		plot.design(xb,fun=mean,main=" ",xlab=" ",ylab=" ")
		for (b in 1:nbfact) {
			boxplot(xb[[indvar]]~xb[[b]],ylab="Intensity")
		}
		## Plot the interactions in function of the number of factors
		if (nbfact>1) {
			y=0
			for (j in 1:(nbfact-1)){ for (u in (j+1):nbfact){
				y=y+1
				xpi=interaction(ods[[fact[j]]],ods[[fact[u]]],sep="*")	
				xp=data.frame(xpf,xpi,xpv);
				indinter=length(fact)+1; names(xp)[indinter]="interaction";		
				indvar=dim(xp)[2]; names(xp)[indvar]=names(ods)[i];
				indinter=length(fact)+1; names(xp)[indinter]="interaction";
				boxplot(xp[[indvar]]~xp[[indinter]],ylab="Intensity")
				interaction.plot(xp[[j]],xp[[u]],xp[[indvar]])	
			}}
		}
	}

}

BoxHistDesign2 = function (ids, fact, firstvar,lastvar=dim(ids)[2]) {
   # ids input dataset; fact type vector containing the index of studied factors; firstvar 1st quantitative variable
   # Function that for a dataset ids plot a histogramm of each quantitative variable with MINIMUM 2 FACTORS
   # Plot a boxplot by factor fact[i]
   # Plot the plo.design of factors
  # If the was a kinetic, put the index in first in the vector index "fact"
   #	lastvar=dim(ids)[2]
   nbfact=length(fact); nbinter=0.5*(nbfact^2)-0.5*nbfact; gco=c(1,2,3,3,4);gli=c(3,3,4,6,6)
   for (i in firstvar:lastvar) {
      if (nbfact < 6) { par(mfrow=c(gli[nbfact],gco[nbfact]))} else { par(mfrow=c(nbfact-1,nbinter-1))}
      # cat("var ",i," -> ",names(ids)[i],"\n")
      # Suppress the missing data 
      ods=ids[!is.na(ids[[i]]),]; lab=names(ods)[i]
      xpf=ods[,fact]
      xpv=ods[,i]
      #		xpi=array(data="inter",nrow=dim(ods)[1],ncol=nbinter)	
      xb=data.frame(xpf,xpv); indvar=dim(xb)[2];lab=names(ods)[i];names(xb)[indvar]=lab
      hist(xb[[indvar]],main=lab,xlab=" ",ylab=" ")
      plot.design(xb,fun=mean,main=" ",xlab=" ",ylab=" ")
      for (b in 1:nbfact) {
         boxplot(xb[[indvar]]~xb[[b]],ylab="Intensity")
      }
    }
}


deter_ioni <- function (aninfo, pm)
{
  # Determine the ionisation of ProbMetab annotation
  # If the difference m/z and proposed mass by probmetab ~1 we use the sign of difference to define the ionisation
  # If adducts diff >>1 we search chains "+" or "2+" or "3+" or "-"... to define the ionisation
  # aninfo : vector resulting from parse (seprator #) from annotation field of ProbMetab
  if (round(abs(as.numeric(aninfo[1]) - pm),0) ==1) {
    if (as.numeric(aninfo[1]) - pm <0) {esi <- "n"} else {esi <- "p"}
  } else 
    if (!is.na(aninfo[4])) {
      anstr <- aninfo[4]
      cat(anstr)
      if ((grepl("]+",anstr,fixed=T)==T) || (grepl("]2+",anstr,fixed=T)==T) || (grepl("]3+",anstr,fixed=T)==T)) { esi <- "p"}
      else 
        if ((grepl("]-",anstr,fixed=T)==T) || (grep("]2-",anstr,fixed=T)==T) || (grep("]3-",anstr,fixed=T)==T)) { esi <- "n"}
      cat(" ioni ",esi,"\n")
    } else
    { esi <- "u"} 
  
  return(esi)
}


##################### fonction VennSignif #####################
## compute common signif variables in a varMedata file with columns of significant variable define with 1 for significant
## And draw a Venn diagram. Able to draw for only 2 or 3 factors.

vennSignif <- function(ids,sigsub,siglab,sigcol,chx=1.5)
   ## ids : input dataframe
   ## nv : number of factor to compute and draw as a Venn diagram : 2 or 3
   ## sigsub : vector of factor subscrip
   ## siglab : vector of names of factors
   ## sigcol : vector of color filling of the circles
{
   nv <- length(sigsub)
   if (nv==2) {
      n1 <- nrow(subset(ids, ids[[sigsub[1]]]==1))
      n2 <- nrow(subset(ids, ids[[sigsub[2]]]==1))
      cat(colnames(ids)[sigsub[1]]," & ",colnames(ids)[sigsub[2]],"\n")
      n12 <- nrow(subset(ids, (ids[[sigsub[1]]]==1 & ids[[sigsub[2]]]==1)))
      
      draw.pairwise.venn(area1 = n1, area2 = n2, cross.area = n12, category = siglab,
                         cex=chx, cat.cex = chx,
                         fill = sigcol[1:2], cat.pos = c(0,0), cat.dist = rep(0.025, 2))
   } else 
      if (nv==3) {
         n1 <- nrow(subset(ids, ids[[sigsub[1]]]==1))
         n2 <- nrow(subset(ids, ids[[sigsub[2]]]==1))
         n3 <- nrow(subset(ids, ids[[sigsub[3]]]==1))
         cat(colnames(ids)[sigsub[1]]," & ",colnames(ids)[sigsub[2]]," & ",colnames(ids)[sigsub[3]],"\n")
         n12 <- nrow(subset(ids, (ids[[sigsub[1]]]==1 & ids[[sigsub[2]]]==1)))
         n13 <- nrow(subset(ids, (ids[[sigsub[1]]]==1 & ids[[sigsub[3]]]==1)))
         n23 <- nrow(subset(ids, (ids[[sigsub[2]]]==1 & ids[[sigsub[3]]]==1)))
         
         n123 <- nrow(subset(ids, (ids[[sigsub[1]]]==1 & ids[[sigsub[2]]]==1 & ids[[sigsub[3]]]==1)))
         
         draw.triple.venn(area1 = n1, area2 = n2, area3 = n3, n12,n23,n13,n123,
                          cex=chx, cat.cex = chx, 
                          category = siglab, lty = "blank", 
                          fill = sigcol[1:3])
         
      } else stop("Only 2 ou 3 factors are supported")
}

crea3files <- function(vmfile, vipfile, splfile, dmfile, outsuf) {
  # This function creates the W4M 3 metadata + matrix files for univariate statistics corresponding with significant VIP 
  # vmfile: variableMetadata file with xcms xtract infos(annotation, ions in rows, columns= mz CAMERA intensities annotations
  # vipfile : file of ions + vip n. These ions must be a subset of ions of vmfile, last column define the signifcance
  # splfile : file of sample(rows) + metadata samples. These samples must be a subset of samples of vmfile
  # dmfile : file of data matrix to be used for univariate test.
  # outsuf : out file common prefix name to create the 3 differents files
  
  ## read variable metadata file
  VM <- read.table(file=vmfile,header=TRUE, sep="\t",stringsAsFactors =FALSE)
  ## order by ions id which must be in the first column
  VM <- VM[order(VM[[1]]),]
  
  VIP <- read.table(file=vipfile,header=TRUE, sep="\t",stringsAsFactors =FALSE)
  VIP <- VIP[order(VIP[[1]]),]
  if (sum(VIP[,1] != VM[,1]) > 0) stop("ions id are different in vip and annot files")
  
  ## selection of signifcant VIP ions, last column of VIP file
  VIPsignif <- VIP[VIP[ncol(VIP)] == 1,] 
  
  ## merge to keep only signifcant ions
  VMs <- merge(x = VM, y=VIPsignif, by.x=1, by.y =1, all.x=FALSE, all.y=TRUE)
  ionsSignif <- VIPsignif[[1]]
  
  ## samples avec uniquement les sampleType ="sample"
  SPL <- read.table(splfile, header=TRUE,sep="\t",stringsAsFactors = FALSE)
  SPLs <- SPL[SPL$sampleType=="sample",]
  splist <- SPLs[[1]]
  
  ## dataMatrix avec les ions signifs et les sampleType ="sample" uniquement et class�s 
  ## read table data matrix
  DM <- read.table(file = dmfile, header = TRUE, sep="\t")
  DM <- DM[order(DM[[1]]),]
  ## on conserve les ions signifs
  DMsignif <- DM[DM[[1]] %in% ionsSignif,]
  DMSS <- data.frame(DMsignif[[1]],DMsignif[,colnames(DMsignif) %in% splist])
  
  ################################## create output files  ##################################
  
  write.table(DMSS,file=paste(outsuf,"DM.txt",sep=""), sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  write.table(SPLs,file=paste(outsuf,"SPL.txt",sep=""),sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  write.table(VMs, file=paste(outsuf,"VM.txt",sep=""), sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
  
  #return(annotAll)
}


moyByFactor <- function(imtx, ispl, subfact, onlyspl=TRUE, oper="mean") {
  
  ## imtx correspond to W4M data matrix : var in rows Id var in the first columns, samples in columns given by vector isample
  ## ispl correspond to W4M samples metadata : samples in rows (same as the columns of imtx data matrix)
  ## The function computes the mean for each level of factor in column #subact, present in the ispl metadata sample dataframe
  ## if onlyspl = TRUE, the average (or other) is computed only on samples defined as "sample" in the sampleType
  ## for each variables of which subscripts given in ivar vector) of the imtx dataframe.
  ## returns a data matrix of means with variales in rows and levels of factor in columns
  
  ## selection of biological sample sampleType=="sample"
  ## in the sampleMetadata ispl
  if (onlyspl) {ispl <- ispl[ispl$sampleType == "sample",]}
  
  ## Selection of biological sample in the dataMatrix by merging with ispl
  varName <- colnames(imtx)[1]
  rownames(imtx) <- imtx[[1]]
  imtx <- imtx[,-1]
  Timtx <- t(imtx)
  Timtx <- data.frame(row.names(Timtx),Timtx)
  ids <- merge(x=ispl, y=Timtx, by.x=1, by.y=1, all.x=TRUE, all.y=FALSE)
  
  
  ## compute mean by levels of factors diet
  
  firstion <- ncol(ispl)+1
  lastion <- ncol(ids)
  nbions <- lastion-firstion+1
  levnames <- levels(as.factor(ispl[[subfact]]))
  nlevels <- length(levnames)
  avInt      <- array(0,dim=c(nlevels,nbions))
  colnames(avInt) <- colnames(ids)[c(firstion:lastion)]
  row.names(avInt) <- levnames
  
  for (i in 1:nbions) {avInt[,i] <- tapply(ids[[(i+firstion-1)]], ids[[subfact]], FUN=oper)}
  
  TavInt <- t(avInt)
  resMoy <- data.frame(row.names(TavInt),TavInt)
  colnames(resMoy)[1] <- varName
  return(resMoy)
  
}


foldchange <- function(ids,reflev,collev=c(2:ncol(ids)), fctype=c("log2","simple")[1]) {
  ## computes fold change between all levels of a factor 
  ## ids : an input dataframe with variables in rows and 
  ## id of variable in 1st column and average of the level of a factor in the other columns
  ## this matrix is computed by the function : moyByFactor
  ## reflev is the reference level which is used as a denominator in fold change computing
  ## collev is a vector of subscript of the different columns of levels. 
  ## ids may contains others columns but after the id+columns of average
  ## So it is possible tu use this function for all the levels in successive runs
  ## fold change by default is "log2 ratio" or "simple" numerator/denominator
  
  
  ireflev <- which(colnames(ids)==reflev)
  #collev <- collev[-which(collev==ireflev)]
  nblev <- length(collev)
  foldch <- data.frame(array(NA,dim=c(nrow(ids),nblev)),stringsAsFactors = FALSE)
  colnames(foldch) <- paste(colnames(ids)[collev],"_vs_",reflev,sep="")
  rownames(foldch) <- ids[[1]]
  for(i in 1:ncol(foldch)) 
    #if (fctype=="simple") foldch[i] <- ids[i+1]/ids[ireflev] else foldch[i] <- log2(ids[i+1]/ids[ireflev])
    if (fctype=="simple") foldch[i] <- ids[collev[i]]/ids[ireflev] else foldch[i] <- log2(ids[collev[i]]/ids[ireflev])
  
  #resids <- cbind(ids,foldch)
  resids <- foldch
  
  return(resids)
  
}

ggplotFC <- function(ids,labDEN,labNUM,sNAME,mtit,stit="log2 Fold change") {
  ## function to ggplot a fold change (labNUM/labDEN) in descending order of absolute value
  ## needs an input dataframe ids (var in rows and FC in columns) with all informations on average and FC
  ## with at least 3 columns :
  # - id in the first column
  # - Name of the variable which is used on the plot
  # - Fold change between labNUM/labDEM
  
  FC2plot <- paste(labNUM,"_vs_",labDEN,sep="")
  sFC2plot <- which(colnames(ids)==FC2plot)
  ids <- ids[,c(1,sNAME,sFC2plot)]
  
  mtitle <- paste(mtit,FC2plot, sep=" ")
  
  ## construct the name of the column to plot and then define his subscrip in log2FC
  sens <- ifelse(ids[[3]] < 0, labDEN, labNUM) 
  ids <- cbind(ids,sens)
  
  ### pour ordonner soit dans le sens decroissant de FC [[3]] soit par ordre alpha de name [[2]]
  ids <- ids[order(abs(ids[[3]])),]
  ## ids <- ids[order(ids[[2]]),]   
  
  colnames(ids)[3] <- "log2FC"
  colnames(ids)[2] <- "name"
  ids$name <- factor(ids$name, levels = ids$name) # convert to factor to retain sorted order in plot.
  
  ## ggplot of the foldchange
  p <- ggplot(ids, aes(x=name, y=log2FC, label= log2FC)) + 
    geom_bar(stat='identity', aes(fill=sens), width=.5)  +
    #geom_text(aes(label = round(pvalue24h3DEM_vs_C,4)), nudge_y = 2) +
    #scale_fill_manual(name="log2FC", labels = c(labDEN, labNUM), values = c(labDEN="green", labNUM="red")) +              
    labs(subtitle=stit,  title= mtitle) + 
    coord_flip()
  
  print(p)
  
}

meanFCggplot <- function(VM,VMsel,DM,SPL,mfact,labDEN,labNUM,sNAME,mtit) {
  
  ## compute means of levels of factor mfact, then computes all the possible FC among levels of factor mfact
  ## add these means and FC values to VM file and save it
  ## Then between the chosen levels labNUM/labDEN for a selection of variable defined by the file VMsel
  ## ggplot the fold change in descending order of the absolute values of FC
  ## sNAME is the colum number defining the label of molecules uses on Y axis
  
  subfact <- which(colnames(SPL)==mfact)
  moyDM <- moyByFactor(imtx=DM, ispl=SPL, subfact, onlyspl=TRUE, oper="mean") 
  moyVM <- merge(x = VM, y =moyDM, by.x=1, by.y=1, all.x=TRUE, all.y=FALSE)
  
  # compute log2 fold change for all level = columns of moyVM corresponding to the means of the different levels
  nlev <- ncol(x = moyDM) - 1
  collev=c((ncol(moyVM)-nlev+1):ncol(moyVM))
  lablev <- colnames(moyVM)[collev]
  
  ## for each level defined as denominator of FC, compute FC and store and write a txt file based on 
  for (l in 1:nlev) {
    FCDEN <- lablev[l]
    FC <- foldchange(moyVM,reflev=FCDEN,collev, fctype="log2")
    #subVM contains all the FC results for each level of the chosen factor mfact
    moyVM <- cbind(moyVM,FC)
    
  }
 
  ## for the mfact factor using 2 levels as Numerator and Denominator performs graphics of the fold change      
  ## prepare dataframe for ploting. If exist 
  if (!is.null(VMsel))  VMsel <- merge(x=VMsel, y=moyVM, by.x=1, by.y=1, all.x=TRUE, all.y=FALSE) else VMsel <- moyVM
  
  ggplotFC(ids=VMsel, labDEN=labDEN,labNUM=labNUM,sNAME,mtit=mtit, stit="log2 Fold change")
  
  return(VMsel)
  
}

ggplotOR <- function(ids,sNAME=2,sOR=3,mtit,stit="Odd Ratio") {
  ## function to ggplot an Odd ratio in descending order 
  ##  needs an input dataframe ids with at least 3 columns :
  # - id in the first column
  # - Name of the variable which is used on the plot (column sNAME) 
  # - odd ratio (column sOR)
  
  ## sens define the color with a threshold for OR=1
  sens <- ifelse(ids[[3]] < 1, yes="OR<1",no="OR>1") 
  ids <- cbind(ids,sens)
  ids <- ids[order(abs(ids[[3]])),]
  colnames(ids)[3] <- "OR"
  colnames(ids)[2] <- "name"
  ids$name <- factor(ids$name, levels = ids$name) # convert to factor to retain sorted order in plot.
  
  ## ggplot of the foldchange
  p <- ggplot(ids, aes(x=name, y=OR, label= OR)) + 
    geom_bar(stat='identity', aes(fill=sens), width=.5)  +
    labs(subtitle=stit,  title= mtit) + 
    coord_flip()
  
  print(p)
  
}


parse_chemFormula <- function(formul, atom2search=c("C","H","O","N","S")) {
   ## function able to count the number of atom (defined by teh vector atom2search)
   ## in a vector formul contaninig only chemical formulas
   ## return a matrix with the number of rows corresponding to formul rows and the number of columns 
   ## corresponding to the number of atom define by atom2search
   ## eliminate blank in the formula
   ## ignore case
   nbfc <- length(formul)
   nba2s <- length(atom2search)
   nbatom <- matrix(data=NA,nrow=nbfc,ncol=nba2s)
   colnames(nbatom) <- atom2search
   
   ## for each formula in vector formul....
   for (cmol in 1:nbfc) {
      # elimination of blank and transform lower to uppercase
      f1 <- toupper(gsub(" ", "", as.character(formul[cmol]), fixed = TRUE))
      ## init of vector position 
      posa2s <- rep(0,nba2s)
      ## detection of the position of the number associated to each atom define by atom2search
      for (a in 1:nba2s) {
         res <- regexpr(atom2search[a],f1, ignore.case = TRUE)
         if (res[1]>0) posa2s[a] <- res[1] + nchar(atom2search[a])
         else posa2s[a] <- 0
      }
      ## detection of number associated to each atom 1 to 99 (none <=> 1)
      for (a in 1:nba2s) {
         i <- posa2s[a]
         if (i >0) { ## i=0 means that atom is not in the formula
            nba <- c(NA,NA); j <-0
            repeat{
               nblu <- substr(f1,i,i)
               if (!is.na(as.numeric(nblu)) & i <= nchar(f1)) {
                  j<-j+1
                  nba[j] <- as.numeric(nblu)
                  i <- i+1
               } else break
            }
            if (j==0) nba[1] <- 1 ## j=0 means atom without number => value=1
            if (is.na(nba[2])) nbatom[cmol,a] <- nba[1] else nbatom[cmol,a] <- (10*nba[1]+nba[2])
         } 
         else nbatom[cmol,a] <-0
      }
   }
   return(nbatom)
}

plot_prince <- function(ids,spl,s_var,s_fact, nbcomp) {
   ## compute the % explained of each factor(spl) on the first ncomp PCA components compute on ids 
   ## ids     : data matrix with variables in rows and samples in columns : type dataframe with the First column = id of ions
   ## s_Var   : substring of quantitative variables to use in ids
   ## spl     : sample metadata with sample in row (idem colmuns of ids) and factors in columns. First column= id of samples
   ## s_fact  : substring of factors tu use in sample Metadata spl : 
   ## nb_comp : number of components tu use as "top" in prince plot 
   
   data <- as.matrix(ids[,-1]) ## on enleve la colonne 1 des identifiants variables
   data <- data[,order(colnames(data))]
   row.names(data) <- ids[[1]] ## ident variable -> colonne 1
   
   row.names(spl) <- spl[[1]]
   spl <- spl[order(row.names(spl)),]
   raw.prince <- prince(as.matrix(data),spl[,s_fact],top=nbcomp)
   prince.plot(raw.prince,Rsquared = T,key=T,cexRow=1,xlab="Principal Components (%total variation)")
   

}

molobsfreq <- function(vmdata,dm, spl,nTimeBL=10, intThreshold=0,freqBy=0,iId=1) {
   ## This function computes the observed frequency of the differents ions of a peaklist
   ## Threshold is defined as nTimeBL*average of non zero blank value (nTimeBL) or an abosulte value (nTimeBL=0 and intThreshold=absolute value of threshold)
   ## for each ion: freq <- 1-(number of samples with intensity<intThreshold / number of samples)
   ## This freq can be computed for all samples or by levels of a factor (freqBy = subscript of the factor)
   ## use the 3 W4M files vmdata : var metadata , dm : dataMatrix, spl : sample metadata
   ## iId : subscript of ID var vmdata
   
   tdm <- t(as.matrix(dm[,-1])) 
   tdm <- data.frame(row.names(tdm),tdm)
   dsf <- merge(x=spl,y=tdm,by.x = 1, by.y=1, all.x=TRUE, all.y = TRUE ) 
   
   firstvar <- ncol(spl)+1
   lastvar  <- ncol(dsf)
   
   ## selection of samples only
   splonly <- dsf[dsf$sampleType=="sample",]
   if (freqBy >0) {
      levby <- as.character(levels(as.factor(splonly[[freqBy]])) )
      nlev  <- nlevels(as.factor(splonly[[freqBy]])) 
   } else
   {
      levby <- "Freq"
      nlev=1
   }
   
   ## extract only blank to compute an average of non zero sample to define the threshold by multiplying this average by nTimeBL
   if (nTimeBL > 0){ 
      thresh <- rep(NA,lastvar - firstvar +1)
      blonly <- dsf[dsf$sampleType=="blank",]
      for(j in firstvar:lastvar) {
         inz <- which(blonly[[j]]>0)
         ## if all blank =0 then threshold is fixed to nTimeBL * intThreshold
         if (length(inz)==0) thresh[j-firstvar+1] <- nTimeBL * intThreshold
            else thresh[j-firstvar+1] <- nTimeBL*mean(blonly[inz,j]) 
      }
   }
   ## if nTimeBL==0 means the same arbitrary threshold for every ions of the data matrix
   if (nTimeBL == 0)  thresh <- rep(intThreshold,lastvar - firstvar +1)
   
   ## init result dataframe 1 line / ion nlev colonnes si by factor + 1 for mean Intensity
   nVar <- nrow(vmdata)
   resFreq <- data.frame(vmdata[,iId],matrix(data=NA, nrow = nVar, ncol = nlev),stringsAsFactors = FALSE)
   resInt  <- data.frame(vmdata[,iId],matrix(data=NA, nrow = nVar, ncol = nlev),stringsAsFactors = FALSE)
   resMin  <- data.frame(vmdata[,iId],matrix(data=NA, nrow = nVar, ncol = nlev),stringsAsFactors = FALSE)
   resMax  <- data.frame(vmdata[,iId],matrix(data=NA, nrow = nVar, ncol = nlev),stringsAsFactors = FALSE)
   
   colnames(resFreq) <- c("ions",paste(levby,"F",sep="."))
   colnames(resInt) <- c("ions",paste(levby,"Med",sep="."))
   colnames(resMin) <- c("ions",paste(levby,"min",sep="."))
   colnames(resMax) <- c("ions",paste(levby,"MAX",sep="."))

   for (i in 1:nVar) {

      ions <- dm[i,1]
      subds <- splonly[,c(1:firstvar-1,firstvar+i-1)]
      ## subscript of current variable (ions) last column of subds
      sv <- ncol(subds)
      
      ## compute number of sample with intensity > intensity Threshold thresh
      if (freqBy > 0) {
         nul1 <- rep(0,nlev)
         int1 <- rep(0,nlev)
         min1 <- rep(0,nlev)
         max1 <- rep(0,nlev)
         for (p in 1:nlev) {
            sublev <- subds[subds[[freqBy]]== levby[p],]
            ## number of sample < threshold
            nul1[p] <- nrow(sublev[sublev[[sv]]< thresh[i],])
            ## frequency : 1 - taux d'inferieur au seuil
            nul1[p] <- (1-(nul1[p]/nrow(subds[subds[[freqBy]]== levby[p],])))
            ## intensity
            int1[p] <- median(sublev[sublev[[sv]]>= thresh[i],sv])
            ## minimum or threshold (computed by blank intensities) if min in sample < threshold
            minp <- min(sublev[sublev[[sv]] >= thresh[i],sv])
            if ( minp == "Inf") min1[p] <- thresh[i]  
            else min1[p] <- min(sublev[sublev[[sv]]>= thresh[i],sv])
            max1[p] <- max(sublev[[sv]])
         } 
      }  else 
      {
         nul1 <- 1-(nrow(subds[subds[[sv]]< thresh[i],])/nrow(subds))
         int1 <- median(subds[subds[[sv]] >= thresh[i],sv])
         if (min(subds[[sv]]) < thresh[i]) min1 <- thresh[i] 
         else min1 <- min <- thresh[i]
         max1 <- max(subds[[sv]])
         
      }
      ## results 
      resFreq[i,1] <- ions
      resFreq[i,c(2:(ncol(resFreq)))] <- nul1
      resInt[i,1] <- ions
      resInt[i,c(2:(ncol(resInt)))] <- int1
      resMin[i,1] <- ions
      resMin[i,c(2:(ncol(resMin)))] <- min1
      resMax[i,1] <- ions
      resMax[i,c(2:(ncol(resMax)))] <- max1
      
      #cat(resFreq[i,1],"\n")
   }
   #resFI <- merge(x=resFreq,y=resInt,by.x=1,by.y=1,all.x=TRUE,all.y=TRUE)
   resFI <- cbind(resFreq, resInt, resMin, resMax)
   return(resFI)
} 

NA_replace <- function(DM,method=c("zero","randmin","halfmin")[2]) {
  ## replace NA values in the dataMatrix by a random value between 0 and min(ions)
  ## data matrix DM is supposed to have ions names in the first column
  ## methods of replacement : 
  ##  - "zero" : NA are replaced by 0
  ##  - "randmin" : NA are replaced by a random value between 0 and min(ions)
  ##  - "halfmin" : NA are replaced by min(ions)/2
  
  rownames(DM) <- DM[[1]]; IDvar <- colnames(DM)[1]
  tDM <- data.frame(t(DM[,-1]), stringAsFactor =FALSE)
  for (i in 1:ncol(tDM)) {
    # iNA subscript of obs with NA 
    iNA <- which(is.na(tDM[[i]]))
    if (method=="randmin") {
      # min value of the current ion is the minimum without NA values
      minval <- min(tDM[-iNA,i])
      tDM[iNA,i] <- runif(n=length(iNA),min=0, max=minval)
    } 
    else if (method=="halfmin") {
      minval <- min(tDM[-iNA,i])
      tDM[iNA,i] <- minval/2
      
    }
    else {
      ## for other method 0 replace NA
      tDM[iNA,i] <- 0
    }
  }
  DMnoNA <- t(tDM)
  DMnoNA <- data.frame(row.names(DMnoNA),DMnoNA,stringsAsFactors = FALSE); colnames(DMnoNA)[1] <- IDvar
  return(DMnoNA)
}

QCBL_filter <- function(VM,DM,spl, threshBL) {
  ## filtration of ions based on QC/BL intensity (must be > thresBL)
  ## this function aimed to filter the xtracted ions () for which QC/BL < threhBL
  
  ## First check if we have BL and QC or Sample 
  splType <- as.character(levels(as.factor(spl$sampleType)))
  spools <- which(splType == "pool")
  ssamp <- which(splType == "sample")
  sbl <- which(splType == "blank")
  itype <- which(colnames(spl)=="sampleType")
  
  ## compute mean by sampleType
  meanByType <- moyByFactor(imtx=DM,ispl = spl,subfact=itype,onlyspl=FALSE, oper="mean")
  
  ## ratio QC/BL average and set a filter 
  QCBL <- meanByType$pool/meanByType$blank
  QCBLfilter <- QCBL>threshBL
  
  ## merge ratio and filter to VM and DM
  filterVM <- data.frame(cbind(meanByType,QCBL,QCBLfilter),stringsAsFactors = FALSE)
  filterVM <- merge(x = VM, y=filterVM, by.x=1, by.y=1, all.x=TRUE, all.y=TRUE)  
  filterVM <- filterVM[filterVM$QCBLfilter==TRUE ,]

  filterDM <- data.frame(cbind(DM,QCBLfilter),stringsAsFactors = FALSE)
  filterDM <- filterDM[filterDM$QCBLfilter==TRUE ,-ncol(filterDM)]
  
  ## a tester stopifnot("error VM and DM don't match",filterVM[[1]] != filterDM[[1]])
  filterData <- list(filterVM,filterDM)
  
  return(filterData)
}

########## extraction function for both neg ou pos ionisation
## pgen : general param common to both ionisation
## pioni: parameters with possibly different values for each ionisation

xcms_xtract <- function(ioniUsed,pgen,pioni,spl,idmanip,repIoni) {
   
   resSuf=paste(idmanip,ioniUsed,sep="") 
   resfile = paste(resSuf,"txt",sep=".")
   resGroup= paste(resSuf,"2ndGrouping","pdf",sep=".")
   resIdent= paste(resSuf,"Ident",sep="_")
   cdffiles <- list.files(repIoni,recursive=TRUE,full.names=TRUE,pattern=fileform,ignore.case=TRUE)
   
   ###### save the cpu time used for each step using system.time
   tpuc <- array(data = NA,dim = c(6,3))
   tpucl <- c("readMS","findChromPeak","1stGroup","adjustRtime","2ndGroup","fillPeaks")
   
   ## data frame of sample data based on sampleMetadata load before
   samples <- sub(basename(cdffiles),pattern = paste(".",fileform,sep=""), replacement = "", fixed = TRUE)
   mfgroup <- spl[[pgen$grIndex]]  ## subscript of the column group
   pd <- data.frame(sample_name = samples, 
                    sample_group=mfgroup,
                    stringsAsFactors = FALSE)
   ## readMSdata  
   tuc <- system.time(raw_data <- readMSData(files = cdffiles, pdata = new("NAnnotatedDataFrame", pd), mode = "onDisk"))
   tpuc[1,] <-tuc[c(1:3)]
   ## find chromatogram peaks xtract from readMS files
   cwp  <- CentWaveParam(ppm = pioni$ppm,
                         peakwidth = pioni$pw, 
                         noise = pioni$noise, 
                         snthresh = pioni$sn,
                         mzdiff = pioni$mzd,
                         prefilter = c(5,pioni$noise)
   )
   tuc <-system.time(xdata <- findChromPeaks(raw_data, param = cwp))
   tpuc[2,] <-tuc[c(1:3)]
   nameGroup <- levels(as.factor(x = xdata$sample_group))
   group_colors <- brewer.pal(n=length(x = nameGroup),"Set3")
   
   ## group and rt correction
   ## 1st grouping.
   pdp1 <- PeakDensityParam(sampleGroups = xdata$sample_group,
                            minFraction = pgen$minf, bw <- pioni$bw1, binSize= pioni$mzw)
   tuc <-  system.time(xdata <- groupChromPeaks(xdata, param = pdp1) )
   tpuc[3,] <-tuc[c(1:3)]
   
   ##### F O R    T E S T
   # x11()
   # ## Define the mz slice.
   # mzr <- c(205.05, 205.1)
   # pdp0 <- PeakDensityParam(sampleGroups = xdata$sample_group,
   #                          minFraction = pgen$minf, 
   #                          bw <- 2, binSize= pioni$mzw)
   # 
   # plotChromPeakDensity(xdata, mz = mzr, param = pdp0, pch = 16, xlim = c(400, 600))
   # plotChromPeakDensity(xdata, mz = mzr, param = pdp1, pch = 16, xlim = c(400, 600))
   
   ##### E N D    T E S T
   
   
   ## retention time correction.
   pgp <- PeakGroupsParam(minFraction = pgen$miss)
   tuc <- system.time(xdata <- adjustRtime(xdata, param = pgp))
   tpuc[4,]<-tuc[c(1:3)]
   
   x11();  plotAdjustedRtime(xdata , col = group_colors)
   
   # example with obiwarp
   # xdata <- adjustRtime(xdata, param = ObiwarpParam(binSize = 0.6))
   
   ## 2nd group after adustRtime
   pdp2 <- PeakDensityParam(sampleGroups = xdata$sample_group,
                            minFraction = pgen$minf, bw <- pioni$bw2, binSize= pioni$mzw)
   tuc <- system.time(xdata <- groupChromPeaks(xdata, param = pdp2))
   tpuc[5,]<-tuc[c(1:3)]
   
   
   ## fillcchromPeaks replace NA values
   tuc <- system.time(xdata <- fillChromPeaks(xdata))
   tpuc[6,]<-tuc[c(1:3)]
   
   ## saving dataframe of result
   save(xdata, file=paste(ioniUsed,"XtractFill.Rdata",sep="")) 
   
   ## list of features
   xf <- featureDefinitions(xdata)
   
   ## CAMERA ####################################################################################
   xs <- as(xdata, 'xcmsSet')
   ## xs <- getxcmsSetObject(xdata)
   #xsa <- xsAnnotate(xs,polarity =ioniUsed)
   #Group after RT value of the xcms grouped peak
   #xsaF <- groupFWHM(xsa, perfwhm=0.6)
   #Verify grouping
   #xsaC <- groupCorr(xsaF)
   #Annotate isotopes, could be done before groupCorr
   #xsaFI <- findIsotopes(xsaC)
   #Annotate adducts
   #xsaFA <- findAdducts(xsaFI, polarity=ioniUsed)
   
   
   ###  To use a the place of the 5 last cod lines
   xsaFA <- annotate(xs, sample=NA, sigma=6, perfwhm=0.6,
                      cor_eic_th=0.75, graphMethod="hcs", pval=0.05, calcCiS=TRUE,
                      calcIso=FALSE, calcCaS=FALSE, maxcharge=3, maxiso=4, minfrac=0.5,
                      ppm=pioni$ppm, mzabs=0.015, quick=FALSE, psg_list=NULL, rules=NULL,
                      polarity=ioniUsed, multiplier=3, max_peaks=100 ,intval="maxo")
   
   
   # Peaklist post CAMERA
   xC <- (getPeaklist(xsaFA))
   ions <- paste(substr(ioniUsed,1,1),round(xC$mz,4),"T",round((xC$rt/60),1),sep="")
   xC <- data.frame(cbind(ions,xC$isotopes,xC$adduct,xC$pcgroup), stringsAsFactors = FALSE)
   colnames(xC)[2:4] <- c("isotopes","adducts","pcgroup")
   
   ## variable metadata
   ncr <- dim(xf)[2]-1
   ions <- paste(substr(ioniUsed,1,1),round(xf$mzmed,4),"T",round((xf$rtmed/60),1),sep="")
   VM <- data.frame(ions,xf[,c(1:ncr)])  
   
   ## merge VM with CAMERA annotation results
   VM <- merge(x=VM, y=xC, by.x=1, by.y=1, all.x=TRUE, all.y=TRUE)
   
   ## data matrix
   ## default missing=NA but can be : missing = "rowmin_half"
   ## values = "into" or "maxo"
   DM <- featureValues(xdata, value = "maxo")
   colnames(DM) <- samples
   DM <- data.frame(ions,DM,stringsAsFactors = FALSE)
   
   ## NA replace 
   DM <- NA_replace(DM,method="randmin")
   
   #############################  F I L T r A T I O N QC/Blank > threshBL #################################
   fildata <- QCBL_filter(VM,DM,spl,threshBL)
   
   ## write on disk
   tpuc <- data.frame(tpucl,tpuc)
   
   allres <- list(VM,DM,tpuc)
   return(allres)
}


computeCV <- function(VM,DMsimca,spl) {
   
   ## compute CV pools and samples using :
   ##  - spl metadataSample file
   ##  - VM variable metadata sample file
   ##  - DMsimca : data matrix simca with ions in column and sample in rows
   
   # subscript of samples and pools
   ispl <- which(DMsimca$sampleType=="sample")
   iQC <- which(DMsimca$sampleType=="pool")
   if (length(ispl)<2 | length(iQC)<1) {stop("pas assez de sample ou de QC")}
   
   
   firstion <- ncol(spl)+1
   lastion <- ncol(DMsimca)
   nbions <- lastion - firstion+1
   labion <- rep(NA,nrow(VM))
   CV <- matrix(0,nrow=nbions, ncol = 3); row.names(CV) <- VM[[1]]; colnames(CV) <- c("CVpools","CVsamples","CVok")
   for (j in firstion:lastion) {
      i <- j-firstion+1
      labion[i] <- colnames(DMsimca)[j]
      ## cv pools
      CV[i,1] <- sd(DMsimca[iQC,j])/mean(DMsimca[iQC,j])  
      ## cv samples
      CV[i,2] <- sd(DMsimca[ispl,j])/mean(DMsimca[ispl,j]) 
      ## if CV pool < 30% and CV pools < CV samples then CV ok =1
      if (!is.na(CV[i,1]) & !is.na(CV[i,2]) & CV[i,1] <0.3 & CV[i,1] < CV[i,2]) CV[i,3]=1 
   }
   CV <- data.frame(labion,CV, stringsAsFactors = FALSE)
   VM <- merge(VM,CV, by.x=1,by.y=1)
   
   return(VM)
}

callDIABLO <- function(Y, data, ncpUser=NA, testkeep, crossValfold = 7, crossValrep=10) {
  ### use mixOmics DIABLO to run multi-block sparsePLS
  ### limited to 2 blocks of data (example Pos and Neg mass spec. ionisation)
  ### Y : qualitative variable to be explained (1 Dimension)
  ### data : a list with the 2 numeric matrix coresponding to the 2 blocks of variables
  ### ncpuser: used to force the number of component. If=NA the algo define the best number of components.
  ### testKeep correspond to the number of variate choosen for each component
  ### crossValdfold : cross validation Mfold
  ### crossValrep   : cross validation number of iteration
  
  design = matrix(0.1, ncol = length(data), nrow = length(data), 
                  dimnames = list(names(data), names(data)))
  diag(design) = 0
  
  ## call misomics block splsda
  
  if (!is.na(ncpUser)) {nc <- ncpUser} else {nc <- (length(levels(Y))-1)}
  sgccda.res = block.splsda(X = data, Y = Y, ncomp = nc, design = design)
  
  set.seed(123) # for reproducibility, only when the `cpus' argument is not used
  # this code takes a couple of min to run
  perf.diablo = perf(sgccda.res, validation = 'Mfold', folds = crossValfold, nrepeat = crossValrep)
  
  #perf.diablo  # lists the different outputs
  # x11()
  # plot(perf.diablo)
  # dev.off()
  
  ## diablo or userdefined number of spls component
  if (!is.na(ncpUser)) ncomp <- ncpUser else ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
  cat("nb comp perf.diablo",ncomp,"\n")
  
  ## tune keepX parameter 
  tune.Zebra = tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                                 test.keepX = testkeep, design = design,
                                 validation = 'Mfold', folds = crossValfold, nrepeat = crossValrep,
                                 cpus = 2, dist = "centroids.dist")
  
  list.keepX = tune.Zebra$choice.keepX

  ## retained model with results of tune
  sgccda.res = block.splsda(X = data, Y = Y, ncomp = ncomp, keepX = list.keepX, design = design)
  print(sgccda.res)
  
  ############################### select for each block and each component, variables selected by sparse ##########################
  nbvartot <- 0
  for (i in 1:length(data)) nbvartot <- nbvartot + ncol(data[[i]])
  labvar <- colnames(data[[1]])
  for (i in 2:length(data)) labvar <- c(labvar,colnames(data[[i]]))
  ## init resu selectVar
  resSel <- data.frame(labvar,array(0,dim = c(nbvartot,ncomp)))
  for (i in 1:length(data))
     for (j in 1:ncomp) {
       varblock <- selectVar(sgccda.res, block = c(1,2), comp = j)[[i]]$name
       resSel[match(varblock,resSel[[1]]),j+1] <- 1
     }
  
  ################################################ diagnostic graphics #############################################################
  
  pdf(file="graphs_evaluation.pdf", width = 11.5, height = 7.8)
  
  for (i in 1:ncomp) plotDiablo(sgccda.res, ncomp = i)
  
  plotIndiv(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'Score plot')
  
  plotArrow(sgccda.res, ind.names = FALSE, legend = TRUE, title = 'DIABLO')
  
  # joli plot mais inutilisable
  circosPlot(sgccda.res, cutoff = 0.8, line = TRUE, color.blocks= c('darkorchid',  'lightgreen'),color.cor = c("chocolate3","grey20"), size.labels = 1.5)
  
  ## reseau de correlation entre block
  network(sgccda.res, blocks = c(1,2), color.node = c('darkorchid', 'lightgreen'), cutoff = 0.7)
  
  ## contribution des var aux composantes
  for (i in 1:ncomp) plotLoadings(sgccda.res, comp = i, contrib = 'max', method = 'median')

  ## cluster heatmap 
  cimDiablo(sgccda.res, margins = c(8, 16),size.legend = 1.2)
  
  ### assess performance of model
  set.seed(123)# for reproducibility, only when the `cpus' argument is not used
  perf.diablo = perf(sgccda.res, validation = 'Mfold', M = crossValfold, nrepeat = crossValrep, dist = 'centroids.dist')
  #perf.diablo  # lists the different outputs
  
  # Performance with Majority vote
  cat("perf Majorityvote ");cat(perf.diablo$MajorityVote.error.rate[[1]]);cat("\n")
  cat("perf Weightedvote ");cat(perf.diablo$WeightedVote.error.rate[[1]]);cat("\n")
 
  
  ## AUC roc curve
  for (i in 1:ncomp) auc.splsda = auroc(sgccda.res, roc.block = data[[i]], roc.comp = 1)
  
  ## FIN plots pdf
  dev.off()

  return(resSel)
  
}

