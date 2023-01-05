####### correct genus names 2022.12.1 #####
	library(data.table)
	library(ape)
	WFO=as.data.frame(fread("classification.txt"))#citation:WFO (2022): World Flora Online. Version 2022.12. Published on the Internet;http://www.worldfloraonline.org. Accessed on: Dec 1, 2022"

	spdis=as.data.frame(fread("Spdat_an_isrm191016.csv"))[,-1]
	genus.dis=colnames(spdis)
	dis.spname=cbind(genus.ori=genus.dis,WFO[match(genus.dis,WFO$scientificName),])
	tapply(dis.spname$genus.ori,dis.spname$taxonomicStatus,length)
	
	tree <- read.tree("processed_aV2/raxml_v2_dated_Av1_final_clean_ann.full.polyRes.tre")
	tre.spname=cbind(genus.ori=tree$tip.label,WFO[match(tree$tip.label,WFO$scientificName),])
	tapply(tre.spname$genus.ori,tre.spname$taxonomicStatus,length)
	
	sym.dis=dis.spname[!dis.spname$taxonomicStatus%in%"ACCEPTED","genus.ori"]
	sym.tre=tre.spname[!tre.spname$taxonomicStatus%in%"ACCEPTED","genus.ori"]
	unmat=sym.dis[sym.dis%in%sym.tre]
	write.csv(unmat,"genuslist.unmat.csv")
	#check names using TRNS
	#ref:https://github.com/ojalaquellueva/TNRSapi/blob/master/tnrs_api_example.R
	# Base URL for TNRS api
	url = "https://tnrsapi.xyz/tnrs_api.php"	# Production - paramo

	# Path and name of input file of taxon names 
	# Comma-delimited CSV file, first column an integer ID, second column the name
	data <- read.csv("genuslist.unmat.csv")
	# Load libraries
	library(httr)		# API requests
	library(jsonlite) # JSON coding/decoding

	# Header for api call
	headers <- list('Accept' = 'application/json', 'Content-Type' = 'application/json', 'charset' = 'UTF-8')
	# Convert the data to JSON
	data_json <- jsonlite::toJSON(unname(data))
	# Set the TNRS options
	sources <- "tropicos,wfo,wcvp,usda"	# Taxonomic sources
	class <- "tropicos"			# Family classification. Options: "tropicos", "wfo"
	mode <- "resolve"			# Processing mode
	matches <- "best"			# Return best match only
	# Convert the options to data frame and then JSON
	opts <- data.frame(c(sources),c(class), c(mode), c(matches))
	names(opts) <- c("sources", "class", "mode", "matches")
	opts_json <-  jsonlite::toJSON(opts)
	opts_json <- gsub('\\[','',opts_json)
	opts_json <- gsub('\\]','',opts_json)

	# Combine the options and data into single JSON object
	input_json <- paste0('{"opts":', opts_json, ',"data":', data_json, '}' )

	# Send the API request
	results_json <- POST(url = url,
					  add_headers('Content-Type' = 'application/json'),
					  add_headers('Accept' = 'application/json'),
					  add_headers('charset' = 'UTF-8'),
					  body = input_json,
					  encode = "json")

	# Convert JSON results to a data frame
	results_raw <- fromJSON(rawToChar(results_json$content)) 
	results <- as.data.frame(results_raw)
	write.csv(results,"genuslist.unmat.tnrs.csv")
	
	# correct distribution data
	library(data.table)	
	WFO=as.data.frame(fread("classification.txt"))
	spdis=as.data.frame(fread("Spdat_an_isrm191016.2cols.csv"))[,-1]
	sym.corr=as.data.frame(fread("genuslist.unmat.tnrs.csv"))[,-1]
	dis.spname=cbind(genus.ori=unique(spdis$sp),WFO[match(unique(spdis$sp),WFO$scientificName),])
	tree <- ape::read.tree("processed_aV2/raxml_v2_dated_Av1_final_clean_ann.full.polyRes.tre")
	tre.spname=cbind(genus.ori=tree$tip.label,WFO[match(tree$tip.label,WFO$scientificName),])
	
	acc.dis=dis.spname[dis.spname$taxonomicStatus%in%"ACCEPTED"&(dis.spname$genus.ori%in%tree$tip.label),
		c("genus.ori","taxonomicStatus","scientificName","scientificNameAuthorship","family","taxonID")]
	acc.dis$taxonID=paste("http://www.worldfloraonline.org/taxon/",acc.dis$taxonID,sep="")
	acc.dis$taxonomicStatus="Accepted"
	sym.dis=sym.corr[sym.corr$Name_submitted%in%unique(spdis$sp),
		c("Name_submitted","Taxonomic_status","Accepted_name","Accepted_name_author","Accepted_family","Accepted_name_url")]
	colnames(acc.dis)=colnames(sym.dis)
	spname.corr=rbind(sym.dis,acc.dis)
	write.csv(spname.corr,"Genus.corr20221226.csv")	
	
	spname.corr=read.csv("Genus.corr20221226.csv")[,-1]
	tapply(spname.corr$Name_submitted,spname.corr$Taxonomic_status,length)
	length(unique(spname.corr$Name_submitted))-length(unique(spname.corr$Accepted_name))
	dis.corr=cbind(spdis,spname.corr[match(spdis$sp,spname.corr$Name_submitted),c("Accepted_name","Accepted_name_author","Accepted_family","Accepted_name_url")])
	dis.corr=unique(dis.corr[,c("dis","Accepted_name","Accepted_name_author","Accepted_family","Accepted_name_url")])
	dis.corr=na.omit(dis.corr[order(dis.corr$Accepted_name),])
	write.csv(dis.corr,"SpDistribution20221208.csv")

#### distribution data matrix ##
spdat<-function(spdis){
		tree.dis.c<-unique(spdis[,c("Adcode99","Species_E1")])		
		units <- unique(tree.dis.c$Adcode99)		
		sprich.tree <- numeric(length(units))
		
		spdisc2spdisgrid <- function(units,  grid.county)
		{	
		grid.sp <- unique(as.numeric(grid.county$Adcode99))	
		#cond <- units %in% grid.sp
		cond <- grid.sp %in% units
		sprich.tree[cond] <- sprich.tree[cond]+1	
		return(as.numeric(cond))
		}
	
		
		tree.dis.grid <- tapply(X=tree.dis.c$Adcode99,INDEX=as.character(tree.dis.c$Species_E1),FUN=spdisc2spdisgrid,grid.county=tree.dis.c)
		sprich.tree[sprich.tree==0] <- NA
		spname=as.character(attributes(tree.dis.grid)$dimnames[[1]])
		
		spdat=c()
		for (i in 1:length(spname))
		{
		sp.dat=tree.dis.grid[[i]]
		spdat=cbind(spdat,sp.dat)
		}
		colnames(spdat)=c(spname)
		rownames(spdat)=c(units)

		return(spdat)		
	}
	
	library(data.table)
	spdis2=as.data.frame(fread("SpDistribution20210106.csv"))[,-1]
	spdat.glb=spdat(spdis2)	
	write.csv(spdat.glb,"Spdat_an_isrm191016.csv")

#### compare phylogeny with/without ITS ###
	library(ape)
	library(data.table)
	#no dating
	tree=read.tree("runfile/runfile/RAxML_bestTree.v2")#phy without ITS by LuoYuan20220930
	tree.ori=read.tree("runfile/raxml_v2.tree")#the ori tree not dating
	#dating
	mat=c("raxml_v2_dated_Av1_final_clean_ann.full.polyRes",#140-210 constraint
			"raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes",#140-150 constraint
			"raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes")#140-259 constraint
	n=1
	tree.ori.date <- read.tree(paste("processed_aV2/",mat[n],".tre",sep=""))
	tree.date <- read.tree(paste("raxml_noITS_dated_Av",n,".tre",sep=""))

	#remove unaccept family from TNRS https://tnrs.biendata.org/
	tpl=read.csv("Genus_synonym.csv")#from Zhiheng 2022.11.1
	#unmat.ge=tpl[is.na(tpl$Accepted_family),"Genus"]
	unmat.ge=tree.ori$tip.label[!tree.ori$tip.label%in%tpl$Genus_E1]
	#remove gym
	gym <- read.csv("gymnosperms.csv", header=F, stringsAsFactors=F)
	tree <- drop.tip(tree, c(gym[,1],unmat.ge))
	tree.ori <- drop.tip(tree.ori, c(gym[,1],unmat.ge))

	tree.date <- drop.tip(tree.date, c(gym[,1],unmat.ge))
	tree.ori.date <- drop.tip(tree.ori.date, c(gym[,1],unmat.ge))
	require(magrittr)
	require(dplyr)

	ori.tip=tree.ori$tip.label
	tip=tree$tip.label
	length(ori.tip[!ori.tip%in%tip])#1176 genus in tree.ori not in tree;
	length(tip[!tip%in%ori.tip])#0 genus in tree not in tree.ori

	ori.tip.date=tree.ori.date$tip.label
	tip.date=tree.date$tip.label
	length(ori.tip.date[!ori.tip.date%in%tip.date])#1727 genus in tree.ori not in tree;
	length(tip.date[!tip.date%in%ori.tip.date])#5 genus in tree not in tree.ori

	#compare pairwised phy dis
	source("read_Phylogenetic_Tree.R")
	di <- phylo.dist(tree.date, tip1=tree.date$tip.label,method = "branch.length")
	di.ori <- phylo.dist(tree.ori.date, tip1=tree.ori.date$tip.label,method = "branch.length")
	di.ori2=di.ori[ori.tip.date[ori.tip.date%in%tip.date],ori.tip.date[ori.tip.date%in%tip.date]]
	di.ori2[is.na(di.ori2)]=0; di.ori2=di.ori2+t(di.ori2)
	di2=di[tip.date[tip.date%in%ori.tip.date],tip.date[tip.date%in%ori.tip.date]]
	di2[is.na(di2)]=0; di2=di2+t(di2)
	di2=di2[rownames(di.ori2),colnames(di.ori2)]
	re=mantel.rtest(as.dist(di2), as.dist(di.ori2), nrepet = 99)
	di.all=list(di2, di.ori2,re)
	save(di.all,file="di.all.Rdata")

	#compare branch length
	get.age=function(tre){
		pos=which(tre$edge[,2]<=Ntip(tre))
		sp.age=tre$edge.length[pos]
		names(sp.age)=tre$tip.label[tre$edge[pos,2]]
		return(sp.age)
	}
	age.ori=get.age(tree.ori.date)
	age=get.age(tree.date)
	age.ori2=na.omit(age.ori[names(age)])
	age2=age[names(age.ori2)]

	## compare missing taxa in each fam without ITS # test monophyletic for each family
	fam.list=unique(tpl$Family);stat=c()
	for (i in fam.list){
	gen.list.ori=tpl[tpl$Family%in%i&tpl$Genus_E1%in%ori.tip,"Genus_E1"]
	gen.list.t=tpl[tpl$Family%in%i&tpl$Genus_E1%in%tip,"Genus_E1"]
	ori.ngenus=length(gen.list.ori)
	noITS.ngenus=length(gen.list.t)
	diff=ori.ngenus-noITS.ngenus
	diff.por=(ori.ngenus-noITS.ngenus)/ori.ngenus
	tmp=data.frame(Family=i,ori.ngenus,noITS.ngenus,diff,diff.por,
		ori.mono=is.monophyletic(tree.ori,gen.list.ori),noITS.mono=is.monophyletic(tree,gen.list.t))
	stat=rbind(stat,tmp)
	}
	stat=stat%>%filter(diff>=0&ori.ngenus>0)#remove gym
	stat%>%nrow()#423
	stat%>%filter(diff==0)%>%nrow()#349,82.5%
	write.csv(stat,"comphy.stat2.csv")#then modify mono mannually based on tree.ori.csv

	#plot tree in a excel
	ori.t=tree.routine(tree.ori)
	ori.t2=cbind(Family=tpl[match(rownames(ori.t),tpl$Genus_E1),"Family"],ori.t)
	write.csv(ori.t2,"tree.ori.csv")

	stat=read.csv("comphy.stat.csv")
	stat%>%filter(ori.mono==FALSE&noITS.mono==TRUE)%>%nrow() #3,1.9%
	stat%>%filter(ori.mono==FALSE)%>%nrow() #20
	stat%>%filter(ori.mono==noITS.mono)%>%nrow() #404,95.3%

	di=get(load("di.all.Rdata"))
	di=get(load("di.all.date2.Rdata"))
	di3=di[[1]][upper.tri(di[[1]])]
	di.ori3=di[[2]][upper.tri(di[[2]])]

	par(mfrow = c(2,2),mar=c(4,5.2,1,2),oma=c(1.5,1.5,0,0.15))
	hist(stat[,"diff.por"],main=NULL,xlab="Proportion of Genus \nwitout ITS in each Family",breaks=50,cex.lab = 1.2)
	mtext("a",side=3,adj=-.2,line=-2,cex=2.5,col="black")

	plot(stat[,c("ori.ngenus","noITS.ngenus")],xlab="Number of Genus in each Family",ylab="Number of Genus in \neach Family(Witout ITS)",cex.lab = 1.2);
	mtext(paste("r = ",round(cor(stat[,c("ori.ngenus","noITS.ngenus")])[1,2],2),sep=""),side=1,adj=0.8,line=-1.5,cex=1.5,col="black",font=3)
	abline(a = 0,b=1, lty = 2,lwd=2,col="red")
	mtext("b",side=3,adj=-.2,line=-2,cex=2.5,col="black")

	smoothScatter(log(age2),log(age.ori2),xlab="Tip branch length (without ITS,Ma,ln-trans)",ylab="Tip branch length (Ma,ln-trans)",cex.lab = 1.2);
	#abline(a = 0,b=1, lty = 2,lwd=2,col="red")
	mtext(paste("r = ",round(cor(log(age2),log(age.ori2)),2),sep=""),side=1,adj=0.8,line=-1.5,cex=1.5,col="black",font=3)
	mtext("c",side=3,adj=-.2,line=-2,cex=2.5,col="black")

	smoothScatter(di3,di.ori3,xlab="Pairwise phylo distance (without ITS,Ma)",ylab="Pairwise phylo distance (Ma)",cex.lab = 1.2);
	#abline(a = 0,b=1, lty = 2,lwd=2,col="red")
	mtext("Mantel r = 0.965",side=1,adj=0.8,line=-1.5,cex=1.5,col="black",font=3)
	mtext("d",side=3,adj=-.2,line=-2,cex=2.5,col="black")

#### caculate beta diversity #####
library(parallel)
## phylogenetic beta diversity function
phyBeta <- function(comm, tree, ancestors=NULL, method=c("sim", "sorrenson")) {
	
	if (length(method) > 1) method <- method[1]

	comm1 <- comm[[1]]
	comm2 <- comm[[2]]

	if (is.null(ancestors)) ancs <- ancestor(tree, tip=tree$tip.label, TRUE)
	else ancs <- ancestors

	a1 <- unique(unlist(ancs[tree$tip.label %in% comm1]))
	a2 <- unique(unlist(ancs[tree$tip.label %in% comm2]))
	shared <- a1[a1 %in% a2]
	uniq.1 <- a1[!a1 %in% a2]
	uniq.2 <- a2[!a2 %in% a1]
	
	aa <- sum(tree$edge.length[which(tree$edge[,2] %in% shared)])
	bb <- sum(tree$edge.length[which(tree$edge[,2] %in% uniq.1)])
	cc <- sum(tree$edge.length[which(tree$edge[,2] %in% uniq.2)])

	if (method == "sim") result <- 1 - aa/(aa + min(bb,cc))
	
	return(result)
	}


## distribution data
d <- read.csv("Spdat_an_isrm191016.csv", header=T)
rownames(d) <- d[,1]
d <- d[,-1]
for (i in 1:ncol(d)) d[,i] <- as.logical(d[,i])
d <- as.matrix(d)

## community pairs
comm <- list()
ii <- 1
for (i in 1:(nrow(d)-1)) {
	for (j in (i+1):nrow(d)) {
		comm1 <- colnames(d)[which(d[i,])]
		comm2 <- colnames(d)[which(d[j,])]
		comm[[ii]] <- list(comm1, comm2)
		ii <- ii + 1
		}
	print(i)
	}


## tree
## molecular trees
mat=c()
for (i in 1:3) mat=c(mat,paste("raxml_v2_dated_Av",i,"_final_clean",sep=""))
## individual trees	
for (i in 1:5) {		
	for (j in 1:3) mat=c(mat,paste("MCC_aV",j,"_",i,"_proc",sep=""),paste("tree",i,".aV",j,".dated",sep=""))
}
#"best" trees
mat=c("raxml_v2_dated_Av1_final_clean_ann.full.polyRes",#140-210 constraint
		"raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes",#140-150 constraint
		"raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes")#140-259 constraint

library(ape)	
for (n in 1:length(mat)){
	rm(tree)
	tree <- read.tree(paste("processed_aV2/",mat[n],".tre",sep=""))
	gym <- read.csv("gymnosperms.csv", header=F, stringsAsFactors=F)
	tree <- drop.tip(tree, gym[,1])

	## ancestors
	ancs <- ancestor(tree=tree, tip=tree$tip.label, include.tip = TRUE)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); date()
	## beta diversity
	#mycl <- makePSOCKcluster(names=6)
	result <- parLapply(cl=mycl, X=comm, fun=phyBeta, tree=tree, ancestors=ancs)

	pBeta.mat <- matrix(NA, nrow=nrow(d), ncol=nrow(d))
	rownames(pBeta.mat) <- colnames(pBeta.mat) <- rownames(d)
	ii <- 1
	for (i in 1:(nrow(d)-1)) {
		for (j in (i+1):nrow(d)) {
			pBeta.mat[i,j] <- result[[ii]]
			ii <- ii + 1
			}
		}
	save(pBeta.mat, file=paste("processed_aV2/psim.",mat[n],".Rdata",sep=""))
	#save(pBeta.mat, file="processed_aV2/psim.noits.Rdata")	
	date()
	stopCluster(mycl)
}

	###calculate pairwise taxonomic distances
	betasp<-function(spdis){
		x=seq(from=6,to=466,by=20)
		result=matrix(data=0,nrow=(dim(spdis)[1]-1),ncol=(dim(spdis)[1]-1))
		colnames(result)=rownames(spdis)[2:dim(spdis)[1]]
		rownames(result)=rownames(spdis)[1:(dim(spdis)[1]-1)]	
		
	   for (i in 1:(dim(spdis)[1]-1)){#
	   print(paste("i=",i))
	   for (j in (i+1) :dim(spdis)[1]){#		   
		   aa=rbind(spdis[i,],spdis[j,])
		   rownames(aa)=c(rownames(spdis)[i],rownames(spdis)[j])
		   aa[,aa[1,]+aa[2,]==0] <- NA
		   bb=na.omit(t(aa))	   
		   a=subset(bb,bb[,1]+bb[,2]==2)
		   ra=dim(a)[1]
		   b=subset(bb,bb[,1]==1&bb[,2]==0)
		   rb=dim(b)[1]
		   c=subset(bb,bb[,1]==0&bb[,2]==1)
		   rc=dim(c)[1]
		   result[i,j-1]=1-(ra/(min(rb,rc)+ra))#i与j的beta   		  
	   }
	   if (i%in%x){write.csv(result,"spbeta0.csv")}
	}		
		return(result)
	}	
	re=betasp(spdis)	
	save(re, file=paste0("processed_aV2/psim_","sp", ".Rdata"))
	
####### clustering analysis######		
	library(dendextend)
	library(recluster)
	betamatrix<-function(beta,as.dist=TRUE){#将beta矩阵结果转换为聚类要求的距离矩阵				
		add=matrix(rep(0,dim(beta)[1]),nrow=dim(beta)[1],ncol=1)
		colnames(add)=rownames(beta)[1]
		rownames(add)=rownames(beta)
		beta.t=cbind(add,beta)

		add.r=matrix(rep(0,dim(beta.t)[2]),nrow=1,ncol=dim(beta.t)[2])
		colnames(add.r)=colnames(beta.t)
		rownames(add.r)=colnames(beta.t)[length(colnames(beta.t))]
		beta.t=rbind(beta.t,add.r)

		beta.c=beta.t+t(beta.t)
		diag(beta.c) <- 0
		if (as.dist==FALSE) return (beta.c) else return(as.dist(beta.c))		
	}
	mat=c("raxml_v2_dated_Av1_final_clean_ann.full.polyRes",#140-210 constraint
		"raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes",#140-150 constraint
		"raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes",#140-259 constraint
		"sp")		
	## evaluate the best cluster methods	
	var.list=c("single","complete","median","mcquitty","average","centroid","ward.D2")
	co2=c()
	for(i in 1:length(mat)){		
		psim1=get(load(paste0("processed_aV2/psim.",mat[n], ".Rdata")))#matrix			
		if (mat%in%"sp") psim2=betamatrix(psim1) else{
			psim1[is.na(psim1)]=0
			psim.t=psim1+t(psim1)
			psim2=as.dist(psim.t)
		}
		por=do.call(cbind,lapply(1:length(var.list),function(j){
			recl <- recluster.region(psim2, method=var.list[j], mincl=4, maxcl=220)
			return(recl[[1]][,"ex.diss"])
		}))
		colnames(por)=var.list
		por2=cbind(cluster=4:10,por)
		save(por2, file=paste0("beta_mat/best_cluster",mat[i], ".Rdata"))		
		co=do.call(c,lapply(1:length(var.list),function(j){
			hc=hclust(psim2,method=var.list[j])
			return(cor(psim2,cophenetic(hc)))
		}))
		names(co)=var.list
		co2=rbind(co2,co)		
	}
	rownames(co2)=mat
	save(co2, file=paste0("beta_mat/coph_index.Rdata"))
		
	#Plot dend,nmds and map	
	library(sp)
	library(maptools)
	library(dendextend)
	library(recluster)
	library(vegan)
		
	plot.map=function(xy,nclust,mycolor,title){				
		if (title%in%"sp"){
		psim1=read.csv("spbeta_an191018.csv",header = TRUE)
		psim=psim1[,-1]		
		colnames(psim)=psim1[,1]
		psim.t=psim+t(psim)
		psim2=as.dist(psim.t)	
		}else{
			psim1=get(load(paste0("processed_aV2/psim.",title, ".Rdata")))#matrix	
			#psim1=get(load("processed_aV2/psim_Carta2021.Rdata"))#matrix	
			#psim1=get(load("processed_aV2/psim.noits.Rdata"))#matrix			
			psim1[is.na(psim1)]=0
			psim.t=psim1+t(psim1)
			psim2=as.dist(psim.t)
		}
				
		hc <- hclust(psim2, method="average")
		dend <- as.dendrogram(hc)
		
		#recl <- recluster.region(psim2, method="average", mincl=4, maxcl=8)
		#nclust=min(recl[[1]][recl[[1]][,"ex.diss"]>=0.8,"k"])
		color_b <- color_branches(dend,k=nclust,col=mycolor)		
		## map UPGMA clustering					
		leaf_col <- leaf_Colors(color_b)
		order <- order.dendrogram(color_b)
		d_leaf_col <- data.frame(order, leaf_col)
		match <- data.frame(row.names(d_leaf_col), leaf_col)
		color <- merge(xy[,1:4], match, by.x="ADCODE99", by.y="row.names.d_leaf_col.", all=F)
		
		shape<-readShapeSpatial("county_islandrm.shp")			
		shape@data=cbind(shape@data,phcol=color[match(shape@data$ADCODE99,color$ADCODE99),"leaf_col"])
		
		# run NMDS		
		nmds <- monoMDS(psim2)		
		x <- nmds$points[,1]
		y <- nmds$points[,2]
		d <- data.frame (row.names(as.matrix(psim2)),x,y)
		col_NMDS <- merge(d,color, by.x="row.names.as.matrix.psim2..", by.y="ADCODE99",sort=F)
		
		layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
		par(mar=c(4,4,1,1),oma=c(2,2,2,2))
		plot(shape, col="gray", border = "white")				
		for (i in 1:length(mycolor)) {
			shp=subset(shape,shape@data$phcol==mycolor[i])		
			plot(shp,col=mycolor[i],border=mycolor[i],add=T)
		}
		# # lon=tapply(X=shape@data$Lon,INDEX=shape@data$phname,FUN=mean)
		# # lat=tapply(X=shape@data$Lat,INDEX=shape@data$phname,FUN=mean)
		# # coor=cbind(lon,lat)
		# # text(coor, labels = rownames(coor), cex =1.2)		
		box()
		title(main=title,outer=TRUE)		
		plot(x, y, col = as.character(col_NMDS$leaf_col),pch=19,xlab="NMDS1",ylab="NMDS2",cex.axis=1)		
		abline(h=0,v=0,lty=2,col="lightgray")	
		#dendcolor=color_b;labels(dendcolor)=rep(NA,420)
		#plot(dendcolor,cex.axis=1.5)
		dendplot=as.phylo(hc)
		plot(dendplot,show.tip.label=F,edge.width=rep(2,Ntip(dendplot)),type="unrooted",use.edge.length = TRUE)
		tiplabels(text=rep("",Ntip(dendplot)),tip=seq(1,Ntip(dendplot),1),cex=0.5,frame="circle",bg=as.character(color[match(dendplot$tip.label,color$ADCODE99),"leaf_col"]))
		add.scale.bar(cex = 1.2, font = 2, col = "red")
		box()		
			
	}
	xy <- read.csv("Geo-isrm.csv",head=T)[,c("ADCODE99","Lon","Lat")]
	xy=cbind(no=1:420,xy)	
	#set color for phylo maps
	#kingdom
	mycolor.kind210259=c("#93AA00","#D39200","#DB72FB","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D")
	mycolor.kind140=c("#93AA00","#D39200","#619CFF","#F8766D","#00C19F","#00BA38","#FF61C3","#DB72FB")
	mycolor.kind.sp=c("#93AA00","#D39200","#F8766D","#619CFF","#00C19F","#00BA38","#FF61C3","#DB72FB")
	for (n in 1:length(mat)){		
		windows()
		plot.map(xy,nclust=8,mycolor.kind.sp,mat[n])
	}
	library(parallel)
## phylogenetic beta diversity function
phyBeta <- function(comm, tree, ancestors=NULL, method=c("sim", "sorrenson")) {
	
	if (length(method) > 1) method <- method[1]

	comm1 <- comm[[1]]
	comm2 <- comm[[2]]

	if (is.null(ancestors)) ancs <- ancestor(tree, tip=tree$tip.label, TRUE)
	else ancs <- ancestors

	a1 <- unique(unlist(ancs[tree$tip.label %in% comm1]))
	a2 <- unique(unlist(ancs[tree$tip.label %in% comm2]))
	shared <- a1[a1 %in% a2]
	uniq.1 <- a1[!a1 %in% a2]
	uniq.2 <- a2[!a2 %in% a1]
	
	aa <- sum(tree$edge.length[which(tree$edge[,2] %in% shared)])
	bb <- sum(tree$edge.length[which(tree$edge[,2] %in% uniq.1)])
	cc <- sum(tree$edge.length[which(tree$edge[,2] %in% uniq.2)])

	if (method == "sim") result <- 1 - aa/(aa + min(bb,cc))
	
	return(result)
	}


## distribution data
d <- read.csv("Spdat_an_isrm191016.csv", header=T)
rownames(d) <- d[,1]
d <- d[,-1]
for (i in 1:ncol(d)) d[,i] <- as.logical(d[,i])
d <- as.matrix(d)

## community pairs
comm <- list()
ii <- 1
for (i in 1:(nrow(d)-1)) {
	for (j in (i+1):nrow(d)) {
		comm1 <- colnames(d)[which(d[i,])]
		comm2 <- colnames(d)[which(d[j,])]
		comm[[ii]] <- list(comm1, comm2)
		ii <- ii + 1
		}
	print(i)
	}


## tree
## molecular trees
mat=c()
for (i in 1:3) mat=c(mat,paste("raxml_v2_dated_Av",i,"_final_clean",sep=""))
## individual trees	
for (i in 1:5) {		
	for (j in 1:3) mat=c(mat,paste("MCC_aV",j,"_",i,"_proc",sep=""),paste("tree",i,".aV",j,".dated",sep=""))
}
#"best" trees
mat=c("raxml_v2_dated_Av1_final_clean_ann.full.polyRes",#140-210 constraint
		"raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes",#140-150 constraint
		"raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes")#140-259 constraint

library(ape)	
for (n in 1:length(mat)){
	rm(tree)
	tree <- read.tree(paste("processed_aV2/",mat[n],".tre",sep=""))
	gym <- read.csv("gymnosperms.csv", header=F, stringsAsFactors=F)
	tree <- drop.tip(tree, gym[,1])

	## ancestors
	ancs <- ancestor(tree=tree, tip=tree$tip.label, include.tip = TRUE)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); date()
	## beta diversity
	#mycl <- makePSOCKcluster(names=6)
	result <- parLapply(cl=mycl, X=comm, fun=phyBeta, tree=tree, ancestors=ancs)

	pBeta.mat <- matrix(NA, nrow=nrow(d), ncol=nrow(d))
	rownames(pBeta.mat) <- colnames(pBeta.mat) <- rownames(d)
	ii <- 1
	for (i in 1:(nrow(d)-1)) {
		for (j in (i+1):nrow(d)) {
			pBeta.mat[i,j] <- result[[ii]]
			ii <- ii + 1
			}
		}
	save(pBeta.mat, file=paste("processed_aV2/psim.",mat[n],".Rdata",sep=""))
	#save(pBeta.mat, file="processed_aV2/psim.noits.Rdata")	
	date()
	stopCluster(mycl)
}

	###calculate pairwise taxonomic distances
	betasp<-function(spdis){
		x=seq(from=6,to=466,by=20)
		result=matrix(data=0,nrow=(dim(spdis)[1]-1),ncol=(dim(spdis)[1]-1))
		colnames(result)=rownames(spdis)[2:dim(spdis)[1]]
		rownames(result)=rownames(spdis)[1:(dim(spdis)[1]-1)]	
		
	   for (i in 1:(dim(spdis)[1]-1)){#
	   print(paste("i=",i))
	   for (j in (i+1) :dim(spdis)[1]){#		   
		   aa=rbind(spdis[i,],spdis[j,])
		   rownames(aa)=c(rownames(spdis)[i],rownames(spdis)[j])
		   aa[,aa[1,]+aa[2,]==0] <- NA
		   bb=na.omit(t(aa))	   
		   a=subset(bb,bb[,1]+bb[,2]==2)
		   ra=dim(a)[1]
		   b=subset(bb,bb[,1]==1&bb[,2]==0)
		   rb=dim(b)[1]
		   c=subset(bb,bb[,1]==0&bb[,2]==1)
		   rc=dim(c)[1]
		   result[i,j-1]=1-(ra/(min(rb,rc)+ra))#i与j的beta   		  
	   }
	   if (i%in%x){write.csv(result,"spbeta0.csv")}
	}		
		return(result)
	}	
	re=betasp(spdis)	
	save(re, file=paste0("processed_aV2/psim_","sp", ".Rdata"))
	
####### clustering analysis ######		
	library(dendextend)
	library(recluster)
	betamatrix<-function(beta,as.dist=TRUE){#将beta矩阵结果转换为聚类要求的距离矩阵				
		add=matrix(rep(0,dim(beta)[1]),nrow=dim(beta)[1],ncol=1)
		colnames(add)=rownames(beta)[1]
		rownames(add)=rownames(beta)
		beta.t=cbind(add,beta)

		add.r=matrix(rep(0,dim(beta.t)[2]),nrow=1,ncol=dim(beta.t)[2])
		colnames(add.r)=colnames(beta.t)
		rownames(add.r)=colnames(beta.t)[length(colnames(beta.t))]
		beta.t=rbind(beta.t,add.r)

		beta.c=beta.t+t(beta.t)
		diag(beta.c) <- 0
		if (as.dist==FALSE) return (beta.c) else return(as.dist(beta.c))		
	}
	mat=c("raxml_v2_dated_Av1_final_clean_ann.full.polyRes",#140-210 constraint
		"raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes",#140-150 constraint
		"raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes",#140-259 constraint
		"sp")		
	## evaluate the best cluster methods	
	var.list=c("single","complete","median","mcquitty","average","centroid","ward.D2")
	co2=c()
	for(i in 1:length(mat)){		
		psim1=get(load(paste0("processed_aV2/psim.",mat[n], ".Rdata")))#matrix			
		if (mat%in%"sp") psim2=betamatrix(psim1) else{
			psim1[is.na(psim1)]=0
			psim.t=psim1+t(psim1)
			psim2=as.dist(psim.t)
		}
		por=do.call(cbind,lapply(1:length(var.list),function(j){
			recl <- recluster.region(psim2, method=var.list[j], mincl=4, maxcl=220)
			return(recl[[1]][,"ex.diss"])
		}))
		colnames(por)=var.list
		por2=cbind(cluster=4:10,por)
		save(por2, file=paste0("beta_mat/best_cluster",mat[i], ".Rdata"))		
		co=do.call(c,lapply(1:length(var.list),function(j){
			hc=hclust(psim2,method=var.list[j])
			return(cor(psim2,cophenetic(hc)))
		}))
		names(co)=var.list
		co2=rbind(co2,co)		
	}
	rownames(co2)=mat
	save(co2, file=paste0("beta_mat/coph_index.Rdata"))
	
	##library(ggpubr)
	prefo=read.csv("Appendix 2-perfomace of methods.csv")#combine all the performance 	
	cpindex=na.omit(prefo[,c(37:41)])	
	method=c("single","complete","median","mcquitty","centroid","ward.D2","average","fuzzyC")
		
	perfrom.plot=function(expvar.n,cpindex.n,method,show.xlab=TRUE){
		title=colnames(cpindex.n)[2]
		colnames(cpindex.n)=c("cluster","index")
		expvar=c();exp.thes=c()
		for (i in 1:length(method)){
			tmp=expvar.n[,c(1,1+i)]
			colnames(tmp)=c("cluster","exp")		
			exp.thes.t=subset(tmp,tmp$exp>=99)[1,"cluster"]
			exp.thes=c(exp.thes,exp.thes.t)
			expvar=rbind(expvar,data.frame(method=method[i],tmp))
		}
		expvar$method=factor(expvar$method,levels=method)
		cpindex.n$cluster=factor(cpindex.n$cluster,levels=method)
		
		main=ggplot(expvar, mapping = aes(x = cluster, y = exp,color=method)) + geom_line(size=1.5)+ theme_bw()+
		theme(axis.title.x = element_blank(),
			axis.title.y = element_blank(),			
			axis.text.y  = element_text(size=20),
			legend.text = element_text(size=20),
			legend.title = element_text(size=20))+
			geom_hline(yintercept=99,col="black",size=1,linetype="longdash",alpha=0.7)+labs(color='Clustering Methods')
			#+ylab("Proportion of total PBD explained by between-cluster PBD")+xlab("Number of clusters")
		for (i in 1:length(method)) {
			if (is.na(exp.thes[i])){
				next
			}else{
				main=main+geom_vline(xintercept=exp.thes[i],col=scales::hue_pal()(8)[i],size=1,linetype="longdash")
			}
		}
		if (show.xlab) {
			main=main+theme(axis.text.x  =element_text(size=20))
		}	else {
			main=main+theme(axis.text.x  =element_blank())
		}
		
		sub <- ggplot(cpindex.n, mapping = aes(x = cluster, y = index,color=cluster)) + geom_point(size=5,show.legend=FALSE)+ theme_bw()+
		scale_color_manual(values=scales::hue_pal()(8)[1:7])+
		theme(axis.title.x = element_blank(),
			axis.title.y = element_text(size=15),
			axis.text.x  = element_text(size=12,angle =30,vjust=1, hjust = 1),
			axis.text.y  = element_text(size=12))+ylab("Co-phenetic correlation")#+ggtitle(title)
		sub$layers <- rev(sub$layers)
		p=main + annotation_custom(ggplotGrob(sub), xmin=50, ymin=5, xmax=220, ymax=65) 
		return(p)
	}
	
	p=list()
	for (i in 1:4){
		expvar.n=prefo[,c((9*(i-1)+1):(9*(i-1)+9))]
		cpindex.n=cpindex[,c(1,(i+1))]
		if (i<=2) {p[[i]]=perfrom.plot(expvar.n,cpindex.n,method,FALSE)} else {p[[i]]=perfrom.plot(expvar.n,cpindex.n,method)}
	}	
	ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],nrow=2,ncol = 2,widths=c(1,1),heights=c(1,1),
		labels=c("a","b","c","d"),label.x=0.1,label.y=1,font.label = list(size = 30),legend="right",common.legend=TRUE)

####### identify floristic realms using the best clustering method######			
	#Plot dend,nmds and map	
	library(sp)
	library(maptools)
	library(dendextend)
	library(recluster)
	library(vegan)
		
	plot.map=function(xy,nclust,mycolor,title){				
		if (title%in%"sp"){
		psim1=read.csv("spbeta_an191018.csv",header = TRUE)
		psim=psim1[,-1]		
		colnames(psim)=psim1[,1]
		psim.t=psim+t(psim)
		psim2=as.dist(psim.t)	
		}else{
			psim1=get(load(paste0("processed_aV2/psim.",title, ".Rdata")))#matrix	
			#psim1=get(load("processed_aV2/psim_Carta2021.Rdata"))#matrix	
			#psim1=get(load("processed_aV2/psim.noits.Rdata"))#matrix			
			psim1[is.na(psim1)]=0
			psim.t=psim1+t(psim1)
			psim2=as.dist(psim.t)
		}
				
		hc <- hclust(psim2, method="average")
		dend <- as.dendrogram(hc)
		
		#recl <- recluster.region(psim2, method="average", mincl=4, maxcl=8)
		#nclust=min(recl[[1]][recl[[1]][,"ex.diss"]>=0.8,"k"])
		color_b <- color_branches(dend,k=nclust,col=mycolor)		
		## map UPGMA clustering					
		leaf_col <- leaf_Colors(color_b)
		order <- order.dendrogram(color_b)
		d_leaf_col <- data.frame(order, leaf_col)
		match <- data.frame(row.names(d_leaf_col), leaf_col)
		color <- merge(xy[,1:4], match, by.x="ADCODE99", by.y="row.names.d_leaf_col.", all=F)
		
		shape<-readShapeSpatial("county_islandrm.shp")			
		shape@data=cbind(shape@data,phcol=color[match(shape@data$ADCODE99,color$ADCODE99),"leaf_col"])
		
		# run NMDS		
		nmds <- monoMDS(psim2)		
		x <- nmds$points[,1]
		y <- nmds$points[,2]
		d <- data.frame (row.names(as.matrix(psim2)),x,y)
		col_NMDS <- merge(d,color, by.x="row.names.as.matrix.psim2..", by.y="ADCODE99",sort=F)
		
		layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
		par(mar=c(4,4,1,1),oma=c(2,2,2,2))
		plot(shape, col="gray", border = "white")				
		for (i in 1:length(mycolor)) {
			shp=subset(shape,shape@data$phcol==mycolor[i])		
			plot(shp,col=mycolor[i],border=mycolor[i],add=T)
		}
		# # lon=tapply(X=shape@data$Lon,INDEX=shape@data$phname,FUN=mean)
		# # lat=tapply(X=shape@data$Lat,INDEX=shape@data$phname,FUN=mean)
		# # coor=cbind(lon,lat)
		# # text(coor, labels = rownames(coor), cex =1.2)		
		box()
		title(main=title,outer=TRUE)		
		plot(x, y, col = as.character(col_NMDS$leaf_col),pch=19,xlab="NMDS1",ylab="NMDS2",cex.axis=1)		
		abline(h=0,v=0,lty=2,col="lightgray")	
		#dendcolor=color_b;labels(dendcolor)=rep(NA,420)
		#plot(dendcolor,cex.axis=1.5)
		dendplot=as.phylo(hc)
		plot(dendplot,show.tip.label=F,edge.width=rep(2,Ntip(dendplot)),type="unrooted",use.edge.length = TRUE)
		tiplabels(text=rep("",Ntip(dendplot)),tip=seq(1,Ntip(dendplot),1),cex=0.5,frame="circle",bg=as.character(color[match(dendplot$tip.label,color$ADCODE99),"leaf_col"]))
		add.scale.bar(cex = 1.2, font = 2, col = "red")
		box()		
			
	}
	xy <- read.csv("Geo-isrm.csv",head=T)[,c("ADCODE99","Lon","Lat")]
	xy=cbind(no=1:420,xy)	
	#set color for phylo maps
	#kingdom
	mycolor.kind=list(c("#93AA00","#D39200","#DB72FB","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D"),
		c("#93AA00","#D39200","#619CFF","#F8766D","#00C19F","#00BA38","#FF61C3","#DB72FB"),
		c("#93AA00","#D39200","#DB72FB","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D"),
		c("#93AA00","#D39200","#F8766D","#619CFF","#00C19F","#00BA38","#FF61C3","#DB72FB"))
	for (n in 1:length(mat)){		
		windows()
		plot.map(xy,nclust=8,mycolor.kind[[n]],mat[n])
	}
	
##compare with carta et al. 2021###
	library(ape)
	tre=read.tree("GBOTB_TS.tre")
	splist=tre$tip.label
	mf=function(x) unlist(strsplit(x,split='_'))[1]
	a=unique(do.call(c,lapply(splist,mf)))#10475 genus (incl. gyo and ferns)
	
	library(data.table)
	spdis=as.data.frame(fread("SpDistribution20191016.csv"))[,-1]
	geoname=read.csv("Geo-isrm.csv")	
	spdis2=subset(spdis,spdis$Species_E1%in%a&spdis$Adcode99%in%geoname$ADCODE99)
	length(unique(spdis2$Species_E1))#9905 genus (angiosperms)
	spdat.glb=spdat(spdis2)	
	d=spdat.glb
	for (i in 1:ncol(d)) d[,i] <- as.logical(d[,i])
	d <- as.matrix(d)

	## community pairs
	comm <- list()
	ii <- 1
	for (i in 1:(nrow(d)-1)) {
		for (j in (i+1):nrow(d)) {
			comm1 <- colnames(d)[which(d[i,]==1)]
			comm2 <- colnames(d)[which(d[j,]==1)]
			comm[[ii]] <- list(comm1, comm2)
			ii <- ii + 1
			}
		print(i)
		}
	
	#caculate phylo beta
	tree <- read.tree("processed_aV2/raxml_v2_dated_Av1_final_clean_ann.full.polyRes.tre")#140-210 constraint
	gym <- read.csv("gymnosperms.csv", header=F, stringsAsFactors=F)
	tree <- drop.tip(tree, gym[,1])

	## ancestors
	ancs <- ancestor(tree=tree, tip=tree$tip.label, include.tip = TRUE)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); date()
	## beta diversity
	#mycl <- makePSOCKcluster(names=6)
	result <- parLapply(cl=mycl, X=comm, fun=phyBeta, tree=tree, ancestors=ancs)
	
	pBeta.mat <- matrix(NA, nrow=nrow(d), ncol=nrow(d))
	rownames(pBeta.mat) <- colnames(pBeta.mat) <- rownames(d)
	ii <- 1
	for (i in 1:(nrow(d)-1)) {
		for (j in (i+1):nrow(d)) {
			pBeta.mat[i,j] <- result[[ii]]
			ii <- ii + 1
			}
		}	
	date()
	stopCluster(mycl)
	save(pBeta.mat, file="processed_aV2/psim_Carta2021.Rdata")
	
	#Plot dend,nmds and map	
	xy <- read.csv("Geo-isrm.csv",head=T)[,c("ADCODE99","Lon","Lat")]
	xy=cbind(no=1:420,xy)
				
	psim=get(load("processed_aV2/psim.Carta2021.Rdata"))
	#psim=get(load("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata"))
	mycolor.kind.cata=c("#93AA00","#D39200","#DB72FB","#FF61C3","#00BA38","#00C19F","#619CFF","#F8766D")
	plot.map(xy,nclust=8,mycolor.kind.cata,"Carta2021")
		
#### compare The maximum clade credibility tree and the individual trees drawn from the posterior distribution ####
	library(sp)
	library(maptools)
	library(dendextend)
	library(recluster)
	library(vegan)	
	xy <- read.csv("Geo-isrm.csv",head=T)[,c("ADCODE99","Lon","Lat")]
	xy=cbind(no=1:420,xy)
	shape0<-readShapeSpatial("county_islandrm.shp")
	mycolor210=c("#93AA00","#D39200","#DB72FB","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D")			
	mycolor150=c("#93AA00","#D39200","#619CFF","#F8766D","#00C19F","#00BA38","#FF61C3","#DB72FB")	
	mycolor150.ind=c("#93AA00","#D39200","#FF61C3","#DB72FB","#F8766D","#619CFF","#00C19F","#00BA38")
	mycolor259.ind=c("#93AA00","#D39200","#DB72FB","#FF61C3","#00BA38","#00C19F","#619CFF","#F8766D")	
	#plot.map(xy,nclust=8,mycolor,colorname.kind,mat.idv[1])
	phykind=function(mat.t,xy,mycolor){				
		psim1=get(load(paste0("processed_aV2/psim.",mat.t, ".Rdata")))#matrix
		psim1[is.na(psim1)]=0
		psim.t=psim1+t(psim1)
		psim2=as.dist(psim.t)
		hc <- hclust(psim2, method="average")
		dend <- as.dendrogram(hc)
		color_b <- color_branches(dend,k=8,col=mycolor)		
		## map UPGMA clustering					
		leaf_col <- leaf_Colors(color_b)
		order <- order.dendrogram(color_b)
		d_leaf_col <- data.frame(order, leaf_col)
		match <- data.frame(row.names(d_leaf_col), leaf_col)
		color <- merge(xy[,1:4], match, by.x="ADCODE99", by.y="row.names.d_leaf_col.", all=F)
		return(color)
	}
	mat.ori=c("raxml_v2_dated_Av1_final_clean_ann.full.polyRes",#140-210 constraint
		"raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes",#140-150 constraint
		"raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes")#140-259 constraint
	mat.idv=c()
	for (j in 1:3) mat.idv=c(mat.idv,paste("MCC_aV",j,"_1_proc",sep=""))
	if (n==2){
		kind.ori=phykind(mat.ori[n],xy,mycolor150)
		kind.idv=phykind(mat.idv[n],xy,mycolor150.ind)
	} else {
		kind.ori=phykind(mat.ori[n],xy,mycolor210)
		if (n==3) kind.idv=phykind(mat.idv[n],xy,mycolor259.ind) else kind.idv=phykind(mat.idv[n],xy,mycolor210)
	}	
	kind.int=cbind(kind.ori,leaf_col.idv=kind.idv[match(kind.ori$ADCODE99,kind.idv$ADCODE99),"leaf_col"])
	
	color=c()
	for (i in 1:420){
		adc=kind.int$ADCODE99[i]
		a=kind.int$leaf_col[i]
		b=kind.int$leaf_col.idv[i]
		if (a%in%b) tmp=data.frame(kind.int[i,1:5],leaf_col2=a) else tmp=data.frame(kind.int[i,1:5],leaf_col2="gray")
		color=rbind(color,tmp)
	}			
	psim1=get(load(paste0("processed_aV2/psim.",mat.ori[n], ".Rdata")))#matrix
	psim1[is.na(psim1)]=0
	psim.t=psim1+t(psim1)
	psim2=as.dist(psim.t)
	hc <- hclust(psim2, method="average")
	dend <- as.dendrogram(hc)
	leaf_no=labels(dend)
	colset=cbind(leaf_no,color[match(leaf_no,color$ADCODE99),])		
	dend2 <- assign_values_to_leaves_edgePar(dend=dend, value = colset[,"leaf_col2"], edgePar = "col")
	labels(dend2)=rep(NA,420)		
	shape=shape0			
	shape@data=cbind(shape@data,color[match(shape@data$ADCODE99,color$ADCODE99),c("leaf_col","leaf_col2")])	
	# run NMDS		
	nmds <- monoMDS(psim2)		
	x <- nmds$points[,1]
	y <- nmds$points[,2]
	d <- data.frame (row.names(as.matrix(psim2)),x,y)
	col_NMDS <- merge(d,color, by.x="row.names.as.matrix.psim2..", by.y="ADCODE99",sort=F)
		
	par(mfrow = c(3,1),mar=c(0.5,4,0.5,1),oma=c(2,2,2,2))
	plot(shape, col="gray", border = "white")				
	mycolor=unique(color$leaf_col2)
	for (i in 1:length(mycolor)) {
		shp=subset(shape,shape@data$leaf_col2==mycolor[i])		
		plot(shp,col=mycolor[i],border=mycolor[i],add=T)
	}		
	title(main=mat.ori[n],outer=TRUE)
	box()		
	plot(dend2)	
	plot(x, y, col = as.character(col_NMDS$leaf_col2),pch=19,xlab="NMDS1",ylab="NMDS2")
	abline(h=0,v=0,lty=2,col="black")		

## caculate sil_witdth and identify uncertain region ###	
	library(cluster)
	library(data.table)
	dis <- as.data.frame(fread("Spdat_an_isrm191016.csv",head = T))
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]		
	adc=dis[,1]
	library(factoextra)	
	library(raster)
	library(ggpubr)
	data.test <- raster("extend.tree/rasters/0.00_temperature.asc",native=TRUE)
	coord=coordinates(data.test)
	rownames(coord)=paste(coord[,1],coord[,2],sep="_")
	cordacd <- read.csv("POINT2.csv")
	cordacd=subset(cordacd,cordacd$ADCODE99>0)
	rownames(cordacd)=paste(cordacd[,1],cordacd[,2],sep="_")	
	mat=c("raxml_v2_dated_Av1_final_clean_ann.full.polyRes",#140-210 constraint
		"raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes",#140-150 constraint
		"raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes")#140-259 constraint
	
	tlt=mat[1]
	psim1=get(load(paste0("processed_aV2/psim.",tlt, ".Rdata")))#matrix	
	psim1[is.na(psim1)]=0
	psim2=psim1+t(psim1)	
	set.seed(123)
	hc <- hclust(as.dist(psim2), method="average")
	plot(hc)
	dend <- rect.hclust(hc,k=8)
	names(dend)=1:8
	fun=function(i,dend) {
		dat=rep(as.numeric(names(dend)[i]),length(dend[[i]]));
		names(dat)=names(dend[[i]]);
		return(dat)
	}
	cluster=do.call(c,lapply(1:8,fun,dend))
	
	#conduct fuzzy c means
	mf=function(i,cluster,adc){
		tmp=matrix(0,420,1);rownames(tmp)=adc
		cls=cluster[cluster==i]
		tmp[names(cls),]=1
		return(tmp)
	}
	iniMem.p=do.call(cbind,lapply(1:8,mf,cluster,adc))
	cm.fanny <- fanny(psim2, k = 8, diss = TRUE, maxit = 8000,iniMem.p=iniMem.p, memb.exp = 1.05)
	plot(cm.fanny)
	
	#sil_witdth	for average clustering
	si=c()
	for (i in 1:420){
		adc.t=adc[i]
		clust.a=cluster[names(cluster)==adc.t]
		adc.a=names(cluster[cluster==clust.a])
		a.i=mean(psim2[which(rownames(psim2)==adc[i]),which(colnames(psim2)%in%adc.a[adc.a!=adc.t])])
		clust.d=unique(cluster)[unique(cluster)!=clust.a]
	
		mf=function(clust.t,adc.t,cluster,psim2){
			clust.adc=names(cluster[cluster==clust.t])
			re=mean(psim2[which(rownames(psim2)==adc.t),which(colnames(psim2)%in%clust.adc)])
			return(re)}
		dis.d=do.call(c,lapply(clust.d,mf,adc.t,cluster,psim2))
		names(dis.d)=clust.d		
		b.i=min(dis.d)
		clust.b=names(dis.d[dis.d==b.i])
		s.i=(b.i-a.i)/max(a.i,b.i)
		si.t=data.frame(adc=adc.t,clust=clust.a,nearest.clust=clust.b,si=s.i)
		si=rbind(si.t,si)
	}
	
	#plot the uncertain region in the map and the sil width
	typ=c("average","fcm")
	si2=cluster::silhouette(cluster,psim2)
	if (typ=="average"){
		si2[,1]=as.numeric(si$clust)
		si2[,2]=as.numeric(si$nearest.clust)
		si2[,3]=as.numeric(si$si)
		rownames(si2)=si$adc
		fcm.value=si[,-1]
		mycolor.kind=data.frame(region=c(0,1,2,3,4,5,6,7,8), 
			col=c("gray","#93AA00","#D39200","#DB72FB","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D"))			
	}
	if (typ=="fcm"){
		si2[,1]=cm.fanny$silinfo$widths[,1]
		si2[,2]=cm.fanny$silinfo$widths[,2]
		si2[,3]=cm.fanny$silinfo$widths[,3]		
		rownames(si2)=rownames(cm.fanny$silinfo$widths)
		fcm.value=cm.fanny$silinfo$widths	
		mycolor.kind=data.frame(region=c(0,1,2,3,4,5,6,7,8), 
			col=c("gray","#D39200","#00BA38","orange","#93AA00","#00C19F","#F8766D","#DB72FB","#FF61C3"))		
	}
	
	sil.map=function(i,thes,fcm.value,cordacd,mycolor.kind){
		fcm.value[fcm.value[,3]<=thes[i],1]=0
		cordacd2=cbind(cordacd,fcm=fcm.value[match(cordacd$ADCODE99,rownames(fcm.value)),1])	
		pC.t4=cbind(coord,cordacd2[match(rownames(coord),rownames(cordacd2)),])		
		data.test[]=pC.t4[,"fcm"]			
		test_spdf <- as(data.test, "SpatialPixelsDataFrame")
		test_df <- as.data.frame(test_spdf)
		colnames(test_df) <- c("value", "x", "y")	
		p.map=ggplot()+geom_raster(data = test_df , aes(x = x, y = y,fill = factor(value))) +
			annotate("text", x =-130, y =-30, label = paste("si >",round(thes[i],2),"\n(",names(thes)[i],")"),size=4) +
			scale_fill_manual(values =as.character(mycolor.kind[match(attributes(factor(test_df$value))$levels,mycolor.kind$region),"col"]))+		
			theme_minimal()+	
			theme(axis.title = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
			panel.grid=element_blank(),plot.title = element_blank(),legend.position="none")
		return(p.map)
	}
	thes=sort(c(0,summary(fcm.value[,3])[-c(1,6)]))
	names(thes)[names(thes)%in%""]=0	
	plotn=lapply(1:length(thes),sil.map,thes,fcm.value,cordacd,mycolor.kind)
	
	fviz=function (sil.obj, ...) 
	{
		df <- as.data.frame(sil.obj[, 1:3], stringsAsFactors = TRUE)
		df <- df[order(df$cluster, -df$sil_width), ]
		if (!is.null(rownames(df))) 
			df$name <- factor(rownames(df), levels = rownames(df))
		else df$name <- as.factor(1:nrow(df))
		df$cluster <- as.factor(df$cluster)
		mapping <- aes_string(x = "name", y = "sil_width", 
			color = "cluster", fill = "cluster")
		p <- ggplot(df, mapping) + geom_bar(stat = "identity") + 
			labs(y = "Silhouette width Si", x = "", title = paste0("Clusters silhouette plot ", 
				"\n Average silhouette width: ", round(mean(df$sil_width), 
					2))) + ggplot2::ylim(c(NA, 1))
		p <- ggpubr::ggpar(p, ...)
		p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+xlab("Geographic standard units (GSUs)")		
		p
	}
		
	p=fviz(si2,ggtheme = theme_bw()+theme(axis.title.y = element_text(size=12),			
			axis.text.y  = element_text(size=15)),
		palette=as.character(mycolor.kind[-1,"col"]),legend = "none")			
	ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],plotn[[5]],
		p,nrow=3,ncol = 2,heights=c(1,1,1),widths=c(1,1),labels=c("a","b","c","d","e","f"))
	
	ggarrange(plotn[[1]],p,nrow=2,ncol = 1,heights=c(1,1),widths=1,labels=c("a","b"))
		
	# plot si values for each realm
	sil.in.realm=function(i,fcm.value,cordacd,mycolor.kind)	{
		fcm.realm=fcm.value
		fcm.realm[fcm.realm[,1]!=i,3]=0
		cordacd2=cbind(cordacd,si=fcm.realm[match(cordacd$ADCODE99,rownames(fcm.realm)),3])	
		pC.t4=cbind(coord,cordacd2[match(rownames(coord),rownames(cordacd2)),])		
		data.test[]=pC.t4[,"si"]			
		test_spdf <- as(data.test, "SpatialPixelsDataFrame")
		test_df <- as.data.frame(test_spdf)
		colnames(test_df) <- c("value", "x", "y")	
		p.map=ggplot()+geom_raster(data = test_df , aes(x = x, y = y,fill = value)) +			
			scale_fill_gradient2(midpoint = 0,low="lightgray",mid="lightgray",high=mycolor.kind[i])+		
			theme_minimal()+	
			theme(axis.title = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
			panel.grid=element_blank(),plot.title = element_blank(),legend.position="none")
		return(p.map)
	}
	plotn=lapply(1:8,sil.in.realm,fcm.value,cordacd,mycolor.kind)
		
	ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],plotn[[5]],plotn[[6]],plotn[[7]],plotn[[8]],
		nrow=4,ncol = 2,heights=c(1,1,1,1),widths=c(1,1),labels=c("a","b","c","d","e","f","g","h"))	
		
### evaluate endemic genera within each realm #####
	library(data.table)
	library(raster)
	library(ggpubr)
	#genus.dis=as.data.frame(fread("SpDistribution20191016.csv"))[,-1]
	genus.dis=as.data.frame(fread("SpDistribution20221208.csv"))[,-1]
	
	#cluster
	psim1=get(load(paste0("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata")))#matrix	
	psim1[is.na(psim1)]=0
	psim2=psim1+t(psim1)	
	hc <- hclust(as.dist(psim2), method="average")
	plot(hc)
	dend <- rect.hclust(hc,k=8)
	names(dend)=1:8
	fun=function(i,dend) {
		dat=rep(as.numeric(names(dend)[i]),length(dend[[i]]));
		names(dat)=names(dend[[i]]);
		return(dat)
	}
	cluster=do.call(c,lapply(1:8,fun,dend))
	
	#raster file
	data.test <- raster("extend.tree/rasters/0.00_temperature.asc",native=TRUE)
	coord=coordinates(data.test)
	rownames(coord)=paste(coord[,1],coord[,2],sep="_")
	cordacd <- read.csv("POINT2.csv")
	cordacd=subset(cordacd,cordacd$ADCODE99>0)
	rownames(cordacd)=paste(cordacd[,1],cordacd[,2],sep="_")	
	
	#plot	
	plot.cluster.rich=function(adc,genus.dis,color,cordacd,data.test,midpoint=50){
		splist=unique(genus.dis[genus.dis$dis%in%adc,"Accepted_name"])
		dis=subset(genus.dis,genus.dis$Accepted_name%in%splist)
		sprich.cluster=tapply(dis$Accepted_name,dis$dis,length)
		sprich.all=tapply(genus.dis$Accepted_name,genus.dis$dis,length)
		sprich=cbind(sprich.all,sprich.cluster[match(names(sprich.all),names(sprich.cluster))])
		sprich[is.na(sprich)]=0
		fcm.value=sprich[,2]/sprich[,1]*100
		cordacd2=cbind(cordacd,fcm=fcm.value[match(cordacd$ADCODE99,names(fcm.value))])	
		pC.t4=cbind(coord,cordacd2[match(rownames(coord),rownames(cordacd2)),])		
		data.test[]=pC.t4[,"fcm"]			
		test_spdf <- as(data.test, "SpatialPixelsDataFrame")
		test_df <- as.data.frame(test_spdf)
		colnames(test_df) <- c("value", "x", "y")	
		p.map=ggplot()+geom_raster(data = test_df , aes(x = x, y = y,fill = value)) + 		
				scale_fill_gradient2(midpoint = midpoint,low="lightgray",mid="lightgray",high=color,name="Proportion \nof Genera")+		
				theme_minimal()+	
				theme(axis.title = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
				panel.grid=element_blank(),plot.title = element_blank(),legend.position="right")
		return(p.map)
	}
	p=list()
	mycolor.kind=c("#93AA00","#D39200","#DB72FB","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D")		
	for (i in 1:8){	
		adc=names(cluster[cluster==i])
		color=mycolor.kind[i]
		#pdf(paste0("processed_aV2/cluster.rich_", i, ".pdf",sep=""), width = 10, height = 6)
		 p[[i]]=plot.cluster.rich(adc,genus.dis,color,cordacd,data.test,30)		
		# #dev.off()
	}
	ggarrange(p[[1]],p[[2]],p[[7]],p[[8]],p[[6]],p[[5]],p[[3]],p[[4]],
	nrow=4,ncol =2,widths=c(1,1,1,1),heights=c(1,1),labels=c("a","b","c","d","e","f","g","h"),label.x=0.03,font.label = list(size = 20))

	
	require(magrittr)
	require(dplyr)
	real.name=data.frame(cluster=c(1,2,7,8,6,5,3,4),realm=c("Sah-Ara","Hol","Pata","NeoTro","Indo-Mal","Afr","Aus","NeZ"))
	dis.realm=cbind(genus.dis,cluster=cluster[match(genus.dis$dis,names(cluster))])%>%
		left_join(real.name,by="cluster")
	genus.stat=dis.realm %>%  group_by(Accepted_name,realm)%>% summarise(n = n()) %>% mutate(freq = n / sum(n))
	max.freq=genus.stat %>%group_by(Accepted_name)%>% summarise(max = max(freq))#min value for the max freq of all genera is 0.23
	summary(genus.stat$freq)
	genus.stat2=genus.stat %>% filter(freq>0.1)
	genus.stat2$dis=rep(0,nrow(genus.stat2))
	genus.stat2[genus.stat2$freq>0.7,"dis"]=1
	bit=function(realm.ls){
		n=length(realm.ls)
		set=0:(2^n-1)
		rst=matrix(0,ncol=n,nrow=2^n)
		for(i in 1:n){
			rst[,i]=ifelse((set-rowSums(rst*rep(c(2^((n-1):0)),each=2^n)))/(2^(n-i))>=1,1,0)
		}
		colnames(rst)=realm.ls
		rst=as.data.frame(rst)
		rst$num.realm=apply(rst,1,sum)
		rst=rst[rst$num.realm>0,]
		rst[order(rst$num.realm),-(n+1)]
	}
	realm.ls=c("Sah-Ara","Hol","Pata","NeoTro","Indo-Mal","Afr","Aus","NeZ")
	selt=bit(realm.ls)
	comb=apply(selt,1,function(x){realm.ls[which(x==1)]})
	comb=lapply(comb,sort)
	stat=vector("list",length=length(comb));names(stat)=lapply(comb,paste,collapse="+")	
	gen.ls=unique(genus.stat2$Accepted_name)
	for(ge in gen.ls){
		tmp=genus.stat2%>% filter(Accepted_name==ge)
		pos=which(names(stat)%in%paste(sort(as.character(tmp$realm)),collapse="+"))
		stat[[pos]]=c(stat[[pos]],ge)
	}
	stat.num=do.call(rbind,lapply(stat,length))%>%as.data.frame()
	colnames(stat.num)="Num.of.genera"
	stat.num$Num.of.realm=do.call(c,lapply(comb,length))
	stat.num$occupied.realm=rownames(stat.num)
	stat.num=stat.num[,c(2,3,1)]
	stat.num=stat.num%>%filter(Num.of.genera>0)
	write.csv(stat.num,"genusstat.csv")
			
### clades in the regions which belong to different kingdoms based on PBD and TBD ####
	phname=read.csv("cluster.average.merge.csv")[,c("adcode_ave2","sp2","phname")]#names of each phylo kingdom
	dis=read.csv("Spdat_an_isrm191016.csv")	
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]	
	taxonmic=read.csv("Gyo-angiosperm.csv")
	apg=read.csv("APGIV-taxonmic.csv")
	taxon=cbind(taxonmic,APG=apg[match(taxonmic$Family,apg$Family_E),"clades"])	
	
	neo=phname[phname[,"phname"]%in%"Neotropical"&phname[,"sp2"]==4,]	
	sp.neo=colSums(subset(spdis,rownames(spdis)%in%neo$adcode_ave2))
	sp.neo=sort(sp.neo[sp.neo>0], decreasing = TRUE)/dim(neo)[1]*100 #the proportion of species range sizes to the whole region
	clades.neo=data.frame(taxon[match(names(sp.neo),taxon$Genus),-1],Genus=names(sp.neo),range.por=sp.neo)
	
	hol=phname[phname[,"phname"]%in%"Holarctic"&phname[,"sp2"]==6,]
	sp.hol=colSums(subset(spdis,rownames(spdis)%in%hol$adcode_ave2))
	sp.hol=sort(sp.hol[sp.hol>0], decreasing = TRUE)/dim(hol)[1]*100 #
	clades.hol=data.frame(taxon[match(names(sp.hol),taxon$Genus),-1],Genus=names(sp.hol),range.por=sp.hol)	
	#clades.neo2=subset(clades.neo,clades.neo$range.por>50)#clades covering more than half of the whole area.
	apg.neo=sort(tapply(clades.neo$Genus,clades.neo$Family,length), decreasing = TRUE)	
	#clades.hol2=subset(clades.hol,clades.hol$range.por>50)#clades covering more than half of the whole area.
	apg.hol=sort(tapply(clades.hol$Genus,clades.hol$Family,length), decreasing = TRUE)
	
	clades.diff=rbind(cbind(region="Hol-IndMal",clades.hol),cbind(region="NT-ChlPat",clades.neo))
	write.csv(clades.diff,"clades.diff.TBD.PBD.csv")
	
###  floristic realm through time ####
	rm(list = ls())
	###ref Evolutionary history of zoogeographical regions surrounding the Tibetan Plateau
	## R function for the calculation of phylogenetic beta diversity at different phylogenetic depth was obtained from Daru et al., 2018.
	# collapsing nodes at a given age in the tree (use mammal tree as an example)
	library(ape)
			
	tre_ori=read.tree("raxml_v2_dated_Av1_final_clean_ann.full.polyRes.tre")#全球属级的树,140-210 constraint
	#tre_ori=read.tree("raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes.tre")#全球属级的树,140-150 constraint
	#tre_ori=read.tree("raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes.tre")#全球属级的树,140-259 constraint
	gym <- read.csv("gymnosperms.csv", header=F, stringsAsFactors=F)
	tre <- drop.tip(tre_ori, gym[,1])
		
	phy.dep=c(5,10,15,20,30,40,50,60,80,100,120,140,160,180)	
	get.tree=function(dep,tre){
		require(ape);require(phytools)
		new.tree <- tre	
		node.age<- dep   ###specify here threshold age for collapsing nodes in the tree	
		goal.length<-max(nodeHeights(tre))-node.age		
		while (round(max(nodeHeights(new.tree)), 4)>round(goal.length,4)){
		  k<-nodeHeights(new.tree)
		  w<-which(k[,2]==max(k[,2]))[1]
		  if(k[w,2]>goal.length){
			z<-k[w,2]-goal.length
			new.tree$edge.length[w]<-new.tree$edge.length[w]-z
			if(new.tree$edge.length[w]<0){
			  new.tree$edge.length[w]<-0
			}#end if
		  }#end if
		}#end while
		return(new.tree)
	}
	
	library(parallel)	
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); date()	
	tree.dep <- parLapply(cl=mycl,phy.dep,get.tree,tre)	
	names(tree.dep)=paste(phy.dep,"Ma")
	save(tree.dep, file="His_trees.Rdata")
	date()	
	stopCluster(mycl)
	
	## R function for calculating pairwise phylobetasim distances
	tree.dep=get(load("His_trees.Rdata"))	
	for (n in 1:length(tree.dep)){
		rm(tree)
		tree <- tree.dep[[n]]	
		dep=names(tree.dep)[n]
		## ancestors
		ancs <- ancestor(tree=tree, tip=tree$tip.label, include.tip = TRUE)
		no_cores <- detectCores() - 1
		mycl <- makePSOCKcluster(no_cores); date()
		## beta diversity
		#mycl <- makePSOCKcluster(names=6)
		result <- parLapply(cl=mycl, X=comm, fun=phyBeta, tree=tree, ancestors=ancs)

		pBeta.mat <- matrix(NA, nrow=nrow(d), ncol=nrow(d))
		rownames(pBeta.mat) <- colnames(pBeta.mat) <- rownames(d)
		ii <- 1
		for (i in 1:(nrow(d)-1)) {
			for (j in (i+1):nrow(d)) {
				pBeta.mat[i,j] <- result[[ii]]
				ii <- ii + 1
				}
			}
		save(pBeta.mat, file=paste("processed_aV2/psim.",dep,"raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata",sep=""))
		date()
		stopCluster(mycl)
	}

	## Figure 2
		# run UPGMA clustering
		# dis=read.csv("Spdat_an_isrm191016.csv")
		# spdis=as.matrix(dis[,2:dim(dis)[2]])
		# rownames(spdis)=dis[,1]	
		# subres <- read.csv("phylobetasimfin_an150.csv")
		# psim <- betamatrix(subres,spdis)
		# save(psim, file=paste0("extend.tree/psim_", 0, "Ma.Rdata"))	
	library(dendextend)
	library(recluster)
	library(vegan)
	library(raster)	
	xy <- read.csv("Geo-isrm.csv",head=T)[,c("ADCODE99","Lon","Lat")]
	xy=cbind(no=1:420,xy)
	phy.dep=c(0,5,10,15,20,30,40,50,60,80,100,120,140,160,180)	
	resList <- readRDS("resList.rds")	
	cordacd <- read.csv("POINT2.csv")
	cordacd=subset(cordacd,cordacd$ADCODE99>0)
	rownames(cordacd)=paste(cordacd[,1],cordacd[,2],sep="_")	
	sublist=resList[match(rownames(cordacd),names(resList))]
	mycolor=vector("list",length(phy.dep));names(mycolor)=paste(phy.dep,"Ma",sep="")
	mycolor[["0Ma"]]=c("#93AA00","#D39200","#DB72FB","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D")	
	mycolor[["5Ma"]]=c("#93AA00","#D39200","#DB72FB","#FF61C3","lightgreen","#619CFF","#F8766D")		
	mycolor[["10Ma"]]=c("#93AA00","#D39200","#DB72FB","#FF61C3","lightgreen","#619CFF","#F8766D")
	mycolor[["15Ma"]]=c("#93AA00","#D39200","#DB72FB","#FF61C3","lightgreen","#619CFF","#F8766D")
	mycolor[["20Ma"]]=c("#D39200","#93AA00","#DB72FB","#FF61C3","#00BA38","#00C19F","#619CFF","#F8766D")	
	mycolor[["30Ma"]]=c("#D39200","#93AA00","#DB72FB","#FF61C3","#00BA38","#00C19F","#619CFF","#F8766D")	
	mycolor[["40Ma"]]=c("#93AA00","#D39200","#DB72FB","#FF61C3","lightgreen","#619CFF","#F8766D")	
	mycolor[["50Ma"]]=c("#93AA00","#D39200","#DB72FB","#FF61C3","lightgreen","#619CFF","#F8766D")	
	mycolor[["60Ma"]]=c("#93AA00","#D39200","#DB72FB","#FF61C3","purple","#619CFF","#F8766D","#00BA38","#00C19F")	
	mycolor[["80Ma"]]=c("red","#D39200","#93AA00","#619CFF","#F8766D","#00C19F","#00BA38","#DB72FB","gray","#FF61C3")	
	mycolor[["100Ma"]]=c("#D39200","#93AA00","lightgreen","#619CFF","#DB72FB","#FF61C3","black","purple","#00BA38","pink")	
	mycolor[["120Ma"]]=c("black","#00BA38","#619CFF","pink","purple","#DB72FB","#FF61C3","darkgray","#D39200","#93AA00")	
	mycolor[["140Ma"]]=c("black","darkgray","#93AA00","#619CFF","pink","purple","#D39200")
	mycolor[["160Ma"]]=c("#00BA38","green","brown","darkgreen","#93AA00","#00C19F","pink","gray","purple","#DB72FB","#FF61C3")
	mycolor[["180Ma"]]=c("red","#F8766D","pink","#D39200","#00BA38","#93AA00")
	#scales::show_col(hue_pal()(18))
	#scales::show_col(c("#93AA00","#D39200","#DB72FB","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D"))
	for (j in (1:length(phy.dep))){
		k <- phy.dep[j]
		if (k==0) psim1=get(load(paste("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata",sep=""))) else {
			psim1=get(load(paste("processed_aV2/psim.",k," Maraxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata",sep="")))
		}			
		psim1[is.na(psim1)]=0
		psim=as.dist(psim1+t(psim1))		
		psim <- sort_dist_mat(psim, by_rows = TRUE, by_cols = TRUE)	
		# evaluate the best cluster		
		# recl <- recluster.region(psim, method="average", mincl=2, maxcl=40)
		# nclust=min(recl[[1]][recl[[1]][,"ex.diss"]>=0.8,"k"])#解释力大于0.8
		nclust=length(mycolor[[j]])
		# color UPGMA clustering
		col <- mycolor[[j]]
		#col <- rainbow(nclust)
		
		hc <- hclust(psim, method="average")
		dend <- as.dendrogram(hc)		
		color_b <- color_branches(dend,k=nclust,col=col)	
		
		# map UPGMA clustering					
		leaf_col <- leaf_Colors(color_b)
		order <- order.dendrogram(color_b)
		d_leaf_col <- data.frame(order, leaf_col)
		match <- data.frame(row.names(d_leaf_col), leaf_col)
		color0 <- merge(xy, match, by.x="ADCODE99", by.y="row.names.d_leaf_col.", all=F)
		phname.id=data.frame(leaf_col=col,phname=seq(1:nclust))
		color=cbind(color0,phname.id=phname.id[match(color0$leaf_col,phname.id$leaf_col),"phname"])
		# run NMDS	
		set.seed(1012)
		nmds <- monoMDS(psim)		
		x <- nmds$points[,1]
		y <- nmds$points[,2]
		d <- data.frame (row.names(as.matrix(psim)),x,y)
		col_NMDS <- merge(d,color, by.x="row.names.as.matrix.psim..", by.y="ADCODE99",sort=F)
		
		# draw plot		
		data.test <- raster::raster(paste0("extend.tree/rasters/",k,".00_temperature.asc",sep=""),native=TRUE)
		coord=coordinates(data.test)
		rownames(coord)=paste(coord[,1],coord[,2],sep="_")					
		pC.t=c()
		for (i in 1:length(sublist)){
			temp=na.omit(sublist[[i]][sublist[[i]]$age==k,c(3,2)])
			pC.t=rbind(pC.t,temp)
		}
		rownames(pC.t)=names(sublist)				
		aa=floor(as.matrix(pC.t)/0.5)*0.5
		aa[which(floor(aa[,1])-aa[,1]==0),1]=aa[which(floor(aa[,1])-aa[,1]==0),1]+0.5
		aa[which(floor(aa[,2])-aa[,2]==0),2]=aa[which(floor(aa[,2])-aa[,2]==0),2]+0.5		
		cordacd2=cbind(cordacd,color[match(cordacd$ADCODE99,color$ADCODE99),c("leaf_col","phname.id")])
		pC.t2=data.frame(adcode=cordacd2$ADCODE99,cordacd2[,c("leaf_col","phname.id")],aa[match(rownames(cordacd2),rownames(aa)),])	
		pC.t3=data.frame(ID=paste(pC.t2[,"lon"],pC.t2[,"lat"],sep="_"),pC.t2)		
		pC.t4=cbind(coord,pC.t3[match(rownames(coord),pC.t3$ID),])
		data.test[]=pC.t4$phname
		windows()
		pdf(paste0("processed_aV2/psim210_", k, "Ma.pdf",sep=""), width = 8, height = 10)
		par(mfrow = c(3,1),mar=c(0.5,4,0.5,1),oma=c(2,2,2,2))		
		raster::plot(data.test,col=col,axes=FALSE,legend=FALSE)
		title(main=paste(k,"Ma"),outer=TRUE)
		dendcolor=color_b;labels(dendcolor)=rep(NA,420)
		plot(dendcolor,cex.axis=2.5)	
		plot(x, y, col = as.character(col_NMDS$leaf_col),pch=19,xlab="NMDS1",ylab="NMDS2",cex.axis=2.5)
		dev.off()
	}
	
	### Figure 3
	library(vegan)
	library(dendextend)
	phy.dep=c(0,5,10,15,20,30,40,50,60,80,100,120,140,160,180)	
	nmds=lapply(1:length(phy.dep),function(i){
		k <- phy.dep[i]
		if (k==0) psim1=get(load(paste("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata",sep=""))) else {
			psim1=get(load(paste("processed_aV2/psim.",k," Maraxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata",sep="")))
		}							
		psim1[is.na(psim1)]=0
		psim=as.dist(psim1+t(psim1))		
		psim <- sort_dist_mat(psim, by_rows = TRUE, by_cols = TRUE)
		set.seed(1012)		
		nmds <- monoMDS(psim)
		return(nmds)
	})	
	vare.proc <-lapply(2:length(phy.dep),function(i){
		set.seed(10)
		vare.proc <- procrustes(nmds[[1]], nmds[[i]])
		if(i==2){tmp=rbind(vare.proc$X,vare.proc$Yrot)} else{
		tmp=vare.proc$Yrot
		}	
		return(tmp)
	})
	vare.proc2=do.call(rbind,vare.proc)
	time <- rep(1:length(phy.dep), each=420)	
	dataf <- data.frame(time,ADCODE99=rownames(vare.proc2),vare.proc2)

	## plot Procrustes result
	#color <- read.csv("02_NMDS_time_col.csv",head=T)
	xy <- read.csv("Geo-isrm.csv",head=T)[,c("ADCODE99","Lon","Lat")]
	psim1=get(load(paste("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata",sep="")))					
	psim1[is.na(psim1)]=0
	psim=as.dist(psim1+t(psim1))		
	psim <- sort_dist_mat(psim, by_rows = TRUE, by_cols = TRUE)	
	col <- mycolor[[1]]
	hc <- hclust(psim, method="average")
	dend <- as.dendrogram(hc)	
	color_b <- color_branches(dend,k=8,col=col)	
	leaf_col <- leaf_Colors(color_b)
	order <- order.dendrogram(color_b)
	d_leaf_col <- data.frame(order, leaf_col)
	match0 <- data.frame(row.names(d_leaf_col), leaf_col)
	color0 <- merge(xy, match0, by.x="ADCODE99", by.y="row.names.d_leaf_col.", all=F)
	color=cbind(color0,region=color0$leaf_col)	
	match <- merge(dataf,color,by="ADCODE99")
	
	plot(NA, NA, xlim=c(min(match$MDS1),max(match$MDS1)), 
	ylim=c(min(match$MDS2),max(match$MDS2)), xlab="NMDS1", ylab="NMDS2")
	abline(h = 0, v = 0, col = "gray60",lty = 3)
	adc=as.character(unique(match$ADCODE99))
	for (i in 1:length(adc)) {
		a <- match[which(match$ADCODE99==adc[i]),c("MDS1","MDS2")]
		col1 <- match[which(match$ADCODE99==adc[i]),]
		col2 <- as.character(col1$region)
		col3 <- adjustcolor(col2, alpha.f = 0.1)
		lines(a,col=col3)
	}
	
	# average gridded species assemblages
	region.seq=c("#FF61C3","#619CFF","#DB72FB","#93AA00","#F8766D","#00C19F","#00BA38","#D39200")
	region <- aggregate(cbind(MDS1,MDS2) ~ region+time, data = match, mean)
	for (i in 1:8) {
		a <- region[which(region$region==region.seq[i]),]		
		line <- data.frame(a$MDS1,a$MDS2)
		lines(line,col="black", lwd=8)
		lines(line,col=adjustcolor(as.character(a$region),alpha.f = 0.9), lwd=4)		
		#arrows <- data.frame(a$MDS1,a$MDS2)[1:2,]
		#arrows(arrows[2,1],arrows[2,2],arrows[1,1],arrows[1,2],col="black", lwd=8,length=0.2,angle=20)
		#arrows(arrows[2,1],arrows[2,2],arrows[1,1],arrows[1,2],col=adjustcolor(as.character(a$region),alpha.f = 0.9), lwd=4,length=0.2,angle=20)
	}
	# add points to each phylogenetic depth
	for (i in 1:length(phy.dep)) {
		b <- region[which(region$time==i),]
		point <- data.frame(b$MDS1,b$MDS2)			
		points(point, col=adjustcolor("white",alpha.f = 0.6), pch=20,cex=1.2)
		points(point, col=adjustcolor(as.character(b$region),alpha.f = 0.6), pch=20,cex=1)		
	}	
	### end	
	
	
####### Node-based analysis of species distributions ######
	library(nodiv)
	library(ape)
	library(sp)
	library(maptools)

	dis=read.csv("Spdat_an_isrm191016.csv")	
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]		
	tre_ori=read.tree("raxml_v2_dated_Av1_final_clean_ann.full.polyRes.tre")#,140-210 constraint
	#tre_ori=read.tree("raxml_v2_dated_Av1_final_clean_b_ann.full.polyRes.tre")#140-150 constraint
	#tre_ori=read.tree("raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes.tre")#140-259 constraint
	d=subset(tre_ori$tip.label,!tre_ori$tip.label%in%colnames(spdis))
	tre <- drop.tip(tre_ori, tip=d)
	regionlist=read.csv("cluster.average.merge.csv")
	projection <- CRS ("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
	shape<-readShapeSpatial("county_islandrm.shp")	
	proj4string(shape) <-"+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs"	
	coord=shape@data[,c(1,3,4)]	
	coords=na.omit(unique(coord[match(rownames(spdis),coord$ADCODE99),]))	
	#plot(shape, axes=TRUE, border="gray")
	region=colnames(regionlist)[8:14]
	
	res=vector("list",length(region))#node-based result
	gnd=vector("list",length(region))#GND values of all nodes
	sos.re=vector("list",length(region))#SOS values for nodes with GND>0.65
	ccoquettes=vector("list",length(region))#data prepared for node based analysis
	for (i in 1:length(region)){
		nodelist1=subset(regionlist,regionlist[[region[i]]]>0)		
		spdis.t0=spdis[match(nodelist1$adcode,as.integer(rownames(spdis))),]
		spdis.t=spdis.t0[,colSums(spdis.t0)>0]					
		d.t=subset(tre$tip.label,!tre$tip.label%in%colnames(spdis.t))
		tre.t <- drop.tip(tre, tip=d.t)	
	
		coords.t=coords[match(rownames(spdis.t),coords$ADCODE99),]
		ccoquettes[[i]]=nodiv_data(tre.t, spdis.t, coords.t, proj4string_in = projection,type="auto",shape = shape)	
		res[[i]] <- Node_analysis(ccoquettes[[i]], 500,"rdtable")#"quasiswap")	
		gnd[[i]]=na.omit(data.frame(node=summary(res[[i]])$nodes,GND=summary(res[[i]])$GND))
		
		gnd[[i]]=gnd[[i]][order(gnd[[i]]$GND,decreasing = TRUE),]
		gfig=as.numeric(as.vector(gnd[[i]][gnd[[i]]$GND>0.65,"node"]))	
		sos=c()
		for (j in 1:length(gfig)) {	
			sos=rbind(sos,SOS(res[[i]],gfig[j]))
		}	
		rownames(sos)=gfig
		sos.t=t(sos)
		sos.re[[i]]=rbind(gnd=gnd[[i]][match(colnames(sos.t),gnd[[i]]$node),"GND"],sos.t)
		
	}
	saveRDS(res,"res.rds")
	saveRDS(gnd,"gnd.rds")
	saveRDS(sos.re,"sos.re.rds")
	saveRDS(ccoquettes,"ccoquettes.rds")
	
## get table of node edge, node label and GND value 
	extractNodeAge <- function(tree, node=NULL) {
		root.dist <- function(tree, method = c("n.node","branch.length"), type=c("tip","node","both"))
		{
		if (length(method) > 1) method <- method[1]
		if (length(type) > 1) type <- type[1]
		
		if (class(tree) != "phylo") stop("The tree is not a phylo tree")
		
		## find the root
		root.label <- unique(tree$edge[,1][!tree$edge[,1] %in% tree$edge[,2]])
		
		tree.edge <- rbind(tree$edge,c(0,root.label))
		ii <- order(tree.edge[,2])
		tree.edge <- tree.edge[ii,]
		
		## the results of root distance
		N.tip <- Ntip(tree)
		N.node <- tree$Nnode
		N.edge <- Nedge(tree)
		if (type == "tip") {
			N.ii <- 1:N.tip
			rd <- numeric(N.tip)
			names(rd) <- tree$tip.label
			}
		else if (type =="node") {
			N.ii <- N.tip + (1:N.node)
			rd <- numeric(N.node)
			names(rd) <- tree$node.label
			}
		else if (type == "both") {
			N.ii <- 1:(N.tip + N.node)
			rd <- numeric(N.tip + N.node)
			names(rd) <- c(tree$tip.label, tree$node.label)
			}
		
		## count the number of nodes from root to tips
		if (method == "n.node")
			{
			for (i in 1:length(N.ii))
				{
				tip <- descendant <- N.ii[i]
				for (j in 1:N.node)
					{
					if (tree.edge[descendant,1] == 0)
						{
						rd[i] <- 0
						break
						}
					if (tree.edge[descendant,1] == root.label) 
						{
						rd[i] <- j
						break
						}
					ancestor <- tree.edge[descendant,1]
					descendant <- tree.edge[ancestor,2]
					}
				}
			}
		else if (method == "branch.length")
			{
			if (is.null(tree$edge.length)) edge.length <- rep(1, times = N.tip + N.node)
			if (!is.null(tree$edge.length))	edge.length <- c(tree$edge.length, 0)[ii]
			
			for (i in 1:length(N.ii))
				{
				tip <- descendant <- N.ii[i]
				for (j in 1:N.node)
					{
					rd[i] <- rd[i] + edge.length[descendant]
					if (tree.edge[descendant,1] == root.label | tree.edge[descendant,1] == 0) 
						{
						#rd[i] <- rd[i] + edge.length[root.label]
						break
						}
					ancestor <- tree.edge[descendant,1]
					descendant <- tree.edge[ancestor,2]
					}
				}		
			}
		return(rd)
		}

		  rd <- root.dist(tree, method = "branch.length", type="both")
		  depth <- max(rd, na.rm=TRUE)     

		  ## find the root
		  root.label <- unique(tree$edge[,1][!tree$edge[,1] %in% tree$edge[,2]])
		  tree.edge <- rbind(tree$edge,c(0,root.label))
		  ii <- order(tree.edge[,2])
		  tree.edge <- tree.edge[ii,]
		  edge.length <- c(tree$edge.length, 0)[ii]  

		  age.crown <- as.numeric(c(depth - rd)[Ntip(tree) + c(1:Nnode(tree))])
		  age.stem <- as.numeric(c(depth - rd + edge.length)[Ntip(tree) + c(1:Nnode(tree))])
		  if (is.null(node)) {
			   result <- data.frame(Node=1:Nnode(tree), crown=age.crown, stem=age.stem)
			   if (!is.null(tree$node.label)) {
					pos <- which(!is.na(tree$node.label) & tree$node.label != "")
					result[pos,1] <- tree$node.label[pos]
					}
		   }
		  else {
			   if (is.character(node)) node.pos <- which(tree$node.label %in% node)
			   else if (is.numeric(node)) node.pos <- which(tree.edge[,2] %in% c(node + Ntip(tree)))
			   result <- data.frame(Node=c(1:Nnode(tree))[tree.edge[node.pos,2]-Ntip(tree)],
									crown=age.crown[tree.edge[node.pos,2]-Ntip(tree)],
									stem=age.stem[tree.edge[node.pos,2]-Ntip(tree)])
			   if (!is.null(tree$node.label)) {
				  pos <- which(!is.na(tree$node.label[tree$edge[node.pos,2]-Ntip(tree)]) & tree$node.label[tree$edge[node.pos,2]-Ntip(tree)] != "")
				  result[pos,1] <- tree$node.label[tree$edge[node.pos,2]-Ntip(tree)][pos]
					 }
				}
				return(result)
		  }

	#match GND and node age
	nodebase=vector("list",length(region))#nodes with GND>0.65 and their edge, age and GND value
	for (i in 1:length(region)){
		tre.t=ccoquettes[[i]]$phylo   
		root.label <- unique(tre.t$edge[,1][!tre.t$edge[,1] %in% tre.t$edge[,2]])			
		tre.edge <- rbind(tre.t$edge,c(0,root.label))
		ii <- order(tre.edge[,2])
		tre.edge <- tre.edge[ii,]
		pos2=which(tre.edge[,2]>Ntip(tre.t))
		sp.age2=tre.edge[pos2,]		
		spgnd=cbind(sp.age2,gnd[[i]][match(sp.age2[,2],gnd[[i]]$node),])
		
		edge.length <- c(tre.t$edge.length, 0)[ii]		
		node.age=extractNodeAge(tre.t)
		re=cbind(spgnd,node.age)
		nage=na.omit(re[,c(1,"node","GND","crown","stem")])
		colnames(nage)=c("edge1","edge2","GND","age","stem")
		nodebase[[i]]=subset(nage,GND>0.65)				
	}
	saveRDS(nodebase,"nodebase.rds")
	
	nodebase=readRDS("nodebase150.rds")
	sos.re=readRDS("sos.re150.rds")
	ccoquettes=readRDS("ccoquettes150.rds")
# R2 of node to  floristic division 
	ss.glm <- function(r.glm){
		r.ss <- summary(r.glm)
		rsq <- 100*(r.ss$null.deviance-r.ss$deviance)/r.ss$null.deviance
		adj.rsq <- 100*(1-(r.ss$deviance/r.ss$df.residual)/(r.ss$null.deviance/r.ss$df.null))
		f.stat <- ((r.ss$null.deviance-r.ss$deviance)/(r.ss$df.null-
		r.ss$df.residual))/(r.ss$deviance/r.ss$df.residual)
		p <- pf(f.stat, r.ss$df.null-r.ss$df.residual, r.ss$df.residual, lower.tail=FALSE)
		return(c(r2=rsq,adj.r2=adj.rsq,p=p))
		}
					
	r2=vector("list",length(region));names(r2)=region;richfin=vector("list",length(region));names(richfin)=region	
	for (i in 1:length(region)){
		spdis=ccoquettes[[i]]$comm	
		tre <- ccoquettes[[i]]$phylo	
		nodelist.t=sos.re[[i]][2:dim(sos.re[[i]])[1],]
		regionlist.t=regionlist[regionlist[,region[i]]>0,c("adcode_ave2",region[i])]
		nodelist=cbind(regionlist.t,nodelist.t[match(regionlist.t$adcode_ave2,rownames(nodelist.t)),])
		node=colnames(nodelist.t)#GND node edge
		
		r2[[i]]=matrix(0,length(node),3)
		rownames(r2[[i]])=node;	colnames(r2[[i]])=c("r2","p","rich")
			
		region.temp=nodelist[[region[i]]]
		for (j in 1:length(node)){			
			node.temp=nodelist[[node[j]]]
			if (length(node.temp)-sum(is.na(nodelist[[node[j]]]))>1){			
			node.aov=glm(node.temp~region.temp)#方差分析模型，glm()与aov()结果一样				
			r2[[i]][j,1]=ss.glm(node.aov)[[1]]		
			r2[[i]][j,2]=ss.glm(node.aov)[[3]]	
			r2[[i]][j,3]=Ntip(extract.clade(tre,node=as.numeric(node[j])))/Ntip(tre)
			}										
		}	
		age=nodebase[[i]]	
		richfin.t=cbind(age[match(rownames(r2[[i]]),age$edge2),],r2[[i]])
		richfin[[i]]=subset(richfin.t,richfin.t$p<0.05)
		#write.csv(richfin[[i]],paste("richfin_",region[i],".csv"))
	}	
	saveRDS(richfin,"richfin259.rds")
	
## clades with high GND value
	library(ape)
	richfin=readRDS("richfin150.rds")
	ccoquettes=readRDS("ccoquettes150.rds")
	taxonmic=read.csv("Gyo-angiosperm.csv")
	apg=read.csv("APGIV-taxonmic.csv")
	taxon=cbind(taxonmic,APG=apg[match(taxonmic$Family,apg$Family_E),"clades"])
	
	clade.list=c()
	for (i in 1:7){
		tre <- ccoquettes[[i]]$phylo
		richfin.t=subset(richfin[[i]],richfin[[i]]$r2>quantile(richfin[[i]]$r2,0.9))
		checklist=c()
		for (j in 1:length(richfin.t$edge2)){
			node=subset(tre$edge,tre$edge[,1]==richfin.t$edge2[j])
			tre.t=extract.clade(tre,node=richfin.t$edge2[j])
			splist=unique(taxon[match(tre.t$tip.label,taxon$Genus),])
			## cladelist.csv
			tmp=data.frame(id=j,richfin.t[j,c("edge2","GND","age","stem","r2","p")],taxa=paste(na.omit(unique(splist$APG)),collapse =","))
			
			## clade.list.xlsx
			tre1=extract.clade(tre,node=node[1,2])
			tre2=extract.clade(tre,node=node[2,2])
			splist1=unique(taxon[match(tre1$tip.label,taxon$Genus),])
			splist2=unique(taxon[match(tre2$tip.label,taxon$Genus),])			
			if (length(na.omit(unique(splist$APG)))>1) {
				a=paste("(",paste(na.omit(unique(splist1$APG)),collapse =","),")",sep="");
				b=paste("(",paste(na.omit(unique(splist2$APG)),collapse =","),")",sep="")
				tmp=data.frame(id=j,richfin.t[j,c("edge2","GND","age","stem","r2","p")],taxa=paste(a,b,sep=" vs. "))
			} else {
				if (length(na.omit(unique(splist$Order)))>1) {
					a=paste("(",paste(na.omit(unique(splist1$Order)),collapse =","),")",sep="");
					b=paste("(",paste(na.omit(unique(splist2$Order)),collapse =","),")",sep="")
					tmp=data.frame(id=j,richfin.t[j,c("edge2","GND","age","stem","r2","p")],taxa=paste(a,b,sep=" vs. "))
				} else {
					if (length(na.omit(unique(splist$Family)))>1) {
						a=paste("(",paste(na.omit(unique(splist1$Family)),collapse =","),")",sep="");
						b=paste("(",paste(na.omit(unique(splist2$Family)),collapse =","),")",sep="")
						tmp=data.frame(id=j,richfin.t[j,c("edge2","GND","age","stem","r2","p")],taxa=paste(a,b,sep=" vs. "))
					} else {						
						tmp=data.frame(id=j,richfin.t[j,c("edge2","GND","age","stem","r2","p")],taxa=paste("clades within ",na.omit(unique(splist$Family)),sep=""))
					}
				}
			}
			
			checklist=rbind(checklist,tmp)
		}
		tmp2=cbind(boundary=names(richfin)[i],checklist)
		clade.list=rbind(clade.list,tmp2)		
	}
	write.csv(clade.list,"cladelist.csv")
	
	#plot clade contribution
	richfin=readRDS("richfin210.rds")
	thes=5
	plotn=lapply(1:length(region),function(i){
		richfin2=cbind(age.cat=round(richfin[[i]]$age/thes)*thes,richfin[[i]])
		age.cata=sort(unique(richfin2$age.cat))
		#age.cata=age.cata[age.cata<=100]
		r2=do.call(rbind,lapply(1:length(age.cata),function(k){
			dat=richfin2[richfin2$age.cat==age.cata[k],]
			if(dim(dat)[1]>1){
				ci.t=t.test(dat$r2)
				re=data.frame(r2=ci.t$est,r2.l=ci.t$conf[1],r2.h=ci.t$conf[2],r2.se=ci.t$std)
			}else {re=data.frame(r2=mean(dat$r2),r2.l=mean(dat$r2),r2.h=mean(dat$r2),r2.se=0)}
		}))
		rownames(r2)=age.cata
		p.r=cbind(time=age.cata,r2)
		plt=ggplot(p.r)+annotate("rect", xmin=200, xmax=173.5, ymin=0, ymax=6,fill= "#619CFF")+
		annotate("rect", xmin=173.5, xmax=163.5, ymin=0, ymax=6,fill= "#00BCDF")+
		annotate("rect", xmin=163.5, xmax=145, ymin=0, ymax=6,fill= "#8AECFF")+
		annotate("rect", xmin=145, xmax=100, ymin=0, ymax=6,fill= "#7CAE00")+
		annotate("rect", xmin=100, xmax=65.5, ymin=0, ymax=6,fill= "#B6DCB6")+		 
		 annotate("rect", xmin=65.5, xmax=55.8,ymin=0, ymax=6,fill= "#D2E9E1")+
		  annotate("rect", xmin=55.8, xmax=33.9, ymin=0, ymax=6,fill= "#FBEDC9")+
		  annotate("rect", xmin=33.9, xmax=23.03,ymin=0, ymax=6,fill="#F8DDA9")+
		  annotate("rect", xmin=23.03, xmax=2.58, ymin=0, ymax=6,fill="#FCB6D0")+		  
		  annotate("rect", xmin=2.58, xmax=0, ymin=0, ymax=6,fill="#D376FF")+		  
		  annotate("text", x =173.5+(200-173.5)/2 , y =3,label = "J1",size=5)+
		  geom_vline(xintercept=173.5,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =163.5+(173.5-163.5)/2 , y =3,label = "J2",size=5)+
		  geom_vline(xintercept=163.5,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =145+(163.5-145)/2 , y =3,label = "J3",size=5)+
		  geom_vline(xintercept=145,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =100+(145-100)/2 , y =3,label = "K1",size=5)+
		  geom_vline(xintercept=100,col="gray",size=0.5,linetype="longdash",alpha=0.5)+ 		  
		  annotate("text", x =65.5+(100-65.5)/2 , y =3,label = "K2",size=5)+
		  geom_vline(xintercept=65.5,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =55.8+(65.5-55.8)/2, y =3, label = "E1",size=5) +
		  geom_vline(xintercept=55.8,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =33.9+(55.8-33.9)/2, y =3, label = "E2",size=5) +
		  geom_vline(xintercept=33.9,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =23.03+(33.9-23.03)/2, y =3, label = "E3",size=5) +
		  geom_vline(xintercept=23.03,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =2.58+(23.03-2.58)/2, y =3, label = "N",size=5)+		  
		  geom_vline(xintercept=2.58,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =2.58/2, y =3, label = "Q",size=5)	
		plt=plt+theme(panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank())
		if (i%in%c(1,5)) {plt=plt+theme(axis.text.y  = element_text(size=20,color = "black"))
			} else {
			plt=plt+theme(axis.text.y  = element_blank())}	
		if (i>4) {plt=plt+theme(axis.text.x  = element_text(angle=30,vjust=1, hjust = 1,size=20,color = "black"))
		} else {plt=plt+theme(axis.text.x  = element_blank())}
					
	plt=plt+		 
		 geom_ribbon(data =p.r,aes(x = time,ymin=r2-r2.se, ymax=r2+r2.se),fill="gray",alpha=0.5)+
		 geom_line(data=p.r,aes(x=time,y=r2),colour="black",show.legend=TRUE,size=1.2,linetype=1)+
		 scale_x_continuous(limits=c(0,200),expand = c(0,0))+
		 scale_y_continuous(expand = c(0,0),breaks=pretty_breaks(5),limits=c(0,100))
	})
	annotate_figure(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], labels = c("a", "b","c","d", "e","f","g"),label.x=0.05,font.label = list(size = 30),
			nrow=2,ncol = 4,widths=c(1,1),heights=c(1,1,1,1),common.legend=TRUE),
			bottom = text_grob("Time before present (Ma)", size=20),left = text_grob("Contribution of clades on floristic division", rot = 90,size=20))

### effect of current climate on floristic division ####
	library(ggplot2)	
	library(gridExtra)		
	library(cowplot)
	library(ggpubr)	
	library(gdm)
	psim1=get(load("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata"))#matrix
	psim1[is.na(psim1)]=0
	phylo=psim1+t(psim1)	
	clim=read.csv("richclim191024.csv")[,c("adcode","Lon","Lat","MAT","MAP")]
	regionlist=read.csv("cluster.average.merge.csv")
	
	node=colnames(regionlist)[8:14]
	region=as.character(unique(regionlist$phname))		
		
	exFormat=vector("list",(length(node)*2+length(region)));names(exFormat)=c(paste("be_",node,sep=""),node,region)
	glm.r=c()
	for (i in 1:(length(node)*2+length(region))){
		if(i <= length(node)){
			envir=cbind(clim,regionlist[match(clim$adcode,regionlist$adcode_ave2),c("phname",node[i])])				
			envt=subset(envir,envir[,node[i]]>0)			
			physim=as.data.frame(phylo[match(envt$adcode,colnames(phylo)),match(envt$adcode,rownames(phylo))])
			Dissim=cbind(adcode=envt$adcode,physim)
			exFormat_pre <- formatsitepair(Dissim,3, XColumn="Lon", YColumn="Lat", predData=envt,siteColumn="adcode")	
			cole=6+(length(colnames(exFormat_pre))-6)/2	
			cole2=length(colnames(exFormat_pre))
			exFormat[[i]]=subset(exFormat_pre,
				  (exFormat_pre[,cole]==1 & exFormat_pre[,cole2]==2)|(exFormat_pre[,cole]==2 & exFormat_pre[,cole2]==1))			
		}
		if (i>length(node)&i<=length(node)*2){
			envir=cbind(clim,regionlist[match(clim$adcode,regionlist$adcode_ave2),c("phname",node[i-length(node)])])			
			envt=subset(envir,envir[,node[i-length(node)]]>0)				
			physim=as.data.frame(phylo[match(envt$adcode,colnames(phylo)),match(envt$adcode,rownames(phylo))])
			Dissim=cbind(adcode=envt$adcode,physim)
			exFormat_pre <- formatsitepair(Dissim,3, XColumn="Lon", YColumn="Lat", predData=envt,siteColumn="adcode")			
			exFormat[[i]]=exFormat_pre
			
		}		
		if (i>length(node)*2){
			envir=cbind(clim,regionlist[match(clim$adcode,regionlist$adcode_ave2),c("phname","phy2150")])
			envt=subset(envir,envir[,"phname"]==region[i-length(node)*2])					
			physim=as.data.frame(phylo[match(envt$adcode,colnames(phylo)),match(envt$adcode,rownames(phylo))])
			Dissim=cbind(adcode=envt$adcode,physim)
			exFormat_pre <- formatsitepair(Dissim,3, XColumn="Lon", YColumn="Lat", predData=envt,siteColumn="adcode")						
			exFormat[[i]]=exFormat_pre
			
		}
		glm.data=data.frame(phylo=exFormat[[i]]$distance,endis=sqrt((exFormat[[i]]$s1.MAT-exFormat[[i]]$s2.MAT)^2+(exFormat[[i]]$s1.MAP-exFormat[[i]]$s2.MAP)^2))
		r.squ <- data.frame(region=names(exFormat)[i],r2=ss.glm(glm(phylo~endis,glm.data,family = gaussian))[1])
		glm.r=rbind(glm.r,r.squ)	
	}
	glm.r[is.na(glm.r[,2]),2]=0
	write.csv(glm.r,"glm.r.csv")
	
	glm.r=read.csv("glm.r.csv")
	regionlist=read.csv("cluster.average.merge.csv")	
	node=colnames(regionlist)[8:14]	
	hp.plot=vector("list",length(node));names(hp.plot)=node
	plotn=vector("list",length(node))
	hrtc=data.frame(node=sort(rep(1:length(node),3)),region=c("be_World","Laurasian","Gondwanan",
				"be_Laurasian","Holarctic","Saharo-Arabian",
				"be_Gondwanan","Neotropic","Indo.Ocean",
				"be_Neotropic","Neotropical","Chlie-Patagonian",
				"be_Indo.Ocean","Palaeotropic","Antarctic",				
				"be_Antarctic","Australian","Neozeylandic",
				"be_Palaeotropic","Indo-Malesian","African"))
				
	for (i in 1:length(node)){
		hrtc.t=subset(hrtc,hrtc$node==i)
		hp.plot[[i]]=subset(glm.r,glm.r$region%in%hrtc.t$region)
		hp.plot[[i]]$region=factor(hp.plot[[i]]$region,levels=c(as.character(hp.plot[[i]]$region)[1],
			as.character(hp.plot[[i]]$region)[2],as.character(hp.plot[[i]]$region)[3]))#调整x轴顺序，使得所有的图都是红色的区域在前面		
		p=ggplot() +
		geom_bar(data=hp.plot[[i]], aes(x=region,y=r2),fill=c("darkgray","#F8766D","#619CFF"), stat="identity",
			   position='stack')+			  
			   theme_bw() + 			
			theme(axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=10),
			axis.text.y  = element_text(size=10))+
			theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
			  panel.background = element_blank(),axis.line = element_line(colour = "black"))
		p=p+ scale_x_discrete(labels = as.factor(c("Between","A","B"))) 		
		p=p+ylim(0,35)
		plotn[[i]]=p
	}	
	
	figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], 
			nrow=2,ncol = 4,widths=c(2,2),heights=c(2,2,2,2),common.legend=TRUE,legend = "top",label.x=0.5,hjust=0,align="hv")	 	
	annotate_figure(figure, bottom = text_grob("Climate Euclidean distance", size=14),left = text_grob("Rsq of GLM (%)", rot = 90,size=14))	

######## effect of plate tectonics and historical climate on floristic division #########			
	
##caculate geographic isolation over time
	library(gdistance)
	dis=read.csv("Spdat_an_isrm191016.csv")	
	resList <- readRDS("resList.rds")	
	
	# coord=read.csv("POINT2.csv")
	# coord=subset(coord,coord$ADCODE99>0)
	# lon=tapply(coord[,1],coord[,3],mean)
	# lat=tapply(coord[,2],coord[,3],mean)
	# point3=data.frame(adcode=names(lon),lon=round(lon,1),lat=round(lat,1))
	# write.csv(point3,"point3.csv")#手动增加ID，经纬度操作见excel文件	

	point3=read.csv("point3.csv")
	point3=point3[match(dis[,1],point3$adcode),]
	sublist=resList[match(point3$id,names(resList))]	

	pC=vector("list",151);names(pC)=seq(150,0,-1)##每一个列表为所有ADCODE中心点经纬度，共200个Ma的列表
	for (i in 1:151){
	for (j in 1:length(sublist)){
	temp=sublist[[j]][i+51,c(3,2)]
	pC[[i]]=rbind(pC[[i]],temp)
	}
	rownames(pC[[i]])=dis[,1]
	}

	load("2_2_rasters_interpolated_land_1d/environmental_variables.RData")
	nmalt <- colnames(environmental_variables[[1]])[-c(1,2)]
	data2load <- paste0("Elevation_",nmalt,"_T.asc")
	nm_res <- which(is.element(as.numeric(colnames(environmental_variables[[1]])[-c(1,2)]),seq(150,0,-1)))
	data2load <- data2load[nm_res]

	require(parallel)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores); 
	Alt.raster.geo <-do.call(cbind,parLapply(cl=mycl,1:length(data2load),function(i,data2load,pC){
		require(gdistance)	
	  data <- aggregate(raster(paste0("2_2_rasters_interpolated_land_1d/elevation_rasters/",data2load[i])),2)
      data <- crop(data,extent(c(-180,180,-90,90)))		  
	  data[!is.na(data)] <-max(data[!is.na(data)])-data[!is.na(data)]+min(data[!is.na(data)])
	  data[is.na(data)] <- 0#min(data[!is.na(data)])/1000 ##ocean
	  # #geodist without ocean isolation
	  # geoDist <- pointDistance(pC[[i]], longlat=TRUE)
	  # geoDist <- as.dist(geoDist)
	  # rook <- matrix(c(	1,1, 1, 1,1,1,1,1,1,
						# 1,1, 1, 1,1,1,1,1,1,
						# 1,1, 1, 1,1,1,1,1,1,
						# 1,1, 1, 1,1,1,1,1,1,
						# 1,1, 1, 1,0,1,1,1,1,
						# 1,1, 1, 1,1,1,1,1,1,
						# 1,1, 1, 1,1,1,1,1,1,
						# 1,1, 1, 1,1,1,1,1,1,
						# 1,1, 1, 1,1,1,1,1,1), ncol=9, byrow=TRUE)
	  tr <- transition(data,mean, directions=16)#directions=rook	 
	  trC <- geoCorrection(tr, "r", scl=TRUE)
	#plot(raster(trC))	
	  cosDist <- costDistance(trC,as.matrix(pC[[i]]))		  
	  return(as.numeric(cosDist))
	  #return(as.numeric(geoDist))
	},data2load,pC))
	colnames(Alt.raster.geo)=paste0("geo",150:0,sep="")
	stopCluster(mycl)
	
	id=dis[,1]
	f <- function(x, y) paste(x,y,sep="_")
	b=as.vector(outer(id,id, f)[lower.tri(outer(id,id, f))])
	rownames(Alt.raster.geo)=b	
	pastdis=data.frame(id=rownames(Alt.raster.geo),Alt.raster.geo)		
	write.csv(pastdis,"pastdis.new.csv")#cost values: ocean=0;land=elevation;direction=16
	
	## hier.part: phlyo beta ~ clim + geo
	library(data.table)
	library(scales)
	library(ggpubr)
	pastdis0=as.data.frame(fread("pastdis.new.csv"))[,-1]
	psim1=get(load("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata"))#matrix
	psim1[is.na(psim1)]=0
	phylo=psim1+t(psim1)	
	beta.sp=get(load("processed_aV2/psim.sp.Rdata"))	
	sp=betamatrix(beta.sp,as.dist=FALSE)	
	pastdis=data.frame(pastdis0,phylo=phylo[lower.tri(phylo)],sp=sp[lower.tri(sp)])
	
	past.clim=as.data.frame(fread("past.env.clim.csv"))		
	regionlist=read.csv("cluster.average.merge.csv")
	region=colnames(regionlist)[8:14]	
	f <- function(x, y) paste(x,y,sep="_")
	hp.re=lapply(1:length(region),function(i){
		#for (i in 1:length(region)){
		region1=regionlist[regionlist[,region[i]]==1,"adcode_ave2"]
		region2=regionlist[regionlist[,region[i]]==2,"adcode_ave2"]	
		adcd <- c(as.vector(outer(region1,region2, f)),as.vector(outer(region2,region1,f)))	
		across.geo=pastdis[pastdis$id%in%adcd,c("id","phylo","sp",paste("geo",150:0,sep=""))]		
		mf=function(x) unlist(strsplit(x,split='_'))
		a=do.call(rbind,lapply(across.geo$id,mf))
		across.geo2=rbind(across.geo,cbind(id=paste(a[,2],a[,1],sep="_"),across.geo[,-1]))
		across.clim=past.clim[past.clim$ID%in%adcd,c("ID",paste("clim_",150:0,sep=""))]
		across.com=cbind(across.geo2,across.clim[match(across.geo2$id,across.clim$ID),])
		across=across.com[!is.na(across.com$ID),-(length(across.geo)+1)]
		# hier.part
		hp=do.call(rbind,lapply(150:0,function(j){
		#for (j in 150:0){
			dat=na.omit(across[,c("phylo","sp",paste("geo",j,sep=""),paste("clim_",j,sep=""))])
			dat=dat[dat$geo>0&dat$clim>0,]			
			colnames(dat)=c("phylo","sp","geo","clim")
			dat2=dat[dat$geo!=Inf,]			
			if(dim(dat)[1]>1){				
				ci.clim=t.test(dat$clim)
				clim=data.frame(time=j,clim.mean=ci.clim$est,clim.l=ci.clim$conf[1],clim.h=ci.clim$conf[2],clim.se=ci.clim$std)						
				if(length(unique(dat2$geo))>1){
					cors=data.frame(cor.geo.clim=cor(scale(dat2[,c("geo","clim")]))[1,2],
						cor_p=ggcorrplot::cor_pmat(scale(dat2[,c("geo","clim")]))[1,2])
					ci.geo=t.test(dat2$geo)
					geo=data.frame(geo.mean=ci.geo$est,geo.l=ci.geo$conf[1],geo.h=ci.geo$conf[2],geo.se=ci.geo$std)
					phylo=hier.part::hier.part(log(dat2[,"phylo"]),scale(dat2[,c("geo","clim")]),gof = "Rsqu",barplot = FALSE)
					sp=hier.part::hier.part(dat2[,"sp"],dat2[,c("geo","clim")],gof = "Rsqu",barplot = FALSE)					
					re=rbind(data.frame(hp=phylo$IJ["geo","I"],hp.type="geo",betatype="phylo",clim,geo,cors),
						data.frame(hp=phylo$IJ["clim","I"],hp.type="clim",betatype="phylo",clim,geo,cors),
						data.frame(hp=phylo$IJ[1,"J"],hp.type="joint",betatype="phylo",clim,geo,cors),
						data.frame(hp=sp$IJ["geo","I"],hp.type="geo",betatype="sp",clim,geo,cors),
						data.frame(hp=sp$IJ["clim","I"],hp.type="clim",betatype="sp",clim,geo,cors),
						data.frame(hp=sp$IJ[1,"J"],hp.type="joint",betatype="sp",clim,geo,cors))
					# using correlation r
					# re=rbind(data.frame(hp=cor(dat2$phylo,dat2$geo),hp.type="geo",betatype="phylo",clim,geo),
						# data.frame(hp=cor(dat$phylo,dat$clim),hp.type="clim",betatype="phylo",clim,geo),						
						# data.frame(hp=cor(dat2$sp,dat2$geo),hp.type="geo",betatype="sp",clim,geo),
						# data.frame(hp=cor(dat$sp,dat$clim),hp.type="clim",betatype="sp",clim,geo))
				}else{re=data.frame(hp=NA,hp.type=NA,betatype=NA,clim,geo.mean=NA,geo.l=NA,geo.h=NA,geo.se=NA,cor.geo.clim=NA,cor_p=NA)}												
			}else{				
				re=data.frame(time=j,hp=NA,hp.type=NA,betatype=NA,
				geo.mean=mean(dat$geo),geo.l=mean(dat$geo),geo.h=mean(dat$geo),geo.se=0,
				clim.mean=mean(dat$clim),clim.l=mean(dat$clim),clim.h=mean(dat$clim),clim.se=0,cor.geo.clim=NA,cor_p=NA)
			}		
		#}			
			return(re)
		 }))
		#}		
		return(hp)
	})	
	names(hp.re)=region
	cor.stat=function(tmp){
		dat=na.omit(unique(tmp[tmp$cor_p<=0.05,"cor.geo.clim"]))
		ci.t=t.test(dat)
		re=data.frame(cor.mean=ci.t$est,cor.l=ci.t$conf[1],cor.h=ci.t$conf[2],cor.se=ci.t$std)
		return(re)
	}
	cor.realm=do.call(rbind,lapply(hp.re,cor.stat))		
	
	rownames(cors)=age.cata
	
	data.p=function(hp.dat,thes){
		if (thes>1){
		p0=cbind(age.cat=round(hp.dat$time/thes)*thes,hp.dat)
		age.cata=sort(unique(p0$age.cat))
		#age.cata=age.cata[age.cata<=100]
		geo=do.call(rbind,lapply(1:length(age.cata),function(k){
		#for (k in 1:length(age.cata)){
		dat=p0[p0$age.cat==age.cata[k],]
		dat2=na.omit(dat)
		if(length(unique(dat2$geo.mean))>1){			
			ci.t=t.test(dat2$geo.mean)
			re=data.frame(time=age.cata[k],geo.mean=ci.t$est,geo.l=ci.t$conf[1],geo.h=ci.t$conf[2],geo.se=ci.t$std)
			}else {re=data.frame(time=age.cata[k],geo.mean=mean(dat2$geo.mean),geo.l=mean(dat2$geo.mean),geo.h=mean(dat2$geo.mean),geo.se=0)}
		#}		
		}))
		rownames(geo)=age.cata
		clim=do.call(rbind,lapply(1:length(age.cata),function(k){
		dat=p0[p0$age.cat==age.cata[k],]
		if(dim(dat)[1]>1){
			ci.t=t.test(dat$clim.mean)
			re=data.frame(clim.mean=ci.t$est,clim.l=ci.t$conf[1],clim.h=ci.t$conf[2],clim.se=ci.t$std)
			}else {re=data.frame(clim.mean=mean(dat$clim.mean),clim.l=mean(dat$clim.mean),clim.h=mean(dat$clim.mean),clim.se=0)}
		}))
		rownames(clim)=age.cata
		p=cbind(geo,clim)	
		
		p2=na.omit(p0[p0$betatype%in%"phylo",])			
		geo.hp=do.call(rbind,lapply(1:length(age.cata),function(k){
			dat=p2[p2$age.cat==age.cata[k]&p2$hp.type=="geo",c("age.cat","hp")]
			if(dim(dat)[1]>1){
				ci.t=t.test(dat$hp)
				re=data.frame(iso.mean=ci.t$est,iso.l=ci.t$conf[1],iso.h=ci.t$conf[2],iso.se=ci.t$std)
			}else {re=data.frame(iso.mean=mean(dat$hp),iso.l=mean(dat$hp),iso.h=mean(dat$hp),iso.se=0)}
		}))
		rownames(geo.hp)=age.cata
		clim.hp=do.call(rbind,lapply(1:length(age.cata),function(k){
			dat=p2[p2$age.cat==age.cata[k]&p2$hp.type=="clim",c("age.cat","hp")]
			if(dim(dat)[1]>1){
				ci.t=t.test(dat$hp)
				re=data.frame(iso.mean=ci.t$est,iso.l=ci.t$conf[1],iso.h=ci.t$conf[2],iso.se=ci.t$std)
			}else {re=data.frame(iso.mean=mean(dat$hp),iso.l=mean(dat$hp),iso.h=mean(dat$hp),iso.se=0)}
		}))
		rownames(clim.hp)=age.cata		
		joint.hp=do.call(rbind,lapply(1:length(age.cata),function(k){
			dat=p2[p2$age.cat==age.cata[k]&p2$hp.type=="joint",c("age.cat","hp")]
			if(dim(dat)[1]>1){
				ci.t=t.test(dat$hp)
				re=data.frame(iso.mean=ci.t$est,iso.l=ci.t$conf[1],iso.h=ci.t$conf[2],iso.se=ci.t$std)
			}else {re=data.frame(iso.mean=mean(dat$hp),iso.l=mean(dat$hp),iso.h=mean(dat$hp),iso.se=0)}
		}))
		rownames(joint.hp)=age.cata
		
		p.iso=rbind(cbind(time=age.cata,geo.hp,iso.typ="geo"),cbind(time=age.cata,clim.hp,iso.typ="clim"),
			cbind(time=age.cata,joint.hp,iso.typ="joint"))
		p.iso$iso.typ=factor(p.iso$iso.typ,levels=c("geo","clim","joint"),ordered=TRUE)
		
		cors=do.call(rbind,lapply(1:length(age.cata),function(k){
			#for (k in 1:length(age.cata)){
			dat=na.omit(unique(p0[p0$cor_p<=0.05&p0$age.cat==age.cata[k],c("age.cat","cor.geo.clim","cor_p")]))
			if(dim(dat)[1]>1){
				ci.t=t.test(dat$cor.geo.clim)
				re=data.frame(cor.mean=ci.t$est,cor.l=ci.t$conf[1],cor.h=ci.t$conf[2],cor.se=ci.t$std)
			}else {re=data.frame(cor.mean=mean(dat$cor.geo.clim),cor.l=mean(dat$cor.geo.clim),cor.h=mean(dat$cor.geo.clim),cor.se=0)}
			#}			
		}))
		rownames(cors)=age.cata	
		
		p.geo=cbind(geo,geo.hp)
		p.clim=cbind(time=age.cata,clim,clim.hp)				
		p.fin=list(p,p.iso,p.geo,p.clim,cors)
		names(p.fin)=c("isolation","hier.r2","geo&r2","clim&r2","cor.geo.clim")
		}
		if (thes==1) {
		p=hp.dat[,c("time","geo.mean","geo.l","geo.h","geo.se","clim.mean","clim.l","clim.h","clim.se")]
		p.iso0=subset(hp.dat[,c("time","hp","hp.type","betatype")],hp.dat$betatype%in%"phylo"&!hp.dat$hp.type%in%"joint")		
		colnames(p.iso0)=c("time","iso.mean","iso.typ","betatype")		
		p.iso=cbind(p.iso0[,-4],iso.l=p.iso0$iso.mean,iso.h=p.iso0$iso.mean,iso.se=0)
		
		roll.mean=function(iso){
			require(zoo)
			x=iso$iso.mean
			names(x)=iso$time
			rollmean=rollmean(x,5)
			return(rollmean)
		}
		clim.roll=subset(p.iso[,c("time","iso.mean","iso.typ")],p.iso$iso.typ%in%"clim")
		geo.roll=subset(p.iso[,c("time","iso.mean","iso.typ")],p.iso$iso.typ%in%"geo")
		p.roll0=rbind(data.frame(time=as.numeric(names(roll.mean(geo.roll))),iso.mean=roll.mean(geo.roll),iso.typ="geo"),
			data.frame(time=as.numeric(names(roll.mean(clim.roll))),iso.mean=roll.mean(clim.roll),iso.typ="clim"))		
		p.roll=cbind(p.roll0,iso.l=p.roll0$iso.mean,iso.h=p.roll0$iso.mean,iso.se=0)
		p.fin=list(p,p.iso,p.roll)
		names(p.fin)=c("isolation","hier.r2","hier.r2.roll")
		}
		return(p.fin)
	}	
	datap=lapply(hp.re,data.p,5)	
	#for (i in region)data.p(hp.re[[i]],5)
	
	#plot HP R2		
	plot.hp=function(pre){
		p.iso=na.omit(pre[["hier.r2"]])	
		p.iso=p.iso[p.iso$time<=80,]
		plt=ggplot(p.iso)+#annotate("rect", xmin=145, xmax=100, ymin=-0.16, ymax=-0.1,fill= "#7CAE00")+
		annotate("rect", xmin=80, xmax=65.5, ymin=-25, ymax=-10,fill= "#B6DCB6")+		 
		 annotate("rect", xmin=65.5, xmax=55.8,ymin=-25, ymax=-10,fill= "#D2E9E1")+
		  annotate("rect", xmin=55.8, xmax=33.9, ymin=-25, ymax=-10,fill= "#FBEDC9")+
		  annotate("rect", xmin=33.9, xmax=23.03,ymin=-25, ymax=-10,fill="#F8DDA9")+
		  annotate("rect", xmin=23.03, xmax=2.58, ymin=-25, ymax=-10,fill="#FCB6D0")+		  
		  annotate("rect", xmin=2.58, xmax=0, ymin=-25, ymax=-10,fill="#D376FF")+		  
		  # annotate("text", x =100+(145-100)/2 , y =-0.13,label = "K1",size=3)+
		  # geom_vline(xintercept=100,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =65.5+(80-65.5)/2 , y = -17.5,label = "K",size=5)+
		  geom_vline(xintercept=65.5,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =55.8+(65.5-55.8)/2, y =-17.5, label = "E1",size=5) +
		  geom_vline(xintercept=55.8,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =33.9+(55.8-33.9)/2, y =-17.5, label = "E2",size=5) +
		  geom_vline(xintercept=33.9,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =23.03+(33.9-23.03)/2, y = -17.5, label = "E3",size=5) +
		  geom_vline(xintercept=23.03,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =2.58+(23.03-2.58)/2, y =-17.5, label = "N",size=5)+		  
		  geom_vline(xintercept=2.58,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =2.58/2, y =-17.5, label = "Q",size=5)	
		plt=plt+theme(panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=20,color = "black"),
			axis.text.y  = element_text(size=20,color = "black"))	
		plt=plt+geom_line(data=na.omit(p.iso),aes(x=time,y=iso.mean*100,colour=iso.typ),size=1.5,show.legend=TRUE)+
			geom_ribbon(data = na.omit(p.iso),aes(x = time,ymin=iso.l*100, ymax=iso.h*100,fill=iso.typ),alpha=0.4)+
			scale_color_manual(values = c("#00BA38","#619CFF","#F8766D"))+
			scale_fill_manual(values = c("#00BA38","#619CFF","#F8766D"))+ ylim(-25,100)					
		return(plt)
	}
	plotn=lapply(datap,plot.hp)
	annotate_figure(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], labels = c("a", "b","c","d", "e","f","g"),label.x=0.1,
			nrow=2,ncol = 4,widths=c(1,1),heights=c(1,1,1,1),common.legend=TRUE),
			bottom = text_grob("Time before present(Ma)", size=20),left = text_grob("Partial R2", rot = 90,size=20))
	
	#plot cor of clim and geo
	plot.cor=function(realm,pre,cor.realm){
		p.iso=na.omit(pre[[realm]][["cor.geo.clim"]])
		p.iso$time=as.numeric(rownames(p.iso))
		p.iso=p.iso[p.iso$time<=80,]		
		plt=ggplot(p.iso)+theme(panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=20,color = "black"),
			axis.text.y  = element_text(size=20,color = "black"))
		plt=plt+annotate("rect", xmin=80, xmax=65.5, ymin=-1, ymax=-0.75,fill= "#B6DCB6")+		 
		 annotate("rect", xmin=65.5, xmax=55.8,ymin=-1, ymax=-0.75,fill= "#D2E9E1")+
		  annotate("rect", xmin=55.8, xmax=33.9, ymin=-1, ymax=-0.75,fill= "#FBEDC9")+
		  annotate("rect", xmin=33.9, xmax=23.03,ymin=-1, ymax=-0.75,fill="#F8DDA9")+
		  annotate("rect", xmin=23.03, xmax=2.58, ymin=-1, ymax=-0.75,fill="#FCB6D0")+		  
		  annotate("rect", xmin=2.58, xmax=0, ymin=-1, ymax=-0.75,fill="#D376FF")+		  
		  # annotate("text", x =100+(145-100)/2 , y =-0.13,label = "K1",size=3)+
		  # geom_vline(xintercept=100,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =65.5+(80-65.5)/2 , y = -0.875,label = "K",size=6)+
		  geom_vline(xintercept=65.5,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =55.8+(65.5-55.8)/2, y =-0.875, label = "E1",size=6) +
		  geom_vline(xintercept=55.8,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =33.9+(55.8-33.9)/2, y =-0.875, label = "E2",size=6) +
		  geom_vline(xintercept=33.9,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =23.03+(33.9-23.03)/2, y = -0.875, label = "E3",size=6) +
		  geom_vline(xintercept=23.03,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =2.58+(23.03-2.58)/2, y =-0.875, label = "N",size=6)+		  
		  geom_vline(xintercept=2.58,col="gray",size=0.5,linetype="longdash",alpha=0.5)+
		  annotate("text", x =2.58/2, y =-0.875, label = "Q",size=6)	+
		  geom_hline(yintercept=0,col="darkred",size=1,linetype="longdash",alpha=0.8)+
		  annotate("text", x =50, y =0.8, label = paste("Cor = ",round(cor.realm[realm,"cor.mean"],2)," ± ",round(cor.realm[realm,"cor.se"],2),sep=""),size=6)
		plt=plt+geom_line(aes(x=time,y=cor.mean),colour="black",size=2,show.legend=TRUE)+
			geom_ribbon(aes(x = time,ymin=cor.l, ymax=cor.h),fill="gray",alpha=0.2)+ylim(-1,1)
		
		}	
	plotn=lapply(names(datap),plot.cor,datap,cor.realm)
	annotate_figure(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], labels = c("a", "b","c","d", "e","f","g"),
			nrow=2,ncol = 4,widths=c(1,1),heights=c(1,1,1,1),common.legend=TRUE),
			bottom = text_grob("Time before present(Ma)", size=20),left = text_grob("Correlation index", rot = 90,size=20))
		
