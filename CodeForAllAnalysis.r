
#谱系多样性

library(ape)
text="((c:25,((d:8,e:8):7,f:15)G2:10)F1:25,(a:10,b:10,g:10)G1:40)N18031:50;"#括号：共同祖先；冒号加数字：枝长；小写字母Tip/物种名;大写字母：科属名；分号：一棵树的结尾；逗号：姊妹种
text="((c:25,((d:8,e:8):7,f:15):10):25,(a:10,b:10,g:10):40):50;"
text="((c:10,d:10):10,(a:15,b:15):5)Node:20;"
tree=read.tree(text=text)#读取谱系数据
plot(tree,show.tip.label=F,show.node.label=F)#显示谱系树，包含物种和科属的名字
axisPhylo()#添加枝长的坐标轴
tree=read.tree("output_tree.tre")

attributes(tre)
tre$edge#物种名和科属名的编号，编号原则从下往上，先编物种名，再从高级到低级编科属名，[,1]表示树的上一级节点所在的编号，[,2]表示本节点所在的编号
tre$tip.label#物种名
tre$node.label#科属名，未添加科属名的显示为空
tree$edge.length#枝长
Ntip(tre)#返回谱系树所包含的物种数

#提取每个物种的枝长
pos=which(tre$edge[,2]<=Ntip(tre))#提取edge中每个物种对应的编号所在的行号
sp.age=tre$edge.length[pos]#返回每个物种的枝长
names(sp.age)=tre$tip.label[tre$edge[pos,2]]


#pos2=subset(tre$edge,tre$edge[,2]<=Ntip(tre))


#提取物种及其上一节点各自的枝长
 pos2=tre$edge[pos,1]#每个物种的上一节点的编号
pos3=c()#每个物种的上一节点的编号所在的行号
 for (i in 1:length(pos2)){
     p=which(tre$edge[,2]==pos2[i])
	 pos3=cbind(pos3,p)}
 sp.age2=tre$edge.length[pos3]#返回每个物种上一节点的枝长

di <- phylo.dist(tre, tip1=tre$tip.label,method = "branch.length")

tre$tip.label
#"c" "g2" "g2" "g2" "g1" "g1" "g1"
tre2 <- tree.backbone(tre, c("c", "g2", "g2", "g2", "g1", "g1", "g1"))

tree.routine(tre,file="e:/tree.csv")

#提取一个树种某几个种建立一个小树
tre3 <- drop.tip(tre, tip=c("a","b", "g"))
plot(tre3,show.tip.label=T,show.node.label=T)


tre3 <- drop.tip(tre, tip=tre$tip.label[15:15035])
plot(tre3)
axisPhylo()#添加枝长的坐标轴

tre2=read.tree("140_210_constrained3_CL0.0001_full_tree_age.tre")
class(tre)
plot(tre,show.tip.label=T,show.node.label=T)#显示谱系树，包含物种和科属的名字
axisPhylo()#添加枝长的坐标轴
attributes(tre)
tre$edge#物种名和科属名的编号，编号原则从下往上，先编物种名，再从高级到低级编科属名，[,1]表示树的上一级节点所在的编号，[,2]表示本节点所在的编号
tre$tip.label#物种名
tre$node.label#科属名，未添加科属名的显示为空
tre$edge.length#枝长
Ntip(tre)#返回谱系树所包含的物种数

#提取每个物种的枝长
pos=which(tre$edge[,2]<=Ntip(tre))#提取edge中每个物种对应的编号所在的行号
sp.age=tre$edge.length[pos]#返回每个物种的枝长
names(sp.age)=tre$tip.label[tre$edge[pos,2]]
write.csv(sp.age,"tre_length.csv")


#谱系beta
library(picante)
help(package="picante")
data(phylocom)
phylosor(phylocom$sample, phylocom$phylo)
 plot(phylocom$phylo)
 head(phylocom$sample)
 #随机化分布矩阵，建立null model
 data(phylocom)
randomizeMatrix(phylocom$sample, null.model="richness")
#随机化分布矩阵，并求算null model的谱系beta
data(phylocom)
phylosor.rnd(phylocom$sample,phylocom$phylo,cstSor=TRUE,null.model="richness",runs=5)
#计算谱系和物种多样性
data(phylocom)
pd(phylocom$sample, phylocom$phylo)
#Calculates MNTD (mean nearest taxon distance) separating taxa in two communities, a measure of phylogenetic beta diversity
data(phylocom)
comdistnt(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=FALSE)
#Calculates MPD (mean pairwise distance) separating taxa in two communities, a measure of phylogenetic beta diversity
data(phylocom)
comdist(phylocom$sample, cophenetic(phylocom$phylo), abundance.weighted=TRUE)


####GTS数据整理 http://www.bgci.org/global_tree_search.php?sec=globaltreesearch
##合并一个文件夹下所有csv文件并增加第一列存储文件名
filelist <- list.files(pattern=".*.csv")
country=sub(pattern = ".csv", replacement = "", filelist)
datafr=c()
for (i in 1:length(filelist)){
   datalist=cbind(country=country[i],read.csv(filelist[i]))
   datafr=rbind(datafr,datalist)
   }
write.csv(datafr,"datafr.csv")
##获得GTS数据与adcode地理单元的对应关系
geo=read.csv("Adcode-GTS.csv")
gts=read.csv("datafr.csv")
gts.adc=cbind(gts,geo[match(gts$country,geo$GTS),])
write.csv(na.omit(gts.adc),"GTSdata.csv")

#补充数据的整理
joint.table2 <- function(data1, index1, data2, index2, case = TRUE, ...)
	{
	if (!is.data.frame(data1))
		data1 <- as.data.frame(data1,stringsAsFactors=F)
	if (!is.data.frame(data2))
		data2 <- as.data.frame(data2,stringsAsFactors=F)	
	
	if (case) index1 <- as.character(data1[[index1]])
	else index1 <- tolower(as.character(data1[[index1]]))
	
	if (class(index2)=="numeric") {
		select <- c(1:(dim(data2)[2]))[-index2]
		select.name <- names(data2)[-index2]
		}
	else if (class(index2)=="character") {
		select <- names(data2)[!names(data2) %in% index2]
		select.name <- select
		}
	else {stop("error1")}

	if (case) index2 <- as.character(data2[[index2]])
	else index2 <- tolower(as.character(data2[[index2]]))
	data2 <- data2[select]
	
	cond <- which((!duplicated(index2)) & (index2 %in% index1) & (!is.na(index2)))
	if (length(cond)>0) {
		index2 <- index2[cond]
		if (length(select)==1) data2 <- as.data.frame(data2[cond,], stringsAsFactors=FALSE)
		else data2 <- data2[cond,]
		}
	names(data2) <- select.name
	rm(cond)

	n <- length(index1)	
	m <- length(index2)
	ii <- order(index1)
	jj <- order(index2)
	index1.1 <- index1[ii]
	index2.1 <- index2[jj]
	
	i <- j <- 1
	pos1 <- pos2 <- numeric(n); pos1[] <- pos2[] <- NA
	while(i <= n & j <= m) {
		if (is.na(index1.1[i])) i <- i + 1
		else if (index1.1[i] < index2.1[j]) i <- i + 1
		else if (index1.1[i] > index2.1[j]) j <- j+1
		else {pos2[i] <- jj[j]; i <- i + 1}
		}
	pos1[ii] <- pos2

	if (all(is.na(pos1))) {stop("no shared row in the two tables")}
	else {
		newdata <- vector(mode="list",length=length(select))
		names(newdata) <- colnames(data2)
		for (i in 1:length(select)) newdata[[i]] <- data2[pos1,i]
		data1 <- data.frame(data1, newdata, stringsAsFactors=FALSE)
		return(data1)
		}
	}

load("All_Sp_table_Acc_161113.R")
#load("RAINBIO.RData")
d.rain=read.csv("RAINBIO.csv")
d2.rain<-read.csv("RAINBIO_not_or_wrong_georef.csv")
d=rbind(d.rain[,2:26],d2.rain[,2:26])
d2 <- joint.table2(d, "species", sp.list.acc[,c("SpName_1", "Taxonomic.status.in.TPL", "ScientificName_1")], "SpName_1")

d.na=subset(d2,is.na(d2$Taxonomic.status.in.TPL))
write.csv(d.na, file="unmatch.csv")
d.ok=subset(d2,!is.na(d2$Taxonomic.status.in.TPL))
write.csv(d.ok, file="dok.csv")

nam=read.csv("results.csv")
d.im=read.csv("unmatch.csv")
d.mo=cbind(d.im,nam[match(d.im$species2,nam$original.search),])
d.mo2=subset(d.mo,!is.na(d.mo$id))
d.ok2<- joint.table2(d.mo2, "spname", sp.list.acc[,c("SpName_1", "Taxonomic.status.in.TPL", "ScientificName_1")], "SpName_1")
write.csv(d.ok2,"dok2.csv")#手动整理了列标题

d.ok2=read.csv("dok2.csv")
d.ok=read.csv("dok.csv")
d.fin=rbind(d.ok[,2:28],d.ok2[,c(2:26,41:42)])
geo=read.csv("Geonames.csv")
re=cbind(d.fin,geo[match(d.fin$country,geo$NAME99),])
#re2=subset(re,(!is.na(re$NAME99))|(!is.na(re$decimalLongitude)))
write.csv(re, file="RAINBIO_yp180505.csv")

#分布数据整理180522
spdis=read.csv("Woody_SpDistribution.csv")
spdel=read.csv("Deleted_SpDistribution.csv")
add=read.csv("Additional_distribution.csv")
addel=read.csv("Deleted_AddDistribution.csv")

temp <- cbind(spdis, spdel[match(spdis$match,spdel$match),c("Adcode99", "Species_E1","Del_Time")])
temp2 <- cbind(add, addel[match(add$match,addel$match),c("Adcode99", "Species_E1","Del_Time")])
spnew=subset(temp[,1:20],is.na(temp$Del_Time))
addnew=subset(temp2[,1:18],is.na(temp2$Del_Time))
write.csv(spnew,"Woody_SpDistribution2.csv")#转换成Woody_SpDistribution.xlsx
write.csv(addnew,"Additional_distribution2.csv")#转换成Additional_distribution.xlsx

spdata=unique(rbind(spnew[,c(1:5,7:18)],addnew[,c(1:4,6:18)]))
write.csv(spdata,"spdis180522.csv")

##校正地名
	geo=read.csv("Link table old to new-yp.csv")#在Link table old to new.csv的基础上，补全了未改动的名字（根据spdis180522.csv的adcode99）
	dis=read.csv("spdis180522.csv")
	dis.cor=cbind(geo[match(dis$Adcode99,geo$Adcode99old),],dis)
	spdis=dis.cor[,c(1:6,8:9,13:24)]	
	write.csv(spdis,"Woody_SpDistribution.csv")#这是截止18年5月22日的分布数据（目前最新），以下分析都基于这个表计算
	
	

#####2019.10.22更新分布数据和制作分布矩阵
	#哈萨克斯坦
	spdisc2spdisgrid_yp <- function(adcode, grid, grid.county)
		{
		area <- grid.county[grid.county$Adcode99 %in% adcode,c("Grid","inter_area","GEOarea")]
		grid.area <- tapply(X=area$inter_area,INDEX=area$Grid,FUN=sum)
		grid.sp <- as.numeric(names(grid.area))
		grid.area <- as.numeric(grid.area)
		
		area.threshold=as.numeric(unique(area$GEOarea,INDEX=area$Grid,FUN=mean))
		spdis.grid <- grid.sp[grid.area/area.threshold >0.5]
		
		cond <- grid %in% spdis.grid
		sprich.tree[cond] <<- sprich.tree[cond]+1
		
		return(as.numeric(cond))
		}

		grid.county <- read.csv("Kaz_geo.csv")
		  tree.dis.c<-read.csv("kaz_dis2.csv")


		#sum(!unique(sp.dis.c$Adcode99) %in% unique(grid.county$Adcode99))
		grid <- unique(grid.county$Grid)
		
		sprich.tree <- numeric(length(grid))
		tree.dis.grid <- tapply(X=tree.dis.c$Adcode99,INDEX=as.character(tree.dis.c$Genus),FUN=spdisc2spdisgrid_yp,
		  grid=grid, grid.county=grid.county)

		sprich.tree[sprich.tree==0] <- NA


	 spname=as.character(attributes(tree.dis.grid)$dimnames[[1]])

	spdat.kaz=c()
	for (i in 1:length(spname))
	{
	sp.dat=tree.dis.grid[[i]]
	spdat.kaz=cbind(spdat.kaz,sp.dat)
	}
	colnames(spdat.kaz)=c(spname)
	rownames(spdat.kaz)=c(grid)

	#东亚
	spdisc2spdisgrid <- function(adcode, grid, grid.county, area.threshold = 8000000000)
		{
		area <- grid.county[grid.county$Adcode99 %in% adcode,c("Grid","Insert_area")]
		grid.area <- tapply(X=area$Insert_area,INDEX=area$Grid,FUN=sum)
		grid.sp <- as.numeric(names(grid.area))
		grid.area <- as.numeric(grid.area)
		
		spdis.grid <- grid.sp[grid.area > area.threshold]
		
		cond <- grid %in% spdis.grid
		sprich.tree[cond] <<- sprich.tree[cond]+1
		
		return(as.numeric(cond))
		}

		grid.county <- read.csv("grid_geo.csv")
		  tree.dis.c<-read.csv("esia.csv")


		#sum(!unique(sp.dis.c$Adcode99) %in% unique(grid.county$Adcode99))
		grid <- unique(grid.county$Grid)
		
		sprich.tree <- numeric(length(grid))
		tree.dis.grid <- tapply(X=tree.dis.c$Adcode99,INDEX=as.character(tree.dis.c$spname),FUN=spdisc2spdisgrid,
		  grid=grid, grid.county=grid.county, area.threshold = 8000000000)

		sprich.tree[sprich.tree==0] <- NA


	 spname=as.character(attributes(tree.dis.grid)$dimnames[[1]])

	spdat.ea=c()
	for (i in 1:length(spname))
	{
	sp.dat=tree.dis.grid[[i]]
	spdat.ea=cbind(spdat.ea,sp.dat)
	}
	colnames(spdat.ea)=c(spname)
	rownames(spdat.ea)=c(grid)

	#将东亚和哈萨克斯坦的分布数据整理成两列并输出
	spdis.t=c()
	for (i in 1:dim(spdat.kaz)[1]){
		for(j in 1:dim(spdat.kaz)[2]){
		if (spdat.kaz[i,j]==1){
		temp=data.frame(adcode=rownames(spdat.kaz)[i],genus=colnames(spdat.kaz)[j])
		spdis.t=rbind(spdis.t,temp)}
		}
	}
	for (i in 1:dim(spdat.ea)[1]){
		for(j in 1:dim(spdat.ea)[2]){
		if (spdat.ea[i,j]==1){
		temp=data.frame(adcode=rownames(spdat.ea)[i],genus=colnames(spdat.ea)[j])
		spdis.t=rbind(spdis.t,temp)}
		}
	}
	colnames(spdis.t)=c("Adcode99","Species_E1")
	write.csv(unique(spdis.t),"spdis_kaz_EA.csv")

	newdis=read.csv("spdis_kaz_EA.csv")[,c("Adcode99","Species_E1")]
	dis.t=unique(read.csv("Woody_SpDistribution.csv")[,c("Adcode99","Species_E1")])
	dis=unique(rbind(newdis,dis.t))
	write.csv(dis,"SpDistribution20191016.csv")#目前最新全球物种分布

# ###整合Plants of the world online(POWO)的introduced数据,校对全球分布数据
	# ##校对powo的地名
	# geoname=read.csv("Geo-isrm.csv")
	# powo.nm=as.character(unique(powo$geoname))
	# powo.nm2=cbind(powo.nm,geoname[match(powo.nm,geoname$Province),],geoname[match(powo.nm,geoname$NAME99),])
	# write.csv(powo.nm2,"geo-isrm-powo.csv")#之后手工校对
	# dim(na.omit(powo.nm2))

	# library(data.table);library(dplyr);library(tidyr)
	# geoname=read.csv("geo-isrm-powo.csv")#powo地名与adcode对应表
	# powo=as.data.frame(fread("PTWDistribution.txt"))[,-1]
	# temp=full_join(powo,geoname,by = c("geoname" = "powo.nm"))
	# powo2=temp[!is.na(temp$ADCODE99),]
	# powo3=separate(powo2,"sp","genus",sep=" ",remove=FALSE)
	# powo.t=read.csv("powo.t.csv")#属名前面有乱码的手动整理的表
	# powo4=rbind(powo.t,powo3)
	# powo5=unique(powo4[,-1])	
	# write.csv(powo5,"powo.csv")#之后去除乱码，得到genus

	# ##将属内各个物种的status归纳为一个
	# library(data.table)
	# powo=as.data.frame(fread("powo.csv"))
	# adcode=unique(powo$ADCODE99)
	# getGenus=function(powo,adc){
		# powo.nz=subset(powo,powo$ADCODE99==adc)
		# mf=function(x){length(unique(x))}
		# intro.nz=tapply(powo.nz$status,powo.nz$genus,mf)
		# powo.nz2=cbind(powo.nz,cal=intro.nz[match(powo.nz$genus,names(intro.nz))])
		# powo.nz2[powo.nz2$cal==2,"status"]="native"
		# powo.nz3=unique(powo.nz2[,c("genus","status","ADCODE99","NAME99")])
		# powo.nz4=cbind(ID=paste(powo.nz3$ADCODE99,powo.nz3$genus,sep="_"),powo.nz3)
		# return(powo.nz4)
	# }
	# re=c()
	# for(i in 1:length(adcode)){
		# tmp=getGenus(powo,adcode[i])
		# re=rbind(tmp,re)
	# }
	# write.csv(re,"powo.ge.csv")

	# ##利用powo剔除全球数据库塔斯马尼亚和新西兰的introduced数据
	# library(data.table)
	# powo.ge=as.data.frame(fread("powo.ge.csv"))[,-1]
	# intro.nz=powo.ge[powo.ge$status%in%"introduced"&powo.ge$ADCODE99%in%c(280,413),]
	# dis=as.data.frame(fread("SpDistribution20191016.csv"))[,-1]#分布数据中有重复数据
	# dis2=cbind(ID=paste(dis$Adcode99,dis$Species_E1,sep="_"),dis)		
	# dis3=dis2[!dis2$ID%in%intro.nz$ID,]
	# an=as.data.frame(fread("GyAn.csv"))	
	# spdis=subset(dis3,dis3$Species_E1%in%an$Genus)#只保留被子植物
	# geo=read.csv("Geo-isrm.csv")#在Geonames.csv(来自county——Albert.shp)的基础上remove the islands(area<2.5 wan km2)，以及所罗门群岛
	# spdis2=subset(spdis,spdis$Adcode99%in%geo$ADCODE99)#去除了岛屿	
	# write.csv(spdis2,"SpDistribution20210106.csv")#目前最新全球物种分布
	
##制作分布矩阵
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

		return(spdat)#全球属的分布数据		
	}
	
	library(data.table)
	spdis2=as.data.frame(fread("SpDistribution20191016.csv"))[,-1]
	spdat.glb=spdat(spdis2)	
	write.csv(spdat.glb,"Spdat_an_isrm191016.csv")##目前最新全球属级分布，14016 genus

######################################################################################################################
######################################################################################################################
####### calculate pairwise betasim distances######	
	#找出物种的所有祖先node所在编号
	ancestor <- function(tree, tip, include.tip = FALSE)
				{
			if (class(tree) != "phylo") stop("The tree is not a phylo tree")
			N.ii <- 1:Ntip(tree)
			N.node <- Nnode(tree)
			
			if (class(tip) == "character") tip <- joint.table1(list(tip=tip),"tip", list(N.ii=N.ii, tip=tree$tip.label),"tip")[,2]
			
			## find the root
			root.label <- unique(tree$edge[,1][!tree$edge[,1] %in% tree$edge[,2]])
			
			tree.edge <- rbind(tree$edge,c(0,root.label))
			ii <- order(tree.edge[,2])
			tree.edge <- tree.edge[ii,]

			anc  <- vector(mode = "list", length = length(tip))
			names(anc) <- tree$tip.label[tip]
			for (i in 1:length(tip))
				{
				descendant <- tip[i]
				anc.t <- numeric(0)
				if (include.tip) anc.t <- descendant
				for (j in 1:N.node)
					{
					if (tree.edge[descendant,1] == 0)
						{
						break
						}
					if (tree.edge[descendant,1] == root.label) 
						{
						anc.t <- c(anc.t, root.label)
						break
						}
					ancestor <- tree.edge[descendant,1]
					descendant <- tree.edge[ancestor,2]
					anc.t <- c(anc.t, ancestor)
					}
				anc[[i]] <- as.numeric(rev(anc.t))
				}
			return(anc)
			}
## to joint table by index1 of data table 1 and index2 of data table 2
		## By Zhiheng Wang
		## join two tables by common rows according to index1 and index2
		joint.table1 <- function(data1, index1, data2, index2, case = TRUE, ...)
			{
			merge.cols <- function(newdata, m,n,i,data2,jj)
				{
				var.class <- class(data2[[select[m]]])

				newdata[[m]] <<- new(var.class)
				length(newdata[[m]]) <<- n

				if (var.class == "factor")
					{			
					f.temp <- as.character(data2[[select[m]]])[jj]
					levels(newdata[[m]]) <<- levels(f.temp)
					newdata[[m]][i[ii]] <<- as.factor(f.temp)#(as.character(data2[[select[m]]][jj]))
					}
				else
					{
					newdata[[m]][i[ii]] <<- data2[[select[m]]][jj]
					}
				}

			if (!is.data.frame(data1))
				data1 <- as.data.frame(data1,stringsAsFactors=F)
			if (!is.data.frame(data2))
				data2 <- as.data.frame(data2,stringsAsFactors=F)	
			
			if (case)
				index1 <- as.character(data1[[index1]])
			else
				index1 <- tolower(as.character(data1[[index1]]))

			n <- length(index1)
			i <- 1:n
			if (length(unique(index1))!=n) stop("duplicates are found in index1")
			
			if (class(index2)=="numeric")
				{
				select <- c(1:dim(data2)[2])[-index2]
				select.name <- names(data2)[-index2]
				}
			else if (class(index2)=="character")
				{
				select <- names(data2)[!names(data2) %in% index2]
				select.name <- select
				}
			else
				{stop("error1")}
			if (case)
				index2 <- as.character(data2[[index2]])
			else
				index2 <- tolower(as.character(data2[[index2]]))
			
			cond <- duplicated(index2)
			if (any(cond,na.rm=T)) {
				index2 <- index2[!cond]
				data2 <- data2[!cond,]
				}
			
			cond <- index2 %in% index1
			index2 <- index2[cond]
			data2 <- as.data.frame(data2[cond,select],stringsAsFactors=F)
			names(data2) <- select.name
			rm(cond)
			
			cond <- index1 %in% index2
			index1 <- index1[cond]
			i <- i[cond]
				
			ii <- order(index1)
			jj <- order(index2)
			if (all(index1[ii]==index2[jj]))
				{
				newdata <- vector(mode="list",length=length(select))
				names(newdata) <- names(data2)
				for (m in 1:length(select)) merge.cols(newdata, m, n, i, data2, jj)
				data1 <- cbind(data1, as.data.frame(newdata, stringsAsFactors=F))
				return(data1)
				}
			else
				{stop("no shared row in the two tables")}

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
	#save(pBeta.mat, file="processed_aV2/psim.ansy.Rdata")	#both ang and gymnosperms
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
			#psim1=get(load("processed_aV2/psim.ansy.Rdata"))#matrix			
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
		par(mfrow = c(3,1),mar=c(0.5,4,0.5,1),oma=c(2,2,2,2))
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
		dendcolor=color_b;labels(dendcolor)=rep(NA,420)
		plot(dendcolor,cex.axis=1.5)				
		plot(x, y, col = as.character(col_NMDS$leaf_col),pch=19,xlab="NMDS1",ylab="NMDS2",cex.axis=1.5)		
		abline(h=0,v=0,lty=2,col="lightgray")				
	}
	xy <- read.csv("Geo-isrm.csv",head=T)[,c("ADCODE99","Lon","Lat")]
	xy=cbind(no=1:420,xy)	
	#set color for phylo maps
	#kingdom
	mycolor.kind=c("#93AA00","#D39200","#DB72FB","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D")
	mycolor.kind.sp=c("#93AA00","#D39200","#F8766D","#619CFF","#00C19F","#00BA38","#DB72FB","#FF61C3")
	colorname.kind=data.frame(color=mycolor.kind,reg=c("Saharo-Arabian","Holarctic","Chlie-Patagonian","Neotropical","Indo-Malesian","African","Australian","Neozeylandic"))	
	mycolor=c("#93AA00","#D39200","#DB72FB","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D")
	for (n in 1:length(mat)){		
		windows()
		plot.map(xy,nclust=8,mycolor.kind.sp,mat[n])
	}
	
	#region
	mycolor.region=c("#95AB74","#95AB74","#95AB74","#CDC78A",#Saharo，将1和2合并到3
	"#FDDE87","#FCDFB3","#D39200",#Holarctic
	"#619CFF", "#619CFF",#Patagonia,合并8和9
	"#E45B41","#F6BBAD",#Neotropical
	"#00BA38","#B8D25F",#Indo-Maly
	"#0DB09D","#79C6B7","#A6D5C3",#African
	"#DB72FB",#Australian
	"#E2318A","#F8BDD1"#NewZealand
	)
	colorname.region=data.frame(color=mycolor.region,reg=c(1,1,1:6,6:16))	
	for (n in 1:length(mat)){
		windows()	
		plot.map(xy,nclust=length(mycolor.region),mycolor.region,colorname.region,mat[n])
	}

##compare with carta et al. 2021
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
	
	#Plot Adcode map
	library(sp)
	library(maptools)
	shape<-readShapeSpatial("county_islandrm.shp")
	coor=cbind(shape@data$Lon,shape@data$Lat)
	rownames(coor)=shape@data$ADCODE99
	windows()
	plot(shape,col="white", border = "gray")
	text(coor, labels = rownames(coor), cex =0.8)	
	
	#Plot dend,nmds and map	
	library(sp)
	library(maptools)
	library(dendextend)
	library(recluster)
	library(vegan)
	xy <- read.csv("Geo-isrm.csv",head=T)[,c("ADCODE99","Lon","Lat")]
	xy=cbind(no=1:420,xy)
				
	psim=get(load("processed_aV2/psim_Carta2021.Rdata"))
	#psim=get(load("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata"))
	
	psim[is.na(psim)]=0
	psim.t=psim+t(psim)
	psim2=as.dist(psim.t)		
	hc <- hclust(psim2, method="average")
	dend <- as.dendrogram(hc)
		#recl <- recluster.region(psim2, method="average", mincl=8, maxcl=32)
		#nclust=min(recl[[1]][recl[[1]][,"ex.diss"]>=0.8,"k"])
	mycolor=rainbow(nclust)
	mycolor=c("#93AA00","#D39200","#DB72FB","#FF61C3","#00BA38","#00C19F","#619CFF","#F8766D")
	color_b <- color_branches(dend,k=nclust,col=mycolor)		
		## map UPGMA clustering					
	leaf_col <- leaf_Colors(color_b)
	order <- order.dendrogram(color_b)
	d_leaf_col <- data.frame(order, leaf_col)
	match <- data.frame(row.names(d_leaf_col), leaf_col)
	color <- merge(xy[,1:4], match, by.x="ADCODE99", by.y="row.names.d_leaf_col.", all=F)
	shape<-readShapeSpatial("county_islandrm.shp")			
	shape@data=cbind(shape@data,phcol=color[match(shape@data$ADCODE99,color$ADCODE99),"leaf_col"])
		#shape@data=cbind(shape@data,phname=colorname[match(shape@data$phcol,colorname$color),"reg"])
		#write.csv(re,"PhyregionForGIS.csv")
		# run NMDS		
	nmds <- monoMDS(psim2)		
	x <- nmds$points[,1]
	y <- nmds$points[,2]
	d <- data.frame (row.names(as.matrix(psim2)),x,y)
	col_NMDS <- merge(d,color, by.x="row.names.as.matrix.psim2..", by.y="ADCODE99",sort=F)
	windows()
	par(mfrow = c(3,1),mar=c(0.5,4,0.5,1),oma=c(2,2,2,2))
	plot(shape, col="gray", border = "white")				
	for (i in 1:nclust) {
	shp=subset(shape,shape@data$phcol==mycolor[i])		
	plot(shp,col=mycolor[i],border=mycolor[i],add=T)
	}
			# lon=tapply(X=shape@data$Lon,INDEX=shape@data$phname,FUN=mean)
			# lat=tapply(X=shape@data$Lat,INDEX=shape@data$phname,FUN=mean)
			# coor=cbind(lon,lat)
			# text(coor, labels = rownames(coor), cex =1.2)		
	box()
	
	dendcolor=color_b;labels(dendcolor)=rep(NA,420)
	plot(dendcolor)	
	plot(x, y, col = as.character(col_NMDS$leaf_col),pch=19,xlab="NMDS1",ylab="NMDS2")
	abline(h=0,v=0,lty=2,col="lightgray")	
	
## 比较 The maximum clade credibility tree 和the individual trees drawn from the posterior distribution
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

###performance of other clustering method
	library(ggpubr)	
	prefo=read.csv("Appendix 2-perfomace of methods.csv")
	cpindex=na.omit(prefo[,c(33:37)])	
	method=c("single","complete","median","mcquitty","centroid","ward.D2","average")
		
	perfrom.plot=function(expvar.n,cpindex.n,method,show.xlab=TRUE){
		title=colnames(cpindex.n)[2]
		colnames(cpindex.n)=c("cluster","index")
		expvar=c();exp.thes=c()
		for (i in 1:length(method)){
			tmp=expvar.n[,c(1,1+i)]
			colnames(tmp)=c("cluster","exp")		
			exp.thes.t=ifelse(max(tmp$exp>99),tmp[tmp$exp>99,][1,1],NA)
			exp.thes=c(exp.thes,exp.thes.t)
			expvar=rbind(expvar,data.frame(method=method[i],tmp))
		}
		expvar$method=factor(expvar$method,levels=method)
		cpindex.n$cluster=factor(cpindex.n$cluster,levels=method)
		
		main=ggplot(expvar, mapping = aes(x = cluster, y = exp,color=method)) + geom_line(size=2)+ theme_bw()+
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
				main=main+geom_vline(xintercept=exp.thes[i],col=hue_pal()(7)[i],size=1,linetype="longdash")
			}
		}
		if (show.xlab) {
			main=main+theme(axis.text.x  =element_text(size=20))
		}	else {
			main=main+theme(axis.text.x  =element_blank())
		}
		
		sub <- ggplot(cpindex.n, mapping = aes(x = cluster, y = index,color=cluster)) + geom_point(size=5,show.legend=FALSE)+ theme_bw()+
		theme(axis.title.x = element_blank(),
			axis.title.y = element_text(size=20),
			axis.text.x  = element_text(size=20,angle =30,vjust=1, hjust = 1),
			axis.text.y  = element_text(size=20))+ylab("Co-phenetic correlation")#+ggtitle(title)
		sub$layers <- rev(sub$layers)
		p=main + annotation_custom(ggplotGrob(sub), xmin=50, ymin=5, xmax=220, ymax=65) 
		return(p)
	}
	
	p=list()
	for (i in 1:4){
		expvar.n=prefo[,c((8*(i-1)+1):(8*(i-1)+8))]
		cpindex.n=cpindex[,c(1,(i+1))]
		if (i<=2) {p[[i]]=perfrom.plot(expvar.n,cpindex.n,method,FALSE)} else {p[[i]]=perfrom.plot(expvar.n,cpindex.n,method)}
	}	
	ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],nrow=2,ncol = 2,widths=c(1,1),heights=c(1,1),
		labels=c("a","b","c","d"),label.x=0.1,label.y=1,font.label = list(size = 30),legend="right",common.legend=TRUE)
		
##计算sil_witdth并确定uncertain region	
#s(i) = (b(i) - a(i)) / max(a(i), b(i)),对于单元格i，a(i)为该单元格与所属类别内所有其他单元格的平均距离；             先获得该单元格与其他cluster内所有单元格平均距离最近的那个cluster，然后b(i)为单元格i与该cluster内所有单元格的平均距离 	
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
	
	# mat=c()
	# for (i in 1:3) mat=c(mat,paste("raxml_v2_dated_Av",i,"_final_clean",sep=""))	
	# for (i in 1:5) {		
	# for (j in 1:3) mat=c(mat,paste("MCC_aV",j,"_",i,"_proc",sep=""),paste("tree",i,".aV",j,".dated",sep=""))	
	# }
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
	
	#基于average方法计算sil_witdth	
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
	
	
## overlay species range of each kingdom
	library(data.table)
	library(raster)
	library(ggpubr)
	genus.dis=as.data.frame(fread("SpDistribution20191016.csv"))[,-1]
	mycolor.kind=c("#93AA00","#D39200","darkblue","#FF61C3","#00C19F","#00BA38","#619CFF","#F8766D")
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
		splist=unique(genus.dis[genus.dis$Adcode99%in%adc,"Species_E1"])
		dis=subset(genus.dis,genus.dis$Species_E1%in%splist)
		sprich.cluster=tapply(dis$Species_E1,dis$Adcode99,length)
		sprich.all=tapply(genus.dis$Species_E1,genus.dis$Adcode99,length)
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
				scale_fill_gradient2(midpoint = midpoint,low="lightgray",mid="lightgray",high=color)+		
				theme_minimal()+	
				theme(axis.title = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
				panel.grid=element_blank(),plot.title = element_blank(),legend.position="none")
		return(p.map)
	}
	p=list()	
	for (i in 1:8){	
		adc=names(cluster[cluster==i])
		color=mycolor.kind[i]		
		#pdf(paste0("processed_aV2/cluster.rich_", i, ".pdf",sep=""), width = 10, height = 6)
		p[[i]]=plot.cluster.rich(adc,genus.dis,color,cordacd,data.test,70)		
		#dev.off()
	}
	ggarrange(p[[1]],p[[2]],p[[7]],p[[8]],p[[6]],p[[5]],p[[3]],p[[4]],
	nrow=4,ncol =2,widths=c(1,1,1,1),heights=c(1,1),labels=c("a","b","c","d","e","f","g","h"),label.x=0.03,font.label = list(size = 20))

##endemic taxa in each phylo region
	#1.combine phylo region,phylo kingdom and adcode
	result=function(re){
	 relist=list()
	 temp=c()
	 class=c()
	 for (i in 1:length(re)){
	   data=as.matrix(re[[i]])
	   relist[[i]]=as.numeric(rownames(data))
	   temp=cbind(adcode=relist[[i]],i)
	   class=rbind(class,temp)	   
	 }
	 rownames(class)=class[,1]
	return(class) 
	}
	
	psim1=get(load("processed_aV2/psim.raxml_v2_dated_Av2_final_clean_ann.full.polyRes.Rdata"))#matrix
	psim1[is.na(psim1)]=0
	psim2=as.dist(psim1+t(psim1))
	hc <- hclust(psim2, method="average")		
	plot(hc)
	re.sim=rect.hclust(hc,k=8)
	resim=result(re.sim)
	re.sim2=rect.hclust(hc,k=19)
	resim2=result(re.sim2)	
	phname=read.csv("cluster.average.merge.csv")[,-1]#names of each phylo kingdom
	phname2=cbind(adcode=resim[,1],phy.region=resim2[match(rownames(resim),rownames(resim2)),2],
		phname[match(rownames(resim),phname$adcode_ave2),c("adcode_ave2","phname")])
	#2.match phylo tree, distribution data to get the genus in each adcode
	taxonmic=read.csv("Gyo-angiosperm.csv")	
	library(data.table)
	dis <- as.data.frame(fread("Spdat_an_isrm191016.csv",head = T))
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]		
	library(picante)
	tre=read.tree("raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes.tre")#140-150 constraint
	phydata <- match.phylo.comm(tre,spdis)	
	occur=unique(as.data.frame(data.table::fread("SpDistribution20191016.csv")))[,-1]
	occur.angio=na.omit(occur[match(occur$Species_E1,colnames(phydata$comm))&match(occur$Adcode99,rownames(phydata$comm)),])
	#3. match genus in each adcode with phylo region,phylo kingdom
	splist=cbind(taxonmic[match(occur.angio$Species_E1,taxonmic$Genus),],phname2[match(occur.angio$Adcode99,phname2$adcode),-3])
	colnames(splist)[6]="phy.kind"
	write.csv(splist,"Genus.dis.csv")#目前最新全球属级分布，属名匹配了phylo tree，分布信息匹配了谱系地理区	
	#4. get endemic genus in each phylo kingdom and phylo region	
	endemic.in.phyregion=function(tx,thes,reg,splist){
		lst=as.character(unique(splist[,tx]))
		endemic=na.omit(do.call(rbind,lapply(1:length(lst),function(i){
			splist.t=splist[splist[,tx]==lst[i],]
			tmp1=tapply(splist.t$adcode,splist.t[,reg],length)
			tmp1=tmp1[!is.na(tmp1)]
			tmp2=sort(tmp1/sum(tmp1),decreasing=T)#每个物种在各个生物地理区占该物种总分布区的比例
			name=ifelse(tmp2[1]>=thes,names(tmp2[1]),NA)
			re=data.frame(region=name,taxa=lst[i])
			return(re)
		}
		)))	
		colnames(endemic)=c(reg,tx)
		rownames(endemic)=NULL
		return(endemic)
	}	
	splist=read.csv("Genus.dis.csv")[,-1]
	#kingdom
	Kingdom=as.character(unique(splist[,"phy.kind"]))
	##order
	kingdom.ord=endemic.in.phyregion("Order",0.8,"phy.kind",splist)	
	kingdom.ord.stat=tapply(kingdom.ord$Order,kingdom.ord$phy.kind,length)[Kingdom]
	names(kingdom.ord.stat)=Kingdom
	write.csv(kingdom.ord,"kingdom.ord.csv")
	##family
	kingdom.fam=endemic.in.phyregion("Family",0.8,"phy.kind",splist)	
	kingdom.fam.stat=tapply(kingdom.fam$Family,kingdom.fam$phy.kind,length)[Kingdom]
	names(kingdom.fam.stat)=Kingdom
	write.csv(kingdom.fam,"kingdom.fam.csv")	
	##genus
	kingdom.genus=endemic.in.phyregion("Genus",0.8,"phy.kind",splist)	
	kingdom.genus.stat=tapply(kingdom.genus$Genus,kingdom.genus$phy.kind,length)[Kingdom]
	names(kingdom.genus.stat)=Kingdom
	write.csv(kingdom.genus,"kingdom.genus.csv")
	kingdom.stat=data.frame(Order=kingdom.ord.stat[Kingdom],Family=kingdom.fam.stat[Kingdom],Genus=kingdom.genus.stat[Kingdom])
	kingdom.stat[is.na(kingdom.stat)]=0
	write.csv(kingdom.stat,"kingdom.stat.csv")
	##Neogenus in Hol
	geo=read.csv("Geo-isrm.csv")
	fcm.dis=cbind(as.data.frame(fcm.value),continent=geo[match(rownames(fcm.value),geo$ADCODE99),"continent"],
		phname=phname[match(rownames(fcm.value),phname$adcode_ave2),"phname"])
	std_error=function(x) sd(x)/sqrt(length(x))
	sil.NAma=fcm.dis[fcm.dis$continent=="NAmerica"&fcm.dis$cluster==3,"sil_width"]
	mean(sil.NAma);std_error(sil.NAma)
	sil.pata=fcm.dis[fcm.dis$phname=="Chlie-Patagonian","sil_width"]
	mean(sil.pata);std_error(sil.pata)
	mean(fcm.dis$sil_width);std_error(fcm.dis$sil_width)
	
	genus.dis=cbind(splist,continent=geo[match(splist$adcode,geo$ADCODE99),"continent"])
	neo.genus = kingdom.genus[kingdom.genus$phy.kind=="Neotropical","Genus"]
	hol.genus = kingdom.genus[kingdom.genus$phy.kind=="Holarctic","Genus"]
	sahro.genus=kingdom.genus[kingdom.genus$phy.kind=="Saharo-Arabian","Genus"]
	pata.genus = kingdom.genus[kingdom.genus$phy.kind=="Chlie-Patagonian","Genus"]
	eur.genus = unique(genus.dis[genus.dis$continent=="Europe","Genus"])
	NAma.genus.dis = unique(genus.dis[genus.dis$continent=="NAmerica"&genus.dis$phy.kind=="Holarctic",])
	NAma.genus0=sort(tapply(NAma.genus.dis$adcode,NAma.genus.dis$Genus,length),decreasing=T)
	NAma.genus=names(NAma.genus0[!(names(NAma.genus0)%in%neo.genus&NAma.genus0<=1)])
	
	Pata.geuns.dis=unique(genus.dis[genus.dis$phy.kind=="Chlie-Patagonian",])
	Pata.genus0=sort(tapply(Pata.geuns.dis$adcode,Pata.geuns.dis$Genus,length),decreasing=T)
	Pata.genus=names(Pata.genus0[!(names(Pata.genus0)%in%neo.genus&Pata.genus0<=1)])
	
	length(eur.genus)
	length(eur.genus[eur.genus%in%sahro.genus])
	length(eur.genus[eur.genus%in%hol.genus])
	
	length(NAma.genus)
	length(NAma.genus[NAma.genus%in%neo.genus])
	length(NAma.genus[NAma.genus%in%hol.genus])
	
	length(Pata.genus)
	length(Pata.genus[Pata.genus%in%neo.genus])
	length(Pata.genus[Pata.genus%in%pata.genus])
	#region
	Region=as.character(unique(splist[,"phy.region"]))
	##order
	region.ord=endemic.in.phyregion("Order",0.8,"phy.region",splist)	
	region.ord.stat=tapply(region.ord$Order,region.ord$phy.region,length)[Region]
	names(region.ord.stat)=Region
	write.csv(region.ord,"region.ord.csv")
	##family
	region.fam=endemic.in.phyregion("Family",0.8,"phy.region",splist)	
	region.fam.stat=tapply(region.fam$Family,region.fam$phy.region,length)[Region]
	names(region.fam.stat)=Region
	write.csv(region.fam,"region.fam.csv")	
	##genus
	region.genus=endemic.in.phyregion("Genus",0.8,"phy.region",splist)	
	region.genus.stat=tapply(region.genus$Genus,region.genus$phy.region,length)[Region]
	names(region.genus.stat)=Region	
	write.csv(region.genus,"region.genus.csv")
	
	region.stat=data.frame(Order=region.ord.stat[Region],Family=region.fam.stat[Region],Genus=region.genus.stat[Region])
	region.stat[is.na(region.stat)]=0
	write.csv(region.stat,"region.stat.csv")	
	
	##clades in the regions which belong to different kingdoms based on BD and TBD
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

####### Node-based analysis,REF: Node-based analysis of species distributions######
# noderich=subtrees(tre,wait=TRUE)#非常耗时			
library(nodiv)
library(ape)
library(sp)
library(maptools)

	dis=read.csv("Spdat_an_isrm191016.csv")	
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]		
	#tre_ori=read.tree("raxml_v2_dated_Av1_final_clean_ann.full.polyRes.tre")#全球属级的树,140-210 constraint
	tre_ori=read.tree("raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes.tre")#全球属级的树,140-150 constraint
	#tre_ori=read.tree("raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes.tre")#全球属级的树,140-259 constraint
	d=subset(tre_ori$tip.label,!tre_ori$tip.label%in%colnames(spdis))
	tre <- drop.tip(tre_ori, tip=d)
	#plot(tre,show.tip.label=T,show.node.label=T)#显示谱系树，包含物种和科属的名字
	#tre$node.label	
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
	
##节点编号与节点科属名以及GND值的对应关系
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

	#匹配GND与node age
	nodebase=vector("list",length(region))#GND>0.65的node的编号、上级编号、年龄和GND值
	for (i in 1:length(region)){
		tre.t=ccoquettes[[i]]$phylo   ####8个高级节点的子树
		root.label <- unique(tre.t$edge[,1][!tre.t$edge[,1] %in% tre.t$edge[,2]])			
		tre.edge <- rbind(tre.t$edge,c(0,root.label))
		ii <- order(tre.edge[,2])
		tre.edge <- tre.edge[ii,]
		pos2=which(tre.edge[,2]>Ntip(tre.t))#提取edge中每个物种对应的编号所在的行号
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
#检验node对bioregion的贡献##R2	
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
		node=colnames(nodelist.t)#GND node节点编号
		
		r2[[i]]=matrix(0,length(node),3)#节点与边界的矩阵，R2
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
	
	#箱线图
	richfin.box=lapply(1:length(richfin),function(j){
		box.data=do.call(rbind,lapply(1:length(richfin[[j]][,"age"]),function(i){	
		  cbind(richfin[[j]][i,],boxg=floor(richfin[[j]][i,"age"]/20)*20)	  
		}))
	})
	
	plotn=vector("list",length(richfin))	
	for (i in 1:length(richfin)){
		p=ggplot(richfin.box[[i]],aes(x=age,y=r2,group=boxg))+stat_boxplot(geom="errorbar",width=0.15,color="#39B600")+geom_boxplot()		
		p=p+ theme(panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=10),
			axis.text.y  = element_text(size=10))+
			ylim(0,100)+xlim(150,0)		
				
		plotn[[i]]=p		
	}	
	
	figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], 
			nrow=3,ncol = 3,widths=c(2,2,2),heights=c(2,2,2),common.legend=TRUE,legend = "top",hjust=0,align="hv")
	
	richfin=readRDS("richfin150.rds")
	#划图	
	library(ggplot2)	
	library(gridExtra)		
	library(cowplot)
	library(ggpubr)	
	noder=c()	
	for (i in 1:length(region)){
		tmp=na.omit(richfin[[i]])		
		tmp=tmp[tmp$r2>0&tmp$p<0.05,]
		tmp2=cbind(tmp,region=region[i])		
		noder=rbind(noder,tmp2)		
		}	
		
	noderp=noder	
	#画GND节点在谱系树的分布图	
	# Create two panels side by side
	layout(matrix(1:9,3,3), widths=rep(1))
	# Set margins and turn all axis labels horizontally (with `las=1`)
	par(mar=rep(.5, 4), oma=rep(3, 4))	
	for (i in 1:length(region)){
		gnd.t=noderp[noderp$region==region[i],c("edge2","GND","r2")]
		gnd.t=gnd.t[order(-gnd.t$GND),]
		tre=ccoquettes[[i]]$phylo		
		plot(tre,show.tip.label=FALSE,edge.width=1.5,no.margin=FALSE)#,type="fan")
		#axisPhylo()	
		for (j in 1:dim(gnd.t)[1]){
			nodelabels(node =as.integer(gnd.t[j,1]), pie = 1,piecol=heat.colors(dim(gnd.t)[1])[j], cex =gnd.t[j,2])	#颜色为GND值，大小为r2
		}
		text(10,50, region[i],cex=2,col="red")		
	}	
	
	##散点图
	library(scales)
	show_col(hue_pal()(30)) #查看前9个默认的颜色梯度
	#mycolor=c("#F8766D","#D39200","#93AA00","#00BA38","#00C19F","#00B9E3","#619CFF","#DB72FB","#FF61C3" )
	leg=c("a","b","c","d","e","f","g")
	plotn=vector("list",length(region))	
	for (i in 1:length(region)){
		p=ggplot(data = noderp[noderp$region==region[i]&noderp$age<=100,],mapping=aes(x=age,y=r2)) 		
		if(i==1) p=p+geom_point(size=2,colour="gray",fill="gray",shape=21)				
		if(i>1)	p=p+geom_point(size=2,colour="gray",fill="gray",shape=21,show.legend=FALSE)		
		p=p+stat_density2d(aes(alpha=..density..), geom="tile",contour=FALSE,show.legend=FALSE)+		
		theme(panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=12),
			axis.text.y  = element_text(size=12))+			
			#scale_x_continuous(limits=c(0,100))+
			scale_y_continuous(expand = c(0,0),breaks=pretty_breaks(5),limits=c(0,100))
			annotate("text",x=150,y=100,label=leg[i],size=8)		
		d=glm(r2~age,data = noderp[noderp$region==region[i],])	
		if(ss.glm(d)[3]<0.05) p=p+geom_smooth(method = "glm",colour="red",fill="#619CFF",size=1)
		if(ss.glm(d)[3]<0.1&ss.glm(d)[3]>0.05) p=p+geom_smooth(method = "glm",colour="red",,fill="#619CFF",linetype="dashed")#marginal sig		
		plotn[[i]]=p		
	}	
	
	figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], 
			nrow=2,ncol = 4,widths=c(2,2),heights=c(2,2,2,2),common.legend=TRUE,legend = "bottom",hjust=0,align="hv",labels = c("a", "b","c","d", "e","f","g"),label.x=0.12)
	annotate_figure(figure, bottom = text_grob("Clade crown age (Ma)", size=14),left = text_grob("Clade contribution on floristic division", rot = 90,size=14))
	
	 #South的极端值
	 install.packages("BiocManager")
	 BiocManager::install("ggtree")
	library(ggtree)
	
	 library(nodiv)
	 res=readRDS("res.rds")	
	 summary(res[[2]])	 
	 tre=res[[2]]$phylo
	 node=subset(tre$edge[,2],tre$edge[,1]==9628)
	 
	 tree <- groupClade(tre, node=c(21, 17))
	ggtree(tree, aes(color=group, linetype=group))

	 plotSOS(res[[2]],9628)	 
	 for (k in 1:length(noderich2[[2]])){
		if (noderich2[[2]][[k]]$node.label[1]==12273)	{	
			tre.t0=noderich2[[i]][[k]]						
			d.t=tre.t0$tip.label					
			}
	 }
	zoom(tre,d.t)	
	text="((c:25,((d:8,e:8):7,f:15)G2:10)F1:25,(a:10,b:10,g:10)G1:40)O1:50;"#括号：共同祖先；冒号加数字：枝长；小写字母Tip/物种名;大写字母：科属名；分号：一棵树的结尾；逗号：姊妹种
	tre=read.tree(text=text)#读取谱系数据
	plot(tre, show.tip.label=F), 
     edge.color = ifelse(tre$edge[,1]<=7 %in% c("a","b","c"),
                      'red','blue'))

##提取高GND值节点对应的物种列表
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
	
#####Fig.2 谱系和环境的关系#####	
	ss.glm <- function(r.glm){
		r.ss <- summary(r.glm)
		rsq <- 100*(r.ss$null.deviance-r.ss$deviance)/r.ss$null.deviance
		adj.rsq <- 100*(1-(r.ss$deviance/r.ss$df.residual)/(r.ss$null.deviance/r.ss$df.null))
		f.stat <- ((r.ss$null.deviance-r.ss$deviance)/(r.ss$df.null-
		r.ss$df.residual))/(r.ss$deviance/r.ss$df.residual)
		p <- pf(f.stat, r.ss$df.null-r.ss$df.residual, r.ss$df.residual, lower.tail=FALSE)
		return(c(r2=rsq,adj.r2=adj.rsq,p=p))
		}
		
	
##气候对物种和谱系影响的比较
	library(ggplot2)	
	library(ggpubr)		
	library(gdm)
	psim1=get(load("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata"))#matrix
	psim1[is.na(psim1)]=0
	phylo=psim1+t(psim1)
	beta.sim=get(load("processed_aV2/psim.sp.Rdata"))
	sp=betamatrix(beta.sim,as.dist=FALSE)
	clim=read.csv("richclim191024.csv")[,c("adcode","Lon","Lat","MAT","MAP")]	
	regionlist=read.csv("cluster.average.sphy.csv")
	
	node=colnames(regionlist)[8:9]
	node.sp=colnames(regionlist)[10:11]
	exFormat=vector("list",(length(node)*2));names(exFormat)=c(paste("ph_",node,sep=""),paste("sp_",node,sep=""))
	glm.r=c()
	for (i in 1:(length(node)*2)){
	##将beta分成between和within region，输出gdm所需的exFormat表
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
			envir=cbind(clim,regionlist[match(clim$adcode,regionlist$adcode_ave2),c("phname",node.sp[i-length(node)])])				
			envt=subset(envir,envir[,node.sp[i-length(node)]]>0)			
			spsim=as.data.frame(sp[match(envt$adcode,colnames(sp)),match(envt$adcode,rownames(sp))])
			Dissim=cbind(adcode=envt$adcode,spsim)
			exFormat_pre <- formatsitepair(Dissim,3, XColumn="Lon", YColumn="Lat", predData=envt,siteColumn="adcode")	
			cole=6+(length(colnames(exFormat_pre))-6)/2	
			cole2=length(colnames(exFormat_pre))
			exFormat[[i]]=subset(exFormat_pre,
				  (exFormat_pre[,cole]==1 & exFormat_pre[,cole2]==2)|(exFormat_pre[,cole]==2 & exFormat_pre[,cole2]==1))
		}		
		
		glm.data=data.frame(beta=exFormat[[i]]$distance,endis=sqrt((exFormat[[i]]$s1.MAT-exFormat[[i]]$s2.MAT)^2+(exFormat[[i]]$s1.MAP-exFormat[[i]]$s2.MAP)^2))
		r.squ <- data.frame(region=names(exFormat)[i],r2=ss.glm(glm(beta~endis,glm.data,family=binomial))[1])
		glm.r=rbind(glm.r,r.squ)	
	}
	glm.r[is.na(glm.r[,2]),2]=0
	
	hp.plot=vector("list",length(node));names(hp.plot)=node
	plotn=vector("list",length(node))
	hrtc=data.frame(node=sort(rep(1:length(node),2)),region=c("ph_Neotropic","sp_Neotropic",
				"ph_Hoarcortal","sp_Hoarcortal"))
				
	for (i in 1:length(node)){
		hrtc.t=subset(hrtc,hrtc$node==i)
		hp.plot[[i]]=subset(glm.r,glm.r$region%in%hrtc.t$region)		
		p=ggplot() +
		geom_bar(data=hp.plot[[i]], aes(x=region,y=r2),fill=c("#F8766D","#619CFF"), stat="identity",
			   position='stack')+			  
			   theme_bw() + 			
			theme(axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=10),
			axis.text.y  = element_text(size=10))+
			theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
			  panel.background = element_blank(),axis.line = element_line(colour = "black"))
		p=p+ scale_x_discrete(labels = as.factor(c("PBD","BD"))) 		
		p=p+ylim(0,20)
		plotn[[i]]=p
	}	
	
	figure=ggarrange(plotn[[1]],plotn[[2]],nrow=1,ncol = 2,widths=c(2,2),heights=2,common.legend=TRUE,legend = "top",label.x=0.5,hjust=0,align="hv")	 	
	annotate_figure(figure, bottom = text_grob("Climate Euclidean distance", size=14),left = text_grob("R2 (%)", rot = 90,size=14))	
	
##Fig.2 谱系和环境的关系
##现代气候对分区边界的影响 GLM
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
	##将beta分成between和within region，输出gdm所需的exFormat表
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

	##version 2：气候对物种和谱系影响的比较Mantel tests warning: non-square matrix
	library(vegan)	
	library(ggplot2)	
	psim1=get(load("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata"))#matrix
	psim1[is.na(psim1)]=0
	phylo=psim1+t(psim1)	
	clim=read.csv("richclim191024.csv")[,c("adcode","Lon","Lat","MAT","MAP")]	
	regionlist=read.csv("cluster.average.sphy.csv")
	phname=c("Saharo-Arabian","Holarctic","Patagonia","Neotropical","African","Oriental","Australian","NewZealand")
	
	env <- clim[ ,c("MAT","MAP")]
	scale.env <- scale(env, center = TRUE, scale = TRUE)
	rownames(scale.env)=clim$adcode
	dist.env <- as.matrix(dist(scale.env, method = 'euclidean'))
	re=c()
	for (i in 1:length(phname)){
		tmp2=c()
		for (j in 1:i){
		adc.i=regionlist[regionlist$phname%in%phname[i],"adcode_ave2"]
		adc.j=regionlist[regionlist$phname%in%phname[j],"adcode_ave2"]
		psim.t=phylo[match(adc.i,rownames(phylo)),match(adc.j,colnames(phylo))]
		env.t=dist.env[match(adc.i,rownames(dist.env)),match(adc.j,colnames(dist.env))]		 
		if (i==j) {	
			if(length(as.dist(psim.t))==1) {tmp=data.frame(x=phname[i],y=phname[j],cor=0,p=0,cor.l=0,cor.h=0);tmp2=rbind(tmp2,tmp);next;} 		
			cor=cor(as.dist(psim.t),as.dist(env.t))
			cor.test=cor.test(as.dist(psim.t),as.dist(env.t))
		} else {
			cor=cor(as.numeric(psim.t),as.numeric(env.t))
			cor.test=cor.test(as.numeric(psim.t),as.numeric(env.t))
		}
		tmp=data.frame(x=phname[i],y=phname[j],cor=cor,p=cor.test$p.value,cor.l=cor.test$conf.int[1],cor.h=cor.test$conf.int[2])
		tmp2=rbind(tmp2,tmp)
		}
	re=rbind(re,tmp2)
	}
	
	#bubble plot	
	re$x=factor(re$x,levels=c("Saharo-Arabian","Holarctic","Patagonia","Neotropical","African","Oriental","Australian","NewZealand"),ordered=TRUE)
	re$y=factor(re$y,levels=c("Saharo-Arabian","Holarctic","Patagonia","Neotropical","African","Oriental","Australian","NewZealand"),ordered=TRUE)	
	ggplot(re, aes(x=x,y=y,fill=cor)) + geom_point(shape=21,stroke =0.5,size=20)+
			geom_abline(intercept=0,slope=1,colour = "darkgray", linetype = "twodash", size = 1)+			
			scale_fill_gradient2(midpoint = 0,low="#c15dd5", mid ="gray",high="#92d050")+	
			scale_x_discrete(labels=c("Sah-Ara","Hol","Chl-Pat","Neotro","Afr","Indo-Mal","Aus","NZ"))+
			scale_y_discrete(labels=c("Sah-Ara","Hol","Chl-Pat","Neotro","Afr","Indo-Mal","Aus","NZ"))+
			#labs(x='Source',y='Sink',fill="Niche evl rate") +
			geom_text(aes(y=y,label=signif(cor,2),hjust=0.5), size=5,color="black",position = position_dodge(width=0.00),check_overlap = FALSE)+ 
			#geom_text(aes(y=y,label=ifelse(p>0.01,"NS",signif(cor,2)),hjust=0.5), size=2,color="black",position = position_dodge(width=0.00),check_overlap = FALSE)+			
			theme(#panel.grid.major.y=element_blank(),
			panel.grid.minor.y=element_blank(),
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),
			axis.text.y = element_text(angle = 90, hjust = 0.5),    
			axis.line.y=element_line(linetype=1,color='black'),
			axis.line.x=element_line(linetype=1,color='black'),
			axis.ticks = element_line(linetype=2,color='black'),
			panel.grid=element_line(linetype=2,color='grey'),
			panel.background = element_blank(),
			legend.background = element_rect(fill = NA),
			legend.text=element_text(face="bold",size=10),
			legend.title=element_text(face="bold",size=12),
			axis.text=element_text(face="bold",size=11.5),
			axis.title=element_text(face="bold",size=11.5),
			legend.position=c(0.1,0.8)
		  )
		  
#####Fig.4 谱系和过去地质历史的关系#####			
	
##1.计算考虑海洋阻隔的地理距离
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

	#加载raster文件
	load("2_2_rasters_interpolated_land_1d/environmental_variables.RData")
	nmalt <- colnames(environmental_variables[[1]])[-c(1,2)]
	data2load <- paste0("Elevation_",nmalt,"_T.asc")
	nm_res <- which(is.element(as.numeric(colnames(environmental_variables[[1]])[-c(1,2)]),seq(150,0,-1)))
	data2load <- data2load[nm_res]

	#计算地理距离，考虑海洋隔离和海拔高差
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
	
	####glm分析各生物地理区地质历史的影响	
	## hier.part: phlyo beta ~ clim + geo
	library(data.table)
	library(scales)
	library(ggpubr)
	pastdis0=as.data.frame(fread("pastdis.new.csv"))[,-1]
	#匹配地理单元与beta多样性
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
		#画下面的五个地质历史时期的矩形
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
		   #画下面的五个地质历史时期的矩形
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

##code graveyard
#plot isolate
	plot.iso=function(pre,iso.typ){						
		p=pre[["isolation"]]
		plt=ggplot(p)
		# #画下面的五个地质历史时期的矩形
		plt=plt+
		 geom_vline(xintercept=80,col="gray",size=0.5,linetype="longdash",alpha=0.5)+		 
		geom_vline(xintercept=65.5,col="gray",size=0.5,linetype="longdash",alpha=0.5)+		
		 geom_vline(xintercept=55.8,col="gray",size=0.5,linetype="longdash",alpha=0.5)+		  
		  geom_vline(xintercept=33.9,col="gray",size=0.5,linetype="longdash",alpha=0.5)+		
		 geom_vline(xintercept=23.03,col="gray",size=0.5,linetype="longdash",alpha=0.5)+		 	  
		  geom_vline(xintercept=2.58,col="gray",size=0.5,linetype="longdash",alpha=0.5)		 
		if (iso.typ=="geo"){
		plt=plt+		
			geom_line(data=na.omit(p),aes(x=time,y=geo.mean),colour="#00BA38",size=1)+
			geom_ribbon(data = na.omit(p),aes(x = time,ymin=geo.l, ymax=geo.h),fill="#00BA38",alpha=0.5)+
			scale_y_continuous(limits=c(0,0.13))
		}
		if (iso.typ=="clim"){
		plt=plt+		
			geom_line(data=p,aes(x=time,y=clim.mean),colour="#619CFF",size=1)+
			geom_ribbon(data = p,aes(x = time,ymin=clim.l, ymax=clim.h),fill="#619CFF",alpha=0.5)+
			scale_y_continuous(limits=c(0.5,2.5))
		}	
		plt=plt+#facet_wrap( ~ phyregion,ncol = 4,scales='free_x')+
			scale_x_continuous(limits=c(0,80))+			
			theme(panel.grid.major =element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(colour = "black"),
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x  = element_text(size=10),
				axis.text.y  = element_text(size=10))		  
		return(plt)
	}
	
	plotn=lapply(datap,plot.iso,"geo")	
	annotate_figure(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], labels = c("a", "b","c","d", "e","f","g"),label.x=0.1,
			nrow=2,ncol = 4,widths=c(1,1),heights=c(1,1,1,1),common.legend=TRUE),
			bottom = text_grob("Time before present(Ma)", size=14),left = text_grob("Geographic isolation", rot = 90,size=14))
	plotn=lapply(datap,plot.iso,"clim")
	annotate_figure(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], labels = c("a", "b","c","d", "e","f","g"),label.x=0.1,
			nrow=2,ncol = 4,widths=c(1,1),heights=c(1,1,1,1),common.legend=TRUE),
			bottom = text_grob("Time before present(Ma)", size=14),left = text_grob("Historical climatic differences", rot = 90,size=14))
			
##各个大区域的地图
	library(scales)
	show_col(color_bio)
	color_bio=c("#F8766D","#D39200","#93AA00","#00BA38","#00C19F","#619CFF","#DB72FB","#FF61C3")
	phname=as.character(unique(shape@data$phname))
	colbio=data.frame(phname,color_bio)
	regionh=read.csv("node2region.csv")
	region2col=cbind(regionh,colbio[match(regionh$region,colbio$phname),])				
	node=as.character(unique(region2col$node))
	#高级分类单元
	for (i in 1:length(node)){
		plot(shape, col="lightgray", border = "lightgray")
		noderange=region2col[region2col$node==node[i],]
		shp.re1=subset(shape,shape@data$phname%in%noderange$phname)
		for (j in 1:dim(noderange)[1]) {
			shp.re2=subset(shp.re1,shp.re1@data$phname%in%noderange$phname[j])
			plot(shp.re2,col=as.character(noderange$color_bio[j]),border=as.character(noderange$color_bio[j]),add=T)
		}
		savePlot(paste("bioregion_",node[i],sep=""), type=c("jpg"),device=dev.cur(),restoreConsole=TRUE)
	}
	#各个 生物地理区 
	for (i in 1:length(phname)){
	plot(shape, col="lightgray", border = "lightgray")
	noderange=region2col[region2col$phname==phname[i],]
	shp.re3=subset(shape,shape@data$phname%in%noderange$phname)
	colr=as.character(unique(region2col[region2col$phname==phname[i],"color_bio"]))
	plot(shp.re3,col=colr,border=colr,add=T)
	savePlot(paste("bio_tip_",phname[i],sep=""), type=c("jpg"),device=dev.cur(),restoreConsole=TRUE)
	}

	#全部生物地理区
	for (i in 1:length(phname)){
	noderange=region2col[region2col$phname==phname[i],]
	shp.re3=subset(shape,shape@data$phname%in%noderange$phname)
	colr=as.character(unique(region2col[region2col$phname==phname[i],"color_bio"]))
	plot(shp.re3,col=colr,border=colr,add=T)
	}
	savePlot("bio_tip_all", type=c("jpg"),device=dev.cur(),restoreConsole=TRUE)	
	#plot geo isolation and their hp R2
	plot.iso.hp=function(pre,iso.typ){
		if (iso.typ=="geo"){
			p=pre[["geo&r2"]]
		}
		if (iso.typ=="clim"){
			p=pre[["clim&r2"]]
		}
		plt=ggplot(p)+scale_x_continuous(limits=c(0,80))
		# #画下面的五个地质历史时期的矩形
		plt=plt+
		 geom_vline(xintercept=100,col="gray",size=0.5,linetype="longdash",alpha=0.5)+		 
		geom_vline(xintercept=65.5,col="gray",size=0.5,linetype="longdash",alpha=0.5)+		
		 geom_vline(xintercept=55.8,col="gray",size=0.5,linetype="longdash",alpha=0.5)+		  
		  geom_vline(xintercept=33.9,col="gray",size=0.5,linetype="longdash",alpha=0.5)+		
		 geom_vline(xintercept=23.03,col="gray",size=0.5,linetype="longdash",alpha=0.5)+		 	  
		  geom_vline(xintercept=2.58,col="gray",size=0.5,linetype="longdash",alpha=0.5)	
		pdat=na.omit(p)
		if (iso.typ=="geo"){		
		plt=plt+geom_line(data=pdat,aes(x=time,y=iso.mean),colour="#00BA38",size=1,show.legend=TRUE)+
			geom_ribbon(data = pdat,aes(x = time,ymin=iso.l, ymax=iso.h),fill="#00BA38",alpha=0.2)+		
			geom_line(data=pdat,aes(x=time,y=rescale(geo.mean,to=c(min(iso.mean),max(iso.mean)))),colour="darkgray",size=1,alpha=0.8,linetype=2)+
			geom_ribbon(data =pdat,aes(x = time,ymin=rescale(geo.l,to=c(min(iso.mean),max(iso.mean))), 
				ymax=rescale(geo.h,to=c(min(iso.mean),max(iso.mean)))),fill="darkgray",alpha=0.5)+
			scale_y_continuous(limits=c(-0.22,0.85),sec.axis = sec_axis( ~rescale(.,c(min(pdat$iso.mean),max(pdat$iso.mean))),breaks=seq(min(pdat$iso.mean),max(pdat$iso.mean), length.out = 5), labels =seq(round(min(pdat$geo.mean),2),round(max(pdat$geo.mean),2), length.out = 5)))
		}
		if (iso.typ=="clim"){		
		plt=plt+geom_line(data=pdat,aes(x=time,y=iso.mean),colour="#619CFF",size=1,show.legend=TRUE)+
			geom_ribbon(data = pdat,aes(x = time,ymin=iso.l, ymax=iso.h),fill="#619CFF",alpha=0.2)+		
			geom_line(data=p,aes(x=time,y=rescale(clim.mean,to=c(min(pdat$iso.mean),max(pdat$iso.mean)))),colour="darkgray",size=1,alpha=0.8,linetype=2)+
			geom_ribbon(data =p,aes(x = time,ymin=rescale(clim.l,to=c(min(pdat$iso.mean),max(pdat$iso.mean))), 
				ymax=rescale(clim.h,to=c(min(pdat$iso.mean),max(pdat$iso.mean)))),fill="darkgray",alpha=0.5)+
			scale_y_continuous(limits=c(-0.22,0.85),sec.axis = sec_axis( ~rescale(.,c(min(pdat$iso.mean),max(pdat$iso.mean))),breaks=seq(min(pdat$iso.mean),max(pdat$iso.mean), length.out = 5), labels =seq(round(min(pdat$clim.mean),1),round(max(pdat$clim.mean),1), length.out = 5)))				
		}	
		plt=plt+#facet_wrap( ~ phyregion,ncol = 4,scales='free_x')+		 
			theme(panel.grid.major =element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(colour = "black"),
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x  = element_text(size=10),
				axis.text.y  = element_text(size=10))		  
		return(plt)
	}
	plotn=lapply(datap,plot.iso.hp,"clim")
	annotate_figure(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], labels = c("a", "b","c","d", "e","f","g"),label.x=0.1,
			nrow=2,ncol = 4,widths=c(1,1),heights=c(1,1,1,1),common.legend=TRUE),
			bottom = text_grob("Time before present(Ma)", size=14),left = text_grob("Contribution of climatic isolation on floristic division", rot = 90,size=14),right=text_grob("Historical climatic isolation", rot = 90,size=14))
	
	plotn=lapply(datap,plot.iso.hp,"geo")
	annotate_figure(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], labels = c("a", "b","c","d", "e","f","g"),label.x=0.1,
			nrow=2,ncol = 4,widths=c(1,1),heights=c(1,1,1,1),common.legend=TRUE),
			bottom = text_grob("Time before present(Ma)", size=14),left = text_grob("Contribution of geographic isolation on floristic division", rot = 90,size=14),right=text_grob("Geographic isolation", rot = 90,size=14))
	
##NRI
	library(data.table)
	dis <- as.data.frame(fread("Spdat_an_isrm191016.csv",head = T))
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]		
	library(picante)
	tre=read.tree("raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes.tre")#140-150 constraint
	phydata <- match.phylo.comm(tre,spdis)
	
	#install.packages("PhyloMeasures")
	library(PhyloMeasures)	
	regionlist=read.csv("cluster.average.merge.csv")
	spdis.t=phydata$comm
	rownames(spdis.t)=regionlist[match(rownames(spdis.t),regionlist$adcode_ave2),"phname"]
	comm.all=t(sapply(by(spdis.t,rownames(spdis.t),colSums),identity))
	comm.all[comm.all>1]=1
	nri <- -1*mpd.query(phydata$phy, comm.all, TRUE)
	names(nri)=rownames(comm.all)
	# African       Australian Chlie-Patagonian        Holarctic 
    # 2.457840        -6.900246         2.862184         4.790133 
   # Indo-Malesian      Neotropical     Neozeylandic   Saharo-Arabian 
    # -9.697798         1.830508        -4.502146         6.154543
			
	###hier.part: SOS ~ clim + geo
	library(ggpubr)
	library(data.table)
	library(scales)
	pastdis2=as.data.frame(fread("pastdis2.csv"))[,-1]	
	colnames(pastdis2)[-c(1:3)]=paste0("geo",150:0,sep="")	
	past.clim=as.data.frame(fread("past.env.clim.csv"))		
	regionlist=read.csv("cluster.average.merge.csv")
	region=colnames(regionlist)[8:14]	
	f <- function(x, y) paste(x,y,sep="_")
	f2 <- function(x,y)abs(x-y)
	mf=function(x) unlist(strsplit(x,split='_'))	
	sos=readRDS("sos.re150.rds")
	richfin=readRDS("richfin150.rds")
	get_cladeiso=function(regionlist,region.p,richfin.p,sos.p,pastdis2,past.clim,method){
		sos.t=sos.p[-1,as.character(richfin.p$edge2)]
		acd1=regionlist[regionlist[,region.p]==1,"adcode_ave2"]
		acd2=regionlist[regionlist[,region.p]==2,"adcode_ave2"]		
		cladedis=do.call(rbind,lapply(1:dim(sos.t)[2],function(j){
			#for (j in 1:dim(sos.t)[2]){
			clade=richfin.p[richfin.p$edge2==colnames(sos.t)[j],]						
			if (method==1){#属于该clade占主要成分的那一支当前分布的不同生物地理区的地理单元之间距离的mean
				dis=na.omit(sos.t[,j])
				region1.t1=dis[dis>0&names(dis)%in%acd1]
				region1.t2=dis[dis<0&names(dis)%in%acd1]
				if(length(region1.t1)>length(region1.t2)) {region1=region1.t1} else {region1=region1.t2}				
				region2.t1=dis[dis>0&names(dis)%in%acd2]
				region2.t2=dis[dis<0&names(dis)%in%acd2]
				if(length(region2.t1)>length(region2.t2)) {region2=region2.t1} else {region2=region2.t2}			
			}
			if (method==2){#属于该clade当前分布的不同生物地理区的地理单元之间距离的mean
				dis=na.omit(sos.t[,j])
				region1=dis[names(dis)%in%acd1]
				region2=dis[names(dis)%in%acd2]
			}
			if (method==3){#不同生物地理区的地理单元之间距离的mean
				dis=sos.t[,j]
				region1=dis[names(dis)%in%acd1]
				region2=dis[names(dis)%in%acd2]
			}
			adcd <- c(as.vector(outer(names(region1),names(region2), f)),as.vector(outer(names(region2),names(region1),f)))			
			across.geo=pastdis2[pastdis2$id%in%adcd,c("id",paste("geo",round(clade$age),sep=""))]
			colnames(across.geo)[2]="geo"
			a=do.call(rbind,lapply(across.geo$id,mf))
			across.geo2=rbind(across.geo,data.frame(id=paste(a[,2],a[,1],sep="_"),geo=across.geo[,-1]))
			across.clim=past.clim[past.clim$ID%in%adcd,c("ID",paste("clim_",round(clade$age),sep=""))]
			across.com=cbind(across.geo2,across.clim[match(across.geo2$id,across.clim$ID),])
			across=across.com[!is.na(across.com$ID),-3]
			
			# hier.part			
			sos.diff<- c(as.vector(outer(region1,region2, f2)),as.vector(outer(region2,region1,f2)))			
			names(sos.diff)=adcd			
			sos.ana=na.omit(cbind(across,sos.diff[across$id]))
			colnames(sos.ana)[2:4]=c("geo","clim","sos")
			#sos.ana=sos.ana[sos.ana$sos>0,]			
			if(dim(sos.ana)[1]==0) 	{
				isolate=data.frame(clade,geo.med=mean(across[,2]),clim.med=mean(na.omit(across[,3])),
					hp.geo=NA,hp.clim=NA,hp.joint=NA)
			}else{
				
				sos.hp=hier.part::hier.part(sos.ana[,"sos"],sos.ana[,c("geo","clim")],gof = "Rsqu",barplot = FALSE)
				isolate=data.frame(clade,geo.med=mean(across[,2]),clim.med=mean(na.omit(across[,3])),
					hp.geo=sos.hp$IJ["geo","I"],hp.clim=sos.hp$IJ["clim","I"],hp.joint=sos.hp$IJ[1,"J"])
			}		
			#}									
			return(isolate)
		}))		
		return(cladedis)
	}
	
	cladedis=lapply(1:length(region),function(i){
		get_cladeiso(regionlist,region[i],richfin[[i]],sos[[i]],pastdis2,past.clim,method=3)
	})
	
	thes=5
	#plot isolate	
	plotn=lapply(1:length(region),function(i){
		#for(i in 1:length(region)){		
		cladedis2=cbind(age.cat=round(cladedis[[i]]$age/thes)*thes,cladedis[[i]])
		age.cata=sort(unique(cladedis2$age.cat))
		r2=do.call(rbind,lapply(1:length(age.cata),function(k){
			dat=cladedis2[cladedis2$age.cat==age.cata[k],]
			if(dim(dat)[1]>1){
				ci.t=t.test(dat$r2)
				re=data.frame(r2=ci.t$est,r2.l=ci.t$conf[1],r2.h=ci.t$conf[2],r2.se=ci.t$std)
			}else {re=data.frame(r2=mean(dat$r2),r2.l=mean(dat$r2),r2.h=mean(dat$r2),r2.se=0)}
		}))
		rownames(r2)=age.cata
		geo.cld=do.call(rbind,lapply(1:length(age.cata),function(k){
			#for(k in 1:length(age.cata)){
			dat=cladedis2[cladedis2$age.cat==age.cata[k],]
			if(dim(dat)[1]>1&length(unique(dat$geo.med))>1){
				ci.t=t.test(dat$geo.med)
				re=data.frame(geo.cld=ci.t$est,geo.cld.l=ci.t$conf[1],geo.cld.h=ci.t$conf[2],geo.cld.se=ci.t$std)
			}else {re=data.frame(geo.cld=mean(dat$geo.med),geo.cld.l=mean(dat$geo.med),geo.cld.h=mean(dat$geo.med),geo.cld.se=0)}
			#}			
		}))
		rownames(geo.cld)=age.cata
		clim.cld=do.call(rbind,lapply(1:length(age.cata),function(k){
			dat=cladedis2[cladedis2$age.cat==age.cata[k],]
			if(dim(dat)[1]>1&length(unique(dat$clim.med))>1){
				ci.t=t.test(dat$clim.med)
				re=data.frame(clim.cld=ci.t$est,clim.cld.l=ci.t$conf[1],clim.cld.h=ci.t$conf[2],clim.cld.se=ci.t$std)
			}else {re=data.frame(clim.cld=mean(dat$clim.med),clim.cld.l=mean(dat$clim.med),clim.cld.h=mean(dat$clim.med),clim.cld.se=0)}
		}))
		rownames(clim.cld)=age.cata
		p0=cbind(time=age.cata,r2,geo.cld,clim.cld)					
		p=p0[p0$time<=101,]
			
		plt=ggplot(p,aes(x=time,y=r2))+		
			geom_line(data=p,aes(x=time,y=rescale(geo.cld,to=c(min(p$r2),max(p$r2)))),colour="#00BA38",size=1,show.legend=TRUE,alpha=0.8,linetype=2)+
			#geom_errorbar(data = p,aes(x = time,ymin=rescale(geo.cld-geo.cld.se,c(max(p$r2),min(p$r2))), ymax=rescale(geo.cld+geo.cld.se,c(max(p$r2),min(p$r2)))),color="blue",size=0.8,width=1,alpha=0.5)+
			geom_ribbon(data = p,aes(x = time,ymin=rescale(geo.cld-geo.cld.se,to=c(min(p$r2),max(p$r2))), ymax=rescale(geo.cld+geo.cld.se,to=c(min(p$r2),max(p$r2)))),fill="#00BA38",alpha=0.08)+			
			geom_line(data=p,aes(x=time,y=rescale(clim.cld,to=c(min(p$r2),max(p$r2)))),colour="#619CFF",size=1,show.legend=TRUE,alpha=0.8,linetype=2)+
			#geom_errorbar(data = p,aes(x = time,ymin=rescale(clim.cld-clim.cld.se,c(max(p$r2),min(p$r2))), ymax=rescale(clim.cld+clim.cld.se,c(max(p$r2),min(p$r2)))),color="red",size=0.8,width=1,alpha=0.5)+
			geom_ribbon(data = p,aes(x = time,ymin=rescale(clim.cld-clim.cld.se,to=c(min(p$r2),max(p$r2))), ymax=rescale(clim.cld+clim.cld.se,to=c(min(p$r2),max(p$r2)))),fill="#619CFF",alpha=0.08)+
			geom_line(colour="black",show.legend=TRUE,size=1.2)+
			#geom_errorbar(aes(ymin=r2-r2.se, ymax=r2+r2.se),size=1.2,width=1,alpha=0.5)+
			geom_ribbon(data = p,aes(x = time,ymin=r2-r2.se, ymax=r2+r2.se),fill="gray",alpha=0.3)+
			scale_y_continuous(breaks=pretty_breaks(5),limits=c(-1.5,100))+		
			xlim(101,15)			
			#画下面的五个地质历史时期的矩形
		 plt=plt+geom_rect(xmin=-100, xmax=-65.5, ymin=-5, ymax=0,fill= "#B6DCB6")+
		  geom_rect(xmin=-65.5, xmax=-55.8,ymin=-5, ymax=0,fill= "#D2E9E1")+
		  geom_rect(xmin=-55.8, xmax=-33.9, ymin=-5, ymax=0,fill= "#FBEDC9")+
		  geom_rect(xmin=-33.9, xmax=-23.03,ymin=-5, ymax=0,fill="#F8DDA9")+
		  geom_rect(xmin=-23.03, xmax=-15, ymin=-5, ymax=0,fill="#FCB6D0")+		 
		   #画下面的五个地质历史时期的矩形
		  annotate("text", x =65.5+(100-65.5)/2 , y = -1.5 ,label = "K2",size=4)+
		  geom_vline(xintercept=65.5,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =55.8+(65.5-55.8)/2, y =-1.5, label = "E1",size=4) +
		  geom_vline(xintercept=55.8,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =33.9+(55.8-33.9)/2, y =-1.5, label = "E2",size=4) +
		  geom_vline(xintercept=33.9,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =23.03+(33.9-23.03)/2, y = -1.5, label = "E3",size=4) +
		  geom_vline(xintercept=23.03,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =15+(23.03-15)/2, y =-1.5, label = "N1",size=4)
		  
		plt=plt+theme(panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=10),
			axis.text.y  = element_text(size=10))
		#}		
		return(plt)
	})
	annotate_figure(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], labels = c("a", "b","c","d", "e","f","g"),label.x=0.1,
			nrow=2,ncol = 4,widths=c(1,1),heights=c(1,1,1,1),common.legend=TRUE),
			bottom = text_grob("Clade crown age(Ma)", size=14),left = text_grob("Contribution of clades(black)", rot = 90,size=14),
			right = text_grob("Isolation(geo:green/clim:blue)", rot = 90,size=14))
	
	#plot SOS hp r2
	plotn=lapply(1:length(region),function(i){
		#for(i in 1:length(region)){
		cladedis2=cbind(age.cat=round(cladedis[[i]]$age/thes)*thes,cladedis[[i]])
		age.cata=sort(unique(cladedis2$age.cat))
		age.cata=age.cata[age.cata<=100]
		r2=do.call(rbind,lapply(1:length(age.cata),function(k){
			dat=cladedis2[cladedis2$age.cat==age.cata[k],]
			if(dim(dat)[1]>1){
				ci.t=t.test(dat$r2)
				re=data.frame(r2=ci.t$est,r2.l=ci.t$conf[1],r2.h=ci.t$conf[2],r2.se=ci.t$std)
			}else {re=data.frame(r2=mean(dat$r2),r2.l=mean(dat$r2),r2.h=mean(dat$r2),r2.se=0)}
		}))
		rownames(r2)=age.cata		
		hp.joint=do.call(rbind,lapply(1:length(age.cata),function(k){
			dat=cladedis2[cladedis2$age.cat==age.cata[k],]
			if(dim(dat)[1]>1){
				ci.t=t.test(dat$hp.joint)
				re=data.frame(hp=ci.t$est,hp.l=ci.t$conf[1],hp.h=ci.t$conf[2],hp.se=ci.t$std)
			}else {re=data.frame(hp=mean(dat$hp.joint),hp.l=mean(dat$hp.joint),hp.h=mean(dat$hp.joint),hp.se=0)}
		}))
		rownames(hp.joint)=age.cata
		hp.geo=do.call(rbind,lapply(1:length(age.cata),function(k){			
			dat=cladedis2[cladedis2$age.cat==age.cata[k],]
			if(dim(dat)[1]>1&length(unique(dat$hp.geo))>1){
				ci.t=t.test(dat$hp.geo)
				re=data.frame(hp=ci.t$est,hp.l=ci.t$conf[1],hp.h=ci.t$conf[2],hp.se=ci.t$std)
			}else {re=data.frame(hp=mean(dat$hp.geo),hp.l=mean(dat$hp.geo),hp.h=mean(dat$hp.geo),hp.se=0)}						
		}))
		rownames(hp.geo)=age.cata
		hp.clim=do.call(rbind,lapply(1:length(age.cata),function(k){
			dat=cladedis2[cladedis2$age.cat==age.cata[k],]
			if(dim(dat)[1]>1&length(unique(dat$hp.clim))>1){
				ci.t=t.test(dat$hp.clim)
				re=data.frame(hp=ci.t$est,hp.l=ci.t$conf[1],hp.h=ci.t$conf[2],hp.se=ci.t$std)
			}else {re=data.frame(hp=mean(dat$hp.clim),hp.l=mean(dat$hp.clim),hp.h=mean(dat$hp.clim),hp.se=0)}
		}))
		rownames(hp.clim)=age.cata
		p.r=cbind(time=age.cata,r2)
		p.iso=rbind(cbind(time=age.cata,hp.joint,iso.typ="joint"),
		cbind(time=age.cata,hp.geo,iso.typ="geo"),
		cbind(time=age.cata,hp.clim,iso.typ="clim"))		
		p.iso$iso.typ=factor(p.iso$iso.typ,levels=c("joint","geo","clim"),ordered=TRUE)		
		plt=ggplot()+	 
		 geom_line(data=p.r,aes(x=time,y=r2),colour="black",show.legend=TRUE,size=1.2)+
		 geom_ribbon(data =p.r,aes(x = time,ymin=r2-r2.se, ymax=r2+r2.se),fill="gray",alpha=0.3)+
		 scale_y_continuous(breaks=pretty_breaks(5),limits=c(-4,100),sec.axis = sec_axis( ~rescale(.,c(min(p.iso$hp),max(p.iso$hp)))))+
		 xlim(100,15)
			#画下面的五个地质历史时期的矩形
		plt=plt+annotate("rect", xmin=100, xmax=65.5, ymin=-4, ymax=0,fill= "#B6DCB6")+		 
		 annotate("rect", xmin=65.5, xmax=55.8,ymin=-4, ymax=0,fill= "#D2E9E1")+
		  annotate("rect", xmin=55.8, xmax=33.9, ymin=-4, ymax=0,fill= "#FBEDC9")+
		  annotate("rect", xmin=33.9, xmax=23.03,ymin=-4, ymax=0,fill="#F8DDA9")+
		  annotate("rect", xmin=23.03, xmax=15, ymin=-4, ymax=0,fill="#FCB6D0")+		 
		   #画下面的五个地质历史时期的矩形
		  annotate("text", x =65.5+(100-65.5)/2 , y = -2 ,label = "K2",size=3)+
		  geom_vline(xintercept=65.5,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =55.8+(65.5-55.8)/2, y =-2, label = "E1",size=3) +
		  geom_vline(xintercept=55.8,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =33.9+(55.8-33.9)/2, y =-2, label = "E2",size=3) +
		  geom_vline(xintercept=33.9,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =23.03+(33.9-23.03)/2, y = -2, label = "E3",size=3) +
		  geom_vline(xintercept=23.03,col="gray",size=1,linetype="longdash",alpha=0.5)+
		  annotate("text", x =15+(23.03-15)/2, y =-2, label = "N1",size=3)
		plt=plt+geom_line(data=p.iso,aes(x=time,y=rescale(hp,to=c(-4,100)),colour=iso.typ),size=1,show.legend=TRUE,alpha=0.8,linetype=2)+
		 geom_ribbon(data = p.iso,aes(x = time,ymin=rescale(hp-hp.se,to=c(-4,100)), ymax=rescale(hp+hp.se,to=c(-4,100)),fill=iso.typ),alpha=0.1)+
		 scale_color_manual(values = c("#F8766D","#00BA38","#619CFF"))+
		 scale_fill_manual(values = c("#F8766D","#00BA38","#619CFF"))		 
		plt=plt+theme(panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=10),
			axis.text.y  = element_text(size=10))		
		return(plt)
	})
	annotate_figure(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], labels = c("a", "b","c","d", "e","f","g"),label.x=0.1,
			nrow=2,ncol = 4,widths=c(1,1),heights=c(1,1,1,1),common.legend=TRUE),
			bottom = text_grob("Clade crown age(Ma)", size=14),left = text_grob("Contribution of clades(black)", rot = 90,size=14),
			right =text_grob("Isolation HP r2 on clade contribution", rot = 90,size=14))
		
	# ####
	# geo=read.csv("Geo-isrm.csv")
	# regionlist2=cbind(geo[match(regionlist$adcode_ave2,geo$ADCODE99),c("NAME99","continent","Chinese_0")],regionlist)
	# Ori_ind <- list(regionlist2[regionlist2$Chinese_0%in%"印度","adcode_ave2"])
	# names(Ori_ind)="Ori_ind"
	# Ori_may <- list(regionlist2[regionlist2$phname%in%"Oriental"&!regionlist2$Chinese_0%in%"印度","adcode_ave2"])
	# names(Ori_may)="Ori_may"
	# Hol_am <- list(regionlist2[regionlist2$phname%in%"Holarctic"&regionlist2$continent%in%"NAmerica","adcode_ave2"])
	# names(Hol_am)="Hol_am"
	# Hol_euras <- list(regionlist2[regionlist2$phname%in%"Holarctic"&regionlist2$continent%in%c("Europe","Asia"),"adcode_ave2"])
	# names(Hol_euras)="Hol_euras"
	# Neo_cam <- list(regionlist2[regionlist2$phname%in%"Neotropical"&regionlist2$continent%in%"NAmerica","adcode_ave2"])
	# names(Neo_cam)="Neo_cam"
	# Neo_sam <- list(regionlist2[regionlist2$phname%in%"Neotropical"&regionlist2$continent%in%"SAmerica","adcode_ave2"])
	# names(Neo_sam)="Neo_sam"
	# Saharo <- list(regionlist2[regionlist2$phname%in%"Saharo-Arabian","adcode_ave2"])
	# names(Saharo)="Saharo"
	# Aus <- list(regionlist2[regionlist2$phname%in%"Australian","adcode_ave2"])
	# names(Aus)="Aus"
	# Nzland <- list(regionlist2[regionlist2$phname%in%"NewZealand","adcode_ave2"])
	# names(Nzland)="Nzland"
	# Afr <- list(regionlist2[regionlist2$phname%in%"African","adcode_ave2"])
	# names(Afr)="Afr"
	# Pata <- list(regionlist2[regionlist2$phname%in%"Patagonia","adcode_ave2"])
	# names(Pata)="Pata"
	
	# get_acros=function(region1,region2,pastdis2){	
		# f <- function(x, y) paste(x,y,sep="_")
		# geotime=paste0("geo",150:0,sep="")
		# lgmr.t2=vector("list",length(region1)*length(region2));
		# names(lgmr.t2)=paste(rep(names(region1),each=length(region2)),names(region2),sep =".")
		# for(g in 1:length(region1)){
			# for (h in 1:length(region2)){
				# adcd <- unique(c(as.vector(outer(region1[[g]],region2[[h]], f)),as.vector(outer(region2[[h]],region1[[g]],f))))
				# across=pastdis2[pastdis2$id%in%adcd,c("phylo",geotime)]
				# #GLM
				# lgmr.t=do.call(rbind,lapply(1:length(geotime),function(j){
					# past.t=across[across[,geotime[j]]>0,c("phylo",geotime[j])]
					# colnames(past.t)=c("phylo","dis")
					# s=boxplot(past.t$dis)
					# past.t2=subset(past.t,!past.t$dis%in%s$out)					
					# #plot(past.t)
					# #thes=quantile(past.t[,1],0.9)
					# #thes2=quantile(past.t[,1],0.1)
					# #past.t2=past.t[past.t$dis<=thes&past.t$dis>=thes2,]			
					# re.geo=glm(phylo~dis,data=past.t2)
					# slp=as.numeric(coef(re.geo)[2])*sd(past.t2[,"dis"])/sd(past.t2[,"phylo"])
					# temp=data.frame(type=paste(names(region1[g]),names(region2[h]),sep="."),
					# time=as.numeric(substr(geotime[j],4,nchar(geotime[j]))),slp=slp,r2=ss.glm(re.geo)[2],p=ss.glm(re.geo)[3],
					# aic=re.geo$aic,r=cor(past.t2$phylo,past.t2$dis),dis=median(past.t2[,"dis"]))	
					# return(temp)
				# }))
				# lgmr.t2[[h+(g-1)*length(region2)]]=lgmr.t[lgmr.t$p<=0.05,]					
			# }
		# }		
		# return(lgmr.t2)
	# }	
	
	# Hol=c(Hol_am,Hol_euras)	
	# North=c(Hol_am,Hol_euras,Saharo)
	# Neo=c(Neo_cam,Neo_sam)
	# NeoPata=c(Neo_cam,Neo_sam,Pata)
	# AusZland=c(Aus,Nzland)
	# Ori=c(Ori_ind,Ori_may)
	# Paleo=c(Afr,Ori_ind,Ori_may)
	# PaleoAus=c(Afr,Ori_ind,Ori_may,Aus,Nzland)
	# South=c(Afr,Ori_ind,Ori_may,Aus,Nzland,Neo_cam,Neo_sam,Pata)	
	
	# region.dis=list(do.call(rbind,get_acros(North,South,pastdis2)),
		# do.call(rbind,get_acros(Hol,Saharo,pastdis2)),
		# do.call(rbind,get_acros(NeoPata,PaleoAus,pastdis2)),
		# do.call(rbind,get_acros(Neo,Pata,pastdis2)),
		# do.call(rbind,get_acros(Paleo,AusZland,pastdis2)),
		# do.call(rbind,get_acros(Aus,Nzland,pastdis2)),
		# do.call(rbind,get_acros(Afr,Ori,pastdis2))
	# )
	
	# ggplot(region.dis[[1]][region.dis[[1]]$slp>0,],aes(x=time,y=slp,colour=type))+geom_line(show.legend=TRUE)
	
	
	##Fig.4 
	library(ggpubr);library(scales)			
	plotn=vector("list",length(region))		
	for (i in 1:length(region)){				
		p=ggplot()+geom_point(data = lgm.plot, aes(x=time, y=r),size=2,fill="#00BD5F",colour="black",shape=21,show.legend=FALSE)+	scale_size_area(max_size = 3.5)			
		# d=glm(slp~time,data = lgm.plot)	
		# if(ss.glm(d)[3]<0.05) p=p+geom_smooth(aes(x=time,y=slp),data=lgm.plot,method = "glm",colour="red",fill="#00BD5F",size=1)
		# if(ss.glm(d)[3]<0.1&ss.glm(d)[3]>0.05) p=p+geom_smooth(aes(x=time,y=slp),data=lgm.plot,method = "glm",colour="red",,fill="#00BD5F",linetype="dashed")#marginal sig		
		p=p+theme(panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=10),
			axis.text.y  = element_text(size=10))+			
		xlim(150,0)
		#p=p+geom_smooth(method = "loess",colour="red",fill="lightgray",size=1)			
		plotn[[i]]=p
		
	}
	annotate_figure(ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], labels = c("a", "b","c","d", "e","f","g"),label.x=0.1,
			nrow=2,ncol = 4,widths=c(1,1),heights=c(1,1,1,1),common.legend=TRUE),
			bottom = text_grob("Time before present (Ma)", size=14),left = text_grob("Contribution of plate tectonics", rot = 90,size=14),
			right = text_grob("Climate/Geographic isolation", rot = 90,size=14))
	
	##Fig.4-extend plaeo maps
	library(raster)
	resList <- readRDS("resList.rds")	
	nmalt=seq(150,0,by=-1)	
	cordacd <- read.csv("POINT2.csv")
	phname=read.csv("cluster.average.merge.csv")	
	cordacd=subset(cordacd,cordacd$ADCODE99>0)
	rownames(cordacd)=paste(cordacd[,1],cordacd[,2],sep="_")	
	sublist=resList[match(rownames(cordacd),names(resList))]
	cordacd2=cbind(cordacd,phname=phname[match(cordacd$ADCODE99,phname$adcode_ave2),"phname"])
	library(parallel)
	no_cores <- detectCores() - 1
	mycl <- makePSOCKcluster(no_cores);
	data.test=parLapply(cl=mycl, X=nmalt,function(i,sublist,cordacd2,phname,nmalt){
		require(raster)
		pC.t=c()
		for (j in 1:length(sublist)){
			temp=sublist[[j]][202-i,c(3,2)]#i+51=150
			pC.t=rbind(pC.t,temp)
		}
		rownames(pC.t)=names(sublist)				
		aa=floor(as.matrix(pC.t)/0.5)*0.5
		aa[which(floor(aa[,1])-aa[,1]==0),1]=aa[which(floor(aa[,1])-aa[,1]==0),1]+0.5
		aa[which(floor(aa[,2])-aa[,2]==0),2]=aa[which(floor(aa[,2])-aa[,2]==0),2]+0.5
		pC.t2=cbind(adcode=cordacd2$ADCODE99,phname=cordacd2$phname,aa[match(rownames(cordacd2),rownames(aa)),])	
		pC.t3=data.frame(ID=paste(pC.t2[,"lon"],pC.t2[,"lat"],sep="_"),pC.t2)	
		
		data.test <- raster(paste0("2_2_rasters_interpolated_land_1d/temperature/global/rasters/",i,".00_temperature.asc",sep=""))
		coord=coordinates(data.test)
		rownames(coord)=paste(coord[,1],coord[,2],sep="_")
		pC.t4=cbind(coord,pC.t3[match(rownames(coord),pC.t3$ID),])
		data.test[]=pC.t4$phname
		return(data.test)
	},sublist,cordacd2,phname,nmalt)	
	names(data.test)=nmalt ##每一个列表为所有ADCODE与古经纬度的对应关系，共150个Ma的列表
	saveRDS(data.test,"paleomaps/paleomaps.rds")
	stopCluster(mycl)
	
	nmalt=seq(150,0,by=-1)	
	data.test=readRDS("paleomaps/paleomaps.rds")
	for (i in 1:length(nmalt)){
		jpeg(file=paste("paleomaps/foo",nmalt[i],".jpg",sep=""))		
		arg <- list(at=c(1:8), labels=c("African","Australian","Holarctic","Neotropic","New Zealand","Oriental","Patagonia","Saharo-Arabian"))
		raster::plot(data.test[[as.character(nmalt[i])]],col=c("#00BFC4","#C77CFF","#CD9600","#F8766D","#FF61CC","#00BE67","#00A9FF","#7CAE00"),main=paste(nmalt[i],"Mya"),
			axis.args=arg,xlab="Lat",ylab="Lon")
		dev.off()
	}	
	
	library(ggpubr)	
	lgmr=read.csv("lgmr10.csv")[,-1]	
	region=unique(as.character(lgmr$region))[-8]
	lgmr2=lgmr[lgmr$region!="all",]
	lgmr2$region=factor(lgmr2$region,levels=c( "World","Laurasian","Gondwanan","Neotropic","Indo.Ocean","Antarctic","Palaeotropic"))
	for (i in 1:length(nmalt)){
	logo <-  image_read(paste("paleomaps/foo",nmalt[i],".jpg",sep=""))
	p=ggplot(lgmr2, aes(x=time, y=geo.r2))+geom_point(size=2,fill="#00BD5F",colour="black",shape=21,show.legend=FALSE)+scale_size_area(max_size = 3.5)+
	facet_wrap( ~ region,ncol = 3,scales='fixed')+
	theme(strip.text= element_blank(),
		panel.grid.major =element_blank(), 
			panel.grid.minor = element_blank(),
			#panel.background = element_blank(),
			axis.line = element_line(colour = "black"),
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=10),
			axis.text.y  = element_text(size=10))+			
		xlim(150,0)+ylim(0,28)+ geom_vline(aes(xintercept=nmalt[i]), colour="#990000", linetype="dashed")
	p2=ggplot()+annotation_raster(logo,-Inf, Inf, -Inf, Inf)
	p3=ggarrange(p,p2,nrow=1,ncol = 2,widths=1,heights=c(1,1))
	ggsave(paste("paleomaps/eft",nmalt[i],".jpg",sep=""),p3)
	}
	##make a gif, slow, alternatively done in PS. ref https://blog.hubspot.com/marketing/add-text-animated-gif-photoshop-tutorial; https://blog.hubspot.com/marketing/how-to-create-animated-gif-quick-tip-ht
	#install.packages("magick")
	nmalt=seq(150,0,by=-1)	
	library(magick)	
	
	eft=vector("list",length(nmalt));names(eft)=nmalt
	for (i in 1:length(nmalt)){
	eft[[i]]=image_read(paste("paleomaps/eft",nmalt[i],".jpg",sep=""))
	}
	saveRDS(eft,"paleomaps/eft.rds")
	
	animation <- image_animate(image_join(eft))	
	image_write(animation, "paleomaps/paleomaps.gif")
	
#####Loic 交流， 未放到文中的分析########			
##2.loic的环境数据，1°网格的气候数据
	library(raster)
	resList <- readRDS("resList.rds")
	nmalt=seq(150,0,by=-1)
	data2load <- paste0(nmalt,".00_temperature.asc",sep="")	
	data.test <- raster(paste0("2_2_rasters_interpolated_land_1d/temperature/global/rasters/",data2load[1]))
	coord=coordinates(data.test)
	rownames(coord)=paste(coord[,1],coord[,2],sep="_")	
	
	Alt.raster <-do.call(cbind,lapply(1:length(data2load),function(i){
	  data <- raster(paste0("2_2_rasters_interpolated_land_1d/temperature/global/rasters/",data2load[i]))	  
	  data[]
	}))
	Alt.raster <- cbind(coord,Alt.raster)
	colnames(Alt.raster)[-c(1,2)] <- nmalt

	data2load.p <- paste0(nmalt,".00_precipitation.asc",sep="")	
	Alt.raster.p <-do.call(cbind,lapply(1:length(data2load.p),function(i){
	  data.p <- raster(paste0("2_2_rasters_interpolated_land_1d/precipitation/global/rasters/",data2load.p[i]))	 
	  data.p[]
	}))
	Alt.raster.p <- cbind(coord,Alt.raster.p)
	colnames(Alt.raster.p)[-c(1,2)] <- nmalt
	
	cordacd <- read.csv("POINT2.csv")
	cordacd=subset(cordacd,cordacd$ADCODE99>0)
	rownames(cordacd)=paste(cordacd[,1],cordacd[,2],sep="_")		
	sublist=resList[match(rownames(cordacd),names(resList))]	
	pC=vector("list",151);names(pC)=seq(150,0,-1)##每一个列表为所有ADCODE与古经纬度的对应关系，共150个Ma的列表
	for (i in 1:151){
		pC.t=c()
		for (j in 1:length(sublist)){
			temp=sublist[[j]][i+51,c(3,2)]
			pC.t=rbind(pC.t,temp)
		}
		rownames(pC.t)=names(sublist)				
		aa=floor(as.matrix(pC.t)/0.5)*0.5
		aa[which(floor(aa[,1])-aa[,1]==0),1]=aa[which(floor(aa[,1])-aa[,1]==0),1]+0.5
		aa[which(floor(aa[,2])-aa[,2]==0),2]=aa[which(floor(aa[,2])-aa[,2]==0),2]+0.5
		pC.t2=cbind(adcode=cordacd$ADCODE99,aa[match(rownames(cordacd),rownames(aa)),])	
		pC[[i]]=data.frame(ID=paste(pC.t2[,"lon"],pC.t2[,"lat"],sep="_"),pC.t2)
	}
	
	past.temp=data.frame(adcode=unique(pC[[1]]$adcode))
	past.prep=data.frame(adcode=unique(pC[[1]]$adcode))
	for (i in 1:151){
		tmp.temp=cbind(pC[[i]],Alt.raster[match(pC[[i]]$ID,rownames(Alt.raster)),i+2])
		tmp0=na.omit(tmp.temp)
		tmp1=tapply(X=tmp0[,5],INDEX=tmp0$adcode,FUN=mean)
		past.temp=cbind(past.temp,tmp1[match(past.temp$adcode,names(tmp1))])
		
		tmp.prep=cbind(pC[[i]],Alt.raster.p[match(pC[[i]]$ID,rownames(Alt.raster.p)),i+2])
		tmp0=na.omit(tmp.prep)
		tmp1=tapply(X=tmp0[,5],INDEX=tmp0$adcode,FUN=mean)
		past.prep=cbind(past.prep,tmp1[match(past.prep$adcode,names(tmp1))])
	}
	colnames(past.temp)=c("ADCODE99",paste("temp",seq(150,0,-1),sep="_"))
	colnames(past.prep)=c("ADCODE99",paste("prep",seq(150,0,-1),sep="_"))
	
	clim=read.csv("richclim191024.csv")
	phname=read.csv("cluster.average.merge.csv")	
	past.clim=cbind(clim[match(past.temp$ADCODE99,clim$adcode),c("adcode","Lon","Lat")],phname[match(past.temp$ADCODE99,phname$adcode_ave2),-c(1:2)],
				past.temp[,c(2:length(past.temp))],past.prep[match(past.temp$ADCODE99,past.prep$ADCODE99),c(2:length(past.prep))])						
	write.csv(past.clim,"past.clim.csv")
	
	past.clim=read.csv("past.clim.csv")
	past.clim=past.clim[,2:length(past.clim)]
	
	#温度标准化，降水*365后标准化
	past.clim2=past.clim[,1:15]
	for (i in 16:dim(past.clim)[2]){
	if (i<167) temp=scale(past.clim[,i]) else temp=scale(past.clim[,i]*365)
	past.clim2=cbind(past.clim2,temp)
	}
	colnames(past.clim2)=colnames(past.clim)
	write.csv(past.clim2,"past.clim2.csv")
	
##3.整合过去气候和地质历史的数据以及beta多样性
	beta.sim.phlo=read.csv("phylobetasimfin_an150.csv")
	dis.sim=read.csv("Spdat_an_isrm191016.csv")	
	phylo=betamatrix2(beta.sim.phlo,dis.sim)
	
	physim=as.data.frame(phylo[match(past.clim2$adcode,colnames(phylo)),match(past.clim2$adcode,rownames(phylo))])
	Dissim.phylo=cbind(adcode=past.clim2$adcode,physim)
	exFormat.all=formatsitepair(Dissim.phylo,3, XColumn="Lon", YColumn="Lat", predData=data.frame(past.clim2,Adcode99=past.clim2$adcode),siteColumn="adcode")
	
	beta_pastenv <-function(exFormat){		
		dis.clim=exFormat[,c("s1.Adcode99","s2.Adcode99")]
		
		dis.temp=c()
		dis.prep=c()		
		for (i in 150:0){
			#温度和降水的欧式距离
			dis.clim.t=sqrt((exFormat[,paste("s1.temp",i,sep="_")]-exFormat[,paste("s2.temp",i,sep="_")])^2+
					(exFormat[,paste("s1.prep",i,sep="_")]-exFormat[,paste("s2.prep",i,sep="_")])^2)
			dis.clim=cbind(dis.clim,dis.clim.t)
			#温度的差异	=|a-b|
			dis.temp.t=abs(exFormat[,paste("s1.temp",i,sep="_")]-exFormat[,paste("s2.temp",i,sep="_")])#/
					#max(abs(exFormat[,paste("s1.temp",i,sep="_")]),abs(exFormat[,paste("s2.temp",i,sep="_")]))
			dis.temp=cbind(dis.temp,dis.temp.t)
			#降水的差异=|a-b|
			dis.prep.t=abs(exFormat[,paste("s1.prep",i,sep="_")]-exFormat[,paste("s2.prep",i,sep="_")])
			dis.prep=cbind(dis.prep,dis.prep.t)			
			}
		colnames(dis.clim)=	c("s1.Adcode99","s2.Adcode99",paste("clim",150:0,sep="_"))	
		colnames(dis.temp)=	c(paste("temp",150:0,sep="_"))
		colnames(dis.prep)=	c(paste("prep",150:0,sep="_"))			
		r.tmp=cbind(dis.clim,dis.temp,dis.prep)		
		return(r.tmp)
		}

	past.climdis=beta_pastenv(exFormat.all)	
	
	exFormat.phy=cbind(ID=paste(past.climdis$s1.Adcode99,past.climdis$s2.Adcode99,sep="_"),past.climdis[,-c(1:2)])
	
	pastdis=read.csv("pastdis2.csv")
	id.t=as.character(pastdis$id)
	id.t2=strsplit(id.t,split="_")
	id=do.call(rbind,lapply(1:length(id.t2),function(i){
	  paste(id.t2[[i]][2],id.t2[[i]][1],sep="_")	  
	}))
	pastdis2=rbind(pastdis[,-1],data.frame(id=as.character(id),pastdis[,-c(1:2)]))
	
	past.env=cbind(exFormat.phy,pastdis2[match(exFormat.phy$ID,pastdis2$id),-1])#地质历史和过去气候与物种beta多样性的对应关系
	write.csv(past.env,"past.env.all.csv")##过去气候和地质历史的数据以及beta多样性, past.env.clim.csv是只保留clim欧式距离的表
	
##4.偏回归分析过去气候和地质历史对物种和谱系beta的相对贡献
	past.env=read.csv("past.env.all.csv")	
	hpr=c();hpr.sp=c()
	for (i in 150:0){
		re=hier.part(past.env[,"phylo"],past.env[,c(paste("X",i,sep=""),paste("clim",i,sep="_"))],gof = "Rsqu",barplot = FALSE)
		temp=c(time=i,geo=re$IJ$I[1],clim=re$IJ$I[2],joint=re$IJ$J[1])
		hpr=rbind(hpr,temp)
		
		re.sp=hier.part(past.env[,"sp"],past.env[,c(paste("X",i,sep=""),paste("clim",i,sep="_"))], gof = "Rsqu",barplot = FALSE)
		temp.sp=c(time=i,geo=re.sp$IJ$I[1],clim=re.sp$IJ$I[2],joint=re.sp$IJ$J[1])
		hpr.sp=rbind(hpr.sp,temp.sp)		
	}	
	
	#作图	
	r.all=rbind(data.frame(betatype="phylo",hpr),data.frame(betatype="sp",hpr.sp))
	r.plot=c()
	for (i in 3:length(colnames(r.all))){
		r.tmp=data.frame(rsqtype=colnames(r.all)[i],r.all[,c(1:2,i)])
		colnames(r.tmp)[4]="rsq"
		r.plot=rbind(r.plot,r.tmp)
	}
	write.csv(r.plot,"r.plot.phylosp.csv")
	r.plot=read.csv("r.plot.phylosp.csv")
	p=ggplot(r.plot, aes(x=time, y=rsq, colour=rsqtype,group=rsqtype)) + geom_line(size=1) + facet_grid(betatype ~ .)+scale_x_reverse()
		
	ggsave("plot.phylosp.pdf",plot=p)
	
##5.生物地理区内和生物地理区间beta与过去气候的关系
	regionlist=read.csv("past.clim2.csv")[,-1]	
	region=as.character(unique(regionlist$phname))		
	
	exFormat.within=c()
	for (i in 1:length(region)){		
		envir=regionlist						
		envt=subset(envir,envir[,"phname"]==region[i])		
		physim=as.data.frame(phylo[match(envt$adcode,colnames(phylo)),match(envt$adcode,rownames(phylo))])
		Dissim=cbind(adcode=envt$adcode,physim)
		exFormat_pre <- formatsitepair(Dissim,3, XColumn="Lon", YColumn="Lat", predData=data.frame(envt,Adcode99=envt$adcode),siteColumn="adcode")						
		exFormat.within=rbind(exFormat.within,exFormat_pre)	
	}
	
	all.temp=cbind(ID=paste(exFormat.all[,"s1.Adcode99"],exFormat.all[,"s2.Adcode99"],sep="_"),exFormat.all)
	within.temp=cbind(ID=paste(exFormat.within[,"s1.Adcode99"],exFormat.within[,"s2.Adcode99"],sep="_"),exFormat.within)
	exFormat.across=all.temp[-match(within.temp$ID,all.temp$ID),-1]
	
	past.climdis.within=beta_pastenv(exFormat.within)	
	past.climdis.across=beta_pastenv(exFormat.across)
	
	exFormat.phy.within=cbind(ID=paste(past.climdis.within$s1.Adcode99,past.climdis.within$s2.Adcode99,sep="_"),past.climdis.within[,-c(1:2)])		 
	past.env.within=cbind(exFormat.phy.within,pastdis2[match(exFormat.phy.within$ID,pastdis2$id),-1])#地质历史和过去气候与物种beta多样性的对应关系
	
	exFormat.phy.across=cbind(ID=paste(past.climdis.across$s1.Adcode99,past.climdis.across$s2.Adcode99,sep="_"),past.climdis.across[,-c(1:2)])		 
	past.env.across=cbind(exFormat.phy.across,pastdis2[match(exFormat.phy.across$ID,pastdis2$id),-1])#地质历史和过去气候与物种beta多样性的对应关系
	write.csv(past.env.within,"past.env.within.csv")
	write.csv(past.env.across,"past.env.across.csv")
	
	library(data.table)
	past.env.within=as.data.frame(fread("past.env.within.csv"))[,-1]	
	past.env.across=as.data.frame(fread("past.env.across.csv"))[,-1]
	hpr.within=c();hpr.across=c()
	for (i in 150:0){
		re=hier.part(past.env.within[,"phylo"],past.env.within[,c(paste("X",i,sep=""),paste("clim",i,sep="_"))], gof = "Rsqu",barplot = FALSE)
		temp.within=c(time=i,geo=re$IJ$I[1],clim=re$IJ$I[2],joint=re$IJ$J[1])
		hpr.within=rbind(hpr.within,temp.within)		
		
		re.sp=hier.part(past.env.across[,"phylo"],past.env.across[,c(paste("X",i,sep=""),paste("clim",i,sep="_"))], gof = "Rsqu",barplot = FALSE)
		temp.across=c(time=i,geo=re.sp$IJ$I[1],clim=re.sp$IJ$I[2],joint=re.sp$IJ$J[1])
		hpr.across=rbind(hpr.across,temp.across)		
	}	
	
	#作图	
	r.all=rbind(data.frame(betatype="within bio-regions",hpr.within),data.frame(betatype="acorss bio-regions",hpr.across))
	r.plot=c()
	for (i in 3:length(colnames(r.all))){
		r.tmp=data.frame(rsqtype=colnames(r.all)[i],r.all[,c(1:2,i)])
		colnames(r.tmp)[4]="rsq"
		r.plot=rbind(r.plot,r.tmp)
	}
	write.csv(r.plot,"r.plot.withacross.csv")
	r.plot=read.csv("r.plot.withacross.csv")
	p=ggplot(r.plot, aes(x=time, y=rsq, colour=rsqtype,group=rsqtype)) + geom_line(size=1) + facet_grid(betatype ~ .)+scale_x_reverse()
	ggsave("plot.withacross.pdf",plot=p)
	
	
##6.偏回归分析过去气候和地质历史对各生物地理区边界的相对贡献	
	regionlist=read.csv("past.clim2.csv")[,-1]		
	node=colnames(regionlist)[9:15]		
	past.env=vector("list",length(node));names(past.env)=paste("be_",node,sep="")	
	for (i in 1:length(node)){	
		envir=regionlist				
		envt=subset(envir,envir[,node[i]]>0)				
		physim=as.data.frame(phylo[match(envt$adcode,colnames(phylo)),match(envt$adcode,rownames(phylo))])
		Dissim=cbind(adcode=envt$adcode,physim)
		exFormat_pre <- formatsitepair(Dissim,3, XColumn="Lon", YColumn="Lat", predData=data.frame(envt,Adcode99=envt$adcode),siteColumn="adcode")			
		
		cole=6+(length(colnames(exFormat_pre))-6)/2	
		cole2=length(colnames(exFormat_pre))
		exFormat=subset(exFormat_pre,
				  (exFormat_pre[,paste("s1.",node[i],sep="")]==1 & exFormat_pre[,paste("s2.",node[i],sep="")]==2)|
					(exFormat_pre[,paste("s1.",node[i],sep="")]==2 & exFormat_pre[,paste("s2.",node[i],sep="")]==1))			
		past.climdis=beta_pastenv(exFormat)	
		exFormat.phy=cbind(ID=paste(past.climdis$s1.Adcode99,past.climdis$s2.Adcode99,sep="_"),past.climdis[,-c(1:2)])		 
		past.env[[i]]=cbind(exFormat.phy,pastdis2[match(exFormat.phy$ID,pastdis2$id),-1])#地质历史和过去气候与物种beta多样性的对应关系	
	}
	saveRDS(past.env,"past.env.RDS")	
	
	
	#hp
	hpr=c()
	for (j in 1:length(node)){
		hpr.t=c()
		for (i in 150:0){
			if (dim(na.omit(past.env[[j]][,c(paste("X",i,sep=""),paste("clim",i,sep="_"))]))[1]==0) temp=c(time=i,geo=0,clim=0,joint=0)
			else {re=hier.part(past.env[[j]][,"phylo"],past.env[[j]][,c(paste("X",i,sep=""),paste("clim",i,sep="_"))], gof = "Rsqu",barplot = FALSE)
			temp=c(time=i,geo=re$IJ$I[1],clim=re$IJ$I[2],joint=re$IJ$J[1])}
			hpr.t=rbind(hpr.t,temp)			
		}
		r.plot=c()
		for (k in 2:4){
		r.tmp=data.frame(rsqtype=colnames(hpr.t)[k],hpr.t[,c(1,k)])
		colnames(r.tmp)[3]="rsq"
		r.plot=rbind(r.plot,r.tmp)		
		}
		hpr=rbind(hpr,cbind(region=node[j],r.plot))
	}
	
	write.csv(hpr,"hpr.csv")
	
	#作图		
	p=ggplot(hpr, aes(x=time, y=rsq, colour=rsqtype,group=rsqtype)) + geom_line(size=1) + facet_wrap(~region)+scale_x_reverse()
	ggsave("plot.bioregion.pdf",plot=p)
	#作图2
	plotn=vector("list",length(node))	
	for (i in 1:length(node)){
		hp.plot=subset(hpr,hpr$region==node[i])		
		p=ggplot(hp.plot, aes(x=time, y=rsq, colour=rsqtype,group=rsqtype)) + geom_line(size=1) +scale_x_reverse()+
			   theme_bw() + 			
			theme(axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=8),
			axis.text.y  = element_text(size=8))+
			theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
			  panel.background = element_blank(),axis.line = element_line(colour = "black"))

		if (i==1){
		p=p+theme(legend.key.width=unit(0.5,'cm'),
			legend.justification = c(0, 1), # pick the upper left corner of the legend box and
			legend.position = c(0, 1), # adjust the position of the corner as relative to axis
			legend.background = element_rect(fill = NA), # transparent legend background
			legend.box = "horizontal", # horizontal arrangement of multiple legends
			legend.spacing.x = unit(0.5, units = "cm"))# horizontal spacing between legends 	
			
		}
		if (i>1) p=p+guides(fill=FALSE)
		plotn[[i]]=p
	}	
	
	figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]], 
			nrow=3,ncol = 3,widths=c(2,2,2),heights=c(2,2,2),common.legend=TRUE,legend = "top",label.x=0.5,hjust=0,align="hv")
	ggsave("plot.bioregion.pdf",plot=figure)	

##8.SEM
	#每个网格的净分化速率	
	dis=read.csv("Spdat_an_isrm191016.csv")
	rownames(dis)=dis[,1]
	dis=dis[,-1]
	region=rep(rownames(dis),dim(dis)[2])
	sp=c()
	for (i in 1:dim(dis)[2]){
	sp=c(sp,rep(colnames(dis)[i],dim(dis)[1]))
	}
	spdis=data.frame(dis=region,sp=sp,pre=as.numeric(as.matrix(dis)))
	spdis=subset(spdis[,1:2],spdis$pre>0)
	write.csv(spdis,"Spdat_an_isrm191016.2cols.csv")
	
	netdiv=read.csv("netdiv.csv")
	spdis=read.csv("Spdat_an_isrm191016.2cols.csv")
	spnet=na.omit(cbind(spdis[,-1],netdiv[match(spdis$sp,netdiv$genus),]))
	geonet=tapply(spnet$Av2_a1_sp_netdiv,spnet$dis,mean)	
	geoext=tapply(spnet$Av2_a1_sp_ext,spnet$dis,mean)
	geoage=tapply(spnet$Av2_a1_sp_age,spnet$dis,mean)
	
	#过去气候的标准差
	past.clim=read.csv("past.clim2.csv")
	temp=apply(past.clim[,17:167],1,function(i) sd(na.omit(i)))
	prep=apply(past.clim[,168:318],1,function(i) sd(na.omit(i)))
	
	t.anc.mean=apply(past.clim[,17:101],1,function(i) mean(na.omit(i)))#66Ma mass extinction
	t.anc.sd=apply(past.clim[,17:101],1,function(i) sd(na.omit(i)))
	t.pre.mean=apply(past.clim[,102:167],1,function(i) mean(na.omit(i)))#66Ma mass extinction
	t.pre.sd=apply(past.clim[,102:167],1,function(i) sd(na.omit(i)))
	
	p.anc.mean=apply(past.clim[,168:252],1,function(i) mean(na.omit(i)))#66Ma mass extinction
	p.anc.sd=apply(past.clim[,168:252],1,function(i) sd(na.omit(i)))
	p.pre.mean=apply(past.clim[,253:318],1,function(i) mean(na.omit(i)))#66Ma mass extinction
	p.pre.sd=apply(past.clim[,253:318],1,function(i) sd(na.omit(i)))	
	
	clim=cbind(past.clim[,2:16],temp,prep,t.anc.mean,t.anc.sd,t.pre.mean,t.pre.sd,p.anc.mean,p.anc.sd,p.pre.mean,p.pre.sd)	
	climnet=cbind(clim,netdiv=geonet[match(clim$adcode,names(geonet))],ext=geoext[match(clim$adcode,names(geoext))],age=geoage[match(clim$adcode,names(geoage))])

	beta.sim.phlo=read.csv("phylobetasimfin_an150.csv")
	dis.sim=read.csv("Spdat_an_isrm191016.csv")	
	phylo=betamatrix2(beta.sim.phlo,dis.sim)
		
	physim=as.data.frame(phylo[match(climnet$adcode,colnames(phylo)),match(climnet$adcode,rownames(phylo))])
	Dissim.phylo=cbind(adcode=climnet$adcode,physim)
	exFormat.all=formatsitepair(Dissim.phylo,3, XColumn="Lon", YColumn="Lat", predData=data.frame(climnet,Adcode99=climnet$adcode),siteColumn="adcode")	
	pastdis.all=read.csv("pastdis2.csv")
	
	get.pastenv<-function(exFormat,pastdis){
		dis.clim=sqrt((exFormat[,"s1.temp"]-exFormat[,"s2.temp"])^2+(exFormat[,"s1.prep"]-exFormat[,"s2.prep"])^2)
		# dis.temp=abs(exFormat[,"s1.temp"]-exFormat[,"s2.temp"])
		# dis.prep=abs(exFormat[,"s1.prep"]-exFormat[,"s2.prep"])
		
		clim.anc.m=sqrt((exFormat[,"s1.t.anc.mean"]-exFormat[,"s2.t.anc.mean"])^2+(exFormat[,"s1.p.anc.mean"]-exFormat[,"s2.p.anc.mean"])^2)
		clim.anc.s=sqrt((exFormat[,"s1.t.anc.sd"]-exFormat[,"s2.t.anc.sd"])^2+(exFormat[,"s1.p.anc.sd"]-exFormat[,"s2.p.anc.sd"])^2)
		clim.pre.m=sqrt((exFormat[,"s1.t.pre.mean"]-exFormat[,"s2.t.pre.mean"])^2+(exFormat[,"s1.p.pre.mean"]-exFormat[,"s2.p.pre.mean"])^2)
		clim.pre.s=sqrt((exFormat[,"s1.t.pre.sd"]-exFormat[,"s2.t.pre.sd"])^2+(exFormat[,"s1.p.pre.sd"]-exFormat[,"s2.p.pre.sd"])^2)
		
		dis.netdiv=abs(exFormat[,"s1.netdiv"]-exFormat[,"s2.netdiv"])
		dis.ext=abs(exFormat[,"s1.ext"]-exFormat[,"s2.ext"])
		dis.age=abs(exFormat[,"s1.age"]-exFormat[,"s2.age"])
		
		past.climdis=cbind(exFormat[,c("s1.Adcode99","s2.Adcode99")],dis.clim,clim.anc.m,clim.anc.s,clim.pre.m,clim.pre.s,dis.netdiv,dis.ext,dis.age)		
		exFormat.phy=cbind(ID=paste(past.climdis$s1.Adcode99,past.climdis$s2.Adcode99,sep="_"),past.climdis[,-c(1:2)])
		
		geo=apply(pastdis[,-c(1:4)],1,sd)
		g.anc.m=apply(pastdis[,c(5:89)],1,mean)
		g.anc.s=apply(pastdis[,c(5:89)],1,sd)
		g.pre.m=apply(pastdis[,c(90:155)],1,mean)
		g.pre.s=apply(pastdis[,c(90:155)],1,sd)
		
		pastdis=cbind(pastdis[,2:4],dis.geo=geo,g.anc.m=g.anc.m,g.anc.s=g.anc.s,g.pre.m=g.pre.m,g.pre.s=g.pre.s)
		id.t=as.character(pastdis$id)
		id.t2=strsplit(id.t,split="_")
		id=do.call(rbind,lapply(1:length(id.t2),function(i){
			paste(id.t2[[i]][2],id.t2[[i]][1],sep="_")	  
		}))
		pastdis2=rbind(pastdis,data.frame(id=as.character(id),pastdis[,-1]))
	
		#sem分析用的表格
		past.env=cbind(exFormat.phy,pastdis2[match(exFormat.phy$ID,pastdis2$id),-1])#地质历史和过去气候与物种beta多样性的对应关系
		return(past.env)
	}	
	
	#sem分析用的表格
	past.env=get.pastenv(exFormat.all,pastdis.all)#地质历史和过去气候与物种beta多样性的对应关系
	
	#write.csv(past.env,"sem.csv")
	
	#sem分析
	sem.data=past.env   
	#install.packages("lavaan", dependencies = TRUE)  # 安装lavaan包
	#install.packages("semPlot")
	library(lavaan)
	library(semPlot)
	
	##可视化中介效应 https://zhuanlan.zhihu.com/p/53206137,https://rpubs.com/cardiomoon/239332
	# install.packages("devtools")
	# devtools::install_github("guhjy/semMediation")	###这个包已经没了
	myModel<- " phylo~b1*dis.netdiv+b2*dis.ext+c1*dis.clim+c2*dis.geo
	dis.netdiv~a1*dis.clim+a2*dis.geo
	dis.ext~a3*dis.clim+a4*dis.geo+d1*dis.netdiv
	ind1:=a1*b1
	ind2:=a2*b1
	ind3:=a3*b2
	ind4:=a4*b2
	secondInd1:=a1*d1*b2+a2*d1*b2
	total1:=c1+c2+a1*b1+a2*b1+a3*b2+a4*b2+a1*d1*b2+a2*d1*b2	
		"		
	fit <- sem(myModel, data = scale(na.omit(sem.data[,-1])))
	summary (fit, rsquare=TRUE, standardized=TRUE, fit.measures=TRUE)
	
	myModel.sp<- " sp~b1*dis.netdiv+b2*dis.ext+c1*dis.clim+c2*dis.geo
	dis.netdiv~a1*dis.clim+a2*dis.geo
	dis.ext~a3*dis.clim+a4*dis.geo+d1*dis.netdiv
	ind1:=a1*b1
	ind2:=a2*b1
	ind3:=a3*b2
	ind4:=a4*b2
	secondInd1:=a1*d1*b2+a2*d1*b2
	total1:=c1+c2+a1*b1+a2*b1+a3*b2+a4*b2+a1*d1*b2+a2*d1*b2	
		"
	fit.sp <- sem(myModel.sp, data = scale(na.omit(sem.data[,-1])))

	par(mfrow=c(1,2))	
	semPaths(fit,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)#标准化的载荷std
	semPaths(fit.sp,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)
	savePlot("semplot",type=c("pdf"),device=dev.cur(),restoreConsole=TRUE)

	##2
	myModel<- " phylo~b1*dis.netdiv+b2*dis.ext+c1*clim.anc.m+c2*clim.pre.m+c3*g.anc.m+c4*g.pre.m
	dis.netdiv~a1*clim.anc.m+a2*clim.pre.m+a3*g.anc.m+a4*g.pre.m
	dis.ext~a5*clim.anc.m+a6*clim.pre.m+a7*g.anc.m+a8*g.pre.m+d1*dis.netdiv
	ind1:=a1*b1
	ind2:=a2*b1
	ind3:=a3*b2
	ind4:=a4*b2
	ind5:=a5*b1
	ind6:=a6*b1
	ind7:=a7*b2
	ind8:=a8*b2
	secondInd1:=a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2
	total1:=c1+c2+c3+c4+a1*b1+a2*b1+a3*b2+a4*b2+a5*b1+a6*b1+a7*b2+a8*b2+a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2	
		"		
	fit <- sem(myModel, data = scale(na.omit(sem.data[,-1])))
	summary (fit, rsquare=TRUE, standardized=TRUE, fit.measures=TRUE)
	
	myModel.sp<- " sp~b1*dis.netdiv+b2*dis.ext+c1*clim.anc.m+c2*clim.pre.m+c3*g.anc.m+c4*g.pre.m
	dis.netdiv~a1*clim.anc.m+a2*clim.pre.m+a3*g.anc.m+a4*g.pre.m
	dis.ext~a5*clim.anc.m+a6*clim.pre.m+a7*g.anc.m+a8*g.pre.m+d1*dis.netdiv
	ind1:=a1*b1
	ind2:=a2*b1
	ind3:=a3*b2
	ind4:=a4*b2
	ind5:=a5*b1
	ind6:=a6*b1
	ind7:=a7*b2
	ind8:=a8*b2
	secondInd1:=a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2
	total1:=c1+c2+c3+c4+a1*b1+a2*b1+a3*b2+a4*b2+a5*b1+a6*b1+a7*b2+a8*b2+a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2	
		"		
	fit.sp <- sem(myModel.sp, data = scale(na.omit(sem.data[,-1])))

	par(mfrow=c(1,2))	
	semPaths(fit,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)#标准化的载荷std
	semPaths(fit.sp,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)
	savePlot("semplotmean",type=c("pdf"),device=dev.cur(),restoreConsole=TRUE)
	
	##3
	myModel<- " phylo~b1*dis.netdiv+b2*dis.ext+c1*clim.anc.s+c2*clim.pre.s+c3*g.anc.s+c4*g.pre.s
	dis.netdiv~a1*clim.anc.s+a2*clim.pre.s+a3*g.anc.s+a4*g.pre.s
	dis.ext~a5*clim.anc.s+a6*clim.pre.s+a7*g.anc.s+a8*g.pre.s+d1*dis.netdiv
	ind1:=a1*b1
	ind2:=a2*b1
	ind3:=a3*b2
	ind4:=a4*b2
	ind5:=a5*b1
	ind6:=a6*b1
	ind7:=a7*b2
	ind8:=a8*b2
	secondInd1:=a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2
	total1:=c1+c2+c3+c4+a1*b1+a2*b1+a3*b2+a4*b2+a5*b1+a6*b1+a7*b2+a8*b2+a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2	
		"		
	fit <- sem(myModel, data = scale(na.omit(sem.data[,-1])))
	summary (fit, rsquare=TRUE, standardized=TRUE, fit.measures=TRUE)
	
	myModel.sp<- " sp~b1*dis.netdiv+b2*dis.ext+c1*clim.anc.s+c2*clim.pre.s+c3*g.anc.s+c4*g.pre.s
	dis.netdiv~a1*clim.anc.s+a2*clim.pre.s+a3*g.anc.s+a4*g.pre.s
	dis.ext~a5*clim.anc.s+a6*clim.pre.s+a7*g.anc.s+a8*g.pre.s+d1*dis.netdiv
	ind1:=a1*b1
	ind2:=a2*b1
	ind3:=a3*b2
	ind4:=a4*b2
	ind5:=a5*b1
	ind6:=a6*b1
	ind7:=a7*b2
	ind8:=a8*b2
	secondInd1:=a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2
	total1:=c1+c2+c3+c4+a1*b1+a2*b1+a3*b2+a4*b2+a5*b1+a6*b1+a7*b2+a8*b2+a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2	
		"		
	fit.sp <- sem(myModel.sp, data = scale(na.omit(sem.data[,-1])))

	par(mfrow=c(1,2))	
	semPaths(fit,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)#标准化的载荷std
	semPaths(fit.sp,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)
	savePlot("semplotsd",type=c("pdf"),device=dev.cur(),restoreConsole=TRUE)
	
	
	#地理区内外的
	region=as.character(unique(climnet$phname))		
	exFormat.within=c()
	for (i in 1:length(region)){		
		envir=climnet						
		envt=subset(envir,envir[,"phname"]==region[i])		
		physim=as.data.frame(phylo[match(envt$adcode,colnames(phylo)),match(envt$adcode,rownames(phylo))])
		Dissim=cbind(adcode=envt$adcode,physim)
		exFormat_pre <- formatsitepair(Dissim,3, XColumn="Lon", YColumn="Lat", predData=data.frame(envt,Adcode99=envt$adcode),siteColumn="adcode")						
		exFormat.within=rbind(exFormat.within,exFormat_pre)	
	}	
	
	past.env.within=get.pastenv(exFormat.within,pastdis.all)
	past.env.across=past.env[-match(past.env.within$ID,past.env$ID),]
	
	myModel<- " phylo~b1*dis.netdiv+b2*dis.ext+c1*dis.clim+c2*dis.geo
	dis.netdiv~a1*dis.clim+a2*dis.geo
	dis.ext~a3*dis.clim+a4*dis.geo+d1*dis.netdiv
	ind1:=a1*b1
	ind2:=a2*b1
	ind3:=a3*b2
	ind4:=a4*b2
	secondInd1:=a1*d1*b2+a2*d1*b2
	total1:=c1+c2+a1*b1+a2*b1+a3*b2+a4*b2+a1*d1*b2+a2*d1*b2	
		"		
	fit.within <- sem(myModel, data = scale(na.omit(past.env.within[,-1])))
	fit.across <- sem(myModel, data = scale(na.omit(past.env.across[,-1])))
	
	par(mfrow=c(1,2))	
	semPaths(fit.within,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)
	semPaths(fit.across,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)
	savePlot("semplot-witacros",type=c("pdf"),device=dev.cur(),restoreConsole=TRUE)
	
	##2
	myModel<- " phylo~b1*dis.netdiv+b2*dis.ext+c1*clim.anc.m+c2*clim.pre.m+c3*g.anc.m+c4*g.pre.m
	dis.netdiv~a1*clim.anc.m+a2*clim.pre.m+a3*g.anc.m+a4*g.pre.m
	dis.ext~a5*clim.anc.m+a6*clim.pre.m+a7*g.anc.m+a8*g.pre.m+d1*dis.netdiv
	ind1:=a1*b1
	ind2:=a2*b1
	ind3:=a3*b2
	ind4:=a4*b2
	ind5:=a5*b1
	ind6:=a6*b1
	ind7:=a7*b2
	ind8:=a8*b2
	secondInd1:=a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2
	total1:=c1+c2+c3+c4+a1*b1+a2*b1+a3*b2+a4*b2+a5*b1+a6*b1+a7*b2+a8*b2+a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2	
		"		
	fit.within <- sem(myModel, data = scale(na.omit(past.env.within[,-1])))
	fit.across <- sem(myModel, data = scale(na.omit(past.env.across[,-1])))
	
	par(mfrow=c(1,2))	
	semPaths(fit.within,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)
	semPaths(fit.across,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)
	savePlot("semplot-witacrosmean",type=c("pdf"),device=dev.cur(),restoreConsole=TRUE)
	
	##3
	myModel<- " phylo~b1*dis.netdiv+b2*dis.ext+c1*clim.anc.s+c2*clim.pre.s+c3*g.anc.s+c4*g.pre.s
	dis.netdiv~a1*clim.anc.s+a2*clim.pre.s+a3*g.anc.s+a4*g.pre.s
	dis.ext~a5*clim.anc.s+a6*clim.pre.s+a7*g.anc.s+a8*g.pre.s+d1*dis.netdiv
	ind1:=a1*b1
	ind2:=a2*b1
	ind3:=a3*b2
	ind4:=a4*b2
	ind5:=a5*b1
	ind6:=a6*b1
	ind7:=a7*b2
	ind8:=a8*b2
	secondInd1:=a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2
	total1:=c1+c2+c3+c4+a1*b1+a2*b1+a3*b2+a4*b2+a5*b1+a6*b1+a7*b2+a8*b2+a1*d1*b2+a2*d1*b2+a3*d1*b2+a4*d1*b2	
		"		
	fit.within <- sem(myModel, data = scale(na.omit(past.env.within[,-1])))
	fit.across <- sem(myModel, data = scale(na.omit(past.env.across[,-1])))
	
	par(mfrow=c(1,2))	
	semPaths(fit.within,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)
	semPaths(fit.across,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex = 1)
	savePlot("semplot-witacrossd",type=c("pdf"),device=dev.cur(),restoreConsole=TRUE)
	
	#各个生物地理区的
	node=colnames(climnet)[9:15]		
	past.env=vector("list",length(node));names(past.env)=paste("be_",node,sep="")	
	for (i in 1:length(node)){	
		envir=climnet				
		envt=subset(envir,envir[,node[i]]>0)				
		physim=as.data.frame(phylo[match(envt$adcode,colnames(phylo)),match(envt$adcode,rownames(phylo))])
		Dissim=cbind(adcode=envt$adcode,physim)
		exFormat_pre <- formatsitepair(Dissim,3, XColumn="Lon", YColumn="Lat", predData=data.frame(envt,Adcode99=envt$adcode),siteColumn="adcode")			
		
		cole=6+(length(colnames(exFormat_pre))-6)/2	
		cole2=length(colnames(exFormat_pre))
		exFormat.t=subset(exFormat_pre,
				  (exFormat_pre[,paste("s1.",node[i],sep="")]==1 & exFormat_pre[,paste("s2.",node[i],sep="")]==2)|
					(exFormat_pre[,paste("s1.",node[i],sep="")]==2 & exFormat_pre[,paste("s2.",node[i],sep="")]==1))			 
		past.env[[i]]=get.pastenv(exFormat.t,pastdis.all)	
	}
	
	par(mfrow=c(2,2))
	for (i in 1:4){	
	fit.t <- sem(myModel, data = scale(na.omit(past.env[[i]][,-1])))
	semPaths(fit.t,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex =1.5)	
	mtext(node[i])
	}
	savePlot("semplot-bioregion1",type=c("pdf"),device=dev.cur(),restoreConsole=TRUE)
	
	par(mfrow=c(2,2))
	for (i in 5:7){	
	fit.t <- sem(myModel, data = scale(na.omit(past.env[[i]][,-1])))
	semPaths(fit.t,what = "std",layout = "circle2",residuals=FALSE,edge.label.cex =1.5)	
	mtext(node[i])
	}
	savePlot("semplot-bioregion2",type=c("pdf"),device=dev.cur(),restoreConsole=TRUE)
	
##7.计算平均成对距离(mean pairwise distance, Dpw),Loic推荐的两个谱系距离算法.pdf,Ecology Letters, (2019) 22: 1126–1135;Phylogenetic turnover patterns consistent with niche conservatism in montane plant species
	source("read_Phylogenetic_Tree.R")
	#tre_ori=read.tree("raxml_v2_dated_Av1_final_clean_ann.full.polyRes.tre")#全球属级的树,140-210 constraint
	tre_ori=read.tree("raxml_v2_dated_Av2_final_clean_b_ann.full.polyRes.tre")#全球属级的树,140-150 constraint
	#tre_ori=read.tree("raxml_v2_dated_Av3_final_clean_b_ann.full.polyRes.tre")#全球属级的树,140-259 constraint
	
	an=read.csv("GyAn.csv")
	d=subset(tre_ori$tip.label,!tre_ori$tip.label%in%an$Genus)
	tre <- drop.tip(tre_ori, tip=d)	
	
	dis=read.csv("Spdat_an_isrm191016.csv")
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]

	#用嵌套循环计算beta，dpw,dnn,
	betaphylo2<-function(tre,spdis){
		x=seq(from=6,to=length(rownames(spdis)),by=20)#防断电失去所有
		result=matrix(data=0,nrow=(dim(spdis)[1]-1),ncol=(dim(spdis)[1]-1))
		colnames(result)=rownames(spdis)[2:dim(spdis)[1]]
		rownames(result)=rownames(spdis)[1:(dim(spdis)[1]-1)]
		
		result2=matrix(data=0,nrow=(dim(spdis)[1]-1),ncol=(dim(spdis)[1]-1))
		colnames(result2)=rownames(spdis)[2:dim(spdis)[1]]
		rownames(result2)=rownames(spdis)[1:(dim(spdis)[1]-1)]
		
		result3=matrix(data=0,nrow=(dim(spdis)[1]-1),ncol=(dim(spdis)[1]-1))
		colnames(result3)=rownames(spdis)[2:dim(spdis)[1]]
		rownames(result3)=rownames(spdis)[1:(dim(spdis)[1]-1)]
		
		result4=matrix(data=0,nrow=(dim(spdis)[1]-1),ncol=(dim(spdis)[1]-1))
		colnames(result4)=rownames(spdis)[2:dim(spdis)[1]]
		rownames(result4)=rownames(spdis)[1:(dim(spdis)[1]-1)]
	for (i in 1:(dim(spdis)[1]-1)){
	   print(paste("i=",i))
	   
	   for (j in (i+1):dim(spdis)[1]){	
		
		   aa=rbind(spdis[i,],spdis[j,])
		   rownames(aa)=c(rownames(spdis)[i],rownames(spdis)[j])
		   aa[,aa[1,]+aa[2,]==0] <- NA
		   bb=na.omit(t(aa))
		   
		   bbb=subset(bb,rownames(bb)%in%tre$tip.label)#去掉了分布列表中全球树没有的属名
		   bbb.i=subset(bbb, bbb[,1]==1)
		   bbb.j=subset(bbb, bbb[,2]==1)
		   
		   ## mean pairwise phylogenetic distance of a community
		   #pd = meanPhyloDist(rownames(bbb), tre, phylo.distance = NULL, method = "branch.length")
		   di <- phylo.dist(tre, tip1=rownames(bbb),method = "branch.length")
		  
		  #填充对角矩阵		
		  di[is.na(di)] <- 0		
		  di.c=di+t(di)
		  dis.i=0;dnn.i=0
		  for (m in 1: dim(bbb.i)[1]){
		  di.t=di.c[rownames(di.c)%in%rownames(bbb.i)[m],colnames(di.c)%in%rownames(bbb.j)]
		  dis.i=dis.i+mean(di.t)		  
		  dnn.i=dnn.i+min(di.t[di.t!=min(di.t)])
		  }	
		  
		  dis.j=0;dnn.j=0
		  for (n in 1: dim(bbb.j)[1]){
		  di.t2=di.c[rownames(di.c)%in%rownames(bbb.j)[n],colnames(di.c)%in%rownames(bbb.i)]
		  dis.j=dis.j+mean(di.t2)		  
		  dnn.j=dnn.j+min(di.t2[di.t2!=min(di.t2)])
		  }
		  
		   result[i,j-1] <- (dis.i+dis.j)/(dim(bbb.i)[1]+dim(bbb.j)[1])
		   result2[i,j-1] <- (dis.i+dis.j)/2
		   
		   result3[i,j-1] <- (dnn.i+dnn.j)/(dim(bbb.i)[1]+dim(bbb.j)[1])
		   result4[i,j-1] <- (dnn.i+dnn.j)/2
		   
		   print(paste("j=",j)) 		  
	   }
		#if (i%in%x){write.csv(result,"phylobeta02.csv")}
	}
	finre=list(dpwspweight=result,dpw=result2,dnnspweight=result3,dnn=result4)
	return(finre)	
	}
	
	dnndpw=betaphylo2(tre,spdis)
	saveRDS(dnndpw,"dnndpw.RDS")	
		

	
	
#######未放到文章中的分析##########
	
##检测各region的谱系集聚情况
	library(ape)
	#install.packages("PhyloMeasures")
	library(PhyloMeasures)
	dis=read.csv("Spdat_an_isrm.csv")	
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]		
	tre_ori=read.tree("140_150_constrained3_CL0.0001_full_tree_age.tre")#全球属级的树
	d=subset(tre_ori$tip.label,!tre_ori$tip.label%in%colnames(spdis))
	tre <- drop.tip(tre_ori, tip=d)	
	spdis=spdis[,match(tre$tip.label,colnames(spdis))]
	
	nodelist=read.csv("sos_node.csv")[,1:13]
	region=colnames(nodelist)[5:13]#谱系beta分区的边界名称（9个）
	
	dis.o=read.csv("Spdat_old.csv")
	spdis.o=as.matrix(dis.o[,2:dim(dis.o)[2]])
	rownames(spdis.o)=dis.o[,1]	
	
	dis.y=read.csv("Spdat_young.csv")
	spdis.y=as.matrix(dis.y[,2:dim(dis.y)[2]])
	rownames(spdis.y)=dis.y[,1]
	
	nri_mat<-function(region,nodelist,spdis,tre){
		all=c()
		for (i in 1:length(region)){		
			nodelist1=subset(nodelist,nodelist[[region[i]]]>0)		
			spdis.t0=spdis[match(nodelist1$adcode,as.integer(rownames(spdis))),]
			spdis.t=spdis.t0[,colSums(spdis.t0)>0]			
			rownames(spdis.t)=nodelist1[[region[i]]]	
			comm=t(sapply(by(spdis.t,rownames(spdis.t),colSums),identity))
			
			d.t=subset(tre$tip.label,!tre$tip.label%in%colnames(spdis.t))
			tre.t <- drop.tip(tre, tip=d.t)							
			nri <- -1*mpd.query(tre.t, comm, TRUE) 		
			
			# asem=colnames(spdis.t)#研究区域的物种名录	
			# u=NRI(asem,tre,method = "branch.length")#south的region1所有物种在谱系上是随机分布的（nri=0.4510512）
			# a=data.frame(region=paste(region[i],1,sep=""),NRI=u)
			
			# nodelist2=subset(nodelist,nodelist[[region[i]]]==2)
			# spdis.t02=spdis.y[match(nodelist2$adcode,as.integer(rownames(spdis.y))),]
			# spdis.t2=spdis.t02[,colSums(spdis.t02)>0]
			# asem2=colnames(spdis.t2)#研究区域的物种名录	
			# u2=NRI(asem2,tre,method = "branch.length")
			# b=data.frame(region=paste(region[i],2,sep=""),NRI=u2)
			# temp=rbind(a,b)
		
			all=rbind(all,nri)
			}
		rownames(all)=region
		colnames(all)=c("region1","region2")
		return(all)
		}
	young=nri_mat(region,nodelist,spdis.y,tre)
	old=nri_mat(region,nodelist,spdis.o,tre)
	all=nri_mat(region,nodelist,spdis,tre)
	re=cbind(all,young,old)
	colnames(re)=c("all1","all2","young1","young2","old1","old2")
	write.csv(re,"nri.csv")
	
	##10个区域总的nri
	spdis.t=spdis
	rownames(spdis.t)=nodelist$phyl
	comm.all=t(sapply(by(spdis.t,rownames(spdis.t),colSums),identity))
	nri <- -1*mpd.query(tre, comm.all, TRUE)
	names(nri)=rownames(comm.all)
	
	spdis.y.t=spdis.y
	rownames(spdis.y.t)=nodelist$phyl
	comm.y=t(sapply(by(spdis.y.t,rownames(spdis.y.t),colSums),identity))
	d.y=subset(tre$tip.label,!tre$tip.label%in%colnames(spdis.y))
	tre.y <- drop.tip(tre, tip=d.y)	
	nri.y <- -1*mpd.query(tre.y, comm.y, TRUE)
	names(nri.y)=rownames(comm.y)
	
	spdis.o.t=spdis.o
	rownames(spdis.o.t)=nodelist$phyl
	comm.o=t(sapply(by(spdis.o.t,rownames(spdis.o.t),colSums),identity))
	d.o=subset(tre$tip.label,!tre$tip.label%in%colnames(spdis.o))
	tre.o <- drop.tip(tre, tip=d.o)	
	nri.o <- -1*mpd.query(tre.o, comm.o, TRUE)
	names(nri.o)=rownames(comm.o)
	
	re2=rbind(nri,nri.y,nri.o)
	rownames(re2)=c("all","young","old")
	write.csv(re2,"nri-10.csv")
	
##探索谱系结构，IC，gamma test,	
	install.packages("CollessLike")
	install.packages("phytools")
	
	library(CollessLike)
	library(ape)	
	library(phytools)	
	dis=read.csv("Spdat_an_isrm.csv")	
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]			
	tre <- read.tree("tre_is_rm.tre")	
	spdis=spdis[,match(tre$tip.label,colnames(spdis))]	
	nodelist=read.csv("sos_node.csv")[,1:13]
	
	#探索9个region小树的谱系结构，IC，gamma test,
	region=colnames(nodelist)[5:13]#谱系beta分区的边界名称（9个）		
	IC=c()
	gama=c()
	gama.re=vector("list",9)
	names(gama.re)=region
	par(mfrow=c(3,3))	
	#node.n=c()
	for (i in 1:length(region)){		
		nodelist1=subset(nodelist,nodelist[[region[i]]]>0)		
		spdis.t0=spdis[match(nodelist1$adcode,as.integer(rownames(spdis))),]
		spdis.t=spdis.t0[,colSums(spdis.t0)>0]			
		d.t=subset(tre$tip.label,!tre$tip.label%in%colnames(spdis.t))
		tre.t <- drop.tip(tre, tip=d.t)	#9个region的小树
		#node.n=c(node.n,Ntip(tre.t)-1)
		
		#IC（Colless index，Colless (1982)，Biol. Rev. (2017), 92, pp. 698–715.）Package ‘CollessLike’ 
		t=cbind(region[i],colless.like.index(tre.t),colless.like.index(tre.t,norm=TRUE))
		IC=rbind(IC,t)
		
		#γ-test of Pybus & Harvey (2000).
		z<-ltt(tre.t,gamma=FALSE,plot=FALSE)
		g<-gammatest(z)
		gama=rbind(gama,g)
		plot(z,log="y",log.lineages=FALSE)		
		gama.re[[i]]=cbind(ltt=z$ltt,times=z$times)
	}
	
	rownames(gama)=region
	write.csv(gama,"gama.csv")
	for (i in 1:length(region)){
	write.csv(gama.re[[i]],paste("gama_",region[i],".csv",sep=""))}
	
	colnames(IC)=c("region","Ic","Ic_norm")
	write.csv(IC,"IC.csv")
	
	#探索10个全区的谱系结构，IC，gamma test
	area=unique(nodelist$phyl)	
	IC=c()
	gama=c()
	gama.re=vector("list",10)
	names(gama.re)=area
	par(mfrow=c(3,3))
	#node.n=c()
	for (i in 1:length(area)){		
		nodelist1=subset(nodelist,nodelist$phyl==area[i])		
		spdis.t0=spdis[match(nodelist1$adcode,as.integer(rownames(spdis))),]
		spdis.t=spdis.t0[,colSums(spdis.t0)>0]			
		d.t=subset(tre$tip.label,!tre$tip.label%in%colnames(spdis.t))
		tre.t <- drop.tip(tre, tip=d.t)	#10个area的小树
		#node.n=c(node.n,Ntip(tre.t)-1)
		t=cbind(area[i],colless.like.index(tre.t),colless.like.index(tre.t,norm=TRUE))
		IC=rbind(IC,t)		
		
		z<-ltt(tre.t,gamma=FALSE,plot=FALSE)
		g<-gammatest(z)
		gama=rbind(gama,g)
		plot(z,log="y",log.lineages=FALSE)		
		gama.re[[i]]=cbind(ltt=z$ltt,times=z$times)
	}
	
	rownames(gama)=area
	write.csv(gama,"a_gama.csv")
	for (i in 1:length(area)){
	write.csv(gama.re[[i]],paste("a_gama_",area[i],".csv",sep=""))}
	
	colnames(IC)=c("area","Ic","Ic_norm")
	write.csv(IC,"a_IC.csv")
	
	#基于ltt（）结果画node随时间的分化曲线(lenage through time)
	library(ggplot2)
	#9 region
	nodelist=read.csv("sos_node.csv")[,1:13]
	region=colnames(nodelist)[5:13]		
	rlt=vector("list",length(region))
	names(rlt)=region	
	for (i in 1:length(region)){	
		temp=read.csv(paste("gama_",region[i],".csv",sep=""))
		n <- max(temp$ltt)
		rlt[[i]]=temp[2:which(temp$ltt==n),]
	}
	
	temp=c();ltt.r=c()
	for (i in 1:length(region)){
		temp=cbind(region=region[i],rlt[[i]])
		ltt.r=rbind(ltt.r,temp)
	}
	p=ggplot(ltt.r, aes(x=times, y=log(ltt), color=region,group=region,linetype=region)) + geom_line(size=1)+theme(legend.key.width=unit(2,'cm'))
	p+theme(legend.key.width=unit(2,'cm'),
	legend.justification = c(0, 1), # pick the upper left corner of the legend box and
    legend.position = c(0, 1), # adjust the position of the corner as relative to axis
    legend.background = element_rect(fill = NA), # transparent legend background
    legend.box = "horizontal", # horizontal arrangement of multiple legends
    legend.spacing.x = unit(0.5, units = "cm"), # horizontal spacing between legends
    panel.grid = element_blank() # eliminate grids
    )
	
	#10 area
	area=c(1:10)
	name=c("NeoArc","EAsia","CAsia","Euro","NAfr","NeoTro","NeoTmp","AfrTro","Aus","IndMal")	
	alt=vector("list",length(area))
	names(alt)=name	
	for (i in 1:length(area)){
		temp=read.csv(paste("a_gama_",area[i],".csv",sep=""))
		n <- max(temp$ltt)
		alt[[i]]=temp[2:which(temp$ltt==n),]
	}
		
	temp=c();ltt.a=c()
	for (i in 1:length(area)){
		temp=cbind(area=names(alt)[i],alt[[i]])
		ltt.a=rbind(ltt.a,temp)
	}
	m=data.frame(hemi=c("Neo-North","Old-North","Old-North","Old-North","Old-North","Neo-South","Neo-South","Old-South","Old-South","Old-South"),name)
	ltt.a2=cbind(m[match(ltt.a$area,m$name),],ltt.a)	
	p=ggplot(ltt.a2, aes(x=times, y=log(ltt), color=area,group=area,linetype=area))+ geom_line(size=1)	
	p+theme(legend.key.width=unit(2,'cm'),
	legend.justification = c(0, 1), # pick the upper left corner of the legend box and
    legend.position = c(0, 1), # adjust the position of the corner as relative to axis
    legend.background = element_rect(fill = NA), # transparent legend background
    legend.box = "horizontal", # horizontal arrangement of multiple legends
    legend.spacing.x = unit(0.5, units = "cm"), # horizontal spacing between legends
    panel.grid = element_blank() # eliminate grids
    )
	
	# my_palette_2 <- colorRampPalette(colors = c("slateblue3","springgreen3","orange","pink2","red4"))(4)
	# p+geom_polygon(data = ltt.a2,aes(x    = times,y    = ltt,fill = hemi)+
	# scale_fill_manual(                   breaks = unique(hemi),
                    # labels = unique(hemi),
                    # values = unique(hemi))
					
	#寻找快速分化的点所在的节点编号，并匹配年龄和gnd值
	dis=read.csv("Spdat_an_isrm.csv")	
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]			
	tre <- read.tree("tre_is_rm.tre")	
	spdis=spdis[,match(tre$tip.label,colnames(spdis))]	
	
	load("D:\\bioregion\\data_analysis\\noderich.RData")	
	nodelist=read.csv("sos_node.csv")
	node=colnames(nodelist)[14:length(colnames(nodelist))]#GND node节点编号
	region=colnames(nodelist)[5:13]#谱系beta分区的边界名称（9个）
	
	IC=matrix(0,length(node),length(region))
	colnames(IC)=region
	rownames(IC)=node
	
	ICn=matrix(0,length(node),length(region))
	colnames(ICn)=region
	rownames(ICn)=node	
	
	ICnor=matrix(0,length(node),length(region))
	colnames(ICnor)=region
	rownames(ICnor)=node

	gamma=matrix(0,length(node),length(region))
	colnames(gamma)=region
	rownames(gamma)=node
	
	gammap=matrix(0,length(node),length(region))
	colnames(gammap)=region
	rownames(gammap)=node
	
	for (i in 1:length(region)){
		nodelist.t=subset(nodelist,nodelist[[region[i]]]>0)
		spdis.t0=spdis[match(nodelist.t$adcode,as.integer(rownames(spdis))),]
		spdis.t=spdis.t0[,colSums(spdis.t0)>0]
		for (j in 1:length(node)){
			n=as.integer(substr(node[j],2,nchar(node[j])))			
			for (k in 1:length(noderich)){
				if (noderich[[k]]$node.label[1]==n)	{	
					tre.t0=noderich[[k]]						
					d.t=subset(tre.t0$tip.label,!tre.t0$tip.label%in%colnames(spdis.t))
					if (length(d.t)!=Ntip(tre.t0)){
						tre.t <- drop.tip(tre.t0, tip=d.t)
						if (Ntip(tre.t)>1){	
							#Ic
							IC[j,i]=colless.like.index(tre.t)
							ICnor[j,i]=colless.like.index(tre.t,norm=TRUE)
							ICn[j,i]=IC[j,i]/(Ntip(tre.t)-1)
							
							#gamma
							z<-ltt(tre.t,gamma=FALSE,plot=FALSE)		
							gamma[j,i]=gammatest(z)[[1]]
							gammap[j,i]=gammatest(z)[[2]]
							
							print(paste(n,"node is at",k,"with",Ntip(tre.t),"tips"))
						}						
					}
					
				}
			}
									
		}
	}
	
	Ic_node=cbind(IC,ICn,ICnor)
	write.csv(Ic_node,"Ic_node.csv")#手动增加Iclog
	
	gammar_node=cbind(gamma,gammap)	
	write.csv(gammar_node,"gammar_node.csv")
	
	#匹配年龄和gnd值
	gamma=read.csv("gammar_node.csv")
	ic=read.csv("Ic_node.csv")
	rich=read.csv("node_r2_p_richpro.csv")
	richfin=cbind(rich,ic[match(rich$node,ic$node),],gamma[match(rich$node,gamma$node),])
	write.csv(richfin,"richfin.csv")	
	
	#筛选gnd贡献大的node并做散点图和小提琴图
	nodesub=list();vars=c();median=c()
	for (i in 1:length(region)){
		r=paste(region[i],"_r2",sep="")
		#ri=paste(region[i],"_rich",sep="")
		p=paste(region[i],"_p",sep="")
		thes.r=5
		#thes.ri=3#阈值类型，3-中位数，4-平均值，5-3/4分位数
		thes.p=0.05#显著性
		nodesub[[i]]=subset(richfin,richfin[[p]]<thes.p #)
			& richfin[[r]]>summary(richfin[[r]])[[thes.r]] )
			#& richfin[[ri]]>summary(richfin[[ri]])[[thes.ri]])#R2和Richness大于阈值的节点
		median.stat=data.frame(region=region[i],min_r2=summary(richfin[[r]])[[thes.r]],
			# min_rich_pro=summary(richfin[[ri]])[[thes.ri]],
			node_num=dim(nodesub[[i]])[1])		
		median=rbind(median,median.stat)		
	}
	
	names(nodesub)=region
	median#显示阈值
	
	#画小提琴图
	#install.packages("vioplot")
	library(vioplot)
	index2=c("nodeage","ic","iclog","gamma")
	width <- 100; height <- 40	
	n.col <- 1; n.row <- length(index2)
	mylayout <- layout(matrix(1:((2+n.col)*(2+n.row)), ncol=2+n.col, byrow=T), 
			width=c(0.7*width/(2+n.col), rep(width/(2+n.col),times = n.col), 0.3*width/(2+n.col)),
			height=c(0.3*height/(2+n.row), rep(height/(2+n.row),times=n.row), 0.7*height/(2+n.row)))
	layout.show((2+n.col)*(2+n.row))
		
	for (i in 1:length(index2)) {	 
		if (i == 1) {
			for (j in 1:(2+n.col)) {par(mar=c(0,0,0,0)); plot.new()}
			}	
			
		par(mar=c(0,0,0,0)); plot.new()			  
		par(mar=c(0.5,0.5,0,0), cex=1, cex.axis=0.8, cex.lab=1.2, mgp=c(2.4,0.3,0), tck=-0.04)
		cex.axis1=1.4;cex1=1.3;tck1=-0.02;cex.lab1=1.8;cex.legend = 1.4;cex.txt=1.4		  
			
			#画图
		if(i == 1){
			vioplot(na.omit(nodesub[[1]][[index2[i]]]),
			na.omit(nodesub[[2]][[index2[i]]]),
			na.omit(nodesub[[3]][[index2[i]]]),
			na.omit(nodesub[[4]][[index2[i]]]),
			na.omit(nodesub[[5]][[index2[i]]]),
			na.omit(nodesub[[6]][[index2[i]]]),
			na.omit(nodesub[[7]][[index2[i]]]),
			na.omit(nodesub[[8]][[index2[i]]]),
			na.omit(nodesub[[9]][[index2[i]]]),
			names=c(rep("",9)),col="gold")			
			abline(v=1.5,col="red",lty = 3,lwd=4)
			abline(v=5.5,col="red",lty = 3,lwd=4)
		}
		if((i >1)& (i< length(index2))){			
			vioplot(na.omit(nodesub[[1]][[paste("All",index2[i],sep="_")]]),
			na.omit(nodesub[[2]][[paste("South",index2[i],sep="_")]]),
			na.omit(nodesub[[3]][[paste("Ind",index2[i],sep="_")]]),
			na.omit(nodesub[[4]][[paste("Aus",index2[i],sep="_")]]),
			na.omit(nodesub[[5]][[paste("Ama",index2[i],sep="_")]]),
			na.omit(nodesub[[6]][[paste("North",index2[i],sep="_")]]),
			na.omit(nodesub[[7]][[paste("AsiaMed",index2[i],sep="_")]]),
			na.omit(nodesub[[8]][[paste("Med",index2[i],sep="_")]]),
			na.omit(nodesub[[9]][[paste("Asia",index2[i],sep="_")]]),
			names=c(rep("",9)),col="gold")			
			abline(v=1.5,col="red",lty = 3,lwd=4)
			abline(v=5.5,col="red",lty = 3,lwd=4) 
		} 		
		if (i == length(index2)) {
			vioplot(na.omit(nodesub[[1]][[paste("All",index2[i],sep="_")]]),
			na.omit(nodesub[[2]][[paste("South",index2[i],sep="_")]]),
			na.omit(nodesub[[3]][[paste("Ind",index2[i],sep="_")]]),
			na.omit(nodesub[[4]][[paste("Aus",index2[i],sep="_")]]),
			na.omit(nodesub[[5]][[paste("Ama",index2[i],sep="_")]]),
			na.omit(nodesub[[6]][[paste("North",index2[i],sep="_")]]),
			na.omit(nodesub[[7]][[paste("AsiaMed",index2[i],sep="_")]]),
			na.omit(nodesub[[8]][[paste("Med",index2[i],sep="_")]]),
			na.omit(nodesub[[9]][[paste("Asia",index2[i],sep="_")]]),
			names=names(nodesub),col="gold")
			abline(v=1.5,col="red",lty = 3,lwd=4)
			abline(v=5.5,col="red",lty = 3,lwd=4) 
			abline(h=-1.645,col="black",lty = 3,lwd=3)#The hypothesis of consistant b and d is rejected when gamma<-1.645
			}
			
		axis(side=2, label=T)
		mtext(side=2, text=index2[i], line=1.6, cex=1.5, las=0)	
		par(mar=c(0,0,0,0)); plot.new()		
	}
		
	#画散点图	
	nodelist=read.csv("sos_node.csv")
	region=colnames(nodelist)[5:13]#谱系beta分区的边界名称（9个）
	index=c("ic","iclog","gamma")
	
	width <- 30; height <- 40
	windows(width=width, height=height)
	n.col <- length(region); n.row <- length(index)
	mylayout <- layout(matrix(1:((2+n.col)*(2+n.row)), ncol=2+n.col, byrow=T), 
			width=c(0.7*width/(2+n.col), rep(width/(2+n.col),times = n.col), 0.3*width/(2+n.col)),
			height=c(0.3*height/(2+n.row), rep(height/(2+n.row),times=n.row), 0.7*height/(2+n.row)))
	layout.show((2+n.col)*(2+n.row))

	for (i in 1:length(index)) {	 
		if (i == 1) {
			for (j in 1:(2+n.col)) {par(mar=c(0,0,0,0)); plot.new()}
			}	
		for (xi in 1:length(region)) {
			name=paste(region[xi],index[i],sep="_")
			x=nodesub[[xi]][["nodeage"]]
			y=nodesub[[xi]][[name]]
			if (xi == 1) {par(mar=c(0,0,0,0)); plot.new()}			  
			par(mar=c(0.5,0.5,0,0), cex=1, cex.axis=0.8, cex.lab=1.2, mgp=c(2.4,0.3,0), tck=-0.04)
			cex.axis1=1.4;cex1=1.3;tck1=-0.02;cex.lab1=1.8;cex.legend = 1.4;cex.txt=1.4		  
			
			#画图
			plot(y~x, data=nodesub[[xi]], pch=19, xlab="", ylab="", cex=1, col='grey', cex.axis=cex.axis1, axes=F); box()			
			m <- glm(y~x, data=nodesub[[xi]])
			a <- seq(min(x), max(x), 0.1)
			lines(x=a, y=predict(m, newdata=list(x=a)), col="red", lty=2, lwd=1) 
			rug(jitter(x),lwd=1.6,col="blue")			
			#画图结束
			
			if (xi == 1) {
				axis(side=2, label=T)
				mtext(side=2, text=index[i], line=1.6, cex=1.5, las=0)				
			}
			
			if (i == length(index)) {
				axis(side=1, lwd = 0.5)
				mtext(side=1, text=paste(region[xi],"_age",sep=""), line=1.6, cex=1, las=0)
			}			
			
			if (xi == length(region)) {par(mar=c(0,0,0,0)); plot.new()}
		}
	}	
	
##环境因子对分区的影响	
#气候因子的相关系数
	library(corrgram)
	richclim=read.csv("richclim191024.csv")
	colin<-cor(richclim[15:27])
	corrgram(colin,order=TRUE,lower.panel=panel.shade,upper.panel=panel.pie,text.panel=panel.txt,main="cor")
	write.csv(colin,"colin.csv")
 #随机森林，Confusion matrix 反映To assess whether the biomes identified could be distinguished using climatic data
	library(randomForest)#/安装随机森林程序包/
	clim=read.csv("Clim.csv")
	region=read.csv("sos_node.csv")[,1:13]
	region2=cbind(region[,c(1,3:13)],phyl=paste("r",region$phyl,sep=""))
	richclim=cbind(region2,clim[match(region2$adcode,clim$ADCODE99),])[,c(1:13,15:27)]
	write.csv(richclim,"richclim.csv")#加了一列phname,温度和降水数据以及ART在源数据基础上除以10，TSN和PSN不变。	
	
	richclim=read.csv("richclim.csv")
	
	rf <- randomForest(phname~.,data=richclim[,14:27],importance=TRUE,proximity=TRUE)
	con=rf$confusion
	imp=round(importance(rf), 2)#the first nclass columns are the class-specific measures computed as mean descrease in accuracy. 
								#The nclass + 1st column is the mean descrease in accuracy over all classes. 
								#The last column is the mean decrease in Gini index. 
	write.csv(con,"confusion.csv")
	write.csv(imp,"importance.csv")	
 #NMDS
	MDSplot(rf,richclim$phname,palette=rainbow(10),pch=as.numeric(richclim$phname))
	#install.packages("ggplot2")
	library(ggplot2)
	rf.mds <- cmdscale(1 - rf$proximity, eig=TRUE)$points
	colnames(rf.mds)=c("MDS_1","MDS_2")
	region=data.frame(region=names(rf$y),num=as.character(rf$y))
	mds=cbind(rf.mds,region[match(rownames(rf.mds),region$region),])
	p=ggplot(mds, aes(x = MDS_1, y = MDS_2, colour = num))
	p+geom_point(size=3,alpha = .7)
	
  #散点图展示气候空间
	# xlist= colnames(richclim)[15:27]
	# for (i in 1:(length(xlist)-1)) {
		# for (j in (i+1):length(xlist)){
			# p=ggplot(richclim, aes(richclim[[xlist[i]]], richclim[[xlist[j]]] , colour = phname))+geom_point(size=3,alpha = .7)+labs(x = xlist[i],y=xlist[j])
			# p
			# name=paste(xlist[i],"_",xlist[j],sep="")
			# ggsave(paste(name, ".jpg", sep = ""),plot=p)	
		# }
	# }
	
	# plotbiomes是一个基于ggplot的可以画出用降雨、温度划分的生物群系的包。
	# https://rdrr.io/github/valentinitnelav/plotbiomes/f/README.md
	# install.packages("devtools")	
	# install_github("valentinitnelav/plotbiomes")	
	library(gridExtra)
	library(devtools)
	library(plotbiomes)
	library(ggplot2)
	library(RColorBrewer)
	#install.packages("ggsignif")
	library(ggsignif)
	#install.packages("cowplot")
	library(cowplot)		
	finalfig=function(richclim,i,ma){
		plotclim <- function(richclim,i){
			reg=colnames(richclim)[4:12]
			richclim2=subset(richclim,richclim[[reg[i]]]!=0)
			Subregion=richclim2[[reg[i]]]
			p=ggplot()+
			 geom_polygon(data = Whittaker_biomes,aes(x= temp_c,y = precp_cm,fill = biome), colour = "gray98",
					   size   = 1,alpha  = 0.5,show.legend = FALSE)+
			 theme(axis.title.x=element_blank(),axis.title.y=element_blank())+
			 geom_point(data = richclim2,aes(x = MAT,y = MAP,colour =Subregion,group =Subregion), shape = 16, size  = 3,show.legend = FALSE)#+
			 # scale_colour_hue(name=reg[i])+			
			 # theme(legend.justification = c(0, 1), # pick the upper left corner of the legend box and
					# legend.position = c(0, 1), # adjust the position of the corner as relative to axis
					# legend.background = element_rect(fill = NA), # transparent legend background
					# legend.box = "horizontal", # horizontal arrangement of multiple legends
					# legend.spacing.x = unit(0.5, units = "cm"), # horizontal spacing between legends
					# panel.grid = element_blank(),
					# legend.title=element_text(size=rel(1)),
					# legend.text=element_text(size=rel(1)))			
			return(p)
		}
		
		boxplotfuc=function(richclim,i,ma,lab=TRUE){	
		  reg=colnames(richclim)[4:12]
		  tmp=richclim0[richclim0[[reg[i]]]!=0,c("MAT","MAP",reg[i])]  
		  compaired=as.character(unique(tmp[[reg[i]]]))
		  p<-ggplot(tmp, aes(tmp[[reg[i]]],tmp[[ma]], fill = tmp[[reg[i]]]))+
		  geom_boxplot(width=0.6,position = position_dodge(0.8),outlier.size = 0,outlier.color = "white",show.legend = FALSE)+  
		  xlab(reg[i])+
		  ylab(ma)+
		  geom_signif(comparisons =list(compaired),vjust = 1.2,map_signif_level = T,test = wilcox.test) +
		  theme_bw()+
		  theme(
			axis.title.x = element_text(face = "bold",size=rel(1)),
			axis.title.y = element_text(face = "bold",size=rel(1)),
			axis.text.x=element_text(size=rel(0.9),face="bold"),
			axis.text.y=element_text(size=rel(0.6),face="bold"),
			axis.line.x = element_line(size = 0.5, colour = "black"),
			axis.line.y = element_line(size = 0.5, colour = "black"),
			panel.border = element_blank(),
			panel.grid = element_blank()
		  )+  guides(color=FALSE)
		  if (lab==FALSE) p=p+ scale_x_discrete(breaks=NULL) 
		 return(p)
		}

		plot0=plotclim(richclim0,i)
		p1=boxplotfuc(richclim0,i,ma[1],FALSE)
		p2=boxplotfuc(richclim0,i,ma[2])
		p=ggdraw() +draw_plot(plot0,0,0,1,1)+draw_plot(p1,0.093,0.65,0.3,0.35)+draw_plot(p2,0.09,0.42,0.3,0.35)
		return(p)
	}
	
	richclim0=read.csv("richclim_plot.csv")#将9个region栏的1和2替换为下级region的名字，如对于All的1，改为North
	plot1=finalfig(richclim0,1,c("MAT","MAP"))
	plot2=finalfig(richclim0,2,c("MAT","MAP"))
	plot3=finalfig(richclim0,3,c("MAT","MAP"))
	plot4=finalfig(richclim0,4,c("MAT","MAP"))
	plot5=finalfig(richclim0,5,c("MAT","MAP"))
	plot6=finalfig(richclim0,6,c("MAT","MAP"))
	plot7=finalfig(richclim0,7,c("MAT","MAP"))
	plot8=finalfig(richclim0,8,c("MAT","MAP"))
	plot9=finalfig(richclim0,9,c("MAT","MAP"))	
	grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,nrow=3)	
	
	#气候与SOS值的关系
	##glm提取R2和p值
	ss.glm <- function(r.glm)
                {
                r.ss <- summary(r.glm)
                rsq <- 100*(r.ss$null.deviance-r.ss$deviance)/r.ss$null.deviance
                adj.rsq <- 100*(1-(r.ss$deviance/r.ss$df.residual)/(r.ss$null.deviance/r.ss$df.null))
                f.stat <- ((r.ss$null.deviance-r.ss$deviance)/(r.ss$df.null-
                        r.ss$df.residual))/(r.ss$deviance/r.ss$df.residual)
                p <- pf(f.stat, r.ss$df.null-r.ss$df.residual, r.ss$df.residual, lower.tail=FALSE)

                return(c(r2=rsq,adj.r2=adj.rsq,p=p))
                }
	
	load("richfin.RData")		
	# nodesub=vector("list",length(region))
	# for (i in 1:length(region)){
		# tmp=na.omit(richfin[[i]])		
		# tmp=tmp[tmp$r2>0&tmp$p<0.05,]
		# tmp2=cbind(tmp,region=region[i])		
		# #noder=rbind(noder,tmp2)
		# nodesub[[i]]=tmp2[tmp2$r2>summary(tmp2$r2)[3],]
		# # noder4=rbind(noder4,tmp2[tmp2$r2>summary(tmp2$r2)[4],])
		# # noder5=rbind(noder5,tmp2[tmp2$r2>summary(tmp2$r2)[5],])#3-中位数；4-mean; 5-3/4分位数
		# }	
	
	clim=read.csv("richclim191024.csv")		
	sos.tar=vector("list",length(region))#各区域环境因子与sos值的列表	
	for (i in 1:length(region)){			
		sos.t=sos.re[[i]][2:dim(sos.re[[i]])[1],]
		sos.tar[[i]]=cbind(clim[match(rownames(sos.t),clim$adcode),c(7:20)],regionlist[match(rownames(sos.t),regionlist$adcode_ave2),c("adcode_ave2",region[i])],sos.t)		
		}
		
	region.r=vector("list",length(region))
	sos.r=vector("list",length(region))	
	names(sos.r)=region
	names(region.r)=region
	for (i in 1:length(region)){	
		#glm r2
		#region与气候因子的关系		
		tmp.r=c()
		for (r in 3:15){
			glm.r<-glm((sos.tar[[i]][,2]-1)~sos.tar[[i]][,r],family=binomial)
			tmp=ss.glm(glm.r)[2:3]			
			tmp.r=rbind(tmp.r,tmp)
		}
		vars.r=data.frame(var=colnames(sos.tar[[i]])[3:15])
		region.r[[i]]=cbind(vars.r,tmp.r)	
			
		#Sos node与气候因子的关系
		tmp.j=c()
		for (j in 16:dim(sos.tar[[i]])[2]){				
				tmp.m=c()
				for (m in 3:15){
					glm.n<-glm(sos.tar[[i]][,j]~sos.tar[[i]][,m])			
					tmp=ss.glm(glm.n)[2:3]
					tmp.m=rbind(tmp.m,tmp)
				}
				vars=data.frame(var=colnames(sos.tar[[i]])[3:15])
				nodename=data.frame(nname=rep(colnames(sos.tar[[i]])[j],dim(tmp.m)[1]))
				tmp2=cbind(vars,nodename,tmp.m)
				tmp.j=rbind(tmp.j,tmp2)						
		}
		sos.r[[i]]=tmp.j
	}	
	save.image("clim-sos.RData")
	
	#画图
	load("clim-sos.RData")	
	library(ggplot2)	
	noder.n=c();tmp=c()
	for (i in 1:length(region)){
		tmp=cbind(sos.r[[i]],region[i])
		colnames(tmp)[5]="region"
		noder.n=rbind(noder.n,tmp)
		}
	noder.n=na.omit(noder.n[noder.n$p<0.05,])	#Ama的节点8686只有2个sos值，因此是NA
	group <- data.frame(var=unique(noder.n$var),group=c(rep("tmp",7),rep("prep",4),rep("other",2)))
	noder2=cbind(noder.n,group=group[match(noder.n$var,group$var),2])
	noder2$var = factor(noder2$var, levels=c("MAT","PET","MDR","TSN","MTWQ","MTCQ","ART","MAP","PSN","MPWQ","MPDQ","Sard","NSoilType"))#调整x轴顺序
	mycolor=c("#D39200","#00B9E3","#F8766D" )
	p=ggplot(data = noder2,mapping = aes(x=var,y=adj.r2))+ 
		geom_boxplot(aes(fill=group),show.legend=FALSE)+
		scale_fill_manual(values = mycolor)+#自定义颜色
		theme(axis.title.x = element_blank(),
		axis.text.x= element_text(angle=45,size=8))+
		geom_hline(aes(yintercept=0),size=1,colour="black", linetype="dashed")+
		facet_wrap(~region,nrow = 3)+
		labs(y="glm R2adj for high GND node")#sos-clim.pdf
	
	noder=c();tmp=c()
	for (i in 1:length(region)){
		tmp=cbind(region.r[[i]],region[i])
		colnames(tmp)[4]="region"
		noder=rbind(noder,tmp)
		}
	#noder[noder$p>0.05,2]=0
	noder=noder[noder$p<0.05,]
	group <- data.frame(var=unique(noder.n$var),group=c(rep("tmp",7),rep("prep",4),rep("other",2)))
	noder2=cbind(noder,group=group[match(noder$var,group$var),2])
	noder2$var = factor(noder2$var, levels=c("MAT","PET","MDR","TSN","MTWQ","MTCQ","ART","MAP","PSN","MPWQ","MPDQ","Sard","NSoilType"))#调整x轴顺序
	mycolor=c("#D39200","#00B9E3","#F8766D" )
	p=ggplot(data = noder2,mapping = aes(x=var,y=adj.r2))+ 
	geom_bar(aes(fill=group),position=position_dodge(), stat="identity",show.legend=FALSE)+	
			scale_fill_manual(values = mycolor)+#自定义颜色
			theme(axis.title.x = element_blank(),
			axis.text.x= element_text(angle=45,size=8))+
			geom_hline(aes(yintercept=0),size=1,colour="black", linetype="dashed")+
		facet_wrap(~region,nrow = 3)+
		labs(y="glm R2adj")	 #region-clim.pdf
	
	##物种矩阵与气候矩阵的回归	
	clim=read.csv("richclim.csv")
	env_texture_sub <- vegdist(sub_env$Soil_texture,method='euclidean')
env_texture_sub_m <- as.matrix(env_texture_sub)
View(env_texture_sub_m)

mantel(Bact_sub_bray,env_Tanomaly_m,permutations=999)



	Bact_sub_bray <- read.csv("Bacteria_subsoil_BrayCurtis.csv")
Bact_sub_bray_c <- as.character()

for (i in 1:nrow(Bact_sub_bray)){
 for (j in 1:ncol(Bact_sub_bray)){
 if (i < j){
Bact_sub_bray_c[(i-1)*ncol(Bact_sub_bray)+j] <- Bact_sub_bray[i,j]
}
}}
View(Bact_sub_bray_c)
write.csv(Bact_sub_bray_c,file="Bact_sub_bray_c.csv")


		
# ##祖先分布区重建，通过以下代码，已经得到分析所需的.trees以及Rasp_dis文件，接下来可以按照A Rough guide to RASP进行分析。
	# 祖先分布区重建的模型可以用于不同生物区系之间的扩散。祖先重构是基于树文件的，然而如果范围扩大到整个被子植物的范围，
	# 其系统发育关系可能已经不是单纯的树结构，而是更加复杂的网络结构。即便假设他们之间的关系就是树，也会由于物种缺失太多导致可信度降低。
	# RASP目前只支持小于512个物种的推演，BioGeoBEARS包实现的是DIVA的最大似然估计版本，虽然并未限制物种数量，
	# 但实际上超过100个物种就要花费十几小时甚至数天才能出结果，更多的物种会导致运算时间呈指数增长。（来自余岩的回复）
	
	##提取每个节点所在的谱系树rasp_tre_xx.tre，并用RASP软件获得RASP_xx_condensed.trees以及有ID的RASP_xx.csv
	library(ape)
	dis=read.csv("Spdat_an_isrm.csv")	
	spdis=as.matrix(dis[,2:dim(dis)[2]])
	rownames(spdis)=dis[,1]		
	tre_ori=read.tree("140_150_constrained3_CL0.0001_full_tree_age.tre")#全球属级的树
	d=subset(tre_ori$tip.label,!tre_ori$tip.label%in%colnames(spdis))
	tre <- drop.tip(tre_ori, tip=d)	
	write.tree(tre,"tre_is_rm.tre")
	
	spdis=spdis[,match(tre$tip.label,colnames(spdis))]	
	nodelist=read.csv("sos_node.csv")[,1:13]
	region=colnames(nodelist)[5:13]#谱系beta分区的边界名称（9个）	
	
		 for (i in 1:length(region)){		
			 nodelist1=subset(nodelist,nodelist[[region[i]]]>0)		
			 spdis.t0=spdis[match(nodelist1$adcode,as.integer(rownames(spdis))),]
			 spdis.t=spdis.t0[,colSums(spdis.t0)>0]			
			
			 d.t=subset(tre$tip.label,!tre$tip.label%in%colnames(spdis.t))
			 tre.t <- drop.tip(tre, tip=d.t)	
			 print(tre.t)
			 write.tree(tre.t,paste("rasp_tre_",region[i],".tre",sep=""))
			 }
	
	##根据RASP_xx.csv及RASP的格式要求补全分布信息,得到rasp_dis_xx.csv
	spdis=unique(read.csv("Woody_SpDistribution.csv")[,c("Adcode99","Species_E1")])
	splist=read.csv("RASP_All.csv",header=F)[,1:2]
	colnames(splist)=c("ID","sp")
	region=read.csv("sos_node.csv")[,c(1:2,5:13)]
	dis=na.omit(cbind(spdis,region[match(spdis$Adcode99,region$adcode),]))#去除了岛屿和被子植物
	dis_rasp=cbind(splist[match(dis$Species_E1,splist$sp),],dis)
	write.csv(dis_rasp,"RASP_dis.csv")#之后将数字替换成大写字母,0替换为NA,ID列改名为All_ID,All_sp	
	
	raspdis<-function(rasp,r){
		rasp.t=na.omit(unique(rasp[,c(paste(r,"_ID",sep=""),paste(r,"_sp",sep=""),r)]))
		rasp_north=cbind(rasp.t,dis=rasp.t[,r])
		levels <- levels(rasp_north[,3])
		levels[length(levels) + 1] <- "AB"
		rasp_north[,4] <- factor(rasp_north[,4], levels = levels)
		for (i in 2:dim(rasp_north)[1])
		{
			if (rasp_north[i-1,1]==rasp_north[i,1]) {
			rasp_north[i-1,4]="AB"
			rasp_north[i,4]="AB"	}		
		}
		return(write.csv(unique(rasp_north[,c(1:2,4)]),paste("rasp_dis_",r,".csv",sep="")))		
	}
	
	# library(BioGeoBEARS)
	# library(ape)
	# install.packages("snow")
	# library(snow)
	# tre=read.tree("rasp_tre_Ama.tre")
	# dis=read.csv("rasp_dis_Ama.csv")
	# bears_2param_DIVA_fast(tre,dis)
	
	
	#匹配其他节点的ID，使分布信息的物种ID与.trees文件中物种ID一致
	rasp.t=read.csv("RASP_dis.csv")
	rasp=rasp.t
	for (i in 2:length(region)){
		r=region[i]
		splist=read.csv(paste("RASP_",r,".csv",sep=""),header=F)[,1:2]
		colnames(splist)=c(paste(r,"_ID",sep=""),paste(r,"_sp",sep=""))
		rasp=cbind(rasp,splist[match(rasp$All_sp,splist[,paste(r,"_sp",sep="")]),])
	}
	write.csv(rasp,"rasp_dis_fin.csv")
	
	for (i in 1:length(region)){
		r=region[i]
		raspdis(rasp,r)
	}
	
##GeoHiSSE,根据diversetree.pdf,模型构建依据GeoSSE2.pdf
	#install.packages("diversitree")
	library(diversitree)	
	geosse.chose<-function(reg){
		tree=read.tree(paste(reg,".newick",sep=""))#将之前生成的各节点.tre文件导入FigTree中，
												#然后Trees > Transform Branches > Select 'Cladogram'，
												#然后File>Export trees,选择newick及两个复选框得到
		dat=read.csv(paste(reg,"-states.csv",sep=""), header=F)#之前生成的rasp_dis_All.csv文件中添加一列，A-1，B-2，AB-0;
		# A species occurring in both biomes
# (widespread species) may diversify and give rise to either one
# endemic plus one widespread daughter species (rate l) or to two
# endemic daughter species, one in each biome (here referred to as
# speciation by biome divergence, rate lTempTrop). Speciation by
# biome divergence can occur if populations belonging to each
# biome experience directional selection in opposite directions
# leading to speciation
		states <- dat[, 4]
		names(states) <- dat[, 2]

		p <- starting.point.geosse(tree)
		
		statecols <- c("AB"="purple", "A"="blue", "B"="red")
		plot(tree, tip.color=statecols[states+1], cex=0.1,type="fan")
		legend("topright", names(statecols), col=statecols, lty=2)
		savePlot(paste(reg,"_plot",sep=""),type=c("pdf"),device=dev.cur(),restoreConsole=TRUE)
		
		model.name=c("full","no.sAB_eq.div_eq.d","no.sAB_eq.ex_eq.d","no.sAB_eq.sp_eq.d","no.sAB_eq.d",
			"no.sAB_eq.div","no.sAB_eq.ex","no.sAB_eq.sp","no.sAB",
			"eq.div_eq.d","eq.ex_eq.d","eq.sp_eq.d","eq.d",
			"eq.div","eq.ex","eq.sp")
			
		lik=vector("list",16)
		names(lik)=model.name

		lik[[1]]=make.geosse(tree,states) #full
		lik[[2]] <- constrain(lik[[1]], sAB ~ 0,sA ~ sB, xA ~ xB,dA ~ dB)#no.sAB_eq.div_eq.d
		lik[[3]] <- constrain(lik[[1]], sAB ~ 0,xA ~ xB,dA ~ dB)#no.sAB_eq.ex_eq.d
		lik[[4]] <- constrain(lik[[1]], sAB ~ 0,sA ~ sB,dA ~ dB)#no.sAB_eq.sp_eq.d
		lik[[5]] <- constrain(lik[[1]], sAB ~ 0,dA ~ dB)#no.sAB_eq.d

		lik[[6]] <- constrain(lik[[1]], sAB ~ 0,sA ~ sB, xA ~ xB)#no.sAB_eq.div
		lik[[7]] <- constrain(lik[[1]], sAB ~ 0,xA ~ xB)#no.sAB_eq.ex
		lik[[8]] <- constrain(lik[[1]], sAB ~ 0,sA ~ sB)#no.sAB_eq.sp
		lik[[9]] <- constrain(lik[[1]], sAB ~ 0)#no.sAB

		lik[[10]] <- constrain(lik[[1]], sA ~ sB, xA ~ xB,dA ~ dB)#eq.div_eq.d
		lik[[11]] <- constrain(lik[[1]], xA ~ xB,dA ~ dB)#eq.ex_eq.d
		lik[[12]] <- constrain(lik[[1]], sA ~ sB,dA ~ dB)#eq.sp_eq.d
		lik[[13]] <- constrain(lik[[1]], dA ~ dB)#eq.d

		lik[[14]] <- constrain(lik[[1]], sA ~ sB, xA ~ xB)#eq.div
		lik[[15]] <- constrain(lik[[1]], xA ~ xB)#eq.ex
		lik[[16]] <- constrain(lik[[1]], sA ~ sB)#eq.sp

		ml=vector("list",16)
		names(ml)=model.name
		ml[[1]]=find.mle(lik[[1]], p)
		p <- coef(ml[[1]])
		para=p

		for (i in 2:16){
			ml[[i]]=find.mle(lik[[i]], p[argnames(lik[[i]])])
			para=rbind(para,coef(ml[[i]], TRUE))
		}
		rownames(para)=model.name		
		write.csv(para,paste("para-",reg,".csv",sep=""))	
		
		# formu="ml[[1]]"
		# for (i in 2:16){
		# formu=paste(formu,paste(names(ml[i]),"=ml[[",i,"]]",sep=""),sep=",")
		# }
		# formu 
		ano=anova(ml[[1]],no.sAB_eq.div_eq.d=ml[[2]],no.sAB_eq.ex_eq.d=ml[[3]],no.sAB_eq.sp_eq.d=ml[[4]],
			no.sAB_eq.d=ml[[5]],no.sAB_eq.div=ml[[6]],no.sAB_eq.ex=ml[[7]],no.sAB_eq.sp=ml[[8]],
			no.sAB=ml[[9]],eq.div_eq.d=ml[[10]],eq.ex_eq.d=ml[[11]],eq.sp_eq.d=ml[[12]],eq.d=ml[[13]],
			eq.div=ml[[14]],eq.ex=ml[[15]],eq.sp=ml[[16]])
		write.csv(ano,paste("ano-",reg,".csv",sep=""))
		
		best.model=rownames(subset(ano,ano$AIC==min(ano$AIC)))
		result=vector("list",4)
		result[[1]]=lik[[best.model]]
		result[[2]]=ml[[best.model]]
		result[[3]]=lik
		result[[4]]=ml
		return(result)
	}
	
	nodelist=read.csv("sos_node.csv")
	region=colnames(nodelist)[5:13]#谱系beta分区的边界名称（9个）
	lik=vector("list",9)
	names(lik)=region
	ml=vector("list",9)
	names(ml)=region
	re.geosse=vector("list",9)
	names(re.geosse)=region
	
	for (i in 1:length(region)){	
		reg=region[i]
		re.geosse[[i]]=geosse.chose(reg)
		lik[[i]]=re.geosse[[i]][[1]]
		ml[[i]]=re.geosse[[i]][[2]]	
		save.image("geosse.RData")
		#load("geosse.RData")
	}
		
	for (i in 1:length(region)){
		lik1=lik[[i]]
		ml1=ml[[i]]
		p <- coef(ml1)
		prior <- make.prior.exponential(1/2)
		set.seed(1)
		tmp <- mcmc(lik1, p, nsteps=100, prior=prior, w=1, print.every=0)
		w <- diff(sapply(tmp[2:(dim(tmp)[2]-1)], quantile, c(0.025, 0.975)))
		#mcmc2 <- mcmc(lik1, p, nsteps=1000, prior=prior, w=w)
		
		# tmp0=read.csv(paste("mcmc-",region[i],".csv",sep=""))
		# tmp=tmp0[,-match(c("nA","nB"),colnames(tmp0))]
		# w <- diff(sapply(tmp[3:(dim(tmp)[2]-1)], quantile, c(0.025, 0.975)))		
		mcmc2 <- mcmc(lik1, p, nsteps=10000, prior=prior, w=w)		
		
		write.csv(mcmc2,paste("mcmc2-",region[i],".csv",sep=""))		
		save.image("geosse2.RData")
	}
	
	#画散点图
	library(diversitree)
	index=c("Speciation","Extinction","Net Diversification","Dispersal")
	nodelist=read.csv("sos_node.csv")
	region=colnames(nodelist)[5:13]#谱系beta分区的边界名称（9个）
	mcmc=vector("list",9)
	names(mcmc)=region
	clnm=vector("list",9)
	names(clnm)=region
	for (i in 1:length(region)){
	mcmc[[i]]=read.csv(paste("mcmc-",region[i],".csv",sep=""))#添加nA=sA-xA;nB=sB-xB
	clnm[[i]]=colnames(mcmc[[i]])
	}	
	rgnm=read.csv("rgnm.csv")[,2:10]
	
	width <- 80; height <- 40
	windows(width=width, height=height)
	n.col <- length(region); n.row <- length(index)
	mylayout <- layout(matrix(1:((2+n.col)*(2+n.row)), ncol=2+n.col, byrow=T), 
			width=c(0.7*width/(2+n.col), rep(width/(2+n.col),times = n.col), 0.3*width/(2+n.col)),
			height=c(0.3*height/(2+n.row), rep(height/(2+n.row),times=n.row), 0.7*height/(2+n.row)))
	layout.show((2+n.col)*(2+n.row))

	for (i in 1:length(index)) {	 
		if (i == 1) {
			for (j in 1:(2+n.col)) {par(mar=c(0,0,0,0)); plot.new()}
			}	
		for (xi in 1:length(region)) {
			
			if (xi == 1) {par(mar=c(0,0,0,0)); plot.new()}			  
			par(mar=c(1.3,1.3,0,0), cex=1, cex.axis=0.7, cex.lab=1.2, mgp=c(2.4,0.3,0), tck=-0.04)
			cex.axis1=1.4;cex1=1.3;tck1=-0.02;cex.lab1=1.8;cex.legend = 1.4;cex.txt=1.4		  
			
			#画图
			var.sp=c();var.ex=c();var.n=c();var.d=c()
			for (j in 3:(dim(mcmc[[xi]])[2]-1)){
			if (substr(colnames(mcmc[[xi]])[j], 1,1)=="s") var.sp=c(var.sp,colnames(mcmc[[xi]])[j])
			if (substr(colnames(mcmc[[xi]])[j], 1,1)=="x") var.ex=c(var.ex,colnames(mcmc[[xi]])[j])
			if (substr(colnames(mcmc[[xi]])[j], 1,1)=="n") var.n=c(var.n,colnames(mcmc[[xi]])[j])
			if (substr(colnames(mcmc[[xi]])[j], 1,1)=="d") var.d=c(var.d,colnames(mcmc[[xi]])[j])
			}
			
			if (i==1){
				rownames(rgnm)=c("sA","sB","sAB")
				if(length(var.sp)==1) {
					if(nchar(var.sp)==2){col="gray";names(var.sp)="Equal"} 
					if(nchar(var.sp)==3){col="blue";names(var.sp)="sAB"} 
					}
				if(length(var.sp)==2) {
					if(nchar(var.sp[2])==2) {col=c("green","red");names(var.sp)=rgnm[match(var.sp,rownames(rgnm)),xi]}
					if(nchar(var.sp[2])==3) {col=c("gray","blue");names(var.sp)=c("Eqal","Both")}
					}
				if(length(var.sp)==3) {col=c("green","red","blue");names(var.sp)=rgnm[match(var.sp,rownames(rgnm)),xi]}
				
				profiles.plot(mcmc[[xi]][var.sp], col.line=col, xlab="", ylab="")		
				legend("top", names(var.sp), col=col, cex=0.8,lty=1,lwd=2,bty="n",x.intersp=0.1,seg.len=1)			
			}
			
			if (i==2){
				rownames(rgnm)=c("xA","xB","sAB")
				if(length(var.ex)==1) {col="gray";names(var.ex)="Equal"}
				if(length(var.ex)==2) {col=c("green","red");names(var.ex)=rgnm[match(var.ex,rownames(rgnm)),xi]}			
				profiles.plot(mcmc[[xi]][var.ex], col.line=col, xlab="", ylab="")
				legend("top", names(var.ex), col=col, cex=0.8,lty=1,lwd=2,bty="n",x.intersp=0.1,seg.len=1)					
			}	
			
			if (i==3){
				if(length(var.n)==1) col="gray"
				if(length(var.n)==2) col=c("green","red")			
				profiles.plot(mcmc[[xi]][var.n], col.line=col, xlab="", ylab="")			
			}	
			if (i==4){
				if(length(var.d)==1) col="gray"
				if(length(var.d)==2) col=c("green","red")			
				profiles.plot(mcmc[[xi]][var.d], col.line=col, xlab="", ylab="")			
			}	
			#画图结束
			
			if (xi == 1) {
				axis(side=2, label=T)
				mtext(side=2, text=index[i], line=1.6, cex=1, las=0)				
			}			
			if (i == length(index)) {
				axis(side=1, lwd = 0.5)
				mtext(side=1, text=region[xi], line=1.6, cex=1, las=0)
			}						
			if (xi == length(region)) {par(mar=c(0,0,0,0)); plot.new()}
			if (xi == 1 & i == 0.5*length(index)) {
                 mtext(side=2, text="Posterior probability density", line=3, cex=1.5)
		     }
		}
	}
	
	
#谱系和环境的关系
	
	#GLM作图
	write.csv(reglm,"reglm.csv")#手动调节变量顺序,改名为rename.csv
	#添加node列	
	rename=read.csv("rename.csv")
	reglm2=cbind(node=rename[match(reglm$region,rename$region),"node"],reglm)
	recoef2=cbind(node=rename[match(recoef$region,rename$region),"node"],recoef)
	rerclim2=cbind(node=rename[match(rerclim$region,rename$region),"node"],rerclim)	
	
	re_exp2=cbind(node=rename[match(re_exp$region,rename$region),"node"],re_exp)
	reispline2=cbind(node=rename[match(reispline$region,rename$region),"node"],reispline)
	regdmcoef2=cbind(node=rename[match(regdmcoef$region,rename$region),"node"],regdmcoef)	
	    
		#将环境变量放到一列
	datglm=c();datrclim=c();datcoef=c()	
	for (i in 1:(dim(reglm2)[2]-2)){
		temp=cbind(reglm2[,1:2],var=colnames(reglm2)[i+2],r2=reglm2[,i+2])	
		datglm=rbind(datglm,temp)
		
		temp2=cbind(rerclim2[,1:3],var=colnames(rerclim2)[i+30],dev=rerclim2[,i+30])
		datrclim=rbind(datrclim,temp2)
		
		temp3=cbind(recoef2[,1:2],var=colnames(recoef2)[i+2],coef=recoef2[,i+2])	
		datcoef=rbind(datcoef,temp3)
	}
	#save.image("gdm-1.RData")
	save.image("gdm-2.RData")
	
	#画图，glm的柱形图
	load("gdm-2.RData")
	library(ggplot2)
	library(ggpubr)
	
	library(scales)
	show_col(hue_pal()(13))
	#9张，以地区为基础的柱形图
	# mycolor=c("#9590FF",#geodis
		# "#EF7F49","#F47B5C","#F8766D","#FB727C","#FE6E8A","#FF6A98","#FF67A4",#Temp
		# "#619CFF","#35A2FF","#00A7FF","#00ACFC",#prep
		# "#00B816")
	mycolor=c("#9590FF","#EF7F49","#FF67A4")	
	dat=datcoef		
	aesx=unique(dat[,c("node","region")])
	plotn=vector("list",9)	
	for (i in 1:9){	 
		subdat=dat[dat$node==i,]
		subaes=aesx[aesx$node==i,]
		subdat$region =factor(subdat$region,levels=as.character(subaes$region))
		 
		p=ggplot(data =subdat ,aes(x = region,y = coef,fill = var)) 		
		if(i==1){
			p=p+geom_bar(stat ="identity",width = 0.6,position = "dodge")+     
				scale_fill_manual(values =mycolor)+ 
				theme(legend.key.width=unit(0.5,'cm'),
				legend.justification = c(0, 1), # pick the upper left corner of the legend box and
				legend.position = c(0.8, 1), # adjust the position of the corner as relative to axis
				legend.background = element_rect(fill = NA), # transparent legend background
				legend.box = "horizontal", # horizontal arrangement of multiple legends
				legend.spacing.x = unit(0.5, units = "cm")# horizontal spacing between legends 	
				)  
		}  		
		if(i>1)	{
			p=p+geom_bar(stat ="identity",width = 0.6,position = "dodge",show.legend=FALSE)+     
				scale_fill_manual(values =mycolor)		
		}
		plotn[[i]]=p+theme(axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=10),
			axis.text.y  = element_text(size=10))+
			expand_limits(y=c(summary(dat[,4])[[1]],summary(dat[,4])[[6]]))
			
			#ggsave(paste("nodescat_",region[i], ".pdf", sep = ""),plot=plotn[[i]])	
	}	
	
	figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]],plotn[[8]],plotn[[9]], 
			nrow=3,ncol = 3,widths=c(2,2,2),heights=c(2,2,2),common.legend=TRUE,legend = "bottom",label.x=0.5,hjust=0,align="hv")
	# figure2=ggarrange(figure,vio,nrow=1,ncol=2,widths=c(8,1),heights=4)		  	
	annotate_figure(figure, bottom = text_grob("plant region", size=14),left = text_grob("coef of glm", rot = 90,size=14))	
	
	#画图，散点图
	mycolor=c("#F564E3","#A3A500","#619CFF")
	dat=datrclim	
	envname=as.character(unique(dat$var))
		
	for (j in 1:length(envname)){
		plotn=vector("list",9)
		for (i in 1:9){	 
			subdat=dat[dat$node==i&dat$var==envname[j],]			
			p=ggplot(data =subdat ,aes(x = dev,y = distance,colour=region))	    
			p=p+geom_point(shape=1,alpha =0.5)			
			
			reg=as.character(unique(subdat$region))
			d1=glm(distance~dev,data = subdat[subdat$region==reg[1],])
			d2=glm(distance~dev,data = subdat[subdat$region==reg[2],])
			d3=glm(distance~dev,data = subdat[subdat$region==reg[3],])
			posx=max(subdat$dev)-min(subdat$dev)			
			p=p+annotate("text",label=c(paste("R2",envname[j]),paste(reg[1],
			 round(ss.glm(d1)[2],2)),paste(reg[2],round(ss.glm(d2)[2],2)),paste(reg[3],round(ss.glm(d3)[2],2))),
			 x=posx*0.5,y=c(1,0.9,0.8,0.7))
			if(ss.glm(d1)[3]<0.05&ss.glm(d1)[2]>10) p=p+geom_smooth(data = subdat[subdat$region==reg[1],], method = "glm",colour="red",size=2)			
			if(ss.glm(d2)[3]<0.05&ss.glm(d2)[2]>10) p=p+geom_smooth(data = subdat[subdat$region==reg[2],], method = "glm",colour="darkgreen",size=2)			
			if(ss.glm(d3)[3]<0.05&ss.glm(d3)[2]>10) p=p+geom_smooth(data = subdat[subdat$region==reg[3],], method = "glm",colour="blue",size=2)
			
			p=p+theme(legend.key.width=unit(0.1,'cm'),
				legend.justification = c(0, 1), # pick the upper left corner of the legend box and
				legend.position = c(0.8, 1), # adjust the position of the corner as relative to axis
				legend.background = element_rect(fill = NA), # transparent legend background
				legend.box = "horizontal", # horizontal arrangement of multiple legends
				legend.spacing.x = unit(0.1, units = "cm")# horizontal spacing between legends 	
				)
			plotn[[i]]=	p+theme(axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x  = element_text(size=10),
				axis.text.y  = element_text(size=10))
			}
		figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]],plotn[[8]],plotn[[9]], 
			nrow=3,ncol = 3,widths=c(2,2,2),heights=c(2,2,2))
	# figure2=ggarrange(figure,vio,nrow=1,ncol=2,widths=c(8,1),heights=4)		  	
	annotate_figure(figure, bottom = text_grob(paste(envname[j],"deviance"), size=14),left = text_grob("phylobeta", rot = 90,size=14))	
	ggsave(paste(envname[j],"_dev", ".jpg", sep = ""),plot=figure)	
	}
	
#画图，GDM
	load("gdm-2.RData")
	library(ggplot2)
	library(ggpubr)	
  ####re_exp2	
	dat=re_exp2		
	aesx=unique(dat[,c("node","region")])
	plotn=vector("list",9)	
	for (i in 1:9){	 
		subdat=dat[dat$node==i,]
		subaes=aesx[aesx$node==i,]
		subdat$region =factor(subdat$region,levels=as.character(subaes$region))
		 
		p=ggplot(data =subdat ,aes(x = region,y = explained,fill="red")) 		
		if(i==1){
			p=p+geom_bar(stat ="identity",width = 0.6,position = "dodge")+    			
				theme(legend.key.width=unit(0.5,'cm'),
				legend.justification = c(0, 1), # pick the upper left corner of the legend box and
				legend.position = c(0.8, 1), # adjust the position of the corner as relative to axis
				legend.background = element_rect(fill = NA), # transparent legend background
				legend.box = "horizontal", # horizontal arrangement of multiple legends
				legend.spacing.x = unit(0.5, units = "cm")# horizontal spacing between legends 	
				)  
		}  		
		if(i>1)	{
			p=p+geom_bar(stat ="identity",width = 0.6,position = "dodge",show.legend=FALSE)     
						
		}
		plotn[[i]]=p+theme(axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=10),
			axis.text.y  = element_text(size=10))+
			expand_limits(y=c(summary(dat[,5])[[1]],summary(dat[,5])[[6]]))
			#ggsave(paste("nodescat_",region[i], ".pdf", sep = ""),plot=plotn[[i]])	
	}	
	
	figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]],plotn[[8]],plotn[[9]], 
			nrow=3,ncol = 3,widths=c(2,2,2),heights=c(2,2,2),common.legend=TRUE,legend = "bottom",label.x=0.5,hjust=0,align="hv")
	# figure2=ggarrange(figure,vio,nrow=1,ncol=2,widths=c(8,1),heights=4)		  	
	annotate_figure(figure, bottom = text_grob("plant region", size=14),left = text_grob("explain power of gdm", rot = 90,size=14))	
	
####regdmcoef2
	mycolor=c("#9590FF",#intercept
		"#00B816",#geodis
		"#EF7F49","#F47B5C","#F8766D","#FB727C","#FE6E8A","#FF6A98","#FF67A4",#Temp
		"#619CFF","#35A2FF","#00A7FF","#00ACFC")#prep
			
	dat=regdmcoef2		
	aesx=unique(dat[,c("node","region")])
	varname=unique(dat[,"var"])
	plotn=vector("list",9)	
	for (i in 1:9){	 
		subdat=dat[dat$node==i,]
		subaes=aesx[aesx$node==i,]		
		subdat$region =factor(subdat$region,levels=as.character(subaes$region))
		subdat$var=factor(subdat$var,levels=as.character(varname))
		
		p=ggplot(data =subdat ,aes(x = region,y = gdmcoef,fill = var)) 		
		if(i==1){
			p=p+geom_bar(stat ="identity",width = 0.6,position = "dodge")+     
				scale_fill_manual(values =mycolor)+ 
				theme(legend.key.width=unit(0.5,'cm'),
				legend.justification = c(0, 1), # pick the upper left corner of the legend box and
				legend.position = c(0.8, 1), # adjust the position of the corner as relative to axis
				legend.background = element_rect(fill = NA), # transparent legend background
				legend.box = "horizontal", # horizontal arrangement of multiple legends
				legend.spacing.x = unit(0.5, units = "cm")# horizontal spacing between legends 	
				)  
		}  		
		if(i>1)	{
			p=p+geom_bar(stat ="identity",width = 0.6,position = "dodge",show.legend=FALSE)+     
				scale_fill_manual(values =mycolor)		
		}
		plotn[[i]]=p+theme(axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x  = element_text(size=10),
			axis.text.y  = element_text(size=10))+
			expand_limits(y=c(summary(dat[,4])[[1]],summary(dat[,4])[[6]]))
			#ggsave(paste("nodescat_",region[i], ".pdf", sep = ""),plot=plotn[[i]])	
	}	
	
	figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]],plotn[[8]],plotn[[9]], 
			nrow=3,ncol = 3,widths=c(2,2,2),heights=c(2,2,2),common.legend=TRUE,legend = "bottom",label.x=0.5,hjust=0,align="hv")
	# figure2=ggarrange(figure,vio,nrow=1,ncol=2,widths=c(8,1),heights=4)		  	
	annotate_figure(figure, bottom = text_grob("plant region", size=14),left = text_grob("coef of gdm", rot = 90,size=14))	
	
####reispline2
	plot(reispline2[,"bio19"], exSplines[[2]][,"bio19"], type="l",
lwd=3, xlab="Winter precipitation (mm)", ylab="Partial Ecological Distance")
	
	mycolor=c("#F564E3","#A3A500","#619CFF")
	dat=reispline2	
	envname=as.character(unique(dat$varname))
		
	for (j in 1:length(envname)){
		plotn=vector("list",9)
		for (i in 1:9){	 
			subdat=dat[dat$node==i&dat$varname==envname[j],]			
			p=ggplot(data =subdat ,aes(x = obs,y = pre,group=region,colour=region))	+geom_line(size=2)			
			p=p+theme(legend.key.width=unit(0.1,'cm'),
				legend.justification = c(0, 1), # pick the upper left corner of the legend box and
				legend.position = c(0, 1), # adjust the position of the corner as relative to axis
				legend.background = element_rect(fill = NA), # transparent legend background
				legend.box = "horizontal", # horizontal arrangement of multiple legends
				legend.spacing.x = unit(0.1, units = "cm")# horizontal spacing between legends 	
				)+expand_limits(y=c(0,1))
			plotn[[i]]=	p+theme(axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				axis.text.x  = element_text(size=10),
				axis.text.y  = element_text(size=10))
			}
		figure=ggarrange(plotn[[1]],plotn[[2]],plotn[[3]],plotn[[4]],
			plotn[[5]],plotn[[6]],plotn[[7]],plotn[[8]],plotn[[9]], 
			nrow=3,ncol = 3,widths=c(2,2,2),heights=c(2,2,2))
	# figure2=ggarrange(figure,vio,nrow=1,ncol=2,widths=c(8,1),heights=4)		  	
	annotate_figure(figure, bottom = text_grob(paste(envname[j]), size=14),left = text_grob("Partial Ecological Distance", rot = 90,size=14))	
	ggsave(paste(envname[j],"_ispline", ".jpg", sep = ""),plot=figure)	
	}
	

## fuzzy c-means https://www.modb.pro/db/137902
	library(cluster)		
	#options("install.lock"=FALSE)
	library(factoextra)
	library(raster)
	library(ggpubr)
	data.test <- raster("extend.tree/rasters/0.00_temperature.asc",native=TRUE)
	coord=coordinates(data.test)
	rownames(coord)=paste(coord[,1],coord[,2],sep="_")
	cordacd <- read.csv("POINT2.csv")
	cordacd=subset(cordacd,cordacd$ADCODE99>0)
	rownames(cordacd)=paste(cordacd[,1],cordacd[,2],sep="_")	
	
	#分类效果评估	
	psim1=get(load("processed_aV2/psim.raxml_v2_dated_Av1_final_clean_ann.full.polyRes.Rdata"))#matrix
	psim1[is.na(psim1)]=0
	psim.t=psim1+t(psim1)
	psim2=as.dist(psim.t)
	k=c(10,15,20,25,30,36,40,45,50,60,70,80)	
	por=c();plt=list()
	for (m in 1:length(k)){
			cm.fanny <- fanny(psim2, k = 8, diss = TRUE, maxit = 8000,memb.exp = 1.1)	
			#作图
			#windows()
			k2=length(unique(cm.fanny$clustering))#最佳分类，即各对象最接近的聚类簇	
			#轮廓图	
			p=fviz_silhouette(cm.fanny,ggtheme = theme_minimal(),palette=rainbow(k2),legend = "none")	
			fcm.value=cm.fanny$clustering
			fcm.value[!names(fcm.value)%in%rownames(subset(cm.fanny$silinfo$widths,cm.fanny$silinfo$widths[,"sil_width"]>0))]=0
			cordacd2=cbind(cordacd,fcm=fcm.value[match(cordacd$ADCODE99,names(fcm.value))])	
			pC.t4=cbind(coord,cordacd2[match(rownames(coord),rownames(cordacd2)),])		
			data.test[]=pC.t4[,"fcm"]			
			test_spdf <- as(data.test, "SpatialPixelsDataFrame")
			test_df <- as.data.frame(test_spdf)
			colnames(test_df) <- c("value", "x", "y")	
			p.map=ggplot()+geom_raster(data = test_df , aes(x = x, y = y,fill = factor(value))) + 		
				scale_fill_manual(name = "Kingdoms",values = c("gray",rainbow(k2)))+		
				ggtitle(paste("Initial cluster:",8))+ 	theme_bw()+	
				theme(axis.title = element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
				panel.grid=element_blank(),plot.title = element_text(hjust = 0.5))			
			plt[[m]]=ggarrange(p.map,p,nrow=2,ncol = 1,widths=2,heights=c(1,1),labels=c("a","b"),label.x=0.02)
	
			#求算组间beta与总beta的比值
			resim=cm.fanny$clustering		
			classnum=unique(cm.fanny$clustering)
			submat=matrix(data=0,nrow=(length(classnum)-1),ncol=(length(classnum)-1))
			colnames(submat)=classnum[2:length(classnum)]
			rownames(submat)=classnum[1:(length(classnum)-1)]			
			for (i in 1:(length(classnum)-1)){
				for (j in (i+1):length(classnum)){
					a=subset(resim,resim==i)
					b=subset(resim,resim==j)
					matemp=subset(psim.t,rownames(psim.t)%in%names(a))
					matemp=matemp[,colnames(psim.t)%in%names(b)]
					submat[i,j-1]=sum(matemp)
				}
			}		
			por.t=data.frame(Desired.clusters=max(classnum),Final.clusters=k,Pbeta=sum(submat)/(sum(psim.t)/2)*100)
			por=rbind(por,por.t)
		}
		write.csv(por,"fuzz.por.csv")
		por=read.csv("fuzz.por.csv")
		p1=ggplot(por, mapping = aes(x = Desired.clusters, y = Pbeta)) + geom_line()
		p2=ggplot(por, mapping = aes(x = Final.clusters, y = Pbeta)) + geom_line()
		ggarrange(p1,p2,nrow=2,ncol = 1,widths=2,heights=c(1,1),labels=c("a","b"),label.x=0.02)
	}

