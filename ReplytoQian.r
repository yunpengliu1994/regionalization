#2023 powo
library(dplyr);library(rWCVP);library(data.table)
genames <- rWCVPdata::wcvp_names%>%filter(taxon_rank=="Genus")
gyn=fread("vascular_genera_ZW.csv")%>%filter(Division=="A")
distributions <- rWCVPdata::wcvp_distributions%>%filter((plant_name_id%in%genames$plant_name_id))
disPOWO=distributions%>%left_join(genames,by="plant_name_id")%>%select(introduced,genus,area_code_l3)%>%distinct()#%>%filter(genus%in%gyn$Genus_E)
intr=disPOWO%>%filter(introduced==1)
dim(disPOWO);dim(intr)
length(unique(disPOWO$genus));length(unique(intr$genus))

#2019 powo
POWO2019=fread("PTWDistribution.txt",header=T)
mf=function(x) strsplit(x," ")[[1]][[1]]
POWO2019$genus=do.call(rbind,lapply(POWO2019$sp,mf))
mf2=function(x) ifelse(x=="native",0,1)
POWO2019$status2=do.call(rbind,lapply(POWO2019$status,mf2))
POWO2019.genus=POWO2019%>%group_by(genus,geoname)%>%summarise(min_in=min(status2))
intr=POWO2019.genus%>%filter(min_in==1)
dim(POWO2019.genus);dim(intr)

#install.packages("rWCVP")
library(data.table);library(dplyr)
library(foreign);
library(rWCVP);
genusdis=fread("APPENDIX S2 GenusDistribution.csv")[,1:2]%>%distinct()
colnames(genusdis)[2]="genus"
genames <- rWCVPdata::wcvp_names%>%filter(taxon_rank=="Genus")#&genus%in%genusdis$Accepted_name)
geo=read.dbf("phyloRealm_WGSRPD/phylo_WGSRPD.dbf")
#since most of GSU are smaller or equal to wcvp geounit (WGSRPD level 3), thus use GSU of our study to caculate proportion
geo$pro=geo$Area_int/geo$Area*100
geo=geo%>%filter(pro>=25)
# #xiaoting's method
# geo=read.csv("GSU_lev3.csv")[,-1]
# colnames(geo)[c(1,5)]=c("LEVEL3_COD","GeoID")
# geoid=unique(read.dbf("phyloRealm_WGSRPD/phylo_WGSRPD.dbf")$GeoID)
# geo=geo%>%filter(GeoID%in%geoid)
disPOWO=rWCVPdata::wcvp_distributions%>%left_join(genames,by="plant_name_id")%>%left_join(geo,by=c("area_code_l3"="LEVEL3_COD"),multiple = "all",relationship ="many-to-many")%>%
	select(GeoID,genus,introduced)%>%na.omit()%>%distinct()%>%group_by(GeoID,genus)%>%summarize(introduced=min(introduced))#in 330records, a geoID matches mutiple units, thus use group_by
#save(disPOWO,file="disPOWO.RData")
#disPOWO=get(load("disPOWO.RData"))
genusdis2=genusdis%>%left_join(disPOWO,by=c("GeoID","genus"))%>%na.omit()
#save(genusdis2,file="disGSU.RData")

#genusdis2=get(load("disGSU.RData"))
intro.stat.geo=genusdis2%>%group_by(GeoID)%>%summarize(intro=mean(introduced))
dim(intro.stat.geo%>%filter(intro==0))#0
dim(intro.stat.geo%>%filter(intro<=0.15))#416

intro.stat.genus=genusdis2%>%group_by(genus)%>%summarize(intro=mean(introduced))
dim(intro.stat.genus%>%filter(intro==0))#10247 out of 11578 (88.5%) genus do not contain intro;#when geo$pro >= 50, 10135
length(unique(genusdis2$genus))
intro=genusdis2%>%filter(introduced==1);dim(intro)#16566,when geo$pro >= 50, intro = 16024

#Qian's method
geo=read.dbf("phyloRealm_WGSRPD/phylo_WGSRPD.dbf")
distributions <- rWCVPdata::wcvp_distributions%>%filter((plant_name_id%in%genames$plant_name_id)&(area_code_l3%in%geo$LEVEL3_COD))
disPOWO=distributions%>%left_join(genames,by="plant_name_id")%>%left_join(geo,by=c("area_code_l3"="LEVEL3_COD"),multiple = "all",relationship ="many-to-many")%>%
	select(ID_int,GeoID,Area,area_code_l3,Area_lv3,plant_name_id,genus,introduced,Area_int)%>%
	distinct()
genusdis2=genusdis%>%left_join(disPOWO,by=c("GeoID","genus"))%>%na.omit()
intro.stat=genusdis2%>%group_by(genus,area_code_l3)%>%summarize(intro=max(introduced))%>%left_join(geo,by=c("area_code_l3"="LEVEL3_COD"),multiple="all",relationship ="many-to-many")
intro.qian=intro.stat%>%filter(intro==1)
length(unique(intro.qian$genus));
dim(intro.qian %>% select(genus,GeoID) %>% distinct())

# library(rgdal);library(ggpubr)
# map <- readOGR(dsn = "phyloRealms/phyloRealms.dbf",stringsAsFactors=FALSE);map1<-fortify(map)
# cordacd<-as.data.frame(map@data)%>%left_join(intro.stat.geo,by="GeoID");
# cordacd[is.na(cordacd)]=0
# cordacd$id<-as.character(0:(length(unique(map1$id))-1))
# map2<-inner_join(map1,cordacd)
# p1=ggplot()+geom_polygon(data=map2,aes(x=long,y=lat,group=group,fill=intro),color="transparent")+ 
		# scale_fill_gradient2(midpoint=0.5,limits=c(0,1),low="blue",mid="white",high="red",name="Proportion of\nnon-native records")+
		# theme(#aspect.ratio=(diff(range(test_df$y))+5)/(diff(range(test_df$x))+10*max(nich.rate)+2),#+75*max(srlist$P)
			  # #text = element_text(size =10),
			  # legend.background=element_rect(fill='transparent'),
			  # panel.border = element_blank(),
			  # axis.title = element_blank(),
			  # axis.text = element_text(face="bold",size =10),
			  # axis.text.x = element_text(angle=0),
			  # axis.ticks = element_blank(),
			  # panel.grid.minor = element_blank(),
			  # legend.text=element_text(face="bold",size=10),
			  # legend.title=element_text(face="bold",size=12),
			  # legend.position = c(0.12,0.38)) 	
# p2=ggplot(intro.stat.geo, aes(x=intro)) + geom_histogram(bins=30)+xlim(-0.1,1)+
	# labs(x="Proportion of non-native species in each GSU",y="Number of GSU")+theme_bw()
# p3=	ggplot(intro.stat.genus, aes(x=intro)) + geom_histogram(bins=30)+xlim(-0.1,1)+
	# labs(x="Proportion of non-native distributions in each genus",y="Number of genera")+theme_bw()
# part1=ggarrange(p2,p3,nrow=1,ncol =2,labels=c("b","c"),font.label = list(size = 20))
# ggarrange(p1,part1,nrow=2,ncol =1,heights=c(1.8,1),labels=c("a",""),font.label = list(size = 20))
genusdis2=genusdis%>%left_join(disPOWO,by=c("GeoID","genus"))
genusdis2[is.na(genusdis2)]=0

geo.stat=geo%>%group_by(LEVEL3_COD)%>%summarise(n=length(GeoID))
geo2=geo%>%left_join(geo.stat,by="LEVEL3_COD")%>%filter(n==1)
powo.add=rWCVPdata::wcvp_distributions%>%left_join(genames,by="plant_name_id")%>%left_join(geo2,by=c("area_code_l3"="LEVEL3_COD"),multiple = "all",relationship ="many-to-many")%>%
	select(GeoID,genus,introduced)%>%na.omit()%>%distinct()%>%group_by(GeoID,genus)%>%summarize(introduced=max(introduced))%>%filter(introduced==0&genus%in%genusdis$genus)
genusdis3=genusdis2%>%filter(introduced==0)%>%rbind(powo.add)%>%distinct()
dim(genusdis3);length(unique(genusdis3$GeoID));length(unique(genusdis3$genus))

native=disPOWO%>%filter(introduced==0&genus%in%c("Bontia","Cocos","Rivina","Synedrella","Tamarindus","Zea"))
genusdis.coret=genusdis3%>%filter(!genus%in%c("Bontia","Cocos","Rivina","Synedrella","Tamarindus","Zea"))%>%
	rbind(native)
spname=fread("APPENDIX S2 GenusDistribution.csv")[,-1]%>%distinct()
genusdis4=cbind(genusdis.coret[,-3],spname[match(genusdis.coret$genus,spname$Accepted_name),])

a=genusdis4%>%filter(genus=="Cocos")
geo.a=geo%>%filter(GeoID%in%a$GeoID)
a=genusdis4%>%filter(genus=="Rivina")
geo.a=geo%>%filter(GeoID%in%a$GeoID)
a=genusdis4%>%filter(GeoID==227);length(unique(a$genus))

library(ape)
tre=read.tree("raxml_v2_dated_Av1_final_clean_ann.full.polyRes.tre")
genusdis4=genusdis4 %>%	filter(genus%in%tre$tip.label)
save(genusdis4,file="disGSU_POWOcorrected.RData")
tre2=keep.tip(tre,tip=unique(genusdis4$genus))
write.tree(tre2,"phylo.powo.tre")

#re-run regionalization
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
	
spdis=genusdis4[,1:2]%>%distinct()
spdis=get(load("disGSU_POWOcorrected.RData"))[,1:2]%>%distinct()#add genus that not in POWO
colnames(spdis)[1:2]=c("Adcode99","Species_E1")
d=spdat(spdis)
save(d,file="Spdat_POWOcorrected.Rdata")

#caculate phylo beta
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

for (i in 1:ncol(d)) d[,i] <- as.logical(d[,i])
d <- as.matrix(d)

## community pairs
comm <- list()
ii <- 1
for (i in 1:(nrow(d)-1)) {
	for (j in (i+1):nrow(d)) {
		comm1 <- colnames(d)[which(as.logical(d[i,]))]
		comm2 <- colnames(d)[which(as.logical(d[j,]))]
		comm[[ii]] <- list(comm1, comm2)
		ii <- ii + 1
		}
	print(i)
	}

tree <- read.tree("phylo.powo.tre")
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
	save(pBeta.mat, file=paste("psim.POWOcorrected.Rdata",sep=""))
	stopCluster(mycl)
	
	#Plot dend,nmds and map	
	library(sp)
	library(maptools)
	library(dendextend)
	xy <- read.csv("Geo-isrm.csv",head=T)[,c("ADCODE99","Lon","Lat")]
	xy=cbind(no=1:420,xy)	
	mycolor=c("#93AA00","#D39200","#619CFF","#F8766D","#00C19F","#00BA38","#DB72FB","#FF61C3")	
		psim1=pBeta.mat#matrix	
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
		#color[color$ADCODE99==107,"leaf_col"]="#93AA00"
		shape<-readShapeSpatial("PhyloRealms/PhyloRealms.shp")			
		shape@data=cbind(shape@data,phcol=color[match(shape@data$GeoID,color$ADCODE99),"leaf_col"])
		plot(shape, col="gray", border = "white")				
		for (i in 1:length(mycolor)) {
			shp=subset(shape,shape@data$phcol==mycolor[i])		
			plot(shp,col=mycolor[i],border=mycolor[i],add=T)
		}
		box()
		
		scales::show_col(c("#93AA00","#D39200","#619CFF","#F8766D","#00BA38","#00C19F","#FF61C3","#DB72FB"))
	realm=data.frame(PhyloRealm=c("Saharo-Arabian","Holarctic","Chile-Patagonian","Neotropical","Indo-Malesian","African","Neovozealandic","Australian"),
		mycolor.kind.powo=c("#93AA00","#D39200","#619CFF","#F8766D","#00BA38","#00C19F","#FF61C3","#DB72FB"))
	shape2=shape
	shape2@data=shape2@data %>% left_join(realm,by=c("phcol"="mycolor.kind.powo")) %>% select(GeoID,NAME99,Lon,Lat,PhyloRealm)
	colnames(shape2@data)[1:2]=c("GeoID","NAME")
	library(rgeos);library(rgdal)
	writeOGR(shape2, dsn = 'PhyloRealms.new', layer = 'PhyloRealms.new', driver = "ESRI Shapefile")
	