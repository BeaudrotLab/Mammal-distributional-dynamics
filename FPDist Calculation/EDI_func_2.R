#EDI=Ecological Distance Index
EDI=function(DF){ 	
SB=as.data.frame(DF)
SB.df=SB[,-1]
SB.cameras=list()

to=dim(SB.df)[2]


for (i in 1:to){ 
com=SB$SpID[SB.df[,i]>0]
	
mat=matrix(nrow=length(com),ncol=length(com))

cnt=0
cntt=0
for (k in unique(com))
{
	cnt=cnt+1
	for(q in unique(com)){
    cntt=cntt+1
mat[cnt,cntt]=FPDist.global[k,q]	
}
dimnames(mat)=list(unique(com),unique(com))
cntt=0
}
SB.cameras[[i]]=mat
}

####EDI_obs
SB.edi=list()
for (h in 1:length(SB.cameras)){
	tm=SB.cameras[[h]]
	div=dim(tm)[1]-1
	tm.edi=rowSums(tm)/div
	tm.edi[is.na(tm.edi)]<-0	
	SB.edi[[h]]=tm.edi		
}

return(SB.edi)
	
}