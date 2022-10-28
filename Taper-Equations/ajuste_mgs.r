rm(list=ls(all.names = T));

cub<-openxlsx::read.xlsx('../dados/cubagem_ufscar_original.xlsx', sheet=1)

#definição da estratificacao
def_estrato<-'ndef'
#def_estrato<-'clidade'
#def_estrato<-'matgen'
#def_estrato<-'matgen + cldap'

cub$id<-cub$idarv
cub$espcasca<-(cub$dicc-cub$disc)/2

if(def_estrato=='matgen + cldap'){
  cub<-cub[order(cub$matgen,cub$dap),]
  f<-function(x){
    bks<-c(-Inf,quantile(x$dap,probs = c(0.25,0.5,0.75)),Inf)
    x$cldap<-as.character(cut(x$dap,breaks = bks,include.lowest = F))
    return(x)
  }
  cub<-Reduce('rbind',lapply(split(cub,cub$matgen),f))
  cub$idestrato<-cmrinvflor::chaveordem(cub[,c('matgen','cldap')])
  desc_estrato<-cub[!duplicated(cub$idestcoef),c('idestcoef','matgen','cldap')]
  desc_estrato$li<- as.numeric(sub("\\[(.+),.*\\]", "\\1",sub("\\((.+),.*", "\\1", desc_estrato$cldap)))
  desc_estrato$ls<- as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", desc_estrato$cldap))
} else if(def_estrato=='ndef'){
  cub$idestrato<-1
}

f<-function(x){
  x$dap<-x$dicc[x$hi==1.3]
  x$vref<-x$x_vicc
  x$vref[1]<-x$hi[1]*pi/40000*x$disc[1]^2
  x<-x[,c('idestrato','id','dap','ht','hi','dicc','disc','espcasca','vref')]
  xp<-x[nrow(x),]
  xp$vref<-pi/120000*xp$disc^2*(xp$ht-xp$hi)
  xp$hi<-xp$ht
  xp$dicc<-0;
  xp$espcasca<-0;
  x<-rbind(x,xp)
  x$disc<-NULL
  return(x)
}

cub<-Reduce('rbind',lapply(split(cub,cub$id),f))
cub<-subset(cub,!id %in% c(18,19))

# openxlsx::write.xlsx(
#   x=cub,
#   file = '../dados/cubagem_ufscar.xlsx',
#   rowNames=F,
#   sheetName="validada"
# )
  
#source('taperMGS.r');

#Método de estimativa do volume de referência
calc_vref<-NULL

curva_media<-T; 
if(curva_media){
  fstat<-function(...) mean(...,trim=0.05) #Função para a medida de posição (media, mediana, etc) do di em cada hi/ht por estrato
}  

com_casca<-F
#pr: predefinição das posições relativas em caso de cubagem em posições absolutas e a opção curva_media==T
pr<-c(seq(0,0.9,0.05),seq(0.91,1,0.01)) 
#pr<-NULL #Posições relativas sem atribuição caso a cubagem já tenha sido realizada em posições relativas 

#curva_media<-F
eafil<-taperMGS(cub, vol_ref='smalian', pr=pr, curva_media=curva_media, fstat=fstat)

saveRDS(list(desc_estrato=desc_estrato, curva_media=curva_media, eafil=eafil), 
                 file.path('..','dados','taper_smalian.rds'))

