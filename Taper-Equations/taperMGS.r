# Criado em:     2021-10-29                                                          
# Modificado em: 2022-10-25                                                          

#' @title Taper function using multiple geometric solids
#' @description Obtenção das função de afilamento considerando a superposição de múltiplos sólidos geométricos.
#' @param  cub 'data frame' \cr
#' ..$ idestrato: Chave primária do estrato\cr
#' ..$ id: Chave primária do fuste\cr
#' ..$ dap: Diâmetro à altura do peito (cm) \cr
#' ..$ ht:  Altura total(m) \cr
#' ..$ hi: Iésima altura de medição do diâmetro ao longo do fuste (m) \cr
#' ..$ dicc: Diâmetro com casca na iésima altura de medição ao longo do fuste (m) \cr 
#' ..$ espcasca: Espessura da casca (cm) 
#' @param vol_ref Método de estimativa do volume de referência para reconstrução do sólido geométrico
#' @param com_casca Estimar volumes com casca (TRUE) ou sem casca (FALSE)
#' @param curva_media Obter as funções para os perfis médios por estrato
#' @param pr Posições relativas fixas para todos os futes.  Caso curva_media==T
#' @param fstat Função para estimativa do perfil médio dos fustes por estrato. Caso curva_media==T
#' @return 'list'\cr 
#' $ fst 'data.frame'\cr
#' ..$ idestrato: Chave primária do estrato\cr
#' ..$ id: Chave primária do fuste\cr
#' ..$ dap: Diâmetro à altura do peito (cm) \cr
#' ..$ ht:  Altura total(m) \cr
#' ..$ eafil:  Expressão da função de afilamento \cr
#' $ cm 'data.frame' (caso curva_media==T)\cr
#' ..$ idestrato: Chave primária do estrato\cr
#' ..$ eafil:  Expressão da função de afilamento \cr
#' @author Cláudio Roberto Thiersch \email{crthiersch@@ufscar.br},\cr 
#' Monica F. B. M. Thiersch \email{monicathiersch@@ufscar.br}
#' @rdname taperMGS
#' @export
taperMGS<-function(cub, calc_vref='smalian', com_casca=T, curva_media=T, pr=NULL, fstat=mean){
  if(!is.null(calc_vref)){
    names(cub)<-c('idestrato','id','dap','ht','hi','dicc','espcasca')
    cub<-cub[order(cub$id,cub$hi),]
    
    #Inserindo hi=ht, hi=0 e hi=1.3 em caso de inexistência
    f<-function(x){
      x<-x[,c('idestrato','id','dap','ht','hi','dicc','espcasca')]
      if(match(x$ht[1],x$hi,nomatch=0)==0){
        x<-rbind(x, data.frame(idestrato=x$idestrato[1], id=x$id[1], dap=x$dap[1],ht=x$ht[1],hi=x$ht[1],dicc=0,espcasca=0))
      }
      if(match(0,x$hi,nomatch=0)==0){
        x<-rbind(x, data.frame(idestrato=x$idestrato[1], id=x$id[1], dap=x$dap[1],ht=x$ht[1],hi=0,dicc=x$dicc[1],espcasca=x$espcasca[1]))
      }
      if(match(1.3,x$hi,nomatch = 0)==0){
        mdist<-abs(1.3-x$hi)
        ii<-match(min(mdist),mdist)
        x<-rbind(x, data.frame(idestrato=x$idestrato[1], id=x$id[1], dap=x$dap[1],ht=x$ht[1],hi=1.3,dicc=x$dap[1],espcasca=x$espcasca[ii]))
      }
      return(x)
    }
    
    if(com_casca){cub$di<-cub$dicc}else{cub$di<-cub$dicc-2*cub$espcasca}
    ii<-cub$di<0; if(sum(ii)>0) cub<-cub[!ii,]
    
    if(calc_vref == 'smalian'){

      pair<-function(x){
        x<-x[order(x$hi),];
        ii<-matrix(sort(c(1:(nrow(x)-1),2:nrow(x))),nrow(x)-1,2,byrow = T)
        return(data.frame(idestrato=x$idestrato[ii[,1]],
                          id=x$id[ii[,1]],
                          dap=x$dap[ii[,1]],ht =x$ht[ii[,1]],
                          hi1=x$hi[ii[,1]],hi2=x$hi[ii[,2]],
                          di1=x$di[ii[,1]],di2=x$di[ii[,2]],
                          stringsAsFactors = F))  
      }
      fsmalian<-function(x) {
        if(x$hi1!=0){
          vsmal<-pi/80000*(x$di2^2+x$di1^2)*(x$hi2-x$hi1)
        }else{
          vsmal<-pi/120000*(x$di2^2)*(x$hi2-x$hi1)
        }
        return(vsmal)
      }
      
      tora<-Reduce('rbind',lapply(split(cub,cub$id),pair))
      tora$vref<-Reduce('c',Map(fsmalian,split(tora,1:nrow(tora))))
    }
  } else{
    #Inserindo hi=0 em caso de inexistência
    f<-function(x){
      x<-x[,c('idestrato','id','dap','ht','hi','dicc','espcasca','vref')]
      if(match(0,x$hi,nomatch=0)==0){
        x<-rbind(
          data.frame(
            idestrato=x$idestrato[1], 
            id=x$id[1], 
            dap=x$dap[1],
            ht=x$ht[1],
            hi=0,
            dicc=x$dicc[1],
            espcasca=x$espcasca[1],
            vref=0
          ),
          x 
        )
      }
      return(x)
    }
    cub<-Reduce('rbind',lapply(split(cub,cub$id),f))
    
    if(com_casca){cub$di<-cub$dicc}else{cub$di<-cub$dicc-2*cub$espcasca}
    ii<-cub$di<0; if(sum(ii)>0) cub<-cub[!ii,]

    cub$hi<-cub$ht-cub$hi
    cub$phi<-cub$hi/cub$ht
    if(is.null(pr)) pr<-sort(unique(cub$phi))

    pair<-function(x){
      x<-x[order(x$hi),];
      ii<-matrix(sort(c(1:(nrow(x)-1),2:nrow(x))),nrow(x)-1,2,byrow = T)
      return(data.frame(idestrato=x$idestrato[ii[,1]],
                        id=x$id[ii[,1]],
                        dap=x$dap[ii[,1]],ht =x$ht[ii[,1]],
                        hi1=x$hi[ii[,1]],hi2=x$hi[ii[,2]],
                        di1=x$di[ii[,1]],di2=x$di[ii[,2]],
                        vref=x$vref[ii[,2]],
                        stringsAsFactors = F))  
    }
    tora<-Reduce('rbind',lapply(split(cub,cub$id),pair))
  }  

  vigs<-function(l,di1,k,r){
    integrate(function(l,di1,k,r) pi*((di1+k*l^r)^2)/40000, lower=0, upper=l, di1, k, r,stop.on.error = FALSE)$value}
  
  fr<-function(s){
    if(!is.list(s)) s<-as.list(s)
    if(s$hi1==0){
      r=1
    }else if(s$di2<=s$di1){
      r=0
    }else{
      fo<-function(x, y){
        l<-y$hi2-y$hi1
        k<-(y$di2-y$di1)/(l^x)
        return(abs(vigs(l,di1=y$di1,k,r=x)-s$vref))
      }
      r<-optimize(f=fo,interval=c(0,10),y=s,maximum=F)$minimum;
    }
    return(r)
  }
  
  create_taper_mgs<-function(x){
    f<-function(x){ 
      ii<-2:nrow(x)
      eafil<-c(paste0('(',x$pdi1[1],'*dap+',x$k[1],'*(abs(hi-',x$phi1[1],'*ht)^',x$r[1],'))*(round(hi/ht,4) >= ',x$phi1[1],' & round(hi/ht,4) <= ',x$phi2[1],')'),
               paste0('(',x$pdi1[ii],'*dap+',x$k[ii],'*(abs(hi-',x$phi1[ii],'*ht)^',x$r[ii],'))*(round(hi/ht,4) > ',x$phi1[ii],' & round(hi/ht,4) <= ',x$phi2[ii],')'))           
    }
    x$eafil<-Reduce('c',by(x,x$id,f))
    
    eafil<-data.frame(idestrato=x$idestrato[1], id=x$id[1], dap=x$dap[1], ht=x$ht[1], eafil=paste0(x$eafil,collapse = '+'), stringsAsFactors = F)
    eafil$eafil<-str2expression(eafil$eafil);
    return(eafil)
  }
  
  tora$pdi1<-round(tora$di1/tora$dap,4)
  tora$phi1<-round(tora$hi1/tora$ht,4)
  tora$phi2<-round(tora$hi2/tora$ht,4)

  tora$r<-apply(tora,1,fr)
  tora$k<-with(tora, (di2-di1)/(hi2-hi1)^r)
  
  ##Verificando a acurária na estimativa dos parâmetros k e r
  # new_tora<-tora[,c('hi1','hi2','ht','di1','k','r')]
  # new_tora$hi2<-new_tora$hi2-new_tora$hi1
  # new_tora$hi1<-0
  # tora$vest<-as.vector(apply(new_tora,1,fvi,ps=NULL,eafil=expression(di1+k*hi^r)))
  # tora$erro<-with(tora,round(vref-vest,5))
  # summary(tora$erro)
  
  fst_eafil<-Reduce('rbind', by(tora, tora$id, create_taper_mgs))

  if(curva_media){
    f<-function(x,pr){
      eafil<-x$eafil
      x$eafil<-NULL
      ncub<-merge(x,data.frame(idestrato=x$idestrato, id=x$id, phi=pr,stringsAsFactors = F))
      ncub$hi<-with(ncub,ht*phi)
      ncub$di<-with(ncub,eval(eafil))
      return(ncub)
    }
    cub<-Reduce('rbind',lapply(split(fst_eafil,1:nrow(fst_eafil)),f,pr))
    
    mcub<-with(subset(cub,!duplicated(id),c('idestrato','dap','ht')),
               aggregate(list(dap=dap,ht=ht),list(idestrato=idestrato, id=idestrato),mean))
    f<-function(x,fstat){
      return(data.frame(
        idestrato=x$idestrato[1],
        phi=x$phi[1],
        di=fstat(x$di),
        stringsAsFactors = F))
    }
    cub<-merge(mcub,Reduce('rbind',lapply(split(cub,cub[,c('idestrato','phi')]),f,fstat)))
    cub$hi<-cub$phi*cub$ht
    cub<-cub[order(cub$dap,cub$phi),]
    row.names(cub)<-NULL;
  }
  
  if(is.null(calc_vref)) calc_vref<-"pre_calc"
  if(calc_vref=='smalian'){
    tora<-Reduce('rbind',lapply(split(cub,cub$id),pair))
    tora$pdi1<-round(tora$di1/tora$dap,4)
    tora$phi1<-round(tora$hi1/tora$ht,4)
    tora$phi2<-round(tora$hi2/tora$ht,4)
    tora$vref<-Reduce('c',Map(fsmalian,split(tora,1:nrow(tora))))
    
    tora$r<-Reduce('c',lapply(split(tora,1:nrow(tora)),fr))
    tora$k<-with(tora, (di2-di1)/(hi2-hi1)^r)
    cm_eafil<-Reduce('rbind',by(tora,tora$id,create_taper_mgs))
  }
  if(curva_media) {
    cm_eafil[,c('id','dap','ht')]<-NULL
    cm_eafil<-cm_eafil[order(cm_eafil$idestrato),]
  }else{
    cm_eafil<-NULL
  }
  return(list(fst=fst_eafil, cm=cm_eafil))
}
