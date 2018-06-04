# export.sav<-function(df,datafile,codefile,sep="\t",dec=",",na="",quote=T,drop.factor=F){
#   
#   if (drop.factor){
#     variable<-NULL
#     for (i in 1:ncol(df)){
#       if (is.factor(df[,i])) variable<-c(variable,i)
#     }
#     df<-df[,-variable]
#   }
#   
#   write.table(df,datafile,row.names=F,col.names=T,na=na,sep=sep,dec=dec,quote=quote)
#   
#   sink(codefile)
#   for(i in 1:ncol(df)){
#     if (!is.null(attr(df[,i],"var.lab"))) cat(paste("VAR LAB ",names(df[i]),"'",attr(df[,i],"var.lab"),"'.",sep=""),fill=T)
#     if (!is.null(attr(df[,i],"val.lab"))){
#       labels<-paste(list.val(df[,i])[,2],"'",list.val(df[,i])[,1],"'",sep="")
#       cat(paste("VAL LAB",names(df)[i]),fill=T)
#       cat(c(labels[-length(labels)],paste(labels[length(labels)],".",sep="")),fill=T)
#     }
#   }
#   sink()
# }

export.sav<-function(df,path){
  require(haven)
  for(i in names(df)){
    if(!is.null(attr(df[,i],"val.lab"))){
      df[,i]<-labelled(df[,i],labels=attr(df[,i],"val.lab"))
      attr(df[,i],"val.lab")<-NULL
    }
    if(!is.null(attr(df[,i],"var.lab"))){
      attr(df[,i],"label")<-attr(df[,i],"var.lab")
      attr(df[,i],"var.lab")<-NULL
    }
  }
  write_sav(df,path)
}

.weighted.var <- function(x,w,std=F) {
  if (missing(w))w<-rep(1,length(x))     
  n<-sum(ifelse(!is.na(x),w,NA),na.rm=T)
  media<-sum(x*w,na.rm=T)/n
  val2<-sum((x^2*w),na.rm=T)
  m2<-(val2)-(2*n*media^2)+(media^2*n)
  var<-m2/(n-1)
  if (std==F) return(var)
  return(sqrt(var))
}


standarized<-function(df,w){
  
  if (missing(w)) w<-rep(1,nrow(df))
  
  df1<-df[complete.cases(w),,drop=F]
  w<-w[complete.cases(w)]
  
  if(ncol(df)>1){
    df2<-sweep(df1,2,apply(df1,2,weighted.mean,w=w,na.rm=T),FUN = "-")
    df3<-sweep(df2,2,apply(df1,2,.weighted.var,w=w,std=T),FUN = "/")
 
    indices<-names(df3)
    copyattr<-function(x){
      v<-df3[,x]  
      attr(v,"var.lab")<-paste(attr(df[,x],"var.lab"),"(Zscore)")
      names(attr(v,"var.lab"))<-paste("z",x,sep="")
      return(v)
    } 
    df3<-data.frame(sapply(indices,copyattr,simplify = F))    
    names(df3)<-paste("z",names(df),sep="")
    rownames(df3)<-rownames(df1)
    return(df3)
  }else{
    df2<-sweep(df1,2,weighted.mean(df1[,1],w=w,na.rm=T),FUN = "-")
    df3<-sweep(df2,2,.weighted.var(df1[,1],w=w,std=T),FUN = "/")
    nombre<-names(df3)
    attr(df3[,1],"var.lab")<-paste(attr(df[,1],"var.lab"),"(Zscore)")
    names(attr(df3[,1],"var.lab"))<-paste("z",nombre,sep="")    
    rownames(df3)<-rownames(df1)
    names(df3)<-paste("z",names(df),sep="")
    return(df3)
  }
  
}

print.tabla<-function (x){
  if (any(c("cross", "ctables", "cmeans", "multiresp") %in% 
            class(x))) {
    if (!"cmeans"%in%class(x)){
      tablap <- ifelse(x == 0, NA, x)
    }else{
      tablap<-x
    }
    tablap <- round(tablap, attr(x, "dec"))
    tablap <- ifelse(is.na(tablap), "", tablap)
    tablap <- cbind(rownames(tablap), tablap)
    rownames(tablap) <- NULL
    colnames(tablap) <- substr(colnames(tablap), 1, attr(x, 
                                                         "ancho"))
    
    cells <- nrow(x)/length(rownames(x)[!rownames(x) %in% 
                                          ""])
    separador <-paste(rep("-",400),collapse="") 
    
    long<-apply(tablap,2,function(k)max(nchar(as.character(k))))
    long<- pmax(long, nchar(as.character(colnames(tablap))))
    long<- ifelse(long < 3, 3, long)
    long[2:length(long)]<-max(long[2:length(long)])#Todas la columnas sean igual de anchas
    
    separadores <- NULL
    for (i in 1:length(long)) {
      separadores <- cbind(separadores, substr(separador, 
                                               1, long[i]))
    }
    tablaprint <- NULL
    saltos <- seq(cells, nrow(tablap), cells)
    for (i in saltos) {
      if (i == min(saltos)) {
        tablaprint <- rbind(tablaprint, tablap[1:i, ], 
                            separadores)
      }
      else {
        tablaprint <- rbind(tablaprint, tablap[(i - (cells - 
                                                       1)):i, ], separadores)
      }
    }
    alinear <- c("l", rep("r", ncol(tablaprint)))
    
    repeat{
      indices<-c(1:ncol(tablaprint))
      indicesp<-nchar(tablaprint[nrow(tablaprint),indices])+1
      indicesp[length(indicesp)]<-indicesp[length(indicesp)]+1
      indicesp<-which(cumsum(indicesp)<options()$width)
      
      indicesp[c(1,length(indicesp))]<-c(1,ncol(tablaprint))
      
      cat(attr(x, "title"))
      print(knitr::kable(tablaprint[, indicesp], format = "markdown", 
                         row.names = F, align = alinear, padding = 0))
      
      indicesp<-c(1,setdiff(indices,indicesp),ncol(tablaprint))
      if (length(indicesp)>2){ 
        cat("\n(Continued)\n\n")
        tablaprint<-tablaprint[,indicesp]
      }else{
        break
      }
    }    
    if (!is.null(attr(x, "resumen"))) 
      cat("\n", attr(x, "resumen"), "\n", fill = T)
    if (any(c("ctables", "multiresp") %in% class(x))) 
      cat("\n")
    if (c("cmeans") %in% class(x)) 
      cat("\n Cell contents:", attr(x, "pct"), "\n", fill = T)
  }
  if (any(c("means", "freq", "desc") %in% class(x))) {
    tablaprint <- x
    attributes(tablaprint) <- NULL
    dim(tablaprint) <- dim(x)
    rownames(tablaprint) <- rownames(x)
    colnames(tablaprint) <- colnames(x)
    cat(attr(x, "title"), "\n", fill = T)
    print(tablaprint, na.print = "")
    if (!is.null(attr(x, "resumen"))) {
      cat("\n", attr(x, "resumen"), "\n", fill = T)
    }
    else {
      cat("\n")
    }
  }
  if ("list" %in% class(x)) {
    for (i in names(x)) {
      if (i == "desc") {
        tablaprint <- x[[i]]
        attributes(tablaprint) <- NULL
        dim(tablaprint) <- dim(x[[i]])
        rownames(tablaprint) <- rownames(x[[i]])
        colnames(tablaprint) <- colnames(x[[i]])
        cat(attr(x[[i]], "title"), "\n", fill = T)
        print(tablaprint, na.print = "")
      }
      else {
        tablaprint <- x[[i]]
        cat("\nPercentiles\n", fill = T)
        print(tablaprint, na.print = "")
        cat("\n", fill = T)
      }
    }
  }
}

freq<-function(x,...){
  UseMethod("freq",x)
}

.frequencies<-function(x,w,order,dec=1,selectcol){
   if (missing(w)) w<-rep(1,length(x))
  tabla<-suppressWarnings(rowsum(w,x,na.rm=T))
  n<-N<-sum(tabla)
  
  if (is.numeric(x) & !is.null(attr(x, "val.lab"))) {
    indices <- as.numeric(rownames(tabla))
    etiquetas <- indices %in% attr(x, "val.lab")
    etiquetas <- sort(c(attr(x, "val.lab"), indices[etiquetas==F]))
    etiquetas <- names(etiquetas[etiquetas %in% indices])
    
    if (is.na(indices[length(indices)])){
      indices[length(indices)]<-""
      etiquetas<-c(etiquetas,"missing")
    }
    rownames(tabla) <- paste(indices, etiquetas)
  }
  rownames(tabla)[is.na(rownames(tabla))]<-" missing"
  tabla<-cbind(tabla,prop.table(tabla)*100)
  if(any(rownames(tabla)==" missing")) {
    tabla<-cbind(tabla,c(prop.table(tabla[-nrow(tabla),1])*100,NA))
    n<-sum(tabla[-nrow(tabla),1])
  }else{
    tabla<-cbind(tabla,prop.table(tabla[,1])*100)
  }

 if (!missing(order)) {
    if (order == "a") 
      tabla <- tabla[order(tabla[, 3]), ]
    if (order == "d") 
      tabla <- tabla[order(-tabla[, 3]), ]
  }
  
  tabla<-cbind(tabla,cumsum(tabla[,3]))
  colnames(tabla)<-c("Frequency","Percent","Valid Pct","Cum Pct")
  tabla[,1]<-round(tabla[,1],0)
  tabla[,2:4]<-round(tabla[,2:4],dec)
  
   if (!missing(selectcol)) 
    tabla <- tabla[, selectcol, drop = F]
  if (is.null(attr(x, "var.lab"))) 
    attr(x, "var.lab") <- ""
  tabla <- structure(tabla, class = c("tabla", "matrix", "freq"))
  attr(tabla, "title") <- paste(names(as.list(attr(x, "var.lab"))), 
                                attr(x, "var.lab"))
  attr(tabla, "resumen") <- paste("Total Cases:", round(N,0), "Valid Cases:", round(n, 0))
  return(tabla)
}


freq.character<-function(x,...) .frequencies(x,...)
freq.factor<-function(x,...) .frequencies(x,...)
freq.numeric<-function(x,...) .frequencies(x,...)
freq.logical<-function(x,...) .frequencies(x,...)
freq.data.frame<-function(x,...)lapply(x,.frequencies,...)


desc<- function (x, ...) {
   UseMethod("desc", x)
 }
 
.descriptives<-function(x,w,ntiles=1,stat=c("Mean","Std.Dev","Minimum","Maximum","Valid.N"),dec=3){
    
  if (is.null(x)) stop("La variable no existe")
  
  if (missing(w))w<-rep(1,length(x))     
  stat<-c("Mean","Median","Mode","Variance","Std.Dev","S.E.Mean","Skewness","Kurtosis","Minimum","Maximum","Range","Sum","Valid.N")[c("Mean","Median","Mode","Variance","Std.Dev","S.E.Mean","Skewness","Kurtosis","Minimum","Maximum","Range","Sum","Valid.N")%in%unique(stat)]
      
  n<-sum(ifelse(!is.na(x),w,NA),na.rm=T)
  media<-sum(x*w,na.rm=T)/n
  
  val2<-sum((x^2*w),na.rm=T)
  m2<-(val2)-(2*n*media^2)+(media^2*n)
  varianza<-m2/(n-1)
  
  desviacion<-sqrt(varianza)
  s.e.media<-desviacion/sqrt(n)
  
  if (any(c("Skewness","Kurtosis")%in%stat)){
  val3<-sum((x^3*w),na.rm=T)
  m3<-(val3)-(3*media*val2)+(3*media^3*n)-(media^3*n)
  asimetria<-(n*m3)/((n-1)*(n-2)*desviacion^3)
    
  val4<-sum(x^4*w,na.rm=T)
  m4<-(val4)-(4*media*val3)+(6*media^2*val2)-(4*media^4*n)+(n*media^4)
  curtosis<-((n*(n-1)*m4)-(3*m2*m2*(n-1)))/((n-1)*(n-2)*(n-3)*desviacion^4)
  }
  
  if (any(c("Minimum","Maximum","Range")%in%stat)){
  minimo<-min(x,na.rm=T)
  maximo<-max(x,na.rm=T)
  rango<-maximo-minimo
  }
  
  suma<-sum(x*w,na.rm=T)
  
  if (any(c("Median","Mode")%in%stat)|ntiles>1){
  percents<-prop.table(xtabs(w~x))*100
  cumpercents<-cumsum(percents)
  mediana<-as.numeric(names(which(cumpercents>=50)[1]))
  moda<-as.numeric(names(which.max(percents)[1]))
  }

tabla<-NULL
if ("Mean"%in%stat) tabla<-cbind(tabla,Mean=round(media,dec))
if ("Mode"%in%stat) tabla<-cbind(tabla,Mode=round(moda,dec))
if ("Median"%in%stat) tabla<-cbind(tabla,Median=round(mediana,dec))
if ("Variance"%in%stat) tabla<-cbind(tabla,Variance=round(varianza,dec))
if ("Std.Dev"%in%stat) tabla<-cbind(tabla,Std.Dev=round(desviacion,dec))
if ("S.E.Mean"%in%stat) tabla<-cbind(tabla,S.E.Mean=round(s.e.media,dec))
if ("Skewness"%in%stat) tabla<-cbind(tabla,Skewness=round(asimetria,dec))
if ("Kurtosis"%in%stat) tabla<-cbind(tabla,Kurtosis=round(curtosis,dec))
if ("Minimum"%in%stat) tabla<-cbind(tabla,Minimum=round(minimo,dec))
if ("Maximum"%in%stat) tabla<-cbind(tabla,Maximum=round(maximo,dec))
if ("Range"%in%stat) tabla<-cbind(tabla,Range=round(rango,dec))
if ("Sum"%in%stat) tabla<-cbind(tabla,Sum=round(suma,dec))
if ("Valid.N"%in%stat) tabla<-cbind(tabla,Valid.N=round(n,0))

 if (ntiles>1){
    indices<-1/ntiles
    indices<-seq(indices,1,indices)*100
    ntil<-NULL
    label<-NULL
    for (i in 1:(length(indices)-1)){
      ntil<-c(ntil,as.numeric(names(which(cumpercents>=indices[i])[1])))
      label<-c(label,paste("To ",as.character(round(indices[i],2)),"%",sep=""))
    }
    names(ntil)<-label
    ntil<-round(ntil,3)
  }
  
  if (is.null(attr(x,"var.lab"))) attr(x,"var.lab")<-""
  rownames(tabla)<-paste(names(as.list(attr(x,"var.lab"))),attr(x,"var.lab"))
      
  tabla<-structure(tabla,class=c("tabla","matrix","desc"))
  attr(tabla,"title")<-"Descriptives"
  if (ntiles>1) {
    lista<-list(desc=tabla,percentiles=ntil)
    lista<-structure(lista,class=c("tabla","list"))
    return(lista)
  }
  return(tabla)
}

desc.logical<-function(x,...) .descriptives(x,...)
desc.numeric<-function(x,...) .descriptives(x,...)
desc.matrix<-function(x,...,stat=c("Mean","Std.Dev","Minimum","Maximum","Valid.N")){
 tabla <- apply(x, 2, .descriptives,...,stat = stat)
  if(length(stat)>1){
  tabla<-t(tabla)
  }else{
  tabla<-t(t(tabla))
  }
  colnames(tabla) <- stat
  tabla
  }
desc.data.frame<-function(x,...){
clase<-sapply(x,class)
if(any(clase%in%c("character","factor"))){
    eliminar<-which(clase%in%c("character","factor"))
    cat("Not numerics variables:\n")
    print(names(clase)[eliminar])
    cat("\n")
    lapply(x[,-eliminar],.descriptives,...)
}else{
lapply(x,.descriptives,...)
}}

means<-function(x,y,w,stat=c("Mean","Std.Dev","Valid.N"),Totrow=T,dec=2,selectrow){
  
  if (is.null(x)) stop("La variable x no existe")
  if (is.null(y)) stop("La variable y no existe")
  
  
  stat<-c("Mean","Variance","Std.Dev","S.E.Mean","Skewness","Kurtosis","Minimum","Maximum","Range","Sum","Valid.N")[c("Mean","Variance","Std.Dev","S.E.Mean","Skewness","Kurtosis","Minimum","Maximum","Range","Sum","Valid.N")%in%unique(stat)]
  
  if ("Valid.N"%in%stat){
  validn<-T
  }else{
  stat<-c(stat,"Valid.N")  
  validn<-F
  }
  
  
  if (missing(w)) w<-rep(1,length(x))
  N<-round(sum(w,na.rm=T),0)
  validos<-complete.cases(cbind(x,y))
  n<-round(sum(w[validos],na.rm=T),0)
  
  values<-sort(unique(y))  
  tabla<-sapply(values,function(k)desc(x[y==k],w[y==k],stat=stat,dec=dec))

  
   tabla<-t(tabla)
   dim(tabla)<-c(length(values),length(stat))
   rownames(tabla)<-values
   colnames(tabla)<-stat

  Total<-desc(x[validos],w[validos],stat=stat)
  rownames(Total)<-"Total"
    
  
  if (is.numeric(x) & !is.null(attr(y, "val.lab"))) {
      indices <- as.numeric(rownames(tabla))
      etiquetas<-indices%in%attr(y,"val.lab")
      etiquetas<-sort(c(attr(y,"val.lab"),indices[etiquetas==F]))
      etiquetas<-names(etiquetas[etiquetas%in%indices])
      rownames(tabla)<-paste(indices,etiquetas)   
    }
        
  tabla<-rbind(tabla,Total)  
  tabla<-round(tabla,dec)
  tabla<-tabla[tabla[,"Valid.N"]!=0,,drop=F]
  if (validn==F) tabla<-tabla[,-ncol(tabla),drop=F]
  if (!Totrow) tabla<-tabla[-nrow(tabla),,drop=F]
  
  if (!missing(selectrow)) tabla<-tabla[selectrow,,drop=F]
  
  if (is.null(attr(x,"var.lab"))) attr(x,"var.lab")<-""
  if (is.null(attr(y,"var.lab"))) attr(y,"var.lab")<-""  
  
  tabla<-structure(tabla,class=c("tabla","matrix","means"))
  attr(tabla,"title")<-paste("Descriptives of",names(as.list(attr(x,"var.lab"))),attr(x,"var.lab"),"by",names(as.list(attr(y,"var.lab"))),attr(y,"var.lab"))
  attr(tabla,"resumen")<-paste("Total Cases =",N,"Missing Cases =",N-n)
  return(tabla)
}

filtro<-function (df, ...) {
  
  bd <- droplevels(subset(df, ...))
  dimbd<-dim(bd)
  rownamesbd<-rownames(bd)
  if (nrow(bd) == 0) 
    stop("La base de datos filtrada no contiene ninguna observacion")
  if (ncol(bd) == 0) 
    stop("La base de datos filtrada no contiene ninguna variable")
  indices <- names(df)[names(df) %in% names(bd)]
  copyattr <- function(x, y) {
    attr(y,"var.lab") <- attr(x,"var.lab")
    attr(y,"val.lab") <- attr(x,"val.lab")
    return(y)
  }
  bd <- as.data.frame(mapply(copyattr, df[, indices,drop=F], bd[,indices,drop=F], 
                               SIMPLIFY = F),stringsAsFactors =F)
  if(any(dim(bd)==dimbd)==F) bd<-t(bd)
  rownames(bd) <- rownamesbd
  return(bd)
}

list.var<-function(str,df){
  var<-names(df)[grep(str,names(df))]
  if (length(var)==0) stop("No hay coincidencias")
  Names<-NULL
  Labels<-NULL
  for (i in 1:length(var)){
    Names<-c(Names,var[i])
    if(is.null(attr(df[,var[i]],"var.lab"))) {
      Labels<-c(Labels," ")
    }else{
      Labels<-c(Labels,attr(df[,var[i]],"var.lab"))
    }
  }
  tabla<-cbind(Names,Labels)
  rownames(tabla)<-grep(str,names(df))
  print(tabla,quote=F)  
}

contar<-function (bbdd, v,thru=F) {

valor <- rep(0,nrow(bbdd))
  
if (length(v)==2 & thru==T)  {
    for (i in 1:ncol(bbdd)) {
      valor[bbdd[, i] %thru% v ] <- valor[bbdd[,i] %thru% v] + 1
    }
    return(valor)
  }else{
    for (i in 1:ncol(bbdd)) {
      valor[bbdd[, i] %in% v & !is.na(bbdd[, i])] <- valor[bbdd[,i] %in% v & !is.na(bbdd[, i])] + 1
    }
    return(valor)
  }
}

recode<-function(x,...){
  
  if (is.null(x)) stop("La variable no existe")
  grupos<-list(...)
  if (!is.list(grupos)) stop("Los grupos han de ser introducidos como objeto list()")
  
  valores<-sapply(grupos,function(k)sum(names(k)=="v"))
  if (any(valores!=1)){
    stop("todos los grupos del list han de contener el valor (v)")
  }else{
    valores<-as.numeric(sapply(grupos,function(k)k[names(k)=="v"]))  
  }  
  
  varlab<-attr(x,"var.lab")  
  vallab<-attr(x,"val.lab")  
  grupos<-lapply(grupos,function(k)k[!names(k)=="v"])
  
  if(length(valores)!=length(grupos)) stop("Las recodificaciones no se han introducido correctamente")
  
  if(!class(x)%in%c("factor","character") & any(sapply(grupos,class)%in%c("factor","character")))
   stop("Esta intentado recodificar variable numerica")
  
  if(class(x)%in%c("factor","character")){
  if(any(sapply(grupos,class)%in%c("numeric","integer"))) stop("Esta intentando recodificar una variable numerica")
  x<-as.character(x)
  } 
  
  variable<-x
  for(i in 1:length(grupos)){
    variable[x%in%grupos[[i]]]<-valores[i]
    if (any(is.na(grupos[[i]]))) variable[is.na(x)]<-valores[i]
  }
  variable<-as.numeric(variable)
  
  if(!is.null(varlab)) {
    nombre<-names(varlab)
    attr(variable,"var.lab")<-varlab
    names(attr(variable,"var.lab"))<-nombre}
  
  if(!is.null(vallab)) attr(variable,"val.lab")<-vallab  
  
  return(variable)
}

compfactor<-function(x){
  
  if (is.null(x)) stop("La variable no existe")
  if (is.null(attr(x,"val.lab"))) stop("La variable no tiene value labels")
  
  indices <- sort(unique(x))
  etiquetas<-indices%in%attr(x,"val.lab")
  etiquetas<-sort(c(attr(x,"val.lab"),indices[etiquetas==F]))
  etiquetas<-names(etiquetas[etiquetas%in%indices])
  
  etiquetas[etiquetas==""]<-indices[etiquetas==""]
      
  variable<-factor(x,levels=sort(unique(x)),labels=etiquetas)  
  attr(variable,"var.lab")<-attr(x,"var.lab")
  return(variable)
}

list.val<-function(x){
if (is.null(x)) stop("La variable no existe")
if (is.null(attr(x, "val.lab"))) stop("La variable no tiene definido ningun nivel")
bbdd<-data.frame(Label=names(attr(x, "val.lab")),Value = attr(x, "val.lab"),stringsAsFactors = T)
rownames(bbdd)<-NULL
return(bbdd)
}

cross<-function (x, y, w, cells = "count", dec = 1, order, ancho = 12, 
          Totrow = T, Totcol = T, selectrow, selectcol) 
{
  
if (is.null(x)) 
    stop("La variable x no existe")
  if (is.null(y)) 
    stop("La variable y no existe")
  if (missing(w)) 
    w <- rep(1, length(x))
  
  cells <- c("count", "row", "col", "tot")[which(c("count", 
                                                   "row", "col", "tot") %in% unique(cells))]
  absolutos <- tapply(w, list(x, y), sum, na.rm = T)
  absolutos[is.na(absolutos)] <- 0
  
  
  absolutos <- absolutos[!(rowSums(absolutos) == 0), !colSums(absolutos) == 
                           0, drop = F]
  if (all(dim(absolutos) == c(0, 0)) == T)  stop("Cross is empty")
  
  
  colpct<-prop.table(absolutos,2)*100
  filapct<-prop.table(absolutos,1)*100
  totpct<-prop.table(absolutos)*100
  
  if (is.numeric(x) & !is.null(attr(x, "val.lab"))) {
    indices <- as.numeric(rownames(absolutos))
    etiquetas <- indices %in% attr(x, "val.lab")
    etiquetas <- sort(c(attr(x, "val.lab"), indices[etiquetas==F]))
    etiquetas <- names(etiquetas[etiquetas %in% indices])
    rownames(absolutos) <- paste(indices, etiquetas)
  }
  
  if (is.numeric(y) & !is.null(attr(y, "val.lab"))) {
    indices <- as.numeric(colnames(absolutos))
    etiquetas <- indices %in% attr(y, "val.lab")
    etiquetas <- sort(c(attr(y, "val.lab"), indices[etiquetas==F]))
    etiquetas <- names(etiquetas[etiquetas %in% indices])
    colnames(absolutos) <- paste(indices, etiquetas)
  }
  
  
if(!missing(order)){
    if(order%in%c("a","d")){
      if(order=="a"){
        orden<-order(rowSums(absolutos))
        absolutos<-absolutos[orden,]
        colpct<-colpct[orden,]
        filapct<-filapct[orden,]
        totpct<-totpct[orden,]
        
      }else{
        orden<-order(-rowSums(tabla))
        absolutos<-absolutos[orden,]
        colpct<-colpct[orden,]
        filapct<-filapct[orden,]
        totpct<-totpct[orden,]
      }  
    }
  }
  
  
  
absolutos <- addmargins(absolutos, FUN = list(Total = sum),quiet = T)
colpct <- addmargins(colpct,1, FUN = list(Total = sum), quiet = T)
filapct <- addmargins(filapct,2, FUN = list(Total = sum),quiet = T)
totpct <- addmargins(totpct, FUN = list(Total = sum), quiet = T)
colpct<-cbind(colpct,Total=totpct[,ncol(totpct)])
filapct<-rbind(filapct,Total=totpct[nrow(totpct),])
absolutos <- round(absolutos, 0)

  tabla <- NULL
  for (i in 1:nrow(absolutos)) {
    tabla <- rbind(tabla, absolutos[i, ], filapct[i, ], 
                   colpct[i, ], totpct[i, ])
  }

rownames(tabla)<-paste(rep(rownames(absolutos), each = 4))
pct<-rep(c("count","row", "col", "tot"),nrow(absolutos))
tabla<-tabla[which(pct %in% cells),,drop=F]
rownames(tabla)[-c(seq(1, length(rownames(tabla)), length(cells)))] <- ""  
  pct <- cells
  pct <- ifelse(pct != "count", paste(pct, "%", sep = ""), 
                pct)
  pct <- paste(pct, collapse = ", ")
  if (is.null(attr(x, "var.lab"))) 
    attr(x, "var.lab") <- ""
  if (is.null(attr(y, "var.lab"))) 
    attr(y, "var.lab") <- ""
  if (Totcol == F) 
    tabla <- tabla[, -ncol(tabla), drop = F]
  if (Totrow == F) 
    tabla <- tabla[-(nrow(tabla):((nrow(tabla) - (length(cells) - 
                                                    1)))), , drop = F]
  if (!missing(selectrow)) 
    tabla <- tabla[selectrow, , drop = F]
  if (!missing(selectcol)) 
    tabla <- tabla[, selectcol, drop = F]
  if (is.null(attr(x, "var.lab"))) 
    attr(x, "var.lab") <- ""
  if (is.null(attr(y, "var.lab"))) 
    attr(y, "var.lab") <- ""
  tabla <- structure(tabla, class = c("tabla", "matrix", "cross"))
  attr(tabla, "formato") <- cumsum(rep(length(cells), nrow(absolutos))[-1])
  attr(tabla, "pct") <- pct
  attr(tabla, "resumen") <- paste("Number of Missing Observations:", 
                                  round(sum(w, na.rm = T), 0) - absolutos[nrow(absolutos), 
                                                                          ncol(absolutos)])
  attr(tabla, "title") <- paste(names(as.list(attr(x, "var.lab"))), 
                                attr(x, "var.lab"), "by", names(as.list(attr(y, "var.lab"))), 
                                attr(y, "var.lab"))
  attr(tabla, "ancho") <- ancho
  attr(tabla, "dec") <- dec
  
  return(tabla)
}


compnumeric<-function(x){
  if (is.null(x)) stop("La variable no existe")    
  if(!is.factor(x)) stop("La variable ha de ser factor")
  etiqueta<-attr(x,"var.lab")
  x<-droplevels(x)
  variable<-as.numeric(x)
  attr(variable,"var.lab")<-etiqueta
  
  value.labs<-c(1:length(levels(x)))
  names(value.labs)<-levels(x)
  attr(variable,"val.lab")<-value.labs
  return(variable)
}

cmeans<-function (x, y, w, dec = 2, stat = "Mean", title = "", ancho = 12, 
  Totcol = T, selectrow, selectcol) 
{
  if (missing(w)) 
    w <- rep(1, length(x[[1]]))
  stat <- c("Mean", "Variance", "Std.Dev", "S.E.Mean", "Skewness", 
    "Kurtosis", "Minimum", "Maximum", "Range", "Sum", "Valid.N")[c("Mean", 
    "Variance", "Std.Dev", "S.E.Mean", "Skewness", "Kurtosis", 
    "Minimum", "Maximum", "Range", "Sum", "Valid.N") %in% 
    unique(stat)]
  if (missing(y)) {
    y <- rep("siny", length(x[[1]]))
  }
  if (is.null(y)) 
    stop("La variable y no existe")
  tabla <- NULL
  etiquetas <- NULL
  for (i in 1:ncol(x)) {
    if (is.null(attr(x[,i], "var.lab"))) 
      attr(x[,i], "var.lab") <- ""
    etiquetas <- c(etiquetas, rep(paste(names(attr(x[,i], "var.lab")), 
    unname(attr(x[,i], "var.lab"))), length(stat)))
    tabla <- cbind(tabla, means(x[,i], y, w, dec = dec, stat = stat))
  }
  tabla <- t(tabla)
  etiquetas[-c(seq(1, length(etiquetas), length(stat)))] <- ""
  tabla <- cbind(` ` = etiquetas, Stat = rownames(tabla), 
    tabla)
  rownames(tabla) <- NULL
  if (length(unique(y) == 1) & all(unique(y) == "siny")) {
    tabla <- tabla[, -3, drop = F]
  }
  pct <- stat
  pct <- paste(pct, collapse = ", ")
  rownames(tabla) <- tabla[, 1]
  tabla <- tabla[, -c(1, 2), drop = F]
  tabla <- apply(tabla, 2, function(x) as.numeric(x))
  if (class(tabla) != "matrix") {
    etiquetasy <- names(tabla)
    dim(tabla) <- c(1, length(tabla))
    colnames(tabla) <- etiquetasy
  }
  rownames(tabla) <- etiquetas
  if (Totcol == F) 
    tabla <- tabla[, -ncol(tabla), drop = F]
  if (!missing(selectrow)) 
    tabla <- tabla[selectrow, , drop = F]
  if (!missing(selectcol)) 
    tabla <- tabla[, selectcol, drop = F]
  if (is.null(attr(y, "var.lab"))) 
    attr(y, "var.lab") <- ""
  tabla <- structure(tabla, class = c("tabla", "matrix", "cmeans"))
  attr(tabla, "formato") <- cumsum(rep(length(stat), nrow(tabla)/length(stat))[-1])
  attr(tabla, "pct") <- pct
  attr(tabla, "title") <- paste(title, "by", names(as.list(attr(y, 
    "var.lab"))), attr(y, "var.lab"))
  attr(tabla, "ancho") <- ancho
  attr(tabla, "dec") <- dec
  return(tabla)
}

multiresp<-function(x,y,w,orden,resp=F,dec=1,cells="count",title="",ancho=12,Totrow=T,Totcol=T,selectrow,selectcol){
  
  
  if (missing(y)) {
    y<-rep("siny",length(x[[1]]))
    }
  
  if (is.null(y)) stop("La variable y no existe")  
  
  if (!is.list(x)) stop("Las variables han de ser introducidas como una lista")
  if (missing(w)) w<-rep(1,length(x[[1]]))
  cells<-c("count","row","col","tot")[which(c("count","row","col","tot")%in%unique(cells))]
  
  clase<-"error"
  if (all(sapply(x,function(x)class(x))=="factor")) clase="factor"
  if (all(sapply(x,function(x)class(x))=="numeric")) clase="numeric"
  
  if (clase=="error") stop("Las variables de la lista han de ser o todas factores o todas numericas")  
  
 
  absolutos<-NULL
  for (i in 1:length(x)){
    absolutos<-rbind(absolutos,tapply(w,list(x[[i]],y),sum,na.rm=T))
  }
  
  absolutos[is.na(absolutos)]<-0
  
  absplit<-split(absolutos,rownames(absolutos))
  absolutos<-t(sapply(absplit,function(x)colSums(matrix(x,ncol=ncol(absolutos)))))
  
  if (length(unique(y)==1) & all(unique(y)=="siny")) absolutos<-t(absolutos)
  
  colnames(absolutos)<-sort(unique(y))
  absolutos<-round(absolutos,0)    
  
  totalx<-rowSums(absolutos)
  totaly<-tapply(w,y,sum,na.rm=T)
  totaly<-ifelse(is.na(totaly),0,totaly)
  
if (is.numeric(x[[1]]) & !is.null(attr(x[[1]],"val.lab"))){
    indices <- as.numeric(rownames(absolutos))
    etiquetasx<-indices%in%attr(x[[1]],"val.lab")
    etiquetasx<-sort(c(attr(x[[1]],"val.lab"),indices[etiquetasx==F]))
    etiquetasx<-names(etiquetasx[etiquetasx%in%indices])
    etiquetasx<-c(paste(indices,etiquetasx),"Total")
} else {
    etiquetasx <- c(rownames(absolutos), "Total")
}
  
  if (is.numeric(y) & !is.null(attr(y,"val.lab"))){
    indices <- as.numeric(colnames(absolutos))
    etiquetasy<-indices%in%attr(y,"val.lab")
    etiquetasy<-sort(c(attr(y,"val.lab"),indices[etiquetasy==F]))
    etiquetasy<-names(etiquetasy[etiquetasy%in%indices])
    etiquetasy<-c("",paste(indices,etiquetasy),"Total")  
  } else {
    etiquetasy <- c("", colnames(absolutos), "Total")
  }
  
  if (!missing(orden)){
    if(orden=="d"){
      absolutos<-absolutos[rev(order(totalx)),]
      etiquetasx<-etiquetasx[c(rev(order(totalx)),length(totalx)+1)]
      totalx<-totalx[rev(order(totalx))]
    }else if(orden=="a"){
      absolutos<-absolutos[order(totalx),]
      etiquetasx<-etiquetasx[c(order(totalx),length(totalx)+1)]
      totalx<-totalx[order(totalx)]
    }
  }
  
if (resp==T){
    filapct<-prop.table(absolutos,1)*100
    colpct<-prop.table(absolutos,2)*100
    totpct<-prop.table(absolutos)*100
    
    totalx<-prop.table(rowSums(absolutos))*100
    totaly<-prop.table(colSums(absolutos))*100
    
    filapct<-cbind(filapct,rowSums(filapct))
    filapct<-rbind(filapct,c(totaly,100))
    
    
    colpct<-rbind(colpct,colSums(colpct))
    colpct<-cbind(colpct,c(totalx,100))
    
    totpct<-cbind(totpct,totalx)
    totpct<-rbind(totpct,c(totaly,100))
  }else{
    filapct<-sweep(absolutos,1,totalx,"/")*100
    colpct<-sweep(absolutos,2,totaly,"/")*100
    totpct<-(absolutos/sum(totaly))*100
    
    filapct<-cbind(filapct,rowSums(filapct))
    filapct<-rbind(filapct,c(prop.table(totaly)*100,100))
    colpct<-rbind(colpct,rep(100,ncol(colpct)))
    colpct<-cbind(colpct,c((totalx/sum(totaly))*100,100))
    
    totpct<-cbind(totpct,c((totalx/sum(totaly))*100))
    totpct<-rbind(totpct,c(prop.table(totaly)*100,100))
}
  
  
  if (resp==T){
    absolutos<-cbind(absolutos,"Total"=rowSums(absolutos))
    absolutos<-rbind(absolutos,"Total"=colSums(absolutos))
  }else{
    absolutos<-cbind(absolutos,"Total"=rowSums(absolutos))
    absolutos<-rbind(absolutos,"Total"=c(totaly,sum(totaly)))
  }
  
  tabla<-NULL
  for (i in 1:nrow(absolutos)){
    tabla<-rbind(tabla,round(absolutos[i,],0),filapct[i,],colpct[i,],totpct[i,])
  }

  tabla<-cbind(rep(etiquetasx,each=4),rep(c("count","row","col","tot"),nrow(absolutos)),tabla)
  
  tabla<-tabla[which(tabla[,2]%in%cells),-2]
  etiquetasx<-tabla[,1]
  etiquetasx[-c(seq(1,length(etiquetasx),length(cells)))]<-""
  tabla[,1]<-etiquetasx
  
  colnames(tabla)<-etiquetasy
  
    
if (length(unique(y)==1) & all(unique(y)=="siny")) {
  tabla<-tabla[,-2,drop=F]
} 


pct<-cells
pct<-ifelse(pct!="count",paste(pct,"%",sep=""),pct)
pct<-paste(pct,collapse=", ")

rownames(tabla)<-tabla[,1]
tabla<-tabla[,-1,drop=F]
tabla<-apply(tabla,2,function(x)as.numeric(x))
rownames(tabla)<-etiquetasx      

if (Totcol==F) tabla<-tabla[,-ncol(tabla),drop=F]
if (Totrow==F) tabla<-tabla[-(nrow(tabla):((nrow(tabla)-(length(cells)-1)))),,drop=F]

if (!missing(selectrow)) tabla<-tabla[selectrow,,drop=F]
if (!missing(selectcol)) tabla<-tabla[,selectcol,drop=F]

if (is.null(attr(y,"var.lab"))) attr(y,"var.lab")<-""

tabla<-structure(tabla,class=c("tabla","matrix","multiresp"))
attr(tabla,"formato")<-cumsum(rep(length(cells),nrow(tabla)/length(cells))[-1])
attr(tabla,"pct")<-pct
attr(tabla,"title")<-paste(title,"by",names(as.list(attr(y,"var.lab"))),attr(y,"var.lab"))
attr(tabla,"ancho")<-ancho
attr(tabla,"dec")<-dec
return(tabla)
}

spssdef<-function (df,to.data.frame=T) { 
  if ("tbl_df"%in%class(df)){
    for (i in names(df)) {
    a <- attr(df[[i]], "label")
    b <- sort(attr(df[[i]], "labels"))
    if(class(df[[i]])=="labelled") df[[i]]<-as.numeric(df[[i]])
    attr(df[[i]], "var.lab") <- a
    if (is.null(a) | length(a)>1) 
      attr(df[[i]], "var.lab") <- ""
    names(attr(df[[i]], "var.lab")) <- toupper(i)
    attr(df[[i]], "val.lab") <- b
}
  }else{
  for (i in names(df)) {
    attr(df[, i], "var.lab") <- attributes(df)$variable.labels[i]
    attr(df[, i], "val.lab") <- sort(attr(df[, i], "value.labels"))
    attr(df[, i], "value.labels") <- NULL
  }
    }
  attr(df, "variable.labels") <- NULL
  names(df)<-tolower(names(df))
  if(to.data.frame==T) df<-as.data.frame(df)
  return(df)
}


.meansby<-function(x,y,w){
  return(tapply(x*w,y,sum,na.rm=T)/tapply(w[complete.cases(x)],y[complete.cases(x)],sum,na.rm=T))
}

k.means<-function (x, centers, w, iter = 10, initial){
  
  noval<-si(nmis(x)==ncol(x),T,F)
  if (missing(w))
    w <- rep(1, nrow(x))

  if (missing(initial)) {
    centroides <-sample(which(complete.cases(x)),centers)
    centroides <-as.matrix(x[centroides,,drop=F])
  }
  else {
    if (!is.matrix(initial))
      stop("Los valores iniciales han de ser una matrix")
    if (ncol(initial) != ncol(x))
      stop("El numero de columnas de initial es diferente a de x")
    if (nrow(initial) != centers)
      stop("El numero de filas de initial tiene que ser igual al numero de clusters (centers)")
    centroides <- initial
  }
  centroidesinit <- centroides
  
  mdist <- apply(centroides, 1, function(k) sqrt(rowSums(sweep(x,
                                                               2, k, FUN = "-")^2,na.rm=T)))
  
  grupofinal <- apply(t(mdist), 2, which.min)
  condicion <- F
  contador <- 0
  
  while (condicion == F) {
    grupoinicial <- grupofinal
    centroides <- apply(x, 2, means, y = grupoinicial, w = w,
                        Totrow = F, stat = "Mean")
    mdist <- apply(centroides, 1, function(k) sqrt(rowSums(sweep(x,
                                                                 2, k, FUN = "-")^2,na.rm=T)))
    
    grupofinal <- apply(t(mdist), 2, which.min)
    condicion <- all(grupoinicial == grupofinal, na.rm = T)
    contador <- contador + 1
    if (contador == iter) {
      warning("Se han alcanzado el numero maximo de iteraciones:",
              iter, call. = F)
      break
    }
  }
  grupofinal[noval]<-NA
  grupofinal[is.na(w)]<-NA
  
  centroides <- apply(x, 2, means, y = grupofinal, w = w,
                      Totrow = F, stat = "Mean",dec=Inf)
  etiquetas <- NULL
  for (i in 1:ncol(x)) {
    if (is.null(attr(x[, i], "var.lab")))
      attr(x[, i], "var.lab") <- ""
    etiquetas <- c(etiquetas, attr(x[, i], "var.lab"))
  }
  colnames(centroides) <- paste(names(x), etiquetas)
  rownames(centroides) <- rownames(c(1:centers), do.NULL = F,
                                   prefix = "Grupo.")
  centroides <- t(centroides)
  colnames(centroidesinit) <- paste(names(x), etiquetas)
  rownames(centroidesinit) <- rownames(c(1:centers), do.NULL = F,
                                       prefix = "Grupo.")
  centroidesinit <- t(centroidesinit)
  
  dim(grupofinal) <- c(length(grupofinal), 1)
  rownames(grupofinal) <- rownames(x)
  
  resultados <- list(Centros_Iniciales = centroidesinit, Grupos = grupofinal,
                     Centros_Finales = centroides, Iterations = contador)
  return(resultados)
}

covar<-function(x,w){
if (ncol(x) < 2) 
    stop("Need array or data.frame with 2 or more columns")
  if (missing(w)) 
    w <- rep(1, nrow(x))
  covarianza <- NULL
  colum <- ncol(x)
  for (i in 1:colum) {
    for (j in 1:colum) {
      colums <- c(i, j)
      v <- complete.cases(x[, colums])
      medias <- apply(x[v, colums], 2, weighted.mean, 
        w = w[v], na.rm = T)
      
      a<-x[v,colums][,1]-medias[1]
      b<-x[v,colums][,2]-medias[2]
      numerador<-sum(a*b*w[v])
      denominador <-(sum(w[v]) - 1)
      covarianza <- c(covarianza, numerador/denominador)
    }
  }
  covarianza <- matrix(covarianza, byrow = T, nrow = colum)
  rownames(covarianza) <- colnames(covarianza) <- colnames(x)
  if (is.data.frame(x)) {
    rownames(covarianza) <- paste(rownames(covarianza), 
      sapply(x, function(k) attr(k, "var.lab")))
    rownames(covarianza) <- sub("NULL", "", rownames(covarianza))
    rownames(covarianza) <- sub("[[:blank:]]+$", "", rownames(covarianza))
  }
  return(covarianza)
}

correl<-function(x,w){
  if (missing(w))  w<-rep(1,nrow(x))
  correlaciones<-NULL
  colum<-ncol(x)
  for(i in 1:colum){
    for(j in 1:colum){
    colums <- c(i, j)
    v <- complete.cases(x[, colums])
    medias <- apply(x[v, colums], 2, weighted.mean, w = w[v],na.rm = T)
    a<-x[v, colums][,1] - medias[1]
    b<-x[v, colums][,2] - medias[2]
    numerador<-sum(a* b*w[v])
    denominador<-sqrt(sum((a^2)*w[v])) * sqrt(sum((b^2)*w[v]))
    correlaciones <- c(correlaciones,numerador/denominador)
    }
  }
  correlaciones<-matrix(correlaciones,byrow = T,nrow=colum)
  rownames(correlaciones)<-colnames(correlaciones)<-colnames(x)
  if(is.data.frame(x)){
    rownames(correlaciones)<-paste(rownames(correlaciones),sapply(x,function(k)attr(k,"var.lab")))
    rownames(correlaciones)<-sub("NULL","",rownames(correlaciones))
    rownames(correlaciones)<-sub("[[:blank:]]+$","",rownames(correlaciones))
  }
  return(correlaciones)
}


.sct<-function(x,w){
  sct<-(x-weighted.mean(x,w=w,na.rm=T))^2
  sct<-sum(sct*w,na.rm=T)
  return(sct)  
}

.scd<-function(x,y,w){
  matriz<-data.frame(x,w,y)
  matriz<-split(matriz,y)
  matriz<-lapply(matriz,function(x)cbind(apply(x,2,function(j)j-weighted.mean(j,w=x[,2],na.rm=T)),x[,2])[,-c(2,3)])
  matriz<-lapply(matriz,function(x)cbind(t(apply(x,1,function(j)j^2)),x[,2])[,-2])
  matriz<-lapply(matriz,function(x)cbind(sweep(x,1,x[,2],FUN="*"),x[,2])[,-2])
  scd<-t(sapply(matriz,function(x)apply(x,2,sum,na.rm=T)))[,-2]
  return(scd)
}

.sce<-function(x,y,w){
  matriz<-cbind(x,y)  
  matriz<-apply(matriz,2,.meansby,y=y,w=w)[,-2]
  matriz<-(matriz-weighted.mean(x,w=w,na.rm=T))^2
  matriz<-matriz*tapply(w[complete.cases(x)],y[complete.cases(x)],sum,na.rm=T)
  sce<-matriz
  return(sce)
}

ctables<-function(x,y,w,niveles,orden,cells="count",title="",ancho=12,dec=1,Totrow=T,Totcol=T,selectrow,selectcol){
  
  cells<-c("count","row","col","tot")[which(c("count","row","col","tot")%in%unique(cells))]
  absolutos<-NULL
  etiquetas<-NULL
  
  if (missing(y)) {
    y<-rep("siny",length(x[[1]]))
  }
  
  if (is.null(y)) stop("La variable y no existe")  
  
  if (!is.list(x)) stop("Las variables han de ser introducidas como una lista")
  if (missing(w)) w<-rep(1,length(x[[1]]))
  
  
  clase<-"error"
  if (all(sapply(x,function(x)class(x))=="factor")) clase="factor"
  if (all(sapply(x,function(x)class(x))=="numeric")) clase="numeric"
  
  
  if (clase=="error") stop("Las variables de la lista han de ser o todas factores o todas numericas")
    
  for (i in 1:length(x)){
    tablita<-tapply(w,list(x[[i]],y),sum,na.rm=T)
    tablita[is.na(tablita)]<-0            
    
    if (is.null(attr(x[[i]],"var.lab"))) attr(x[[i]],"var.lab")<-""
    
    if (is.null(attr(x[[i]],"val.lab"))) {
    etiqueta<-paste(names(attr(x[[i]],"var.lab")),attr(x[[i]],"var.lab"),rownames(tablita))
    }else{
      indices <- as.numeric(rownames(tablita))
      etiqueta<-indices%in%attr(x[[i]],"val.lab")
      etiqueta<-sort(c(attr(x[[i]],"val.lab"),indices[etiqueta==F]))
      etiqueta<-names(etiqueta[etiqueta%in%indices])
      etiqueta<-paste(names(attr(x[[i]],"var.lab")),attr(x[[i]],"var.lab"),etiqueta)#si introduzco el objeto indices,agrego el valor aletiquetado      
}
absolutos<-rbind(absolutos,tablita)
etiquetas<-c(etiquetas,etiqueta)
}

  if (!missing(orden) & !missing(niveles)){
    totalx<-prop.table(rowSums(absolutos))*100
    if(orden=="d"){
      absolutos<-absolutos[rev(order(totalx)),,drop=F]
      etiquetas<-etiquetas[rev(order(totalx))]
    }
    if(orden=="a"){
      absolutos<-absolutos[order(totalx),,drop=F]
      etiquetas<-etiquetas[order(totalx)]
    }
  }
  
etiquetas<-c(etiquetas,"Total")

    
if (is.numeric(y) & !is.null(attr(y,"val.lab"))){
  
  indices <- as.numeric(colnames(absolutos))
  etiquetasy<-indices%in%attr(y,"val.lab")
  etiquetasy<-sort(c(attr(y,"val.lab"),indices[etiquetasy==F]))
  etiquetasy<-names(etiquetasy[etiquetasy%in%indices])
  etiquetasy<-paste(indices,etiquetasy)   
  etiquetasy<-substr(etiquetasy,1,12)    
  
  }else{
  etiquetasy<-substr(colnames(absolutos),1,12)
}
    
  colnames(absolutos)<-etiquetasy
  
  eliminary<-colSums(absolutos)==0
    
  absolutos<-rbind(absolutos,tapply(w,y,sum,na.rm=T))
  filapct<-addmargins(prop.table(absolutos,1)*100,2)
  absolutos<-cbind(absolutos,Total=rowSums(absolutos))
  colpct<-apply(absolutos,2,function(x)x/x[length(x)])*100
  totpct<-(absolutos/absolutos[nrow(absolutos),ncol(absolutos)])*100
  
  if (clase=="numeric"){
    absolutos<-cbind(as.numeric(rownames(absolutos)),absolutos)
    filapct<-cbind(as.numeric(rownames(filapct)),filapct)
    colpct<-cbind(as.numeric(rownames(colpct)),colpct)
    totpct<-cbind(as.numeric(rownames(totpct)),totpct)
    
      tabla<-NULL
    for (i in 1:nrow(absolutos)){
      tabla<-rbind(tabla,round(absolutos[i,],0),filapct[i,],colpct[i,],totpct[i,])  
    }
  
  
  tabla<-cbind(rep(etiquetas,each=4),rep(c("count","row","col","tot"),nrow(absolutos)),tabla)
    
    if (!missing(niveles)){
      tabla<-tabla[which(tabla[,2]%in%cells & tabla[,3]%in%unique(niveles) | tabla[,2]%in%cells & is.na(tabla[,3])),-c(2,3)]
      etiquetas<-tabla[,1]
      etiquetas[-c(seq(1,length(etiquetas),length(cells)))]<-""
      tabla[,1]<-etiquetas
    }else{
      tabla<-tabla[which(tabla[,2]%in%cells),-c(2,3)]
      etiquetas<-ifelse(duplicated(tabla[,1]),"",tabla[,1])
      tabla[,1]<-etiquetas
    }
  }else{
    absolutos<-cbind(rownames(absolutos),ifelse(absolutos==0,"",round(absolutos,0)))
    filapct<-cbind(rownames(filapct),ifelse(filapct==0,"",filapct))
    colpct<-cbind(rownames(colpct),ifelse(colpct==0,"",colpct))
    totpct<-cbind(rownames(totpct),ifelse(totpct==0,"",totpct))
    tabla<-NULL
    for (i in 1:nrow(absolutos)){
      tabla<-rbind(tabla,absolutos[i,],filapct[i,],colpct[i,],totpct[i,])  
    }
      
    tabla<-cbind(rep(etiquetas,each=4),rep(c("count","row","col","tot"),nrow(absolutos)),tabla)
        
    if (!missing(niveles)){
      tabla<-tabla[which(tabla[,2]%in%cells & tabla[,3]%in%unique(niveles) | tabla[,2]%in%cells & tabla[,3]==""),-c(2,3)]
      etiquetas<-tabla[,1]
      etiquetas[-c(seq(1,length(etiquetas),length(cells)))]<-""
      tabla[,1]<-etiquetas
      
    }else{
      tabla<-tabla[which(tabla[,2]%in%cells),-c(2,3)];tabla
      etiquetas<-ifelse(duplicated(tabla[,1]),"",tabla[,1])
      tabla[,1]<-etiquetas
    }
  }
  
    
    
  if (length(unique(y)==1) & all(unique(y)=="siny")) {
    tabla<-tabla[,-2,drop=F]
  }
  
  
  pct<-cells
  pct<-ifelse(pct!="count",paste(pct,"%",sep=""),pct)
  pct<-paste(pct,collapse=", ")
  
  
  rownames(tabla)<-tabla[,1]
  tabla<-tabla[,-1,drop=F]
  tabla<-apply(tabla,2,function(x)as.numeric(x))
  rownames(tabla)<-etiquetas      
  tabla<-tabla[,!eliminary,drop=F]
  
  if (Totcol==F) tabla<-tabla[,-ncol(tabla),drop=F]
  if (Totrow==F) tabla<-tabla[-(nrow(tabla):((nrow(tabla)-(length(cells)-1)))),,drop=F]

  if (!missing(selectrow)) tabla<-tabla[selectrow,,drop=F]
  if (!missing(selectcol)) tabla<-tabla[,selectcol,drop=F]
  
  if (is.null(attr(y,"var.lab"))) attr(y,"var.lab")<-""
  
  tabla<-structure(tabla,class=c("tabla","matrix","ctables"))
  attr(tabla,"formato")<-cumsum(rep(length(cells),nrow(tabla)/length(cells))[-1])
  attr(tabla,"pct")<-pct
  attr(tabla,"title")<-paste(title,"by",names(as.list(attr(y,"var.lab"))),attr(y,"var.lab"))
  attr(tabla,"ancho")<-ancho
  attr(tabla,"dec")<-dec
  return(tabla)
}

fusion<-function(df1,df2,...){
  df<-merge(x=df1,y=df2,...)  
  
  indices<-names(df)
  copyattr<-function(x){
    v<-df[,x]  
    if(x%in%names(df1)){
    attributes(v)<-attributes(df1[,x])
    return(v)
    }else{
    attributes(v)<-attributes(df2[,x])  
    return(v)
    }
  } 
  df<-data.frame(sapply(indices,copyattr,simplify = F))
  return(df)
}

comul<-function(x,ncp,sufix){

  require(FactoMineR)
  if (class(x) != "list") 
    stop("Las tablas han de ser introducidas como una list()")
  if (length(x) < 2) 
    stop("Para una sola tabla haz un Analisis de Correspondencia Simple")
  if (length(apply(sapply(x, dim), 1, function(k) unique(k))) != 
        2) 
    stop("Las tablas han de tener las misma dimensiones")
  etiquetas <- sapply(x, colnames)
  if (sum(apply(etiquetas, 1, function(k) length(unique(k)))) != 
        ncol(x[[1]])) 
    stop("Los nombres de las columnas no coinciden")
  etiquetas <- sapply(x, rownames)
  if (sum(apply(etiquetas, 1, function(k) length(unique(k)))) != 
        nrow(x[[1]])) 
    stop("Los nombres de las filas no coinciden")
  if (missing(sufix)) {
    sufix <- c(1:length(x))
  }
  else {
    if (length(sufix) != length(x)) 
      stop("Introduzca un numero de sufijos igual al de tabla que contiene la lista")
  }
  
  x<-lapply(x,as.data.frame)
  
  tablam <- x[[1]]
  for (k in 2:length(x)) tablam <- tablam + x[[k]]
  tablam <- tablam/length(x)
  if (missing(ncp)) 
    ncp <- min(dim(tablam))
  
  ncolt<-ncol(tablam)
  nrowt<-nrow(tablam)
    
  for (i in 1:length(x)) {
    colnames(x[[i]])<-paste(colnames(x[[i]]),sufix[i],sep="_")
    tablam <- cbind(tablam,x[[i]])
  }
  
  resultado <- CA(tablam, ncp = ncp, col.sup=c((ncolt+1):(ncol(tablam))),graph = F)
  
  tablam<-tablam[,-c(((ncolt+1):(ncol(tablam)))),drop=F]
  
  x<-lapply(x,function(k){
    colnames(k)<-colnames(tablam)
    return(k)})
  
  
  for (i in 1:length(x)) {
    rownames(x[[i]])<-paste(rownames(x[[i]]),sufix[i],sep="_")
    tablam <- rbind(tablam,x[[i]])
  }
  
  resultado$row.sup <-CA(tablam, ncp = ncp, row.sup = ((nrowt+1):(nrow(tablam))), graph = F)$row.sup
  resultado<-structure(resultado,class=c("list","comul"))
  return(resultado)
}

plot.comul<-function(x,dim=c(1,2),draw=c("col.sup","row.sup"),select){

X <- factor(c(rep("col", nrow(x$col$coord)), 
              rep("row",nrow(x$row$coord)),
              rep("col.sup", nrow(x$col.sup$coord)), 
              rep("row.sup", nrow(x$row.sup$coord))),
              levels = c("col","row", "col.sup", "row.sup"))

df<-data.frame(rbind(x$col$coord,x$row$coord,x$col.sup$coord,x$row.sup$coord),Puntos=X)

limx<-c(min(df[,dim[1]]),max(df[,dim[1]]))
limy<-c(min(df[,dim[2]]),max(df[,dim[2]]))


if (!missing(select)) 
  df <- df[unlist(as.vector((sapply(select, function(x) grep(x,rownames(df), ignore.case = T))))), ]
  
df<-df[df$Puntos%in%draw,c(dim,ncol(df)),drop=F]
  
plot(df[,-ncol(df),drop=F],xlim=limx,ylim=limy,cex=0,cex.axis=0.6,cex.lab=0.6)
text(df[,-ncol(df),drop=F],rownames(df),cex=0.6,col=rainbow(nlevels(df$Puntos),v=0.6)[as.numeric(df$Puntos)])
abline(h = 0, v = 0, lty = 3)

}

print.comul<-function(x,...){
  print.CA(structure(x,class=c("CA","list","comul")),...)
}
summary.comul<-function(x,...){
  summary.CA(structure(x,class=c("CA","list","comul")),...)
}

ponderar<-function (variables, pesos, vp, dif = 1, iter = 100, N) {
    if (!is.list(variables)) 
        stop("Las variables han de ser introducidas como una lista")
    if (!is.list(pesos)) 
        stop("Los pesos han de ser introducidos como una lista")
    if (length(variables) != length(pesos)) 
        stop("El numero de variables es diferente al de pesos")
    if (missing(vp)) 
        vp <- rep(1, length(variables[[1]]))
    
    contador <- 0
    posicion <- 1:length(variables)
    while (length(posicion) != 0) {
        contador <- contador + 1
        for (v in posicion) {
            vari <- variables[[v]]
            muestra <- tapply(vp, vari, sum, na.rm = T)
            muestra<-prop.table(muestra)*100
            pes <- pesos[[v]]
            valores <- sort(unique(vari))
            for (i in 1:length(valores)) {
                vp[vari == valores[i]] <- vp[vari == valores[i]] * (pes[i]/muestra[i])
            }
        }
        final <- lapply(variables, function(k) 
          prop.table(tapply(vp,k, sum, na.rm = T))*100)
        condicion <- NULL
        for (i in 1:length(variables)) {
          condicion <- c(condicion, all(abs(round(final[[i]],1) - pesos[[i]]) <= dif))
        }
        posicion <- which(condicion == F)
        if (contador == iter) 
            break
    }
if (!missing(N)) vp <- vp*(N/sum(vp,na.rm=T))
cat("Total iteraciones alcanzadas: ",contador,"\n\n")
print(lapply(variables,freq,w=vp))
return(vp)
}                        
                        
varlab<-gtools::defmacro(df,var,val,
                         expr={
                           if(!as.character(substitute(var))%in%names(df)) stop("La variable no existe")
                           attr(df$var,"var.lab")<-val
                           names(attr(df$var,"var.lab"))<-as.character(substitute(var))
                         })

vallab<-gtools::defmacro(df,var,val,
                         expr={
                           if(!as.character(substitute(var))%in%names(df)) stop("La variable no existe")
                           if (!is.numeric(val)) stop("El vector ha de ser numerico")
                           if (length(names(val)) == 0) stop("La variable def se ha introducido sin etiquetas")
                           if (length(unique(val)) != length(val)) stop("valores repetidos")
                           attr(df$var,"val.lab")<-sort(val)
                         })

addvallab<-gtools::defmacro(df,var,val,
                            expr={
                              if(!as.character(substitute(var))%in%names(df)) stop("La variable no existe")
                              if (!is.numeric(val)) stop("El vector ha de ser numerico")
                              if (length(names(val)) == 0) stop("La variable def se ha introducido sin etiquetas")
                              if (length(unique(val)) != length(val)) stop("valores repetidos")
                              attr(df$var, "val.lab")<-sort(c(sort(val),attr(df$var,"val.lab"))[!duplicated(c(sort(val),attr(df$var,"val.lab")))])
                              #if (length(a[a%in%bal])>0) names(a)[a[a%in%bal]]<-names(bal)[bal%in%a[a%in%bal]]
                              #attr(df$var,"val.lab")<-sort(c(a,bal[!bal%in%a[a%in%bal]]))
                            })
sumar<-function(...){
  x<-cbind(...)
  suma<-rowSums(cbind(x), na.rm =T)
  esna<-rowSums(!is.na(x))
  suma[esna==0]<-NA
  return(suma)
}

si<-function(test,yes,no){
  atributos<-attributes(no)    
  test[is.na(test)] <- F
  valor<-ifelse(test, yes, no)
  attributes(valor)<-atributos
  return(valor)
}

nvalid<-function(x){
rowSums(!is.na(x))
}

nmis<-function(x){
rowSums(is.na(x))
}


cruce<-function (x, y, w, cells = "count", dec = 1, order, ancho = 12, 
          Totrow = T, Totcol = T, selectrow, selectcol) {
if (is.null(x)) 
  stop("La variable x no existe")
if (is.null(y)) 
  stop("La variable y no existe")
if (missing(w)) 
  w <- rep(1, length(x))
cells <- c("count", "row", "col", "tot")[which(c("count", 
                                                 "row", "col", "tot") %in% unique(cells))]
comb<-expand.grid(sort(unique(x)),sort(unique(y)))
comb$Freq<-apply(comb,1,function(k) sum(w[x==k[1] & y==k[2]],na.rm=T))
absolutos<-tapply(comb$Freq,list(comb$Var1,comb$Var2),sum)
absolutos <- absolutos[!(rowSums(absolutos) == 0), !colSums(absolutos) == 
                         0, drop = F]
if (all(dim(absolutos) == c(0, 0)) == T) 
  stop("Cross is empty")
colpct <- prop.table(absolutos, 2) * 100
filapct <- prop.table(absolutos, 1) * 100
totpct <- prop.table(absolutos) * 100
if (is.numeric(x) & !is.null(attr(x, "val.lab"))) {
  indices <- as.numeric(rownames(absolutos))
  etiquetas <- indices %in% attr(x, "val.lab")
  etiquetas <- sort(c(attr(x, "val.lab"), indices[etiquetas == 
                                                    F]))
  etiquetas <- names(etiquetas[etiquetas %in% indices])
  rownames(absolutos) <- paste(indices, etiquetas)
}
if (is.numeric(y) & !is.null(attr(y, "val.lab"))) {
  indices <- as.numeric(colnames(absolutos))
  etiquetas <- indices %in% attr(y, "val.lab")
  etiquetas <- sort(c(attr(y, "val.lab"), indices[etiquetas == 
                                                    F]))
  etiquetas <- names(etiquetas[etiquetas %in% indices])
  colnames(absolutos) <- paste(indices, etiquetas)
}
if (!missing(order)) {
  if (order %in% c("a", "d")) {
    if (order == "a") {
      orden <- order(rowSums(absolutos))
      absolutos <- absolutos[orden, ]
      colpct <- colpct[orden, ]
      filapct <- filapct[orden, ]
      totpct <- totpct[orden, ]
    }
    else {
      orden <- order(-rowSums(tabla))
      absolutos <- absolutos[orden, ]
      colpct <- colpct[orden, ]
      filapct <- filapct[orden, ]
      totpct <- totpct[orden, ]
    }
  }
}
absolutos <- addmargins(absolutos, FUN = list(Total = sum), 
                        quiet = T)
colpct <- addmargins(colpct, 1, FUN = list(Total = sum), 
                     quiet = T)
filapct <- addmargins(filapct, 2, FUN = list(Total = sum), 
                      quiet = T)
totpct <- addmargins(totpct, FUN = list(Total = sum), quiet = T)
colpct <- cbind(colpct, Total = totpct[, ncol(totpct)])
filapct <- rbind(filapct, Total = totpct[nrow(totpct), ])
absolutos <- round(absolutos, 0)
tabla <- NULL
for (i in 1:nrow(absolutos)) {
  tabla <- rbind(tabla, absolutos[i, ], filapct[i, ], colpct[i, 
                                                             ], totpct[i, ])
}
rownames(tabla) <- paste(rep(rownames(absolutos), each = 4))
pct <- rep(c("count", "row", "col", "tot"), nrow(absolutos))
tabla <- tabla[which(pct %in% cells), , drop = F]
rownames(tabla)[-c(seq(1, length(rownames(tabla)), length(cells)))] <- ""
pct <- cells
pct <- ifelse(pct != "count", paste(pct, "%", sep = ""), 
              pct)
pct <- paste(pct, collapse = ", ")
if (is.null(attr(x, "var.lab"))) 
  attr(x, "var.lab") <- ""
if (is.null(attr(y, "var.lab"))) 
  attr(y, "var.lab") <- ""
if (Totcol == F) 
  tabla <- tabla[, -ncol(tabla), drop = F]
if (Totrow == F) 
  tabla <- tabla[-(nrow(tabla):((nrow(tabla) - (length(cells) - 
                                                  1)))), , drop = F]
if (!missing(selectrow)) 
  tabla <- tabla[selectrow, , drop = F]
if (!missing(selectcol)) 
  tabla <- tabla[, selectcol, drop = F]
if (is.null(attr(x, "var.lab"))) 
  attr(x, "var.lab") <- ""
if (is.null(attr(y, "var.lab"))) 
  attr(y, "var.lab") <- ""
tabla <- structure(tabla, class = c("tabla", "matrix", "cross"))
attr(tabla, "formato") <- cumsum(rep(length(cells), nrow(absolutos))[-1])
attr(tabla, "pct") <- pct
attr(tabla, "resumen") <- paste("Number of Missing Observations:", 
                                round(sum(w, na.rm = T), 0) - absolutos[nrow(absolutos), 
                                                                        ncol(absolutos)])
attr(tabla, "title") <- paste(names(as.list(attr(x, "var.lab"))), 
                              attr(x, "var.lab"), "by", names(as.list(attr(y, "var.lab"))), 
                              attr(y, "var.lab"))
attr(tabla, "ancho") <- ancho
attr(tabla, "dec") <- dec
return(tabla)
}
