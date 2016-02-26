library(ReporteRs)
options('ReporteRs-default-font'='Arial')
#options('ReporteRs-fontsize'=14, 'ReporteRs-default-font'='Arial')
#doc = declareTitlesStyles(doc,stylenames = c("Titre1", "Titre2", "Titre3" , "Titre4", "Titre5", "Titre6", "Titre7", "Titre8", "Titre9","Normal" ) )

echo<-function(txt){
  print(txt)
  doc<<-addParagraph(doc,pot(txt,textProperties(color='black', font.size =9))) 
  doc<<-addParagraph(doc,pot("",textProperties(color='black', font.size =9))) 
}

newpage<-function(...){
  doc<<-addPageBreak(doc,...)
}

tabu<-function(tabla,title,anchos,...){

if (any(c("cross","ctables","cmeans","multiresp")%in%class(tabla))) {  

tablaprint<-tabla
if (any(tabla==0))tablaprint[tabla[,]==0]<-NA
tablaprint<-round(tablaprint,attr(tabla,"dec"))
tablaprint<-cbind(rownames(tablaprint),tablaprint)[,-1,drop=F]
if (any(is.na(tablaprint))) tablaprint[is.na(tablaprint)]<-""

if(!c("cmeans")%in%class(tabla)){
  cells<-attr(tabla,"pct")
  cells<-unlist(strsplit(cells, ", ", fixed = TRUE))
  
  if ("Count"%in%cells) {
    tablaprint[-c(seq(1, length(tablaprint), length(cells))),]<-paste0(tablaprint[-c(seq(1, length(tablaprint), length(cells))),],"%")
    tablaprint[tablaprint=="%"]<-""
  }else{
    tablaprint[tablaprint!=""]<-paste0(tablaprint[tablaprint!=""],"%")
  }
}

}else{
tablaprint<-tabla
}

tb<-FlexTable(tablaprint,add.rownames = T,header.columns = F,
                #header.par.props = parProperties(text.align="center",padding = 0),
                #header.text.props = textProperties( color = "black",font.size = 9), 
                body.text.props = textProperties(font.size = 8))
                
tb[,1]<-parProperties(text.align = 'left', padding=0)
tb[,2:tb$numcol]<-parProperties(text.align='center',padding=0)

tb<-setFlexTableBorders(tb,
                          inner.vertical = borderProperties(color="black", style="solid",width = 3 ),
                          inner.horizontal = borderNone(),
                          outer.vertical = borderProperties( color = "black", style = "solid", width = 3 ),
                          outer.horizontal = borderProperties( color = "black", style = "solid", width = 3 ))

if (missing(title)){
tb<-addHeaderRow( tb, value =attr(tabla,"title"),colspan=tb$numcol,
                  text.properties=textBold(font.size=9),
                  parProperties(text.align = 'center'),
                  cellProperties(border.top.style="none",border.bottom.style="none",border.left.style="none",border.right.style="none"))
}else{
  tb<-addHeaderRow( tb, value =title,colspan=tb$numcol,
                    text.properties=textBold(font.size=9),
                    parProperties(text.align = 'center'),
                    cellProperties(border.top.style="none",border.bottom.style="none",border.left.style="none",border.right.style="none"))
}

tb<-addHeaderRow(tb,value=c("",colnames(tabla)),
                 text.properties=textNormal(font.size=9),
                 par.properties=parProperties(text.align = 'center'),
                 cell.properties = cellProperties(border.width=3))
                                  

if (!is.null(attr(tabla,"resumen")))
  tb<-addFooterRow(tb,value=attr(tabla,"resumen"),colspan=tb$numcol,
                   textProperties(font.size = 6),
                   parProperties(text.align = 'left'),
                   cellProperties(border.top.style="none",border.bottom.style="none",border.left.style="none",border.right.style="none"))


if (!is.null(attr(tabla,"pct")))
 tb<-addFooterRow(tb,value=paste("Cell contents:",attr(tabla,"pct")),colspan=tb$numcol,
                   textProperties(font.size = 6),
                   parProperties(text.align = 'left'),
                   cellProperties(border.top.style="none",border.bottom.style="none",border.left.style="none",border.right.style="none"))

if (!is.null(attr(tabla,"formato")))
  tb[attr(tabla,"formato"),1:tb$numcol,side='bottom']<-borderProperties(style='solid',width=3)                   
if (!missing(anchos)) 
  tb<-setFlexTableWidths( tb, widths=anchos)

doc<<-addFlexTable(doc,tb)
doc<<-addParagraph(doc,pot("",textProperties(color='black', font.size =9))) 
}
