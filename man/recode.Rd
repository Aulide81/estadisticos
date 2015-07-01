\name{recode}
\alias{recode}

\title{
Recodificar Variables
}
\description{
Variable que permite recodificar variables numerica/character/factor a otra variable numerica, que aprovechara la definicion de la variable original
}
\usage{
recode(x,...)
}

\arguments{
  \item{x}{
Variable a recodificar
}
  \item{...}{
Las agrupaciones de la variable x, que formaran los grupos de la nueva variables
}
}


\author{
Emilio Arenas
}

\examples{
## Not run: 
#var<-c(rep(1,5),rep(2,3),rep(9,6))#creamos una variable
#attr(var,"var.lab")<-c(var="VARIABLE ORIGINAL")#definicion de la variable

##Definimos los niveles de la variable
#attr(var,"val.lab")<-c("Cod.1"=1,"Cod.2"=2,"Cod.3"=3)

#varr<-recode(var,c(1,2,v=1),c(9,v=2))
            
#attr(varr,"var.lab")<-c(varr="VARIABLE RECODIFICADA")#definicion de la variable
##Definimos los niveles de la variable
#attr(varr,"val.lab")<-c("Cod.1y2"=1,"Cod.9"=2)
#cross(var,varr,cells=c("count","col","tot","row"))
## End(Not run)
}

\keyword{recode}