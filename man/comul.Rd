\name{comul}
\alias{comul}

\title{
Analisis de correspondencias simple de K-tablas
}
\description{
Funcion que permite realizar un analisis de imagen con diferentes numeros de tablas y determinar como se han desplazado los diferentes items analizados
}
\usage{
comul(x, ncp, sufix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Lista que contenga las diferentes tablas para analizar. Estas tablas han de ser de identicas dimensiones y han de tener identicos rownames( ) e identicos colnames( )
}
  \item{ncp}{
Numero de ejes a conservar. Por omision el maximo posible
}
  \item{sufix}{
Vector de misma longitud que la lista x, que contenga los sufijos para nombrar los rownames( ) y colnames( )
}
}

\value{
Devuelve un objeto del tipo CA
}
\references{
Package FactoMineR
}
\author{
Emilio Arenas
}


\keyword{comul}