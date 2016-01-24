\name{freq}
\alias{freq}
\alias{freq.numeric}
\alias{freq.character}
\alias{freq.factor}
\alias{freq.data.frame}

\title{
Tabla de Frecuencias
}
\description{
Funcion que implementa una tabla de frecuencias con un output similar al del SPSS
}
\usage{
freq(x,w,order,dec=1,selectcol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Variable o data.frame
}

  \item{w}{
Variable de ponderacion de los casos
}
  
  \item{dec}{
Numero de decimales que se desean mostrar
}
  \item{order}{
Valor logico, order="a" o order="d", para ordenar segun columna del Total de moda ascendete o descendente respectivamente
}
  \item{selectcol}{
Valor o vector numerico que indica las filas a retener y por lo tanto que se mostraran
}
}

\value{
Devuelve una matrix de la clase tabla
}

\author{
Emilio Arenas
}

\keyword{freq}
