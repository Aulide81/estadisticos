\name{cross}
\alias{cross}

\title{
CrossTable
}
\description{
Funcion que implementa una tabulacion cruzada con un output similar al del SPSS
}
\usage{
cross(x, y, w, cells = "count", dec = 1, order, ancho = 12, Totrow = T, Totcol = T, selectrow, selectcol)
}

\arguments{
  \item{x}{
Variable fila
}
  \item{y}{
Variable columna
}
  \item{w}{
Variable de ponderacion de los casos
}
  \item{cells}{
Vector string de estadisticos a mostrar, cells=c("count","row","col","tot")
}
  \item{dec}{
Numero de decimales que se desean mostrar
}
  \item{order}{
Valor string, orden="a" o orden="d", para ordenar segun columna del Total de moda ascendete o descendente respectivamente
}
  \item{ancho}{
Valor numerico que indica el ancho de las columnas, por defecto ancho=12 
}
  \item{Totrow}{
Valor logico que indica si deseamos mostrar el Total fila, por defecto Totrow=T
}
  \item{Totcol}{
Valor logico que indica si deseamos mostrar el Total columna, por defecto Totcol=T
}
  \item{selectrow}{
Valor o vector numerico que indica las columnas a retener y por lo tanto que se mostraran
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

\keyword{cross}
