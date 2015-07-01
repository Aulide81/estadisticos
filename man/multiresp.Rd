\name{multiresp}
\alias{multiresp}

\title{
Tabla de Respuesta Multiple
}
\description{
Funcion que obtiene una tabla de respuestas multiple con un output similiar al de SPSS
}
\usage{
multiresp(x, y, w, orden, resp = F, dec = 1, cells = "count", title = "", ancho = 12, Totrow = T, Totcol = T, selectrow, selectcol)
}

\arguments{
  \item{x}{
Lista de variables fila
}
  \item{y}{
Variable cabecera
}
  \item{w}{
Variable de ponderacion
}
  \item{orden}{
Valor string, order="a" o order="d", para ordenar segun columna del Total de moda ascendete o descendente respectivamente
}
  \item{resp}{
Valor logico. Indica si deseamos que los porcentajes se realicen sobre el total citaciones o no. Por defecto resp=F
}
  \item{dec}{
Numero de decimales que se desean mostrar
}
  \item{cells}{
Vector string de estadisticos a mostrar, cells=c("count","row","col","tot
}
  \item{title}{
Titulo de la tabla
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

\keyword{multiresp}