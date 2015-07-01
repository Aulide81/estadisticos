\name{ctables}
\alias{ctables}

\title{
Tabla de Perfiles / Tabla de variables dicotomicas
}
\description{
Funcion que proporciona una tabla de perfiles o una tabla de medias con output similar al de SPSS
}
\usage{
ctables(x, y, w, niveles, orden, cells = "count", title = "", ancho = 12, dec = 1, Totrow = T, Totcol = T, selectrow, selectcol)
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
  \item{niveles}{
Valor numerico que indica el nivel de las variables que seran mostrados en las tablas.Por omision se mostraran todos los niveles creando una tabla de perfiles.
}
  \item{orden}{
Valor string, orden="a" o orden="d", para ordenar segun columna del Total de moda ascendete o descendente respectivamente
}
  \item{cells}{
Vector string de estadisticos a mostrar, cells=c("count","row","col","tot")
}
  \item{title}{
Titulo de la tabla
}
  \item{ancho}{
Valor numerico que indica el ancho de las columnas, por defecto ancho=12 
}
  \item{dec}{
Numero de decimales que se desean mostrar
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

\keyword{ctables}