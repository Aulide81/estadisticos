\name{cmeans}
\alias{cmeans}

\title{
Tabla de medias
}
\description{
Funcion que proporciona tablas con diferentes estadisticos para variables cuantitativas con un output similar al del SPSS.
}
\usage{
cmeans(x, y, w, dec = 2, stat = "Mean", title = "", ancho = 12, Totcol = T, selectrow, selectcol)
}

\arguments{
  \item{x}{
lista de variables cuantitativas que se desea analizar
}
  \item{y}{
Variables de cabecera por la que se desea obtener los estadisticos de la lista de variables x
}
  \item{w}{
Variable de ponderacion de los casos
}
  \item{dec}{
Numero de decimales que se desean mostrar. Por defecto 2
}
  \item{stat}{
vector string que contenga los estadisticos a mostrar, stat=c("Sum","Mean","Std.Dev","Variance", "Minimum","Maximum", "Valid.N")
}
  \item{title}{
Titulo de la tabla
}
  \item{ancho}{
Ancho de las columnas. Por defecto ancho=12
}
  \item{Totcol}{
Valor logico que indica si se desea mostrar el total columna. Por defecto Totcol=T
}
  \item{selectrow}{
Valor o vector numerico que indica las columnas a retener y por lo tanto que se mostraran
}
  \item{selectcol}{
Valor o vector numerico que indica las filas a retener y por lo tanto que se mostraran
}
}

\value{
Se obtiene una matriz de clase tabla
}

\author{
Emilio Arenas
}

\keyword{cmeans}