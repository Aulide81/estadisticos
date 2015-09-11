\name{plot.comul}
\alias{plot.comul}

\title{
Grafico del analisis de correspondencias simple de K-tablas
}
\description{
Funcion que permite graficar los resultado del analisis de correspondencias simple de K-tablas
}
\usage{
plot(x,dim=c(1,2),draw=c("col.sup","row.sup"),select)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Objeto de la clase comul
}
  \item{dim}{
Vector numerico que expresa las dimensiones a graficar. Por defecto las dos primeras
}
  \item{draw}{
Vector caracter que indica los cuatro elementos a graficar, draw=c("col","row","col.sup","row.sup"). "col"= Columnas de la matriz de compromiso, "row"= Folumnas de la matriz de compromiso, "col.sup"= Columnas de las tablas suplementarias,
"row.sup"= Filas de las tablas suplementarias. Por defecto dibujara las filas y las columnas suplementarias.
}
  \item{select}{
Vector caracter que permite dibujar solo las filas que se deseen mediante el nombre de estas
}
}

\value{
Devuelve un mapa de imagen
}

\author{
Emilio Arenas
}

\keyword{plot.comul}
