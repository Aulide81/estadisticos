\name{k.means}
\alias{k.means}

\title{
K-Means
}
\description{
K -means es un metodo de agrupamiento, que tiene como objetivo la particion de un conjunto n observaciones en k grupos en el que cada observacion pertenece al grupo mas cercano a la media.
}
\usage{
k.means(x, centers, w, iter = 10, initial, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data.frame solo con valores numericos
}
  \item{centers}{
Valor numerico que indica el numerode grupos que deseamos obtener
}
  \item{w}{
Variable de ponderacion
}
  \item{iter}{
Numero de iteraciones. Por defecto 10
}
  \item{initial}{
Matrix numerica con los centros iniciales. Por omision, el algoritmo elegira uno al azar
}
  \item{seed}{
Valor numerico para fijar la semilla de probabilidad
}
}

\value{
\item{Grupos}{Vector numeric que contiene la clasificacion para cada observacion}
\item{Anova}{Tabla de analisis de Varianza}
}

\author{
Emilio Arenas
}

\keyword{k.means}