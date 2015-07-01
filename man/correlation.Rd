\name{correlation}
\alias{correlation}

\title{
Matriz de correlaciones
}
\description{
Esta funcion calcula la matriz de correlaciones
}
\usage{
correlation(x, y, w, dec = 3)
}

\arguments{
  \item{x}{
Objeto del tipo vector,matrix o data.frame que contenga variables numericas
}
  \item{y}{
En caso de que x sea un vector numerico, y debe ser otro vector numerico. En casos de que x sea una matrix o data.frame con ncol>1, este argumento debe ser omitido
}
  \item{w}{
Variable de ponderacion de los casos
}
  \item{dec}{
Numero de decimales que se desean mostrar
}
}


\author{
Emilio Arenas
}

\keyword{correlation}