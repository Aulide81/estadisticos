\name{export.sav}
\alias{export.sav}

\title{
Exportar de Rdata a SAV
}
\description{
Funcion que permite pasar archivos Rdata de variables numericas definidas a SAV
}
\usage{
export.sav(df, datafile, codefile, sep = "\t", dec = ",", na = "", quote = T, drop.factor = F)
}

\arguments{
  \item{df}{
Data.frame con variables numericas definidas
}
  \item{datafile}{
Nombre y ubicacion del archivo de datos
}
  \item{codefile}{
Nombre y ubicacion del archivo de definicion
}
  \item{sep}{
Delimitacion de las variables del archivo datafile
}
  \item{dec}{
Representacion del signo decimal
}
  \item{na}{
Codificacion de los valora ausentes (missings)
}
  \item{quote}{
Valor logico.Los factores o variables string se grabaran entrecomillados? Por defecto quote=T
}
  \item{drop.factor}{
Valor logico. Desea no exportar las variables string o factor? Por defecto drop.factor=F
}
}

\author{
Emilio Arenas
}

\keyword{export.sav}