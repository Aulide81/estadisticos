\name{desc}
\alias{desc}

\title{
Descriptivos
}
\description{
Funcion para el analisis univariado
}
\usage{
desc(x, w, ntiles = 1, stat = c("Mean", "Std.Dev", "Minimum", "Maximum", "Valid.N"), dec = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Variable numerica
}
  \item{w}{
Variable de ponderacion
}
  \item{ntiles}{
Valor numeric que indica los percentiles que se desean extraer
}
  \item{stat}{
Vector string que indica los estadisticos que se mostraran, stat=c(
"Mean","Median","Mode","Variance","Std.Dev","S.E.Mean","Skewness","Kurtosis","Minimum","Maximum","Range","Sum","Valid.N") 
}
  \item{dec}{
Valor numerico que indica los decimales que se desean mostrar
}
}

\value{
Devuelve una matrix de la clase tabla
}

\author{
Emilio Arenas
}

\keyword{desc}