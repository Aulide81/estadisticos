\name{means}
\alias{means}

\title{
Descriptivos por grupos
}
\description{
Funcion que permite obtener los estadisticos descriptivos de una variable numerica por una varible categorica
}
\usage{
means(x, y, w, stat = c("Mean", "Std.Dev", "Valid.N"), Totrow = T, dec = 2, selectrow)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Varible numerica
}
  \item{y}{
Variable categorica
}
  \item{w}{
Variable de ponderacion
}
  \item{stat}{
Vector string de los estadisticos a mostrar, stat=c("Sum", "Mean", "Std.Dev", "Variance","Minimum", "Maximum", "Valid.N")
}
  \item{Totrow}{
Valor logico que indica si deseamos mostrar el total fila. Por defecto Totrow=F
}
  \item{dec}{
Valor numerico que indica los decimales que deseamos que muestre
}
  \item{selectrow}{
Valor o vector numerico que indica las columnas a retener y por lo tanto que se mostraran
}
}

\value{
Devuelve un matrix de la clase table
}

\author{
Emilio Arenas
}

\keyword{means}