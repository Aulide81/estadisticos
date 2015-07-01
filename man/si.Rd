\name{si}
\alias{si}

\title{
Seleccion condicional
}
\description{
Comprueba si se cumple una condicion y devuelve un valor si se evalua como TRUE y otro si se evalua como FALSE
}
\usage{
si(test,yes,no)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test}{
Prueba logica que se evalua como TRUE o FALSE
}
  \item{yes}{
Valor que se devolvera si test resulta TRUE
}
  \item{no}{
Valor que se devolvera si test resulta FALSE
}
}

\value{
Devuelve un vector de la misma longitud que test
}
\references{
Funcion ifelse
}
\author{
Emilio Arenas
}

\keyword{si}