\name{compfactor}
\alias{compfactor}

\title{
Transformar variable numerica categorica en factor
}
\description{
Funcion que transforma una variable numerica en factor, aprobechando las etiquetas de valor contenidas en attr(x,"val.lab")
}
\usage{
compfactor(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Variable numerica que disponga de etiquetas en attr(x,"val.lab")
}
}
\value{
Devuelve una variable factor, cuyos niveles son los attr(x,"val.lab")
}

\author{
Emilio Arenas
}

\keyword{compfactor}