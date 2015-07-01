\name{fusion}
\alias{fusion}

\title{
Fusionar variables
}
\description{
Funcion que permite fusionar una o mas variable a otro grupo de variables, como pueda ser un data.frame
}
\usage{
fusion(df1, df2, ...)
}

\arguments{
  \item{df1}{
Data.frame anfrition. Sus definiciones prevalecen sobre df2
}
  \item{df2}{
Data.frame invitada
}
  \item{\dots}{
argumentos basado en la funcion merge( )
}
}

\value{
Devuelve un data.frame
}

\author{
Emilio Arenas
}

\keyword{fusion}