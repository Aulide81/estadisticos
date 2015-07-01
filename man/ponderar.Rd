\name{ponderar}
\alias{ponderar}

\title{
Ponderar Casos
}
\description{
Funcion que permite obtener los coeficientes de ponderacion para cada caso. 
}
\usage{
ponderar(variables, pesos, dif = 1, iter = 100, N)
}

\arguments{
  \item{variables}{
Lista de variables a ponderar
}
  \item{pesos}{
Lista con los pesos de las variables a ponderar
}
  \item{dif}{
Valor numerico que permite que el resultado difiera como maximo (en valor absoluto) al objetivo señalado en pesos. Por defecto dif=1
}
  \item{iter}{
Valor numerico. Iteraciones que la funcion realizara como máximo. Por defecto iter=100
}
  \item{N}{
Valor numerico al que se igualara la suma de los coeficientes de ponderacion
}
}
\value{
Se obtiene una variable numerica
}

\author{
Emilio Arenas
}

\keyword{ponderar}