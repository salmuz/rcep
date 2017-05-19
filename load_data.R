## Jeu de données : biscuit non cuit 

## Echantillon pour la regression simple
ozone_simple<- read.csv("data/ozone_simple.csv", header = TRUE, sep = ";")
ozone_long<- read.csv("data/ozone_long.csv", header = TRUE, sep = ";")

## Charger l'ehantillon
Xbrut.app <- matrix(scan("data/nirc.asc"), ncol=700, byrow=T)
Ybrut.app <- matrix(scan("data/labc.asc"), ncol=4, byrow=T)
Xbrut.val <- matrix(scan("data/nirp.asc"), ncol=700, byrow=T)
Ybrut.val <- matrix(scan("data/labp.asc"), ncol=4, byrow=T)

#Modele Univarie
Y.unit <- 2

## Echantillon pour trouver les estimateurs du modele
cookie.app <- cbind.data.frame(Ybrut.app[,Y.unit], Xbrut.app)
names(cookie.app) <- c("sucres", paste("X", 1:ncol(Xbrut.app), sep=""))

## Echantillon pour validar el modele 
cookie.val <- cbind.data.frame(Ybrut.val[,Y.unit], Xbrut.val)
names(cookie.val) <- c("sucres", paste("X", 1:ncol(Xbrut.val), sep=""))
