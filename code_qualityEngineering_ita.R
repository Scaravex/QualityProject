#file nella stessa cartella di origine
library(mvtnorm) 
library(MASS)
library(car)#qualche libreria per rendere i grafici migliori

#i dati sono stati precedentemente filtrati con le regole spiegate sul documento.
setwd("C:\\Users\\MSC\\Desktop\\quality engineering\\Project") ##Settare la cartella dove viene inserito il file

Dati<-read.table("Dataset1.csv",header=T,dec=".",sep=",")
head(Dati)

M = colMeans(Dati) #vettore delle medie--> sarà necessario più avanti per standardizzare i dati in ingresso
S = var(Dati)      #matrice di varianza-covarianza
R = cor(Dati)      #Matrice di correlazione
n = dim(Dati)[1]
k = dim(Dati)[2]
alpha = 0.01

shapiro.test(Dati[,1])  #verifico per ogni variabile
#dovrei verificare ogni direzione, es. somma dei primi due meno il terzo
shapiro.test(Dati[,2] + Dati[,1] - Dati[,3]) #non basta --> voglio avere un test per la normalità multivariata

#voglio che sia normale in ogni direzione --> uso MC Shapiro Test, che sfrutta la logica FDR( Folk Discovery Rate)
mcshapiro.test <- function
                  (X, devstmax = 0.01, sim = ceiling(1/(4*devstmax^2)))
    { 
    library(mvnormtest)
    n   <- dim(X)[1]
    p   <- dim(X)[2]
    mu  <- rep(0,p)
    sig <- diag(p)
    W   <- NULL
    for(i in 1:sim) {Xsim <- rmvnorm(n, mu, sig)
                     W<- c(W, mshapiro.test(t(Xsim))$stat) 
                     }# mshapiro.test: calcola statistica min(W.a) x l'input
    Wmin   <- mshapiro.test(t(X))$stat# min(W.a) per il campione dato
    pvalue <- sum(W < Wmin)/sim       # proporzione di min(W.a) e estremi di Wmin osservato
    devst  <- sqrt(pvalue*(1-pvalue)/sim)
    list(Wmin = as.vector(Wmin), pvalue = pvalue, devst = devst, sim = sim)
    }

mcshapiro.test(Dati)  #Normalità anche multivariata, per cui possiamo usare le usuali 
#carte di controllo senza fare assunzioni sulla dimensione del campione di simulazione in ingresso
Dati <- as.matrix(Dati)

x11()
boxplot(Dati)  #qui analizzo invece la variabilità
pairs(Dati, col='orange',pch=19)  #alcune variabili sono correlate, altre di meno. In generale ho un po' di 
#multicollinearità: non fare model selection e riduzione dimensionale, ma sfruttare le Componenti Principali
# V=sqrt(diag(S))*diag(9)  #A=solve(V)%*%(Dati-M)    --> se volessimo standardizzare a mano: ora usiamo solve, per i prossimi arrivi useremo dev.st e media memorizzati nel training set

Dati.std = scale(Dati)
R = var(Dati.std) #che è uguale a R di prima (ho riscalato tutti le variabili per renderle confrontabili)

x11()  #disegniamo la matrice di correlazione, più è chiaro(bianco) e più le variabili sono correlate
image(R)

# PRINCIPAL COMPONENT ANALYSIS
Dati.std <- data.frame(Dati.std) #dati standardizzati: ora analizziamo quante selezionarne
pc.Dati  <- princomp(Dati.std, scores=T)
summary(pc.Dati)
pc.Dati$sd^2 / sum(pc.Dati$sd^2) # proporzione di varianza spiegata     
scores.Dati  <- pc.Dati$scores # scores (dati proiettati)
load.Dati    <- pc.Dati$loadings  # loadings (coefficienti della combinazione lineare delle variabili, definisce ciascuna componente)
# 0.458665+0.2420976+0.1362624

x11()     #mostro ora perchè conviene scegiere le prime tre
layout(matrix(c(2,3,1,3),2,byrow=T))
barplot(pc.Dati$sdev^2, las=2, main='Componenti principali', ylim=c(0,6), ylab='Variances',col='cornflowerblue')
abline(h=c(1,1), lwd=1, col='red')
barplot(sapply(Dati.std,sd)^2, las=2, main='Variabili originarie', ylim=c(0,3), ylab='Variances',col='lightblue')
plot(cumsum(pc.Dati$sdev^2)/sum(pc.Dati$sde^2), type='b', axes=F, xlab='numero di componenti', ylab='contributo alla varianza totale', ylim=c(0,1),lwd=3,col='forestgreen')
abline(h=c(1,0.9), lwd=2, col='blue')
box()
axis(2,at=0:10/10,labels=0:10/10)
axis(1,at=1:ncol(Dati),labels=1:ncol(Dati),las=2)


x.mean = colMeans(scores.Dati)
x.cov  = var(scores.Dati)
p      = 3
qchisq(1 - alpha, p)  #DEVO TROVARE I QUANTILI GIUSTI:CHIQUADRO --> MA IL CAMPIONE è PICCOLO... (SOTTO I 30 USO FISHER)
cfr.fisher <- ((n^2-1)*p/(n*(n-p)))*qf(1-alpha, p, n-p)  #usando la chisquare sottostimerei:uso il metodo esatto di Fisher
#Identifico il tipo di regione di confidenza di interesse: Regione di previsione(regione ellissoidale) 
eigen(x.cov)$vectors# Direzione degli assi principali dell'ellisse:
r <- sqrt(cfr.fisher)# Lunghezza dei semiassi principali dell'ellisse:
raggio = r * sqrt(eigen(x.cov)$values)

E   = eigen(var(scores.Dati))$values  #autovalori
#costruisco una carta Tquadro sulle prime 3 componenti principali
T2A = scores.Dati[,1]^2 /E[1] +scores.Dati[,2]^2 /E[2]+scores.Dati[,3]^2 /E[3]
T2A

x11()
plot(T2A, type="b", pch=16,ylim=c(0,18),lwd=1)
abline(h=cfr.fisher, lty=1, col='red',lwd=3)

cfr.fisher

#Piano CP 1 - CP2
x11() #DISEGNO LE COMPONENTI PRINCIPALI 
v=c(1,2)
plot(scores.Dati[,1],scores.Dati[,2],type="n",xlab="pc1",
     ylab="pc2", asp=1,xlim=c(-8,8),ylim=c(-6,6))
text(scores.Dati[,1],scores.Dati[,2], cex=0.7)
ellipse(x.mean[v], x.cov[v,v], radius=sqrt(cfr.fisher),col = 'cornflowerblue')

#Piano CP 1 - CP3
x11() #DISEGNO LE COMPONENTI PRINCIPALI
v=c(1,3)
plot(scores.Dati[,1],scores.Dati[,3],type="n",xlab="pc1",ylab="pc3", asp=1,xlim=c(-8,8),ylim=c(-5,5))
text(scores.Dati[,1],scores.Dati[,3], cex=0.7)
ellipse(x.mean[v], x.cov[v,v], radius=sqrt(cfr.fisher), col = 'cornflowerblue')


#Piano CP 2 - CP3
x11()
v=c(2,3)
plot(scores.Dati[,2],scores.Dati[,3],type="n",xlab="pc2",ylab="pc3", asp=1,xlim=c(-8,8),ylim=c(-5,5))
text(scores.Dati[,2],scores.Dati[,3], cex=0.7)
ellipse(x.mean[v], x.cov[v,v], radius=sqrt(cfr.fisher), col = 'cornflowerblue')

library(rgl)   #grafico a tre dimensioni delle prime tre componenti principali
open3d()                    # apre un nuovo device
points3d(scores.Dati[,c(1,2,3)], asp=1, size=4, col="tomato")  # disegna i punti
axes3d()  
v = c(1,2,3)
plot3d(ellipse3d(x.cov[v,v], centre=x.mean[v], level= 9.9/10), alpha=0.15, add = TRUE,col="lightblue") # aggiunge l'ellissoide

#andiamo ora a monitorare la variabilità restante nelle altre 6 componenti principali

cfr.fisherB <- ((n^2-1)*6/(n*(n-6)))*qf(1-alpha,6,n-6)
T2B = scores.Dati[,4]^2 /E[4] +scores.Dati[,5]^2 /E[5]+scores.Dati[,6]^2 /E[6]+scores.Dati[,7]^2 /E[7]+scores.Dati[,8]^2 /E[8]+scores.Dati[,9]^2 /E[9]

x11()
plot(T2B, type="b", pch=16, ylim=c(0,30), lwd=1)
abline(h=cfr.fisherB, lty=1, col='red', lwd=3)



#terzo metodo proposto da Alt(1985): lavoro sulle medie di sottocampioni

head(scores.Dati) #terzo metodo che verifichiamo è quello di cercare di raggruppare i dati in triplette e vedere se l'analisi trova meglio fuori controllo sistematici
Datitriple = scores.Dati[1:10,1:3]
Datitriple[1,1:3]  = (scores.Dati[1,1:3]+scores.Dati[2,1:3]+scores.Dati[3,1:3])/3
Datitriple[2,1:3]  = (scores.Dati[4,1:3]+scores.Dati[5,1:3]+scores.Dati[6,1:3])/3
Datitriple[3,1:3]  = (scores.Dati[7,1:3]+scores.Dati[8,1:3]+scores.Dati[9,1:3])/3
Datitriple[4,1:3]  = (scores.Dati[10,1:3]+scores.Dati[11,1:3]+scores.Dati[12,1:3])/3
Datitriple[5,1:3]  = (scores.Dati[13,1:3]+scores.Dati[14,1:3]+scores.Dati[15,1:3])/3
Datitriple[6,1:3]  = (scores.Dati[16,1:3]+scores.Dati[17,1:3]+scores.Dati[18,1:3])/3
Datitriple[7,1:3]  = (scores.Dati[19,1:3]+scores.Dati[20,1:3]+scores.Dati[21,1:3])/3
Datitriple[8,1:3]  = (scores.Dati[22,1:3]+scores.Dati[23,1:3]+scores.Dati[24,1:3])/3
Datitriple[9,1:3]  = (scores.Dati[25,1:3]+scores.Dati[26,1:3]+scores.Dati[27,1:3])/3
Datitriple[10,1:3] = (scores.Dati[28,1:3]+scores.Dati[29,1:3]+scores.Dati[30,1:3])/3

head(Datitriple)
ET = eigen(var(Datitriple))$values  #autovalori
#ora abbiamo, per usare la nmotazione di Alt(1985), n=3, m=10, p=3
T2Triple = Datitriple[,1]^2 /ET[1] +Datitriple[,2]^2 /ET[2]+Datitriple[,3]^2 /ET[3]
#confrontiamo questi valori con il target, dato da:
cfr.fisherTrip<-3*(10-1)*(3-1)/(30-10-3+1)*qf(1-alpha,3,30-10-3+1) #alt divide tra creazione carta e osservazioni future
cfr.fisherTripFuture<-3*(10+1)*(3-1)/(30-10-3+1)*qf(1-alpha,3,30-10-3+1)

x11()
plot(T2Triple,type="b", pch=16,ylim=c(0,20),lwd=1)
abline(h=cfr.fisherTrip, lty=2, col='lightblue')
abline(h=cfr.fisherTripFuture, lty=1, col='red',lwd=3)

#per finire mostriamo una breve analisi di come abbiamo scelto di usare per gli intervalli di confidenza unitari
#le variabili posizione, Area ed Entropia: clusterizziamo le variabili con il metodo di Ward
 
Dati2<-read.table("Dataset2.csv",header=T,dec=".",sep=",")
head(Dati2)

M2 = colMeans(Dati2) #vettore delle medie--> sarà necessario più avanti per standardizzare i dati in ingresso
S2 = var(Dati2)      #matrice di varianza-covarianza

Dati2<-as.matrix(Dati2)
Dati2.std = Dati2
for(i in 1:9)  
  Dati2.std[,i] = (Dati2.std[,i]-M[i])/(sqrt(S[i,i])) #Ho standardizzato i dati
#con media e deviazione standard calcolata prima
L1 = as.matrix(load.Dati[,1])
L2 = as.matrix(load.Dati[,2])
L3 = as.matrix(load.Dati[,3])

scores2.1 = Dati2.std%*%L1
scores2.2 = Dati2.std%*%L2
scores2.3 = Dati2.std%*%L3

T2test = scores2.1^2 /E[1] +scores2.2^2 /E[2]+scores2.3^2 /E[3]
x11()
plot(T2test, type="b", pch=16, ylim=c(0,55), lwd=1)
abline(h=cfr.fisher, lty=1, col='red', lwd=3)


#Piano CP 1 - CP2
x11() #DISEGNO LE COMPONENTI PRINCIPALI 
v = c(1,2)
plot(scores2.1, scores2.2, type="n", xlab="pc1", ylab="pc2", 
     asp=1, xlim=c(-8,8),ylim=c(-6,6))
text(scores2.1, scores2.2, col='red',cex=0.8)
text(scores.Dati[,1],scores.Dati[,2], cex=0.6)
ellipse(x.mean[v], x.cov[v,v], radius=sqrt(cfr.fisher),col = 'cornflowerblue')

#Piano CP 1- CP3
x11() #DISEGNO LE COMPONENTI PRINCIPALI
v = c(1,3)
plot(scores2.1,scores2.3,type="n",xlab="pc1",ylab="pc3", asp=1,xlim=c(-8,8),ylim=c(-8,8))
text(scores2.1,scores2.3, col='red',cex=0.8)
text(scores.Dati[,1],scores.Dati[,3], cex=0.6)
ellipse(x.mean[v], x.cov[v,v], radius=sqrt(cfr.fisher), col = 'cornflowerblue')

#Piano CP 2 - CP3
x11()
v = c(2,3)
plot(scores2.2,scores2.3,type="n",xlab="pc2",ylab="pc3", asp=1,xlim=c(-4,4),ylim=c(-8,8))
text(scores2.2,scores2.3, ,col='red', cex=0.8)
text(scores.Dati[,2],scores.Dati[,3], cex=0.6)
ellipse(x.mean[v], x.cov[v,v], radius=sqrt(cfr.fisher), col = 'cornflowerblue')

library(rgl)   #grafico a tre dimensioni delle prime tre componenti principali
open3d()                    # apre un nuovo device
points3d(scores.Dati[,c(1,2,3)], asp=1, size=4,col="black")  # disegna i punti
points3d(scores2.1,scores2.2,scores2.3, asp=1, size=4,col="red")  # disegna i punti
axes3d()  
v = c(1,2,3)
plot3d(ellipse3d(x.cov[v,v], centre=x.mean[v], level= 9.9/10), alpha=0.15, add = TRUE,col="lightblue") # aggiunge l'ellissoide

L4 = as.matrix(load.Dati[,4])
L5 = as.matrix(load.Dati[,5])
L6 = as.matrix(load.Dati[,6])
L7 = as.matrix(load.Dati[,7])
L8 = as.matrix(load.Dati[,8])
L9 = as.matrix(load.Dati[,9])

scores2.4 = Dati2.std%*%L4
scores2.5 = Dati2.std%*%L5
scores2.6 = Dati2.std%*%L6
scores2.7 = Dati2.std%*%L7
scores2.8 = Dati2.std%*%L8
scores2.9 = Dati2.std%*%L9

T2Btest = scores2.4^2 /E[4] + scores2.5^2 /E[5] + scores2.6^2 /E[6]+
          scores2.7^2 /E[7] + scores2.8^2 /E[8] + scores2.9^2 /E[9]
#i il T2 sui residui è gigante: infatti ci servivano per costruire i T2 sul campione iniziale e essere certi
#di non costruirli con dati che erano già di per sé degli outliers.


Datitriple2      = scores.Dati[1:9,1:3]
Datitriple2[1,1] = (scores2.1[1]+scores2.1[2]+scores2.1[3])/3
Datitriple2[1,2] = (scores2.2[1]+scores2.2[2]+scores2.2[3])/3
Datitriple2[1,3] = (scores2.3[1]+scores2.3[2]+scores2.3[3])/3
Datitriple2[2,1] = (scores2.1[4]+scores2.1[5]+scores2.1[6])/3
Datitriple2[2,2] = (scores2.2[4]+scores2.2[5]+scores2.2[6])/3
Datitriple2[2,3] = (scores2.3[4]+scores2.3[5]+scores2.3[6])/3
Datitriple2[3,1] = (scores2.1[7]+scores2.1[8]+scores2.1[9])/3
Datitriple2[3,2] = (scores2.2[7]+scores2.2[8]+scores2.2[9])/3
Datitriple2[3,3] = (scores2.3[7]+scores2.3[8]+scores2.3[9])/3
Datitriple2[4,1] = (scores2.1[10]+scores2.1[11]+scores2.1[12])/3 
Datitriple2[4,2] = (scores2.2[10]+scores2.2[11]+scores2.2[12])/3
Datitriple2[4,3] = (scores2.3[10]+scores2.3[11]+scores2.3[12])/3
Datitriple2[5,1] = (scores2.1[13]+scores2.1[14]+scores2.1[15])/3
Datitriple2[5,2] = (scores2.2[13]+scores2.2[14]+scores2.2[15])/3
Datitriple2[5,3] = (scores2.3[13]+scores2.3[14]+scores2.3[15])/3
Datitriple2[6,1] = (scores2.1[16]+scores2.1[17]+scores2.1[18])/3
Datitriple2[6,2] = (scores2.2[16]+scores2.2[17]+scores2.2[18])/3
Datitriple2[6,3] = (scores2.3[16]+scores2.3[17]+scores2.3[18])/3
Datitriple2[7,1] = (scores2.1[19]+scores2.1[20]+scores2.1[21])/3
Datitriple2[7,2] = (scores2.2[19]+scores2.2[20]+scores2.2[21])/3
Datitriple2[7,3] = ( scores2.3[19]+scores2.3[20]+scores2.3[21])/3
Datitriple2[8,1] = (scores2.1[22]+scores2.1[23]+scores2.1[24])/3
Datitriple2[8,2] = (scores2.2[22]+scores2.2[23]+scores2.2[24])/3
Datitriple2[8,3] = (scores2.3[22]+scores2.3[23]+scores2.3[24])/3
Datitriple2[9,1] = (scores2.1[25]+scores2.1[26]+scores2.1[27])/3
Datitriple2[9,2] = (scores2.2[25]+scores2.2[26]+scores2.2[27])/3 
Datitriple2[9,3] = (scores2.3[25]+scores2.3[26]+scores2.3[27])/3

  
Datitriple[2,1:3]  = (scores.Dati[4,1:3]+scores.Dati[5,1:3]+scores.Dati[6,1:3])/3
Datitriple[3,1:3]  = (scores.Dati[7,1:3]+scores.Dati[8,1:3]+scores.Dati[9,1:3])/3
Datitriple[4,1:3]  = (scores.Dati[10,1:3]+scores.Dati[11,1:3]+scores.Dati[12,1:3])/3
Datitriple[5,1:3]  = (scores.Dati[13,1:3]+scores.Dati[14,1:3]+scores.Dati[15,1:3])/3
Datitriple[6,1:3]  = (scores.Dati[16,1:3]+scores.Dati[17,1:3]+scores.Dati[18,1:3])/3
Datitriple[7,1:3]  = (scores.Dati[19,1:3]+scores.Dati[20,1:3]+scores.Dati[21,1:3])/3
Datitriple[8,1:3]  = (scores.Dati[22,1:3]+scores.Dati[23,1:3]+scores.Dati[24,1:3])/3
Datitriple[9,1:3]  = (scores.Dati[25,1:3]+scores.Dati[26,1:3]+scores.Dati[27,1:3])/3
Datitriple[10,1:3] = (scores.Dati[28,1:3]+scores.Dati[29,1:3]+scores.Dati[30,1:3])/3
head(Datitriple2)


#ora abbiamo, per usare la notazione di Alt(1985), n=3, m=9, p=3
T2Triple2 = Datitriple2[,1]^2 /ET[1] + Datitriple2[,2]^2 /ET[2]+
            Datitriple2[,3]^2 /ET[3]
#confrontiamo questi valori con il target, dato da:
cfr.fisherTripFuture

x11()
plot(T2Triple2,type="b", pch=16,ylim=c(0,90),lwd=1)
abline(h=cfr.fisherTripFuture, lty=1, col='red',lwd=5)


Mdiff=colMeans(Dati) - colMeans(Dati2)
n1 = 30
n2 = 27
Sp = ((n1-1)*S + (n2-1)*S2)/(n1+n2-2) #media pesata per i gradi di libertà tra le due covarianze

image(cor(Dati))
image(cor(Dati2)) #le due matrici di varianza-covarianza sono inidicativamente simili
mcshapiro.test(Dati )$pvalue #normalità dati verificata
mcshapiro.test(Dati2)$pvalue

ntot = 1/(1/n1+1/n2)
cfr.t <- qt(1-alpha/(2*k),ntot)
IC.BF <- cbind( Mdiff - cfr.t*sqrt(diag(Sp)/ntot) , Mdiff, Mdiff + cfr.t*sqrt(diag(Sp)/ntot) )
IC.BF # percepisco cambiamenti strong nel trattamento 2, leggeri nel trattamento 1
#intervalli di confidenza alla Bonferroni
colnames(IC.BF) <- c("inf","mean","sup")
for (i in 1:9)   
  print(paste('Rifiuto H0 per a',i,': ', !(0>IC.BF[i,1] & 0<IC.BF[i,3]),sep=''))

z = T2test
lambda = 0.4
z[1] = lambda*T2test[1] + (1-lambda)*mean(T2A)
for (i in 2:27) 
  z[i]=lambda*T2test[i]*(1-(1-lambda)^(2*i))+(1-lambda)*z[i-1]



squadro = T2test
media   = mean(scores.Dati[,1])
squadro[1] = var(scores.Dati[,1])
for (i in 2:27) 
  squadro[i] = lambda*(scores2.1[i]-media)^2 + (1-lambda)*squadro[i-1]

squadro
vi = (2-lambda)/lambda
qchisq(1-alpha/2, vi)
qchisq(alpha/2, vi)

UCL = sqrt(var(scores.Dati[,1]))*sqrt(qchisq(1-alpha/2, vi)/vi)
LCL = sqrt(var(scores.Dati[,1]))*sqrt(qchisq(alpha/2  , vi)/vi)
sqrt(squadro)

x11()
plot(sqrt(squadro), pch=19, ylim=c(1.3,3))
abline(h=LCL, lty=2, lwd=4, col='red')
abline(h=UCL, lty=2, lwd=4, col='red')
