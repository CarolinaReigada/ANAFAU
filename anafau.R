
# abundancia: total de individuos coletados + CI a 5% probabilidade

# muito abundante: no. individuos > do que o limite superior do CI da media
# abundante: no. individuos situados entre o limite superior do CI
# comum : no individuos entre os limites inferior e superior do CI da media comum
# disperso: no. individuos entre o limite inferior do CI da media
# raro: no. individuos < do que o limite inferior do CI da media


# frequencia relativa (%): (no. individuos /total de individuos coletados)*100
# indice da frequencia é estimado pela distribuição dos valores em histograma
# e comparados com base no CI, com os resultados classificados em:
# pouco frequente: freq. <= CI inferior,
# frequente: IC inf <= freq < media
# muito frequente:  media <= freq < CI superior
# super frequente: freq >=IC superior

# Dominancia:
# superdominante e dominante; indice maior que o limite da dominancia
# não dominante: freq. menor

# Constancia: (no. coletas contendo a sp/no.total de spp coletadas)*100
# são classificados em:
# acidental- Z: (<25% das amostras);
# adicional -Y: (25%<=x<50% das amostras);
#constante -W: (>=50% das amostras)


####################### ANALISE FAUNISTICA: ANAFAU  ##########

n_spp=15
n_coleta=12


d=read.table("dados.txt",header=T)
d
attach(d)
summary(d)

N=sum(Total_indiv)
N


### estimativa da media e desvio padrao amostra e desvio padrao media dados total

#media
media=mean(Total_indiv)
media

#desvio padrao amostral
desv_pad=sd(Total_indiv)
desv_pad

#desvio padrao media
dpm=desv_pad/(sqrt(n_spp))
dpm


###### ANALISE DE RESIDUO PARA TIRAR AS ESPECIES DISCREPANTES ####

e=c()
newd=d

#PARA CADA ESPECIE:

for(i in 1:n_spp){
  
  e[i]= (newd[i,1]- media)/desv_pad
  newd[i,1]=newd[i,1]
  
   if(-2<=e[i]&e[i]>=2) {
    newd[i,1:2] = NA
  } 
 
}

e
newd
newd=na.omit(newd)
newd
colnames(newd)=c("T_ind","T_col")
attach(newd)

#### calculo de media desv pad apos remocao dos "super"

Nnew=sum(T_ind)
Nnew
spp_new=length(T_ind)

#media
media=mean(T_ind)
media

#desvio padrao amostral
desv_pad=sd(T_ind)
desv_pad

#desvio padrao media
dpm=desv_pad/(sqrt(spp_new)) 
dpm


## calculo Intervalo de confianca (IC)

LS5= ceiling(media+2*(dpm))
LI5= ceiling(media-2*(dpm))
LS1= ceiling(media+2.58*(dpm))
LI1= ceiling(media-2.58*(dpm))


resultado_abundancia_5=c(LI5,media,LS5)
resultado_abundancia_5

resultado_abundancia_1=c(LI1,media,LS1)
resultado_abundancia_1


########## 1. ABUNDANCIA #############

attach(d)
abund=c()

for(i in 1:n_spp){
  
  
  if(Total_indiv[i]<LI1) {
    abund[i] = "r"
  }
  
  if(LI1<=Total_indiv[i]&Total_indiv[i]<LI5 ) {
    abund[i] = "d"
  }
  
  if(LI5<=Total_indiv[i]&Total_indiv[i]<LS5 ) {
    abund[i] = "c"
  }
  
  if(LS5<=Total_indiv[i]&Total_indiv[i]<LS1 ) {
    abund[i] = "a"
  }
  
  if(Total_indiv[i]>=LS1 ) {
    abund[i] = "ma"
  }
    
    if(-2<=e[i]&e[i]>=2) {
    abund[i] = "sa"
    } 
   
}

abund

################# 2. Frequencia ####################


##Classificacao frequencia 

freq_d=(Total_indiv/sum(Total_indiv))*100 # dados d
freq_d
freq=(T_ind/(sum(T_ind)))*100 # newd com remocao de dados discrepantes
freq

#media
media=mean(freq)
media

#desvio padrao amostral
desv_pad=sd(freq)
desv_pad

#desvio padrao media
dpm=desv_pad/(sqrt(spp_new))
dpm

## calculo Intervalo de confianca (IC)

LS5= ceiling(media+2*(dpm))
LI5= ceiling(media-2*(dpm))

resultado_freq=c(LI5,media,LS5)
resultado_freq


r_freq=c()

for(i in 1:n_spp){
 
#  if(freq_d[i]<=LI5) {
#    r_freq[i] = "pf"
#  }
  
#  if(freq_d[i]>LI5 & freq[i]<LS5 ) {
#   r_freq[i] = "f"
#  }
  
#  if(freq_d[i]>=LS5 ) {
#   r_freq[i] = "mf"
#  }
  
#  if(-2<=e[i]&e[i]>=2) {
#    r_freq[i] = "sf"
 # } 
  
  if(abund[i]=="ma") {
    r_freq[i] = "mf"
  }
  
  if(abund[i]=="a" ) {
    r_freq[i] = "mf"
  }
  
  if(abund[i]=="c") {
    r_freq[i] = "f"
  }
  
  if(abund[i]=="d") {
    r_freq[i] = "pf"
  }
  
  if(abund[i]=="r" ) {
    r_freq[i] = "pf"
  } 
  
  if(abund[i]=="sa" ) {
    r_freq[i] = "sf"
  } 
  
}

r_freq

################## 3. constancia ##################

const=(Total_coletas/n_coleta)*100
const

# classificação:
const

it=length(const)
r_const=rep("Y",it) #acessorias

for(i in 1:it){
  
  if(const[i]<25) {
    r_const[i] = "Z"    #acidentais
  }
  
  
  if(const[i]>=50 ) {
    r_const[i] = "W"  #constante
    }
  
}

r_const

####################### 4. dominancia  #########################

#### segundo kato et al
attach(d)
K=c()
r_dom=rep("ND",n_spp) #nao dominante

for(i in 1:n_spp){
  
  K=Total_indiv[i]
 
  n1=2   #2*(K+1), sendo K=0
  n2= 2*(N+1)    # 2*(N-K+1), sendo K=0 condicao de Kato et al
  F0=qf(0.95,n1,n2)
  #F0 = valor obtido atraves da tabela de distribuicao
  # F ao nivel de 5% de probabilidade, com valores des graus
  # de liberdade iguais a n1 e n2
  
  # limite superior
  LS=((n1*F0)/(n2+(n1*F0)))*100
  
  
  ni1=2*(N-K+1)
  ni2=2*(K+1)
  Fi0=qf(0.95,ni1,ni2)

 # limite inferior
 LI=(1-((ni1*Fi0)/(ni2+ni1*Fi0)))*100
  
  if(LI>LS) {
    r_dom[i] = "D"    #dominante
  }
  
 if(-2<=e[i]&e[i]>=2) {
   r_dom[i] = "SD"
 } 
 
}

r_dom


####################### TABELA FINAL #########################

especie=seq(1:n_spp)

resultados=cbind(especie,d,r_dom,abund,r_freq,r_const)
resultados

## predominantes e indicadores ecológicos:

for(i in 1:n_spp){
 
  if(r_dom[i]=="ND"|r_dom[i]=="D"){
    if(abund[i]=="sa"|abund[i]=="ma"|abund[i]=="a"){
      if(r_freq[i]=="mf"|r_freq[i]=="f"){
        if(r_const[i]=="W"){
          print(paste("especie predominante:", especie[i]))
          
        }
      }
    }
  }
  if(i==n_spp){
    print(paste("Número total de individuos:",N))
    print(paste("Número de espécies:",n_spp))
    print(paste("Número total de coletas:",n_coleta))
   
  } 
  
  
}

