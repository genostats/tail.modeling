# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 09:14:24 2018

@author: Marion
"""

"""Fonctions utiles pour estimer la p-valeur"""
import numpy as np
import statsmodels.api as sm
from numpy import sklearn.linear_model
import sklearn
from sklearn.linear_model import LinearRegression

"""calcul du nombre de stats simulees superieures a Zobs""" 
def nb_exc(x0,z):
    M = len([x for x in z if x >= x0])
    return M
    
   
"""critere de queue lourde a appliquer a Zsim[[nbsim]]"""
def critereLourd(z,N):
  z = np.sort(z)
  p = [ x / N for x in range(1,N+1)] #proba empiriques p(X<Z)
  d = [z[len([x for x in p if x<=0.97])-1] -  z[len([x for x in p if x<=0.98])-1], #indices commencent à zero en python
               z[len([x for x in p if x<=0.98])-1] - z[len([x for x in p if x<=0.99])-1], 
                 z[len([x for x in p if x<=0.99])-1] - z[len([x for x in p if x<=0.999])-1]] #diff des ecarts des quantiles empiriques
  R = (d[2]-d[3])/(d[1]-d[0])

  if R > 6 : 
    return True  ##queue lourde, on prend le log des statistiques de test
  else:
    return False
   
""" estimateur PWM pour les parametres de la GPD"""
def PWM(z):
  n = len(z)
  p = [ (x - 0.35) / n for x in range(1,n+1)]
  t = np.mean( [(1-x)*y for x,y in zip(p,z)] )
  mu = np.mean(z)
  k = mu/(mu-2*t)-2
  a = 2*mu*t/(mu-2*t)
  return {'a' : a, 'k' : k} #dictionnaire

"""Fonction de repartition de Pareto"""
def FGPD(z,zexc):
  coeff = PWM(zexc)
  a = coeff['a']
  k = coeff['k']
  if k!= 0:
    return [1-(1-k*x/a)**(1/k) for x in z]
  else:
    return [1-np.exp(-x/a) for x in z]

"""calcul de Pgpd"""
def PGPD(x0,zexc,y,seuil):
  M = nb_exc(x0,y)
  if M >= 10:
    return M/len(y)
  else:
    return [len(zexc)/len(y)*(1-x) for x in FGPD(x0-seuil,zexc)]

"""le modele lineaire """
def PML(z,Ntot,N):
  pp = [x / Ntot for x in range(1,N+1)]
  p = pp[::-1] #inverser la liste
  # le modele ~ lineaire 
  # estimation de P( Z > Zobs) par regression lineaire
  if critereLourd(Zsim[[nbsim]],end) == True :
    coeffs = lm(  log(-log(p)) ~ log(log(z)) )$coefficients
    return {'p' : exp(-exp( sum( coeffs*c(1, log(log(Zobs))) ) )), 
                'inter': coeffs[[1]], 
                'pente' : coeffs[[2]]))}
  else:
    coeffs = lm(  log(-log(p)) ~ log(z) )$coefficients
    return{'p' : exp(-exp( sum( coeffs*c(1, log(Zobs)) ) )), 
                'inter': coeffs[[1]], 
                'pente' : coeffs[[2]]))}





























