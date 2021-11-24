import numpy as np

def Dielectric(distance):
    denominator=1-7.7839*np.exp(-0.3153*distance)

    Er=(86.9525/denominator)-8.5225  # This is a constant called Melher-Solajer Dielectric
    return(Er)


def Eelec(distance, ch1, ch2):

    Eelc=332.16*((ch1*ch2)/(Dielectric(distance)*distance)) #here we get the electreostatics energy
    return(Eelc)

def Evdw(distance, sig1, sig2, eps1, eps2):
    
    cof1=((((sig2)**6)*(sig1)**6)/((distance)**12))
    cof2=((((sig2)**3)*(sig1)**3)/((distance)**6))

    Evdw=(4*np.sqrt(eps1*eps2))*(cof1-cof2)# here we get van der waals energy
    return(Evdw)

def Esolvation(sig, ASA):
    R = sig*float(ASA)
    return(R)