from Jakobi import Jakobi, Direktna_kinematika
import numpy as np 
import time

def ik(Tc,Tg,theta_trenutni):

    """Tc in Tg predstavljata matriki ki vsebujeta translacijo in rotacijo
    theta_trenutni pa predstavljajo kot rotacije posameznega sklepa. Podani so v obliki numpy matrike"""

    while x >= 10**-2:

        x = Tg-Tc

        dx = 1/10*x

        J = Jakobi(theta_trenutni[0],theta_trenutni[1],theta_trenutni[2],theta_trenutni[3],theta_trenutni[4],theta_trenutni[5],theta_trenutni[6])
        Psevdoinverz_J = np.linalg.pinv(J)

        dq = Psevdoinverz_J*dx

        theta_trenutni += theta_trenutni+dq
        Tc = Direktna_kinematika(theta_trenutni[0],theta_trenutni[1],theta_trenutni[2],theta_trenutni[3],theta_trenutni[4],theta_trenutni[5],theta_trenutni[6])
        x = Tg-Tc

    return theta_trenutni
    