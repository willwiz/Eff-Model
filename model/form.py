########################################
# This is defines the form of the model 
# the stress, response functions etc. 
########################################

# The parameters of the models are listed as followed
# . c0
# . c11 for y1^2
# . c22 for y2^2
# . c33 for y3^2
# . c12 for y1-y2
# . c13 for y1-y3^2
# . c23 for y2-y3^2

# The input of the model are:
# . y1 for gamma_1 which is the log strain along the fiber direction
# . y2 for gamma_2 which is the log strain orthogonal to the fiber direction
# . y3 for gamma_3 which is the shear angle


########################################
# Packages that are needed
########################################
import numpy as np
from numpy.linalg import inv

########################################
# Define the Kinematics
########################################


########################################
# Define preferred orientation

def tensor(F):
    return np.array([[F[0], F[1]], [F[2], F[3]]])

def cross_2D(M):
    return np.array([-M[1], M[0]])

def calc_m(F,M):
    tF = tensor(F)
    vec = tF.dot(M)
    m = vec/np.sqrt(vec.dot(vec))
    return m

def calc_s(F,M):
    m = calc_m(F,M)
    return cross_2D(m)


########################################
# Define invariants

def lambda_M(F,M):
    tF = tensor(F)
    m = calc_m(F,M)
    return m.dot(tF).dot(M)

def lambda_S(F,M):
    tF = tensor(F)
    S = cross_2D(M)
    s = calc_s(F,M)
    return s.dot(tF).dot(S)

def phi(F,M):
    tF = tensor(F)
    S = cross_2D(M)
    m = calc_m(F,M)
    return m.dot(tF).dot(S) / lambda_M(F,M)


def gamma_1(F,M):
    return np.log(lambda_M(F,M))

def gamma_2(F,M):
    return np.log(lambda_S(F,M))

def gamma_3(F,M):
    return phi(F,M)

########################################
# Define the response functions

def W1_fromdat(T,F,M):
    m = calc_m(F,M)
    return m.dot(T).dot(m)

def W2_fromdat(T,F,M):
    s = calc_s(F,M)
    return s.dot(T).dot(s)

def W3_fromdat(T,F,M):
    s = calc_s(F,M)
    return m.dot(T).dot(s) * lambda_M(F,M) / lambda_S(F,M)





########################################
# first we define the strain energy function

def model_Q(c, F, M):
    y1 = gamma_1(F,M)
    y2 = gamma_2(F,M)
    y3 = gamma_3(F,M)
    return np.array([c[1]*y1*y1, c[2]*y2*y2, c[3]*y3*y3, c[4]*y1*y2, c[5]*y1*y3, c[6]*y2*y3])

def model_Q_s(c, F, M, ymax):
    y1 = gamma_1(F,M)
    y2 = gamma_2(F,M)
    y3 = gamma_3(F,M)
    return np.array([ c[1]*(y1*y1 - ymax[0]*ymax[0]), c[2]*(y2*y2 - ymax[1]*ymax[1]), 
                 c[3]*(y3*y3 - ymax[2]*ymax[2]), c[4]*(y1*y2 - ymax[0]*ymax[1]), 
                 c[5]*(y1*y3 - ymax[0]*ymax[2]), c[6]*(y2*y3 - ymax[1]*ymax[2]) ])


    
def model_Q1(c, F, M):
    y1 = gamma_1(F,M)
    y2 = gamma_2(F,M)
    y3 = gamma_3(F,M)
    return np.array([ 2*c[1]*y1, 0, 0, c[4]*y2, c[5]*y3, 0 ])    
    
def model_Q2(c, F, M):
    y1 = gamma_1(F,M)
    y2 = gamma_2(F,M)
    y3 = gamma_3(F,M)
    return np.array([ 0, 2*c[2]*y2, 0, c[4]*y1, 0, c[6]*y3 ])
    
def model_Q3(c, F, M):
    y1 = gamma_1(F,M)
    y2 = gamma_2(F,M)
    y3 = gamma_3(F,M)
    return np.array([ 0, 0, 2*c[3]*y3, 0, c[5]*y1, c[6]*y2 ])




def model_Psi(c, F, M):
    Q = model_Q(c, tF)
    ans = c[0] * ( (np.exp(Q)).sum() - 1.0 )
#    ans = c[0] * ( np.exp(Q.sum()) - 1.0 )
    return ans




def model_W1(c, F, M):
    Q = model_Q(c, F, M)
    Q1 = model_Q1(c, F, M)
    ans = c[0] * (Q1 * np.exp(Q)).sum()
#    ans = c[0] * (Q1.sum() * np.exp(Q.sum()))
    return ans

def model_W2(c, F, M):
    Q = model_Q(c, F, M)
    Q2 = model_Q2(c, F, M)
    ans = c[0] * (Q2 * np.exp(Q)).sum()
#    ans = c[0] * (Q2.sum() * np.exp(Q.sum()))
    return ans

def model_W3(c, F, M):
    Q = model_Q(c, F, M)
    Q3 = model_Q3(c, F, M)
    ans = c[0] * (Q3 * np.exp(Q)).sum()
#    ans = c[0] * (Q3.sum() * np.exp(Q.sum()))
    return ans



def model_W1_s(c, F, M, ymax):
    Q = model_Q_s(c, F, M, ymax)
    Q1 = model_Q1(c, F, M)
    ans = c[0] * (Q1 * np.exp(Q)).sum()
#    ans = c[0] * (Q1.sum() * np.exp(Q.sum()))
    return ans

def model_W2_s(c, F, M, ymax):
    Q = model_Q_s(c, F, M, ymax)
    Q2 = model_Q2(c, F, M)
    ans = c[0] * (Q2 * np.exp(Q)).sum()
#    ans = c[0] * (Q2.sum() * np.exp(Q.sum()))
    return ans

def model_W3_s(c, F, M, ymax):
    Q = model_Q_s(c, F, M, ymax)
    Q3 = model_Q3(c, F, M)
    ans = c[0] * (Q3 * np.exp(Q)).sum()
#    ans = c[0] * (Q3.sum() * np.exp(Q.sum()))
    return ans

    
########################################
# first we define the strain energy function
    
def T_basis(F, M):
    m = calc_m(F,M)
    s = calc_s(F,M)
    [np.array([m[0]*m[0], m[0]*m[1], m[1]*m[0], m[1]*m[1]]),
     np.array([s[0]*s[0], m[0]*m[1], m[1]*m[0], m[1]*m[1]]),
     np.array([ 2*m[0]*s[0], m[0]*s[1] + s[0]*m[1], m[1]*s[0] + s[1]*m[0], 2*m[1]*s[1] ]) * ( lambda_S(F,M) / lambda_M(F,M) )
    ]
    
def S_basis(F, M):
    
    m = inv(tensor(F)).dot(calc_m(F,M))
    s = inv(tensor(F)).dot(calc_s(F,M))
    [np.array([m[0]*m[0], m[0]*m[1], m[1]*m[0], m[1]*m[1]]),
     np.array([s[0]*s[0], m[0]*m[1], m[1]*m[0], m[1]*m[1]]),
     np.array([ 2*m[0]*s[0], m[0]*s[1] + s[0]*m[1], m[1]*s[0] + s[1]*m[0], 2*m[1]*s[1] ]) * ( lambda_S(F,M) / lambda_M(F,M) )
    ]
    
    
    
    

