########################################
# This is defines the form of the model 
# the stress, response functions etc. 

########################################
# The following defines how the parameters are vectorized
#   The parameters of the model corresponds to the 
# . exponents of the invariants for each term

# How the index are determined is based on the maximum poly
# . -nomial used, which is the top_degree


#### This get the vector index from the exponents
def vectorize(top_degree, i, j, k):
    n = (top_degree//2 + 1) 
    # Note that this is divided by 2 because only 
    # . even exponent k is allowed
    n2 = (top_degree + 1) * (n)
    return n2*i + n*j + k//2

#### This recovers the exponenets from the vector index
def find_index(top_degree, y):
    n = (top_degree//2 + 1)
    n2 = (top_degree + 1) * (n)
    return y//n2, (y%n2)//n, 2*(y%n)




########################################
# The following gives the strain energy density given a set of 
# . parameters and vectorized index (equivalent of the exponenets)

# g1 is gamma1 the first Hencky strain invariant corresponding 
# . to the stretch in the preferred direction
# g2 is the second invariant
# g1 is gamma1 
# g2 is gamma2
# g3 is gamma3

#### Gives the strain energy contribution for this terms
def Psi(g1, g2, g3, c, y): 
    # c is the model parameter values
    # y is the exponent index
    i,j,k = find_index(poly_order, y) # recovers the exponents
    return c*(g1**i)*(g2**j)*(g3**k) 

#### Gives the total strain energy contribution given a set of terms
def strainenergy(g1, g2, g3, const, degrees):
    
    # const is a set of model parameter values
    # y is a set of corresponding parameter index
    
    result = 0.0
    for y in range(const.size):
        result += Psi(g1, g2, g3, c, y)
    return result



########################################
# The following gives defines the strain energy function derivatives
# . Note here that dWX is equivalent of dW/dX, dWXY is equivalent of 
# . d^2W/dX/dY etc.

# y is the parameters index of the term

#### The following are the first derivatives other known as the 
# . response functions
def dW1(g1, g2, g3, y):
    i,j,k = find_index(poly_order, y)
    if i > 0:
        return (i*g1**(i-1)) * (g2**j) * (g3**k)
    return 0

def dW2(g1, g2, g3, degrees):
    i,j,k = find_index(poly_order, degrees)
    if j > 0:
        return (g1**i) * (j*g2**(j-1)) * (g3**k)
    return 0
    
def dW3(g1, g2, g3, degrees):
    i,j,k = find_index(poly_order, degrees)
    if k > 0:
        return (g1**i) * (g2**j) * (k*g3**(k-1))
    return 0


#### The following are the second derivatives 
def dW11(g1, g2, g3, degrees):
    i,j,k = find_index(poly_order, degrees)
    if i > 1:
        return (i * (i-1) * g1**(i-2)) * (g2**j) * (g3**k)
    return 0
    
def dW22(g1, g2, g3, degrees):
    i,j,k = find_index(poly_order, degrees)
    if j > 1:
        return (g1**i) * (j * (j-1) * g2**(j-2)) * (g3**k)
    return 0
    
def dW33(g1, g2, g3, degrees):
    i,j,k = find_index(poly_order, degrees)
    if k > 1:
        return (g1**i)*(g2**j)*(k*(k-1)*g3**(k-2))
    return 0

    
#### The following are the cross derivatives 
def dW12(g1, g2, g3, degrees):
    i,j,k = find_index(poly_order, degrees)
    if i > 0 and j > 0:
        return (i*g1**(i-1))*(j*g2**(j-1))*(g3**k)
    return 0
    
def dW13(g1, g2, g3, degrees):
    i,j,k = find_index(poly_order, degrees)
    if i > 0 and k > 0:
        return (i*g1**(i-1))*(g2**j)*(k*g3**(k-1))
    return 0
    
def dW23(g1, g2, g3, degrees):
    i,j,k = find_index(poly_order, degrees)
    if j > 0 and k > 0:
        return (g1**i)*(j*g2**(j-1))*(k*g3**(k-1))
    return 0




########################################
# The following gives defines the response functions which is a 
# . sum of each terms of the model with respect to the first derivatives

# const is a set of model parameter values
# degrees is a set of corresponding parameter index

def responseW1(g1, g2, g3, const, degrees):
    result = 0.0
    for y in range(const.size):
        result += const[y]*dW1(g1,g2,g3,degrees[y])
    return result


def responseW2(g1, g2, g3, const, degrees):
    result = 0.0
    for y in range(const.size):
        result += const[y]*dW2(g1,g2,g3,degrees[y])
    return result


def responseW3(g1, g2, g3, const, degrees):
    result = 0.0
    for y in range(const.size):
        result += const[y]*dW3(g1,g2,g3,degrees[y])
    return result
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

