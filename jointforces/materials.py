def linear(stiffness):
    """
    Defines the stiffness for a linear elastic material.
    
    Args:
        stiffness(float) : Bulk modulus of the material in Pascal for SAENO Simulation (see [Steinwachs,2015])
    """
    return {'K_0': stiffness, 'D_0': 1e30, 'L_S': 1e30, 'D_S': 1e30}


def custom(K_0, D_0, L_S, D_S):
    """
    Defines the material properties for a custom nonlinear material.
    
    Args:
        K_0(float) : Bulk modulus of the material in Pascal for SAENO Simulation (see [Steinwachs,2015])
        D_0(float) : Buckling coefficient of the fibers for SAENO Simulation (see [Steinwachs,2015])
        L_S(float) : Onset of strain stiffening of the fibers (see [Steinwachs,2015])
        D_S(float) : Strain stiffening coefficient of the fibers (see [Steinwachs,2015])   
    """
    return {'K_0': K_0, 'D_0': D_0, 'L_S': L_S, 'D_S': D_S}


# The following material parameters are taken from Steinwachs et al. (2016).
# See https://www.nature.com/articles/nmeth.3685 for detailed protocols to
# replicate the bioploymer gels - 

# 1:1 mixture of Collagen G and Collagen R (0.6mg/ml)
collagen06 = {'K_0': 447, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}

# 1:1 mixture of Collagen G and Collagen R (1.2mg/ml)
collagen12 = {'K_0': 1645, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}

# 1:1 mixture of Collagen G and Collagen R (2.4mg/ml)
collagen24 = {'K_0': 5208, 'D_0': 0.0008, 'L_S': 0.0075, 'D_S': 0.033}

# Fibrin gel (4.0mg/ml)
fibrin40 = {'K_0': 2091, 'D_0': 0.002, 'L_S': 1e30, 'D_S': 1e30}

# Matrigel (10mg/ml)
matrigel100 = {'K_0': 2364, 'D_0': 1e30, 'L_S': 1e30, 'D_S': 1e30}
